#include <faunus/faunus.h>
using namespace Faunus;
typedef Space<Geometry::Cuboid> Tspace;

namespace Faunus {
  namespace Potential {
    class IonIonSP3 : public Coulomb {
      private:
        string _brief() { return "Coulomb SP3"; }
        double rc1, rc1i, rc2, _lB;
      public:
        IonIonSP3(Tmjson &j, const string &sec="coulomb") : Coulomb(j,sec) { 
          name += " SP3"; 
          _lB = Coulomb(j,sec).bjerrumLength();
          rc1 = j[sec]["cutoff"] | pc::infty;
          rc1i = 1.0/rc1;
          rc2 = rc1*rc1;
        }

        template<class Tparticle>
          double operator()(const Tparticle &a, const Tparticle &b, double r2) const {
            if(r2 < rc2) {
              double r1 = sqrt(r2);
              double q = r1 * rc1i;
              double q5 = _powi<5>(q);
              //return _lB*(a.charge*b.charge/r1)*(1.0 - 1.75*q + 5.25*q5 - 7.0*q5*q + 2.5*q5*q*q);
              return _lB*(a.charge*b.charge/r1)*(1.0 - 1.75*q + q5*(5.25 + q*(2.5*q - 7.0) ) );
            }
            return 0.0;
          }

        template<class Tparticle>
          double operator()(const Tparticle &a, const Tparticle &b, const Point &r) const {
            double r2 = r.squaredNorm();
            return operator()(a,b,r2);
          }

        string info(char w) {
          using namespace textio;
          std::ostringstream o;
          o << Coulomb::info(w)
            << pad(SUB,w,"Cutoff") << rc1 << " "+angstrom << endl;
          return o.str();
        }
    };
  }
}

template<class Tspace, class Tenergy, class Tpvec>
double systemEnergy2(Tspace &spc, Tenergy &pot, const Tpvec &p) {
  pot.setSpace(spc); // ensure pot geometry is in sync with spc
  double u = pot.external(p);
  for (auto g : spc.groupList())
    u += pot.g_external(p, *g);// + pot.g_internal(p, *g);
  for (int i=0; i<(int)spc.groupList().size()-1; i++)
    for (int j=i+1; j<(int)spc.groupList().size(); j++)
      u += pot.g2g(p, *spc.groupList()[i], *spc.groupList()[j]);
  return u;
}


using namespace Potential;

typedef CombinedPairPotential<IonIonSP3, HardSphere> Tpairpot;
//typedef CombinedPairPotential<IonIonSP3, LennardJones> Tpairpot;
//typedef HardSphere Tpairpot;
//typedef R12Repulsion Tpairpot;

int main() {
  Faunus::MPI::MPIController mpi;               //init MPI

  InputMap in( textio::prefix+"temper.json" );   // read input file
  Tspace spc( in );                             //create simulation space

  Energy::Nonbonded<Tspace,Tpairpot> pot(in); 

  spc.load( textio::prefix+"state"); // load old config. from disk (if any)

  Move::ParallelTempering<Tspace> pt(pot,spc,in,mpi);//temper move
  Move::Propagator<Tspace> mv( in, pot, spc );
  Histogram<double> rdf(0.2), rdf_hshs(0.2);

  pt.setEnergyFunction(
      systemEnergy2<Tspace,Energy::Energybase<Tspace>,Tspace::ParticleVector> );

  EnergyDrift sys;                              // class for tracking system energy drifts
  sys.init(systemEnergy2(spc,pot,spc.p));// store initial total system energy

  mpi.cout << spc.info() << pot.info();

  vector<int> cations, anions, hs;
  for (size_t i=0; i<spc.p.size(); i++) {
    if (atom[spc.p[i].id].name=="POS")
      cations.push_back(i);
    if (atom[spc.p[i].id].name=="NEG")
      anions.push_back(i);
    if (atom[spc.p[i].id].name=="HS")
      hs.push_back(i);
  }

  MCLoop loop( in );                            //handle mc loops
  while ( loop[0] ) {                           //start markov chain
    while ( loop[1] ) {                        
      sys += mv.move();
    }                                          
    // sample rdf
    for (auto i : cations)
      for (auto j : anions) {
        double r = spc.geo.dist( spc.p[i], spc.p[j] );
        if (r < spc.geo.len_half.x())
          rdf(r)++;
      }
    for (auto i : hs)
      for (auto j : hs)
        if (j>i) {
          double r = spc.geo.dist( spc.p[i], spc.p[j] );
          if (r < spc.geo.len_half.x())
            rdf_hshs(r)++;
        }

    sys += pt.move();
    mpi.cout << loop.timing();                  //print progress
  }                                            

  sys.checkDrift( systemEnergy2(spc,pot,spc.p) ); // calc. energy drift

  spc.save( textio::prefix+"state");
  rdf.save( textio::prefix+"cation-anion.rdf");
  rdf_hshs.save( textio::prefix+"hs-hs.rdf");
  FormatPQR::save(textio::prefix+"confout.pqr", spc.p, spc.geo.len);

  mpi.cout << pt.info() << mv.info() << sys.info();

  double u=0;
  int cnt=0;
  Energy::Nonbonded<Tspace, HardSphere> poths(in);
  for (int i=0; i<(int)spc.groupList().size()-1; i++)
    for (int j=i+1; j<(int)spc.groupList().size(); j++) {
      double _u = poths.g2g(spc.p, *spc.groupList()[i], *spc.groupList()[j]);
      if (_u==pc::infty)
        cnt++;
      u += _u;
    }
  mpi.cout << "cnt = " << cnt << endl;

}
