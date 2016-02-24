#include <faunus/faunus.h>

using namespace Faunus;

typedef Space<Geometry::Cuboidslit> Tspace;

// Split large oligomer into smaller ones
class Oligomers {
  private:
    int n; // residues per peptide
  public:
    Tspace::ParticleVector v;
    // @param file AAM file with loooong oligomer
    // @param peplen Length of each peptide/monomer
    Oligomers(string file, int peplen) : n(peplen) {
      FormatAAM::load(file, v);
      if (v.size() % n != 0 )
        std::cerr << "aam NOT file ok...N=" << v.size() << " is not divisible w. " << n << "\n";

      //set ntr and ctr
      for (int p=0; p<int(v.size())/n; p++) {
        v.at( p*n+0  )  = atom["NTR"];
        v.at( p*n+n-1 ) = atom["CTR"];
      }
      //set charges
      for (auto &i : v) {
        if (  atom[i.id].name=="ASP" ) i.charge = -1;
        if (  atom[i.id].name=="GLU" ) i.charge = -1;
        if (  atom[i.id].name=="CTR" ) i.charge = -1;
        if (  atom[i.id].name=="NTR" ) i.charge = 1;
        if (  atom[i.id].name=="LYS" ) i.charge = 1;
      }
    }

    // return maximum n-mer that can be extracted
    int maxMer() { return int(v.size())/n; }

    // @brief Extract a suboligomer
    // @param begpep Number of first peptide
    // @param numpep Total number of peptides to extract
    void extract(int begpep, int numpep, Tspace::ParticleVector &dest) {
      assert (size_t(numpep*n) <= v.size());
      int j=0;
      int currpep=begpep;
      dest.resize( numpep*n );
      while (j < numpep ) {
        cout << currpep << " ";
        for (int i=0; i<n; i++) {
          dest.at(j*n+i) = v.at(currpep*n+i);
        }
        if ((currpep/5)%2 == 0) 
          currpep+=5;
        else {
          if ((currpep+1)%5 == 0)
            currpep+=1;
          else
            currpep-=4;
        }
        j++;
      }
      cout << endl;
      /* for (int i=begpep*n; i<(begpep+numpep)*n; i++)
         dest.at(j++) = v.at(i);
         */
    }
};

namespace Faunus {
  namespace Move {
    template<class Tspace>
      class Swapper : public Movebase<Tspace> {
        private:
          vector<unsigned int> hist;
          typename Tspace::ParticleVector v;
          typedef Movebase<Tspace> base;
          using base::spc;
          using base::pot;

          string _info() {
            std::ostringstream o;
            o << "Max mer     = " << oligo.maxMer() << endl;
            o << "Max atoms   = " << oligo.v.size() << endl;
            return o.str();
          }

          void _trialMove() {
            //int numpep = 1 + rand() % oligo.maxMer();
            int numpep = 1 + rand() % 32;
            oligo.extract(0,numpep,v);

            // find random location {and orientation..todo}
            Point a;
            spc->geo.randompos(a);
            Geometry::cm2origo(spc->geo, v);
            Geometry::translate(spc->geo, v, a);

            spc->trial.resize( v.size() );
            for (size_t i=0; i<v.size(); i++)
              spc->trial[i] = v[i];
          }

          double _energyChange() {
            double uold=0, unew=0;
            for (size_t i=0; i<spc->trial.size(); i++)
              if (spc->geo.collision(spc->trial[i]))
                return pc::infty;
              else unew += pot->i_external(spc->trial,i);
            for (size_t i=0; i<spc->p.size(); i++) {
              assert( spc->geo.collision(spc->p[i])==false && "accepted config collides w. container" );
              uold += pot->i_external(spc->p,i);
            }
            return unew-uold;
          } 

          void _acceptMove() {
            spc->p.resize( spc->trial.size()   );
            for (size_t i=0; i<spc->trial.size(); i++)
              spc->p[i] = spc->trial[i];
            for (auto i : spc->groupList())
              i->resize(spc->p.size());
            hist.at( spc->p.size() / 26 -1 )++;
          }

          void _rejectMove() {
            spc->trial.resize( spc->p.size()  );
            for (size_t i=0; i<spc->trial.size(); i++)
              spc->trial[i] = spc->p[i];
            hist.at( spc->p.size() / 26 -1)++;
          } 

          Oligomers oligo;

        public:
          Swapper(InputMap &in, Energy::Energybase<Tspace> &e, Tspace &s, string pfx="swapper_")
            : base(e,s,pfx), oligo(in.get<string>("oligomer",""), 26)
            {
              base::title="Oligomer Swap";
              hist.resize( oligo.maxMer(), 0 );
            }

          ~Swapper() {
            std::ofstream f(base::prefix+"merhist.dat");
            int cnt=1;
            for (auto i : hist)
              f << cnt++ << " " << i << "\n"; 
          }
      };
  }//namespace
}//namespace

int main() {
  InputMap mcp("r2s.input");         // Open input parameter file
  MCLoop loop(mcp);                  // handle mc loops
  EnergyDrift sys;                   // track system energy drifts
  Tspace spc(mcp);                   // Simulation space (all particles and group info)

  spc.reserve(6000);

  Energy::ExternalPotential<Tspace,Potential::GouyChapman<double,true> > pot(mcp);
  pot.expot.setSurfPositionZ( &spc.geo.len_half.z() ); // Pass position of GC surface

  // Load and add polymer to Space
  int numpep = mcp.get<int>("numpep", 1);
  string file = mcp.get<string>("oligomer", "");
  Tspace::ParticleVector v;
  Oligomers oligo(file, 26);
  oligo.extract(0,numpep,v);

  /*
  for (int i=1 ; i<32; i++) {
    oligo.extract(0,i,v);
    FormatPQR::save("numpep" + std::to_string(i) + ".pqr", v);
  }
  return 0;
*/
  Geometry::FindSpace().find(spc.geo,spc.p,v);// find empty spot in particle vector
  Group pol = spc.insert(v);                  // Insert into Space and return matching group
  pol.name="polymer";                         // Give polymer arbitrary name
  spc.enroll(pol);                            // Enroll polymer in Space

  for (auto &i : spc.p)
    cout << atom[i.id].name << " " << i.charge << "\n";

  // MC moves
  Move::TranslateRotate<Tspace> gmv(mcp,pot,spc);
  Move::Swapper<Tspace> swp(mcp,pot,spc);

  Analysis::LineDistribution<> surfmapall, cmdist;    // monomer-surface histogram
  spc.load("state");                                  // Load start configuration, if any
  sys.init( Energy::systemEnergy(spc,pot,spc.p) );    // Store total system energy

  cout << spc.info() + pot.info() + pol.info()
    + swp.info() + textio::header("MC Simulation Begins!");

  while ( loop.macroCnt() ) {  // Markov chain 
    while ( loop.microCnt() ) {
      int i=slp_global.rand() % 1;
      switch (i) {
        case 0: // translate and rotate polymer
          gmv.setGroup(pol);
          sys += gmv.move(); 
          cmdist( pot.expot.surfDist( pol.cm ) )++;
          break;
        case 10:
          double rnd = slp_global(); // [0:1[
          if (rnd>0.9)
            sys += swp.move();
          break;
      }

      double rnd = slp_global(); // [0:1[
      if (rnd<0.05)
        for (auto i : pol)
          surfmapall( pot.expot.surfDist( spc.p[i] ) )++;  // sample monomer distribution

    } // end of micro loop

    sys.checkDrift( Energy::systemEnergy(spc,pot,spc.p) ); // re-calc system energy and detect drift
    cout << loop.timing();                                 // print timing and ETA information


    spc.save("state");               // save final state of simulation (positions etc)
    surfmapall.save("surfall.dat");  // save monomer-surface distribution
    cmdist.save("cmdist.dat");
    FormatPQR::save("confout.pqr", spc.p);  // save PQR file

  } // end of macro loop

  cout << loop.info() + sys.info() + gmv.info() + spc.info() + swp.info() ;
}
