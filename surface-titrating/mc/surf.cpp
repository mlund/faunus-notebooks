#include <faunus/faunus.h>
#include <faunus/ewald.h>
using namespace Faunus;
using namespace Faunus::Potential;

typedef Space<Geometry::Cuboidslit,PointParticle> Tspace;
typedef CombinedPairPotential<Coulomb,LennardJonesLB> Tpairpot;

int main() {
  InputMap mcp("surf.json");                    // open user input file
  Tspace spc(mcp);                              // simulation space
  auto pot =
    Energy::Nonbonded<Tspace,Tpairpot>(mcp)     // hamiltonian
    + Energy::EquilibriumEnergy<Tspace>(mcp);

  pot.first.pairpot.second.customParameters(mcp["customlj"]);

  pot.setSpace(spc);                            // share space w. hamiltonian

  spc.load("state",Tspace::RESIZE);             // load old config. from disk (if any)

  // set initial binding state of the lipid layer
  auto g = spc.findMolecules("lipid");
  size_t N = g.size();
  if (N==1) {
    double f = mcp["initial_degree_of_protonation"] | -1.0;
    if (f>=0) {
      cout << "Setting protonated lipid fraction to " << f << endl;
      for (int i=0; i<g[0]->size(); i++)
        if (i<f*g[0]->size())
          spc.p[i] = atom["COOH"];
        else
          spc.p[i] = atom["COO"];
      spc.trial = spc.p;
      spc.initTracker();
    }
  }

  Average<double> Zlipid;
  Histogram<double> hist_anion(0.1);
  Histogram<double> hist_cation(0.1);
  Move::Propagator<Tspace> mv(mcp,pot,spc);
  EnergyDrift sys;                              // class for tracking system energy drifts
  sys.init(Energy::systemEnergy(spc,pot,spc.p));// store initial total system energy
  auto salt = spc.findMolecules("salt");

  cout << atom.info() + spc.info() + pot.info() + "\n";

  MCLoop loop(mcp);                             // class for handling mc loops
  while ( loop[0] ) {
    while ( loop[1] ) {
      sys+=mv.move();                           // move!

      // sample ion distributions
      for (auto i : *salt[0])
        if (spc.p[i].charge<0)
          hist_anion( spc.p[i].z() )++;
        else
          if (spc.p[i].charge>0)
            hist_cation( spc.p[i].z() )++;

      Zlipid += netCharge(spc.p, *g[0]);        // average lipid charge 
    }                                           // end of micro loop
    sys.checkDrift(Energy::systemEnergy(spc,pot,spc.p)); // calc. energy drift
    cout << loop.timing();
  }                                             // end of macro loop

  cout << loop.info() << sys.info() << mv.info();

  hist_anion.save("zhist_anion.dat");
  hist_cation.save("zhist_cation.dat");
  FormatPQR::save("confout.pqr", spc.p, spc.geo.len);
  spc.save("state");                            // final simulation state

  cout << "System net charge        = " << netCharge( spc.p.begin(), spc.p.end() ) << endl;
  cout << "Average lipid charge     = " << Zlipid.avg() << endl; 
  cout << "Avg. deg. of protonation = " << (g[0]->size() + Zlipid.avg())/g[0]->size() << endl; 
}
