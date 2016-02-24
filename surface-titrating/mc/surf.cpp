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
  Energy::Nonbonded<Tspace,Tpairpot>(mcp)       // hamiltonian
  + Energy::EquilibriumEnergy<Tspace>(mcp);

  pot.first.pairpot.second.customParameters(mcp["customlj"]);

  pot.setSpace(spc);                            // share space w. hamiltonian

  spc.load("state",Tspace::RESIZE);             // load old config. from disk (if any)

  // Two different Widom analysis methods
  double lB = 7.1;
  Analysis::LineDistribution<> rdf_ab(0.5);      // 0.1 angstrom resolution

  Move::Propagator<Tspace> mv(mcp,pot,spc);

  EnergyDrift sys;                              // class for tracking system energy drifts
  sys.init(Energy::systemEnergy(spc,pot,spc.p));// store initial total system energy

  cout << atom.info() + spc.info() + pot.info() + "\n";

  MCLoop loop(mcp);                             // class for handling mc loops
  while ( loop[0] ) {
    while ( loop[1] ) {
      sys+=mv.move();                           // move!

      if (slump() < 0.10) {}

    }                                           // end of micro loop
    sys.checkDrift(Energy::systemEnergy(spc,pot,spc.p)); // calc. energy drift
    cout << loop.timing();
  }                                             // end of macro loop

  cout << "System net charge = " << netCharge( spc.p.begin(), spc.p.end() ) << endl;
  cout << loop.info() << sys.info() << mv.info();

  FormatPQR::save("confout.pqr", spc.p, spc.geo.len);
  spc.save("state");                            // final simulation state
}
