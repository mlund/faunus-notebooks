#include <faunus/faunus.h>

using namespace Faunus;
using namespace Faunus::Potential;

typedef Space<Geometry::Sphere> Tspace;
typedef CombinedPairPotential<Coulomb, LennardJonesLB> Tpairpot;

int main() {
  InputMap mcp("titrate.json");
  Tspace spc(mcp);
  auto pot = Energy::Nonbonded<Tspace,Tpairpot>(mcp)
    + Energy::EquilibriumEnergy<Tspace>(mcp);

  pot.setSpace(spc);
  spc.load("state",Tspace::RESIZE);

  Move::Propagator<Tspace> mv(mcp,pot,spc);
  Analysis::CombinedAnalysis analysis(mcp,pot,spc);

  EnergyDrift sys;
  sys.init( Energy::systemEnergy(spc,pot,spc.p) );

  cout << atom.info() + spc.info() + pot.info()
    + textio::header("MC Simulation Begins!");

  MCLoop loop(mcp);
  while ( loop[0] ) {
    while ( loop[1] ) {
      sys+=mv.move();
      analysis.sample();
    } // end of micro loop

    sys.checkDrift( Energy::systemEnergy(spc,pot,spc.p) );
    cout << loop.timing();

  } // end of macro loop

  cout << loop.info() + sys.info() + mv.info() + analysis.info();

  FormatPQR().save("confout.pqr", spc.p);
  spc.save("state");
}
