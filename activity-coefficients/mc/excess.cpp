#include <faunus/faunus.h>
#include <faunus/ewald.h>
using namespace Faunus;
using namespace Faunus::Potential;

typedef Space<Geometry::Cuboid,PointParticle> Tspace;
typedef CombinedPairPotential<Coulomb, LennardJonesLB> Tpairpot;

int main() {
  InputMap mcp("excess.json");                  // open user input file
  Tspace spc(mcp);                              // simulation space

  auto pot = Energy::Nonbonded<Tspace,Tpairpot>(mcp)
    + Energy::EquilibriumEnergy<Tspace>(mcp);
  pot.setSpace(spc);                            // share space w. hamiltonian

  spc.load("state",Tspace::RESIZE);             // load old config. from disk (if any)

  // Two different Widom analysis methods
  double lB = pot.first.pairpot.first.bjerrumLength();// get bjerrum length

  Analysis::CombinedAnalysis analyze(mcp,pot,spc);
  Move::Propagator<Tspace> mv(mcp,pot,spc);

  cout << atom.info() + spc.info() + pot.info() + "\n";

  MCLoop loop(mcp);                             // class for handling mc loops
  while ( loop[0] ) {
    while ( loop[1] ) {

      mv.move();                           // move!
      analyze.sample();
  
    }                                           // end of micro loop
    cout << loop.timing();
  }                                             // end of macro loop

  FormatGRO gro;
  gro.len=spc.geo.len.x();
  gro.save("confout.gro", spc.p);

  cout << loop.info() + mv.info();// + widom1.info();
  
  std::ofstream o("move_out.json");
  o << std::setw(4) << mv.json() << endl;
}
