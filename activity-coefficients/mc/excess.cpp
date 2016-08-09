#include <faunus/faunus.h>
#include <faunus/ewald.h>
using namespace Faunus;
using namespace Faunus::Potential;

typedef Space<Geometry::Cuboid,PointParticle> Tspace;
typedef CombinedPairPotential<Coulomb, LennardJonesLB> Tpairpot;

int main() {
  InputMap mcp("excess.json");                  // open user input file
  Tspace spc(mcp);                              // simulation space
  Energy::Nonbonded<Tspace,Tpairpot> pot(mcp);  // hamiltonian
  pot.setSpace(spc);                            // share space w. hamiltonian

  spc.load("state",Tspace::RESIZE);             // load old config. from disk (if any)

  // Two different Widom analysis methods
  double lB = pot.pairpot.first.bjerrumLength();// get bjerrum length
  Analysis::Widom<PointParticle> widom1;        // widom analysis (I)
  widom1.add(spc.p);

  Move::Propagator<Tspace> mv(mcp,pot,spc);

  EnergyDrift sys;                              // class for tracking system energy drifts
  sys.init(Energy::systemEnergy(spc,pot,spc.p));// store initial total system energy

  cout << atom.info() + spc.info() + pot.info() + "\n";

  MCLoop loop(mcp);                             // class for handling mc loops
  while ( loop[0] ) {
    while ( loop[1] ) {
      sys+=mv.move();                           // move!
      widom1.sample(spc,pot,1);
    }                                           // end of micro loop
    sys.checkDrift(Energy::systemEnergy(spc,pot,spc.p)); // calc. energy drift
    cout << loop.timing();
  }                                             // end of macro loop

  FormatPQR::save("confout.pqr", spc.p);        // PQR snapshot for VMD etc.
  spc.save("state");                            // final simulation state

  cout << loop.info() + sys.info() + mv.info() + widom1.info();
  
  std::ofstream o("move_out.json");
  o << std::setw(4) << mv.json() << endl;
}
