#include <faunus/faunus.h>
#include <faunus/ewald.h>
using namespace Faunus;
using namespace Faunus::Potential;

typedef Space<Geometry::Cylinder,PointParticle> Tspace;

int main() {
  InputMap mcp("gctit.json");                   // open user input file
  Tspace spc(mcp);                              // simulation space
  auto pot =
  Energy::Nonbonded<Tspace,CoulombLJ>(mcp)      // hamiltonian
  + Energy::EquilibriumEnergy<Tspace>(mcp)
  + Energy::MassCenterConstrain<Tspace>(mcp, spc);

  pot.setSpace(spc);                            // share space w. hamiltonian

  spc.load("state",Tspace::RESIZE);             // load old config. from disk (if any)

  Analysis::LineDistribution<> rdf_ab(0.1);      // 0.1 angstrom resolution
  Analysis::ChargeMultipole mpol;
  Average<double> QQ, Q1, Q2;

  Move::Propagator<Tspace> mv(mcp,pot,spc);

  EnergyDrift sys;                              // class for tracking system energy drifts
  sys.init(Energy::systemEnergy(spc,pot,spc.p));// store initial total system energy

  cout << atom.info() + spc.info() + pot.info() + "\n";

  auto g1 = spc.findMolecules( "protein1" ).at(0);
  auto g2 = spc.findMolecules( "protein2" ).at(0);

  MCLoop loop(mcp);                             // class for handling mc loops
  while ( loop[0] ) {
    while ( loop[1] ) {
      sys+=mv.move();                           // move!
      rdf_ab( spc.geo.dist( g1->cm, g2->cm ) )++; 
      if (slump()<0.5) {
        mpol.sample(*g1, spc);
        mpol.sample(*g2, spc);
        double z1=netCharge( spc.p, *g1 );
        double z2=netCharge( spc.p, *g2 );
        Q1 += z1;
        Q2 += z2;
        QQ += z1*z2;
      }
    }                                           // end of micro loop
    sys.checkDrift(Energy::systemEnergy(spc,pot,spc.p)); // calc. energy drift
    cout << loop.timing();
    rdf_ab.save("rdf.dat");                     // g(r) - not normalized!
  }                                             // end of macro loop

  UnitTest test(mcp);                           // class for unit testing

  cout << loop.info() + sys.info() + mv.info() + test.info() + mpol.info();
  cout << "Avg. charge <Q1> <Q2> <Q1Q2>   = " << Q1 << " " << Q2 << " " << QQ << endl;
  cout << "Charge product <Q1><Q2>-<Q1Q2> = " << Q1*Q2-QQ << endl;

  FormatPQR::save("confout.pqr", spc.p);        // PQR snapshot for VMD etc.
  spc.save("state");                            // final simulation state
}
