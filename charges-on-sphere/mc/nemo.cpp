#include <faunus/faunus.h>
using namespace Faunus;
using namespace Faunus::Potential;

typedef Coulomb Tpairpot; 
typedef Geometry::SphereSurface Tgeometry; 
typedef Space<Tgeometry> Tspace;

int main() {
  InputMap mcp("nemo.json"); 
  MCLoop loop(mcp);                
  EnergyDrift sys;                 
  Tspace spc(mcp);

  auto pot = Energy::Nonbonded<Tspace,Tpairpot>(mcp);

  Move::Propagator<Tspace> mv(mcp,pot,spc);
  Analysis::RadialDistribution<> rdf(0.1);     
  //FormatXTC xtc(1000);
  
  vector<double> energy;;

  //spc.load("state");                          
  sys.init( Energy::systemEnergy(spc,pot,spc.p)  );

  cout << atom.info() + spc.info() + pot.info();

  while ( loop[0] ) {  // Markov chain 
    while ( loop[1] ) {
      sys += mv.move();
      rdf.sample( spc, atom["Na"].id, atom["Na"].id );
      //xtc.setbox( 10.0 );
      //xtc.save( "traj.xtc", spc.p );
      
      //energy.push_back(Energy::systemEnergy(spc,pot,spc.p));
    } 

    sys.checkDrift(Energy::systemEnergy(spc,pot,spc.p)); 
    cout << loop.timing();

  }

  FormatPQR::save("confout.pqr", spc.p); 
  rdf.save("rdf.dat");              
  spc.save("state");               

  // print information
  cout << loop.info() + sys.info() + mv.info();
  
  //string file = "energy.dat";
  //std::ofstream f(file.c_str());
  //if (f)
  //  for (unsigned int i = 0; i < energy.size(); i++)
  //    f << std::left << std::setw(10) << energy.at(i) << endl;

  return 0;
}
