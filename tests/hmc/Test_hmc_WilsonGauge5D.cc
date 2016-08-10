    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_hmc_WilsonGauge.cc

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: neo <cossu@post.kek.jp>
Author: paboyle <paboyle@ph.ed.ac.uk>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

    See the full license in the file "LICENSE" in the top level distribution directory
    *************************************************************************************/
    /*  END LEGAL */
#include <Grid/Grid.h>

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

namespace Grid { 
  namespace QCD { 


typedef GaugeImplTypes<vComplex,  Nc,5> GimplTypes5DR;
typedef GaugeImplTypes<vComplexF, Nc,5> GimplTypes5DF;
typedef GaugeImplTypes<vComplexD, Nc,5> GimplTypes5DD;

typedef PeriodicGaugeImpl<GimplTypes5DR> PeriodicGimpl5DR; // Real.. whichever prec
typedef PeriodicGaugeImpl<GimplTypes5DF> PeriodicGimpl5DF; // Float
typedef PeriodicGaugeImpl<GimplTypes5DD> PeriodicGimpl5DD; // Double


typedef WilsonGaugeAction<PeriodicGimpl5DR>        WilsonGaugeAction5DR;
typedef WilsonGaugeAction<PeriodicGimpl5DF>        WilsonGaugeAction5DF;
typedef WilsonGaugeAction<PeriodicGimpl5DD>        WilsonGaugeAction5DD;

class HmcRunner : public NerscHmcRunnerTemplate<PeriodicGimpl5DR> {
public:

  void BuildTheAction (int argc, char **argv)
  {
    const int Ndim=5;

    UGrid   =  new GridCartesian(GridDefaultLatt(),GridDefaultSimd(Ndim,vComplex::Nsimd()),GridDefaultMpi());

    // 5d grid expected
    assert(UGrid->_ndimension ==5);

    // temporarily need a gauge field
    typedef typename WilsonGaugeAction5DR::GaugeLinkField GaugeLinkField;
    typedef typename WilsonGaugeAction5DR::GaugeField GaugeField;

    GaugeField U(UGrid);

    // Gauge action
    WilsonGaugeAction5DR Waction(5.6);

    //Collect actions
    ActionLevel<GaugeField> Level1(1);
    Level1.push_back(&Waction);
    TheAction.push_back(Level1);

    Run(argc,argv);
  };

};

}}

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  int threads = GridThread::GetThreads();
  std::cout<<GridLogMessage << "Grid is setup to use "<<threads<<" threads"<<std::endl;

  HmcRunner TheHMC;
  
  TheHMC.BuildTheAction(argc,argv);

}

