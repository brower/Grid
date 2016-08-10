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

#if 1
typedef OpenGaugeImpl<GimplTypes5DR> Gimpl5DR; // Real.. whichever prec
typedef OpenGaugeImpl<GimplTypes5DF> Gimpl5DF; // Float
typedef OpenGaugeImpl<GimplTypes5DD> Gimpl5DD; // Double


typedef WilsonGaugeActionNonUniform<Gimpl5DR>        WilsonGaugeAction5DR;
typedef WilsonGaugeActionNonUniform<Gimpl5DF>        WilsonGaugeAction5DF;
typedef WilsonGaugeActionNonUniform<Gimpl5DD>        WilsonGaugeAction5DD;
#else
typedef PeriodicGaugeImpl<GimplTypes5DR> Gimpl5DR; // Real.. whichever prec
typedef PeriodicGaugeImpl<GimplTypes5DF> Gimpl5DF; // Float
typedef PeriodicGaugeImpl<GimplTypes5DD> Gimpl5DD; // Double


typedef WilsonGaugeActionNonUniform<Gimpl5DR>        WilsonGaugeAction5DR;
typedef WilsonGaugeActionNonUniform<Gimpl5DF>        WilsonGaugeAction5DF;
typedef WilsonGaugeActionNonUniform<Gimpl5DD>        WilsonGaugeAction5DD;
#endif

class HmcRunner : public NerscHmcRunnerTemplate<Gimpl5DR> {
public:

  void BuildTheAction (int argc, char **argv)
  {
    const int Ndim=5;
    
    // 5d grid expected
    UGrid   =  new GridCartesian(GridDefaultLatt(),GridDefaultSimd(Ndim,vComplex::Nsimd()),GridDefaultMpi());
    int Ntau= UGrid->_fdimensions[Ndim-1];
    assert(UGrid->_ndimension ==Ndim);

    // temporarily need a gauge field
    typedef typename WilsonGaugeAction5DR::GaugeLinkField GaugeLinkField;
    typedef typename WilsonGaugeAction5DR::GaugeField GaugeField;

    GaugeField U(UGrid);

    // Gauge action
    std::cout << "Ntau is "<<Ntau;

    std::vector<RealD> beta_4(Ntau,5.6);
    std::vector<RealD> beta_tau(Ntau,6.0);
    beta_tau[Ntau-1]=0.0;

    WilsonGaugeAction5DR Waction(beta_4,beta_tau,UGrid);

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

