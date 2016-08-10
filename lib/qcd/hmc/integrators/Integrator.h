    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/hmc/integrators/Integrator.h

    Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
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
//--------------------------------------------------------------------
/*! @file Integrator.h
 * @brief Classes for the Molecular Dynamics integrator
 *
 * @author Guido Cossu
 * Time-stamp: <2015-07-30 16:21:29 neo>
 */
//--------------------------------------------------------------------

#ifndef INTEGRATOR_INCLUDED
#define INTEGRATOR_INCLUDED

//class Observer;

#include <memory>

 namespace Grid{
 	namespace QCD{

 		struct IntegratorParameters{

 			int Nexp;
      int MDsteps;  // number of outer steps
      RealD trajL;  // trajectory length 
      RealD stepsize;

      IntegratorParameters(int MDsteps_, 
      	RealD trajL_=1.0,
      	int Nexp_=12):
      Nexp(Nexp_),
      MDsteps(MDsteps_),
      trajL(trajL_),
      stepsize(trajL/MDsteps)
      {
	  // empty body constructor
      };

  };

    /*! @brief Class for Molecular Dynamics management */   
  template<class Gimpl, class SmearingPolicy>
  class Integrator : public Gimpl{
  public:
    INHERIT_GIMPL_TYPES(Gimpl);

  protected:

  	typedef IntegratorParameters ParameterType;

  	IntegratorParameters Params;

  	const ActionSet<GaugeField> as;

      int levels;              //
      double t_U;              // Track time passing on each level and for U and for P
      std::vector<double> t_P; //

      GaugeField P;

      SmearingPolicy &Smearer;
      
      // Should match any legal (SU(n)) gauge field
      // Need to use this template to match Ncol to pass to SU<N> class
      void generate_momenta(GaugeField &P,GridParallelRNG& pRNG){
	GaugeLinkField Pmu(P._grid);
	Pmu = zero;
	for(int mu=0;mu<Ndim;mu++){
	  SU<Nrepresentation>::GaussianLieAlgebraMatrix(pRNG, Pmu);
	  PokeIndex<LorentzIndex>(P, Pmu, mu);
	}
      }


      //ObserverList observers; // not yet
      //      typedef std::vector<Observer*> ObserverList;
      //      void register_observers();
      //      void notify_observers();

  void update_P(GaugeField&U, int level, double ep){
  	t_P[level]+=ep;
  	update_P(P,U,level,ep);

  	std::cout<<GridLogIntegrator<<"["<<level<<"] P " << " dt "<< ep <<" : t_P "<< t_P[level] <<std::endl;
  }

  void update_P(GaugeField &Mom,GaugeField&U, int level,double ep){

    // input U actually not used... 
    for(int a=0; a<as[level].actions.size(); ++a){
      GaugeField force(U._grid);
      GaugeField& Us = Smearer.get_U(as[level].actions.at(a)->is_smeared);
      as[level].actions.at(a)->deriv(Us,force); // deriv should NOT include Ta

      std::cout<< GridLogIntegrator << "Smearing (on/off): "<<as[level].actions.at(a)->is_smeared <<std::endl;
      if (as[level].actions.at(a)->is_smeared) Smearer.smeared_force(force);
      force = Ta(force);
      std::cout<< GridLogIntegrator << "Force average: "<< norm2(force)/(U._grid->gSites()) <<std::endl;
      Mom -= force*ep;
    }
    Gimpl::EnforceMomentumBCs(Mom);
  }

  void update_U(GaugeField&U, double ep){
    update_U(P,U,ep);
    
    t_U+=ep;
    int fl = levels-1;
    std::cout<< GridLogIntegrator <<"   "<<"["<<fl<<"] U " << " dt "<< ep <<" : t_U "<< t_U <<std::endl;
  }
  void update_U(GaugeField &Mom, GaugeField&U, double ep){
    //rewrite exponential to deal automatically  with the lorentz index?
    //	GaugeLinkField Umu(U._grid);
    //	GaugeLinkField Pmu(U._grid);
    for (int mu = 0; mu < Ndim; mu++){
      auto Umu=PeekIndex<LorentzIndex>(U, mu);
      auto Pmu=PeekIndex<LorentzIndex>(Mom, mu);
      Umu = expMat(Pmu, ep, Params.Nexp)*Umu;
      ProjectOnGroup(Umu);
      PokeIndex<LorentzIndex>(U, Umu, mu);
    }
    Gimpl::EnforceLinkBCs(U);
    // Update the smeared fields, can be implemented as observer
    Smearer.set_GaugeField(U);
  }

  virtual void step (GaugeField& U,int level, int first,int last)=0;

public:

	Integrator(GridBase* grid, 
		IntegratorParameters Par,
		ActionSet<GaugeField> & Aset,
		SmearingPolicy &Sm):
	Params(Par),
	as(Aset),
	P(grid),
	levels(Aset.size()),
	Smearer(Sm)
	{
		t_P.resize(levels,0.0);
		t_U=0.0;
	// initialization of smearer delegated outside of Integrator
	};

	virtual ~Integrator(){}

      //Initialization of momenta and actions
	void refresh(GaugeField& U,GridParallelRNG &pRNG){
	  std::cout<<GridLogIntegrator<< "Integrator refresh\n";
	  generate_momenta(P,pRNG);
	  Gimpl::EnforceMomentumBCs(P);
	  for(int level=0; level< as.size(); ++level){
	    for(int actionID=0; actionID<as[level].actions.size(); ++actionID){
	      // get gauge field from the SmearingPolicy and
	      // based on the boolean is_smeared in actionID
	      GaugeField& Us = Smearer.get_U(as[level].actions.at(actionID)->is_smeared);
	      as[level].actions.at(actionID)->refresh(Us, pRNG);
	    }
	  }
	}

      // Calculate action
	RealD S(GaugeField& U){// here also U not used

		LatticeComplex Hloc(U._grid);	Hloc = zero;
	// Momenta
		for (int mu=0; mu <Ndim; mu++){
			auto Pmu = PeekIndex<LorentzIndex>(P, mu);
			Hloc -= trace(Pmu*Pmu);
		}
		Complex Hsum = sum(Hloc);

		RealD H = Hsum.real();
		RealD Hterm;
		std::cout<<GridLogMessage << "Momentum action H_p = "<< H << "\n";

	// Actions
		for(int level=0; level<as.size(); ++level){
			for(int actionID=0; actionID<as[level].actions.size(); ++actionID){
	    // get gauge field from the SmearingPolicy and
	    // based on the boolean is_smeared in actionID
				GaugeField& Us = Smearer.get_U(as[level].actions.at(actionID)->is_smeared);
				Hterm = as[level].actions.at(actionID)->S(Us);
				std::cout<<GridLogMessage << "S Level "<<level<<" term "<<actionID<<" H = "<<Hterm<<std::endl;
				H += Hterm;
			}
		}

		return H;
	}

	void integrate(GaugeField& U){

	// reset the clocks
	  t_U=0;
	  for(int level=0; level<as.size(); ++level){
	    t_P[level]=0;
	  }	
	  
	  for(int step=0; step< Params.MDsteps; ++step){   // MD step
	    int first_step = (step==0);
	    int  last_step = (step==Params.MDsteps-1);
	    this->step(U,0,first_step,last_step);
	  }

	  // Check the clocks all match on all levels
	  for(int level=0; level<as.size(); ++level){
	    assert(fabs(t_U - t_P[level])<1.0e-6); // must be the same
	    std::cout<<GridLogIntegrator<<" times["<<level<<"]= "<<t_P[level]<< " " << t_U <<std::endl;
	  }	
	  
	// and that we indeed got to the end of the trajectory
	  assert(fabs(t_U-Params.trajL) < 1.0e-6);
	  
	  
	}
};

}
}
#endif//INTEGRATOR_INCLUDED
