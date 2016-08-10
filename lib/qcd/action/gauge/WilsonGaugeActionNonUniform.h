    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/action/gauge/WilsonGaugeAction.h

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
#ifndef QCD_WILSON_GAUGE_ACTION_NONUNIFORM_H
#define QCD_WILSON_GAUGE_ACTION_NONUNIFORM_H

namespace Grid{
  namespace QCD{
    
    ////////////////////////////////////////////////////////////////////////
    // Wilson Gauge Action .. should I template the Nc etc..
    ////////////////////////////////////////////////////////////////////////
    template<class Gimpl>
    class WilsonGaugeActionNonUniform : public Action<typename Gimpl::GaugeField> {
    public:

      INHERIT_GIMPL_TYPES(Gimpl);

    private:
      std::vector<RealD> bs;
      std::vector<RealD> bt;
      GridBase *grid;
      LatticeComplex beta_tau_sft;
      LatticeComplex beta_tau;
      LatticeComplex beta_space;
    public:

      WilsonGaugeActionNonUniform(std::vector<RealD> _bs, std::vector<RealD> _bt,GridBase *_grid) : 
        bs(_bs), bt(_bt) ,
	beta_tau(_grid),
	beta_tau_sft(_grid),
	beta_space(_grid),
	grid(_grid)
      {
	assert(grid->_ndimension == Ndim );
	LatticeComplex temp(grid);

	Lattice<iScalar<vInteger> > coor(grid);
	
	LatticeCoordinate(coor,Ndim-1);

	int Ntau = grid->_fdimensions[Ndim-1];
	assert(bt.size()==Ntau);
	assert(bs.size()==Ntau);
	//	assert(bt[Ntau-1]==0.0);

	beta_tau=zero;
	for(int tau=0;tau<Ntau;tau++) {
	  temp = bt[tau];
	  beta_tau = where(coor==tau,temp,beta_tau);
	}


	beta_tau_sft = Cshift(beta_tau, Ndim-1, -1);

	beta_space=zero;
	for(int tau=0;tau<Ntau;tau++) {
	  temp = bs[tau];
	  beta_space = where(coor==tau,temp,beta_space);
	}
	
      };
      
      virtual void refresh(const GaugeField &U, GridParallelRNG& pRNG) {}; // noop as no pseudoferms
      
      virtual RealD S(const GaugeField &Umu) {

	conformable(grid,Umu._grid);

	std::vector<GaugeLinkField> U(Ndim,grid);

	for (int mu = 0; mu < Ndim; mu++) {
	  U[mu] = PeekIndex<LorentzIndex>(Umu, mu);
	}

	LatticeComplex dirPlaq(grid);
	LatticeComplex Plaq(grid);
	
	RealD OneOnNc = 1.0/Nrepresentation;

	/////////////
	// Lower dim plaquettes
	/////////////
	Plaq = zero;
	for (int mu = 1; mu < Ndim-1; mu++) {
	  for (int nu = 0; nu < mu; nu++) {
	    WilsonLoops<Gimpl>::traceDirPlaquette(dirPlaq, U, mu, nu);
	    Plaq = Plaq + (1.0 - dirPlaq*OneOnNc)*beta_space;
	  }
	}

	/////////////
	// Taumporal 
	/////////////
	{
	  int mu=Ndim-1;
	  for (int nu = 0; nu < mu; nu++) {
	    WilsonLoops<Gimpl>::traceDirPlaquette(dirPlaq, U, mu, nu);
	    Plaq = Plaq + (1.0 - dirPlaq*OneOnNc)*beta_tau;
	  }
	}

	TComplex Tp = sum(Plaq);
	Complex p   = TensorRemove(Tp);
	RealD action = p.real();
	return action;

      }

      virtual void deriv(const GaugeField &U,GaugeField & dSdU) {

	conformable(grid,U._grid);

	RealD factor = 0.5/RealD(Nrepresentation);

	GaugeLinkField Umu(grid);
	GaugeLinkField dSdU_mu(grid);
	GaugeLinkField staple(grid);

	for (int mu=0; mu < Ndim; mu++){

	  Umu = PeekIndex<LorentzIndex>(U,mu);
	  dSdU_mu = zero;

	  for(int nu=0;nu<Ndim;nu++) {

	    if (nu != mu) {

	      if (  (mu < (Ndim-1)) && (nu < (Ndim-1)) ) {

		// Spacelike case apply beta space
		WilsonLoops<Gimpl>::Staple(staple,U,mu,nu);
		staple = staple * beta_space;
		dSdU_mu += staple;

	      } else if ( mu ==  (Ndim-1) ) {

		// nu space; mu time link
		assert ( nu < (Ndim-1));
		assert ( mu == (Ndim-1));

		// mu==tau dir link deriv, nu spatial
		WilsonLoops<Gimpl>::Staple(staple,U,mu,nu);
		staple = staple * beta_tau;
		dSdU_mu += staple;


	      } else {

		assert ( mu < (Ndim-1));
		assert ( nu == (Ndim-1));

		// nu time; mu space link

		// staple forwards in tau
		WilsonLoops<Gimpl>::StapleUpper(staple,U,mu,nu);
		staple = staple * beta_tau;
		dSdU_mu += staple;

		// staple backwards in tau
		WilsonLoops<Gimpl>::StapleLower(staple,U,mu,nu);
		staple = staple * beta_tau_sft;
		dSdU_mu += staple;
	      }
	    }
	  }

	  dSdU_mu = Ta(Umu*dSdU_mu)*factor;
	  PokeIndex<LorentzIndex>(dSdU, dSdU_mu, mu);
	}
      };
    };
    
  }
}

#endif
