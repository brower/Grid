/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/action/gauge/GaugeImpl.h

Copyright (C) 2015

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

See the full license in the file "LICENSE" in the top level distribution
directory
*************************************************************************************/
/*  END LEGAL */
#ifndef GRID_QCD_GAUGE_IMPL_H
#define GRID_QCD_GAUGE_IMPL_H

namespace Grid {

namespace QCD {

////////////////////////////////////////////////////////////////////////
// Implementation dependent gauge types
////////////////////////////////////////////////////////////////////////

template <class Gimpl> class WilsonLoops;

#define INHERIT_GIMPL_TYPES(GImpl)                                             \
  typedef typename GImpl::Simd Simd;                                           \
  typedef typename GImpl::GaugeLinkField GaugeLinkField;                       \
  typedef typename GImpl::GaugeField GaugeField;                               \
  typedef typename GImpl::SiteGaugeField SiteGaugeField;                       \
  typedef typename GImpl::SiteGaugeLink SiteGaugeLink; \
  static const int Ndim = GImpl::Ndim;\
  static const int Nrepresentation = GImpl::Nrepresentation;

//
template <class S, int Ncc= Nc, int Ndd=Nd> class GaugeImplTypes {
public:
  typedef S Simd;

  static const int Ndim   = Ndd;
  static const int Nrepresentation=Ncc;

  template <typename vtype>
  using iImplGaugeLink = iScalar<iScalar<iMatrix<vtype, Nrepresentation> > >;
  template <typename vtype>
  using iImplGaugeField = iVector<iScalar<iMatrix<vtype, Nrepresentation> >, Ndim>;

  typedef iImplGaugeLink<Simd>  SiteGaugeLink;
  typedef iImplGaugeField<Simd> SiteGaugeField;

  typedef Lattice<SiteGaugeLink> GaugeLinkField; // bit ugly naming; polarised
                                                 // gauge field, lorentz... all
                                                 // ugly
  typedef Lattice<SiteGaugeField> GaugeField;

  // Move this elsewhere?
  static inline void AddGaugeLink(GaugeField &U, GaugeLinkField &W,int mu) { // U[mu] += W
    PARALLEL_FOR_LOOP
    for (auto ss = 0; ss < U._grid->oSites(); ss++) {
      U._odata[ss]._internal[mu] =
          U._odata[ss]._internal[mu] + W._odata[ss]._internal;
    }
  }
};

// Composition with smeared link, bc's etc.. probably need multiple inheritance
// Variable precision "S" and variable Nc
template <class GimplTypes> class PeriodicGaugeImpl : public GimplTypes {
public:
  INHERIT_GIMPL_TYPES(GimplTypes);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Support needed for the assembly of loops including all boundary condition
  // effects such as conjugate bcs
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////

  template <class covariant>
  static inline Lattice<covariant>
  CovShiftForward(const GaugeLinkField &Link, int mu,
                  const Lattice<covariant> &field) {
    return PeriodicBC::CovShiftForward(Link, mu, field);
  }

  template <class covariant>
  static inline Lattice<covariant>
  CovShiftBackward(const GaugeLinkField &Link, int mu,
                   const Lattice<covariant> &field) {
    return PeriodicBC::CovShiftBackward(Link, mu, field);
  }
  static inline GaugeLinkField
  CovShiftIdentityBackward(const GaugeLinkField &Link, int mu) {
    return Cshift(adj(Link), mu, -1);
  }
  static inline GaugeLinkField
  CovShiftIdentityForward(const GaugeLinkField &Link, int mu) {
    return Link;
  }
  static inline GaugeLinkField ShiftStaple(const GaugeLinkField &Link, int mu) {
    return Cshift(Link, mu, 1);
  }

  static inline bool isPeriodicGaugeField(void) { return true; }
  static void EnforceLinkBCs(GaugeField & U) { };
  static void EnforceMomentumBCs(GaugeField & U) { };
};

// Composition with smeared link, bc's etc.. probably need multiple inheritance
// Variable precision "S" and variable Nc
template <class GimplTypes> class ConjugateGaugeImpl : public GimplTypes {
public:
  INHERIT_GIMPL_TYPES(GimplTypes);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Support needed for the assembly of loops including all boundary condition
  // effects such as Gparity.
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  template <class covariant>
  static Lattice<covariant> CovShiftForward(const GaugeLinkField &Link, int mu,
                                            const Lattice<covariant> &field) {
    return ConjugateBC::CovShiftForward(Link, mu, field);
  }

  template <class covariant>
  static Lattice<covariant> CovShiftBackward(const GaugeLinkField &Link, int mu,
                                             const Lattice<covariant> &field) {
    return ConjugateBC::CovShiftBackward(Link, mu, field);
  }

  static inline GaugeLinkField
  CovShiftIdentityBackward(const GaugeLinkField &Link, int mu) {
    GridBase *grid = Link._grid;
    int Lmu = grid->GlobalDimensions()[mu] - 1;

    Lattice<iScalar<vInteger> > coor(grid);
    LatticeCoordinate(coor, mu);

    GaugeLinkField tmp(grid);
    tmp = adj(Link);
    tmp = where(coor == Lmu, conjugate(tmp), tmp);
    return Cshift(tmp, mu, -1); // moves towards positive mu
  }
  static inline GaugeLinkField
  CovShiftIdentityForward(const GaugeLinkField &Link, int mu) {
    return Link;
  }

  static inline GaugeLinkField ShiftStaple(const GaugeLinkField &Link, int mu) {
    GridBase *grid = Link._grid;
    int Lmu = grid->GlobalDimensions()[mu] - 1;

    Lattice<iScalar<vInteger> > coor(grid);
    LatticeCoordinate(coor, mu);

    GaugeLinkField tmp(grid);
    tmp = Cshift(Link, mu, 1);
    tmp = where(coor == Lmu, conjugate(tmp), tmp);
    return tmp;
  }

  static inline bool isPeriodicGaugeField(void) { return false; }
  static void EnforceLinkBCs(GaugeField & U) { };
  static void EnforceMomentumBCs(GaugeField & P) { };
};

// Composition with smeared link, bc's etc.. probably need multiple inheritance
// Variable precision "S" and variable Nc
template <class GimplTypes> class OpenGaugeImpl : public GimplTypes {
public:

  INHERIT_GIMPL_TYPES(GimplTypes);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Support needed for the assembly of loops including all boundary condition
  // effects such as conjugate bcs
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  static  void EnforceLinkBCs(GaugeField & U) { 
    Lattice<iScalar<vInteger> > coor(U._grid);
    LatticeCoordinate(coor,Ndim-1);
    
    int Ntau = U._grid->_fdimensions[Ndim-1];

    GaugeLinkField Utau(U._grid);
    GaugeLinkField zz(U._grid); zz=zero;

    Utau = PeekIndex<LorentzIndex>(U, Ndim-1);
    Utau = where(coor==(Ntau-1),zz,Utau);
    PokeIndex<LorentzIndex>(U,Utau,Ndim-1);
  };
  static void EnforceMomentumBCs(GaugeField & P) { 
    EnforceLinkBCs(P);
  };
  template <class covariant>
  static inline Lattice<covariant>
  CovShiftForward(const GaugeLinkField &Link, int mu,
                  const Lattice<covariant> &field) {
    return PeriodicBC::CovShiftForward(Link, mu, field);
  }

  template <class covariant>
  static inline Lattice<covariant>
  CovShiftBackward(const GaugeLinkField &Link, int mu,
                   const Lattice<covariant> &field) {
    return PeriodicBC::CovShiftBackward(Link, mu, field);
  }
  static inline GaugeLinkField
  CovShiftIdentityBackward(const GaugeLinkField &Link, int mu) {
    return Cshift(adj(Link), mu, -1);
  }
  static inline GaugeLinkField
  CovShiftIdentityForward(const GaugeLinkField &Link, int mu) {
    return Link;
  }
  static inline GaugeLinkField ShiftStaple(const GaugeLinkField &Link, int mu) {
    return Cshift(Link, mu, 1);
  }

  static inline bool isPeriodicGaugeField(void) { return true; }
};

typedef GaugeImplTypes<vComplex, Nc> GimplTypesR;
typedef GaugeImplTypes<vComplexF, Nc> GimplTypesF;
typedef GaugeImplTypes<vComplexD, Nc> GimplTypesD;

typedef PeriodicGaugeImpl<GimplTypesR> PeriodicGimplR; // Real.. whichever prec
typedef PeriodicGaugeImpl<GimplTypesF> PeriodicGimplF; // Float
typedef PeriodicGaugeImpl<GimplTypesD> PeriodicGimplD; // Double

typedef ConjugateGaugeImpl<GimplTypesR>
    ConjugateGimplR; // Real.. whichever prec
typedef ConjugateGaugeImpl<GimplTypesF> ConjugateGimplF; // Float
typedef ConjugateGaugeImpl<GimplTypesD> ConjugateGimplD; // Double
}
}

#endif
