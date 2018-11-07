#ifndef ALGOIM_BLITZINC_HPP
#define ALGOIM_BLITZINC_HPP

// Management of Algoim's use of blitz++, mainly to inject blitz::TinyVector into the Algoim namespace
// and to enable vector expressions when using QD for double-double or quadruple-double precision.

#include <blitz/array.h>
// #include <blitz/tinyvec-et.h> // comment out if using blitz-0.10
#include "algoim_real.hpp"

namespace Algoim
{
    // Algoim makes heavy use of blitz::TinyVector, but try not to pollute global namespace by using all of blitz
    using ::blitz::TinyVector;
}

// std::equal_to for blitz::TinyVector
namespace std
{
    template<typename T, int N>
    struct equal_to<blitz::TinyVector<T,N>>
    {
        bool operator()(const blitz::TinyVector<T,N>& x, const blitz::TinyVector<T,N>& y) const
        {
            return blitz::all(x == y);
        }
    };
}

#ifdef ALGOIM_HPREAL

// Enable blitz vector expressions for Real

BZ_NAMESPACE(blitz)

// Should be fine to leave this out, but specifying it could inform blitz that a trivial
// constructor is available.
//BZDECLNUMTRAIT(Real, Real, Real, Real, Real);

// Use this to inform blitz that all types can be promoted to Real as it has more precision
//BZ_DECLARE_PRECISION(Real,710)

// The following has been copied from vecbops.cc, isolating those routines containing
// 'long double' and replaced with Real.

// Vector<P_numtype1> + Real
template<class P_numtype1>
inline
_bz_VecExpr<_bz_VecExprOp<VectorIterConst<P_numtype1>, 
      _bz_VecExprConstant<Real>,
      _bz_Add<P_numtype1, Real > > >
operator+(const Vector<P_numtype1>& d1, 
      Real d2)
{
    typedef _bz_VecExprOp<VectorIterConst<P_numtype1>, 
      _bz_VecExprConstant<Real>, 
      _bz_Add<P_numtype1, Real> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(d1.beginFast(), 
      _bz_VecExprConstant<Real>(d2)));
}

// _bz_VecExpr<P_expr1> + Real
template<class P_expr1>
inline
_bz_VecExpr<_bz_VecExprOp<_bz_VecExpr<P_expr1>, 
      _bz_VecExprConstant<Real>,
      _bz_Add<typename P_expr1::T_numtype, Real > > >
operator+(_bz_VecExpr<P_expr1> d1, 
      Real d2)
{
    typedef _bz_VecExprOp<_bz_VecExpr<P_expr1>, 
      _bz_VecExprConstant<Real>, 
      _bz_Add<typename P_expr1::T_numtype, Real> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(d1, 
      _bz_VecExprConstant<Real>(d2)));
}

// VectorPick<P_numtype1> + Real
template<class P_numtype1>
inline
_bz_VecExpr<_bz_VecExprOp<VectorPickIterConst<P_numtype1>, 
      _bz_VecExprConstant<Real>,
      _bz_Add<P_numtype1, Real > > >
operator+(const VectorPick<P_numtype1>& d1, 
      Real d2)
{
    typedef _bz_VecExprOp<VectorPickIterConst<P_numtype1>, 
      _bz_VecExprConstant<Real>, 
      _bz_Add<P_numtype1, Real> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(d1.beginFast(), 
      _bz_VecExprConstant<Real>(d2)));
}

// Range + Real

inline
_bz_VecExpr<_bz_VecExprOp<Range, 
      _bz_VecExprConstant<Real>,
      _bz_Add<int, Real > > >
operator+(Range d1, 
      Real d2)
{
    typedef _bz_VecExprOp<Range, 
      _bz_VecExprConstant<Real>, 
      _bz_Add<int, Real> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(d1, 
      _bz_VecExprConstant<Real>(d2)));
}

// TinyVector<P_numtype1, N_length1> + Real
template<class P_numtype1, int N_length1>
inline
_bz_VecExpr<_bz_VecExprOp<TinyVectorIterConst<P_numtype1, N_length1>, 
      _bz_VecExprConstant<Real>,
      _bz_Add<P_numtype1, Real > > >
operator+(const TinyVector<P_numtype1, N_length1>& d1, 
      Real d2)
{
    typedef _bz_VecExprOp<TinyVectorIterConst<P_numtype1, N_length1>, 
      _bz_VecExprConstant<Real>, 
      _bz_Add<P_numtype1, Real> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(d1.beginFast(), 
      _bz_VecExprConstant<Real>(d2)));
}

// Real + Vector<P_numtype2>
template<class P_numtype2>
inline
_bz_VecExpr<_bz_VecExprOp<_bz_VecExprConstant<Real>, 
      VectorIterConst<P_numtype2>,
      _bz_Add<Real, P_numtype2 > > >
operator+(Real d1, 
      const Vector<P_numtype2>& d2)
{
    typedef _bz_VecExprOp<_bz_VecExprConstant<Real>, 
      VectorIterConst<P_numtype2>, 
      _bz_Add<Real, P_numtype2> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(_bz_VecExprConstant<Real>(d1), 
      d2.beginFast()));
}

// Real + _bz_VecExpr<P_expr2>
template<class P_expr2>
inline
_bz_VecExpr<_bz_VecExprOp<_bz_VecExprConstant<Real>, 
      _bz_VecExpr<P_expr2>,
      _bz_Add<Real, typename P_expr2::T_numtype > > >
operator+(Real d1, 
      _bz_VecExpr<P_expr2> d2)
{
    typedef _bz_VecExprOp<_bz_VecExprConstant<Real>, 
      _bz_VecExpr<P_expr2>, 
      _bz_Add<Real, typename P_expr2::T_numtype> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(_bz_VecExprConstant<Real>(d1), 
      d2));
}

// Real + VectorPick<P_numtype2>
template<class P_numtype2>
inline
_bz_VecExpr<_bz_VecExprOp<_bz_VecExprConstant<Real>, 
      VectorPickIterConst<P_numtype2>,
      _bz_Add<Real, P_numtype2 > > >
operator+(Real d1, 
      const VectorPick<P_numtype2>& d2)
{
    typedef _bz_VecExprOp<_bz_VecExprConstant<Real>, 
      VectorPickIterConst<P_numtype2>, 
      _bz_Add<Real, P_numtype2> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(_bz_VecExprConstant<Real>(d1), 
      d2.beginFast()));
}

// Real + Range

inline
_bz_VecExpr<_bz_VecExprOp<_bz_VecExprConstant<Real>, 
      Range,
      _bz_Add<Real, int > > >
operator+(Real d1, 
      Range d2)
{
    typedef _bz_VecExprOp<_bz_VecExprConstant<Real>, 
      Range, 
      _bz_Add<Real, int> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(_bz_VecExprConstant<Real>(d1), 
      d2));
}

// Real + TinyVector<P_numtype2, N_length2>
template<class P_numtype2, int N_length2>
inline
_bz_VecExpr<_bz_VecExprOp<_bz_VecExprConstant<Real>, 
      TinyVectorIterConst<P_numtype2, N_length2>,
      _bz_Add<Real, P_numtype2 > > >
operator+(Real d1, 
      const TinyVector<P_numtype2, N_length2>& d2)
{
    typedef _bz_VecExprOp<_bz_VecExprConstant<Real>, 
      TinyVectorIterConst<P_numtype2, N_length2>, 
      _bz_Add<Real, P_numtype2> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(_bz_VecExprConstant<Real>(d1), 
      d2.beginFast()));
}

// Vector<P_numtype1> - Real
template<class P_numtype1>
inline
_bz_VecExpr<_bz_VecExprOp<VectorIterConst<P_numtype1>, 
      _bz_VecExprConstant<Real>,
      _bz_Subtract<P_numtype1, Real > > >
operator-(const Vector<P_numtype1>& d1, 
      Real d2)
{
    typedef _bz_VecExprOp<VectorIterConst<P_numtype1>, 
      _bz_VecExprConstant<Real>, 
      _bz_Subtract<P_numtype1, Real> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(d1.beginFast(), 
      _bz_VecExprConstant<Real>(d2)));
}

// _bz_VecExpr<P_expr1> - Real
template<class P_expr1>
inline
_bz_VecExpr<_bz_VecExprOp<_bz_VecExpr<P_expr1>, 
      _bz_VecExprConstant<Real>,
      _bz_Subtract<typename P_expr1::T_numtype, Real > > >
operator-(_bz_VecExpr<P_expr1> d1, 
      Real d2)
{
    typedef _bz_VecExprOp<_bz_VecExpr<P_expr1>, 
      _bz_VecExprConstant<Real>, 
      _bz_Subtract<typename P_expr1::T_numtype, Real> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(d1, 
      _bz_VecExprConstant<Real>(d2)));
}

// VectorPick<P_numtype1> - Real
template<class P_numtype1>
inline
_bz_VecExpr<_bz_VecExprOp<VectorPickIterConst<P_numtype1>, 
      _bz_VecExprConstant<Real>,
      _bz_Subtract<P_numtype1, Real > > >
operator-(const VectorPick<P_numtype1>& d1, 
      Real d2)
{
    typedef _bz_VecExprOp<VectorPickIterConst<P_numtype1>, 
      _bz_VecExprConstant<Real>, 
      _bz_Subtract<P_numtype1, Real> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(d1.beginFast(), 
      _bz_VecExprConstant<Real>(d2)));
}

// Range - Real

inline
_bz_VecExpr<_bz_VecExprOp<Range, 
      _bz_VecExprConstant<Real>,
      _bz_Subtract<int, Real > > >
operator-(Range d1, 
      Real d2)
{
    typedef _bz_VecExprOp<Range, 
      _bz_VecExprConstant<Real>, 
      _bz_Subtract<int, Real> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(d1, 
      _bz_VecExprConstant<Real>(d2)));
}

// TinyVector<P_numtype1, N_length1> - Real
template<class P_numtype1, int N_length1>
inline
_bz_VecExpr<_bz_VecExprOp<TinyVectorIterConst<P_numtype1, N_length1>, 
      _bz_VecExprConstant<Real>,
      _bz_Subtract<P_numtype1, Real > > >
operator-(const TinyVector<P_numtype1, N_length1>& d1, 
      Real d2)
{
    typedef _bz_VecExprOp<TinyVectorIterConst<P_numtype1, N_length1>, 
      _bz_VecExprConstant<Real>, 
      _bz_Subtract<P_numtype1, Real> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(d1.beginFast(), 
      _bz_VecExprConstant<Real>(d2)));
}

// Real - Vector<P_numtype2>
template<class P_numtype2>
inline
_bz_VecExpr<_bz_VecExprOp<_bz_VecExprConstant<Real>, 
      VectorIterConst<P_numtype2>,
      _bz_Subtract<Real, P_numtype2 > > >
operator-(Real d1, 
      const Vector<P_numtype2>& d2)
{
    typedef _bz_VecExprOp<_bz_VecExprConstant<Real>, 
      VectorIterConst<P_numtype2>, 
      _bz_Subtract<Real, P_numtype2> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(_bz_VecExprConstant<Real>(d1), 
      d2.beginFast()));
}

// Real - _bz_VecExpr<P_expr2>
template<class P_expr2>
inline
_bz_VecExpr<_bz_VecExprOp<_bz_VecExprConstant<Real>, 
      _bz_VecExpr<P_expr2>,
      _bz_Subtract<Real, typename P_expr2::T_numtype > > >
operator-(Real d1, 
      _bz_VecExpr<P_expr2> d2)
{
    typedef _bz_VecExprOp<_bz_VecExprConstant<Real>, 
      _bz_VecExpr<P_expr2>, 
      _bz_Subtract<Real, typename P_expr2::T_numtype> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(_bz_VecExprConstant<Real>(d1), 
      d2));
}

// Real - VectorPick<P_numtype2>
template<class P_numtype2>
inline
_bz_VecExpr<_bz_VecExprOp<_bz_VecExprConstant<Real>, 
      VectorPickIterConst<P_numtype2>,
      _bz_Subtract<Real, P_numtype2 > > >
operator-(Real d1, 
      const VectorPick<P_numtype2>& d2)
{
    typedef _bz_VecExprOp<_bz_VecExprConstant<Real>, 
      VectorPickIterConst<P_numtype2>, 
      _bz_Subtract<Real, P_numtype2> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(_bz_VecExprConstant<Real>(d1), 
      d2.beginFast()));
}

// Real - Range

inline
_bz_VecExpr<_bz_VecExprOp<_bz_VecExprConstant<Real>, 
      Range,
      _bz_Subtract<Real, int > > >
operator-(Real d1, 
      Range d2)
{
    typedef _bz_VecExprOp<_bz_VecExprConstant<Real>, 
      Range, 
      _bz_Subtract<Real, int> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(_bz_VecExprConstant<Real>(d1), 
      d2));
}

// Real - TinyVector<P_numtype2, N_length2>
template<class P_numtype2, int N_length2>
inline
_bz_VecExpr<_bz_VecExprOp<_bz_VecExprConstant<Real>, 
      TinyVectorIterConst<P_numtype2, N_length2>,
      _bz_Subtract<Real, P_numtype2 > > >
operator-(Real d1, 
      const TinyVector<P_numtype2, N_length2>& d2)
{
    typedef _bz_VecExprOp<_bz_VecExprConstant<Real>, 
      TinyVectorIterConst<P_numtype2, N_length2>, 
      _bz_Subtract<Real, P_numtype2> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(_bz_VecExprConstant<Real>(d1), 
      d2.beginFast()));
}

// Vector<P_numtype1> * Real
template<class P_numtype1>
inline
_bz_VecExpr<_bz_VecExprOp<VectorIterConst<P_numtype1>, 
      _bz_VecExprConstant<Real>,
      _bz_Multiply<P_numtype1, Real > > >
operator*(const Vector<P_numtype1>& d1, 
      Real d2)
{
    typedef _bz_VecExprOp<VectorIterConst<P_numtype1>, 
      _bz_VecExprConstant<Real>, 
      _bz_Multiply<P_numtype1, Real> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(d1.beginFast(), 
      _bz_VecExprConstant<Real>(d2)));
}

// _bz_VecExpr<P_expr1> * Real
template<class P_expr1>
inline
_bz_VecExpr<_bz_VecExprOp<_bz_VecExpr<P_expr1>, 
      _bz_VecExprConstant<Real>,
      _bz_Multiply<typename P_expr1::T_numtype, Real > > >
operator*(_bz_VecExpr<P_expr1> d1, 
      Real d2)
{
    typedef _bz_VecExprOp<_bz_VecExpr<P_expr1>, 
      _bz_VecExprConstant<Real>, 
      _bz_Multiply<typename P_expr1::T_numtype, Real> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(d1, 
      _bz_VecExprConstant<Real>(d2)));
}

// VectorPick<P_numtype1> * Real
template<class P_numtype1>
inline
_bz_VecExpr<_bz_VecExprOp<VectorPickIterConst<P_numtype1>, 
      _bz_VecExprConstant<Real>,
      _bz_Multiply<P_numtype1, Real > > >
operator*(const VectorPick<P_numtype1>& d1, 
      Real d2)
{
    typedef _bz_VecExprOp<VectorPickIterConst<P_numtype1>, 
      _bz_VecExprConstant<Real>, 
      _bz_Multiply<P_numtype1, Real> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(d1.beginFast(), 
      _bz_VecExprConstant<Real>(d2)));
}

// Range * Real

inline
_bz_VecExpr<_bz_VecExprOp<Range, 
      _bz_VecExprConstant<Real>,
      _bz_Multiply<int, Real > > >
operator*(Range d1, 
      Real d2)
{
    typedef _bz_VecExprOp<Range, 
      _bz_VecExprConstant<Real>, 
      _bz_Multiply<int, Real> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(d1, 
      _bz_VecExprConstant<Real>(d2)));
}

// TinyVector<P_numtype1, N_length1> * Real
template<class P_numtype1, int N_length1>
inline
_bz_VecExpr<_bz_VecExprOp<TinyVectorIterConst<P_numtype1, N_length1>, 
      _bz_VecExprConstant<Real>,
      _bz_Multiply<P_numtype1, Real > > >
operator*(const TinyVector<P_numtype1, N_length1>& d1, 
      Real d2)
{
    typedef _bz_VecExprOp<TinyVectorIterConst<P_numtype1, N_length1>, 
      _bz_VecExprConstant<Real>, 
      _bz_Multiply<P_numtype1, Real> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(d1.beginFast(), 
      _bz_VecExprConstant<Real>(d2)));
}

// Real * Vector<P_numtype2>
template<class P_numtype2>
inline
_bz_VecExpr<_bz_VecExprOp<_bz_VecExprConstant<Real>, 
      VectorIterConst<P_numtype2>,
      _bz_Multiply<Real, P_numtype2 > > >
operator*(Real d1, 
      const Vector<P_numtype2>& d2)
{
    typedef _bz_VecExprOp<_bz_VecExprConstant<Real>, 
      VectorIterConst<P_numtype2>, 
      _bz_Multiply<Real, P_numtype2> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(_bz_VecExprConstant<Real>(d1), 
      d2.beginFast()));
}

// Real * _bz_VecExpr<P_expr2>
template<class P_expr2>
inline
_bz_VecExpr<_bz_VecExprOp<_bz_VecExprConstant<Real>, 
      _bz_VecExpr<P_expr2>,
      _bz_Multiply<Real, typename P_expr2::T_numtype > > >
operator*(Real d1, 
      _bz_VecExpr<P_expr2> d2)
{
    typedef _bz_VecExprOp<_bz_VecExprConstant<Real>, 
      _bz_VecExpr<P_expr2>, 
      _bz_Multiply<Real, typename P_expr2::T_numtype> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(_bz_VecExprConstant<Real>(d1), 
      d2));
}

// Real * VectorPick<P_numtype2>
template<class P_numtype2>
inline
_bz_VecExpr<_bz_VecExprOp<_bz_VecExprConstant<Real>, 
      VectorPickIterConst<P_numtype2>,
      _bz_Multiply<Real, P_numtype2 > > >
operator*(Real d1, 
      const VectorPick<P_numtype2>& d2)
{
    typedef _bz_VecExprOp<_bz_VecExprConstant<Real>, 
      VectorPickIterConst<P_numtype2>, 
      _bz_Multiply<Real, P_numtype2> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(_bz_VecExprConstant<Real>(d1), 
      d2.beginFast()));
}

// Real * Range

inline
_bz_VecExpr<_bz_VecExprOp<_bz_VecExprConstant<Real>, 
      Range,
      _bz_Multiply<Real, int > > >
operator*(Real d1, 
      Range d2)
{
    typedef _bz_VecExprOp<_bz_VecExprConstant<Real>, 
      Range, 
      _bz_Multiply<Real, int> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(_bz_VecExprConstant<Real>(d1), 
      d2));
}

// Real * TinyVector<P_numtype2, N_length2>
template<class P_numtype2, int N_length2>
inline
_bz_VecExpr<_bz_VecExprOp<_bz_VecExprConstant<Real>, 
      TinyVectorIterConst<P_numtype2, N_length2>,
      _bz_Multiply<Real, P_numtype2 > > >
operator*(Real d1, 
      const TinyVector<P_numtype2, N_length2>& d2)
{
    typedef _bz_VecExprOp<_bz_VecExprConstant<Real>, 
      TinyVectorIterConst<P_numtype2, N_length2>, 
      _bz_Multiply<Real, P_numtype2> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(_bz_VecExprConstant<Real>(d1), 
      d2.beginFast()));
}

// Vector<P_numtype1> / Real
template<class P_numtype1>
inline
_bz_VecExpr<_bz_VecExprOp<VectorIterConst<P_numtype1>, 
      _bz_VecExprConstant<Real>,
      _bz_Divide<P_numtype1, Real > > >
operator/(const Vector<P_numtype1>& d1, 
      Real d2)
{
    typedef _bz_VecExprOp<VectorIterConst<P_numtype1>, 
      _bz_VecExprConstant<Real>, 
      _bz_Divide<P_numtype1, Real> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(d1.beginFast(), 
      _bz_VecExprConstant<Real>(d2)));
}

// _bz_VecExpr<P_expr1> / Real
template<class P_expr1>
inline
_bz_VecExpr<_bz_VecExprOp<_bz_VecExpr<P_expr1>, 
      _bz_VecExprConstant<Real>,
      _bz_Divide<typename P_expr1::T_numtype, Real > > >
operator/(_bz_VecExpr<P_expr1> d1, 
      Real d2)
{
    typedef _bz_VecExprOp<_bz_VecExpr<P_expr1>, 
      _bz_VecExprConstant<Real>, 
      _bz_Divide<typename P_expr1::T_numtype, Real> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(d1, 
      _bz_VecExprConstant<Real>(d2)));
}

// VectorPick<P_numtype1> / Real
template<class P_numtype1>
inline
_bz_VecExpr<_bz_VecExprOp<VectorPickIterConst<P_numtype1>, 
      _bz_VecExprConstant<Real>,
      _bz_Divide<P_numtype1, Real > > >
operator/(const VectorPick<P_numtype1>& d1, 
      Real d2)
{
    typedef _bz_VecExprOp<VectorPickIterConst<P_numtype1>, 
      _bz_VecExprConstant<Real>, 
      _bz_Divide<P_numtype1, Real> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(d1.beginFast(), 
      _bz_VecExprConstant<Real>(d2)));
}

// Range / Real

inline
_bz_VecExpr<_bz_VecExprOp<Range, 
      _bz_VecExprConstant<Real>,
      _bz_Divide<int, Real > > >
operator/(Range d1, 
      Real d2)
{
    typedef _bz_VecExprOp<Range, 
      _bz_VecExprConstant<Real>, 
      _bz_Divide<int, Real> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(d1, 
      _bz_VecExprConstant<Real>(d2)));
}

// TinyVector<P_numtype1, N_length1> / Real
template<class P_numtype1, int N_length1>
inline
_bz_VecExpr<_bz_VecExprOp<TinyVectorIterConst<P_numtype1, N_length1>, 
      _bz_VecExprConstant<Real>,
      _bz_Divide<P_numtype1, Real > > >
operator/(const TinyVector<P_numtype1, N_length1>& d1, 
      Real d2)
{
    typedef _bz_VecExprOp<TinyVectorIterConst<P_numtype1, N_length1>, 
      _bz_VecExprConstant<Real>, 
      _bz_Divide<P_numtype1, Real> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(d1.beginFast(), 
      _bz_VecExprConstant<Real>(d2)));
}

// Real / Vector<P_numtype2>
template<class P_numtype2>
inline
_bz_VecExpr<_bz_VecExprOp<_bz_VecExprConstant<Real>, 
      VectorIterConst<P_numtype2>,
      _bz_Divide<Real, P_numtype2 > > >
operator/(Real d1, 
      const Vector<P_numtype2>& d2)
{
    typedef _bz_VecExprOp<_bz_VecExprConstant<Real>, 
      VectorIterConst<P_numtype2>, 
      _bz_Divide<Real, P_numtype2> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(_bz_VecExprConstant<Real>(d1), 
      d2.beginFast()));
}

// Real / _bz_VecExpr<P_expr2>
template<class P_expr2>
inline
_bz_VecExpr<_bz_VecExprOp<_bz_VecExprConstant<Real>, 
      _bz_VecExpr<P_expr2>,
      _bz_Divide<Real, typename P_expr2::T_numtype > > >
operator/(Real d1, 
      _bz_VecExpr<P_expr2> d2)
{
    typedef _bz_VecExprOp<_bz_VecExprConstant<Real>, 
      _bz_VecExpr<P_expr2>, 
      _bz_Divide<Real, typename P_expr2::T_numtype> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(_bz_VecExprConstant<Real>(d1), 
      d2));
}

// Real / VectorPick<P_numtype2>
template<class P_numtype2>
inline
_bz_VecExpr<_bz_VecExprOp<_bz_VecExprConstant<Real>, 
      VectorPickIterConst<P_numtype2>,
      _bz_Divide<Real, P_numtype2 > > >
operator/(Real d1, 
      const VectorPick<P_numtype2>& d2)
{
    typedef _bz_VecExprOp<_bz_VecExprConstant<Real>, 
      VectorPickIterConst<P_numtype2>, 
      _bz_Divide<Real, P_numtype2> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(_bz_VecExprConstant<Real>(d1), 
      d2.beginFast()));
}

// Real / Range

inline
_bz_VecExpr<_bz_VecExprOp<_bz_VecExprConstant<Real>, 
      Range,
      _bz_Divide<Real, int > > >
operator/(Real d1, 
      Range d2)
{
    typedef _bz_VecExprOp<_bz_VecExprConstant<Real>, 
      Range, 
      _bz_Divide<Real, int> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(_bz_VecExprConstant<Real>(d1), 
      d2));
}

// Real / TinyVector<P_numtype2, N_length2>
template<class P_numtype2, int N_length2>
inline
_bz_VecExpr<_bz_VecExprOp<_bz_VecExprConstant<Real>, 
      TinyVectorIterConst<P_numtype2, N_length2>,
      _bz_Divide<Real, P_numtype2 > > >
operator/(Real d1, 
      const TinyVector<P_numtype2, N_length2>& d2)
{
    typedef _bz_VecExprOp<_bz_VecExprConstant<Real>, 
      TinyVectorIterConst<P_numtype2, N_length2>, 
      _bz_Divide<Real, P_numtype2> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(_bz_VecExprConstant<Real>(d1), 
      d2.beginFast()));
}

// Vector<P_numtype1> > Real
template<class P_numtype1>
inline
_bz_VecExpr<_bz_VecExprOp<VectorIterConst<P_numtype1>, 
      _bz_VecExprConstant<Real>,
      _bz_Greater<P_numtype1, Real > > >
operator>(const Vector<P_numtype1>& d1, 
      Real d2)
{
    typedef _bz_VecExprOp<VectorIterConst<P_numtype1>, 
      _bz_VecExprConstant<Real>, 
      _bz_Greater<P_numtype1, Real> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(d1.beginFast(), 
      _bz_VecExprConstant<Real>(d2)));
}

// _bz_VecExpr<P_expr1> > Real
template<class P_expr1>
inline
_bz_VecExpr<_bz_VecExprOp<_bz_VecExpr<P_expr1>, 
      _bz_VecExprConstant<Real>,
      _bz_Greater<typename P_expr1::T_numtype, Real > > >
operator>(_bz_VecExpr<P_expr1> d1, 
      Real d2)
{
    typedef _bz_VecExprOp<_bz_VecExpr<P_expr1>, 
      _bz_VecExprConstant<Real>, 
      _bz_Greater<typename P_expr1::T_numtype, Real> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(d1, 
      _bz_VecExprConstant<Real>(d2)));
}

// VectorPick<P_numtype1> > Real
template<class P_numtype1>
inline
_bz_VecExpr<_bz_VecExprOp<VectorPickIterConst<P_numtype1>, 
      _bz_VecExprConstant<Real>,
      _bz_Greater<P_numtype1, Real > > >
operator>(const VectorPick<P_numtype1>& d1, 
      Real d2)
{
    typedef _bz_VecExprOp<VectorPickIterConst<P_numtype1>, 
      _bz_VecExprConstant<Real>, 
      _bz_Greater<P_numtype1, Real> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(d1.beginFast(), 
      _bz_VecExprConstant<Real>(d2)));
}

// Range > Real

inline
_bz_VecExpr<_bz_VecExprOp<Range, 
      _bz_VecExprConstant<Real>,
      _bz_Greater<int, Real > > >
operator>(Range d1, 
      Real d2)
{
    typedef _bz_VecExprOp<Range, 
      _bz_VecExprConstant<Real>, 
      _bz_Greater<int, Real> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(d1, 
      _bz_VecExprConstant<Real>(d2)));
}

// TinyVector<P_numtype1, N_length1> > Real
template<class P_numtype1, int N_length1>
inline
_bz_VecExpr<_bz_VecExprOp<TinyVectorIterConst<P_numtype1, N_length1>, 
      _bz_VecExprConstant<Real>,
      _bz_Greater<P_numtype1, Real > > >
operator>(const TinyVector<P_numtype1, N_length1>& d1, 
      Real d2)
{
    typedef _bz_VecExprOp<TinyVectorIterConst<P_numtype1, N_length1>, 
      _bz_VecExprConstant<Real>, 
      _bz_Greater<P_numtype1, Real> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(d1.beginFast(), 
      _bz_VecExprConstant<Real>(d2)));
}

// Real > Vector<P_numtype2>
template<class P_numtype2>
inline
_bz_VecExpr<_bz_VecExprOp<_bz_VecExprConstant<Real>, 
      VectorIterConst<P_numtype2>,
      _bz_Greater<Real, P_numtype2 > > >
operator>(Real d1, 
      const Vector<P_numtype2>& d2)
{
    typedef _bz_VecExprOp<_bz_VecExprConstant<Real>, 
      VectorIterConst<P_numtype2>, 
      _bz_Greater<Real, P_numtype2> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(_bz_VecExprConstant<Real>(d1), 
      d2.beginFast()));
}

// Real > _bz_VecExpr<P_expr2>
template<class P_expr2>
inline
_bz_VecExpr<_bz_VecExprOp<_bz_VecExprConstant<Real>, 
      _bz_VecExpr<P_expr2>,
      _bz_Greater<Real, typename P_expr2::T_numtype > > >
operator>(Real d1, 
      _bz_VecExpr<P_expr2> d2)
{
    typedef _bz_VecExprOp<_bz_VecExprConstant<Real>, 
      _bz_VecExpr<P_expr2>, 
      _bz_Greater<Real, typename P_expr2::T_numtype> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(_bz_VecExprConstant<Real>(d1), 
      d2));
}

// Real > VectorPick<P_numtype2>
template<class P_numtype2>
inline
_bz_VecExpr<_bz_VecExprOp<_bz_VecExprConstant<Real>, 
      VectorPickIterConst<P_numtype2>,
      _bz_Greater<Real, P_numtype2 > > >
operator>(Real d1, 
      const VectorPick<P_numtype2>& d2)
{
    typedef _bz_VecExprOp<_bz_VecExprConstant<Real>, 
      VectorPickIterConst<P_numtype2>, 
      _bz_Greater<Real, P_numtype2> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(_bz_VecExprConstant<Real>(d1), 
      d2.beginFast()));
}

// Real > Range

inline
_bz_VecExpr<_bz_VecExprOp<_bz_VecExprConstant<Real>, 
      Range,
      _bz_Greater<Real, int > > >
operator>(Real d1, 
      Range d2)
{
    typedef _bz_VecExprOp<_bz_VecExprConstant<Real>, 
      Range, 
      _bz_Greater<Real, int> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(_bz_VecExprConstant<Real>(d1), 
      d2));
}

// Real > TinyVector<P_numtype2, N_length2>
template<class P_numtype2, int N_length2>
inline
_bz_VecExpr<_bz_VecExprOp<_bz_VecExprConstant<Real>, 
      TinyVectorIterConst<P_numtype2, N_length2>,
      _bz_Greater<Real, P_numtype2 > > >
operator>(Real d1, 
      const TinyVector<P_numtype2, N_length2>& d2)
{
    typedef _bz_VecExprOp<_bz_VecExprConstant<Real>, 
      TinyVectorIterConst<P_numtype2, N_length2>, 
      _bz_Greater<Real, P_numtype2> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(_bz_VecExprConstant<Real>(d1), 
      d2.beginFast()));
}

// Vector<P_numtype1> < Real
template<class P_numtype1>
inline
_bz_VecExpr<_bz_VecExprOp<VectorIterConst<P_numtype1>, 
      _bz_VecExprConstant<Real>,
      _bz_Less<P_numtype1, Real > > >
operator<(const Vector<P_numtype1>& d1, 
      Real d2)
{
    typedef _bz_VecExprOp<VectorIterConst<P_numtype1>, 
      _bz_VecExprConstant<Real>, 
      _bz_Less<P_numtype1, Real> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(d1.beginFast(), 
      _bz_VecExprConstant<Real>(d2)));
}

// _bz_VecExpr<P_expr1> < Real
template<class P_expr1>
inline
_bz_VecExpr<_bz_VecExprOp<_bz_VecExpr<P_expr1>, 
      _bz_VecExprConstant<Real>,
      _bz_Less<typename P_expr1::T_numtype, Real > > >
operator<(_bz_VecExpr<P_expr1> d1, 
      Real d2)
{
    typedef _bz_VecExprOp<_bz_VecExpr<P_expr1>, 
      _bz_VecExprConstant<Real>, 
      _bz_Less<typename P_expr1::T_numtype, Real> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(d1, 
      _bz_VecExprConstant<Real>(d2)));
}

// VectorPick<P_numtype1> < Real
template<class P_numtype1>
inline
_bz_VecExpr<_bz_VecExprOp<VectorPickIterConst<P_numtype1>, 
      _bz_VecExprConstant<Real>,
      _bz_Less<P_numtype1, Real > > >
operator<(const VectorPick<P_numtype1>& d1, 
      Real d2)
{
    typedef _bz_VecExprOp<VectorPickIterConst<P_numtype1>, 
      _bz_VecExprConstant<Real>, 
      _bz_Less<P_numtype1, Real> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(d1.beginFast(), 
      _bz_VecExprConstant<Real>(d2)));
}

// Range < Real

inline
_bz_VecExpr<_bz_VecExprOp<Range, 
      _bz_VecExprConstant<Real>,
      _bz_Less<int, Real > > >
operator<(Range d1, 
      Real d2)
{
    typedef _bz_VecExprOp<Range, 
      _bz_VecExprConstant<Real>, 
      _bz_Less<int, Real> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(d1, 
      _bz_VecExprConstant<Real>(d2)));
}

// TinyVector<P_numtype1, N_length1> < Real
template<class P_numtype1, int N_length1>
inline
_bz_VecExpr<_bz_VecExprOp<TinyVectorIterConst<P_numtype1, N_length1>, 
      _bz_VecExprConstant<Real>,
      _bz_Less<P_numtype1, Real > > >
operator<(const TinyVector<P_numtype1, N_length1>& d1, 
      Real d2)
{
    typedef _bz_VecExprOp<TinyVectorIterConst<P_numtype1, N_length1>, 
      _bz_VecExprConstant<Real>, 
      _bz_Less<P_numtype1, Real> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(d1.beginFast(), 
      _bz_VecExprConstant<Real>(d2)));
}

// Real < Vector<P_numtype2>
template<class P_numtype2>
inline
_bz_VecExpr<_bz_VecExprOp<_bz_VecExprConstant<Real>, 
      VectorIterConst<P_numtype2>,
      _bz_Less<Real, P_numtype2 > > >
operator<(Real d1, 
      const Vector<P_numtype2>& d2)
{
    typedef _bz_VecExprOp<_bz_VecExprConstant<Real>, 
      VectorIterConst<P_numtype2>, 
      _bz_Less<Real, P_numtype2> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(_bz_VecExprConstant<Real>(d1), 
      d2.beginFast()));
}

// Real < _bz_VecExpr<P_expr2>
template<class P_expr2>
inline
_bz_VecExpr<_bz_VecExprOp<_bz_VecExprConstant<Real>, 
      _bz_VecExpr<P_expr2>,
      _bz_Less<Real, typename P_expr2::T_numtype > > >
operator<(Real d1, 
      _bz_VecExpr<P_expr2> d2)
{
    typedef _bz_VecExprOp<_bz_VecExprConstant<Real>, 
      _bz_VecExpr<P_expr2>, 
      _bz_Less<Real, typename P_expr2::T_numtype> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(_bz_VecExprConstant<Real>(d1), 
      d2));
}

// Real < VectorPick<P_numtype2>
template<class P_numtype2>
inline
_bz_VecExpr<_bz_VecExprOp<_bz_VecExprConstant<Real>, 
      VectorPickIterConst<P_numtype2>,
      _bz_Less<Real, P_numtype2 > > >
operator<(Real d1, 
      const VectorPick<P_numtype2>& d2)
{
    typedef _bz_VecExprOp<_bz_VecExprConstant<Real>, 
      VectorPickIterConst<P_numtype2>, 
      _bz_Less<Real, P_numtype2> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(_bz_VecExprConstant<Real>(d1), 
      d2.beginFast()));
}

// Real < Range

inline
_bz_VecExpr<_bz_VecExprOp<_bz_VecExprConstant<Real>, 
      Range,
      _bz_Less<Real, int > > >
operator<(Real d1, 
      Range d2)
{
    typedef _bz_VecExprOp<_bz_VecExprConstant<Real>, 
      Range, 
      _bz_Less<Real, int> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(_bz_VecExprConstant<Real>(d1), 
      d2));
}

// Real < TinyVector<P_numtype2, N_length2>
template<class P_numtype2, int N_length2>
inline
_bz_VecExpr<_bz_VecExprOp<_bz_VecExprConstant<Real>, 
      TinyVectorIterConst<P_numtype2, N_length2>,
      _bz_Less<Real, P_numtype2 > > >
operator<(Real d1, 
      const TinyVector<P_numtype2, N_length2>& d2)
{
    typedef _bz_VecExprOp<_bz_VecExprConstant<Real>, 
      TinyVectorIterConst<P_numtype2, N_length2>, 
      _bz_Less<Real, P_numtype2> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(_bz_VecExprConstant<Real>(d1), 
      d2.beginFast()));
}

// Vector<P_numtype1> >= Real
template<class P_numtype1>
inline
_bz_VecExpr<_bz_VecExprOp<VectorIterConst<P_numtype1>, 
      _bz_VecExprConstant<Real>,
      _bz_GreaterOrEqual<P_numtype1, Real > > >
operator>=(const Vector<P_numtype1>& d1, 
      Real d2)
{
    typedef _bz_VecExprOp<VectorIterConst<P_numtype1>, 
      _bz_VecExprConstant<Real>, 
      _bz_GreaterOrEqual<P_numtype1, Real> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(d1.beginFast(), 
      _bz_VecExprConstant<Real>(d2)));
}

// _bz_VecExpr<P_expr1> >= Real
template<class P_expr1>
inline
_bz_VecExpr<_bz_VecExprOp<_bz_VecExpr<P_expr1>, 
      _bz_VecExprConstant<Real>,
      _bz_GreaterOrEqual<typename P_expr1::T_numtype, Real > > >
operator>=(_bz_VecExpr<P_expr1> d1, 
      Real d2)
{
    typedef _bz_VecExprOp<_bz_VecExpr<P_expr1>, 
      _bz_VecExprConstant<Real>, 
      _bz_GreaterOrEqual<typename P_expr1::T_numtype, Real> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(d1, 
      _bz_VecExprConstant<Real>(d2)));
}

// VectorPick<P_numtype1> >= Real
template<class P_numtype1>
inline
_bz_VecExpr<_bz_VecExprOp<VectorPickIterConst<P_numtype1>, 
      _bz_VecExprConstant<Real>,
      _bz_GreaterOrEqual<P_numtype1, Real > > >
operator>=(const VectorPick<P_numtype1>& d1, 
      Real d2)
{
    typedef _bz_VecExprOp<VectorPickIterConst<P_numtype1>, 
      _bz_VecExprConstant<Real>, 
      _bz_GreaterOrEqual<P_numtype1, Real> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(d1.beginFast(), 
      _bz_VecExprConstant<Real>(d2)));
}

// Range >= Real

inline
_bz_VecExpr<_bz_VecExprOp<Range, 
      _bz_VecExprConstant<Real>,
      _bz_GreaterOrEqual<int, Real > > >
operator>=(Range d1, 
      Real d2)
{
    typedef _bz_VecExprOp<Range, 
      _bz_VecExprConstant<Real>, 
      _bz_GreaterOrEqual<int, Real> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(d1, 
      _bz_VecExprConstant<Real>(d2)));
}

// TinyVector<P_numtype1, N_length1> >= Real
template<class P_numtype1, int N_length1>
inline
_bz_VecExpr<_bz_VecExprOp<TinyVectorIterConst<P_numtype1, N_length1>, 
      _bz_VecExprConstant<Real>,
      _bz_GreaterOrEqual<P_numtype1, Real > > >
operator>=(const TinyVector<P_numtype1, N_length1>& d1, 
      Real d2)
{
    typedef _bz_VecExprOp<TinyVectorIterConst<P_numtype1, N_length1>, 
      _bz_VecExprConstant<Real>, 
      _bz_GreaterOrEqual<P_numtype1, Real> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(d1.beginFast(), 
      _bz_VecExprConstant<Real>(d2)));
}

// Real >= Vector<P_numtype2>
template<class P_numtype2>
inline
_bz_VecExpr<_bz_VecExprOp<_bz_VecExprConstant<Real>, 
      VectorIterConst<P_numtype2>,
      _bz_GreaterOrEqual<Real, P_numtype2 > > >
operator>=(Real d1, 
      const Vector<P_numtype2>& d2)
{
    typedef _bz_VecExprOp<_bz_VecExprConstant<Real>, 
      VectorIterConst<P_numtype2>, 
      _bz_GreaterOrEqual<Real, P_numtype2> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(_bz_VecExprConstant<Real>(d1), 
      d2.beginFast()));
}

// Real >= _bz_VecExpr<P_expr2>
template<class P_expr2>
inline
_bz_VecExpr<_bz_VecExprOp<_bz_VecExprConstant<Real>, 
      _bz_VecExpr<P_expr2>,
      _bz_GreaterOrEqual<Real, typename P_expr2::T_numtype > > >
operator>=(Real d1, 
      _bz_VecExpr<P_expr2> d2)
{
    typedef _bz_VecExprOp<_bz_VecExprConstant<Real>, 
      _bz_VecExpr<P_expr2>, 
      _bz_GreaterOrEqual<Real, typename P_expr2::T_numtype> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(_bz_VecExprConstant<Real>(d1), 
      d2));
}

// Real >= VectorPick<P_numtype2>
template<class P_numtype2>
inline
_bz_VecExpr<_bz_VecExprOp<_bz_VecExprConstant<Real>, 
      VectorPickIterConst<P_numtype2>,
      _bz_GreaterOrEqual<Real, P_numtype2 > > >
operator>=(Real d1, 
      const VectorPick<P_numtype2>& d2)
{
    typedef _bz_VecExprOp<_bz_VecExprConstant<Real>, 
      VectorPickIterConst<P_numtype2>, 
      _bz_GreaterOrEqual<Real, P_numtype2> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(_bz_VecExprConstant<Real>(d1), 
      d2.beginFast()));
}

// Real >= Range

inline
_bz_VecExpr<_bz_VecExprOp<_bz_VecExprConstant<Real>, 
      Range,
      _bz_GreaterOrEqual<Real, int > > >
operator>=(Real d1, 
      Range d2)
{
    typedef _bz_VecExprOp<_bz_VecExprConstant<Real>, 
      Range, 
      _bz_GreaterOrEqual<Real, int> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(_bz_VecExprConstant<Real>(d1), 
      d2));
}

// Real >= TinyVector<P_numtype2, N_length2>
template<class P_numtype2, int N_length2>
inline
_bz_VecExpr<_bz_VecExprOp<_bz_VecExprConstant<Real>, 
      TinyVectorIterConst<P_numtype2, N_length2>,
      _bz_GreaterOrEqual<Real, P_numtype2 > > >
operator>=(Real d1, 
      const TinyVector<P_numtype2, N_length2>& d2)
{
    typedef _bz_VecExprOp<_bz_VecExprConstant<Real>, 
      TinyVectorIterConst<P_numtype2, N_length2>, 
      _bz_GreaterOrEqual<Real, P_numtype2> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(_bz_VecExprConstant<Real>(d1), 
      d2.beginFast()));
}

// Vector<P_numtype1> <= Real
template<class P_numtype1>
inline
_bz_VecExpr<_bz_VecExprOp<VectorIterConst<P_numtype1>, 
      _bz_VecExprConstant<Real>,
      _bz_LessOrEqual<P_numtype1, Real > > >
operator<=(const Vector<P_numtype1>& d1, 
      Real d2)
{
    typedef _bz_VecExprOp<VectorIterConst<P_numtype1>, 
      _bz_VecExprConstant<Real>, 
      _bz_LessOrEqual<P_numtype1, Real> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(d1.beginFast(), 
      _bz_VecExprConstant<Real>(d2)));
}

// _bz_VecExpr<P_expr1> <= Real
template<class P_expr1>
inline
_bz_VecExpr<_bz_VecExprOp<_bz_VecExpr<P_expr1>, 
      _bz_VecExprConstant<Real>,
      _bz_LessOrEqual<typename P_expr1::T_numtype, Real > > >
operator<=(_bz_VecExpr<P_expr1> d1, 
      Real d2)
{
    typedef _bz_VecExprOp<_bz_VecExpr<P_expr1>, 
      _bz_VecExprConstant<Real>, 
      _bz_LessOrEqual<typename P_expr1::T_numtype, Real> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(d1, 
      _bz_VecExprConstant<Real>(d2)));
}

// VectorPick<P_numtype1> <= Real
template<class P_numtype1>
inline
_bz_VecExpr<_bz_VecExprOp<VectorPickIterConst<P_numtype1>, 
      _bz_VecExprConstant<Real>,
      _bz_LessOrEqual<P_numtype1, Real > > >
operator<=(const VectorPick<P_numtype1>& d1, 
      Real d2)
{
    typedef _bz_VecExprOp<VectorPickIterConst<P_numtype1>, 
      _bz_VecExprConstant<Real>, 
      _bz_LessOrEqual<P_numtype1, Real> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(d1.beginFast(), 
      _bz_VecExprConstant<Real>(d2)));
}

// Range <= Real

inline
_bz_VecExpr<_bz_VecExprOp<Range, 
      _bz_VecExprConstant<Real>,
      _bz_LessOrEqual<int, Real > > >
operator<=(Range d1, 
      Real d2)
{
    typedef _bz_VecExprOp<Range, 
      _bz_VecExprConstant<Real>, 
      _bz_LessOrEqual<int, Real> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(d1, 
      _bz_VecExprConstant<Real>(d2)));
}

// TinyVector<P_numtype1, N_length1> <= Real
template<class P_numtype1, int N_length1>
inline
_bz_VecExpr<_bz_VecExprOp<TinyVectorIterConst<P_numtype1, N_length1>, 
      _bz_VecExprConstant<Real>,
      _bz_LessOrEqual<P_numtype1, Real > > >
operator<=(const TinyVector<P_numtype1, N_length1>& d1, 
      Real d2)
{
    typedef _bz_VecExprOp<TinyVectorIterConst<P_numtype1, N_length1>, 
      _bz_VecExprConstant<Real>, 
      _bz_LessOrEqual<P_numtype1, Real> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(d1.beginFast(), 
      _bz_VecExprConstant<Real>(d2)));
}

// Real <= Vector<P_numtype2>
template<class P_numtype2>
inline
_bz_VecExpr<_bz_VecExprOp<_bz_VecExprConstant<Real>, 
      VectorIterConst<P_numtype2>,
      _bz_LessOrEqual<Real, P_numtype2 > > >
operator<=(Real d1, 
      const Vector<P_numtype2>& d2)
{
    typedef _bz_VecExprOp<_bz_VecExprConstant<Real>, 
      VectorIterConst<P_numtype2>, 
      _bz_LessOrEqual<Real, P_numtype2> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(_bz_VecExprConstant<Real>(d1), 
      d2.beginFast()));
}

// Real <= _bz_VecExpr<P_expr2>
template<class P_expr2>
inline
_bz_VecExpr<_bz_VecExprOp<_bz_VecExprConstant<Real>, 
      _bz_VecExpr<P_expr2>,
      _bz_LessOrEqual<Real, typename P_expr2::T_numtype > > >
operator<=(Real d1, 
      _bz_VecExpr<P_expr2> d2)
{
    typedef _bz_VecExprOp<_bz_VecExprConstant<Real>, 
      _bz_VecExpr<P_expr2>, 
      _bz_LessOrEqual<Real, typename P_expr2::T_numtype> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(_bz_VecExprConstant<Real>(d1), 
      d2));
}

// Real <= VectorPick<P_numtype2>
template<class P_numtype2>
inline
_bz_VecExpr<_bz_VecExprOp<_bz_VecExprConstant<Real>, 
      VectorPickIterConst<P_numtype2>,
      _bz_LessOrEqual<Real, P_numtype2 > > >
operator<=(Real d1, 
      const VectorPick<P_numtype2>& d2)
{
    typedef _bz_VecExprOp<_bz_VecExprConstant<Real>, 
      VectorPickIterConst<P_numtype2>, 
      _bz_LessOrEqual<Real, P_numtype2> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(_bz_VecExprConstant<Real>(d1), 
      d2.beginFast()));
}

// Real <= Range

inline
_bz_VecExpr<_bz_VecExprOp<_bz_VecExprConstant<Real>, 
      Range,
      _bz_LessOrEqual<Real, int > > >
operator<=(Real d1, 
      Range d2)
{
    typedef _bz_VecExprOp<_bz_VecExprConstant<Real>, 
      Range, 
      _bz_LessOrEqual<Real, int> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(_bz_VecExprConstant<Real>(d1), 
      d2));
}

// Real <= TinyVector<P_numtype2, N_length2>
template<class P_numtype2, int N_length2>
inline
_bz_VecExpr<_bz_VecExprOp<_bz_VecExprConstant<Real>, 
      TinyVectorIterConst<P_numtype2, N_length2>,
      _bz_LessOrEqual<Real, P_numtype2 > > >
operator<=(Real d1, 
      const TinyVector<P_numtype2, N_length2>& d2)
{
    typedef _bz_VecExprOp<_bz_VecExprConstant<Real>, 
      TinyVectorIterConst<P_numtype2, N_length2>, 
      _bz_LessOrEqual<Real, P_numtype2> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(_bz_VecExprConstant<Real>(d1), 
      d2.beginFast()));
}

// Vector<P_numtype1> == Real
template<class P_numtype1>
inline
_bz_VecExpr<_bz_VecExprOp<VectorIterConst<P_numtype1>, 
      _bz_VecExprConstant<Real>,
      _bz_Equal<P_numtype1, Real > > >
operator==(const Vector<P_numtype1>& d1, 
      Real d2)
{
    typedef _bz_VecExprOp<VectorIterConst<P_numtype1>, 
      _bz_VecExprConstant<Real>, 
      _bz_Equal<P_numtype1, Real> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(d1.beginFast(), 
      _bz_VecExprConstant<Real>(d2)));
}

// _bz_VecExpr<P_expr1> == Real
template<class P_expr1>
inline
_bz_VecExpr<_bz_VecExprOp<_bz_VecExpr<P_expr1>, 
      _bz_VecExprConstant<Real>,
      _bz_Equal<typename P_expr1::T_numtype, Real > > >
operator==(_bz_VecExpr<P_expr1> d1, 
      Real d2)
{
    typedef _bz_VecExprOp<_bz_VecExpr<P_expr1>, 
      _bz_VecExprConstant<Real>, 
      _bz_Equal<typename P_expr1::T_numtype, Real> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(d1, 
      _bz_VecExprConstant<Real>(d2)));
}

// VectorPick<P_numtype1> == Real
template<class P_numtype1>
inline
_bz_VecExpr<_bz_VecExprOp<VectorPickIterConst<P_numtype1>, 
      _bz_VecExprConstant<Real>,
      _bz_Equal<P_numtype1, Real > > >
operator==(const VectorPick<P_numtype1>& d1, 
      Real d2)
{
    typedef _bz_VecExprOp<VectorPickIterConst<P_numtype1>, 
      _bz_VecExprConstant<Real>, 
      _bz_Equal<P_numtype1, Real> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(d1.beginFast(), 
      _bz_VecExprConstant<Real>(d2)));
}

// Range == Real

inline
_bz_VecExpr<_bz_VecExprOp<Range, 
      _bz_VecExprConstant<Real>,
      _bz_Equal<int, Real > > >
operator==(Range d1, 
      Real d2)
{
    typedef _bz_VecExprOp<Range, 
      _bz_VecExprConstant<Real>, 
      _bz_Equal<int, Real> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(d1, 
      _bz_VecExprConstant<Real>(d2)));
}

// TinyVector<P_numtype1, N_length1> == Real
template<class P_numtype1, int N_length1>
inline
_bz_VecExpr<_bz_VecExprOp<TinyVectorIterConst<P_numtype1, N_length1>, 
      _bz_VecExprConstant<Real>,
      _bz_Equal<P_numtype1, Real > > >
operator==(const TinyVector<P_numtype1, N_length1>& d1, 
      Real d2)
{
    typedef _bz_VecExprOp<TinyVectorIterConst<P_numtype1, N_length1>, 
      _bz_VecExprConstant<Real>, 
      _bz_Equal<P_numtype1, Real> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(d1.beginFast(), 
      _bz_VecExprConstant<Real>(d2)));
}

// Real == Vector<P_numtype2>
template<class P_numtype2>
inline
_bz_VecExpr<_bz_VecExprOp<_bz_VecExprConstant<Real>, 
      VectorIterConst<P_numtype2>,
      _bz_Equal<Real, P_numtype2 > > >
operator==(Real d1, 
      const Vector<P_numtype2>& d2)
{
    typedef _bz_VecExprOp<_bz_VecExprConstant<Real>, 
      VectorIterConst<P_numtype2>, 
      _bz_Equal<Real, P_numtype2> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(_bz_VecExprConstant<Real>(d1), 
      d2.beginFast()));
}

// Real == _bz_VecExpr<P_expr2>
template<class P_expr2>
inline
_bz_VecExpr<_bz_VecExprOp<_bz_VecExprConstant<Real>, 
      _bz_VecExpr<P_expr2>,
      _bz_Equal<Real, typename P_expr2::T_numtype > > >
operator==(Real d1, 
      _bz_VecExpr<P_expr2> d2)
{
    typedef _bz_VecExprOp<_bz_VecExprConstant<Real>, 
      _bz_VecExpr<P_expr2>, 
      _bz_Equal<Real, typename P_expr2::T_numtype> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(_bz_VecExprConstant<Real>(d1), 
      d2));
}

// Real == VectorPick<P_numtype2>
template<class P_numtype2>
inline
_bz_VecExpr<_bz_VecExprOp<_bz_VecExprConstant<Real>, 
      VectorPickIterConst<P_numtype2>,
      _bz_Equal<Real, P_numtype2 > > >
operator==(Real d1, 
      const VectorPick<P_numtype2>& d2)
{
    typedef _bz_VecExprOp<_bz_VecExprConstant<Real>, 
      VectorPickIterConst<P_numtype2>, 
      _bz_Equal<Real, P_numtype2> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(_bz_VecExprConstant<Real>(d1), 
      d2.beginFast()));
}

// Real == Range

inline
_bz_VecExpr<_bz_VecExprOp<_bz_VecExprConstant<Real>, 
      Range,
      _bz_Equal<Real, int > > >
operator==(Real d1, 
      Range d2)
{
    typedef _bz_VecExprOp<_bz_VecExprConstant<Real>, 
      Range, 
      _bz_Equal<Real, int> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(_bz_VecExprConstant<Real>(d1), 
      d2));
}

// Real == TinyVector<P_numtype2, N_length2>
template<class P_numtype2, int N_length2>
inline
_bz_VecExpr<_bz_VecExprOp<_bz_VecExprConstant<Real>, 
      TinyVectorIterConst<P_numtype2, N_length2>,
      _bz_Equal<Real, P_numtype2 > > >
operator==(Real d1, 
      const TinyVector<P_numtype2, N_length2>& d2)
{
    typedef _bz_VecExprOp<_bz_VecExprConstant<Real>, 
      TinyVectorIterConst<P_numtype2, N_length2>, 
      _bz_Equal<Real, P_numtype2> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(_bz_VecExprConstant<Real>(d1), 
      d2.beginFast()));
}

// Vector<P_numtype1> != Real
template<class P_numtype1>
inline
_bz_VecExpr<_bz_VecExprOp<VectorIterConst<P_numtype1>, 
      _bz_VecExprConstant<Real>,
      _bz_NotEqual<P_numtype1, Real > > >
operator!=(const Vector<P_numtype1>& d1, 
      Real d2)
{
    typedef _bz_VecExprOp<VectorIterConst<P_numtype1>, 
      _bz_VecExprConstant<Real>, 
      _bz_NotEqual<P_numtype1, Real> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(d1.beginFast(), 
      _bz_VecExprConstant<Real>(d2)));
}

// _bz_VecExpr<P_expr1> != Real
template<class P_expr1>
inline
_bz_VecExpr<_bz_VecExprOp<_bz_VecExpr<P_expr1>, 
      _bz_VecExprConstant<Real>,
      _bz_NotEqual<typename P_expr1::T_numtype, Real > > >
operator!=(_bz_VecExpr<P_expr1> d1, 
      Real d2)
{
    typedef _bz_VecExprOp<_bz_VecExpr<P_expr1>, 
      _bz_VecExprConstant<Real>, 
      _bz_NotEqual<typename P_expr1::T_numtype, Real> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(d1, 
      _bz_VecExprConstant<Real>(d2)));
}

// VectorPick<P_numtype1> != Real
template<class P_numtype1>
inline
_bz_VecExpr<_bz_VecExprOp<VectorPickIterConst<P_numtype1>, 
      _bz_VecExprConstant<Real>,
      _bz_NotEqual<P_numtype1, Real > > >
operator!=(const VectorPick<P_numtype1>& d1, 
      Real d2)
{
    typedef _bz_VecExprOp<VectorPickIterConst<P_numtype1>, 
      _bz_VecExprConstant<Real>, 
      _bz_NotEqual<P_numtype1, Real> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(d1.beginFast(), 
      _bz_VecExprConstant<Real>(d2)));
}

// Range != Real

inline
_bz_VecExpr<_bz_VecExprOp<Range, 
      _bz_VecExprConstant<Real>,
      _bz_NotEqual<int, Real > > >
operator!=(Range d1, 
      Real d2)
{
    typedef _bz_VecExprOp<Range, 
      _bz_VecExprConstant<Real>, 
      _bz_NotEqual<int, Real> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(d1, 
      _bz_VecExprConstant<Real>(d2)));
}

// TinyVector<P_numtype1, N_length1> != Real
template<class P_numtype1, int N_length1>
inline
_bz_VecExpr<_bz_VecExprOp<TinyVectorIterConst<P_numtype1, N_length1>, 
      _bz_VecExprConstant<Real>,
      _bz_NotEqual<P_numtype1, Real > > >
operator!=(const TinyVector<P_numtype1, N_length1>& d1, 
      Real d2)
{
    typedef _bz_VecExprOp<TinyVectorIterConst<P_numtype1, N_length1>, 
      _bz_VecExprConstant<Real>, 
      _bz_NotEqual<P_numtype1, Real> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(d1.beginFast(), 
      _bz_VecExprConstant<Real>(d2)));
}

// Real != Vector<P_numtype2>
template<class P_numtype2>
inline
_bz_VecExpr<_bz_VecExprOp<_bz_VecExprConstant<Real>, 
      VectorIterConst<P_numtype2>,
      _bz_NotEqual<Real, P_numtype2 > > >
operator!=(Real d1, 
      const Vector<P_numtype2>& d2)
{
    typedef _bz_VecExprOp<_bz_VecExprConstant<Real>, 
      VectorIterConst<P_numtype2>, 
      _bz_NotEqual<Real, P_numtype2> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(_bz_VecExprConstant<Real>(d1), 
      d2.beginFast()));
}

// Real != _bz_VecExpr<P_expr2>
template<class P_expr2>
inline
_bz_VecExpr<_bz_VecExprOp<_bz_VecExprConstant<Real>, 
      _bz_VecExpr<P_expr2>,
      _bz_NotEqual<Real, typename P_expr2::T_numtype > > >
operator!=(Real d1, 
      _bz_VecExpr<P_expr2> d2)
{
    typedef _bz_VecExprOp<_bz_VecExprConstant<Real>, 
      _bz_VecExpr<P_expr2>, 
      _bz_NotEqual<Real, typename P_expr2::T_numtype> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(_bz_VecExprConstant<Real>(d1), 
      d2));
}

// Real != VectorPick<P_numtype2>
template<class P_numtype2>
inline
_bz_VecExpr<_bz_VecExprOp<_bz_VecExprConstant<Real>, 
      VectorPickIterConst<P_numtype2>,
      _bz_NotEqual<Real, P_numtype2 > > >
operator!=(Real d1, 
      const VectorPick<P_numtype2>& d2)
{
    typedef _bz_VecExprOp<_bz_VecExprConstant<Real>, 
      VectorPickIterConst<P_numtype2>, 
      _bz_NotEqual<Real, P_numtype2> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(_bz_VecExprConstant<Real>(d1), 
      d2.beginFast()));
}

// Real != Range

inline
_bz_VecExpr<_bz_VecExprOp<_bz_VecExprConstant<Real>, 
      Range,
      _bz_NotEqual<Real, int > > >
operator!=(Real d1, 
      Range d2)
{
    typedef _bz_VecExprOp<_bz_VecExprConstant<Real>, 
      Range, 
      _bz_NotEqual<Real, int> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(_bz_VecExprConstant<Real>(d1), 
      d2));
}

// Real != TinyVector<P_numtype2, N_length2>
template<class P_numtype2, int N_length2>
inline
_bz_VecExpr<_bz_VecExprOp<_bz_VecExprConstant<Real>, 
      TinyVectorIterConst<P_numtype2, N_length2>,
      _bz_NotEqual<Real, P_numtype2 > > >
operator!=(Real d1, 
      const TinyVector<P_numtype2, N_length2>& d2)
{
    typedef _bz_VecExprOp<_bz_VecExprConstant<Real>, 
      TinyVectorIterConst<P_numtype2, N_length2>, 
      _bz_NotEqual<Real, P_numtype2> > T_expr;

    return _bz_VecExpr<T_expr>(T_expr(_bz_VecExprConstant<Real>(d1), 
      d2.beginFast()));
}

// Expression templates for blitz::Array
BZ_DECLARE_ARRAY_ET_SCALAR_OPS(Real)

BZ_NAMESPACE_END
#endif

#endif
