#ifndef FIELD_OPS_H
#define FIELD_OPS_H

#include "AMReX_MultiFab.H"

namespace amr_wind {
namespace field_ops {

template <typename T1, typename T2>
inline void
add(T1& dst,
    const T2& src,
    int srccomp,
    int dstcomp,
    int numcomp,
    const amrex::IntVect& nghost)
{
    const int nlevels = dst.repo().num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        amrex::MultiFab::Add(
            dst(lev), src(lev), srccomp, dstcomp, numcomp, nghost);
    }
}

template <typename T1, typename T2>
inline void
add(T1& dst, const T2& src, int srccomp, int dstcomp, int numcomp, int nghost)
{
    add(dst, src, srccomp, dstcomp, numcomp, amrex::IntVect(nghost));
}

template <typename T1, typename T2>
inline void copy(
    T1& dst,
    const T2& src,
    int srccomp,
    int dstcomp,
    int numcomp,
    const amrex::IntVect& nghost)
{
    const int nlevels = dst.repo().num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        amrex::MultiFab::Copy(
            dst(lev), src(lev), srccomp, dstcomp, numcomp, nghost);
    }
}

template <typename T1, typename T2>
inline void
copy(T1& dst, const T2& src, int srccomp, int dstcomp, int numcomp, int nghost)
{
    copy(dst, src, srccomp, dstcomp, numcomp, amrex::IntVect(nghost));
}

template <typename T1, typename T2>
inline void saxpy(
    T1& dst,
    amrex::Real a,
    const T2& src,
    int srccomp,
    int dstcomp,
    int numcomp,
    const amrex::IntVect& nghost)
{
    const int nlevels = dst.repo().num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        amrex::MultiFab::Saxpy(
            dst(lev), a, src(lev), srccomp, dstcomp, numcomp, nghost);
    }
}

template <typename T1, typename T2>
inline void saxpy(
    T1& dst,
    amrex::Real a,
    const T2& src,
    int srccomp,
    int dstcomp,
    int numcomp,
    int nghost)
{
    saxpy(dst, a, src, srccomp, dstcomp, numcomp, amrex::IntVect(nghost));
}

template <typename T1, typename T2>
inline void xpay(
    T1& dst,
    amrex::Real a,
    const T2& src,
    int srccomp,
    int dstcomp,
    int numcomp,
    const amrex::IntVect& nghost)
{
    const int nlevels = dst.repo().num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        amrex::MultiFab::Xpay(
            dst(lev), a, src(lev), srccomp, dstcomp, numcomp, nghost);
    }
}

template <typename T1, typename T2>
inline void xpay(
    T1& dst,
    amrex::Real a,
    const T2& src,
    int srccomp,
    int dstcomp,
    int numcomp,
    int nghost)
{
    xpay(dst, a, src, srccomp, dstcomp, numcomp, amrex::IntVect(nghost));
}

template <typename T1, typename T2, typename T3>
inline void lincomb(
    T1& dst,
    amrex::Real a,
    const T2& x,
    int xcomp,
    amrex::Real b,
    const T3& y,
    int ycomp,
    int dstcomp,
    int numcomp,
    const amrex::IntVect& nghost)
{
    const int nlevels = dst.repo().num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        amrex::MultiFab::LinComb(
            dst(lev), a, x(lev), xcomp, b, y(lev), ycomp, dstcomp, numcomp,
            nghost);
    }
}

template <typename T1, typename T2, typename T3>
inline void lincomb(
    T1& dst,
    amrex::Real a,
    const T2& x,
    int xcomp,
    amrex::Real b,
    const T3& y,
    int ycomp,
    int dstcomp,
    int numcomp,
    int nghost)
{
    lincomb(dst, a, x, xcomp, b, y, ycomp, dstcomp, numcomp, amrex::IntVect(nghost));
}

} // namespace field_ops
} // namespace amr_wind

#endif /* FIELD_OPS_H */
