#ifndef FIELD_OPS_H
#define FIELD_OPS_H

#include "amr-wind/core/Field.H"
#include "AMReX_MultiFab.H"

/**
 *  \defgroup field_ops Field operations
 *
 *  Operators for manipulating fields and computing derived quantities.
 *
 *  This group documents functions that provide various operations defined on
 *  the field. These operations can be grouped into two major categories: simple
 *  linear algebra operations (e.g., saxpy), and discretization operators such
 *  as gradient, laplacian computations.
 *
 *  \ingroup fields
 */

namespace amr_wind {
namespace field_ops {

/** Add two fields \f$y = x + y\f$
 *  \ingroup field_ops
 *
 *
 *  \tparam T1 Field or ScratchField
 *  \tparam T2 Field or ScratchField
 *  \param [out] dst Field that is updated
 *  \param [in] src Field to be added
 *  \param [in] srccomp Starting component index of source field
 *  \param [in] dstcomp Starting component index of destination field
 *  \param [in] numcomp Number of components to be updated
 *  \param [in] nghost Number of ghost cells to be updated
 */
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

/** Add two fields \f$y = x + y\f$
 *  \ingroup field_ops
 *
 *  \tparam T1 Field or ScratchField
 *  \tparam T2 Field or ScratchField
 *  \param [out] dst Field that is updated
 *  \param [in] src Field to be added
 *  \param [in] srccomp Starting component index of source field
 *  \param [in] dstcomp Starting component index of destination field
 *  \param [in] numcomp Number of components to be updated
 *  \param [in] nghost Number of ghost cells to be updated
 */
template <typename T1, typename T2>
inline void
add(T1& dst, const T2& src, int srccomp, int dstcomp, int numcomp, int nghost)
{
    add(dst, src, srccomp, dstcomp, numcomp, amrex::IntVect(nghost));
}

/** Copy source field to destination field
 *  \ingroup field_ops
 *
 *  \tparam T1 Field or ScratchField
 *  \tparam T2 Field or ScratchField
 *  \param [out] dst Field that is updated
 *  \param [in] src Field to be added
 *  \param [in] srccomp Starting component index of source field
 *  \param [in] dstcomp Starting component index of destination field
 *  \param [in] numcomp Number of components to be updated
 *  \param [in] nghost Number of ghost cells to be updated
 */
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

/** Copy source field to destination field
 *  \ingroup field_ops
 *
 *  \tparam T1 Field or ScratchField
 *  \tparam T2 Field or ScratchField
 *  \param [out] dst Field that is updated
 *  \param [in] src Field to be added
 *  \param [in] srccomp Starting component index of source field
 *  \param [in] dstcomp Starting component index of destination field
 *  \param [in] numcomp Number of components to be updated
 *  \param [in] nghost Number of ghost cells to be updated
 */
template <typename T1, typename T2>
inline void
copy(T1& dst, const T2& src, int srccomp, int dstcomp, int numcomp, int nghost)
{
    copy(dst, src, srccomp, dstcomp, numcomp, amrex::IntVect(nghost));
}

/** Perform operation \f$y = a x + y\f$
 *  \ingroup field_ops
 *
 *  \tparam T1 Field or ScratchField
 *  \tparam T2 Field or ScratchField
 *  \param [out] dst Field that is updated
 *  \param [in] a  Coefficient
 *  \param [in] src Field to be added
 *  \param [in] srccomp Starting component index of source field
 *  \param [in] dstcomp Starting component index of destination field
 *  \param [in] numcomp Number of components to be updated
 *  \param [in] nghost Number of ghost cells to be updated
 */
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

/** Perform operation \f$y = a x + y\f$
 *  \ingroup field_ops
 *
 *  \tparam T1 Field or ScratchField
 *  \tparam T2 Field or ScratchField
 *  \param [out] dst Field that is updated
 *  \param [in] a  Coefficient
 *  \param [in] src Field to be added
 *  \param [in] srccomp Starting component index of source field
 *  \param [in] dstcomp Starting component index of destination field
 *  \param [in] numcomp Number of components to be updated
 *  \param [in] nghost Number of ghost cells to be updated
 */
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

/** Perform operation \f$y = x + a y\f$
 *  \ingroup field_ops
 *
 *  \tparam T1 Field or ScratchField
 *  \tparam T2 Field or ScratchField
 *  \param [out] dst Field that is updated
 *  \param [in] a  Coefficient
 *  \param [in] src Field to be added
 *  \param [in] srccomp Starting component index of source field
 *  \param [in] dstcomp Starting component index of destination field
 *  \param [in] numcomp Number of components to be updated
 *  \param [in] nghost Number of ghost cells to be updated
 */
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

/** Perform operation \f$y = x + a y\f$
 *  \ingroup field_ops
 *
 *  \tparam T1 Field or ScratchField
 *  \tparam T2 Field or ScratchField
 *  \param [out] dst Field that is updated
 *  \param [in] a  Coefficient
 *  \param [in] src Field to be added
 *  \param [in] srccomp Starting component index of source field
 *  \param [in] dstcomp Starting component index of destination field
 *  \param [in] numcomp Number of components to be updated
 *  \param [in] nghost Number of ghost cells to be updated
 */
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

/** Perform operation \f$z = a x + b y\f$
 *  \ingroup field_ops
 *
 *  \tparam T1 Field or ScratchField
 *  \tparam T2 Field or ScratchField
 *  \param [out] dst Field that is updated
 *  \param [in] a  Coefficient
 *  \param [in] x Field to be added
 *  \param [in] xcomp Starting component index of `x` field
 *  \param [in] b  Coefficient
 *  \param [in] y Field to be added
 *  \param [in] ycomp Starting component index of `y` field
 *  \param [in] dstcomp Starting component index of destination field
 *  \param [in] numcomp Number of components to be updated
 *  \param [in] nghost Number of ghost cells to be updated
 */
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

/** Perform operation \f$z = a x + b y\f$
 *  \ingroup field_ops
 *
 *  \tparam T1 Field or ScratchField
 *  \tparam T2 Field or ScratchField
 *  \param [out] dst Field that is updated
 *  \param [in] a  Coefficient
 *  \param [in] x Field to be added
 *  \param [in] xcomp Starting component index of `x` field
 *  \param [in] b  Coefficient
 *  \param [in] y Field to be added
 *  \param [in] ycomp Starting component index of `y` field
 *  \param [in] dstcomp Starting component index of destination field
 *  \param [in] numcomp Number of components to be updated
 *  \param [in] nghost Number of ghost cells to be updated
 */
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

/** Set the lower bound for a given scalar field   
 * \ingroup field_ops
 * 
 * \param [in] min_value Minimum value for a given field
 * \param [in] icomp Component to set the minimum bound
 */
template<typename FType>
inline void lower_bound(
    FType& field,
    const amrex::Real min_value,
    const int icomp = 0)
{
  const auto& repo = field.repo();
  const int nlevels = repo.num_active_levels();
  for (int lev=0; lev < nlevels; ++lev) {
    
    for (amrex::MFIter mfi(field(lev)); mfi.isValid(); ++mfi) {
      const auto& bx = mfi.tilebox();
      const auto& field_arr = field(lev).array(mfi);

      amrex::ParallelFor(
        bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          field_arr(i,j,k,icomp) = amrex::max(min_value, field_arr(i,j,k,icomp));
      });
    }
  }

}

} // namespace field_ops
} // namespace amr_wind

#endif /* FIELD_OPS_H */
