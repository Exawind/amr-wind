#ifndef FILTER_H
#define FILTER_H

#include "amr-wind/fvm/fvm_utils.H"
#include "AMReX_Array4.H"
#include "AMReX_Geometry.H"
#include "AMReX_MFIter.H"

namespace amr_wind {
namespace fvm {

/** filter operator
 *  \ingroup fvm
 */
template <typename FTypeIn, typename FTypeOut>
struct Filter
{
    /**
     *  \param filterphi The filter field \f$\widetilde{ \mathbf{\phi}}\f$
     *  \param phi The input file \f$\phi\f$
     */
    Filter(FTypeOut& filterphi, const FTypeIn& phi)
        : m_filterphi(filterphi), m_phi(phi)
    {
        AMREX_ALWAYS_ASSERT(
            filterphi.num_comp() == AMREX_SPACEDIM * phi.num_comp());
    }

    template <typename Stencil>
    void apply(const int lev, const amrex::MFIter& mfi) const
    {
        const int ncomp = m_phi.num_comp();
        const auto& geom = m_phi.repo().mesh().Geom(lev);
        const auto& filterphi_arr = m_filterphi(lev).array(mfi);
        const auto& phi_arr = m_phi(lev).const_array(mfi);

        const auto& bx_in = mfi.tilebox();
        const auto& bx = Stencil::box(bx_in, geom);
        if (bx.isEmpty()) return;

        amrex::ParallelFor(
            bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                for (int icomp = 0; icomp < ncomp; icomp++) {
                    amrex::Real fp1 = Stencil::f00;
                    amrex::Real f = Stencil::f01;
                    amrex::Real fm1 = Stencil::f02;
                    filterphi_arr(i, j, k, icomp * AMREX_SPACEDIM + 0) =
                        (fp1 * phi_arr(i + 1, j, k, icomp) +
                         f * phi_arr(i, j, k, icomp) +
                         fm1 * phi_arr(i - 1, j, k, icomp));

                    fp1 = Stencil::f10;
                    f = Stencil::f11;
                    fm1 = Stencil::f12;
                    filterphi_arr(i, j, k, icomp * AMREX_SPACEDIM + 1) =
                        (fp1 * phi_arr(i, j + 1, k, icomp) +
                         f * phi_arr(i, j, k, icomp) +
                         fm1 * phi_arr(i, j - 1, k, icomp));

                    fp1 = Stencil::f20;
                    f = Stencil::f21;
                    fm1 = Stencil::f22;
                    filterphi_arr(i, j, k, icomp * AMREX_SPACEDIM + 2) =
                        (fp1 * phi_arr(i, j, k + 1, icomp) +
                         f * phi_arr(i, j, k, icomp) +
                         fm1 * phi_arr(i, j, k - 1, icomp));
                }
            });
    }

    FTypeOut& m_filterphi;
    const FTypeIn& m_phi;
};

/** Compute the filtered value of a given field
 *  \ingroup fvm
 *
 *  \param filterphi [inout] Field where the filtered term is populated
 *  \param phi [in] Field whose filter is computed
 */
template <typename FTypeIn, typename FTypeOut>
inline void filter(FTypeOut& filterphi, const FTypeIn& phi)
{
    BL_PROFILE("amr-wind::fvm::filter");
    Filter<FTypeIn, FTypeOut> filtering(filterphi, phi);
    impl::apply(filtering, phi);
}

/** Compute the filter of a given field and return it as ScratchField
 *  \ingroup fvm
 *
 *  \param phi [in] Field whose filter is computed
 */
template <typename FType>
inline std::unique_ptr<ScratchField> filter(const FType& phi)
{
    const std::string gname = phi.name() + "_filter";
    auto filterphi =
        phi.repo().create_scratch_field(gname, phi.num_comp() * AMREX_SPACEDIM);
    filter(*filterphi, phi);
    return filterphi;
}

} // namespace fvm
} // namespace amr_wind

#endif /* FILTER_H */
