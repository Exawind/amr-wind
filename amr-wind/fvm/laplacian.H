#ifndef LAPLACIAN_H
#define LAPLACIAN_H

#include "amr-wind/fvm/fvm_utils.H"

namespace amr_wind {
namespace fvm {

template <typename FTypeIn, typename FTypeOut>
struct Laplacian
{
    Laplacian(FTypeOut& lphi, const FTypeIn& phi) : m_lapphi(lphi), m_phi(phi)
    {
        AMREX_ALWAYS_ASSERT(m_lapphi.num_comp() == m_phi.num_comp());
    }

    template <typename Stencil>
    void apply(const int lev, const amrex::MFIter& mfi) const
    {
        const int ncomp = m_phi.num_comp();
        const auto& geom = m_phi.repo().mesh().Geom(lev);
        const auto& idx = geom.InvCellSizeArray();
        const auto& lapphi = m_lapphi(lev).array(mfi);
        const auto& phi = m_phi(lev).const_array(mfi);

        const auto& bx_in = mfi.tilebox();
        const auto& bx = Stencil::box(bx_in, geom);
        if (bx.isEmpty()) return;

        amrex::ParallelFor(
            bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                amrex::Real sp1, s, sm1, d2phidx2, d2phidy2, d2phidz2;
                for (int icomp = 0; icomp < ncomp; icomp++) {
                    sp1 = Stencil::s00;
                    s = Stencil::s01;
                    sm1 = Stencil::s02;
                    d2phidx2 = (sp1 * phi(i + 1, j, k, icomp) +
                                s * phi(i, j, k, icomp) +
                                sm1 * phi(i - 1, j, k, icomp)) *
                               idx[0] * idx[0];
                    sp1 = Stencil::s10;
                    s = Stencil::s11;
                    sm1 = Stencil::s12;
                    d2phidy2 = (sp1 * phi(i, j + 1, k, icomp) +
                                s * phi(i, j, k, icomp) +
                                sm1 * phi(i, j - 1, k, icomp)) *
                               idx[1] * idx[1];
                    sp1 = Stencil::s20;
                    s = Stencil::s21;
                    sm1 = Stencil::s22;
                    d2phidz2 = (sp1 * phi(i, j, k + 1, icomp) +
                                s * phi(i, j, k, icomp) +
                                sm1 * phi(i, j, k - 1, icomp)) *
                               idx[2] * idx[2];
                    lapphi(i, j, k, icomp) = d2phidx2 + d2phidy2 + d2phidz2;
                }
            });
    }

    FTypeOut& m_lapphi;
    const FTypeIn& m_phi;
};

} // namespace fvm
} // namespace amr_wind

#endif /* LAPLACIAN_H */
