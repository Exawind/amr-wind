#ifndef CURVATURE_H
#define CURVATURE_H

#include "amr-wind/fvm/fvm_utils.H"

namespace amr_wind {
namespace fvm {

template <typename FTypeIn, typename FTypeOut>
struct Curvature
{
    Curvature(FTypeOut& curphi, const FTypeIn& phi)
        : m_curphi(curphi), m_phi(phi)
    {
        AMREX_ALWAYS_ASSERT(m_phi.num_comp() == AMREX_SPACEDIM);
        AMREX_ALWAYS_ASSERT(m_curphi.num_comp() == 1);
    }

    template <typename Stencil>
    void apply(const int lev, const amrex::MFIter& mfi) const
    {
        const auto& geom = m_phi.repo().mesh().Geom(lev);
        const auto& idx = geom.InvCellSizeArray();
        const auto& curphi = m_curphi(lev).array(mfi);
        const auto& phi = m_phi(lev).const_array(mfi);

        const auto& bx_in = mfi.tilebox();
        const auto& bx = Stencil::box(bx_in, geom);
        if (bx.isEmpty()) return;

        amrex::ParallelFor(
            bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                amrex::Real cp1, c, cm1;

                cp1 = Stencil::c00;
                c = Stencil::c01;
                cm1 = Stencil::c02;
                const amrex::Real phix = phi(i, j, k, 0);
                const amrex::Real phixx =
                    (cp1 * phi(i + 1, j, k, 0) + c * phi(i, j, k, 0) +
                     cm1 * phi(i - 1, j, k, 0)) *
                    idx[0];

                cp1 = Stencil::c10;
                c = Stencil::c11;
                cm1 = Stencil::c12;
                const amrex::Real phiy = phi(i, j, k, 1);
                const amrex::Real phixy =
                    (cp1 * phi(i, j + 1, k, 0) + c * phi(i, j, k, 0) +
                     cm1 * phi(i, j - 1, k, 0)) *
                    idx[1];
                const amrex::Real phiyy =
                    (cp1 * phi(i, j + 1, k, 1) + c * phi(i, j, k, 1) +
                     cm1 * phi(i, j - 1, k, 1)) *
                    idx[1];

                cp1 = Stencil::c20;
                c = Stencil::c21;
                cm1 = Stencil::c22;
                const amrex::Real phiz = phi(i, j, k, 2);
                const amrex::Real phixz =
                    (cp1 * phi(i, j, k + 1, 0) + c * phi(i, j, k, 0) +
                     cm1 * phi(i, j, k - 1, 0)) *
                    idx[2];
                const amrex::Real phiyz =
                    (cp1 * phi(i, j, k + 1, 1) + c * phi(i, j, k, 1) +
                     cm1 * phi(i, j, k - 1, 1)) *
                    idx[2];
                const amrex::Real phizz =
                    (cp1 * phi(i, j, k + 1, 2) + c * phi(i, j, k, 2) +
                     cm1 * phi(i, j, k - 1, 2)) *
                    idx[2];

                curphi(i, j, k) =
                    -(phix * phix * phiyy - 2. * phix * phiy * phixy +
                      phiy * phiy * phixx + phix * phix * phizz -
                      2. * phix * phiz * phixz + phiz * phiz * phixx +
                      phiy * phiy * phizz - 2. * phiy * phiz * phiyz +
                      phiz * phiz * phiyy) /
                    std::pow(phix * phix + phiy * phiy + phiz * phiz, 1.5);
            });
    }

    FTypeOut& m_curphi;
    const FTypeIn& m_phi;
};

template <typename FTypeIn, typename FTypeOut>
inline void curvature(FTypeOut& curphi, const FTypeIn& phi)
{
    Curvature<FTypeIn, FTypeOut> cur(curphi, phi);
    impl::apply(cur, phi);
}

template <typename FType>
inline std::unique_ptr<ScratchField> curvature(const FType& phi)
{
    const std::string gname = phi.name() + "_curvature";
    auto curphi = phi.repo().create_scratch_field(gname, 1);
    curvature(*curphi, phi);
    return curphi;
}

} // namespace fvm
} // namespace amr_wind

#endif /* CURVATURE_H */
