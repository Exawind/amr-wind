#ifndef VORTICITY_H
#define VORTICITY_H

#include "amr-wind/fvm/fvm_utils.H"

namespace amr_wind {
namespace fvm {

/** Vorticity operator
 *  \ingroup fvm
 */
template <typename FTypeIn, typename FTypeOut>
struct Vorticity
{
    Vorticity(FTypeOut& vortphi, const FTypeIn& phi) : m_vortphi(vortphi), m_phi(phi)
    {
        AMREX_ALWAYS_ASSERT(AMREX_SPACEDIM == m_phi.num_comp());
    }

    template <typename Stencil>
    void apply(const int lev, const amrex::MFIter& mfi) const
    {
        const auto& geom = m_phi.repo().mesh().Geom(lev);
        const auto& idx = geom.InvCellSizeArray();
        const auto& vortphi = m_vortphi(lev).array(mfi);
        const auto& phi = m_phi(lev).const_array(mfi);

        const auto& bx_in = mfi.tilebox();
        const auto& bx = Stencil::box(bx_in, geom);
        if (bx.isEmpty()) return;

        amrex::ParallelFor(
            bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                amrex::Real cp1, c, cm1, ux, uy, uz, vx, vy, vz, wx, wy, wz;
                cp1 = Stencil::c00;
                c = Stencil::c01;
                cm1 = Stencil::c02;

                ux = (cp1 * phi(i + 1, j, k, 0) + c * phi(i, j, k, 0) +
                      cm1 * phi(i - 1, j, k, 0)) *
                     idx[0];
                vx = (cp1 * phi(i + 1, j, k, 1) + c * phi(i, j, k, 1) +
                      cm1 * phi(i - 1, j, k, 1)) *
                     idx[0];
                wx = (cp1 * phi(i + 1, j, k, 2) + c * phi(i, j, k, 2) +
                      cm1 * phi(i - 1, j, k, 2)) *
                     idx[0];

                cp1 = Stencil::c10;
                c = Stencil::c11;
                cm1 = Stencil::c12;

                uy = (cp1 * phi(i, j + 1, k, 0) + c * phi(i, j, k, 0) +
                      cm1 * phi(i, j - 1, k, 0)) *
                     idx[1];
                vy = (cp1 * phi(i, j + 1, k, 1) + c * phi(i, j, k, 1) +
                      cm1 * phi(i, j - 1, k, 1)) *
                     idx[1];
                wy = (cp1 * phi(i, j + 1, k, 2) + c * phi(i, j, k, 2) +
                      cm1 * phi(i, j - 1, k, 2)) *
                     idx[1];

                cp1 = Stencil::c20;
                c = Stencil::c21;
                cm1 = Stencil::c22;

                uz = (cp1 * phi(i, j, k + 1, 0) + c * phi(i, j, k, 0) +
                      cm1 * phi(i, j, k - 1, 0)) *
                     idx[2];
                vz = (cp1 * phi(i, j, k + 1, 1) + c * phi(i, j, k, 1) +
                      cm1 * phi(i, j, k - 1, 1)) *
                     idx[2];
                wz = (cp1 * phi(i, j, k + 1, 2) + c * phi(i, j, k, 2) +
                      cm1 * phi(i, j, k - 1, 2)) *
                     idx[2];

                vortphi(i, j, k) = sqrt(std::pow(uy - vx, 2)
                                       + std::pow(vz - wy, 2)
                                       + std::pow(wx - uz, 2));
            });
    }

    FTypeOut& m_vortphi;
    const FTypeIn& m_phi;
};

/** Compute the magnitude of vorticity
 *  \ingroup fvm
 *
 *  \param vortphi [out] Field where the vorticity magnitude is populated
 *  \param phi [in] The velocity vector field
 */
template <typename FTypeIn, typename FTypeOut>
inline void vorticity(FTypeOut& vortphi, const FTypeIn& phi)
{
    BL_PROFILE("amr-wind::fvm::vorticity");
    Vorticity<FTypeIn, FTypeOut> vort(vortphi, phi);
    impl::apply(vort, phi);
}

/** Compute the magnitude of vorticity return as a ScratchField
 *  \ingroup fvm
 *
 *  \param phi [in] The velocity vector field
 */
template <typename FType>
inline std::unique_ptr<ScratchField> vorticity(const FType& phi)
{
    const std::string gname = phi.name() + "_vorticity";
    auto vortphi = phi.repo().create_scratch_field(gname, 1);
    vorticity(*vortphi, phi);
    return vortphi;
}

} // namespace fvm
} // namespace amr_wind

#endif /* VORTICITY_H */
