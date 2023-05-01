#ifndef STRAINRATE_H
#define STRAINRATE_H

#include "amr-wind/fvm/fvm_utils.H"

namespace amr_wind::fvm {

/** Strain rate operator
 *  \ingroup fvm
 */
template <typename FTypeIn, typename FTypeOut>
struct StrainRate
{
    StrainRate(FTypeOut& strphi, const FTypeIn& phi)
        : m_strphi(strphi), m_phi(phi)
    {
        AMREX_ALWAYS_ASSERT(AMREX_SPACEDIM == m_phi.num_comp());
    }

    template <typename Stencil>
    void apply(const int lev, const amrex::MFIter& mfi) const
    {
        const auto& geom = m_phi.repo().mesh().Geom(lev);
        const auto& idx = geom.InvCellSizeArray();
        const auto& strphi = m_strphi(lev).array(mfi);
        const auto& phi = m_phi(lev).const_array(mfi);

        const auto& bx_in = mfi.tilebox();
        const auto& bx = Stencil::box(bx_in, geom);
        if (bx.isEmpty()) {
            return;
        }
        const auto& bctype = m_phi.bc_type();

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

                // --- Check for wall boundary BCs ---
                amrex::Real cp2 = 0.0, cm2 = 0.0;
                // xlo --> modifies gradients in x
                int dir = 0;
                amrex::Orientation oxlo(dir, amrex::Orientation::low);
                if (bctype[oxlo] == BC::wall_model &&
                    bx.smallEnd(dir) == geom.Domain().smallEnd(dir) &&
                    bx.bigEnd(dir) == geom.Domain().smallEnd(dir)) {
                    cp2 = Stencil::oc00;
                    cp1 = Stencil::oc01;
                    c = Stencil::oc02;
                    ux = (cp2 * phi(i + 2, j, k, 0) +
                          cp1 * phi(i + 1, j, k, 0) + c * phi(i, j, k, 0)) *
                         idx[0];
                    vx = (cp2 * phi(i + 2, j, k, 1) +
                          cp1 * phi(i + 1, j, k, 1) + c * phi(i, j, k, 1)) *
                         idx[0];
                    wx = (cp2 * phi(i + 2, j, k, 2) +
                          cp1 * phi(i + 1, j, k, 2) + c * phi(i, j, k, 2)) *
                         idx[0];
                }
                // xhi --> modifies gradients in x
                amrex::Orientation oxhi(dir, amrex::Orientation::high);
                if (bctype[oxhi] == BC::wall_model &&
                    bx.bigEnd(dir) == geom.Domain().bigEnd(dir) &&
                    bx.smallEnd(dir) == geom.Domain().bigEnd(dir)) {
                    c = Stencil::oc00;
                    cm1 = Stencil::oc01;
                    cm2 = Stencil::oc02;
                    ux = (c * phi(i, j, k, 0) + cm1 * phi(i - 1, j, k, 0) +
                          cm2 * phi(i - 2, j, k, 0)) *
                         idx[0];
                    vx = (c * phi(i, j, k, 1) + cm1 * phi(i - 1, j, k, 1) +
                          cm2 * phi(i - 2, j, k, 1)) *
                         idx[0];
                    wx = (c * phi(i, j, k, 2) + cm1 * phi(i - 1, j, k, 2) +
                          cm2 * phi(i - 2, j, k, 2)) *
                         idx[0];
                }
                // ylo --> modifies gradients in x
                dir = 1;
                amrex::Orientation oylo(dir, amrex::Orientation::low);
                if (bctype[oylo] == BC::wall_model &&
                    bx.smallEnd(dir) == geom.Domain().smallEnd(dir) &&
                    bx.bigEnd(dir) == geom.Domain().smallEnd(dir)) {
                    cp2 = Stencil::oc10;
                    cp1 = Stencil::oc11;
                    c = Stencil::oc12;
                    uy = (cp2 * phi(i, j + 2, k, 0) +
                          cp1 * phi(i, j + 1, k, 0) + c * phi(i, j, k, 0)) *
                         idx[1];
                    vy = (cp2 * phi(i, j + 2, k, 1) +
                          cp1 * phi(i, j + 1, k, 1) + c * phi(i, j, k, 1)) *
                         idx[1];
                    wy = (cp2 * phi(i, j + 2, k, 2) +
                          cp1 * phi(i, j + 1, k, 2) + c * phi(i, j, k, 2)) *
                         idx[1];
                }
                // yhi --> modifies gradients in y
                amrex::Orientation oyhi(dir, amrex::Orientation::high);
                if (bctype[oyhi] == BC::wall_model &&
                    bx.bigEnd(dir) == geom.Domain().bigEnd(dir) &&
                    bx.smallEnd(dir) == geom.Domain().bigEnd(dir)) {
                    c = Stencil::oc10;
                    cm1 = Stencil::oc11;
                    cm2 = Stencil::oc12;
                    uy = (c * phi(i, j, k, 0) + cm1 * phi(i, j - 1, k, 0) +
                          cm2 * phi(i, j - 2, k, 0)) *
                         idx[1];
                    vy = (c * phi(i, j, k, 1) + cm1 * phi(i, j - 1, k, 1) +
                          cm2 * phi(i, j - 2, k, 1)) *
                         idx[1];
                    wy = (c * phi(i, j, k, 2) + cm1 * phi(i, j - 1, k, 2) +
                          cm2 * phi(i, j - 2, k, 2)) *
                         idx[1];
                }
                // zlo --> modifies gradients in z
                dir = 2;
                amrex::Orientation ozlo(dir, amrex::Orientation::low);
                if (bctype[ozlo] == BC::wall_model &&
                    bx.smallEnd(dir) == geom.Domain().smallEnd(dir) &&
                    bx.bigEnd(dir) == geom.Domain().smallEnd(dir)) {
                    cp2 = Stencil::oc20;
                    cp1 = Stencil::oc21;
                    c = Stencil::oc22;
                    uz = (cp2 * phi(i, j, k + 2, 0) +
                          cp1 * phi(i, j, k + 1, 0) + c * phi(i, j, k, 0)) *
                         idx[2];
                    vz = (cp2 * phi(i, j, k + 2, 1) +
                          cp1 * phi(i, j, k + 1, 1) + c * phi(i, j, k, 1)) *
                         idx[2];
                    wz = (cp2 * phi(i, j, k + 2, 2) +
                          cp1 * phi(i, j, k + 1, 2) + c * phi(i, j, k, 2)) *
                         idx[2];
                }
                // zhi --> modifies gradients in z
                amrex::Orientation ozhi(dir, amrex::Orientation::high);
                if (bctype[ozhi] == BC::wall_model &&
                    bx.bigEnd(dir) == geom.Domain().bigEnd(dir) &&
                    bx.smallEnd(dir) == geom.Domain().bigEnd(dir)) {
                    c = Stencil::oc20;
                    cm1 = Stencil::oc21;
                    cm2 = Stencil::oc22;
                    uz = (c * phi(i, j, k, 0) + cm1 * phi(i, j, k - 1, 0) +
                          cm2 * phi(i, j, k - 2, 0)) *
                         idx[2];
                    vz = (c * phi(i, j, k, 1) + cm1 * phi(i, j, k - 1, 1) +
                          cm2 * phi(i, j, k - 2, 1)) *
                         idx[2];
                    wz = (c * phi(i, j, k, 2) + cm1 * phi(i, j, k - 1, 2) +
                          cm2 * phi(i, j, k - 2, 2)) *
                         idx[2];
                }

                // --- Check for slip wall BCs ---
                // xlo --> modifies gradients in x
                dir = 0;
                if (bctype[oxlo] == BC::slip_wall &&
                    bx.smallEnd(dir) == geom.Domain().smallEnd(dir) &&
                    bx.bigEnd(dir) == geom.Domain().smallEnd(dir)) {
                    cp1 = Stencil::nc00;
                    c = Stencil::nc01;
                    cm1 = Stencil::nc02;
                    vx = (cp1 * phi(i + 1, j, k, 1) + c * phi(i, j, k, 1)) *
                             idx[0] +
                         cm1 * phi(i - 1, j, k, 1);
                    wx = (cp1 * phi(i + 1, j, k, 2) + c * phi(i, j, k, 2)) *
                             idx[0] +
                         cm1 * phi(i - 1, j, k, 2);
                }
                // xhi --> modifies gradients in x
                if (bctype[oxhi] == BC::slip_wall &&
                    bx.bigEnd(dir) == geom.Domain().bigEnd(dir) &&
                    bx.smallEnd(dir) == geom.Domain().bigEnd(dir)) {
                    cp1 = Stencil::nc00;
                    c = Stencil::nc01;
                    cm1 = Stencil::nc02;
                    vx = cp1 * phi(i + 1, j, k, 1) +
                         (c * phi(i, j, k, 1) + cm1 * phi(i - 1, j, k, 1)) *
                             idx[0];
                    wx = cp1 * phi(i + 1, j, k, 2) +
                         (c * phi(i, j, k, 2) + cm1 * phi(i - 1, j, k, 2)) *
                             idx[0];
                }
                // ylo --> modifies gradients in x
                dir = 1;
                if (bctype[oylo] == BC::slip_wall &&
                    bx.smallEnd(dir) == geom.Domain().smallEnd(dir) &&
                    bx.bigEnd(dir) == geom.Domain().smallEnd(dir)) {
                    cp1 = Stencil::nc10;
                    c = Stencil::nc11;
                    cm1 = Stencil::nc12;
                    uy = (cp1 * phi(i, j + 1, k, 0) + c * phi(i, j, k, 0)) *
                             idx[1] +
                         cm1 * phi(i, j - 1, k, 0);
                    wy = (cp1 * phi(i, j + 1, k, 2) + c * phi(i, j, k, 2)) *
                             idx[1] +
                         cm1 * phi(i, j - 1, k, 2);
                }
                // yhi --> modifies gradients in y
                if (bctype[oyhi] == BC::slip_wall &&
                    bx.bigEnd(dir) == geom.Domain().bigEnd(dir) &&
                    bx.smallEnd(dir) == geom.Domain().bigEnd(dir)) {
                    cp1 = Stencil::nc10;
                    c = Stencil::nc11;
                    cm1 = Stencil::nc12;
                    uy = cp1 * phi(i, j + 1, k, 0) +
                         (c * phi(i, j, k, 0) + cm1 * phi(i, j - 1, k, 0)) *
                             idx[1];
                    wy = cp1 * phi(i, j + 1, k, 2) +
                         (c * phi(i, j, k, 2) + cm1 * phi(i, j - 1, k, 2)) *
                             idx[1];
                }
                // zlo --> modifies gradients in z
                dir = 2;
                if (bctype[ozlo] == BC::slip_wall &&
                    bx.smallEnd(dir) == geom.Domain().smallEnd(dir) &&
                    bx.bigEnd(dir) == geom.Domain().smallEnd(dir)) {
                    cp1 = Stencil::nc20;
                    c = Stencil::nc21;
                    cm1 = Stencil::nc22;
                    uz = (cp1 * phi(i, j, k + 1, 0) + c * phi(i, j, k, 0)) *
                             idx[2] +
                         cm1 * phi(i, j, k - 1, 0);
                    vz = (cp1 * phi(i, j, k + 1, 1) + c * phi(i, j, k, 1)) *
                             idx[2] +
                         cm1 * phi(i, j, k - 1, 1);
                }
                // zhi --> modifies gradients in z
                if (bctype[ozhi] == BC::slip_wall &&
                    bx.bigEnd(dir) == geom.Domain().bigEnd(dir) &&
                    bx.smallEnd(dir) == geom.Domain().bigEnd(dir)) {
                    cp1 = Stencil::nc20;
                    c = Stencil::nc21;
                    cm1 = Stencil::nc22;
                    uz = cp1 * phi(i, j, k + 1, 0) +
                         (c * phi(i, j, k, 0) + cm1 * phi(i, j, k - 1, 0)) *
                             idx[2];
                    vz = cp1 * phi(i, j, k + 1, 1) +
                         (c * phi(i, j, k, 1) + cm1 * phi(i, j, k - 1, 1)) *
                             idx[2];
                }

                strphi(i, j, k) = sqrt(
                    2.0 * std::pow(ux, 2) + 2.0 * std::pow(vy, 2) +
                    2.0 * std::pow(wz, 2) + std::pow(uy + vx, 2) +
                    std::pow(vz + wy, 2) + std::pow(wx + uz, 2));
            });
    }

    FTypeOut& m_strphi;
    const FTypeIn& m_phi;
};

/** Compute the magnitude of strain rate
 *  \ingroup fvm
 *
 *  \param strphi [out] Field where the strain rate magnitude is populated
 *  \param phi [in] The velocity vector field
 */
template <typename FTypeIn, typename FTypeOut>
inline void strainrate(FTypeOut& strphi, const FTypeIn& phi)
{
    BL_PROFILE("amr-wind::fvm::strainrate");
    StrainRate<FTypeIn, FTypeOut> str(strphi, phi);
    impl::apply(str, phi);
}

/** Compute the magnitude of strain rate return as a ScratchField
 *  \ingroup fvm
 *
 *  \param phi [in] The velocity vector field
 */
template <typename FType>
inline std::unique_ptr<ScratchField> strainrate(const FType& phi)
{
    const std::string gname = phi.name() + "_strainrate";
    auto strphi = phi.repo().create_scratch_field(gname, 1);
    strainrate(*strphi, phi);
    return strphi;
}

} // namespace amr_wind::fvm

#endif /* STRAINRATE_H */
