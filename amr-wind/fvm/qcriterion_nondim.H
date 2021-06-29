#ifndef QCRITERIONNONDIM_H
#define QCRITERIONNONDIM_H

#include "amr-wind/fvm/fvm_utils.H"

namespace amr_wind {
namespace fvm {

/** Qcriterion operator
 *  \ingroup fvm
 */
template <typename FTypeIn, typename FTypeOut>
struct QcriterionNondim
{
    QcriterionNondim(FTypeOut& qcritnomdimphi, const FTypeIn& phi)
        : m_qcritnomdimphi(qcritnomdimphi), m_phi(phi)
    {
        AMREX_ALWAYS_ASSERT(AMREX_SPACEDIM == m_phi.num_comp());
    }

    template <typename Stencil>
    void apply(const int lev, const amrex::MFIter& mfi) const
    {
        const auto& geom = m_phi.repo().mesh().Geom(lev);
        const auto& idx = geom.InvCellSizeArray();
        const auto& qcritnomdimphi = m_qcritnomdimphi(lev).array(mfi);
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

                const amrex::Real S2 =
                    2.0 * std::pow(ux, 2) + 2.0 * std::pow(vy, 2) +
                    2.0 * std::pow(wz, 2) + std::pow(uy + vx, 2) +
                    std::pow(vz + wy, 2) + std::pow(wx + uz, 2);

                const amrex::Real W2 = std::pow(uy - vx, 2) +
                                       std::pow(vz - wy, 2) +
                                       std::pow(wx - uz, 2);

                qcritnomdimphi(i, j, k) = 0.5 * (W2 / amrex::max(S2, 1.0e-12) - 1.0);
            });
    }

    FTypeOut& m_qcritnomdimphi;
    const FTypeIn& m_phi;
};

/** Compute the non-dimensional q-criterion
 *  \ingroup fvm
 *
 *  \param qcritnomdimphi [out] Field where the Q-criterion is populated
 *  \param phi [in] The velocity vector field
 */
template <typename FTypeIn, typename FTypeOut>
inline void q_criterion_nondim(FTypeOut& qcritnomdimphi, const FTypeIn& phi)
{
    BL_PROFILE("amr-wind::fvm::q_criterion_nondim");
    QcriterionNondim<FTypeIn, FTypeOut> qcriterionnondim(qcritnomdimphi, phi);
    impl::apply(qcriterionnondim, phi);
}

/** Compute the non-dimensional Q-criterion return as a ScratchField
 *  \ingroup fvm
 *
 *  \param phi [in] The velocity vector field
 */
template <typename FType>
inline std::unique_ptr<ScratchField> q_criterion_nondim(const FType& phi)
{
    const std::string gname = phi.name() + "_q_criterion_nondim";
    auto qcritnomdimphi = phi.repo().create_scratch_field(gname, 1);
    q_criterion_nondim(*qcritnomdimphi, phi);
    return qcritnomdimphi;
}

} // namespace fvm
} // namespace amr_wind

#endif /* QCRITERIONNONDIM_H */
