
#include "amr-wind/boundary_conditions/MassInflowOutflowBC.H"

namespace amr_wind {
namespace {

} // namespace

MassInflowOutflowBC::MassInflowOutflowBC(Field& field, amrex::Orientation ori)
    : m_field(field), m_ori(ori)
{}

void MassInflowOutflowBC::operator()(Field& /*field*/, const FieldState /*rho_state*/)
{
    const auto& repo = m_field.repo();
    const auto& velocity = repo.get_field("velocity");
    // const auto bcvals = m_field.bc_values_device();
    const int ncomp = m_field.num_comp();
    // const int idx = static_cast<int>(m_ori);
    const int idim = m_ori.coordDir();
    const auto islow = m_ori.isLow();
    const auto ishigh = m_ori.isHigh();
    const int nlevels = m_field.repo().num_active_levels();
    const amrex::IntVect iv_dir = {idim == 0, idim == 1, idim == 2};

    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& domain = repo.mesh().Geom(lev).Domain();

        amrex::MFItInfo mfi_info{};
        if (amrex::Gpu::notInLaunchRegion()) {
            mfi_info.SetDynamic(true);
        }
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (amrex::MFIter mfi(m_field(lev), mfi_info); mfi.isValid(); ++mfi) {
            auto bx = mfi.validbox();
            bx.grow({!iv_dir[0], !iv_dir[1], !iv_dir[2]});
            const auto& bc_a = m_field(lev).array(mfi);
            const auto& vel = velocity(lev).array(mfi);

            if (islow && (bx.smallEnd(idim) == domain.smallEnd(idim))) {
                amrex::ParallelFor(
                    amrex::bdryLo(bx, idim),
                    [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        const amrex::IntVect iv = {i, j, k};
                        const amrex::IntVect ivm = iv - iv_dir;
                        if (vel(ivm[0], ivm[1], ivm[2], idim) < 0) {
                            for (int n = 0; n < ncomp; n++) {
                                bc_a(ivm[0], ivm[1], ivm[2], n) = bc_a(i, j, k, n);
                            }
                        }
                    });
            }

            if (ishigh && (bx.bigEnd(idim) == domain.bigEnd(idim))) {
                amrex::ParallelFor(
                    amrex::bdryHi(bx, idim),
                    [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        const amrex::IntVect iv = {i, j, k};
                        const amrex::IntVect ivm = iv - iv_dir;
                        if (vel(i, j, k, idim) > 0) {
                            for (int n = 0; n < ncomp; n++) {
                                bc_a(i, j, k, n) = bc_a(ivm[0], ivm[1], ivm[2], n);
                            }
                        }
                    });
            }
        }
    }
}

} // namespace amr_wind