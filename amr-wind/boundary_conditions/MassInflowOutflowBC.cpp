
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
    //const auto bcvals = m_field.bc_values_device();
    const int ncomp = m_field.num_comp();
    const int idx = static_cast<int>(m_ori);
    const int idim = m_ori.coordDir();
    const auto islow = m_ori.isLow();
    const auto ishigh = m_ori.isHigh();
    const int nlevels = m_field.repo().num_active_levels();
    const bool ib = (idim == 0), jb = (idim == 1), kb = (idim == 2);

    amrex::Print() << "******* applying MIO fill-patch at orientation: " << idx << std::endl;
    amrex::Print() << nlevels << " level(s)" << std::endl << std::endl;
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
            const auto& bx = mfi.validbox();
            const auto& bc_a = m_field(lev).array(mfi);
            const auto& vel = velocity(lev).array(mfi);

            if (islow && (bx.smallEnd(idim) == domain.smallEnd(idim))) {
                amrex::ParallelFor(
                    amrex::bdryLo(bx, idim), ncomp,
                    [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                    //amrex::Print() << i << " " << j << " " << k << std::endl;
                    //amrex::Print() << i-ib << " " << j-jb << " " << k-kb << std::endl;
                    //amrex::Print() << vel(i-ib, j-jb, k-kb, idim) << " " << bc_a(i-ib, j-jb, k-kb, n) << " " << bc_a(i, j, k, n) << std::endl << std::endl;
                        if (vel(i-ib, j-jb, k-kb, idim) < 0) {
                                bc_a(i-ib, j-jb, k-kb, n) = bc_a(i, j, k, n);
                        }
                    });
            }

            if (ishigh && (bx.bigEnd(idim) == domain.bigEnd(idim))) {
                amrex::ParallelFor(
                    amrex::bdryHi(bx, idim), ncomp,
                    [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                        if (vel(i, j, k, idim) > 0) {
                            bc_a(i, j, k, n) = bc_a(i-ib, j-jb, k-kb, n);
                        }
                    });
            }
        }
    }
}

} // namespace amr_wind