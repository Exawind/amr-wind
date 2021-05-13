#include "amr-wind/equation_systems/icns/source_terms/IBForcing.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/wind_energy/actuator/Actuator.H"

#include "AMReX_Gpu.H"

namespace amr_wind {
namespace pde {
namespace icns {

IBForcing::IBForcing(const CFDSim& sim)
    : m_ib_src(sim.repo().get_field("ib_src_term"))
{
    if (!sim.physics_manager().contains("IB")) {
        amrex::Abort("IBForcing requires IB physics to be active");
    }
}

IBForcing::~IBForcing() = default;

void IBForcing::operator()(
    const int lev,
    const amrex::MFIter& mfi,
    const amrex::Box& bx,
    const FieldState,
    const amrex::Array4<amrex::Real>& src_term) const
{
    const auto varr = m_ib_src(lev).const_array(mfi);

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        src_term(i, j, k, 0) += varr(i, j, k, 0);
        src_term(i, j, k, 1) += varr(i, j, k, 1);
        src_term(i, j, k, 2) += varr(i, j, k, 2);
    });
}

} // namespace icns
} // namespace pde
} // namespace amr_wind
