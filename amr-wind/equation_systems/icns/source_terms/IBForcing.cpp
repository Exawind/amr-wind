#include "amr-wind/equation_systems/icns/source_terms/IBForcing.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/immersed_boundary/IB.H"

#include "AMReX_Gpu.H"

namespace amr_wind {
namespace pde {
namespace icns {

IBForcing::IBForcing(const CFDSim& sim)
    : m_ib_src(sim.repo().get_field("ib_src_term"))
{}

IBForcing::~IBForcing() = default;

void IBForcing::operator()(
    const int lev,
    const amrex::MFIter& mfi,
    const amrex::Box& bx,
    const FieldState /*fstate*/,
    const amrex::Array4<amrex::Real>& src_term)
{
    const auto varr = m_ib_src(lev).const_array(mfi);

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        for (int ii = 0; ii < 3; ii++) {
            src_term(i, j, k, ii) += varr(i, j, k, ii);
        }
    });
}

} // namespace icns
} // namespace pde
} // namespace amr_wind