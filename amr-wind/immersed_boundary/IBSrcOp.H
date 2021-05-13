#ifndef IBSRCOP_H
#define IBSRCOP_H

#include "amr-wind/immersed_boundary/IBTypes.H"
#include "amr-wind/immersed_boundary/IBOps.H"
#include "amr-wind/immersed_boundary/IBUtils.H"
#include "amr-wind/core/FieldRepo.H"

namespace amr_wind {
namespace ib {
namespace ops {

template <typename GeomTrait>
class IBSrcOp<GeomTrait>
{
private:
    typename GeomTrait::DataType& m_data;
    Field& m_ib_src;

    DeviceVecList m_pos;
    DeviceVecList m_forcing_vel;
    DeviceVecList m_intrp_vel;

    void copy_to_device();

public:
    IBSrcOp(typename GeomTrait::DataType& data)
        : m_data(data), m_ib_src(m_data.sim().repo().get_field("ib_src_term"))
    {}

    void initialize();

    void setup_op() { copy_to_device(); }

    void operator()(
        const int lev, const amrex::MFIter& mfi, const amrex::Geometry& geom);
};

template <typename GeomTrait>
void IBSrcOp<GeomTrait>::initialize()
{
    const auto& grid = m_data.grid();
    m_pos.resize(grid.pos.size());
    m_forcing_vel.resize(grid.forcing_vel.size());
    m_intrp_vel.resize(grid.intrp_vel.size());
}

template <typename GeomTrait>
void IBSrcOp<GeomTrait>::copy_to_device()
{
    const auto& grid = m_data.grid();

    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, grid.pos.begin(), grid.pos.end(),
        m_pos.begin());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, grid.forcing_vel.begin(),
        grid.forcing_vel.end(), m_forcing_vel.begin());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, grid.intrp_vel.begin(), grid.intrp_vel.end(),
        m_intrp_vel.begin());
}

template <typename GeomTrait>
void IBSrcOp<GeomTrait>::operator()(
    const int lev, const amrex::MFIter& mfi, const amrex::Geometry& geom)
{
    const std::string fname = GeomTrait::identifier();
    BL_PROFILE("amr-wind::IBSrcOp<" + fname + ">");

    const amrex::Real dt = m_data.sim().time().deltaT();

    const auto& bx = mfi.tilebox();
    const auto& problo = geom.ProbLoArray();
    const auto& dx = geom.CellSizeArray();
    const amrex::Real dV = dx[0] * dx[1] * dx[2];

    // extract relevant fields
    const auto* pos = m_pos.data();
    const auto* forcing_vel = m_forcing_vel.data();
    const auto* intrp_vel = m_intrp_vel.data();
    const auto& sarr = m_ib_src(lev).array(mfi);

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        const vs::Vector ijk_x{
            problo[0] + (i + 0.5) * dx[0], problo[1] + (j + 0.5) * dx[1],
            problo[2] + (k + 0.5) * dx[2]};

        // loop over Lagrangian markers
        amrex::Real src_force[AMREX_SPACEDIM]{0.0, 0.0, 0.0};
        for (unsigned ip = 0; ip < m_pos.size(); ++ip) {
            const auto lag_force = (forcing_vel[ip] - intrp_vel[ip]) / dt;
            const auto dist_x = ijk_x - pos[ip];

            const auto dirac_delta_fac =
                utils::dirac_delta(dist_x, {dx[0], dx[1], dx[2]});

            src_force[0] += lag_force[0] * dirac_delta_fac * dV;
            src_force[1] += lag_force[1] * dirac_delta_fac * dV;
            src_force[2] += lag_force[2] * dirac_delta_fac * dV;
        }
        sarr(i, j, k, 0) += src_force[0];
        sarr(i, j, k, 1) += src_force[1];
        sarr(i, j, k, 2) += src_force[2];
    });
}

} // namespace ops
} // namespace ib
} // namespace amr_wind

#endif /* IBSRCOP_H */
