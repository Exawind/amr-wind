#include "amr-wind/immersed_boundary/bluff_body/bluff_body_ops.H"
#include "amr-wind/core/MultiParser.H"
#include "amr-wind/fvm/gradient.H"
#include "amr-wind/utilities/ncutils/nc_interface.H"
#include "amr-wind/utilities/io_utils.H"

// Used for mms
#include "amr-wind/physics/ConvectingTaylorVortex.H"

#include "AMReX_MultiFabUtil.H"
#include "AMReX_ParmParse.H"

namespace amr_wind {
namespace ib {
namespace bluff_body {

void read_inputs(
    BluffBodyBaseData& wdata, IBInfo&, const ::amr_wind::utils::MultiParser& pp)
{
    pp.query("has_wall_model", wdata.has_wall_model);
    pp.query("is_moving", wdata.is_moving);
    pp.query("is_mms", wdata.is_mms);
    pp.queryarr("vel_bc", wdata.vel_bc);
}

void init_data_structures(BluffBodyBaseData&) {}

void apply_mms_vel(CFDSim& sim)
{
    const int nlevels = sim.repo().num_active_levels();
    auto& mask_cell = sim.repo().get_int_field("mask_cell");
    auto& velocity = sim.repo().get_field("velocity");
    auto& m_conv_taylor_green =
        sim.physics_manager().get<ctv::ConvectingTaylorVortex>();

    const amrex::Real u0 = m_conv_taylor_green.get_u0();
    const amrex::Real v0 = m_conv_taylor_green.get_v0();
    const amrex::Real omega = m_conv_taylor_green.get_omega();

    amrex::Real t = sim.time().new_time();
    auto& geom = sim.mesh().Geom();

    for (int lev = 0; lev < nlevels; ++lev) {

        const auto& dx = geom[lev].CellSizeArray();
        const auto& problo = geom[lev].ProbLoArray();

        for (amrex::MFIter mfi(mask_cell(lev)); mfi.isValid(); ++mfi) {
            const auto& bx = mfi.growntilebox();
            auto epsilon_cell = mask_cell(lev).array(mfi);
            auto varr = velocity(lev).array(mfi);
            amrex::ParallelFor(
                bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    const amrex::Real x = problo[0] + (i + 0.5) * dx[0];
                    const amrex::Real y = problo[1] + (j + 0.5) * dx[1];
                    if (epsilon_cell(i, j, k) == 0) {
                        varr(i, j, k, 0) =
                            u0 - std::cos(utils::pi() * (x - u0 * t)) *
                                     std::sin(utils::pi() * (y - v0 * t)) *
                                     std::exp(-2.0 * omega * t);
                        varr(i, j, k, 1) =
                            v0 + std::sin(utils::pi() * (x - u0 * t)) *
                                     std::cos(utils::pi() * (y - v0 * t)) *
                                     std::exp(-2.0 * omega * t);
                        varr(i, j, k, 2) = 0.0;
                    }
                });
        }
    }
}

void apply_dirichlet_vel(CFDSim& sim, amrex::Vector<amrex::Real>& vel_bc)
{
    const int nlevels = sim.repo().num_active_levels();
    auto& mask_cell = sim.repo().get_int_field("mask_cell");
    auto& velocity = sim.repo().get_field("velocity");

    for (int lev = 0; lev < nlevels; ++lev) {
        for (amrex::MFIter mfi(mask_cell(lev)); mfi.isValid(); ++mfi) {
            const auto& bx = mfi.growntilebox();
            auto epsilon_cell = mask_cell(lev).array(mfi);
            auto varr = velocity(lev).array(mfi);
            amrex::ParallelFor(
                bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    if (epsilon_cell(i, j, k) == 0) {
                        varr(i, j, k, 0) = vel_bc[0];
                        varr(i, j, k, 1) = vel_bc[1];
                        varr(i, j, k, 2) = vel_bc[2];
                    }
                });
        }
    }
}

void prepare_netcdf_file(
    const std::string& ncfile,
    const BluffBodyBaseData& meta,
    const IBInfo& info)
{
    amrex::ignore_unused(ncfile, meta, info);
}

void write_netcdf(
    const std::string& ncfile,
    const BluffBodyBaseData& meta,
    const IBInfo& info,
    const amrex::Real time)
{
    amrex::ignore_unused(ncfile, meta, info, time);
}

} // namespace bluff_body

namespace ops {

    template <typename GeomTrait>
    void ComputeForceOp< GeomTrait, typename std::enable_if<
        std::is_base_of<BluffBodyType, GeomTrait>::value>::type>::operator()(
            typename GeomTrait::DataType& data)
    {
        auto& sim = data.sim();
        auto geom = sim.mesh().Geom();
        auto nlevels = sim.repo().num_active_levels();
        amrex::Real dt = sim.time().deltaT();

        // get relevant fields
        const auto& velocity = sim.repo().get_field("velocity");
        const auto& vel_np1  = velocity.state(FieldState::NP1);
        const auto& vel_n    = velocity.state(FieldState::N);
        const auto& vel_nm1  = velocity.state(FieldState::NM1);
        auto grad_vel = sim.repo().create_scratch_field(9);
        fvm::gradient(*grad_vel, vel_np1); // arranged row-wise

        const auto& pressure = sim.repo().get_field("p");
        auto p_cc = sim.repo().create_scratch_field(1);

        const auto& density = sim.repo().get_field("density");
        const auto& viscosity = sim.repo().get_field("velocity_mueff");

        const auto& mask_cell = sim.repo().get_int_field("mask_cell");

        amrex::RealBox bounding_box = data().info().bound_box;
        const amrex::Real* bb_lo = bounding_box.lo();
        const amrex::Real* bb_hi = bounding_box.hi();

        // forces on the immersed boundary as evaluated in all directions
        amrex::Real fx, fy, fz = 0.0;

        for (int lev = 0; lev < nlevels; ++lev) {
            const auto& problo = geom[lev].ProbLoArray();
            const auto& dx = geom[lev].CellSizeArray();

            // average pressure from node to cell centers
            amrex::average_node_to_cellcenter((*p_cc)(lev), 0, pressure(lev), 0, 1);

            for (amrex::MFIter mfi(velocity(lev)); mfi.isValid(); ++mfi) {

                // only process this box if it intersects with bounding box
                const auto& bx = mfi.tilebox();
                amrex::RealBox rbx(bx, dx, problo);
                if (!(rbx.intersects(bounding_box))) {
                    continue;
                }

                amrex::Array4<amrex::Real const> const& unp1   = vel_np1[lev]->const_array(mfi);
                amrex::Array4<amrex::Real const> const& un     = vel_np1[lev]->const_array(mfi);
                amrex::Array4<amrex::Real const> const& unm1   = vel_np1[lev]->const_array(mfi);
                amrex::Array4<amrex::Real const> const& grad_u = grad_vel[lev]->const_array(mfi);
                amrex::Array4<amrex::Real const> const& p      = p_cc[lev]->const_array(mfi);
                amrex::Array4<amrex::Real const> const& rho    = density[lev]->const_array(mfi);
                amrex::Array4<amrex::Real const> const& mu     = viscosity[lev]->const_array(mfi);
                amrex::Array4<amrex::Real const> const& mask   = mask_cell[lev]->const_array(mfi);

                amrex::ParallelFor(
                    bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        const amrex::Real x = problo[0] + (i + 0.5) * dx[0];
                        const amrex::Real y = problo[1] + (j + 0.5) * dx[1];
                        const amrex::Real z = problo[2] + (k + 0.5) * dx[2];

                        // compute forces only if cell is inside the bounding box
                        if((x-dx[0]/2 >= bb_lo[0]) && (x+dx[0]/2 <= bb_hi[0]) &&
                           (y-dx[1]/2 >= bb_lo[1]) && (y+dx[1]/2 <= bb_hi[1]) &&
                           (z-dx[2]/2 >= bb_lo[2]) && (z+dx[2]/2 <= bb_hi[2])) {

                            // add volume terms
                            const amrex::Real dv = dx[0]*dx[1]*dx[2];
                            const amrex::Real f_v_x =
                                    (1.5*unp1(i,j,k,0) - 2*un(i,j,k,0) + 0.5*unm1(i,j,k,0)) * dv/dt;
                            const amrex::Real f_v_y =
                                    (1.5*unp1(i,j,k,1) - 2*un(i,j,k,1) + 0.5*unm1(i,j,k,1)) * dv/dt;
                            const amrex::Real f_v_z =
                                    (1.5*unp1(i,j,k,2) - 2*un(i,j,k,2) + 0.5*unm1(i,j,k,2)) * dv/dt;

                            // add surface terms
                            amrex::Real f_s_x, f_s_y, f_s_z = 0.0;

                            if(((bb_lo[0] >= x-dx[0]/2) && (bb_lo[0] <= x+dx[0]/2)) ||
                               ((bb_hi[0] >= x-dx[0]/2) && (bb_hi[0] <= x+dx[0]/2))) { // y-z plane
                                const amrex::Real da =
                                    (bb_hi[0] <= x+dx[0]/2) ? (dx[1]*dx[2]) : (-dx[1]*dx[2]);

                                f_s_x += rho(i,j,k) * unp1(i,j,k,0) * unp1(i,j,k,0) * da;
                                f_s_y += rho(i,j,k) * unp1(i,j,k,1) * unp1(i,j,k,0) * da;
                                f_s_z += rho(i,j,k) * unp1(i,j,k,2) * unp1(i,j,k,0) * da;

                                f_s_x += rho(i,j,k) * (-p(i,j,k)
                                                    + 2*mu(i,j,k) * grad_u(i,j,k,0)) * da;
                                f_s_y += rho(i,j,k) * (-p(i,j,k)
                                                    + mu(i,j,k) * (grad_u(i,j,k,1) + grad_u(i,j,k,3))) * da;
                                f_s_z += rho(i,j,k) * (-p(i,j,k)
                                                    + mu(i,j,k) * (grad_u(i,j,k,2) + grad_u(i,j,k,6))) * da;
                            }
                            else if(((bb_lo[1] >= y-dx[1]/2) && (bb_lo[1] <= y+dx[1]/2))
                                    ((bb_hi[1] >= y-dx[1]/2) && (bb_hi[1] <= y+dx[1]/2))) { // x-z plane
                                const amrex::Real da =
                                    (bb_hi[1] <= y+dx[1]/2) ? (dx[0]*dx[2]) : (-dx[0]*dx[2]);

                                f_s_x += rho(i,j,k) * unp1(i,j,k,0) * unp1(i,j,k,1) * da;
                                f_s_y += rho(i,j,k) * unp1(i,j,k,1) * unp1(i,j,k,1) * da;
                                f_s_z += rho(i,j,k) * unp1(i,j,k,2) * unp1(i,j,k,1) * da;

                                f_s_x += rho(i,j,k) * (-p(i,j,k)
                                                    + mu(i,j,k) * (grad_u(i,j,k,1) + grad_u(i,j,k,3))) * da;
                                f_s_y += rho(i,j,k) * (-p(i,j,k)
                                                    + 2*mu(i,j,k) * grad_u(i,j,k,4)) * da;
                                f_s_z += rho(i,j,k) * (-p(i,j,k)
                                                    + mu(i,j,k) * (grad_u(i,j,k,5) + grad_u(i,j,k,7))) * da;

                            }
                            else if(((bb_lo[2] >= z-dx[2]/2) && (bb_lo[2] <= z+dx[2]/2))
                                    ((bb_hi[2] >= z-dx[2]/2) && (bb_hi[2] <= z+dx[2]/2))) { // x-y plane
                                const amrex::Real da =
                                    (bb_hi[2] <= z+dx[2]/2) ? (dx[0]*dx[1]) : (-dx[0]*dx[1]);

                                f_s_x += rho(i,j,k) * unp1(i,j,k,0) * unp1(i,j,k,2) * da;
                                f_s_y += rho(i,j,k) * unp1(i,j,k,1) * unp1(i,j,k,2) * da;
                                f_s_z += rho(i,j,k) * unp1(i,j,k,2) * unp1(i,j,k,2) * da;

                                f_s_x += rho(i,j,k) * (-p(i,j,k)
                                                    + mu(i,j,k) * (grad_u(i,j,k,2) + grad_u(i,j,k,6))) * da;
                                f_s_y += rho(i,j,k) * (-p(i,j,k)
                                                    + mu(i,j,k) * (grad_u(i,j,k,5) + grad_u(i,j,k,7))) * da;
                                f_s_z += rho(i,j,k) * (-p(i,j,k)
                                                    + 2*mu(i,j,k) * grad_u(i,j,k,8)) * da;

                            }
                        }
                    });
            }
        }
    };

} // namespace ops

} // namespace ib
} // namespace amr_wind
