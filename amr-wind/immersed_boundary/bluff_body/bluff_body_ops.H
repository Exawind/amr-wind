#ifndef BLUFF_BODY_OPS_H
#define BLUFF_BODY_OPS_H

#include "amr-wind/immersed_boundary/bluff_body/BluffBody.H"
#include "amr-wind/immersed_boundary/IBTypes.H"
#include "amr-wind/immersed_boundary/IBOps.H"
#include "amr-wind/core/MultiParser.H"
#include "amr-wind/fvm/gradient.H"

#include "AMReX_MultiFabUtil.H"

namespace amr_wind {
namespace ib {
namespace bluff_body {

/** Parse common data for modeling non-blade bluff bodies via an immersed
 * bounary
 */
void read_inputs(
    BluffBodyBaseData&, IBInfo&, const ::amr_wind::utils::MultiParser&);

/** Initialize core data structures when modeling the IB geometry
 */
void init_data_structures(CFDSim&);

/** Set the velocity inside the IB based on a manufactured solution
 */
void apply_mms_vel(CFDSim&);

/** Set the velocity inside the IB based on a dirichlet BC
 */
void apply_dirichlet_vel(CFDSim&, amrex::Vector<amrex::Real>& vel_bc);

void prepare_netcdf_file(
    const std::string&, const BluffBodyBaseData&, const IBInfo&);

void write_netcdf(
    const std::string&,
    const BluffBodyBaseData&,
    const IBInfo&,
    const amrex::Real);

} // namespace bluff_body

namespace ops {

template <typename GeomTrait>
struct UseDefaultOp
{
    static constexpr bool update_pos = true;
    static constexpr bool update_vel = true;
    static constexpr bool compute_force = true;
    static constexpr bool process_outputs = true;
};

template <typename GeomTrait>
struct UpdatePosOp<
    GeomTrait,
    typename std::enable_if<
        std::is_base_of<BluffBodyType, GeomTrait>::value>::type>
{
    void operator()(typename GeomTrait::DataType&) {}
};

template <typename GeomTrait>
struct UpdateVelOp<
    GeomTrait,
    typename std::enable_if<
        std::is_base_of<BluffBodyType, GeomTrait>::value>::type>
{
    void operator()(typename GeomTrait::DataType& data)
    {
        BL_PROFILE("amr-wind::ib::IB::update_velocities");

        auto& wdata = data.meta();
        auto& sim = data.sim();

        if (wdata.is_mms) {
            bluff_body::apply_mms_vel(sim);
        } else {
            bluff_body::apply_dirichlet_vel(sim, wdata.vel_bc);
        }
    }
};

template <typename GeomTrait>
struct ComputeForceOp<
    GeomTrait,
    typename std::enable_if<
        std::is_base_of<BluffBodyType, GeomTrait>::value>::type>
{
    void operator()(typename GeomTrait::DataType& data)
    {
        auto& info = data.info();
        auto& wdata = data.meta();
        
        auto& sim = data.sim();
        auto geom = sim.mesh().Geom();
        auto nlevels = sim.repo().num_active_levels();
        amrex::Real dt = sim.time().deltaT();

        // get relevant fields
        const auto& velocity = sim.repo().get_field("velocity");
        const auto& vel_np1 = velocity.state(FieldState::NP1);
        const auto& vel_n = velocity.state(FieldState::N);
        auto grad_vel = sim.repo().create_scratch_field(9);
        fvm::gradient(*grad_vel, vel_np1); // arranged row-wise

        const auto& pressure = sim.repo().get_field("p");
        auto p_cc = sim.repo().create_scratch_field(1);

        // scratch field for cell-centered force on immersed boundary
        auto f_on_ib = sim.repo().create_scratch_field(3);

        const auto& density = sim.repo().get_field("density");
        const auto& viscosity = sim.repo().get_field("velocity_mueff");

        // TODO: We should be using a different field for ib_mask because
        // current field is shared with overset iblank, and this framework
        // may break if overset and immersed boundary are used together
        const auto& mask_cell = sim.repo().get_int_field("mask_cell");

        // get ib bounding box information
        amrex::RealBox bounding_box = info.bound_box;
        const amrex::Real* bb_lo = bounding_box.lo();
        const amrex::Real* bb_hi = bounding_box.hi();

        // zero out forces
        wdata.frc = {{0.0, 0.0, 0.0}};

        for (int lev = 0; lev < nlevels; ++lev) {
            // extract geometry information at this level
            const amrex::Real* problo = geom[lev].ProbLo();
            const amrex::Real* dx = geom[lev].CellSize();

            // average pressure from node to cell centers
            amrex::average_node_to_cellcenter(
                (*p_cc)(lev), 0, pressure(lev), 0, 1);

            for (amrex::MFIter mfi(velocity(lev)); mfi.isValid(); ++mfi) {
                // only process this box if it intersects with bounding box
                const auto& bx = mfi.tilebox();
                amrex::RealBox rbx(bx, dx, problo);
                if (!(rbx.intersects(bounding_box))) {
                    continue;
                }

                amrex::Array4<amrex::Real> const& f_arr =
                        (*f_on_ib)(lev).array(mfi);
                amrex::Array4<amrex::Real const> const& grad_u =
                        (*grad_vel)(lev).array(mfi);
                amrex::Array4<int const> const& ib_mask =
                        mask_cell(lev).const_array(mfi);
                amrex::Array4<amrex::Real const> const& mu =
                        viscosity(lev).const_array(mfi);
                amrex::Array4<amrex::Real const> const& p =
                        (*p_cc)(lev).const_array(mfi);
                amrex::Array4<amrex::Real const> const& rho =
                        density(lev).const_array(mfi);
                amrex::Array4<amrex::Real const> const& unp1 =
                        vel_np1(lev).const_array(mfi);
                amrex::Array4<amrex::Real const> const& un =
                        vel_n(lev).const_array(mfi);

                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    const amrex::Real x = problo[0] + (i + 0.5) * dx[0];
                    const amrex::Real y = problo[1] + (j + 0.5) * dx[1];
                    const amrex::Real z = problo[2] + (k + 0.5) * dx[2];

                    // compute forces only if cell is inside the bounding
                    // box Newton's law -> rho*du/dt * dV + rho * u*(u.n) *
                    // dA = tau.n * dA - F F is -ve above because it is the
                    // force acting on the immersed boundary
                    if ((x - dx[0] / 2 >= bb_lo[0]) &&
                        (x + dx[0] / 2 <= bb_hi[0]) &&
                        (y - dx[1] / 2 >= bb_lo[1]) &&
                        (y + dx[1] / 2 <= bb_hi[1]) &&
                        (z - dx[2] / 2 >= bb_lo[2]) &&
                        (z + dx[2] / 2 <= bb_hi[2])) {

                        // add volume terms -> rho*du/dt * dV
                        const amrex::Real dv = dx[0] * dx[1] * dx[2];
                        const amrex::Real f_vol_x =
                                (unp1(i, j, k, 0) - un(i, j, k, 0)) * dv / dt;
                        const amrex::Real f_vol_y =
                                (unp1(i, j, k, 1) - un(i, j, k, 1)) * dv / dt;
                        const amrex::Real f_vol_z =
                                (unp1(i, j, k, 2) - un(i, j, k, 2)) * dv / dt;

                        // initialize surface forces
                        amrex::Real f_vel_x = 0.0, f_vel_y = 0.0, f_vel_z = 0.0;
                        amrex::Real f_tau_x = 0.0, f_tau_y = 0.0, f_tau_z = 0.0;

                        if (((bb_lo[0] >= x - dx[0] / 2) &&
                             (bb_lo[0] <= x + dx[0] / 2)) ||
                            ((bb_hi[0] >= x - dx[0] / 2) &&
                             (bb_hi[0] <= x + dx[0] / 2))) { // y-z plane
                            // determine dS while accounting for direction
                            // of normal
                            const amrex::Real da = (bb_hi[0] <= x + dx[0] / 2) ?
                                                   (dx[1] * dx[2]) : (-dx[1] * dx[2]);

                            // add advection surface force -> rho * u*(u.n)
                            // * dA
                            f_vel_x = rho(i, j, k) * unp1(i, j, k, 0) * unp1(i, j, k, 0) * da;
                            f_vel_y = rho(i, j, k) * unp1(i, j, k, 1) * unp1(i, j, k, 0) * da;
                            f_vel_z = rho(i, j, k) * unp1(i, j, k, 2) * unp1(i, j, k, 0) * da;

                            // add viscous surface force -> tau.n * dA
                            f_tau_x = rho(i, j, k)
                                    * (-p(i, j, k) + 2 * mu(i, j, k) * grad_u(i, j, k, 0))
                                    * da;
                            f_tau_y = rho(i, j, k)
                                    * (-p(i, j, k) + mu(i, j, k) * (grad_u(i, j, k, 1) +
                                                                    grad_u(i, j, k, 3)))
                                    * da;
                            f_tau_z = rho(i, j, k)
                                    * (-p(i, j, k) + mu(i, j, k) * (grad_u(i, j, k, 2) +
                                                                    grad_u(i, j, k, 6)))
                                    * da;
                        } else if (((bb_lo[1] >= y - dx[1] / 2) &&
                                    (bb_lo[1] <= y + dx[1] / 2)) ||
                                   ((bb_hi[1] >= y - dx[1] / 2) &&
                                    (bb_hi[1] <= y + dx[1] / 2))) { // x-z plane
                            // determine dS while accounting for direction
                            // of normal
                            const amrex::Real da = (bb_hi[1] <= y + dx[1] / 2) ?
                                (dx[0] * dx[2]) : (-dx[0] * dx[2]);

                            // add advection surface force -> rho *u*(u.n) *
                            // dA
                            f_vel_x = rho(i, j, k) * unp1(i, j, k, 0) * unp1(i, j, k, 1) * da;
                            f_vel_y = rho(i, j, k) * unp1(i, j, k, 1) * unp1(i, j, k, 1) * da;
                            f_vel_z = rho(i, j, k) * unp1(i, j, k, 2) * unp1(i, j, k, 1) * da;

                            // add viscous surface force -> tau.n * dA
                            f_tau_x = rho(i, j, k)
                                    * (-p(i, j, k) + mu(i, j, k) * (grad_u(i, j, k, 1) +
                                                                    grad_u(i, j, k, 3)))
                                    * da;
                            f_tau_y = rho(i, j, k)
                                    * (-p(i, j, k) + 2 * mu(i, j, k) * grad_u(i, j, k, 4))
                                    * da;
                            f_tau_z = rho(i, j, k)
                                    * (-p(i, j, k) + mu(i, j, k) * (grad_u(i, j, k, 5) +
                                                                    grad_u(i, j, k, 7)))
                                    * da;

                        } else if (((bb_lo[2] >= z - dx[2] / 2) &&
                                    (bb_lo[2] <= z + dx[2] / 2)) ||
                                   ((bb_hi[2] >= z - dx[2] / 2) &&
                                    (bb_hi[2] <= z + dx[2] / 2))) { // x-y plane
                            // determine dS while accounting for direction
                            // of normal
                            const amrex::Real da = (bb_hi[2] <= z + dx[2] / 2) ?
                                (dx[0] * dx[1]) : (-dx[0] * dx[1]);

                            // add advection surface force -> rho * u*(u.n)
                            // * dA
                            f_vel_x = rho(i, j, k) * unp1(i, j, k, 0) * unp1(i, j, k, 2) * da;
                            f_vel_y = rho(i, j, k) * unp1(i, j, k, 1) * unp1(i, j, k, 2) * da;
                            f_vel_z = rho(i, j, k) * unp1(i, j, k, 2) * unp1(i, j, k, 2) * da;

                            // add viscous surface force -> tau.n * dA
                            f_tau_x = rho(i, j, k)
                                    * (-p(i, j, k) + mu(i, j, k) * (grad_u(i, j, k, 2) +
                                                                    grad_u(i, j, k, 6)))
                                    * da;
                            f_tau_y = rho(i, j, k)
                                    * (-p(i, j, k) + mu(i, j, k) * (grad_u(i, j, k, 5) +
                                                                    grad_u(i, j, k, 7)))
                                    * da;
                            f_tau_z = rho(i, j, k)
                                    * (-p(i, j, k) + 2 * mu(i, j, k) * grad_u(i, j, k, 8))
                                    * da;
                        }
                        // fill force scratch field
                        f_arr(i, j, k, 0) = (f_tau_x - f_vel_x - f_vol_x) * ib_mask(i, j, k);
                        f_arr(i, j, k, 1) = (f_tau_y - f_vel_y - f_vol_y) * ib_mask(i, j, k);
                        f_arr(i, j, k, 2) = (f_tau_z - f_vel_z - f_vol_z) * ib_mask(i, j, k);
                    }
                });
            }

            // assign level mask to make sure that integration doesn't
            // double count
            amrex::iMultiFab level_mask;
            if (lev < nlevels - 1) {
                level_mask = makeFineMask(
                    sim.mesh().boxArray(lev), sim.mesh().DistributionMap(lev),
                    sim.mesh().boxArray(lev + 1), amrex::IntVect(2), 1, 0);
            } else {
                level_mask.define(
                    sim.mesh().boxArray(lev), sim.mesh().DistributionMap(lev),
                    1, 0, amrex::MFInfo());
                level_mask.setVal(1);
            }

            // sum across all GPUs for current level
            //          const auto& f_lev = (*f_on_ib)(lev);
            wdata.frc[0] += amrex::ReduceSum(
                (*f_on_ib)(lev), level_mask, 0,
                [=] AMREX_GPU_HOST_DEVICE(
                    amrex::Box const& bx,
                    amrex::Array4<amrex::Real const> const& f_arr,
                    amrex::Array4<int const> const& mask_arr) -> amrex::Real {
                amrex::Real force = 0.0;
                amrex::Loop(bx, [=, &force](int i, int j, int k) noexcept {
                    force += f_arr(i, j, k, 0) * mask_arr(i, j, k);
                });
                return force;
            });

            wdata.frc[1] += amrex::ReduceSum(
                (*f_on_ib)(lev), level_mask, 0,
                [=] AMREX_GPU_HOST_DEVICE(
                    amrex::Box const& bx,
                    amrex::Array4<amrex::Real const> const& f_arr,
                    amrex::Array4<int const> const& mask_arr) -> amrex::Real {
                amrex::Real force = 0.0;
                amrex::Loop(bx, [=, &force](int i, int j, int k) noexcept {
                    force += f_arr(i, j, k, 1) * mask_arr(i, j, k);
                });
                return force;
            });

            wdata.frc[2] += amrex::ReduceSum(
                (*f_on_ib)(lev), level_mask, 0,
                [=] AMREX_GPU_HOST_DEVICE(
                    amrex::Box const& bx,
                    amrex::Array4<amrex::Real const> const& f_arr,
                    amrex::Array4<int const> const& mask_arr) -> amrex::Real {
                amrex::Real force = 0.0;
                amrex::Loop(bx, [=, &force](int i, int j, int k) noexcept {
                    force += f_arr(i, j, k, 2) * mask_arr(i, j, k);
                });
                return force;
            });
        }

        // sum across all MPI ranks
        amrex::ParallelDescriptor::ReduceRealSum(wdata.frc[0]);
        amrex::ParallelDescriptor::ReduceRealSum(wdata.frc[1]);
        amrex::ParallelDescriptor::ReduceRealSum(wdata.frc[2]);

        amrex::Print() << "Forces on immersed boundary : " << info.label
                       << " --> Fx, Fy, Fz = " << wdata.frc[0] << ", "
                       << wdata.frc[1] << ", " << wdata.frc[2] << std::endl;
    }
};

template <typename GeomTrait>
struct ProcessOutputsOp<
    GeomTrait,
    typename std::enable_if<
        std::is_base_of<BluffBodyType, GeomTrait>::value>::type>
{
private:
    typename GeomTrait::DataType& m_data;

    std::string m_out_dir{""};
    std::string m_nc_filename{""};
    int m_out_freq{10};

public:
    explicit ProcessOutputsOp(typename GeomTrait::DataType& data) : m_data(data)
    {}

    void read_io_options(const ::amr_wind::utils::MultiParser& pp)
    {
        pp.query("output_frequency", m_out_freq);
    }

    void prepare_outputs(const std::string& out_dir)
    {
        m_nc_filename = out_dir + "/" + m_data.info().label + ".nc";
        bluff_body::prepare_netcdf_file(
            m_nc_filename, m_data.meta(), m_data.info());
    }

    void write_outputs()
    {
        const auto& time = m_data.sim().time();
        const int tidx = time.time_index();
        if (tidx % m_out_freq != 0) return;

        bluff_body::write_netcdf(
            m_nc_filename, m_data.meta(), m_data.info(), time.new_time());
    }
};

} // namespace ops

} // namespace ib
} // namespace amr_wind

#endif /* BLUFF_BODY_OPS_H */
