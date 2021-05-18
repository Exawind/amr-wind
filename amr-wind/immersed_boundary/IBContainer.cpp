#include "amr-wind/immersed_boundary/IBContainer.H"
#include "amr-wind/immersed_boundary/IB.H"
#include "amr-wind/immersed_boundary/IBUtils.H"
#include "amr-wind/core/gpu_utils.H"
#include "amr-wind/core/Field.H"

#include "AMReX_Scan.H"

#include <algorithm>

namespace amr_wind {
namespace ib {

IBCloud::IBCloud(const int nobjects)
    : num_pts(nobjects, 0), global_id(nobjects, -1), num_objects(nobjects)
{}

IBContainer::IBContainer(amrex::AmrCore& mesh, const int num_objects)
    : amrex::AmrParticleContainer<
          NumPStructReal,
          NumPStructInt,
          NumPArrayReal,
          NumPArrayInt>(&mesh)
    , m_mesh(mesh)
    , m_data(num_objects)
    , m_proc_pos(amrex::ParallelDescriptor::NProcs(), vs::Vector::zero())
    , m_pos_device(amrex::ParallelDescriptor::NProcs(), vs::Vector::zero())
    , m_proc_offsets(amrex::ParallelDescriptor::NProcs() + 1, 0)
    , m_proc_offsets_device(amrex::ParallelDescriptor::NProcs() + 1)
{}

/** Allocate memory and initialize the particles within the container
 *
 *  This method is only called once during the simulation. It allocates the
 *  arrays for holding the position vector and velocity data on host memory and
 *  also initializes corresponding particles within the first available particle
 *  tile within the container. It is expected that the immersed boundary manager
 * instance has already populated the number of points per turbine before
 * invoking this method.
 */
void IBContainer::initialize_container()
{
    BL_PROFILE("amr-wind::ib::IBContainer::initialize_container");

    compute_local_coordinates();

    // Initialize global data arrays
    const int total_pts =
        std::accumulate(m_data.num_pts.begin(), m_data.num_pts.end(), 0);
    m_data.position.resize(total_pts);
    m_data.velocity.resize(total_pts);

    {
        const int nproc = amrex::ParallelDescriptor::NProcs();
        amrex::Vector<int> pts_per_proc(nproc, 0);
#ifdef AMREX_USE_MPI
        MPI_Allgather(
            &total_pts, 1, MPI_INT, pts_per_proc.data(), 1, MPI_INT,
            amrex::ParallelDescriptor::Communicator());
#else
        pts_per_proc[0] = total_pts;
#endif
        m_proc_offsets[0] = 0;
        for (int i = 1; i <= nproc; ++i) {
            m_proc_offsets[i] = m_proc_offsets[i - 1] + pts_per_proc[i - 1];
        }
        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, m_proc_offsets.begin(),
            m_proc_offsets.end(), m_proc_offsets_device.begin());
    }

    initialize_particles(total_pts);
}

void IBContainer::initialize_particles(const int total_pts)
{
    // Initialize particle container data structures.
    //
    // We assign all particles into the first available container within this
    // MPI rank and let redistribute take care of scattering the particles into
    // the respective MPI rank containing the cell enclosing this particle. The
    // actual redistribution happens after position vectors are updated in
    // update_positions.

    // query the particle ID from the container. We should always be starting
    // from 1.
    ParticleType::NextID(1u);
    const auto id_start = ParticleType::NextID();
    AMREX_ALWAYS_ASSERT(id_start == 1u);
    const int iproc = amrex::ParallelDescriptor::MyProc();

    // Flag indicating if a tile was found where all particles were deposited.
    bool assigned = false;
    const int nlevels = m_mesh.finestLevel() + 1;
    for (int lev = 0; ((lev < nlevels) && !assigned); ++lev) {
        for (auto mfi = MakeMFIter(lev); (mfi.isValid() && !assigned); ++mfi) {
            auto& ptile = GetParticles(
                lev)[std::make_pair(mfi.index(), mfi.LocalTileIndex())];
            AMREX_ASSERT(ptile.size() == 0);
            ptile.resize(total_pts);
            auto* pstruct = ptile.GetArrayOfStructs()().data();

            amrex::ParallelFor(
                total_pts, [=] AMREX_GPU_DEVICE(const int ip) noexcept {
                    auto& pp = pstruct[ip];

                    pp.id() = id_start + ip;
                    pp.cpu() = iproc;
                    pp.idata(0) = ip;
                });
            assigned = true;
        }
    }

    // Indicate that we have initialized the containers and remaining methods
    // are safe to use
    m_container_initialized = true;
    m_is_scattered = false;
}

void IBContainer::reset_container()
{
    const int nlevels = m_mesh.finestLevel() + 1;
    for (int lev = 0; lev < nlevels; ++lev) {
        for (ParIterType pti(*this, lev); pti.isValid(); ++pti) {
            const int np = pti.numParticles();
            auto* pstruct = pti.GetArrayOfStructs()().data();

            amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(const int ip) noexcept {
                auto& pp = pstruct[ip];
                pp.id() = -1;
            });
        }
    }
    Redistribute();

    const int total_pts = m_data.velocity.size();
    initialize_particles(total_pts);
}

/** Update position vectors of the particles within a container based on the
 *  data provided by immersed boundary instances.
 *
 *  This method assumes that the particles have been restored to their
 *  originating MPI ranks after any scattering. This happens during
 *  initialization (after AcutatorContainer::initialize_particles) and after
 *  IBContainer::sample_velocities.
 *
 *  After executing this method, the particles have been scattered across the
 *  domain such that each particle is contained within a tile that has the cell
 *  enclosing the particle.
 *
 */
void IBContainer::update_positions()
{
    BL_PROFILE("amr-wind::ib::IBContainer::update_positions");
    AMREX_ALWAYS_ASSERT(m_container_initialized && !m_is_scattered);

    const auto dpos = gpu::device_view(m_data.position);
    const auto dptr = dpos.data();
    const int nlevels = m_mesh.finestLevel() + 1;
    for (int lev = 0; lev < nlevels; ++lev) {
        for (ParIterType pti(*this, lev); pti.isValid(); ++pti) {
            const int np = pti.numParticles();
            auto* pstruct = pti.GetArrayOfStructs()().data();

            amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(const int ip) noexcept {
                auto& pp = pstruct[ip];
                const auto idx = pp.idata(0);

                auto& pvec = dptr[idx];
                for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                    pp.pos(n) = pvec[n];
                }
            });
        }
    }

    // Scatter particles to appropriate MPI ranks
    Redistribute();

    // Indicate that it is safe to sample velocities
    m_is_scattered = true;
}

/** Interpolate the velocity field using a trilinear interpolation
 *
 *  This method performs three tasks:
 *    - Sample the velocity field and interpolate onto particle location
 *    - Restore the particles back to their original MPI rank
 *    - Copy data from particles into the buffer used by immersed boundary
 * instances
 *
 *  After invocation if this method, the particles are back in their original
 *  rank and the position vectors can be updated safely.
 */
void IBContainer::sample_velocities(const Field& vel)
{
    BL_PROFILE("amr-wind::ib::IBContainer::sample_velocities");
    AMREX_ALWAYS_ASSERT(m_container_initialized && m_is_scattered);

    // Sample velocity field
    interpolate_velocities(vel);

    // Recall particles to the MPI ranks that contains their corresponding
    // turbines
    // Redistribute();

    // Populate the velocity buffer that all immersed boundary instances can
    // access
    populate_vel_buffer();

    // Indicate that the particles have been restored to their original MPI rank
    m_is_scattered = false;
}

/** Helper method for IBContainer::sample_velocities
 *
 *  Loops over the particle tiles and copies velocity data from particles to the
 *  velocity array.
 */
void IBContainer::populate_vel_buffer()
{
    BL_PROFILE("amr-wind::ib::IBContainer::populate_vel_buffer");
    amrex::Vector<vs::Vector> velh(m_proc_offsets.back());
    amrex::Gpu::DeviceVector<vs::Vector> vel(
        m_proc_offsets.back(), vs::Vector::zero());
    auto* varr = vel.data();
    auto* offsets = m_proc_offsets_device.data();
    const int nlevels = m_mesh.finestLevel() + 1;
    for (int lev = 0; lev < nlevels; ++lev) {
        for (ParIterType pti(*this, lev); pti.isValid(); ++pti) {
            const int np = pti.numParticles();
            auto* pstruct = pti.GetArrayOfStructs()().data();

            amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(const int ip) noexcept {
                auto& pp = pstruct[ip];
                const auto iproc = pp.cpu();
                const auto idx = offsets[iproc] + pp.idata(0);

                auto& vvel = varr[idx];
                for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                    vvel[n] = pp.rdata(n);
                }
            });
        }
    }

    amrex::Gpu::copy(
        amrex::Gpu::deviceToHost, vel.begin(), vel.end(), velh.begin());
#ifdef AMREX_USE_MPI
    MPI_Allreduce(
        MPI_IN_PLACE, &(velh[0][0]), m_proc_offsets.back() * AMREX_SPACEDIM,
        MPI_DOUBLE, MPI_SUM, amrex::ParallelDescriptor::Communicator());
#endif
    {
        auto& vel_arr = m_data.velocity;
        const int npts = vel_arr.size();
        const int ioff = m_proc_offsets[amrex::ParallelDescriptor::MyProc()];
        for (int i = 0; i < npts; ++i) {
            vel_arr[i] = velh[ioff + i];
        }
    }
}

/** Helper method for IBContainer::interpolate_velocities
 *
 */
void IBContainer::interpolate_velocities(const Field& vel)
{
    BL_PROFILE("amr-wind::ib::IBContainer::interpolate_velocities");

    auto* dptr = m_pos_device.data();
    const int nlevels = m_mesh.finestLevel() + 1;
    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& geom = m_mesh.Geom(lev);
        const auto dx = geom.CellSizeArray();
        const auto dxi = geom.InvCellSizeArray();
        const auto plo = geom.ProbLoArray();
        const amrex::Real dV = dx[0] * dx[1] * dx[2];

        for (ParIterType pti(*this, lev); pti.isValid(); ++pti) {
            const int np = pti.numParticles();
            auto* pstruct = pti.GetArrayOfStructs()().data();
            const auto varr = vel(lev).const_array(pti);

            amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(const int ip) noexcept {
                auto& pp = pstruct[ip];
                // Determine offsets within the containing cell
                const amrex::Real x =
                    (pp.pos(0) - plo[0] - 0.5 * dx[0]) * dxi[0];
                const amrex::Real y =
                    (pp.pos(1) - plo[1] - 0.5 * dx[1]) * dxi[1];
                const amrex::Real z =
                    (pp.pos(2) - plo[2] - 0.5 * dx[2]) * dxi[2];

                // Index of the low corner
                const int i = static_cast<int>(amrex::Math::floor(x));
                const int j = static_cast<int>(amrex::Math::floor(y));
                const int k = static_cast<int>(amrex::Math::floor(z));

                const int iproc = pp.cpu();
                amrex::Real interp_vel = 0.0;
                for (int ic = 0; ic < AMREX_SPACEDIM; ++ic) {
                    // Interpolating from five neighbouring points
                    // clang-format off
                    for (int ii = -2; ii <= 2; ++ii) {
                        for (int jj = -2; jj <= 2; ++jj) {
                            for (int kk = -2; kk <= 2; ++kk) {
                                const vs::Vector deltaX{
                                    plo[0] + (i + ii + 0.5) * dx[0] - pp.pos(0),
                                    plo[1] + (j + jj + 0.5) * dx[1] - pp.pos(1),
                                    plo[2] + (k + kk + 0.5) * dx[2] - pp.pos(2)};

                                interp_vel +=
                                    utils::dirac_delta(deltaX, {dx[0], dx[1], dx[2]}) *
                                    varr(i + ii, j + jj, k + kk, ic) * dV;
                            }
                        }
                    }

                    pp.rdata(ic) = interp_vel;
                    // clang-format on
                    // Reset position vectors so that the
                    // particles return back to the MPI ranks
                    // with the turbines upon redistribution
                    pp.pos(ic) = dptr[iproc][ic];
                }
            });
        }
    }
} // namespace ib

/** Determine position vector of a point within each MPI rank
 *
 *  Loops over the boxArray and determines the first patch that belongs
 * to the current MPI rank. Uses that to generate a position vector that
 * is known to exist in a given rank. Setting the position vectors of
 * the particles to this know location will recall particles belonging
 * to this rank during Redistribute.
 */
void IBContainer::compute_local_coordinates()
{
    BL_PROFILE("amr-wind::ib::IBContainer::compute_local_coordinates");
    const int nprocs = amrex::ParallelDescriptor::NProcs();
    const int iproc = amrex::ParallelDescriptor::MyProc();

    // Reset position vectors to zero (required for parallel reduce sum)
    m_proc_pos.assign(nprocs, vs::Vector::zero());

    // Flag indicating whether a point within the domain belonging to
    // this rank has been found
    bool assigned = false;
    const int nlevels = m_mesh.finestLevel() + 1;
    for (int lev = 0; ((lev < nlevels) && !assigned); ++lev) {
        const auto& ba = m_mesh.boxArray(lev);
        const auto& dm = m_mesh.DistributionMap(lev);

        const int nbx = dm.size();
        for (int i = 0; (i < nbx) && !assigned; ++i) {
            if (dm[i] != iproc) continue;

            const auto& geom = m_mesh.Geom(lev);
            const auto& dx = geom.CellSize();
            const auto& problo = geom.ProbLo();
            const auto& bx = ba[i];
            const int* lo = bx.loVect();

            auto& pvec = m_proc_pos[iproc];
            pvec.x() = problo[0] + (lo[0] + 0.5) * dx[0];
            pvec.y() = problo[1] + (lo[1] + 0.5) * dx[1];
            pvec.z() = problo[2] + (lo[2] + 0.5) * dx[2];

            // Indicate that we have found a point and it is safe to
            // exit the loop
            assigned = true;
        }
    }

    // Share position vectors with every process
    amrex::ParallelDescriptor::ReduceRealSum(
        &(m_proc_pos[0].x()), m_proc_pos.size() * vs::Vector::ncomp);
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, m_proc_pos.begin(), m_proc_pos.end(),
        m_pos_device.begin());
}

void IBContainer::post_regrid_actions() { compute_local_coordinates(); }

} // namespace ib
} // namespace amr_wind
