
#include "amr-wind/utilities/sampling/SamplingContainer.H"
#include "amr-wind/utilities/sampling/SamplerBase.H"
#include "amr-wind/core/Field.H"

namespace amr_wind {
namespace sampling {

namespace {

/** Interpolate a field to the sampling locations
 *
 *  \param np Number of particles in the container
 *  \param ic Component of the field to be interpolated
 *  \param pvec Vector containing particle info
 *  \param pavec Array information for the real component data
 *  \param farr Array of field data for this multifab
 *  \param dxi Inverse cell size array
 *  \param dx Cell size array
 *  \param offset Offsets for cell/node/face fields
 */
void sample_field(
    const int np,
    const int ic,
    SamplingContainer::ParticleVector& pvec,
    SamplingContainer::RealVector& pavec,
    const amrex::Array4<const amrex::Real>& farr,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& problo,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dxi,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dx,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& offset)
{
    BL_PROFILE("amr-wind::SamplingContainer::sample_impl");

    auto* pstruct = pvec.data();
    auto* parr = &(pavec[0]);

    amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(int ip) noexcept {
        auto& p = pstruct[ip];
        // Determine offsets within the containing cell
        const amrex::Real x =
            (p.pos(0) - problo[0] - offset[0] * dx[0]) * dxi[0];
        const amrex::Real y =
            (p.pos(1) - problo[1] - offset[1] * dx[1]) * dxi[1];
        const amrex::Real z =
            (p.pos(2) - problo[2] - offset[2] * dx[2]) * dxi[2];

        // Index of the low corner
        const int i = static_cast<int>(amrex::Math::floor(x));
        const int j = static_cast<int>(amrex::Math::floor(y));
        const int k = static_cast<int>(amrex::Math::floor(z));

        // Interpolation weights in each direction (linear basis)
        const amrex::Real wx_hi = (x - i);
        const amrex::Real wy_hi = (y - j);
        const amrex::Real wz_hi = (z - k);

        const amrex::Real wx_lo = 1.0 - wx_hi;
        const amrex::Real wy_lo = 1.0 - wy_hi;
        const amrex::Real wz_lo = 1.0 - wz_hi;

        parr[ip] = wx_lo * wy_lo * wz_lo * farr(i, j, k, ic) +
                   wx_lo * wy_lo * wz_hi * farr(i, j, k + 1, ic) +
                   wx_lo * wy_hi * wz_lo * farr(i, j + 1, k, ic) +
                   wx_lo * wy_hi * wz_hi * farr(i, j + 1, k + 1, ic) +
                   wx_hi * wy_lo * wz_lo * farr(i + 1, j, k, ic) +
                   wx_hi * wy_lo * wz_hi * farr(i + 1, j, k + 1, ic) +
                   wx_hi * wy_hi * wz_lo * farr(i + 1, j + 1, k, ic) +
                   wx_hi * wy_hi * wz_hi * farr(i + 1, j + 1, k + 1, ic);
    });
}

void sample_field(
    const int np,
    SamplingContainer::ParticleVector& pvec,
    SamplingContainer::RealVector& pavec,
    SamplingContainer::IntVector& piavec,
    const int nf,
    const amrex::Array4<const amrex::Real>& farr,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& problo,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dxi,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dx,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& offset)
{
    BL_PROFILE("amr-wind::SamplingContainer::sample_field_iso");
    const int ic = 0;

    auto* pstruct = pvec.data();
    auto* parr = &(pavec[0]);
    auto* piarr = &(piavec[0]);

    amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(int ip) noexcept {
        auto& p = pstruct[ip];
        // Check if current particle is concerned with current field
        if (p.idata(IIx::sid) != nf) return;
        // Check if current particle has no discernible valid range
        if (piarr[ip] == -3) return;

        // Determine offsets within the containing cell
        const amrex::Real x =
            (p.pos(0) - problo[0] - offset[0] * dx[0]) * dxi[0];
        const amrex::Real y =
            (p.pos(1) - problo[1] - offset[1] * dx[1]) * dxi[1];
        const amrex::Real z =
            (p.pos(2) - problo[2] - offset[2] * dx[2]) * dxi[2];

        // Index of the low corner
        const int i = static_cast<int>(amrex::Math::floor(x));
        const int j = static_cast<int>(amrex::Math::floor(y));
        const int k = static_cast<int>(amrex::Math::floor(z));

        // Interpolation weights in each direction (linear basis)
        const amrex::Real wx_hi = (x - i);
        const amrex::Real wy_hi = (y - j);
        const amrex::Real wz_hi = (z - k);

        const amrex::Real wx_lo = 1.0 - wx_hi;
        const amrex::Real wy_lo = 1.0 - wy_hi;
        const amrex::Real wz_lo = 1.0 - wz_hi;

        parr[ip] = wx_lo * wy_lo * wz_lo * farr(i, j, k, ic) +
                   wx_lo * wy_lo * wz_hi * farr(i, j, k + 1, ic) +
                   wx_lo * wy_hi * wz_lo * farr(i, j + 1, k, ic) +
                   wx_lo * wy_hi * wz_hi * farr(i, j + 1, k + 1, ic) +
                   wx_hi * wy_lo * wz_lo * farr(i + 1, j, k, ic) +
                   wx_hi * wy_lo * wz_hi * farr(i + 1, j, k + 1, ic) +
                   wx_hi * wy_hi * wz_lo * farr(i + 1, j + 1, k, ic) +
                   wx_hi * wy_hi * wz_hi * farr(i + 1, j + 1, k + 1, ic);
    });
}

void init_bounds(
    const int np,
    SamplingContainer::ParticleVector& pvec,
    amrex::Array<amrex::Real*, AMREX_SPACEDIM>& plvec,
    amrex::Array<amrex::Real*, AMREX_SPACEDIM>& prvec,
    amrex::Array<amrex::Real*, AMREX_SPACEDIM>& povec,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& problo,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& probhi,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dx)
{
    BL_PROFILE("amr-wind::SamplingContainer::init_bounds");

    auto* pstruct = pvec.data();

    amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(int ip) noexcept {
        auto& p = pstruct[ip];
        amrex::Real dxmag = dx[0];
        amrex::Real scale = 0.0;
        for (int n = 0; n < AMREX_SPACEDIM; ++n) {
            // Store current location as left bound
            (plvec[n])[ip] = p.pos(n);
            // Calculate the magnitude of dx
            dxmag = std::min(dx[n], dxmag);
            // Calculate the max dimension of domain to use with epsilon
            scale = std::max(scale, probhi[n] - problo[n]);
        }
        // Step along orientation vector until bounds are exceeded
        bool flag = false;
        amrex::Real dist = 0.0;
        int nn = 0;
        while (!flag) {
            for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                // Start at left bound
                if (nn == 0) (prvec[n])[ip] = (plvec[n])[ip];
                // Increment with dx size
                (prvec[n])[ip] += (povec[n])[ip] * dxmag;
                // Check bounds
                if ((prvec[n])[ip] >= probhi[n] || (prvec[n])[ip] < problo[n]) {
                    // Flag to indicate finished
                    flag = true;
                    // Distance that bounds have been exceeded, along vector
                    dist = std::max(
                        dist, std::max(
                                  (prvec[n])[ip] - probhi[n],
                                  problo[n] - (prvec[n])[ip]) /
                                  (povec[n])[ip]);
                }
            }
            // After bounds have been exceeded, remove latest contribution
            // to stay in domain
            if (flag) {
                for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                    (prvec[n])[ip] -=
                        (povec[n])[ip] *
                        (dist +
                         scale * std::numeric_limits<amrex::Real>::epsilon());
                }
            }
            ++nn;
        }
    });
}

void reset_iflag(const int np, SamplingContainer::IntVector& piarr)
{
    BL_PROFILE("amr-wind::SamplingContainer::reset_iflag");

    auto* iflag = &(piarr[0]);
    amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(int ip) noexcept {
        // Reset integer values so that bounds can be checked
        iflag[ip] = 0;
    });
}

void iso_fields(
    const int lev,
    const int np,
    const amrex::Vector<Field*> fields,
    SamplingContainer::ParIterType& pti,
    SamplingContainer::ParticleVector& pvec,
    SamplingContainer::RealVector& parr,
    SamplingContainer::IntVector& piarr,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& plo,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dxi,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dx)
{
    BL_PROFILE("amr-wind::SamplingContainer::iso_fields");
    int fidx = 0;
    for (const auto* fld : fields) {
        const auto farr = (*fld)(lev).const_array(pti);

        switch (fld->field_location()) {
        case FieldLoc::NODE: {
            amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> offset{
                {0.0, 0.0, 0.0}};
            sample_field(
                np, pvec, parr, piarr, fidx, farr, plo, dxi, dx, offset);
            break;
        }

        case FieldLoc::CELL: {
            amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> offset{
                {0.5, 0.5, 0.5}};
            sample_field(
                np, pvec, parr, piarr, fidx, farr, plo, dxi, dx, offset);
            break;
        }

        case FieldLoc::XFACE: {
            amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> offset{
                {0.0, 0.5, 0.5}};
            sample_field(
                np, pvec, parr, piarr, fidx, farr, plo, dxi, dx, offset);
            break;
        }

        case FieldLoc::YFACE: {
            amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> offset{
                {0.5, 0.0, 0.5}};
            sample_field(
                np, pvec, parr, piarr, fidx, farr, plo, dxi, dx, offset);
            break;
        }

        case FieldLoc::ZFACE: {
            amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> offset{
                {0.5, 0.5, 0.0}};
            sample_field(
                np, pvec, parr, piarr, fidx, farr, plo, dxi, dx, offset);
            break;
        }
        }
        // Increment field counter
        ++fidx;
    }
}

void update_position(
    const int np,
    SamplingContainer::ParticleVector& pvec,
    amrex::Array<amrex::Real*, AMREX_SPACEDIM>& posvec)
{
    BL_PROFILE("amr-wind::SamplingContainer::update_position");

    auto* pstruct = pvec.data();

    amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(int ip) noexcept {
        auto& p = pstruct[ip];
        for (int n = 0; n < AMREX_SPACEDIM; ++n) {
            // Copy
            p.pos(n) = (posvec[n])[ip];
        }
    });
}

void update_position(
    const int np,
    SamplingContainer::ParticleVector& pvec,
    amrex::Array<amrex::Real*, AMREX_SPACEDIM>& posvecl,
    amrex::Array<amrex::Real*, AMREX_SPACEDIM>& posvecr)
{
    BL_PROFILE("amr-wind::SamplingContainer::update_position_mid");

    auto* pstruct = pvec.data();

    amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(int ip) noexcept {
        auto& p = pstruct[ip];
        for (int n = 0; n < AMREX_SPACEDIM; ++n) {
            // Copy
            p.pos(n) = 0.5 * ((posvecl[n])[ip] + (posvecr[n])[ip]);
        }
    });
}

// Check sign change, get new position, update flag
void pre_bisect_work(
    const int np,
    SamplingContainer::ParticleVector& pvec,
    SamplingContainer::RealVector& ptarr,
    SamplingContainer::RealVector& plarr,
    SamplingContainer::RealVector& prarr,
    amrex::Array<amrex::Real*, AMREX_SPACEDIM>& plvec,
    amrex::Array<amrex::Real*, AMREX_SPACEDIM>& prvec,
    amrex::Array<amrex::Real*, AMREX_SPACEDIM>& povec,
    SamplingContainer::IntVector& piarr,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& problo,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& probhi,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dx,
    const int ct,
    bool& flag)
{
    BL_PROFILE("amr-wind::SamplingContainer::pre_bisect_work");

    // Pointers
    auto* pstruct = pvec.data();
    auto* target = &(ptarr[0]);
    auto* lval = &(plarr[0]);
    auto* rval = &(prarr[0]);
    auto* iflag = &(piarr[0]);

    int loopsum = 0;
    amrex::ParallelFor(np, [=, &loopsum] AMREX_GPU_DEVICE(int ip) noexcept {
        int not_finished = 1;
        // Check for flags indicating done
        if (iflag[ip] == -3 || iflag[ip] == 1) {
            not_finished = 0;
            loopsum += not_finished;
            return;
        }
        // Check for sign change
        if ((lval[ip] - target[ip]) * (rval[ip] - target[ip]) <= 0) {
            // Sign change is present, signify that work is done
            iflag[ip] = 1;
            not_finished = 0;
            loopsum += not_finished;
            return;
        }

        // Take a dx step along orientation vector
        auto& p = pstruct[ip];
        amrex::Real dxmag = dx[0];
        for (int n = 0; n < AMREX_SPACEDIM; ++n) {
            // Calculate the magnitude of dx
            dxmag = std::min(dx[n], dxmag);
        }

        // Left or right depends on modulus of loop counter
        switch (ct % 2) {
        case 0: {
            // Check if bounds have been exceeded in this direction
            if (iflag[ip] == -1) {
                loopsum += not_finished;
                return;
            }
            for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                (plvec[n])[ip] -= (povec[n])[ip] * dxmag;
                if ((plvec[n])[ip] > probhi[n] || (plvec[n])[ip] < problo[n]) {
                    // If domain has been exceeded in other direction
                    // this sampling point is hopeless (-3)
                    if (iflag[ip] == -2) {
                        iflag[ip] = -3;
                        not_finished = 0;
                    } else {
                        // Flag to indicate bounds exceeded
                        iflag[ip] = -1;
                    }
                }
            }
            if (iflag[ip] == -1 || iflag[ip] == -3) {
                // Make sure point remains in domain by removing last increment
                for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                    (plvec[n])[ip] += (povec[n])[ip] * dxmag;
                }
            }
            // Position of particle assigned to left position
            for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                p.pos(n) = (plvec[n])[ip];
            }
        }
        case 1: {
            // Check if bounds have been exceeded in this direction
            if (iflag[ip] == -2) {
                loopsum += not_finished;
                return;
            }
            for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                (prvec[n])[ip] += (povec[n])[ip] * dxmag;
                if ((prvec[n])[ip] > probhi[n] || (prvec[n])[ip] < problo[n]) {
                    // If domain has been exceeded in other direction
                    // this sampling point is hopeless (-3)
                    if (iflag[ip] == -1) {
                        iflag[ip] = -3;
                        not_finished = 0;
                    } else {
                        // Flag to indicate bounds exceeded
                        iflag[ip] = -2;
                    }
                }
            }
            if (iflag[ip] == -2 || iflag[ip] == -3) {
                // Make sure point remains in domain by removing last increment
                for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                    (prvec[n])[ip] -= (povec[n])[ip] * dxmag;
                }
            }
            // Position of particle assigned to right position
            for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                p.pos(n) = (prvec[n])[ip];
            }
        }
        }
        loopsum += not_finished;
        return;
    });
    if (loopsum == 0) flag = true;
}

// Update flag, send to new middle position
void bisect_work(
    const int np,
    SamplingContainer::ParticleVector& pvec,
    SamplingContainer::RealVector& pcarr,
    SamplingContainer::RealVector& ptarr,
    SamplingContainer::RealVector& plarr,
    SamplingContainer::RealVector& prarr,
    amrex::Array<amrex::Real*, AMREX_SPACEDIM>& plvec,
    amrex::Array<amrex::Real*, AMREX_SPACEDIM>& prvec,
    SamplingContainer::IntVector& piarr,
    const amrex::Real& tol,
    bool& flag)
{
    BL_PROFILE("amr-wind::SamplingContainer::bisect_work");

    // Pointers
    auto* pstruct = pvec.data();
    auto* current = &(pcarr[0]);
    auto* target = &(ptarr[0]);
    auto* lval = &(plarr[0]);
    auto* rval = &(prarr[0]);
    auto* iflag = &(piarr[0]);

    int loopsum = 0;
    amrex::ParallelFor(np, [=, &loopsum] AMREX_GPU_DEVICE(int ip) noexcept {
        int not_finished = 1;
        // Check for flags indicating no solution for current probe
        if (iflag[ip] == -3) {
            not_finished = 0;
            loopsum += not_finished;
            return;
        }
        // Check if tolerance is satisfied for target value
        if (std::abs(current[ip] - target[ip]) < tol) {
            not_finished = 0;
            loopsum += not_finished;
            return;
        }
        auto& p = pstruct[ip];
        // Use sign change to determine new middle, reassign bounds
        if ((current[ip] - target[ip]) * (lval[ip] - target[ip]) <= 0) {
            // Sign change is in left interval
            // (or could be that left or middle is at target)
            rval[ip] = current[ip];
            for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                // Right position moves to current
                (prvec[n])[ip] = p.pos(n);
                // Particle is moved to new middle
                p.pos(n) = 0.5 * ((plvec[n])[ip] + (prvec[n])[ip]);
            }
        } else {
            // Sign change is in right interval
            lval[ip] = current[ip];
            for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                // Left position moves to current
                (plvec[n])[ip] = p.pos(n);
                // Particle is moved to new middle
                p.pos(n) = 0.5 * ((plvec[n])[ip] + (prvec[n])[ip]);
            }
        }
        loopsum += not_finished;
        return;
    });
    if (loopsum == 0) flag = true;
}
} // namespace

void SamplingContainer::setup_container(
    const int num_real_components, const int num_int_components)
{
    BL_PROFILE("amr-wind::SamplingContainer::setup");
    const bool communicate_comp = true;
    for (int i = 0; i < num_real_components; ++i) AddRealComp(communicate_comp);
    for (int i = 0; i < num_int_components; ++i) AddIntComp(communicate_comp);

    const int nlevels = m_mesh.finestLevel() + 1;
    for (int lev = 0; lev < nlevels; ++lev) {
        for (amrex::MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi) {
            DefineAndReturnParticleTile(lev, mfi.index(), mfi.LocalTileIndex());
        }
    }
}

void SamplingContainer::initialize_particles(
    const amrex::Vector<std::unique_ptr<SamplerBase>>& samplers)
{
    BL_PROFILE("amr-wind::SamplingContainer::initialize");

    // We will assign all particles to the first box in level 0 and let
    // redistribute scatter it to the appropriate rank and box.
    const int lev = 0;
    const int iproc = amrex::ParallelDescriptor::MyProc();
    const int owner = ParticleDistributionMap(lev)[0];

    // Let only the MPI rank owning the first box do the work
    if (owner != iproc) return;

    int num_particles = 0;
    for (auto& probes : samplers) num_particles += probes->num_points();
    m_total_particles = num_particles;

    const int grid_id = 0;
    const int tile_id = 0;
    // Setup has already processed runtime array information, safe to do
    // GetParticles here.
    auto& ptile = GetParticles(lev)[std::make_pair(grid_id, tile_id)];
    ptile.resize(num_particles);

    int pidx = 0;
    const int nextid = ParticleType::NextID();
    auto* pstruct = ptile.GetArrayOfStructs()().data();
    SamplerBase::SampleLocType locs;
    for (auto& probe : samplers) {
        probe->sampling_locations(locs);
        const int npts = locs.size();
        const auto probe_id = probe->id();
        amrex::Gpu::DeviceVector<amrex::Real> dlocs(npts * AMREX_SPACEDIM);
        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, locs.begin(), locs.end(), dlocs.begin());
        const auto* dpos = dlocs.data();

        amrex::ParallelFor(npts, [=] AMREX_GPU_DEVICE(const int ip) noexcept {
            const auto uid = pidx + ip;
            auto& pp = pstruct[uid];
            pp.id() = nextid + uid;
            pp.cpu() = iproc;

            for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                pp.pos(n) = dpos[ip * AMREX_SPACEDIM + n];
            }
            pp.idata(IIx::uid) = uid;
            pp.idata(IIx::sid) = probe_id;
            pp.idata(IIx::nid) = ip;
        });
        amrex::Gpu::streamSynchronize();
        pidx += npts;
    }

    AMREX_ALWAYS_ASSERT(pidx == num_particles);
}

void SamplingContainer::initialize_particles(
    const amrex::Vector<std::unique_ptr<SamplerBase>>& samplers,
    const amrex::Vector<amrex::Real>& field_vals)
{
    BL_PROFILE("amr-wind::SamplingContainer::initialize_iso");

    // We will assign all particles to the first box in level 0 and let
    // redistribute scatter it to the appropriate rank and box.
    const int lev = 0;
    const int iproc = amrex::ParallelDescriptor::MyProc();
    const int owner = ParticleDistributionMap(lev)[0];

    // Let only the MPI rank owning the first box do the work
    if (owner != iproc) return;

    int num_particles = 0;
    for (auto& probes : samplers) num_particles += probes->num_points();
    m_total_particles = num_particles;

    const int grid_id = 0;
    const int tile_id = 0;
    // Setup has already processed runtime array information, safe to do
    // GetParticles here.
    auto& ptile = GetParticles(lev)[std::make_pair(grid_id, tile_id)];
    ptile.resize(num_particles);

    int pidx = 0;
    const int nextid = ParticleType::NextID();
    auto* pstruct = ptile.GetArrayOfStructs()().data();
    // Get data access for target values and orientations
    auto* pvalues = &((ptile.GetStructOfArrays().GetRealData(1))[0]);
    amrex::Array<decltype(pvalues), AMREX_SPACEDIM> porients;
    // Data access for lone integer array
    auto* pints = &((ptile.GetStructOfArrays().GetIntData(0))[0]);
    // First index where orientation is stored
    int roffset = 4 + 2 * AMREX_SPACEDIM;
    for (int n = roffset; n < roffset + AMREX_SPACEDIM; ++n) {
        porients[n - roffset] =
            &((ptile.GetStructOfArrays().GetRealData(n))[0]);
    }
    SamplerBase::SampleLocType locs, oris;
    for (auto& probe : samplers) {
        probe->sampling_locations(locs);
        probe->sampling_orientations(oris);
        const int npts = locs.size();
        const auto probe_id = probe->id();
        amrex::Gpu::DeviceVector<amrex::Real> dlocs(npts * AMREX_SPACEDIM);
        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, locs.begin(), locs.end(), dlocs.begin());
        const auto* dpos = dlocs.data();
        amrex::Gpu::DeviceVector<amrex::Real> doris(npts * AMREX_SPACEDIM);
        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, oris.begin(), oris.end(), doris.begin());
        const auto* dors = doris.data();

        amrex::ParallelFor(npts, [=] AMREX_GPU_DEVICE(const int ip) noexcept {
            const auto uid = pidx + ip;
            auto& pp = pstruct[uid];
            pp.id() = nextid + uid;
            pp.cpu() = iproc;

            for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                pp.pos(n) = dpos[ip * AMREX_SPACEDIM + n];
                // Data members - vector
                (porients[n])[uid] = dors[ip * AMREX_SPACEDIM + n];
            }
            pp.idata(IIx::uid) = uid;
            pp.idata(IIx::sid) = probe_id;
            pp.idata(IIx::nid) = ip;
            // Data members - scalar
            pvalues[uid] = field_vals[probe_id];
            pints[uid] = 0;
        });
        amrex::Gpu::streamSynchronize();
        pidx += npts;
    }

    AMREX_ALWAYS_ASSERT(pidx == num_particles);
}

void SamplingContainer::iso_bounds_pos()
{
    BL_PROFILE("amr-wind::SamplingContainer::iso_bounds_pos");

    const int nlevels = m_mesh.finestLevel() + 1;

    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& geom = m_mesh.Geom(lev);
        const auto dx = geom.CellSizeArray();
        const auto plo = geom.ProbLoArray();
        const auto phi = geom.ProbHiArray();

        for (ParIterType pti(*this, lev); pti.isValid(); ++pti) {
            const int np = pti.numParticles();
            auto& pvec = pti.GetArrayOfStructs()();

            // Set up real array components
            amrex::Array<amrex::Real*, AMREX_SPACEDIM> pllocs;
            amrex::Array<amrex::Real*, AMREX_SPACEDIM> prlocs;
            amrex::Array<amrex::Real*, AMREX_SPACEDIM> porients;

            // First index where left location is stored, among real components
            int roffset = 4;
            for (int n = roffset; n < roffset + AMREX_SPACEDIM; ++n) {
                pllocs[n - roffset] =
                    &(pti.GetStructOfArrays().GetRealData(n))[0];
                int nn = n + AMREX_SPACEDIM;
                prlocs[n - roffset] =
                    &(pti.GetStructOfArrays().GetRealData(nn))[0];
                nn += AMREX_SPACEDIM;
                porients[n - roffset] =
                    &(pti.GetStructOfArrays().GetRealData(nn))[0];
            }

            // Get left and right positions
            init_bounds(np, pvec, pllocs, prlocs, porients, plo, phi, dx);
        }
    }
}

void SamplingContainer::iso_bounds_val(const amrex::Vector<Field*> fields)
{
    BL_PROFILE("amr-wind::SamplingContainer::iso_bounds_val");

    const int nlevels = m_mesh.finestLevel() + 1;

    // Loop for sending particles to current left location
    for (int lev = 0; lev < nlevels; ++lev) {

        for (ParIterType pti(*this, lev); pti.isValid(); ++pti) {
            const int np = pti.numParticles();
            auto& pvec = pti.GetArrayOfStructs()();

            // Get left location data
            amrex::Array<amrex::Real*, AMREX_SPACEDIM> pllocs;
            // First index where left location is stored, among real components
            int roffset = 4;
            for (int n = roffset; n < roffset + AMREX_SPACEDIM; ++n) {
                pllocs[n - roffset] =
                    &(pti.GetStructOfArrays().GetRealData(n))[0];
            }
            // Update location of particles to current left
            update_position(np, pvec, pllocs);
        }
    }
    // Redistribute, since position has changed
    this->Redistribute();

    // Loop to measure left value and send to right
    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& geom = m_mesh.Geom(lev);
        const auto dx = geom.CellSizeArray();
        const auto dxi = geom.InvCellSizeArray();
        const auto plo = geom.ProbLoArray();

        for (ParIterType pti(*this, lev); pti.isValid(); ++pti) {
            const int np = pti.numParticles();
            auto& pvec = pti.GetArrayOfStructs()();

            // Set up real and int array components
            auto& plvals = pti.GetStructOfArrays().GetRealData(2);
            auto& pints = pti.GetStructOfArrays().GetIntData(0);
            amrex::Array<amrex::Real*, AMREX_SPACEDIM> prlocs;

            // First index where right location is stored, among real components
            int roffset = 4 + AMREX_SPACEDIM;
            for (int n = roffset; n < roffset + AMREX_SPACEDIM; ++n) {
                prlocs[n - roffset] =
                    &(pti.GetStructOfArrays().GetRealData(n))[0];
            }

            // Reset integer flag
            reset_iflag(np, pints);
            // Get current value and set as left value
            iso_fields(lev, np, fields, pti, pvec, plvals, pints, plo, dxi, dx);
            // Update location of particles to current right value
            update_position(np, pvec, prlocs);
        }
    }
    // Redistribute, since position has changed
    this->Redistribute();

    // Loop to measure right value
    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& geom = m_mesh.Geom(lev);
        const auto dx = geom.CellSizeArray();
        const auto dxi = geom.InvCellSizeArray();
        const auto plo = geom.ProbLoArray();

        for (ParIterType pti(*this, lev); pti.isValid(); ++pti) {
            const int np = pti.numParticles();
            auto& pvec = pti.GetArrayOfStructs()();

            // Set up real components
            auto& prvals = pti.GetStructOfArrays().GetRealData(3);
            // Integer flag (should be all 0's at this point)
            auto& pints = pti.GetStructOfArrays().GetIntData(0);

            // Get current value and set as right value
            iso_fields(lev, np, fields, pti, pvec, prvals, pints, plo, dxi, dx);
        }
    }
}

void SamplingContainer::interpolate_fields(const amrex::Vector<Field*> fields)
{
    BL_PROFILE("amr-wind::SamplingContainer::interpolate");

    const int nlevels = m_mesh.finestLevel() + 1;

    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& geom = m_mesh.Geom(lev);
        const auto dx = geom.CellSizeArray();
        const auto dxi = geom.InvCellSizeArray();
        const auto plo = geom.ProbLoArray();

        for (ParIterType pti(*this, lev); pti.isValid(); ++pti) {
            const int np = pti.numParticles();
            auto& pvec = pti.GetArrayOfStructs()();

            int fidx = 0;
            for (const auto* fld : fields) {
                const auto farr = (*fld)(lev).const_array(pti);
                for (int ic = 0; ic < fld->num_comp(); ++ic) {
                    auto& parr = pti.GetStructOfArrays().GetRealData(fidx++);

                    switch (fld->field_location()) {
                    case FieldLoc::NODE: {
                        amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> offset{
                            {0.0, 0.0, 0.0}};
                        sample_field(
                            np, ic, pvec, parr, farr, plo, dxi, dx, offset);
                        break;
                    }

                    case FieldLoc::CELL: {
                        amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> offset{
                            {0.5, 0.5, 0.5}};
                        sample_field(
                            np, ic, pvec, parr, farr, plo, dxi, dx, offset);
                        break;
                    }

                    case FieldLoc::XFACE: {
                        amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> offset{
                            {0.0, 0.5, 0.5}};
                        sample_field(
                            np, ic, pvec, parr, farr, plo, dxi, dx, offset);
                        break;
                    }

                    case FieldLoc::YFACE: {
                        amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> offset{
                            {0.5, 0.0, 0.5}};
                        sample_field(
                            np, ic, pvec, parr, farr, plo, dxi, dx, offset);
                        break;
                    }

                    case FieldLoc::ZFACE: {
                        amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> offset{
                            {0.5, 0.5, 0.0}};
                        sample_field(
                            np, ic, pvec, parr, farr, plo, dxi, dx, offset);
                        break;
                    }
                    }
                }
            }
        }
    }
}

void SamplingContainer::iso_relocate(const amrex::Vector<Field*> fields)
{
    BL_PROFILE("amr-wind::SamplingContainer::iso_relocate");

    const int nlevels = m_mesh.finestLevel() + 1;
    amrex::Array<amrex::Real*, AMREX_SPACEDIM> pllocs;
    amrex::Array<amrex::Real*, AMREX_SPACEDIM> prlocs;
    amrex::Array<amrex::Real*, AMREX_SPACEDIM> porients;

    //! Pre-bisection loop - for adjusting bounds
    bool flag = false;
    int ct = 0;
    while (!flag) {
        flag = true;
        for (int lev = 0; lev < nlevels; ++lev) {
            const auto& geom = m_mesh.Geom(lev);
            const auto dx = geom.CellSizeArray();
            const auto plo = geom.ProbLoArray();
            const auto phi = geom.ProbHiArray();

            for (ParIterType pti(*this, lev); pti.isValid(); ++pti) {
                const int np = pti.numParticles();
                auto& pvec = pti.GetArrayOfStructs()();
                auto& ptvals = pti.GetStructOfArrays().GetRealData(1);
                auto& plvals = pti.GetStructOfArrays().GetRealData(2);
                auto& prvals = pti.GetStructOfArrays().GetRealData(3);
                auto& pints = pti.GetStructOfArrays().GetIntData(0);

                int roffset = 4;
                for (int n = roffset; n < roffset + AMREX_SPACEDIM; ++n) {
                    pllocs[n - roffset] =
                        &(pti.GetStructOfArrays().GetRealData(n))[0];
                    int nn = n + AMREX_SPACEDIM;
                    prlocs[n - roffset] =
                        &(pti.GetStructOfArrays().GetRealData(nn))[0];
                    nn += AMREX_SPACEDIM;
                    porients[n - roffset] =
                        &(pti.GetStructOfArrays().GetRealData(nn))[0];
                }
                bool out = false;
                // Check sign change, get new position, update flag
                pre_bisect_work(
                    np, pvec, ptvals, plvals, prvals, pllocs, prlocs, porients,
                    pints, plo, phi, dx, ct, out);
                if (!out) flag = false;
            }
            ++ct;
        }
        // With particles having moved, go to new position
        this->Redistribute();
        // Get value at new position
        for (int lev = 0; lev < nlevels; ++lev) {
            const auto& geom = m_mesh.Geom(lev);
            const auto dx = geom.CellSizeArray();
            const auto dxi = geom.InvCellSizeArray();
            const auto plo = geom.ProbLoArray();

            for (ParIterType pti(*this, lev); pti.isValid(); ++pti) {
                const int np = pti.numParticles();
                auto& pvec = pti.GetArrayOfStructs()();
                auto& pints = pti.GetStructOfArrays().GetIntData(0);

                // Use left or right based on modulus
                auto& parr = pti.GetStructOfArrays().GetRealData(2 + ct % 2);
                // Get value at current location
                iso_fields(
                    lev, np, fields, pti, pvec, parr, pints, plo, dxi, dx);
            }
        }
    }

    //! In prep for main loop
    for (int lev = 0; lev < nlevels; ++lev) {
        for (ParIterType pti(*this, lev); pti.isValid(); ++pti) {
            const int np = pti.numParticles();
            auto& pvec = pti.GetArrayOfStructs()();

            int roffset = 4;
            for (int n = roffset; n < roffset + AMREX_SPACEDIM; ++n) {
                pllocs[n - roffset] =
                    &(pti.GetStructOfArrays().GetRealData(n))[0];
                int nn = n + AMREX_SPACEDIM;
                prlocs[n - roffset] =
                    &(pti.GetStructOfArrays().GetRealData(nn))[0];
            }

            // Send particle to center position
            update_position(np, pvec, pllocs, prlocs);
        }
    }
    // With particles having moved, go to new position
    this->Redistribute();

    //! Bisection loop
    flag = false;
    while (!flag) {
        flag = true;
        for (int lev = 0; lev < nlevels; ++lev) {
            const auto& geom = m_mesh.Geom(lev);
            const auto dx = geom.CellSizeArray();
            const auto dxi = geom.InvCellSizeArray();
            const auto plo = geom.ProbLoArray();

            for (ParIterType pti(*this, lev); pti.isValid(); ++pti) {
                const int np = pti.numParticles();
                auto& pvec = pti.GetArrayOfStructs()();
                auto& pcvals = pti.GetStructOfArrays().GetRealData(0);
                auto& ptvals = pti.GetStructOfArrays().GetRealData(1);
                auto& plvals = pti.GetStructOfArrays().GetRealData(2);
                auto& prvals = pti.GetStructOfArrays().GetRealData(3);
                auto& pints = pti.GetStructOfArrays().GetIntData(0);

                int roffset = 4;
                for (int n = roffset; n < roffset + AMREX_SPACEDIM; ++n) {
                    pllocs[n - roffset] =
                        &(pti.GetStructOfArrays().GetRealData(n))[0];
                    int nn = n + AMREX_SPACEDIM;
                    prlocs[n - roffset] =
                        &(pti.GetStructOfArrays().GetRealData(nn))[0];
                }
                bool out = false;
                // Get middle value
                iso_fields(
                    lev, np, fields, pti, pvec, pcvals, pints, plo, dxi, dx);
                // Update flag, send to new middle position
                bisect_work(
                    np, pvec, pcvals, ptvals, plvals, prvals, pllocs, prlocs,
                    pints, m_tol, out);
                if (!out) flag = false;
            }
        }
        // With particles having moved, go to new position
        this->Redistribute();
    }
}

void SamplingContainer::populate_buffer(std::vector<double>& buf)
{
    BL_PROFILE("amr-wind::SamplingContainer::populate_buffer");

    amrex::Gpu::DeviceVector<double> dbuf(buf.size(), 0.0);
    auto* dbuf_ptr = dbuf.data();
    const int nlevels = m_mesh.finestLevel() + 1;
    for (int lev = 0; lev < nlevels; ++lev) {
        for (int fid = 0; fid < NumRuntimeRealComps(); ++fid) {
            const int offset = fid * num_sampling_particles();
            for (ParIterType pti(*this, lev); pti.isValid(); ++pti) {
                const int np = pti.numParticles();
                auto* pstruct = pti.GetArrayOfStructs()().data();
                auto* parr = &pti.GetStructOfArrays().GetRealData(fid)[0];

                amrex::ParallelFor(
                    np, [=] AMREX_GPU_DEVICE(const int ip) noexcept {
                        auto& pp = pstruct[ip];
                        const int pidx = pp.idata(IIx::uid);
                        const int ii = offset + pidx;
                        dbuf_ptr[ii] = parr[ip];
                    });
            }
        }
    }

    amrex::Gpu::copy(
        amrex::Gpu::deviceToHost, dbuf.begin(), dbuf.end(), buf.begin());
    amrex::ParallelDescriptor::ReduceRealSum(
        buf.data(), buf.size(), amrex::ParallelDescriptor::IOProcessorNumber());
}

void SamplingContainer::populate_buffer(std::vector<int>& buf)
{
    BL_PROFILE("amr-wind::SamplingContainer::populate_buffer_int");

    amrex::Gpu::DeviceVector<int> dbuf(buf.size(), 0.0);
    auto* dbuf_ptr = dbuf.data();
    const int nlevels = m_mesh.finestLevel() + 1;
    for (int lev = 0; lev < nlevels; ++lev) {
        for (int fid = 0; fid < NumRuntimeIntComps(); ++fid) {
            const int offset = fid * num_sampling_particles();
            for (ParIterType pti(*this, lev); pti.isValid(); ++pti) {
                const int np = pti.numParticles();
                auto* pstruct = pti.GetArrayOfStructs()().data();
                auto* parr = &pti.GetStructOfArrays().GetIntData(fid)[0];

                amrex::ParallelFor(
                    np, [=] AMREX_GPU_DEVICE(const int ip) noexcept {
                        auto& pp = pstruct[ip];
                        const int pidx = pp.idata(IIx::uid);
                        const int ii = offset + pidx;
                        dbuf_ptr[ii] = parr[ip];
                    });
            }
        }
    }

    amrex::Gpu::copy(
        amrex::Gpu::deviceToHost, dbuf.begin(), dbuf.end(), buf.begin());
    amrex::ParallelDescriptor::ReduceIntSum(
        buf.data(), buf.size(), amrex::ParallelDescriptor::IOProcessorNumber());
}

} // namespace sampling
} // namespace amr_wind
