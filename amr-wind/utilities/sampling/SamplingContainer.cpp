
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
        const int i = std::floor(x);
        const int j = std::floor(y);
        const int k = std::floor(z);

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
    auto* pstruct = ptile.GetArrayOfStructs()().data();
    SamplerBase::SampleLocType locs;
    for (auto& probe : samplers) {
        probe->sampling_locations(locs);

        const int npts = locs.size();
        for (int ip = 0; ip < npts; ++ip) {
            auto& pp = pstruct[pidx];
            pp.id() = ParticleType::NextID();
            pp.cpu() = iproc;

            pp.pos(0) = locs[ip][0];
            pp.pos(1) = locs[ip][1];
            pp.pos(2) = locs[ip][2];

            pp.idata(IIx::uid) = pidx;
            pp.idata(IIx::sid) = probe->id();
            pp.idata(IIx::nid) = ip;

            ++pidx;
        }
    }

    AMREX_ALWAYS_ASSERT(pidx == num_particles);
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

void SamplingContainer::populate_buffer(std::vector<double>& buf)
{
    BL_PROFILE("amr-wind::SamplingContainer::populate_buffer");

    const int nlevels = m_mesh.finestLevel() + 1;
    for (int lev=0; lev < nlevels; ++lev) {
        for (int fid=0; fid < NumRuntimeRealComps(); ++fid) {
            const int offset = fid * num_sampling_particles();
            for (ParIterType pti(*this, lev); pti.isValid(); ++pti) {
                const int np = pti.numParticles();
                auto* pstruct = pti.GetArrayOfStructs()().data();
                auto* parr = &pti.GetStructOfArrays().GetRealData(fid)[0];

                for (int ip=0; ip < np; ++ip) {
                    auto& pp = pstruct[ip];
                    const int pidx = pp.idata(IIx::uid);
                    const int ii = offset + pidx;

                    buf[ii] = parr[ip];
                }
            }
        }
    }
    amrex::ParallelDescriptor::ReduceRealSum(
        buf.data(), buf.size(), amrex::ParallelDescriptor::IOProcessorNumber());
}

} // namespace sampling
} // namespace amr_wind
