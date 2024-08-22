
#include "amr-wind/utilities/sampling/SamplingContainer.H"
#include "amr-wind/utilities/sampling/SamplerBase.H"
#include "amr-wind/core/Field.H"

namespace amr_wind::sampling {

void SamplingContainer::setup_container(
    const int num_real_components, const int num_int_components)
{
    BL_PROFILE("amr-wind::SamplingContainer::setup");
    const bool communicate_comp = true;
    for (int i = 0; i < num_real_components; ++i) {
        AddRealComp(communicate_comp);
    }
    for (int i = 0; i < num_int_components; ++i) {
        AddIntComp(communicate_comp);
    }

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
    if (owner != iproc) {
        return;
    }

    long num_particles = 0;
    for (const auto& probes : samplers) {
        num_particles += probes->num_points();
    }
    m_total_particles = num_particles;

    const int grid_id = 0;
    const int tile_id = 0;
    // Setup has already processed runtime array information, safe to do
    // GetParticles here.
    auto& ptile = GetParticles(lev)[std::make_pair(grid_id, tile_id)];
    ptile.resize(num_particles);

    int pidx = 0;
    const int nextid = static_cast<int>(ParticleType::NextID());
    auto* pstruct = ptile.GetArrayOfStructs()().data();
    SamplerBase::SampleLocType locs;
    for (const auto& probe : samplers) {
        probe->sampling_locations(locs);
        const int npts = static_cast<int>(locs.size());
        const auto probe_id = probe->id();
        amrex::Gpu::DeviceVector<amrex::Array<amrex::Real, AMREX_SPACEDIM>>
            dlocs(npts);
        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, locs.begin(), locs.end(), dlocs.begin());
        const auto* dpos = dlocs.data();

        amrex::ParallelFor(npts, [=] AMREX_GPU_DEVICE(const int ip) noexcept {
            const auto uid = pidx + ip;
            auto& pp = pstruct[uid];
            pp.id() = nextid + uid;
            pp.cpu() = iproc;

            for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                pp.pos(n) = dpos[ip][n];
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

void SamplingContainer::interpolate_derived_fields(
    const DerivedQtyMgr& derived_mgr, const FieldRepo& repo, const int scomp)
{
    BL_PROFILE("amr-wind::SamplingContainer::interpolate_derived_fields");

    auto outfield = repo.create_scratch_field(derived_mgr.num_comp(), 1);
    derived_mgr(*outfield, 0);

    const int nlevels = m_mesh.finestLevel() + 1;

    for (int lev = 0; lev < nlevels; ++lev) {
        for (ParIterType pti(*this, lev); pti.isValid(); ++pti) {
            const auto farr = (*outfield)(lev).const_array(pti);
            interpolate(
                pti, farr, lev, outfield->field_location(),
                outfield->num_comp(), scomp);
        }
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
            const long offset = fid * num_sampling_particles();
            for (ParIterType pti(*this, lev); pti.isValid(); ++pti) {
                const int np = pti.numParticles();
                auto* pstruct = pti.GetArrayOfStructs()().data();
                auto* parr = pti.GetStructOfArrays().GetRealData(fid).data();

                amrex::ParallelFor(
                    np, [=] AMREX_GPU_DEVICE(const int ip) noexcept {
                        auto& pp = pstruct[ip];
                        const int pidx = pp.idata(IIx::uid);
                        const long ii = offset + pidx;
                        dbuf_ptr[ii] = parr[ip];
                    });
            }
        }
    }

    amrex::Gpu::copy(
        amrex::Gpu::deviceToHost, dbuf.begin(), dbuf.end(), buf.begin());
    amrex::ParallelDescriptor::ReduceRealSum(
        buf.data(), static_cast<int>(buf.size()),
        amrex::ParallelDescriptor::IOProcessorNumber());
}

} // namespace amr_wind::sampling
