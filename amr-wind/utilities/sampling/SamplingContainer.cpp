
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

    const int lev = 0;
    const auto iproc = amrex::ParallelDescriptor::MyProc();

    m_total_particles = 0;
    for (const auto& probes : samplers) {
        m_total_particles += probes->num_points();
    }

    const int nprobes = static_cast<int>(samplers.size());
    amrex::iMultiFab particle_counts(
        ParticleBoxArray(lev), ParticleDistributionMap(lev), nprobes, 0,
        amrex::MFInfo());
    amrex::iMultiFab offsets(
        ParticleBoxArray(lev), ParticleDistributionMap(lev), nprobes, 0,
        amrex::MFInfo());
    particle_counts.setVal(0);
    offsets.setVal(0);

    amrex::Vector<amrex::Box> owned_boxes;
    for (auto const idx : particle_counts.IndexArray()) {
        owned_boxes.push_back(particle_counts.boxArray()[idx]);
    }

    const auto& dxinv = m_mesh.Geom(lev).InvCellSizeArray();
    const auto& plo = m_mesh.Geom(lev).ProbLoArray();
    for (int iprobe = 0; iprobe < nprobes; iprobe++) {
        const auto& probe = samplers[iprobe];
        SampleLocType sample_locs;
        probe->sampling_locations(sample_locs, owned_boxes);
        const auto& locs = sample_locs.locations();
        const int npts = static_cast<int>(locs.size());
        amrex::Gpu::DeviceVector<amrex::RealVect> dlocs(npts);

        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, locs.begin(), locs.end(), dlocs.begin());
        const auto* dpos = dlocs.data();

        // don't use openmp
        for (amrex::MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi) {
            const amrex::Box& box = mfi.tilebox();

            // count the number of particles in each cell
            const auto& np_arr = particle_counts[mfi].array();
            amrex::ParallelFor(
                dlocs.size(), [=] AMREX_GPU_DEVICE(long ip) noexcept {
                    const amrex::RealVect pos(
                        AMREX_D_DECL(dpos[ip][0], dpos[ip][1], dpos[ip][2]));

                    const amrex::IntVect div(AMREX_D_DECL(
                        static_cast<int>(
                            amrex::Math::floor((pos[0] - plo[0]) * dxinv[0])),
                        static_cast<int>(
                            amrex::Math::floor((pos[1] - plo[1]) * dxinv[1])),
                        static_cast<int>(
                            amrex::Math::floor((pos[2] - plo[2]) * dxinv[2]))));
                    if (box.contains(div)) {
                        amrex::Gpu::Atomic::AddNoRet(&np_arr(div, iprobe), 1);
                    }
                });
        }
        AMREX_ALWAYS_ASSERT(particle_counts.sum(iprobe) == npts);
    }

    // compute the offsets and size the particle tiles
    // don't use openmp
    for (amrex::MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi) {
        const amrex::Box& box = mfi.tilebox();
        const auto ncells = static_cast<int>(box.numPts());

        const auto& np_arr = particle_counts[mfi].const_array();
        int* p_offsets = offsets[mfi].dataPtr();
        const auto np_box = amrex::Scan::PrefixSum<int>(
            ncells,
            [=] AMREX_GPU_DEVICE(int i) -> int {
                const auto iv = box.atOffset(i);
                int total_np = 0;
                for (int iprobe = 0; iprobe < nprobes; iprobe++) {
                    total_np += np_arr(iv, iprobe);
                }
                return total_np;
            },
            [=] AMREX_GPU_DEVICE(int i, const int& xi) { p_offsets[i] = xi; },
            amrex::Scan::Type::exclusive);

        const auto& offsets_arr = offsets[mfi].array();
        amrex::ParallelFor(
            box, [=] AMREX_GPU_DEVICE(
                     int i, int j, int AMREX_D_PICK(, , k)) noexcept {
                const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
                for (int iprobe = 1; iprobe < nprobes; iprobe++) {
                    offsets_arr(iv, iprobe) =
                        offsets_arr(iv, iprobe - 1) + np_arr(iv, iprobe - 1);
                }
            });

        const int grid_id = mfi.index();
        const int tile_id = mfi.LocalTileIndex();
        auto& ptile = GetParticles(lev)[std::make_pair(grid_id, tile_id)];
        ptile.resize(np_box);

        AMREX_ASSERT(
            np_box ==
            particle_counts[mfi].sum<amrex::RunOn::Device>(box, 0, nprobes));
    }

    for (int iprobe = 0; iprobe < nprobes; iprobe++) {
        const auto& probe = samplers[iprobe];
        const auto probe_id = probe->id();
        SampleLocType sample_locs;
        probe->sampling_locations(sample_locs, owned_boxes);
        const auto& locs = sample_locs.locations();
        const int npts = static_cast<int>(locs.size());
        amrex::Gpu::DeviceVector<amrex::RealVect> dlocs(npts);
        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, locs.begin(), locs.end(), dlocs.begin());
        const auto* p_dlocs = dlocs.data();
        const auto& ids = sample_locs.ids();
        amrex::Gpu::DeviceVector<amrex::Long> dids(npts);
        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, ids.begin(), ids.end(), dids.begin());
        const auto* p_dids = dids.data();

        // don't use openmp
        for (amrex::MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi) {
            const amrex::Box& box = mfi.tilebox();
            const auto& np_arr = particle_counts[mfi].const_array();
            const auto& offset_arr = offsets[mfi].const_array();

            const int grid_id = mfi.index();
            const int tile_id = mfi.LocalTileIndex();
            auto& ptile = GetParticles(lev)[std::make_pair(grid_id, tile_id)];
            if (ptile.numParticles() == 0) {
                continue;
            }

            const int np_probe =
                particle_counts[mfi].sum<amrex::RunOn::Device>(box, iprobe, 1);
            if (np_probe == 0) {
                continue;
            }

            const amrex::Long nextid = ParticleType::NextID();
            ParticleType::NextID(nextid + np_probe);

            auto* pstruct = ptile.GetArrayOfStructs()().data();
            amrex::ParallelFor(
                box, [=] AMREX_GPU_DEVICE(
                         int i, int j, int AMREX_D_PICK(, , k)) noexcept {
                    const amrex::IntVect iv(AMREX_D_DECL(i, j, k));

                    if (np_arr(iv, iprobe) > 0) {
                        const int start = offset_arr(iv, iprobe);
                        int n = start;
                        for (int ip = 0; ip < npts; ip++) {
                            const amrex::RealVect loc(AMREX_D_DECL(
                                p_dlocs[ip][0], p_dlocs[ip][1],
                                p_dlocs[ip][2]));

                            const amrex::IntVect div(AMREX_D_DECL(
                                static_cast<int>(amrex::Math::floor(
                                    (loc[0] - plo[0]) * dxinv[0])),
                                static_cast<int>(amrex::Math::floor(
                                    (loc[1] - plo[1]) * dxinv[1])),
                                static_cast<int>(amrex::Math::floor(
                                    (loc[2] - plo[2]) * dxinv[2]))));
                            if (div == iv) {
                                auto& pp = pstruct[n];
                                pp.id() = nextid + n;
                                pp.cpu() = iproc;
                                for (int idim = 0; idim < AMREX_SPACEDIM;
                                     ++idim) {
                                    pp.pos(idim) = loc[idim];
                                }
                                pp.idata(IIx::uid) = ip + npts * iprobe;
                                pp.idata(IIx::sid) = probe_id;
                                pp.idata(IIx::nid) =
                                    static_cast<int>(p_dids[ip]);
                                n++;
                            }
                        }
                        AMREX_ALWAYS_ASSERT(np_arr(iv, iprobe) == (n - start));
                    }
                });
            amrex::Gpu::streamSynchronize();
        }
    }
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
