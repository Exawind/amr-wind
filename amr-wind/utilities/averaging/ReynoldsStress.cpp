#include "amr-wind/utilities/averaging/ReynoldsStress.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/core/Field.H"
#include "amr-wind/core/FieldRepo.H"
#include "amr-wind/utilities/IOManager.H"

namespace amr_wind::averaging {
namespace {

const Field& get_field_or_error(const FieldRepo& repo, const std::string& fname)
{
    if (!repo.field_exists(fname)) {
        amrex::Abort("ReynoldsStress: Cannot find field: " + fname);
    }

    return repo.get_field(fname);
}

} // namespace

ReynoldsStress::ReynoldsStress(CFDSim& sim, const std::string& fname)
    : m_field(get_field_or_error(sim.repo(), "velocity"))
    , m_average(get_field_or_error(sim.repo(), "velocity_mean"))
    , m_stress(sim.repo().declare_field(
          "velocity_stress",
          6, // number of components of the reynolds stress tensor
          1, // Ghost cells
          1,
          m_field.field_location()))
    , m_re_stress(sim.repo().declare_field(
          "velocity_reynolds_stress",
          6, // number of components of the reynolds stress tensor
          1, // Ghost cells
          1,
          m_field.field_location()))
{
    if (fname != "velocity") {
        amrex::Abort("ReynoldsStress only implemented for velocity field");
    }

    // Register default fillpatch operations
    m_stress.set_default_fillpatch_bc(sim.time());
    m_re_stress.set_default_fillpatch_bc(sim.time());

    // Do coarse/fine interpolations upon regrid
    m_stress.fillpatch_on_regrid() = true;

    // Register average field with the IO manager
    auto& iomgr = sim.io_manager();
    iomgr.register_io_var(m_stress.name());
    iomgr.register_io_var(m_re_stress.name());
}

const std::string& ReynoldsStress::average_field_name()
{
    return m_re_stress.name();
}

void ReynoldsStress::operator()(
    const SimTime& time,
    const amrex::Real filter_width,
    const amrex::Real elapsed_time)
{
    const amrex::Real dt = time.deltaT();
    const amrex::Real filter =
        amrex::max(amrex::min(filter_width, elapsed_time), dt);
    const amrex::Real factor = amrex::max<amrex::Real>(filter - dt, 0.0);

    const int ncomp = m_field.num_comp();
    const int nlevels = m_field.repo().num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {

        const auto& ffab = m_field(lev);
        const auto& afab = m_average(lev);
        auto& sfab = m_stress(lev);
        auto& rfab = m_re_stress(lev);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (amrex::MFIter mfi(ffab, amrex::TilingIfNotGPU()); mfi.isValid();
             ++mfi) {
            const auto& bx = mfi.tilebox();
            const auto& fldarr = ffab.const_array(mfi);
            const auto& avgarr = afab.array(mfi);
            const auto& stressarr = sfab.array(mfi);
            const auto& restressarr = rfab.array(mfi);

            amrex::ParallelFor(
                bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    // The tensor index
                    int mn = 0;
                    for (int n = 0; n < ncomp; ++n) {
                        for (int m = n; m < ncomp; ++m) {
                            // AB
                            const amrex::Real fval2 =
                                fldarr(i, j, k, m) * fldarr(i, j, k, n);
                            // <A><B>
                            const amrex::Real aval2 =
                                avgarr(i, j, k, m) * avgarr(i, j, k, n);
                            // The current value
                            const amrex::Real avg = stressarr(i, j, k, mn);
                            // The stress <AB>
                            stressarr(i, j, k, mn) =
                                (avg * factor + fval2 * dt) / filter;
                            // The Reynolds stress <ab>
                            restressarr(i, j, k, mn) =
                                stressarr(i, j, k, mn) - aval2;
                            ++mn;
                        }
                    }
                });
        }
    }

    m_stress.fillpatch(time.new_time());
    m_re_stress.fillpatch(time.new_time());
}

} // namespace amr_wind::averaging
