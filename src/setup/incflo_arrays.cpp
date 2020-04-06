#include <incflo.H>
#include "FieldRepo.H"

using namespace amrex;

LevelData::LevelData(int lev, amr_wind::FieldRepo& repo)
    : velocity(repo.get_field("velocity")(lev))
    , velocity_o(repo.get_field("velocity", amr_wind::FieldState::Old)(lev))
    , density(repo.get_field("density")(lev))
    , density_o(repo.get_field("density", amr_wind::FieldState::Old)(lev))
    , tracer(repo.get_field("tracer")(lev))
    , tracer_o(repo.get_field("tracer", amr_wind::FieldState::Old)(lev))
    , gp(repo.get_field("gp")(lev))
    , p(repo.get_field("p")(lev))
    , conv_velocity(repo.get_field("conv_velocity")(lev))
    , conv_velocity_o(repo.get_field("conv_velocity", amr_wind::FieldState::Old)(lev))
    , conv_density(repo.get_field("conv_density")(lev))
    , conv_density_o(repo.get_field("conv_density", amr_wind::FieldState::Old)(lev))
    , conv_tracer(repo.get_field("conv_tracer")(lev))
    , conv_tracer_o(repo.get_field("conv_tracer", amr_wind::FieldState::Old)(lev))
    , divtau(repo.get_field("divtau")(lev))
    , divtau_o(repo.get_field("divtau", amr_wind::FieldState::Old)(lev))
    , laps(repo.get_field("laps")(lev))
    , laps_o(repo.get_field("laps", amr_wind::FieldState::Old)(lev))
{}

#if 0
LevelData::LevelData (amrex::BoxArray const& ba,
                              amrex::DistributionMapping const& dm,
                              amrex::FabFactory<FArrayBox> const& fact,
                              int ntrac, int ng_state,
                              bool use_godunov, bool implicit_diffusion,
                              bool advect_tracer)
    : velocity  (ba, dm, AMREX_SPACEDIM, ng_state, MFInfo(), fact),
      velocity_o(ba, dm, AMREX_SPACEDIM, ng_state, MFInfo(), fact),
      density   (ba, dm, 1             , ng_state, MFInfo(), fact),
      density_o (ba, dm, 1             , ng_state, MFInfo(), fact),
      tracer    (ba, dm, ntrac         , ng_state, MFInfo(), fact),
      tracer_o  (ba, dm, ntrac         , ng_state, MFInfo(), fact),
      gp        (ba, dm, AMREX_SPACEDIM, 0       , MFInfo(), fact),
      p         (amrex::convert(ba,IntVect::TheNodeVector()),
                     dm, 1             , 0 , MFInfo(), fact),
      conv_velocity_o(ba, dm, AMREX_SPACEDIM, 0, MFInfo(), fact),
      conv_density_o (ba, dm, 1             , 0, MFInfo(), fact),
      conv_tracer_o  (ba, dm, ntrac         , 0, MFInfo(), fact)
{
    if (use_godunov) {
        divtau_o.define(ba, dm, AMREX_SPACEDIM, 0, MFInfo(), fact);
        if (advect_tracer) {
            laps_o.define(ba, dm, ntrac, 0, MFInfo(), fact);
        }
    } else {
        conv_velocity.define(ba, dm, AMREX_SPACEDIM, 0, MFInfo(), fact);
        conv_density.define (ba, dm, 1             , 0, MFInfo(), fact);
        conv_tracer.define (ba, dm, ntrac         , 0, MFInfo(), fact);
        if (!implicit_diffusion) {
            divtau.define  (ba, dm, AMREX_SPACEDIM, 0, MFInfo(), fact);
            divtau_o.define(ba, dm, AMREX_SPACEDIM, 0, MFInfo(), fact);
            if (advect_tracer) {
                laps.define  (ba, dm, ntrac, 0, MFInfo(), fact);
                laps_o.define(ba, dm, ntrac, 0, MFInfo(), fact);
            }
        }
    }

}
#endif

// Resize all arrays when instance of incflo class is constructed.
// This is only done at the very start of the simulation. 
void incflo::ResizeArrays ()
{
    // Time holders for fillpatch stuff
    m_t_new.resize(max_level + 1);
    m_t_old.resize(max_level + 1);

    m_leveldata.resize(max_level+1);

    m_factory.resize(max_level+1);
}
