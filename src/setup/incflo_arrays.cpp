#include <incflo.H>

using namespace amrex;

incflo::LevelData::LevelData (amrex::BoxArray const& ba,
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
        if (!implicit_diffusion) 
        {
            divtau.define  (ba, dm, AMREX_SPACEDIM, 0, MFInfo(), fact);
            divtau_o.define(ba, dm, AMREX_SPACEDIM, 0, MFInfo(), fact);
            if (advect_tracer) {
                laps.define  (ba, dm, ntrac, 0, MFInfo(), fact);
                laps_o.define(ba, dm, ntrac, 0, MFInfo(), fact);
            }
        }
    }
}

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
