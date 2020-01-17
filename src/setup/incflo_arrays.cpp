#include <incflo.H>

using namespace amrex;

incflo::LevelData::LevelData (amrex::BoxArray const& ba,
                              amrex::DistributionMapping const& dm,
                              amrex::FabFactory<FArrayBox> const& fact,
                              int ntrac, int ng_state, int ng_force)
    : velocity  (ba, dm, AMREX_SPACEDIM, ng_state, MFInfo(), fact),
      velocity_o(ba, dm, AMREX_SPACEDIM, ng_state, MFInfo(), fact),
      density   (ba, dm, 1             , ng_state, MFInfo(), fact),
      density_o (ba, dm, 1             , ng_state, MFInfo(), fact),
      tracer    (ba, dm, ntrac         , ng_state, MFInfo(), fact),
      tracer_o  (ba, dm, ntrac         , ng_state, MFInfo(), fact),
      gp        (ba, dm, AMREX_SPACEDIM, 0       , MFInfo(), fact),
      vel_forces(ba, dm, AMREX_SPACEDIM, ng_force, MFInfo(), fact),
      tra_forces(ba, dm, ntrac         , ng_force, MFInfo(), fact),
      p         (amrex::convert(ba,IntVect::TheNodeVector()),
                     dm, 1             , 0 , MFInfo(), fact),
      conv_velocity  (ba, dm, AMREX_SPACEDIM, 0, MFInfo(), fact),
      conv_velocity_o(ba, dm, AMREX_SPACEDIM, 0, MFInfo(), fact),
      conv_density   (ba, dm, 1             , 0, MFInfo(), fact),
      conv_density_o (ba, dm, 1             , 0, MFInfo(), fact),
      conv_tracer    (ba, dm, ntrac         , 0, MFInfo(), fact),
      conv_tracer_o  (ba, dm, ntrac         , 0, MFInfo(), fact)
{
    // xxxxx TODO we probably do not need the new conv_* for godunov
}

// Resize all arrays when instance of incflo class is constructed.
// This is only done at the very start of the simulation. 
void incflo::ResizeArrays ()
{
    // Time holders for fillpatch stuff
    t_new.resize(max_level + 1);
    t_old.resize(max_level + 1);

    m_leveldata.resize(max_level+1);

    m_factory.resize(max_level+1);
}
