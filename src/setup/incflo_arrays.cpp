#include <incflo.H>

using namespace amrex;

incflo::LevelData::LevelData (amrex::BoxArray const& ba,
                              amrex::DistributionMapping const& dm,
                              amrex::FabFactory<FArrayBox> const& fact)
    : velocity  (ba, dm, AMREX_SPACEDIM, incflo::nghost, MFInfo(), fact),
      velocity_o(ba, dm, AMREX_SPACEDIM, incflo::nghost, MFInfo(), fact),
      density   (ba, dm, 1             , incflo::nghost, MFInfo(), fact),
      density_o (ba, dm, 1             , incflo::nghost, MFInfo(), fact),
      tracer    (ba, dm, incflo::ntrac , incflo::nghost, MFInfo(), fact),
      tracer_o  (ba, dm, ntrac         , incflo::nghost, MFInfo(), fact),
      gp        (ba, dm, AMREX_SPACEDIM, incflo::nghost, MFInfo(), fact),
      p         (amrex::convert(ba,IntVect::TheNodeVector()),
                     dm, 1             , incflo::nghost, MFInfo(), fact),
      conv_velocity  (ba, dm, AMREX_SPACEDIM, 0, MFInfo(), fact),
      conv_velocity_o(ba, dm, AMREX_SPACEDIM, 0, MFInfo(), fact),
      conv_density   (ba, dm, 1             , 0, MFInfo(), fact),
      conv_density_o (ba, dm, 1             , 0, MFInfo(), fact),
      conv_tracer    (ba, dm, incflo::ntrac , 0, MFInfo(), fact),
      conv_tracer_o  (ba, dm, incflo::ntrac , 0, MFInfo(), fact)
{}

void incflo::AllocateArrays (int lev)
{
    // xxxxx remove?
}

void incflo::RegridArrays(int lev)
{
    // xxxxx remove?
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

void incflo::MakeBCArrays()
{
}

