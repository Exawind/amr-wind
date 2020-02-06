#include <incflo.H>

using namespace amrex;

// Make a new level using provided BoxArray and DistributionMapping and
// fill with interpolated coarse level data.
// overrides the pure virtual function in AmrCore
void incflo::MakeNewLevelFromCoarse (int lev,
                                     Real time,
                                     const BoxArray& ba,
                                     const DistributionMapping& dm)
{
    BL_PROFILE("incflo::MakeNewLevelFromCoarse()");

    if (m_verbose > 0) {
        amrex::Print() << "Making new level " << lev << " from coarse" << std::endl;
    }

#ifdef AMREX_USE_EB
    std::unique_ptr<FabFactory<FArrayBox> > new_fact = makeEBFabFactory(geom[lev], ba, dm,
                                                                        {nghost_eb_basic(),
                                                                         nghost_eb_volume(),
                                                                         nghost_eb_full()},
                                                                        EBSupport::full);
#else
    std::unique_ptr<FabFactory<FArrayBox> > new_fact(new FArrayBoxFactory());
#endif
    std::unique_ptr<LevelData> new_leveldata
        (new LevelData(ba, dm, *new_fact, m_ntrac, nghost_state(),
                       m_use_godunov,
                       m_diff_type==DiffusionType::Implicit,
                       m_advect_tracer));
    fillcoarsepatch_velocity(lev, time, new_leveldata->velocity, 0);
    fillcoarsepatch_density(lev, time, new_leveldata->density, 0);
    if (m_ntrac > 0) {
        fillcoarsepatch_tracer(lev, time, new_leveldata->tracer, 0);
    }
    fillcoarsepatch_gradp(lev, time, new_leveldata->gp, 0);
    new_leveldata->p.setVal(0.0);

    m_leveldata[lev] = std::move(new_leveldata);
    m_factory[lev] = std::move(new_fact);

    m_diffusion_tensor_op.reset();
    m_diffusion_scalar_op.reset();
}

// Remake an existing level using provided BoxArray and DistributionMapping and
// fill with existing fine and coarse data.
// overrides the pure virtual function in AmrCore
void incflo::RemakeLevel (int lev, Real time, const BoxArray& ba,
			 const DistributionMapping& dm)
{
    BL_PROFILE("incflo::RemakeLevel()");

    if (m_verbose > 0) {
        amrex::Print() << "Remaking level " << lev << std::endl;
    }

#ifdef AMREX_USE_EB
    std::unique_ptr<FabFactory<FArrayBox> > new_fact = makeEBFabFactory(geom[lev], ba, dm,
                                                                        {nghost_eb_basic(),
                                                                         nghost_eb_volume(),
                                                                         nghost_eb_full()},
                                                                        EBSupport::full);
#else
    std::unique_ptr<FabFactory<FArrayBox> > new_fact(new FArrayBoxFactory());
#endif
    std::unique_ptr<LevelData> new_leveldata
        (new LevelData(ba, dm, *new_fact, m_ntrac, nghost_state(),
                       m_use_godunov,
                       m_diff_type==DiffusionType::Implicit,
                       m_advect_tracer));
    fillpatch_velocity(lev, time, new_leveldata->velocity, 0);
    fillpatch_density(lev, time, new_leveldata->density, 0);
    if (m_ntrac > 0) {
        fillpatch_tracer(lev, time, new_leveldata->tracer, 0);
    }
    fillpatch_gradp(lev, time, new_leveldata->gp, 0);
    new_leveldata->p.setVal(0.0);

    m_leveldata[lev] = std::move(new_leveldata);
    m_factory[lev] = std::move(new_fact);

    m_diffusion_tensor_op.reset();
    m_diffusion_scalar_op.reset();
}

// Delete level data
// overrides the pure virtual function in AmrCore
void incflo::ClearLevel (int lev)
{
    BL_PROFILE("incflo::ClearLevel()");
    m_leveldata[lev].reset();
    m_factory[lev].reset();
    m_diffusion_tensor_op.reset();
    m_diffusion_scalar_op.reset();
}
