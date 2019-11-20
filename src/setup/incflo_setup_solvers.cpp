#include <incflo_proj_F.H>
#include <diffusion_F.H>
#include <incflo.H>

void
incflo::incflo_init_solvers ()
{
    BL_PROFILE("incflo::incflo_init_solvers");

    int bc_lo[3], bc_hi[3];
    Box domain(geom[0].Domain());

    //
    // First the nodal projection
    //
    set_ppe_bcs(bc_lo, bc_hi,
                domain.loVect(), domain.hiVect(),
                &nghost,
                bc_ilo[0]->dataPtr(), bc_ihi[0]->dataPtr(),
                bc_jlo[0]->dataPtr(), bc_jhi[0]->dataPtr(),
                bc_klo[0]->dataPtr(), bc_khi[0]->dataPtr());

    ppe_lobc = {(LinOpBCType)bc_lo[0], (LinOpBCType)bc_lo[1], (LinOpBCType)bc_lo[2]};
    ppe_hibc = {(LinOpBCType)bc_hi[0], (LinOpBCType)bc_hi[1], (LinOpBCType)bc_hi[2]};

    LPInfo info;
    info.setMaxCoarseningLevel(100);

#ifdef AMREX_USE_EB
    nodal_projector.reset(new NodalProjector(geom, grids, dmap, ppe_lobc, ppe_hibc,
                                             GetVecOfConstPtrs(ebfactory), info));
#else
    nodal_projector.reset(new NodalProjector(geom, grids, dmap, ppe_lobc, ppe_hibc, info));
#endif

    //
    // Now the diffusion solver
    //
    set_vel_diff_bc( bc_lo, bc_hi,
                    domain.loVect(), domain.hiVect(),
                    &nghost,
                    bc_ilo[0]->dataPtr(), bc_ihi[0]->dataPtr(),
                    bc_jlo[0]->dataPtr(), bc_jhi[0]->dataPtr(),
                    bc_klo[0]->dataPtr(), bc_khi[0]->dataPtr());


    vel_diff_lobc = {(LinOpBCType)bc_lo[0], (LinOpBCType)bc_lo[1], (LinOpBCType)bc_lo[2]};
    vel_diff_hibc = {(LinOpBCType)bc_hi[0], (LinOpBCType)bc_hi[1], (LinOpBCType)bc_hi[2]};

    //
    // Now the diffusion solver
    //
    set_scal_diff_bc( bc_lo, bc_hi,
                     domain.loVect(), domain.hiVect(),
                     &nghost,
                     bc_ilo[0]->dataPtr(), bc_ihi[0]->dataPtr(),
                     bc_jlo[0]->dataPtr(), bc_jhi[0]->dataPtr(),
                     bc_klo[0]->dataPtr(), bc_khi[0]->dataPtr());

    scal_diff_lobc = {(LinOpBCType)bc_lo[0], (LinOpBCType)bc_lo[1], (LinOpBCType)bc_lo[2]};
    scal_diff_hibc = {(LinOpBCType)bc_hi[0], (LinOpBCType)bc_hi[1], (LinOpBCType)bc_hi[2]};

#ifdef AMREX_USE_EB
    diffusion_op.reset(new DiffusionOp(this, &ebfactory, vel_diff_lobc,  vel_diff_hibc,
                                                        scal_diff_lobc, scal_diff_hibc, nghost,probtype));
#else
    diffusion_op.reset(new DiffusionOp(this,             vel_diff_lobc,  vel_diff_hibc,
                                                        scal_diff_lobc, scal_diff_hibc, nghost,probtype));
#endif
}

void
incflo::incflo_setup_solvers ()
{
    BL_PROFILE("incflo::incflo_setup_solvers");

    int bc_lo[3], bc_hi[3];
    Box domain(geom[0].Domain());

    //
    // First the nodal projection
    //
    set_ppe_bcs(bc_lo, bc_hi,
                domain.loVect(), domain.hiVect(),
                &nghost,
                bc_ilo[0]->dataPtr(), bc_ihi[0]->dataPtr(),
                bc_jlo[0]->dataPtr(), bc_jhi[0]->dataPtr(),
                bc_klo[0]->dataPtr(), bc_khi[0]->dataPtr());

    ppe_lobc = {(LinOpBCType)bc_lo[0], (LinOpBCType)bc_lo[1], (LinOpBCType)bc_lo[2]};
    ppe_hibc = {(LinOpBCType)bc_hi[0], (LinOpBCType)bc_hi[1], (LinOpBCType)bc_hi[2]};
    LPInfo info;
    info.setMaxCoarseningLevel(100);

#ifdef AMREX_USE_EB
    nodal_projector.reset(new NodalProjector(geom, grids, dmap, ppe_lobc, ppe_hibc,
                                             GetVecOfConstPtrs(ebfactory), info));
#else
    nodal_projector.reset(new NodalProjector(geom, grids, dmap, ppe_lobc, ppe_hibc, info));
#endif


    //
    // Now the diffusion solver
    //

#ifdef AMREX_USE_EB
    diffusion_op->setup(this, &ebfactory);
#else
    diffusion_op->setup(this);
#endif
}
