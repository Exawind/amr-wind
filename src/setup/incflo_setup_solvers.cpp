#include <incflo_proj_F.H>
#include <diffusion_F.H>
#include <incflo.H>

#include <NodalProjection.H>
#include <DiffusionOp.H>

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

    nodal_projector.reset(new NodalProjection(geom, grids, dmap, ppe_lobc, ppe_hibc,
                                              GetVecOfConstPtrs(ebfactory)));

    //
    // Now the diffusion solver
    //
    set_diff_bc( bc_lo, bc_hi,
                 domain.loVect(), domain.hiVect(),
                 &nghost,
                 bc_ilo[0]->dataPtr(), bc_ihi[0]->dataPtr(),
                 bc_jlo[0]->dataPtr(), bc_jhi[0]->dataPtr(),
                 bc_klo[0]->dataPtr(), bc_khi[0]->dataPtr());

    diff_lobc = {(LinOpBCType)bc_lo[0], (LinOpBCType)bc_lo[1], (LinOpBCType)bc_lo[2]};
    diff_hibc = {(LinOpBCType)bc_hi[0], (LinOpBCType)bc_hi[1], (LinOpBCType)bc_hi[2]};

    diffusion_op.reset(new DiffusionOp(this, &ebfactory, diff_lobc, diff_hibc, nghost));
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

    nodal_projector.reset(new NodalProjection(geom, grids, dmap, ppe_lobc, ppe_hibc,
                                              GetVecOfConstPtrs(ebfactory)));

    //
    // Now the diffusion solver
    //
    set_diff_bc( bc_lo, bc_hi,
                 domain.loVect(), domain.hiVect(),
                 &nghost,
                 bc_ilo[0]->dataPtr(), bc_ihi[0]->dataPtr(),
                 bc_jlo[0]->dataPtr(), bc_jhi[0]->dataPtr(),
                 bc_klo[0]->dataPtr(), bc_khi[0]->dataPtr());

    diff_lobc = {(LinOpBCType)bc_lo[0], (LinOpBCType)bc_lo[1], (LinOpBCType)bc_lo[2]};
    diff_hibc = {(LinOpBCType)bc_hi[0], (LinOpBCType)bc_hi[1], (LinOpBCType)bc_hi[2]};

    diffusion_op->setup(this, &ebfactory);
}
