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

    Vector<Geometry> lgeom {geom.begin(), geom.begin()+finest_level+1};
    Vector<BoxArray> lgrids{grids.begin(), grids.begin()+finest_level+1};
    Vector<DistributionMapping> ldmap{dmap.begin(), dmap.begin()+finest_level+1};
    Vector<EBFArrayBoxFactory const*> fact = GetVecOfConstPtrs(ebfactory);
    Vector<EBFArrayBoxFactory const*> lfact{fact.begin(), fact.begin()+finest_level+1};

    nodal_projector.reset(new NodalProjector(lgeom, lgrids, ldmap, ppe_lobc, ppe_hibc,
                                             lfact));

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

    diffusion_op.reset(new DiffusionOp(this, &ebfactory, vel_diff_lobc,  vel_diff_hibc,
                                                        scal_diff_lobc, scal_diff_hibc, nghost));
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

    nodal_projector.reset(new NodalProjector(geom, grids, dmap, ppe_lobc, ppe_hibc,
                                             GetVecOfConstPtrs(ebfactory)));


    //
    // Now the diffusion solver
    //

    diffusion_op->setup(this, &ebfactory);
}
