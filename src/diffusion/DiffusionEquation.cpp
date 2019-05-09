#include <AMReX_EBFArrayBox.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Vector.H>

#include <DiffusionEquation.H>
#include <diffusion_F.H>
#include <constants.H>

using namespace amrex;

//
// Constructor:
// We set up everything which doesn't change between timesteps here
//
DiffusionEquation::DiffusionEquation(AmrCore* _amrcore,
                                     Vector<std::unique_ptr<EBFArrayBoxFactory>>* _ebfactory,
                                     Vector<std::unique_ptr<IArrayBox>>& bc_ilo,
                                     Vector<std::unique_ptr<IArrayBox>>& bc_ihi,
                                     Vector<std::unique_ptr<IArrayBox>>& bc_jlo,
                                     Vector<std::unique_ptr<IArrayBox>>& bc_jhi,
                                     Vector<std::unique_ptr<IArrayBox>>& bc_klo,
                                     Vector<std::unique_ptr<IArrayBox>>& bc_khi,
                                     int _nghost, Real _cyl_speed)
{
    // Get inputs from ParmParse
	readParameters();

    if(verbose > 0)
    {
        amrex::Print() << "Constructing DiffusionEquation class" << std::endl;
    }

    // Set AmrCore and ebfactory based on input, fetch some data needed in constructor
    amrcore = _amrcore;
    ebfactory = _ebfactory;
    nghost = _nghost;
    Vector<Geometry> geom = amrcore->Geom();
    Vector<BoxArray> grids = amrcore->boxArray();
    Vector<DistributionMapping> dmap = amrcore->DistributionMap();
    int max_level = amrcore->maxLevel();

    // Cylinder speed
    cyl_speed = _cyl_speed;

    // Whole domain
    Box domain(geom[0].Domain());

    // The boundary conditions need only be set at level 0
    set_diff_bc(bc_lo, bc_hi,
                domain.loVect(), domain.hiVect(),
				&nghost,
                bc_ilo[0]->dataPtr(), bc_ihi[0]->dataPtr(),
                bc_jlo[0]->dataPtr(), bc_jhi[0]->dataPtr(),
                bc_klo[0]->dataPtr(), bc_khi[0]->dataPtr());

    // Resize and reset data
    b.resize(max_level + 1);
    phi.resize(max_level + 1);
    rhs.resize(max_level + 1);
    ueb.resize(max_level + 1);
    veb.resize(max_level + 1);
    for(int lev = 0; lev <= max_level; lev++)
    {
        for(int dir = 0; dir < 3; dir++)
        {
            BoxArray edge_ba = grids[lev];
            edge_ba.surroundingNodes(dir);
            b[lev][dir].reset(new MultiFab(edge_ba, dmap[lev], 1, nghost,
                                           MFInfo(), *(*ebfactory)[lev]));
        }
        phi[lev].reset(new MultiFab(grids[lev], dmap[lev], 1, nghost,
                                    MFInfo(), *(*ebfactory)[lev]));
        rhs[lev].reset(new MultiFab(grids[lev], dmap[lev], 1, nghost,
                                    MFInfo(), *(*ebfactory)[lev]));
        ueb[lev].reset(new MultiFab(grids[lev], dmap[lev], 1, nghost,
                                    MFInfo(), *(*ebfactory)[lev]));
        veb[lev].reset(new MultiFab(grids[lev], dmap[lev], 1, nghost,
                                    MFInfo(), *(*ebfactory)[lev]));
    }

    // Fill the Dirichlet values on the EB surface
    for(int lev = 0; lev <= max_level; lev++)
    {
        // Get EB normal vector
        const amrex::MultiCutFab*                 bndrynormal;
        bndrynormal = &((*ebfactory)[lev] -> getBndryNormal());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for(MFIter mfi(*ueb[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            // Tilebox
            Box bx = mfi.tilebox();

            // This is to check efficiently if this tile contains any eb stuff
            const EBFArrayBox& ueb_fab = static_cast<EBFArrayBox const&>((*ueb[lev])[mfi]);
            const EBCellFlagFab& flags = ueb_fab.getEBCellFlagFab();

            if (flags.getType(bx) == FabType::covered || flags.getType(bx) == FabType::regular)
            {
                (*ueb[lev])[mfi].setVal(0.0, bx);
                (*veb[lev])[mfi].setVal(0.0, bx);
            }
            else
            {
                const auto& ueb_arr = ueb[lev]->array(mfi);
                const auto& veb_arr = veb[lev]->array(mfi);
                const auto& nrm_fab = bndrynormal->array(mfi);

                for(int i = bx.smallEnd(0); i <= bx.bigEnd(0); i++)
                for(int j = bx.smallEnd(1); j <= bx.bigEnd(1); j++)
                for(int k = bx.smallEnd(2); k <= bx.bigEnd(2); k++)
                {
                    Real theta = atan2(-nrm_fab(i,j,k,1), -nrm_fab(i,j,k,0));
                    ueb_arr(i,j,k) =   cyl_speed * sin(theta);
                    veb_arr(i,j,k) = - cyl_speed * cos(theta);
                }
            }
        }
    }

	// Define the matrix.
	LPInfo info;
    info.setMaxCoarseningLevel(mg_max_coarsening_level);
    matrix.define(geom, grids, dmap, info, GetVecOfConstPtrs(*ebfactory));

    // It is essential that we set MaxOrder to 2 if we want to use the standard
    // phi(i)-phi(i-1) approximation for the gradient at Dirichlet boundaries.
    // The solver's default order is 3 and this uses three points for the gradient.
	matrix.setMaxOrder(2);

	// LinOpBCType Definitions are in amrex/Src/Boundary/AMReX_LO_BCTYPES.H
	matrix.setDomainBC({(LinOpBCType) bc_lo[0], (LinOpBCType) bc_lo[1], (LinOpBCType) bc_lo[2]},
					   {(LinOpBCType) bc_hi[0], (LinOpBCType) bc_hi[1], (LinOpBCType) bc_hi[2]});
}

DiffusionEquation::~DiffusionEquation()
{
}

void DiffusionEquation::readParameters()
{
    ParmParse pp("diffusion");

    pp.query("verbose", verbose);
    pp.query("mg_verbose", mg_verbose);
    pp.query("mg_cg_verbose", mg_cg_verbose);
    pp.query("mg_max_iter", mg_max_iter);
    pp.query("mg_cg_maxiter", mg_cg_maxiter);
    pp.query("mg_max_fmg_iter", mg_max_fmg_iter);
    pp.query("mg_max_coarsening_level", mg_max_coarsening_level);
    pp.query("mg_rtol", mg_rtol);
    pp.query("mg_atol", mg_atol);
    pp.query("bottom_solver_type", bottom_solver_type);
}

void DiffusionEquation::updateInternals(AmrCore* amrcore_in,
                                        Vector<std::unique_ptr<EBFArrayBoxFactory>>* ebfactory_in)
{
    // This must be implemented when we want dynamic meshing
    //
    amrex::Print() << "ERROR: DiffusionEquation::updateInternals() not yet implemented" << std::endl;
    amrex::Abort();
}

//
// Solve the matrix equation
//
void DiffusionEquation::solve(Vector<std::unique_ptr<MultiFab>>& vel,
                              const Vector<std::unique_ptr<MultiFab>>& ro,
                              const Vector<std::unique_ptr<MultiFab>>& eta,
                              Real dt)
{
	BL_PROFILE("DiffusionEquation::solve");

    // Update the coefficients of the matrix going into the solve based on the current state of the
    // simulation. Recall that the relevant matrix is
    //
    //      alpha a - beta div ( b grad )   <--->   rho - dt div ( eta grad )
    //
    // So the constants and variable coefficients are:
    //
    //      alpha: 1
    //      beta: dt
    //      a: ro
    //      b: eta

    // Set alpha and beta
    matrix.setScalars(1.0, dt);

    for(int lev = 0; lev <= amrcore->finestLevel(); lev++)
    {
        // Compute the spatially varying b coefficients (on faces) to equal the apparent viscosity
        average_cellcenter_to_face(GetArrOfPtrs(b[lev]), *eta[lev], amrcore->Geom(lev));
        for(int dir = 0; dir < 3; dir++)
        {
            b[lev][dir]->FillBoundary(amrcore->Geom(lev).periodicity());
        }
        
        // This sets the coefficients
        matrix.setACoeffs(lev, (*ro[lev]));
        matrix.setBCoeffs(lev, GetArrOfConstPtrs(b[lev])); 
    }

    if(verbose > 0)
    {
        amrex::Print() << "Diffusing velocity..." << std::endl; 
    }

    // Loop over the velocity components
    for(int dir = 0; dir < 3; dir++)
    {
        for(int lev = 0; lev <= amrcore->finestLevel(); lev++)
        {
            // Set the right hand side to equal rho
            rhs[lev]->copy(*ro[lev], 0, 0, 1, nghost, nghost);

            // Multiply rhs by vel(dir) to get momentum
            // Note that vel holds the updated velocity:
            //
            //      u_old + dt ( - u grad u + div ( eta (grad u)^T ) / rho - grad p / rho + gravity )
            //
            MultiFab::Multiply((*rhs[lev]), (*vel[lev]), dir, 0, 1, nghost);

            // By this point we must have filled the Dirichlet values of phi stored in ghost cells
            phi[lev]->copy(*vel[lev], dir, 0, 1, nghost, nghost);
            phi[lev]->FillBoundary(amrcore->Geom(lev).periodicity());
            matrix.setLevelBC(lev, GetVecOfConstPtrs(phi)[lev]);

            // This sets the coefficient on the wall and defines the wall as a Dirichlet bc
            if(cyl_speed > 0.0 && dir == 0)
            {
                matrix.setEBDirichlet(lev, *ueb[lev], *eta[lev]);
            }
            else if(cyl_speed > 0.0 && dir == 1)
            {
                matrix.setEBDirichlet(lev, *veb[lev], *eta[lev]);
            }
            else
            {
                matrix.setEBHomogDirichlet(lev, *eta[lev]);
            }
        }

        MLMG solver(matrix);
        setSolverSettings(solver);
        solver.solve(GetVecOfPtrs(phi), GetVecOfConstPtrs(rhs), mg_rtol, mg_atol);

        for(int lev = 0; lev <= amrcore->finestLevel(); lev++)
        {
            phi[lev]->FillBoundary(amrcore->Geom(lev).periodicity());
            vel[lev]->copy(*phi[lev], 0, dir, 1, nghost, nghost);
        }

        if(verbose > 0)
        {
            amrex::Print() << " done!" << std::endl;
        }

    }
}

//
// Set the user-supplied settings for the MLMG solver
// (this must be done every time step, since MLMG is created after updating matrix
//
void DiffusionEquation::setSolverSettings(MLMG& solver)
{
    // The default bottom solver is BiCG
    if(bottom_solver_type == "smoother")
    {
       solver.setBottomSolver(MLMG::BottomSolver::smoother);
    }
    else if(bottom_solver_type == "hypre")
    {
       solver.setBottomSolver(MLMG::BottomSolver::hypre);
    }

    // Maximum iterations for MultiGrid / ConjugateGradients
	solver.setMaxIter(mg_max_iter);
	solver.setMaxFmgIter(mg_max_fmg_iter);
	solver.setCGMaxIter(mg_cg_maxiter);

    // Verbosity for MultiGrid / ConjugateGradients
	solver.setVerbose(mg_verbose);
	solver.setCGVerbose(mg_cg_verbose);

	// This ensures that ghost cells of phi are correctly filled when returned from the solver
	solver.setFinalFillBC(true);
}

