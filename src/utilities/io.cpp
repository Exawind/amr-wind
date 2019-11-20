#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_buildInfo.H>

#include <incflo.H>

namespace { const std::string level_prefix{"Level_"}; }

void GotoNextLine(std::istream& is)
{
	constexpr std::streamsize bl_ignore_max{100000};
	is.ignore(bl_ignore_max, '\n');
}

void incflo::WriteHeader(
	const std::string& name, bool is_checkpoint) const
{
	if(ParallelDescriptor::IOProcessor())
	{
		std::string HeaderFileName(name + "/Header");
		VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
		std::ofstream HeaderFile;

		HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

		HeaderFile.open(HeaderFileName.c_str(),
						std::ofstream::out | std::ofstream::trunc | std::ofstream::binary);

		if(!HeaderFile.good())
			amrex::FileOpenFailed(HeaderFileName);

		HeaderFile.precision(17);

		if(is_checkpoint)
			HeaderFile << "Checkpoint version: 1\n";
		else
			HeaderFile << "HyperCLaw-V1.1\n";

		HeaderFile << finest_level << "\n";

		// Time stepping controls
		HeaderFile << nstep << "\n";
		HeaderFile << dt << "\n";
		HeaderFile << cur_time << "\n";

		// Geometry
		for(int i = 0; i < BL_SPACEDIM; ++i)
			HeaderFile << Geom(0).ProbLo(i) << ' ';
		HeaderFile << '\n';

		for(int i = 0; i < BL_SPACEDIM; ++i)
			HeaderFile << Geom(0).ProbHi(i) << ' ';
		HeaderFile << '\n';

		// BoxArray
		for(int lev = 0; lev <= finest_level; ++lev)
		{
			boxArray(lev).writeOn(HeaderFile);
			HeaderFile << '\n';
		}
	}
}

void incflo::WriteCheckPointFile() const
{
	BL_PROFILE("incflo::WriteCheckPointFile()");

	const std::string& checkpointname = amrex::Concatenate(check_file, nstep);

    amrex::Print() << "\n\t Writing checkpoint " << checkpointname << std::endl;

	amrex::PreBuildDirectorHierarchy(checkpointname, level_prefix, finest_level + 1, true);

    bool is_checkpoint = true;
	WriteHeader(checkpointname, is_checkpoint);
	WriteJobInfo(checkpointname);

	for(int lev = 0; lev <= finest_level; ++lev)
	{

		// This writes all three velocity components
		VisMF::Write((*vel[lev]),
			     amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, vecVarsName[0]));

		// This writes all three pressure gradient components
		VisMF::Write((*gp[lev]),
			      amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, vecVarsName[3]));

		// Write scalar variables
		for(int i = 0; i < chkscalarVars.size(); i++)
		{
		   VisMF::Write(*((*chkscalarVars[i])[lev]),
				 amrex::MultiFabFileFullPrefix(
				 lev, checkpointname, level_prefix, chkscaVarsName[i]));
		}
	}
}

void incflo::ReadCheckpointFile()
{
	BL_PROFILE("incflo::ReadCheckpointFile()");

	amrex::Print() << "Restarting from checkpoint " << restart_file << std::endl;

	Real prob_lo[BL_SPACEDIM];
	Real prob_hi[BL_SPACEDIM];

	/***************************************************************************
     * Load header: set up problem domain (including BoxArray)                 *
     *              allocate incflo memory (incflo::AllocateArrays)            *
     *              (by calling MakeNewLevelFromScratch)
     ***************************************************************************/

    std::string File(restart_file + "/Header");

    VisMF::IO_Buffer io_buffer(VisMF::GetIOBufferSize());

    Vector<char> fileCharPtr;
    ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
    std::string fileCharPtrString(fileCharPtr.dataPtr());
    std::istringstream is(fileCharPtrString, std::istringstream::in);

    std::string line, word;

    // Start reading from checkpoint file 
    
    // Title line
    std::getline(is, line);

    // Finest level
    is >> finest_level;
    GotoNextLine(is);

    // Step count
    is >> nstep;
    GotoNextLine(is);

    // Time step size
    is >> dt;
    GotoNextLine(is);

    // Current time
    is >> cur_time;
    GotoNextLine(is);

    // Low coordinates of domain bounding box
    std::getline(is, line);
    {
        std::istringstream lis(line);
        int i = 0;
        while(lis >> word)
        {
            prob_lo[i++] = std::stod(word);
        }
    }

    // High coordinates of domain bounding box
    std::getline(is, line);
    {
        std::istringstream lis(line);
        int i = 0;
        while(lis >> word)
        {
            prob_hi[i++] = std::stod(word);
        }
    }

    // Set up problem domain
    RealBox rb(prob_lo, prob_hi);
    Geometry::ResetDefaultProbDomain(rb);
    for (int lev = 0; lev <= max_level; ++lev) {
        SetGeometry(lev, Geometry(Geom(lev).Domain(), rb, Geom(lev).CoordInt(),
                                  Geom(lev).isPeriodic()));
    }

    for(int lev = 0; lev <= finest_level; ++lev)
    {
        // read in level 'lev' BoxArray from Header
        BoxArray ba;
        ba.readFrom(is);
        GotoNextLine(is);

        // Create distribution mapping
        DistributionMapping dm{ba, ParallelDescriptor::NProcs()};

        MakeNewLevelFromScratch(lev, cur_time, ba, dm);
    }

	/***************************************************************************
     * Load fluid data                                                         *
     ***************************************************************************/

	// Load the field data
	for(int lev = 0; lev <= finest_level; ++lev)
	{
	    // Read velocity and pressure gradients
	    MultiFab mf_vel;
	    VisMF::Read(mf_vel, MultiFabFileFullPrefix(lev, restart_file, level_prefix, "velx"));
            MultiFab::Copy((*vel[lev]), mf_vel, 0, 0, AMREX_SPACEDIM, 0);

	    MultiFab mf_gp;
	    VisMF::Read(mf_gp, MultiFabFileFullPrefix(lev, restart_file, level_prefix, "gpx"));
            MultiFab::Copy((*gp[lev]), mf_gp, 0, 0, AMREX_SPACEDIM, 0);

	    // Read scalar variables
	    for (int i = 0; i < chkscalarVars.size(); i++)
	    {
                // If we have created the walls using the domain boundary conditions and not
                //    by creating them from implicit functions, then the implicit_functions mf
                //    will be empty.  We don't want to fail when reading so we allow the code
                //    to read it in an empty multifab just for this one.
                int allow_empty_mf = 0;
                if (chkscaVarsName[i] == "implicit_functions") allow_empty_mf = 1;

                MultiFab mf;
                VisMF::Read(mf, amrex::MultiFabFileFullPrefix(lev, restart_file, level_prefix,
                                                              chkscaVarsName[i]), nullptr,
                                                              ParallelDescriptor::IOProcessorNumber(), 
                                                              allow_empty_mf);
                MultiFab::Copy(*((*chkscalarVars[i])[lev]), mf, 0, 0, 1, 0);
	    }
	}

	amrex::Print() << "Restart complete" << std::endl;
}

void incflo::WriteJobInfo(const std::string& path) const
{
	if(ParallelDescriptor::IOProcessor())
	{
		// job_info file with details about the run
		std::ofstream jobInfoFile;
		std::string FullPathJobInfoFile = path;
		std::string PrettyLine =
			"===============================================================================\n";

		FullPathJobInfoFile += "/incflo_job_info";
		jobInfoFile.open(FullPathJobInfoFile.c_str(), std::ios::out);

		// job information
		jobInfoFile << PrettyLine;
		jobInfoFile << " incflo Job Information\n";
		jobInfoFile << PrettyLine;

		jobInfoFile << "number of MPI processes: " << ParallelDescriptor::NProcs() << "\n";
#ifdef _OPENMP
		jobInfoFile << "number of threads:       " << omp_get_max_threads() << "\n";
#endif

		jobInfoFile << "\n\n";

		// build information
		jobInfoFile << PrettyLine;
		jobInfoFile << " Build Information\n";
		jobInfoFile << PrettyLine;

		jobInfoFile << "build date:    " << buildInfoGetBuildDate() << "\n";
		jobInfoFile << "build machine: " << buildInfoGetBuildMachine() << "\n";
		jobInfoFile << "build dir:     " << buildInfoGetBuildDir() << "\n";
		jobInfoFile << "AMReX dir:     " << buildInfoGetAMReXDir() << "\n";

		jobInfoFile << "\n";

		jobInfoFile << "COMP:          " << buildInfoGetComp() << "\n";
		jobInfoFile << "COMP version:  " << buildInfoGetCompVersion() << "\n";
		jobInfoFile << "FCOMP:         " << buildInfoGetFcomp() << "\n";
		jobInfoFile << "FCOMP version: " << buildInfoGetFcompVersion() << "\n";

		jobInfoFile << "\n";

		const char* githash1 = buildInfoGetGitHash(1);
		const char* githash2 = buildInfoGetGitHash(2);
		if(strlen(githash1) > 0)
		{
			jobInfoFile << "incflo git hash: " << githash1 << "\n";
		}
		if(strlen(githash2) > 0)
		{
			jobInfoFile << "AMReX git hash: " << githash2 << "\n";
		}

		jobInfoFile << "\n\n";

		// grid information
		jobInfoFile << PrettyLine;
		jobInfoFile << " Grid Information\n";
		jobInfoFile << PrettyLine;

		for(int i = 0; i <= finest_level; i++)
		{
			jobInfoFile << " level: " << i << "\n";
			jobInfoFile << "   number of boxes = " << grids[i].size() << "\n";
			jobInfoFile << "   maximum zones   = ";
			for(int dir = 0; dir < BL_SPACEDIM; dir++)
			{
				jobInfoFile << geom[i].Domain().length(dir) << " ";
			}
			jobInfoFile << "\n\n";
		}

		jobInfoFile << "\n\n";

		// runtime parameters
		jobInfoFile << PrettyLine;
		jobInfoFile << " Inputs File Parameters\n";
		jobInfoFile << PrettyLine;

		ParmParse::dumpTable(jobInfoFile, true);

		jobInfoFile.close();
	}
}

void incflo::WritePlotFile() 
{
	BL_PROFILE("incflo::WritePlotFile()");

    if (plt_divu      ) ComputeDivU(cur_time);
    if (plt_strainrate) ComputeStrainrate(cur_time);
    if (plt_eta) {
        ComputeViscosity(eta,cur_time);
        add_eddy_viscosity(eta,eta_tracer,cur_time);
    }
    if (plt_vort      ) ComputeVorticity(cur_time);

	const std::string& plotfilename = amrex::Concatenate(plot_file, nstep);

	amrex::Print() << "  Writing plotfile " << plotfilename << " at time " << cur_time << std::endl;

	const int ngrow = 0;

	Vector<std::unique_ptr<MultiFab>> mf(finest_level + 1);

        Vector<std::string> pltscaVarsName;

        // First just set all of the names once (not once per level)

        // Velocity components
        if (plt_velx == 1)
            pltscaVarsName.push_back("velx");
        if (plt_vely == 1)
            pltscaVarsName.push_back("vely");
        if (plt_velz == 1)
            pltscaVarsName.push_back("velz");

        // Pressure gradient components
        if (plt_gpx == 1)
            pltscaVarsName.push_back("gpx");
        if (plt_gpy == 1)
            pltscaVarsName.push_back("gpy");
        if (plt_gpz == 1)
            pltscaVarsName.push_back("gpz");

        // Density
        if (plt_rho == 1)
            pltscaVarsName.push_back("density");

        // Tracers
        if (plt_tracer == 1)
        {
            pltscaVarsName.push_back("tracer0");
            if (ntrac > 1)
               pltscaVarsName.push_back("tracer1");
            if (ntrac > 2)
               pltscaVarsName.push_back("tracer2");
            if (ntrac > 3)
               amrex::Error("Haven't made names for ntrac > 3 yet");
        }


        // Pressure
        if(plt_p == 1)
            pltscaVarsName.push_back("p");

        // Apparent viscosity
        if(plt_eta == 1) 
            pltscaVarsName.push_back("eta");

        // Vorticity
        if(plt_vort == 1) 
            pltscaVarsName.push_back("vort");

        // Magnitude of the rate-of-strain tensor 
        if(plt_strainrate == 1) 
            pltscaVarsName.push_back("strainrate");

        // Magnitude of the stress tensor 
        if(plt_stress == 1) 
            pltscaVarsName.push_back("stress");

        // Divergence of velocity field
        if(plt_divu == 1) 
            pltscaVarsName.push_back("divu");

        // Cut cell volume fraction
        if(plt_vfrac == 1) 
            pltscaVarsName.push_back("vfrac");

        // Now fill the data at every level

	for(int lev = 0; lev <= finest_level; ++lev)
	{
            // Multifab to hold all the variables -- there can be only one!!!!
    	    const int ncomp = pltVarCount;
#ifdef AMREX_USE_EB
	    mf[lev].reset(new MultiFab(grids[lev], dmap[lev], ncomp, ngrow, MFInfo(), *ebfactory[lev]));
#else
	    mf[lev].reset(new MultiFab(grids[lev], dmap[lev], ncomp, ngrow, MFInfo()));
#endif

            int lc = 0;

            // Velocity components
            if (plt_velx == 1)
            {
                MultiFab::Copy(*mf[lev], (*vel[lev]), 0, lc, 1, 0);
                lc += 1;
            }
            if (plt_vely == 1)
            {
                MultiFab::Copy(*mf[lev], (*vel[lev]), 1, lc, 1, 0);
                lc += 1;
            }
            if (plt_velz == 1)
            {
                MultiFab::Copy(*mf[lev], (*vel[lev]), 2, lc, 1, 0);
                lc += 1;
            }

            // Pressure gradient components
            if (plt_gpx == 1)
            {
                MultiFab::Copy(*mf[lev], (*gp[lev]), 0, lc, 1, 0);
                lc += 1;
            }
            if (plt_gpy == 1)
            {
                MultiFab::Copy(*mf[lev], (*gp[lev]), 1, lc, 1, 0);
                lc += 1;
            }
            if (plt_gpz == 1)
            {
                MultiFab::Copy(*mf[lev], (*gp[lev]), 2, lc, 1, 0);
                lc += 1;
            }

            // Density
            if (plt_rho == 1)
            {
                MultiFab::Copy(*mf[lev], (*density[lev]), 0, lc, 1, 0);
                lc += 1;
            }

            // Tracer
            if(plt_tracer == 1) 
            {
                MultiFab::Copy(*mf[lev], (*tracer[lev]), 0, lc, tracer[lev]->nComp(), 0);
                lc += tracer[lev]->nComp();
            }

            // Pressure
            if(plt_p == 1)
            {
                MultiFab p_nd(p[lev]->boxArray(), dmap[lev], 1, 0);
                p_nd.setVal(0.0);
                MultiFab::Copy(p_nd, (*p[lev]), 0, 0, 1, 0);
                MultiFab::Add(p_nd, (*p0[lev]), 0, 0, 1, 0);
                amrex::average_node_to_cellcenter(*mf[lev], lc, p_nd, 0, 1);
    
                lc += 1;
            }

            // Apparent viscosity
            if (plt_eta == 1)
            {
                MultiFab::Copy(*mf[lev], (*eta[lev]), 0, lc, 1, 0);
                lc += 1;
            }

            // Vorticity
            if(plt_vort == 1) 
            {
                MultiFab::Copy(*mf[lev], (*vort[lev]), 0, lc, 1, 0);
    
                lc += 1;
            }

            // Magnitude of the rate-of-strain tensor 
            if(plt_strainrate == 1) 
            {
                MultiFab::Copy(*mf[lev], (*strainrate[lev]), 0, lc, 1, 0);
    
                lc += 1;
            }

            // Magnitude of the stress tensor 
            if(plt_stress == 1) 
            {
                MultiFab::Copy(*mf[lev], (*strainrate[lev]), 0, lc, 1, 0);
                MultiFab::Multiply(*mf[lev], (*eta[lev]), 0, lc, 1, 0);

                lc += 1;
            }
    
            // Divergence of velocity field
            if(plt_divu == 1) 
            {
                amrex::average_node_to_cellcenter(*mf[lev], lc, (*divu[lev]), 0, 1);
                lc += 1;
            }

            // Cut cell volume fractino
            if(plt_vfrac == 1) 
            {
#ifdef AMREX_USE_EB
                if (ebfactory[lev]) 
                {
                    MultiFab::Copy(*mf[lev], ebfactory[lev]->getVolFrac(), 0, lc, 1, 0);
                }
                else
#endif
                {
                    mf[lev]->setVal(1.0, lc, 1.0);
                }

                lc += 1;
            }

#ifdef AMREX_USE_EB
            // Zero out all the values in covered cells
            EB_set_covered(*mf[lev], 0.0);
#endif
        }

        // This needs to be defined in order to use amrex::WriteMultiLevelPlotfile, 
        // but will never change unless we use subcycling. 
        // If we do use subcycling, this should be a incflo class member. 
        Vector<int> istep(finest_level + 1, 1);

        // Write the plotfile
        amrex::WriteMultiLevelPlotfile(plotfilename, finest_level + 1, GetVecOfConstPtrs(mf), 
                                   pltscaVarsName, Geom(), cur_time, istep, refRatio());

	WriteJobInfo(plotfilename);
}





void incflo::set_mfab_spatial_averaging_quantities(MultiFab &mfab, int lev, FArrayBox &avg_fab, int axis)
{

    if(axis!=2) amrex::Abort("not implemented for other index yet\n");
    
    AMREX_ASSERT(mfab->nComp() == fab->nComp());
    AMREX_ASSERT(mfab->nComp() == 19); // must match function that is calling this
    
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(mfab, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box bx = mfi.tilebox();

        const auto& mfab_arr = mfab.array(mfi);
        const auto& vel_arr = vel[lev]->array(mfi);
        const auto& tracer_arr = tracer[lev]->array(mfi);
        const auto& eta_arr = eta[lev]->array(mfi);
        const auto& den_arr = density[lev]->array(mfi);
        const auto& avg_fab_arr = avg_fab.array();

         // No cut cells in tile + 1-cell witdh halo -> use non-eb routine
        AMREX_FOR_3D(bx, i, j, k,
        {
            // velocities
            mfab_arr(i,j,k,0) = vel_arr(i,j,k,0);
            mfab_arr(i,j,k,1) = vel_arr(i,j,k,1);
            mfab_arr(i,j,k,2) = vel_arr(i,j,k,2);

            // potential temperature
            mfab_arr(i,j,k,3) = tracer_arr(i,j,k,0);
            
            // fluctuations
            const Real up = vel_arr(i,j,k,0) - avg_fab_arr(0,0,k,0);
            const Real vp = vel_arr(i,j,k,1) - avg_fab_arr(0,0,k,1);
            const Real wp = vel_arr(i,j,k,2) - avg_fab_arr(0,0,k,2);
            const Real Tp = tracer_arr(i,j,k,0) - avg_fab_arr(0,0,k,3);
            
            mfab_arr(i,j,k,4) = up*up;
            mfab_arr(i,j,k,5) = up*vp;
            mfab_arr(i,j,k,6) = up*wp;
            mfab_arr(i,j,k,7) = vp*vp;
            mfab_arr(i,j,k,8) = vp*wp;
            mfab_arr(i,j,k,9) = wp*wp;
            
            mfab_arr(i,j,k,10) = wp*up*up;
            mfab_arr(i,j,k,11) = wp*up*vp;
            mfab_arr(i,j,k,12) = wp*up*wp;
            mfab_arr(i,j,k,13) = wp*vp*vp;
            mfab_arr(i,j,k,14) = wp*vp*wp;
            mfab_arr(i,j,k,15) = wp*wp*wp;
            
            mfab_arr(i,j,k,16) = Tp*up;
            mfab_arr(i,j,k,17) = Tp*vp;
            mfab_arr(i,j,k,18) = Tp*wp;
            
            // nu+nut = (mu+mut)/rho
            mfab_arr(i,j,k,19) = eta_arr(i,j,k)/den_arr(i,j,k);

         });
    }
    
}



void incflo::spatially_average_quantities_down(bool plot)
{
    if(finest_level > 0) {
        amrex::Abort("commented out the fine to coarse averaging uncomment to test \n");
    }
    if(finest_level > 1) {
        amrex::Abort("average down only works for 2 levels for now need to add for loops \n");
    }

    const int ncomp = 20;// must match number of quantities in set_mfab_spatial_averaging_quantities(...)
    const int ngrow = 0;

    int crse_lev = 0;
    int fine_lev = finest_level;

    // create a 3D box with 1D indices
    const Box& domain = geom[0].Domain();
    IntVect dom_lo(domain.loVect());
    IntVect dom_hi(domain.hiVect());
    
    
    const int axis = 2;// fixme if you change this the loop below will break

    
    const int s = domain.smallEnd(axis);
    const int b = domain.bigEnd(axis);

    IntVect small(0,0,0);
    IntVect big(0,0,0);

    small.setVal(axis,s);
    big.setVal(axis,b);

    // have every processor make a 1d box
    Box bx1d = Box(small,big);

    // have every processor make a 1d fab
    FArrayBox fab(bx1d,ncomp);
    fab.setVal(0.0);

    // count number of cells in plane
    Real nxny = 1.0;
    for(int i=0;i<AMREX_SPACEDIM;++i){
        if(i!=axis) nxny *= (dom_hi[i]-dom_lo[i]+1);
    }
         
    
     // need to put this in to a for loop over levels starting from finest and working down
    MultiFab crse_mfab(grids[crse_lev], dmap[crse_lev], ncomp, ngrow);
//    MultiFab fine_mfab(grids[fine_lev], dmap[fine_lev], ncomp, ngrow);

    set_mfab_spatial_averaging_quantities(crse_mfab,crse_lev,fab,axis);
//    set_mfab_spatial_averaging_quantities(fine_mfab,fine_lev,fab,axis);
    
//    const int ratio = 2;
//    average_down(fine_mfab, crse_mfab, geom[fine_lev], geom[crse_lev], 0, ncomp, ratio);

    
    // average the fab
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(crse_mfab, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box bx = mfi.tilebox();

        const auto& crse_arr = crse_mfab.array(mfi);
        const auto& fab_arr = fab.array();
        // ncomp is used here but it could actually be 4 depending on the average down function
        int nvar = 4;
        switch (axis) {
            case 0:
                AMREX_FOR_4D(bx, nvar, i, j, k, n, {fab_arr(i,0,0,n) += crse_arr(i,j,k,n)/nxny;} );
                break;
            case 1:
                AMREX_FOR_4D(bx, nvar, i, j, k, n, {fab_arr(0,j,0,n) += crse_arr(i,j,k,n)/nxny;} );
                break;
            case 2:
                AMREX_FOR_4D(bx, nvar, i, j, k, n, {fab_arr(0,0,k,n) += crse_arr(i,j,k,n)/nxny;} );
                break;
            default:
                amrex::Abort("index value in average quantities is garbage \n");
                break;
        }
    }

    
    //fixme put in a loop like above
    // sum all fabs together
    ParallelDescriptor::ReduceRealSum(fab.dataPtr(0),fab.size());

    set_mfab_spatial_averaging_quantities(crse_mfab,crse_lev,fab,axis);
//    set_mfab_spatial_averaging_quantities(fine_mfab,fine_lev,fab,axis);

//    average_down(fine_mfab, crse_mfab, geom[fine_lev], geom[crse_lev], 0, ncomp, ratio);


    fab.setVal(0.0);
    
    // average the fab again now that velocity and temperature are averaged!
    #ifdef _OPENMP
    #pragma omp parallel if (Gpu::notInLaunchRegion())
    #endif
        for (MFIter mfi(crse_mfab, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Box bx = mfi.tilebox();

            const auto& crse_arr = crse_mfab.array(mfi);
            const auto& fab_arr = fab.array();

            switch (axis) {
                case 0:
                    AMREX_FOR_4D(bx, ncomp, i, j, k, n, {fab_arr(i,0,0,n) += crse_arr(i,j,k,n)/nxny;} );
                    break;
                case 1:
                    AMREX_FOR_4D(bx, ncomp, i, j, k, n, {fab_arr(0,j,0,n) += crse_arr(i,j,k,n)/nxny;} );
                    break;
                case 2:
                    AMREX_FOR_4D(bx, ncomp, i, j, k, n, {fab_arr(0,0,k,n) += crse_arr(i,j,k,n)/nxny;} );
                    break;
                default:
                    amrex::Abort("index value in average quantities is garbage \n");
                    break;
            }
        }
    
    // sum all fabs together
    ParallelDescriptor::ReduceRealSum(fab.dataPtr(0),fab.size());
  
    
    // fixme piggy backing on this function but this belongs somewhere else
    // maybe keep this 1d average in global memory and then make these functions separate?
    AMREX_ASSERT(axis==2);
    const auto& fab_arr = fab.array();
    
    for(int i=s; i <= b; ++i){
        const Real z = geom[0].ProbLo(axis) + (i+0.5)*geom[0].CellSize(axis);
        if(z > 0.0) {
            vx_mean_ground = fab_arr(0,0,i,0);
            vy_mean_ground = fab_arr(0,0,i,1);
            nu_mean_ground = fab_arr(0,0,i,19);
            z_ground = z;
           
            amrex::Print() << "z: " << z_ground << " vx_mean_ground: " << vx_mean_ground << " vy_mean_ground: " << vy_mean_ground << std::endl;
            
            // fixme circular dependency so need to hack this for now
            if(nstep == -1) nu_mean_ground = 1.0;
            amrex::Print() << "nu mean ground: " << nu_mean_ground << std::endl;

            break;
        }
    }
    
    // simple shear stress model for neutral BL
    // apply as an inhomogeneous Neumann BC
    const Real uh = sqrt(pow(vx_mean_ground,2) + pow(vy_mean_ground,2)) + 1.0e-12;
    utau = kappa*uh/log10(z_ground/surface_roughness_z0);
    amrex::Print() << "utau: " << utau << std::endl;
    
    for(int i=s; i <= b; ++i){
        const Real z = geom[0].ProbLo(axis) + (i+0.5)*geom[0].CellSize(axis);
        if(z > abl_forcing_height) {
            vx_mean = fab_arr(0,0,i,0);
            vy_mean = fab_arr(0,0,i,1);
            amrex::Print() << "z: " << z << " vx_mean: " << vx_mean << " vy_mean: " << vy_mean << std::endl;
            break;
        }
    }
    
    
    
    
    if(!plot) return;
    
    // begin collecting and output values
    
     // make a box array with a single box
    BoxArray ba1d(bx1d);

    DistributionMapping dm1d;// {ba1d}; // could use the boxarray but not sure if always guaranteed that proc=0 owns it
    // only processor 0 owns the box
    Vector<int> pmap {0};
    dm1d.define(pmap);

     // create a multiFAB using one box on one proc
    MultiFab mfab1d(ba1d,dm1d,ncomp,ngrow);


     // fixme there has to be a way to copy a fab into a mfab right??
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(mfab1d, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box bx = mfi.tilebox();

        const auto& fab_arr1 = mfab1d.array(mfi);
        const auto& fab_arr2 = fab.array();

        AMREX_FOR_4D(bx, ncomp, i, j, k, n,
        {
            fab_arr1(0,0,k,n) = fab_arr2(0,0,k,n);
        });
    }
//
//    if(ParallelDescriptor::IOProcessorNumber() == 0){
//        FArrayBox& fabby = mfab[0];
//        Real *fptr1 = fabby.dataPtr(0);
//        Real *fptr2 = fab.dataPtr(0);
//        for(int i=0; i < fab.size(); ++i){
//            fptr1[i] = fptr2[i];
//        }
//    }

    std::string line_plot_file{"line_plot"};

    const std::string& plotfilename = amrex::Concatenate(line_plot_file, nstep);

    amrex::Print() << "  Writing line plot file " << plotfilename << " at time " << cur_time << std::endl;

    Vector<std::string> pltscaVarsName;

    pltscaVarsName.push_back("<u>");
    pltscaVarsName.push_back("<v>");
    pltscaVarsName.push_back("<w>");
    pltscaVarsName.push_back("<T>");
    
    pltscaVarsName.push_back("<u'u'>");
    pltscaVarsName.push_back("<u'v'>");
    pltscaVarsName.push_back("<u'w'>");
    pltscaVarsName.push_back("<v'v'>");
    pltscaVarsName.push_back("<v'w'>");
    pltscaVarsName.push_back("<w'w'>");
    
    pltscaVarsName.push_back("<w'u'u'>");
    pltscaVarsName.push_back("<w'u'v'>");
    pltscaVarsName.push_back("<w'u'w'>");
    pltscaVarsName.push_back("<w'v'v'>");
    pltscaVarsName.push_back("<w'v'w'>");
    pltscaVarsName.push_back("<w'w'w'>");
    
    pltscaVarsName.push_back("<T'u'>");
    pltscaVarsName.push_back("<T'v'>");
    pltscaVarsName.push_back("<T'w'>");
    
    //fixme todo add Rij, qj, nu_SGS https://a2ehfm-ecp.slack.com/archives/C3V26K34G/p1522245519000450
    
    int level_step = 0;
    amrex::WriteSingleLevelPlotfile(plotfilename, mfab1d, pltscaVarsName, geom[0], cur_time, level_step);
 
    WriteJobInfo(plotfilename);


 }
