#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_MultiFabUtil.H>
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
		VisMF::Write(
			(*vel[lev]),
			amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, vecVarsName[0]));

		// This writes all three pressure gradient components
		VisMF::Write(
			(*gp[lev]),
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
        vel[lev]->copy(mf_vel, 0, 0, AMREX_SPACEDIM, 0, 0);

		MultiFab mf_gp;
		VisMF::Read(mf_gp, MultiFabFileFullPrefix(lev, restart_file, level_prefix, "gpx"));
        gp[lev]->copy(mf_gp, 0, 0, AMREX_SPACEDIM, 0, 0);

		// Read scalar variables
		for(int i = 0; i < chkscalarVars.size(); i++)
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
                        ParallelDescriptor::IOProcessorNumber(), allow_empty_mf);

            (*chkscalarVars[i])[lev]->copy(mf, 0, 0, 1, 0, 0);
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

void incflo::WritePlotFile() const
{

	BL_PROFILE("incflo::WritePlotFile()");

	const std::string& plotfilename = amrex::Concatenate(plot_file, nstep);

	amrex::Print() << "  Writing plotfile " << plotfilename << std::endl;

	const int ngrow = 0;

	Vector<std::unique_ptr<MultiFab>> mf(finest_level + 1);

        Vector<std::string> pltscaVarsName;

        // First just set all of the names once (not once per level)

        // Velocity components
        if (plt_vel == 1)
        {
            pltscaVarsName.push_back("velx");
            pltscaVarsName.push_back("vely");
            pltscaVarsName.push_back("velz");
        }

        // Pressure gradient components
        if (plt_gradp == 1)
        {
            pltscaVarsName.push_back("gpx");
            pltscaVarsName.push_back("gpy");
            pltscaVarsName.push_back("gpz");
        }

        // Density
        if(plt_rho == 1) 
            pltscaVarsName.push_back("ro");

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
	    mf[lev].reset(new MultiFab(grids[lev], dmap[lev], ncomp, ngrow, MFInfo(), *ebfactory[lev]));

            int lc = 0;

            // Velocity components
            if(plt_vel == 1)
            {
                MultiFab::Copy(*mf[lev], (*vel[lev]), 0, lc  , 1, 0);
                MultiFab::Copy(*mf[lev], (*vel[lev]), 1, lc+1, 1, 0);
                MultiFab::Copy(*mf[lev], (*vel[lev]), 2, lc+2, 1, 0);

                lc += AMREX_SPACEDIM;
            }

            // Pressure gradient components
            if(plt_gradp == 1)
            {
                MultiFab::Copy(*mf[lev], (*gp[lev]), 0, lc  , 1, 0);
                MultiFab::Copy(*mf[lev], (*gp[lev]), 1, lc+1, 1, 0);
                MultiFab::Copy(*mf[lev], (*gp[lev]), 2, lc+2, 1, 0);
    
                lc += AMREX_SPACEDIM;
            }

            // Density
            if(plt_rho == 1) 
            {
                MultiFab::Copy(*mf[lev], (*ro[lev]), 0, lc, 1, 0);
    
                lc += 1;
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
            if(plt_eta == 1) 
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
                MultiFab::Copy(*mf[lev], (*divu[lev]), 0, lc, 1, 0);
    
                lc += 1;
            }

            // Cut cell volume fractino
            if(plt_vfrac == 1) 
            {
                if (ebfactory[lev]) 
                {
                    MultiFab::Copy(*mf[lev], ebfactory[lev]->getVolFrac(), 0, lc, 1, 0);
                }
                else
                {
                    mf[lev]->setVal(1.0, lc, 1.0);
                }

                lc += 1;
            }

            // Zero out all the values in covered cells
            EB_set_covered(*mf[lev], 0.0);
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
