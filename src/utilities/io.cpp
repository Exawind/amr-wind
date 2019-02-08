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
			HeaderFile << Geometry::ProbLo(i) << ' ';
		HeaderFile << '\n';

		for(int i = 0; i < BL_SPACEDIM; ++i)
			HeaderFile << Geometry::ProbHi(i) << ' ';
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
    Geometry::ProbDomain(RealBox(prob_lo, prob_hi));

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
        vel[lev]->copy(mf_vel, 0, 0, 3, 0, 0);

		MultiFab mf_gp;
		VisMF::Read(mf_gp, MultiFabFileFullPrefix(lev, restart_file, level_prefix, "gpx"));
        gp[lev]->copy(mf_gp, 0, 0, 3, 0, 0);

		// Read scalar variables
		for(int i = 0; i < chkscalarVars.size(); i++)
		{
			MultiFab mf;
            VisMF::Read(mf, MultiFabFileFullPrefix(lev, restart_file, 
                                                   level_prefix, chkscaVarsName[i]));
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

	for(int lev = 0; lev <= finest_level; ++lev)
	{
		// the "+1" here is for volfrac
		const int ncomp = vecVarsName.size() + pltscalarVars.size() + 1;
		mf[lev].reset(new MultiFab(grids[lev], dmap[lev], ncomp, ngrow, MFInfo(), *ebfactory[lev]));

        for(int i = 0; i < 3; i++)
        {
            // Velocity components
            MultiFab::Copy(*mf[lev], (*vel[lev]), i, i, 1, 0);

            // Pressure gradient components
            MultiFab::Copy(*mf[lev], (*gp[lev]), i, i+3, 1, 0);
            
            // Multiply by volume fraction to get smoother looks in EB cells
            // (but otherwise, the results are more accurate inside the contour vfrac=0.5)
            // if(ebfactory[lev])
            // {
            //     MultiFab::Multiply(*mf[lev], ebfactory[lev]->getVolFrac(), 0, i, 1, 0);
            //     MultiFab::Multiply(*mf[lev], ebfactory[lev]->getVolFrac(), 0, i+3, 1, 0);
            // }
        }

		// Scalar variables
		int dcomp = vecVarsName.size();
		for(int i = 0; i < pltscalarVars.size(); i++)
		{
			if(pltscaVarsName[i] == "p")
			{
                MultiFab p_nd(p[lev]->boxArray(), dmap[lev], 1, 0);
                p_nd.setVal(0.0);
                MultiFab::Copy(p_nd, (*p[lev]), 0, 0, 1, 0);
                MultiFab::Add(p_nd, (*p0[lev]), 0, 0, 1, 0);
                amrex::average_node_to_cellcenter(*mf[lev], dcomp, p_nd, 0, 1);
			}
			else if(pltscaVarsName[i] == "divu")
			{
				amrex::average_node_to_cellcenter(
					*mf[lev], dcomp, *(*pltscalarVars[i])[lev].get(), 0, 1);
			}
			else if(pltscaVarsName[i] == "strainrate")
			{
				MultiFab::Copy(*mf[lev], (*strainrate[lev]), 0, dcomp, 1, 0);
			}
			else if(pltscaVarsName[i] == "stress")
			{
				MultiFab::Copy(*mf[lev], (*strainrate[lev]), 0, dcomp, 1, 0);
				MultiFab::Multiply(*mf[lev], (*eta[lev]), 0, dcomp, 1, 0);
			}
			else if(pltscaVarsName[i] == "vort")
			{
				MultiFab::Copy(*mf[lev], (*vort[lev]), 0, dcomp, 1, 0);
			}
			else
			{
				MultiFab::Copy(*mf[lev], *((*pltscalarVars[i])[lev].get()), 0, dcomp, 1, 0);
			}

			dcomp++;
		}

		if(ebfactory[lev])
		{
			MultiFab::Copy(*mf[lev], ebfactory[lev]->getVolFrac(), 0, dcomp, 1, 0);
		}
		else
		{
			mf[lev]->setVal(1.0, dcomp, 1, 0);
		}

        EB_set_covered(*mf[lev], 0.0);
	}

	// Concatenate scalar and vector var names
	Vector<std::string> names;
	names.insert(names.end(), vecVarsName.begin(), vecVarsName.end());
	names.insert(names.end(), pltscaVarsName.begin(), pltscaVarsName.end());

    // This needs to be defined in order to use amrex::WriteMultiLevelPlotfile, 
    // but will never change unless we use subcycling. 
    // If we do use subcycling, this should be a incflo class member. 
    Vector<int> istep(finest_level + 1, 1);

    // Write the plotfile
    amrex::WriteMultiLevelPlotfile(plotfilename, finest_level + 1, GetVecOfConstPtrs(mf), names, 
                                   Geom(), cur_time, istep, refRatio());

	WriteJobInfo(plotfilename);
}
