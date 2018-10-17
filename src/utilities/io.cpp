#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>

#include <AMReX_AmrCore.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_VectorIO.H> // amrex::[read,write]IntData(array_of_ints)
#include <AMReX_VisMF.H> // amrex::VisMF::Write(MultiFab)

#include <AMReX_buildInfo.H>

#include <incflo.H>

namespace
{
const std::string level_prefix{"Level_"};
}

// This function initializes the attributes vecVarsName,
//                                          pltscalarVars, pltscaVarsName,
//                                          chkscalarVars, chkscaVarsName.
// If new variables need to be added to the output/checkpoint, simply add them
// here and the IO routines will automatically take care of them.
void incflo::InitIOData()
{
	// Define the list of vector variables on faces that need to be written
	// to plotfile/checkfile.
	vecVarsName = {"velx", "vely", "velz", "gpx", "gpy", "gpz"};

	// Define the list of scalar variables at cell centers that need to be
	// written to plotfile/checkfile. "volfrac" MUST always be last without any
	// mf associated to it!!!
	pltscaVarsName = {"p", "ro", "eta", "strainrate", "stress", "vort", "divu", "volfrac"};
	pltscalarVars = {&p, &ro, &eta, &strainrate, &strainrate, &vort, &divu};

	chkscaVarsName = {"p", "ro", "eta"};
	chkscalarVars = {&p, &ro, &eta};
}

void incflo::WritePlotHeader(const std::string& name, int nstep, Real dt, Real time) const
{
	bool is_checkpoint = 0;
	WriteHeader(name, nstep, dt, time, is_checkpoint);
}

void incflo::WriteCheckHeader(const std::string& name, int nstep, Real dt, Real time) const
{
	bool is_checkpoint = 1;
	WriteHeader(name, nstep, dt, time, is_checkpoint);
}

void incflo::WriteHeader(
	const std::string& name, int nstep, Real dt, Real time, bool is_checkpoint) const
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

		const int nlevels = finestLevel() + 1;
		HeaderFile << nlevels << "\n";

		// Time stepping controls
		HeaderFile << nstep << "\n";
		HeaderFile << dt << "\n";
		HeaderFile << time << "\n";

		// Geometry
		for(int i = 0; i < BL_SPACEDIM; ++i)
			HeaderFile << Geometry::ProbLo(i) << ' ';
		HeaderFile << '\n';

		for(int i = 0; i < BL_SPACEDIM; ++i)
			HeaderFile << Geometry::ProbHi(i) << ' ';
		HeaderFile << '\n';

		// BoxArray
		for(int lev = 0; lev < nlevels; ++lev)
		{
			boxArray(lev).writeOn(HeaderFile);
			HeaderFile << '\n';
		}
	}
}

void incflo::WriteCheckPointFile(std::string& check_file, int nstep, Real dt, Real time) const
{
	BL_PROFILE("incflo::WriteCheckPointFile()");

	const std::string& checkpointname = amrex::Concatenate(check_file, nstep);

    amrex::Print() << "\n\t Writing checkpoint " << checkpointname << std::endl;

	const int nlevels = finestLevel() + 1;
	amrex::PreBuildDirectorHierarchy(checkpointname, level_prefix, nlevels, true);

	WriteCheckHeader(checkpointname, nstep, dt, time);

	WriteJobInfo(checkpointname);

	for(int lev = 0; lev < nlevels; ++lev)
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

void incflo::Restart(
	std::string& restart_file, int* nstep, Real* dt, Real* time, IntVect& Nrep)
{
	BL_PROFILE("incflo::Restart()");

	amrex::Print() << "  Restarting from checkpoint " << restart_file << std::endl;

	Real prob_lo[BL_SPACEDIM];
	Real prob_hi[BL_SPACEDIM];

	/***************************************************************************
     * Load header: set up problem domain (including BoxArray)                 *
     *              allocate incflo memory (incflo::AllocateArrays)    *
     ***************************************************************************/

	{
		std::string File(restart_file + "/Header");

		VisMF::IO_Buffer io_buffer(VisMF::GetIOBufferSize());

		Vector<char> fileCharPtr;
		ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
		std::string fileCharPtrString(fileCharPtr.dataPtr());
		std::istringstream is(fileCharPtrString, std::istringstream::in);

		std::string line, word;

		std::getline(is, line);

		int nlevs;
		int int_tmp;
		Real real_tmp;

		is >> nlevs;
		GotoNextLine(is);
		// finest_level = nlevs-1;

		// Time stepping controls
		is >> int_tmp;
		*nstep = int_tmp;
		GotoNextLine(is);

		is >> real_tmp;
		*dt = real_tmp;
		GotoNextLine(is);

		is >> real_tmp;
		*time = real_tmp;
		GotoNextLine(is);

		std::getline(is, line);
		{
			std::istringstream lis(line);
			int i = 0;
			while(lis >> word)
			{
				prob_lo[i++] = std::stod(word);
			}
		}

		std::getline(is, line);
		{
			std::istringstream lis(line);
			int i = 0;
			while(lis >> word)
			{
				prob_hi[i++] = std::stod(word);
			}
		}

		if(Nrep != IntVect::TheUnitVector())
		{
			for(int d = 0; d < BL_SPACEDIM; d++)
			{
				prob_lo[d] = Nrep[d] * prob_lo[d];
				prob_hi[d] = Nrep[d] * prob_hi[d];
			}
		}
		Geometry::ProbDomain(RealBox(prob_lo, prob_hi));

		for(int lev = 0; lev < nlevs; ++lev)
		{
			BoxArray orig_ba, ba;
			orig_ba.readFrom(is);
			GotoNextLine(is);

			SetBoxArray(lev, orig_ba);
			DistributionMapping orig_dm{orig_ba, ParallelDescriptor::NProcs()};
			SetDistributionMap(lev, orig_dm);

			Box orig_domain(orig_ba.minimalBox());

			if(Nrep != IntVect::TheUnitVector())
			{
				amrex::Print() << " OLD BA had " << orig_ba.size() << " GRIDS " << std::endl;
				amrex::Print() << " OLD Domain" << orig_domain << std::endl;
			}

			BoxList bl;
			for(int nb = 0; nb < orig_ba.size(); nb++)
			{
				for(int k = 0; k < Nrep[2]; k++)
				{
					for(int j = 0; j < Nrep[1]; j++)
					{
						for(int i = 0; i < Nrep[0]; i++)
						{
							Box b(orig_ba[nb]);
							IntVect lo(b.smallEnd());
							IntVect hi(b.bigEnd());
							IntVect shift_vec(i * orig_domain.length(0),
											  j * orig_domain.length(1),
											  k * orig_domain.length(2));
							b.shift(shift_vec);
							bl.push_back(b);
						}
					}
				}
			}
			ba.define(bl);

			if(Nrep != IntVect::TheUnitVector())
			{
				Box new_domain(ba.minimalBox());
				geom[lev].Domain(new_domain);

				amrex::Print() << " NEW BA has " << ba.size() << " GRIDS " << std::endl;
				amrex::Print() << " NEW Domain" << geom[0].Domain() << std::endl;

				DistributionMapping dm{ba, ParallelDescriptor::NProcs()};
				ReMakeNewLevelFromScratch(lev, ba, dm);
			}

			// This needs is needed before initializing level MultiFabs: ebfactories should
			// not change after the eb-dependent MultiFabs are allocated.
			make_eb_geometry();

			// Allocate the fluid data, NOTE: this depends on the ebfactories.
			AllocateArrays(lev);
		}
	}

	amrex::Print() << "  Finished reading header" << std::endl;

	/***************************************************************************
     * Load fluid data                                                         *
     ***************************************************************************/

	// Load the field data
	for(int lev = 0, nlevs = finestLevel() + 1; lev < nlevs; ++lev)
	{
		// Read velocity and pressure gradients
		MultiFab mf_vel;
		VisMF::Read(mf_vel, amrex::MultiFabFileFullPrefix(lev, restart_file, level_prefix, "velx"));

		MultiFab mf_gp;
		VisMF::Read(mf_gp, amrex::MultiFabFileFullPrefix(lev, restart_file, level_prefix, "gpx"));

		if(Nrep == IntVect::TheUnitVector())
		{
			// Simply copy mf_vel into vel, mf_gp intp gp
			vel[lev]->copy(mf_vel, 0, 0, 3, 0, 0);
			gp[lev]->copy(mf_gp, 0, 0, 3, 0, 0);
		}
		else
		{

			if(mf_vel.boxArray().size() > 1)
				amrex::Abort("Replication only works if one initial grid");

			mf_vel.FillBoundary(geom[lev].periodicity());
			mf_gp.FillBoundary(geom[lev].periodicity());

			FArrayBox single_fab_vel(mf_vel.boxArray()[0], 3);
			FArrayBox single_fab_gp(mf_gp.boxArray()[0], 3);

			// Copy and replicate mf into velocity
			for(MFIter mfi(*vel[lev]); mfi.isValid(); ++mfi)
			{
				int ib = mfi.index();
				(*vel[lev])[ib].copy(single_fab_vel, single_fab_vel.box(), 0, mfi.validbox(), 0, 3);
				(*gp[lev])[ib].copy(single_fab_gp, single_fab_gp.box(), 0, mfi.validbox(), 0, 3);
			}
		}

		// Read scalar variables
		for(int i = 0; i < chkscalarVars.size(); i++)
		{
			MultiFab mf;
			VisMF::Read(
				mf,
				amrex::MultiFabFileFullPrefix(lev, restart_file, level_prefix, chkscaVarsName[i]));

			if(Nrep == IntVect::TheUnitVector())
			{
				amrex::Print() << "  - loading scalar data: " << chkscaVarsName[i] << std::endl;

				// Copy mf into chkscalarVars
				if(chkscaVarsName[i] == "level-set")
				{
					// The level-set data is special, and because we want
					// to access it even without a fluid present, it is
					// loaded below.
				}
				else
				{
					(*chkscalarVars[i])[lev]->copy(mf, 0, 0, 1, 0, 0);
				}
			}
			else
			{

				if(mf.boxArray().size() > 1)
					amrex::Abort("Replication only works if one initial grid");

				mf.FillBoundary(geom[lev].periodicity());

				FArrayBox single_fab(mf.boxArray()[0], 1);
				mf.copyTo(single_fab);

				// Copy and replicate mf into chkscalarVars
				for(MFIter mfi(*(*chkscalarVars[i])[lev]); mfi.isValid(); ++mfi)
				{
					int ib = mfi.index();
					(*(*chkscalarVars[i])[lev])[ib].copy(
						single_fab, single_fab.box(), 0, mfi.validbox(), 0, 1);
				}
			}
		}
	}
	amrex::Print() << "  Finished reading fluid data" << std::endl;

	for(int lev = 0, nlevs = finestLevel() + 1; lev < nlevs; ++lev)
	{
        fill_mf_bc(lev, *p[lev]);
        fill_mf_bc(lev, *p_o[lev]);

		fill_mf_bc(lev, *ro[lev]);
		fill_mf_bc(lev, *ro_o[lev]);

		fill_mf_bc(lev, *eta[lev]);

		// Fill the bc's just in case
		vel[lev]->FillBoundary(geom[lev].periodicity());
		vel_o[lev]->FillBoundary(geom[lev].periodicity());

		// used in load balancing
		if(load_balance_type == "KnapSack")
		{
			fluid_cost[lev].reset(new MultiFab(grids[lev], dmap[lev], 1, 0));
			fluid_cost[lev]->setVal(0.0);
		}
	}
	amrex::Print() << "  Done with incflo::Restart " << std::endl;
}

void incflo::GotoNextLine(std::istream& is)
{
	constexpr std::streamsize bl_ignore_max{100000};
	is.ignore(bl_ignore_max, '\n');
}

void incflo::WriteJobInfo(const std::string& dir) const
{
	if(ParallelDescriptor::IOProcessor())
	{
		// job_info file with details about the run
		std::ofstream jobInfoFile;
		std::string FullPathJobInfoFile = dir;
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
			for(int n = 0; n < BL_SPACEDIM; n++)
			{
				jobInfoFile << geom[i].Domain().length(n) << " ";
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

void incflo::WritePlotFile(std::string& plot_file, int nstep, Real dt, Real time) const
{
	BL_PROFILE("incflo::WritePlotFile()");

	const std::string& plotfilename = amrex::Concatenate(plot_file, nstep);

	amrex::Print() << "  Writing plotfile " << plotfilename << std::endl;

	const int ngrow = 0;

	Vector<std::unique_ptr<MultiFab>> mf(nlev);

	for(int lev = 0; lev < nlev; ++lev)
	{
		// the "+1" here is for volfrac
		const int ncomp = vecVarsName.size() + pltscalarVars.size() + 1;
		mf[lev].reset(new MultiFab(grids[lev], dmap[lev], ncomp, ngrow));

        for(int i = 0; i < 3; i++)
        {
            // Velocity components
            MultiFab::Copy(*mf[lev], (*vel[lev]), i, i, 1, 0);

            // Pressure gradient components
            MultiFab::Copy(*mf[lev], (*gp[lev]), i, i+3, 1, 0);
            
            // Multiply by volume fraction to get proper results in EB cells
            if(ebfactory[lev])
            {
                MultiFab::Multiply(*mf[lev], ebfactory[lev]->getVolFrac(), 0, i, 1, 0);
                MultiFab::Multiply(*mf[lev], ebfactory[lev]->getVolFrac(), 0, i+3, 1, 0);
            }
        }

		// Scalar variables
		int dcomp = vecVarsName.size();
		for(int i = 0; i < pltscalarVars.size(); i++)
		{
			if(pltscaVarsName[i] == "p")
			{
                MultiFab::Copy(*mf[lev], (*p[lev]), 0, dcomp, 1, 0);
                MultiFab::Add(*mf[lev], (*p0[lev]), 0, dcomp, 1, 0);
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
         // Multiply by volume fraction to get proper results in EB cells
         if(ebfactory[lev])
         {
             MultiFab::Multiply(*mf[lev], ebfactory[lev]->getVolFrac(), 0, dcomp, 1, 0);
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
	}

	Vector<const MultiFab*> mf2(nlev);

	for(int lev = 0; lev < nlev; ++lev)
	{
		mf2[lev] = mf[lev].get();
	}

	// Concatenate scalar and vector var names
	Vector<std::string> names;
	names.insert(names.end(), vecVarsName.begin(), vecVarsName.end());
	names.insert(names.end(), pltscaVarsName.begin(), pltscaVarsName.end());

    amrex::WriteMultiLevelPlotfile(plotfilename, nlev, mf2, names, Geom(), time, istep, refRatio());

	WriteJobInfo(plotfilename);
}
