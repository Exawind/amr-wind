#include <fstream>
#include <iomanip>

#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>
#include <AMReX_VisMF.H>
#include <AMReX_iMultiFab.H>

#include <incflo_level.H>
#include <derive_F.H>
#include <io_F.H>
#include <setup_F.H>

// Declare and initialise variables
int verbose = -1;

Real stop_time = -1.0;
int max_step = -1;
bool steady_state = false;

int check_int = -1;
int last_chk = -1;
std::string check_file{"chk"};
std::string restart_file{""};

int plot_int = -1;
int last_plt = -1;
std::string plot_file{"plt"};
bool write_eb_surface = false;

int repl_x = 1;
int repl_y = 1;
int repl_z = 1;

int regrid_int = -1;

void ReadParameters();


int main(int argc, char* argv[])
{
	// AMReX will now read the inputs file and the command line arguments, but the
	//        command line arguments are in incflo-format so it will just ignore them.
	amrex::Initialize(argc, argv);
	BL_PROFILE_VAR("main()", pmain)
	BL_PROFILE_REGION_START("incflo::main()");

	// Issue an error if AMR input file is not given
	if(argc < 2)
		amrex::Abort("AMReX input file missing");

	// Setting format to NATIVE rather than default of NATIVE_32
	FArrayBox::setFormat(FABio::FAB_NATIVE);

	// Copy arguments into incflo -- note that the first argument is now the name of the
	//      inputs file to be read by AMReX, so we only pass the arguments after that
	for(int i = 2; i < argc; i++)
	{
		int nlen = strlen(argv[i]);

		// If-statement avoids passing the name of the incflo input file if it is
		// specified on the command line or any AMReX command.
		if((strstr(argv[i], "input_file") == NULL) && (strstr(argv[i], "amr") == NULL) &&
		   (strstr(argv[i], "incflo") == NULL))
			incflo_add_argument(argv[i], &nlen);
	}

    // Start timing the initialisation process
	Real strt_time = ParallelDescriptor::second();

    // Read parameters from ParmParse
	ReadParameters();

    // Time and time step counters
	Real time = 0.0L;
	int nstep = 0; 

	// Default AMR level = 0
	int lev = 0;

	// Define dt here in main (dt is not a member of incflo_level)
	Real dt = -1;

	// Default constructor. Note inheritance: incflo_level : AmrCore : AmrMesh
	//                                                                  |
	//  => Geometry is constructed here:  (constructs Geometry) --------+
	incflo_level my_incflo;

	// Initialize internals from ParmParse database
	my_incflo.InitParams();

	// Initialize memory for data-array internals
	my_incflo.ResizeArrays();

	// Initialize derived internals
	my_incflo.Init(lev, time);

	// Either init from scratch or from the checkpoint file
	int restart_flag = 0;
	if(restart_file.empty())
	{
		// NOTE: this also builds ebfactories and level-set
		my_incflo.InitLevelData(lev, time);
	}
	else
	{
		restart_flag = 1;

		// NOTE: 1) this also builds ebfactories and level-set 
        //       2) this can change the grids (during replication)
		IntVect Nrep(repl_x, repl_y, repl_z);
		my_incflo.Restart(restart_file, &nstep, &dt, &time, Nrep);
	}


	// Regrid
	my_incflo.Regrid(lev);

    // Post-initialisation step
	my_incflo.PostInit(lev, dt, time, nstep, restart_flag, stop_time, steady_state);

	// Write out EB sruface
	if(write_eb_surface)
		my_incflo.WriteEBSurface(lev);

	Real end_init = ParallelDescriptor::second() - strt_time;
	ParallelDescriptor::ReduceRealMax(end_init, ParallelDescriptor::IOProcessorNumber());

    amrex::Print() << "Time spent in init      " << end_init << std::endl;

	int finish = 0;

	// Initialize prev_dt here; it will be re-defined by call to evolve_fluid
	Real prev_dt = dt;

	// We automatically write checkpoint and plotfiles with the initial data
	//    if plot_int > 0
	if(restart_file.empty() && plot_int > 0)
	{
		my_incflo.incflo_compute_vort(lev);
		my_incflo.WritePlotFile(plot_file, nstep, dt, time);
	}

	// We automatically write checkpoint files with the initial data
	//    if check_int > 0
	if(restart_file.empty() && check_int > 0)
	{
		my_incflo.WriteCheckPointFile(check_file, nstep, dt, time);
		last_chk = nstep;
	}

	bool do_not_evolve =
		!steady_state && ((max_step == 0) || ((stop_time >= 0.) && (time > stop_time)) ||
						  ((stop_time <= 0.) && (max_step <= 0)));

	{ // Start profiling solve here

		BL_PROFILE("incflo_solve");
		BL_PROFILE_REGION("incflo_solve");

		if(!do_not_evolve)
		{
			while(finish == 0)
			{
				Real strt_step = ParallelDescriptor::second();

				if(!steady_state && regrid_int > -1 && nstep % regrid_int == 0)
                    my_incflo.Regrid(lev);

				my_incflo.Advance(lev, nstep, steady_state, dt, prev_dt, time, stop_time);

				Real end_step = ParallelDescriptor::second() - strt_step;
				ParallelDescriptor::ReduceRealMax(end_step,
												  ParallelDescriptor::IOProcessorNumber());
			    amrex::Print() << "Time per step " << end_step << std::endl;

				if(!steady_state)
				{
					time += prev_dt;
					nstep++;

					if((plot_int > 0) && (nstep % plot_int == 0))
					{
						my_incflo.incflo_compute_vort(lev);
						my_incflo.WritePlotFile(plot_file, nstep, dt, time);
						last_plt = nstep;
					}

					if((check_int > 0) && (nstep % check_int == 0))
					{
						my_incflo.WriteCheckPointFile(check_file, nstep, dt, time);
						last_chk = nstep;
					}
				}

				// Mechanism to terminate incflo normally.
				do_not_evolve =
					steady_state || (((stop_time >= 0.) && (time + 0.1 * dt >= stop_time)) ||
									 (max_step >= 0 && nstep >= max_step));
				if(do_not_evolve)
					finish = 1;
			}
		}
	}

	if(steady_state)
		nstep = 1;

	// Dump plotfile at the final time
	if(check_int > 0 && nstep != last_chk)
		my_incflo.WriteCheckPointFile(check_file, nstep, dt, time);
	if(plot_int > 0 && nstep != last_plt)
    {
        my_incflo.incflo_compute_vort(lev);
		my_incflo.WritePlotFile(plot_file, nstep, dt, time);
    }

	Real end_time = ParallelDescriptor::second() - strt_time;
	ParallelDescriptor::ReduceRealMax(end_time, ParallelDescriptor::IOProcessorNumber());

    amrex::Print() << "Time spent in main      " << end_time << std::endl;
    amrex::Print() << "Time spent in main-init " << end_time - end_init << std::endl;

	BL_PROFILE_REGION_STOP("incflo::main()");
	BL_PROFILE_VAR_STOP(pmain);

	amrex::Finalize();
	return 0;
}


void ReadParameters()
{
	{
		ParmParse pp("amr");

		pp.query("verbose", verbose);

		pp.query("stop_time", stop_time);
		pp.query("max_step", max_step);
		pp.query("steady_state", steady_state);

		pp.query("check_file", check_file);
		pp.query("check_int", check_int);
		pp.query("restart", restart_file);

		pp.query("plot_file", plot_file);
		pp.query("plot_int", plot_int);
		pp.query("write_eb_surface", write_eb_surface);

		pp.query("repl_x", repl_x);
		pp.query("repl_y", repl_y);
		pp.query("repl_z", repl_z);
		pp.query("regrid_int", regrid_int);
	}
}

