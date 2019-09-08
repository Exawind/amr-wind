#include <incflo.H>
#include <AMReX_buildInfo.H>

// TODO: find better way of making the fillpatch stuff work
void set_ptr_to_incflo(incflo& my_incflo);
 
void writeBuildInfo();

int main(int argc, char* argv[])
{

    // check to see if it contains --describe
    if (argc >= 2) {
        for (auto i = 1; i < argc; i++) {
            if (std::string(argv[i]) == "--describe") {
                writeBuildInfo();
                return 0;
            }
        }
    }

    amrex::Initialize(argc, argv);
    { /* These braces are necessary to ensure amrex::Finalize() can be called without explicitly
        deleting all the incflo member MultiFabs */

    BL_PROFILE_VAR("main()", pmain)

    // Issue an error if input file is not given 
    if(argc < 2) amrex::Abort("Input file must be given as command-line argument.");

    // Write out the incflo git hash (the AMReX git hash is already written)
    const char* githash_incflo = buildInfoGetGitHash(1);
    amrex::Print() << "incflo git hash: " << githash_incflo << "\n";

    // Start timing the program
    Real start_time = ParallelDescriptor::second();

    // Default constructor. Note inheritance: incflo : AmrCore : AmrMesh. 
    incflo my_incflo;

    // Get boundary conditions from inputs file
    my_incflo.GetInputBCs();

    // Set global static pointer to incflo object. Used by fillpatch utility
    set_ptr_to_incflo(my_incflo);

    // Initialize data, parameters, arrays and derived internals
    my_incflo.InitData();

    // Time spent on initialization
    Real init_time = ParallelDescriptor::second() - start_time;
    ParallelDescriptor::ReduceRealMax(init_time, ParallelDescriptor::IOProcessorNumber());

    // Evolve system to final time 
    my_incflo.Evolve();

    // Time spent in total 
    Real end_time = ParallelDescriptor::second() - start_time;
    ParallelDescriptor::ReduceRealMax(end_time, ParallelDescriptor::IOProcessorNumber());

    // Print timing results
    amrex::Print() << "Time spent in InitData():    " << init_time << std::endl;
    amrex::Print() << "Time spent in Evolve():      " << end_time - init_time << std::endl;

    BL_PROFILE_VAR_STOP(pmain);
    }
	amrex::Finalize();
}
