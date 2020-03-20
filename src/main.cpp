#include <incflo.H>
#include <AMReX_buildInfo.H>

void writeBuildInfo();

using namespace amrex;

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

        BL_PROFILE("main()");

        // Issue an error if input file is not given 
        if(argc < 2) amrex::Abort("Input file must be given as command-line argument.");

        // Write out the incflo git hash (the AMReX git hash is already written)
        const char* githash_incflo = buildInfoGetGitHash(1);
        amrex::Print() << "incflo git hash: " << githash_incflo << "\n";

        // Start timing the program
        Real start_time = ParallelDescriptor::second();

        // Default constructor. Note inheritance: incflo : AmrCore : AmrMesh. 
        incflo my_incflo;

        // Initialize data, parameters, arrays and derived internals
        my_incflo.InitData();

        // Time spent on initialization
        Real init_time = ParallelDescriptor::second() - start_time;

        // Evolve system to final time
        my_incflo.Evolve();

        // Time spent in total
        Real end_time = ParallelDescriptor::second() - start_time;

        ParallelDescriptor::ReduceRealMax(init_time, ParallelDescriptor::IOProcessorNumber());
        ParallelDescriptor::ReduceRealMax(end_time, ParallelDescriptor::IOProcessorNumber());

        // Print timing results
        amrex::Print() << "Time spent in InitData():    " << init_time << std::endl;
        amrex::Print() << "Time spent in Evolve():      " << end_time - init_time << std::endl;
    }
    amrex::Finalize();
}
