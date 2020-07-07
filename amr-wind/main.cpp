#include "amr-wind/incflo.H"
#include "amr-wind/utilities/console_io.H"

int main(int argc, char* argv[])
{
#ifdef AMREX_USE_MPI
    MPI_Init(&argc, &argv);
#endif

    // check to see if it contains --describe
    if (argc >= 2) {
        for (auto i = 1; i < argc; i++) {
            if (std::string(argv[i]) == "--describe") {
                amr_wind::io::print_banner(std::cout);
                return 0;
            }
        }
    }
    amr_wind::io::print_banner(std::cout);

    amrex::Initialize(argc, argv, true, MPI_COMM_WORLD);
    { /* These braces are necessary to ensure amrex::Finalize() can be called without explicitly
        deleting all the incflo member MultiFabs */

        BL_PROFILE("amr-wind::main()");

        // Issue an error if input file is not given 
        if(argc < 2) amrex::Abort("Input file must be given as command-line argument.");

        // Start timing the program
        amrex::Real start_time = amrex::ParallelDescriptor::second();
        amrex::Print() << "Initializing AMR-Wind ..." << std::endl;

        // Default constructor. Note inheritance: incflo : AmrCore : AmrMesh. 
        incflo my_incflo;

        // Initialize data, parameters, arrays and derived internals
        my_incflo.InitData();

        // Time spent on initialization
        amrex::Real init_time = amrex::ParallelDescriptor::second() - start_time;
        amrex::Print() << "Initialization successful. Time elapsed = " << init_time << std::endl;

        // Evolve system to final time
        my_incflo.Evolve();

        // Time spent in total
        amrex::Real end_time = amrex::ParallelDescriptor::second() - start_time;

        amrex::ParallelDescriptor::ReduceRealMax(
            init_time, amrex::ParallelDescriptor::IOProcessorNumber());
        amrex::ParallelDescriptor::ReduceRealMax(
            end_time, amrex::ParallelDescriptor::IOProcessorNumber());

        // Print timing results
        amrex::Print() << "Time spent in InitData():    " << init_time << std::endl;
        amrex::Print() << "Time spent in Evolve():      " << end_time - init_time << std::endl;
    }
    amrex::Finalize();
#ifdef AMREX_USE_MPI
    MPI_Finalize();
#endif
}
