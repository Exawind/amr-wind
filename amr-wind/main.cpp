#include "amr-wind/incflo.H"
#include "amr-wind/utilities/console_io.H"

#include "AMReX_FileSystem.H"

int main(int argc, char* argv[])
{
#ifdef AMREX_USE_MPI
    MPI_Init(&argc, &argv);
#endif

    if (argc >= 2) {
        // Look for "-h" or "--help" flag and print usage
        for (auto i = 1; i < argc; i++) {
            const std::string param(argv[i]);
            if ((param == "--help") || (param == "-h")) {
                amr_wind::io::print_banner(MPI_COMM_WORLD, std::cout);
                amr_wind::io::print_usage(MPI_COMM_WORLD, std::cout);
                return 0;
            }
        }
    } else if (argc < 2) {
        // Print usage and exit with error code if no input file was provided.
        amr_wind::io::print_usage(MPI_COMM_WORLD, std::cout);
        amr_wind::io::print_error(
            MPI_COMM_WORLD, "No input file provided. Exiting!!");
        return 1;
    }
    if (!amrex::FileSystem::Exists(std::string(argv[1]))) {
        // Print usage and exit with error code if we cannot find the input file
        amr_wind::io::print_usage(MPI_COMM_WORLD, std::cout);
        amr_wind::io::print_error(
            MPI_COMM_WORLD, "Input file does not exist = " +
                                std::string(argv[1]) + ". Exiting!!");
        return 1;
    }

    amr_wind::io::print_banner(MPI_COMM_WORLD, std::cout);
    amrex::Initialize(argc, argv, true, MPI_COMM_WORLD, []() {
        amrex::ParmParse pp("amrex");
        // Set the defaults so that we throw an exception instead of attempting
        // to generate backtrace files. However, if the user has explicitly set
        // these options in their input files respect those settings.
        if (!pp.contains("throw_exception")) pp.add("throw_exception", 1);
        if (!pp.contains("signal_handling")) pp.add("signal_handling", 0);
    });

    { /* These braces are necessary to ensure amrex::Finalize() can be called
        without explicitly deleting all the incflo member MultiFabs */

        BL_PROFILE("amr-wind::main()");

        // Start timing the program
        amrex::Real start_time = amrex::ParallelDescriptor::second();
        amrex::Print() << "Initializing AMR-Wind ..." << std::endl;

        // Default constructor. Note inheritance: incflo : AmrCore : AmrMesh.
        incflo my_incflo;

        // Initialize data, parameters, arrays and derived internals
        my_incflo.InitData();

        // Time spent on initialization
        amrex::Real init_time =
            amrex::ParallelDescriptor::second() - start_time;
        amrex::Print() << "Initialization successful. Time elapsed = "
                       << init_time << std::endl;

        // Evolve system to final time
        my_incflo.Evolve();

        // Time spent in total
        amrex::Real end_time = amrex::ParallelDescriptor::second() - start_time;

        amrex::ParallelDescriptor::ReduceRealMax(
            init_time, amrex::ParallelDescriptor::IOProcessorNumber());
        amrex::ParallelDescriptor::ReduceRealMax(
            end_time, amrex::ParallelDescriptor::IOProcessorNumber());

        // Print timing results
        amrex::Print() << "Time spent in InitData():    " << init_time
                       << std::endl;
        amrex::Print() << "Time spent in Evolve():      "
                       << end_time - init_time << std::endl;
    }
    amrex::Finalize();
#ifdef AMREX_USE_MPI
    MPI_Finalize();
#endif

    return 0;
}
