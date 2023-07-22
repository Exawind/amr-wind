#include "CoarsenCheckpt.H"
#include "amr-wind/utilities/console_io.H"

int main(int argc, char* argv[])
{
#ifdef AMREX_USE_MPI
    MPI_Init(&argc, &argv);
#endif

    amr_wind::io::print_banner(MPI_COMM_WORLD, std::cout);

    amrex::Initialize(argc, argv, true, MPI_COMM_WORLD, []() {
        {
            amrex::ParmParse pp("amrex");
            // Set the defaults so that we throw an exception instead of
            // attempting to generate backtrace files. However, if the user has
            // explicitly set these options in their input files respect those
            // settings.
            if (!pp.contains("throw_exception")) pp.add("throw_exception", 1);
            if (!pp.contains("signal_handling")) pp.add("signal_handling", 0);
        }
        {
            // Preprocess input file
            amrex::ParmParse pp("amr");
            int max_lev_input = 0;
            pp.get("max_level", max_lev_input);
            pp.add("max_level", max_lev_input + 1);
            // Getting the right number of levels is important for geometry init
        }
    });

    {
        BL_PROFILE("coarsen-chkpt::main");
        amr_wind::tools::CoarsenCheckpt obj;
        obj.run_utility();
    }

    amrex::Finalize();

#ifdef AMREX_USE_MPI
    MPI_Finalize();
#endif

    return 0;
}
