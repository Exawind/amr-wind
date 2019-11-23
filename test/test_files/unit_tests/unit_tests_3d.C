#include <gtest/gtest.h>
//#include <incflo.H>
//#include <mpi.h>

int gl_argc = 0;
char** gl_argv = 0;

int main(int argc, char **argv)
{
    //MPI_Init(&argc, &argv);

    //amrex::Initialize(argc, argv);
    //{
    int returnVal = 0;

    {
      testing::InitGoogleTest(&argc, argv);

      gl_argc = argc;
      gl_argv = argv;

      returnVal = RUN_ALL_TESTS();
    }

    return returnVal;
    //}
    //amrex::Finalize();
}

