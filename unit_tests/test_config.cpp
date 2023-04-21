/** \file test_config.cpp
 *
 *  Tests various configurations for GPU builds
 */

#include "gtest/gtest.h"
#include "amr-wind/AMRWindVersion.H"
#include "amr-wind/utilities/console_io.H"
#include "AMReX_ccse-mpi.H"
#include "AMReX_ParallelDescriptor.H"
#include "AMReX_Print.H"
#include "AMReX_Gpu.H"

namespace amrex {
const char* buildInfoGetGitHash(int i);
}

namespace amr_wind_tests {

TEST(Configuration, Build)
{
    const std::string dirty_tag =
        (amr_wind::version::amr_wind_dirty_repo == "DIRTY")
            ? ("-" + amr_wind::version::amr_wind_dirty_repo)
            : "";
    const std::string awind_git_sha =
        amr_wind::version::amr_wind_git_sha + dirty_tag;
    const char* amrex_git = amrex::buildInfoGetGitHash(2);
    amrex::Print() << "AMR-Wind SHA = " << awind_git_sha
                   << "\nAMReX    SHA = " << amrex_git << std::endl;
}

TEST(Configuration, MPI)
{
#ifdef AMREX_USE_MPI
    int nprocs = amrex::ParallelDescriptor::NProcs();

    amrex::Print() << "MPI configuration: " << nprocs << " ranks" << std::endl;
    char mpi_lib_ver[MPI_MAX_LIBRARY_VERSION_STRING];
    int len;
    MPI_Get_library_version(mpi_lib_ver, &len);
    amrex::Print() << mpi_lib_ver << std::endl;
#else
    amrex::Print() << "AMR-Wind not built with MPI support." << std::endl;
    GTEST_SKIP();
#endif
}

TEST(Configuration, GPU)
{
#ifdef AMREX_USE_GPU
#ifdef AMREX_USE_CUDA
    amrex::Print() << "GPU backend: CUDA" << std::endl;
#if defined(CUDA_VERSION)
    amrex::Print() << "CUDA_VERSION: " << CUDA_VERSION << " "
                   << CUDA_VERSION / 1000 << "." << (CUDA_VERSION % 1000) / 10
                   << std::endl;
#endif
#elif defined(AMREX_USE_HIP)
    amrex::Print() << "GPU backend: HIP" << std::endl;
#elif defined(AMREX_USE_SYCL)
    amrex::Print() << "GPU backend: SYCL" << std::endl;
#endif

    using Dev = amrex::Gpu::Device;
    const int myrank = amrex::ParallelDescriptor::MyProc();
    std::stringstream ss;
    // clang-format off
    ss << "[" << myrank << "] " << Dev::deviceId()
       << ": " << Dev::deviceName() << "\n"
       << "    Warp size          : " << Dev::warp_size << "\n"
       << "    Global memory      : "
       << (static_cast<double>(Dev::totalGlobalMem()) / (1 << 30)) << "GB\n"
       << "    Shared mem/ block  : "
       << (Dev::sharedMemPerBlock() / (1 << 10)) << "KB\n"
       << "    Max. threads/block : " << Dev::maxThreadsPerBlock()
       << " (" << Dev::maxThreadsPerBlock(0) << ", "
       << Dev::maxThreadsPerBlock(1) << ", " << Dev::maxThreadsPerBlock(2) << ")\n"
       << "    Max. blocks/grid   : (" << Dev::maxBlocksPerGrid(0)
       << ", " << Dev::maxBlocksPerGrid(1) << ", " << Dev::maxBlocksPerGrid(2) << ")\n"
       << std::endl;
    // clang-format on
    amrex::OutStream() << ss.str();
#else
    amrex::Print() << "AMR-Wind not built with GPU support" << std::endl;
    GTEST_SKIP();
#endif
}

TEST(Configuration, TPLs)
{
    if (amrex::ParallelDescriptor::IOProcessor()) {
        amr_wind::io::print_tpls(std::cout);
    }
}

} // namespace amr_wind_tests
