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

TEST(Configuration, CUDA)
{
#ifdef AMREX_USE_CUDA
#if defined(CUDA_VERSION)
    amrex::Print() << "CUDA configuration: "
                   << "CUDA_VERSION: " << CUDA_VERSION << " "
                   << CUDA_VERSION / 1000 << "." << (CUDA_VERSION % 1000) / 10
                   << std::endl;
#endif
    cudaError_t error;
    int ndevices;
    cudaDeviceProp dev;

    error = cudaGetDeviceCount(&ndevices);
    if (error != cudaSuccess) {
        ADD_FAILURE() << cudaGetErrorString(error);
        return;
    } else {
        std::cout << "Num. devices = " << ndevices << std::endl;
    }

    const int myrank = amrex::ParallelDescriptor::MyProc();
    const int rankDevice = amrex::Gpu::Device::deviceId();
    error = cudaGetDeviceProperties(&dev, rankDevice);
    if (error != cudaSuccess) {
        ADD_FAILURE() << cudaGetErrorString(error);
        return;
    }
    char busid[512];
    cudaDeviceGetPCIBusId(busid, 512, rankDevice);
    std::cout << "[" << myrank << "] " << rankDevice << ": " << dev.name
              << " CC: " << dev.major << "." << dev.minor << " ID: " << busid
              << " GM: "
              << (static_cast<double>(dev.totalGlobalMem) / (1 << 30)) << "GB"
              << " ShMem/Blk: " << (dev.sharedMemPerBlock / (1 << 10)) << "KB"
              << std::endl;
#else
    amrex::Print() << "AMR-Wind not build with CUDA support" << std::endl;
    GTEST_SKIP();
#endif
}

TEST(Configuration, TPLs)
{
    if (amrex::ParallelDescriptor::IOProcessor())
        amr_wind::io::print_tpls(std::cout);
}

} // namespace amr_wind_tests
