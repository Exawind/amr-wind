#include <chrono>
#include <ctime>
#include "amr-wind/utilities/console_io.H"
#include "amr-wind/AMRWindVersion.H"
#include "AMReX.H"

namespace amrex {
const char* buildInfoGetBuildDate();
const char* buildInfoGetComp();
const char* buildInfoGetGitHash(int i);
const char* buildInfoGetCompVersion();
} // namespace amrex

namespace amr_wind {
namespace io {

namespace {
const std::string dbl_line = std::string(78, '=') + "\n";
const std::string dash_line = "\n" + std::string(78, '-') + "\n";
} // namespace

void print_banner(MPI_Comm comm, std::ostream& out)
{
#ifdef AMREX_USE_MPI
    int irank = 0;
    int num_ranks = 1;
    MPI_Comm_size(comm, &num_ranks);
    MPI_Comm_rank(comm, &irank);

    // Only root process does the printing
    if (irank != 0) return;
#endif

    auto exec_time = std::chrono::system_clock::now();
    auto exect = std::chrono::system_clock::to_time_t(exec_time);
    const std::string dirty_tag = (version::amr_wind_dirty_repo == "DIRTY")
                                      ? ("-" + version::amr_wind_dirty_repo)
                                      : "";
    const std::string awind_version = version::amr_wind_version + dirty_tag;
    const std::string awind_git_sha = version::amr_wind_git_sha + dirty_tag;

    // clang-format off
    out << dbl_line
        << "                AMR-Wind (https://github.com/exawind/amr-wind)"
        << std::endl << std::endl
        << "  AMR-Wind version :: " << awind_version << std::endl
        << "  AMR-Wind Git SHA :: " << awind_git_sha << std::endl
        << "  AMReX version    :: " << amrex::Version() << std::endl << std::endl
        << "  Exec. time       :: " << std::ctime(&exect)
        << "  Build time       :: " << amrex::buildInfoGetBuildDate() << std::endl
        << "  C++ compiler     :: " << amrex::buildInfoGetComp()
        << " " << amrex::buildInfoGetCompVersion() << std::endl << std::endl
        << "  MPI              :: "
#ifdef AMREX_USE_MPI
        << "ON    (Num. ranks = " << num_ranks << ")" << std::endl
#else
        << "OFF " << std::endl
#endif
        << "  GPU              :: "
#ifdef AMREX_USE_GPU
        << "ON   "
#if defined(AMREX_USE_CUDA)
        << "(Backend: CUDA)"
#elif defined(AMREX_USE_HIP)
        << "(Backend: HIP)"
#elif defined(AMREX_USE_DPCPP)
        << "(Backend: SYCL)"
#endif
        << std::endl
#else
        << "OFF" << std::endl
#endif
        << "  OpenMP           :: "
#ifdef _OPENMP
        << "ON    (Num. threads = " << omp_get_max_threads() << ")" << std::endl
#else
        << "OFF" << std::endl << std::endl
#endif
        << "           This software is released under the BSD 3-clause license.           "
        << std::endl
        << " See https://github.com/Exawind/amr-wind/blob/development/LICENSE for details. "
        << dash_line << std::endl;
    // clang-format on
}

void print_mlmg_header(const std::string& key)
{
    const int name_width = 26;
    amrex::Print() << "\n" << key << std::endl;
    amrex::Print() << "  " << std::setw(name_width) << std::left << "System"
                   << std::setw(6) << std::right << "Iters" << std::setw(22)
                   << std::right << "Initial residual" << std::setw(22)
                   << std::right << "Final residual" << std::endl
                   << "  "
                      "--------------------------------------------------------"
                      "--------------------"
                   << std::endl;
}

void print_mlmg_info(const std::string& solve_name, const amrex::MLMG& mlmg)
{
    const int name_width = 26;
    amrex::Print() << "  " << std::setw(name_width) << std::left << solve_name
                   << std::setw(6) << std::right << mlmg.getNumIters()
                   << std::setw(22) << std::right << mlmg.getInitResidual()
                   << std::setw(22) << std::right << mlmg.getFinalResidual()
                   << std::endl;
}

} // namespace io
} // namespace amr_wind
