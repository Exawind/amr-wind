#include <chrono>
#include <ctime>
#include "console_io.H"
#include "AMReX.H"
#include "AMReX_buildInfo.H"

namespace amr_wind {
namespace io {

namespace {
const std::string dbl_line = std::string(78, '=') + "\n";
const std::string dash_line = "\n" + std::string(78, '-') + "\n";
} // namespace

void print_banner(std::ostream& out)
{
    int irank = 0;
#ifdef AMREX_USE_MPI
    int num_ranks = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);
    MPI_Comm_rank(MPI_COMM_WORLD, &irank);
#endif
    if (irank != 0) return;

    auto exec_time = std::chrono::system_clock::now();
    auto exect = std::chrono::system_clock::to_time_t(exec_time);

    const char* git_hash = amrex::buildInfoGetGitHash(1);
    const char* amrex_hash = amrex::buildInfoGetGitHash(2);
    const std::string git_str =
        (strlen(git_hash) > 0) ? std::string(git_hash) : std::string("Unknown");
    const std::string amrex_git_str =
        (strlen(amrex_hash) > 0) ? std::string(amrex_hash) : std::string("Unknown");

    // clang-format off
    out << dbl_line
        << "                AMR-Wind (https://github.com/exawind/amr-wind)"
        << std::endl << std::endl
        << "  AMR-Wind Git SHA :: " << git_hash << std::endl
        << "  AMReX version    :: " << amrex::Version() << " ( "  << amrex_git_str << " )" << std::endl << std::endl
        << "  Exec. date       :: " << std::ctime(&exect)
        << "  Build date       :: " << amrex::buildInfoGetBuildDate() << std::endl
        << "  C++ compiler     :: " << amrex::buildInfoGetComp() << " " << amrex::buildInfoGetCompVersion() << std::endl << std::endl
        << "  MPI              :: "
#ifdef AMREX_USE_MPI
        << "ON    (Num. ranks = " << num_ranks << ")" << std::endl
#else
        << "OFF " << std::endl
#endif
        << "  GPU              :: "
#ifdef AMREX_USE_GPU
        << "ON" << std::endl
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
                   << std::setw(6) << std::right << "Iters"
                   << std::setw(22) << std::right << "Initial residual"
                   << std::setw(22) << std::right << "Final residual" << std::endl
                   << "  ----------------------------------------------------------------------------"
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
