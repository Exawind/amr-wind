#include <chrono>
#include <ctime>
#include "amr-wind/utilities/console_io.H"
#include "amr-wind/AMRWindVersion.H"
#include "AMReX.H"

#ifdef AMR_WIND_USE_NETCDF
#include "netcdf.h"
#ifdef NC_HAVE_META_H
#include "netcdf_meta.h"
#endif
#endif
#ifdef AMR_WIND_USE_MASA
#include "masa.h"
#endif
#ifdef AMREX_USE_HYPRE
#include "HYPRE_config.h"
#endif
#ifdef AMR_WIND_USE_ASCENT
#include "ascent_config.h"
#endif

namespace amrex {
const char* buildInfoGetBuildDate();
const char* buildInfoGetComp();
const char* buildInfoGetGitHash(int i);
const char* buildInfoGetCompVersion();
} // namespace amrex

namespace amr_wind::io {

namespace {
const std::string dbl_line = std::string(78, '=') + "\n";
const std::string dash_line = "\n" + std::string(78, '-') + "\n";
} // namespace

void print_usage(MPI_Comm comm, std::ostream& out)
{
#ifdef AMREX_USE_MPI
    int irank = 0;
    int num_ranks = 1;
    MPI_Comm_size(comm, &num_ranks);
    MPI_Comm_rank(comm, &irank);

    // Only root process does the printing
    if (irank != 0) return;
#else
    amrex::ignore_unused(comm);
#endif

    out << R"doc(Usage:
    amr_wind <input_file> [param=value] [param=value] ...

Required:
    input_file   : Input file with simulation settings

Optional:
    param=value  : Overrides for parameters during runtime
)doc" << std::endl;
}

void print_error(MPI_Comm comm, const std::string& msg)
{
#ifdef AMREX_USE_MPI
    int irank = 0;
    int num_ranks = 1;
    MPI_Comm_size(comm, &num_ranks);
    MPI_Comm_rank(comm, &irank);

    // Only root process does the printing
    if (irank != 0) return;
#else
    amrex::ignore_unused(comm);
#endif

    std::cout << "ERROR: " << msg << std::endl;
}

void print_banner(MPI_Comm comm, std::ostream& out)
{
#ifdef AMREX_USE_MPI
    int irank = 0;
    int num_ranks = 1;
    MPI_Comm_size(comm, &num_ranks);
    MPI_Comm_rank(comm, &irank);

    // Only root process does the printing
    if (irank != 0) return;
#else
    amrex::ignore_unused(comm);
#endif

    auto etime = std::chrono::system_clock::now();
    auto etimet = std::chrono::system_clock::to_time_t(etime);
    char time_buf[64];
    ctime_r(&etimet, time_buf);
    const std::string tstamp(time_buf);

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
        << "  Exec. time       :: " << tstamp
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
        << "ON    "
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
#ifdef AMREX_USE_OMP
        << "ON    (Num. threads = " << omp_get_max_threads() << ")" << std::endl
#else
        << "OFF" << std::endl
#endif
        << std::endl;

    print_tpls(out);

    out << "           This software is released under the BSD 3-clause license.           "
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

void print_tpls(std::ostream& out)
{
    amrex::Vector<std::string> tpls;

#ifdef AMR_WIND_USE_NETCDF
    tpls.push_back(std::string("NetCDF    ") + NC_VERSION);
#endif
#ifdef AMREX_USE_HYPRE
    tpls.push_back(std::string("HYPRE     ") + HYPRE_RELEASE_VERSION);
#endif
#ifdef AMR_WIND_USE_OPENFAST
    tpls.push_back(std::string("OpenFAST  "));
#endif
#ifdef AMR_WIND_USE_MASA
    tpls.push_back(std::string("MASA      ") + MASA_LIB_VERSION);
#endif
#ifdef AMR_WIND_USE_ASCENT
    tpls.push_back(std::string("ASCENT    ") + ASCENT_VERSION);
#endif

    if (!tpls.empty()) {
        out << "  Enabled third-party libraries: ";
        for (const auto& val : tpls) {
            out << "\n    " << val;
        }
        out << std::endl << std::endl;
    } else {
        out << "  No additional third-party libraries enabled" << std::endl
            << std::endl;
    }
}

} // namespace amr_wind::io
