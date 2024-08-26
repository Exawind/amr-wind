#include <chrono>
#include <ctime>
#include "amr-wind/utilities/console_io.H"
#include "amr-wind/AMRWindVersion.H"
#include "AMReX.H"
#include "AMReX_OpenMP.H"
#include "amr-wind/CFDSim.H"

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
    if (irank != 0) {
        return;
    }
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
    if (irank != 0) {
        return;
    }
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
    if (irank != 0) {
        return;
    }
#else
    amrex::ignore_unused(comm);
#endif

    auto etime = std::chrono::system_clock::now();
    auto etimet = std::chrono::system_clock::to_time_t(etime);
    amrex::Array<char, 64> time_buf;
    ctime_r(&etimet, time_buf.begin());
    const std::string tstamp(time_buf.begin());

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
#elif defined(AMREX_USE_SYCL)
        << "(Backend: SYCL)"
#endif
        << std::endl
#else
        << "OFF" << std::endl
#endif
        << "  OpenMP           :: "
#ifdef AMREX_USE_OMP
        << "ON    (max threads = " << amrex::OpenMP::get_max_threads()
        << ")" << std::endl
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

void print_nonlinear_residual(
    CFDSim& sim, ScratchField& vel_diff, ScratchField& vel_star)
{
    // using namespace amrex;
    const int nlevels = sim.repo().num_active_levels();
    const auto& geom = sim.mesh().Geom();
    const auto& mesh = sim.mesh();

    const auto& velocity_new = sim.pde_manager().icns().fields().field;
    // const auto& vel_np1_old = sim.repo().get_field("vel_np1_old");
    // auto& vel_diff = sim.repo().get_field("vel_diff");

    for (int lev = 0; lev < nlevels; ++lev) {

        amrex::iMultiFab level_mask;
        if (lev < nlevels - 1) {
            level_mask = makeFineMask(
                mesh.boxArray(lev), mesh.DistributionMap(lev),
                mesh.boxArray(lev + 1), mesh.refRatio(lev), 1, 0);
        } else {
            level_mask.define(
                mesh.boxArray(lev), mesh.DistributionMap(lev), 1, 0,
                amrex::MFInfo());
            level_mask.setVal(1);
        }

        const auto& vnew = velocity_new(lev);
        // const auto& velstar = vel_np1_old(lev);
        // auto& veldiff = vel_diff(lev);
        const auto& velstar = vel_star(lev);
        auto& veldiff = vel_diff(lev);

        for (amrex::MFIter mfi(velocity_new(lev)); mfi.isValid(); ++mfi) {
            const auto bx = mfi.tilebox();
            const auto& velstar_arr = velstar.const_array(mfi);
            const auto& veldiff_arr = veldiff.array(mfi);
            const auto& velnew_arr = vnew.const_array(mfi);
            const auto& levelmask_arr = level_mask.const_array(mfi);

            amrex::ParallelFor(
                bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    veldiff_arr(i, j, k, 0) =
                        (velnew_arr(i, j, k, 0) - velstar_arr(i, j, k, 0)) *
                        levelmask_arr(i, j, k);
                    veldiff_arr(i, j, k, 1) =
                        (velnew_arr(i, j, k, 1) - velstar_arr(i, j, k, 1)) *
                        levelmask_arr(i, j, k);
                    veldiff_arr(i, j, k, 2) =
                        (velnew_arr(i, j, k, 2) - velstar_arr(i, j, k, 2)) *
                        levelmask_arr(i, j, k);
                });
        }
    }

    amrex::Real rms_ucell = 0.0;
    amrex::Real rms_vcell = 0.0;
    amrex::Real rms_wcell = 0.0;

    for (int lev = 0; lev < nlevels; ++lev) {
        rms_ucell += vel_diff(lev).norm2(0);
        rms_vcell += vel_diff(lev).norm2(1);
        rms_wcell += vel_diff(lev).norm2(2);
    }
    amrex::Print() << "Non-linear residual for u velocity " << rms_ucell
                   << std::endl;
    amrex::Print() << "Non-linear residual for v velocity " << rms_vcell
                   << std::endl;
    amrex::Print() << "Non-linear residual for w velocity " << rms_wcell
                   << std::endl;
}

} // namespace amr_wind::io
