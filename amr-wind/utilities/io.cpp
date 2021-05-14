#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>
#include "amr-wind/incflo.H"
#include "amr-wind/core/Physics.H"
#include "amr-wind/utilities/console_io.H"
#include "amr-wind/utilities/IOManager.H"

using namespace amrex;

namespace {
const std::string level_prefix{"Level_"};
}

void GotoNextLine(std::istream& is)
{
    constexpr std::streamsize bl_ignore_max{100000};
    is.ignore(bl_ignore_max, '\n');
}

void incflo::ReadCheckpointFile()
{
    BL_PROFILE("amr-wind::incflo::ReadCheckpointFile()");

    const std::string& restart_file = m_sim.io_manager().restart_file();
    amrex::Print() << "Restarting from checkpoint " << restart_file
                   << std::endl;

    Real prob_lo[BL_SPACEDIM];
    Real prob_hi[BL_SPACEDIM];

    /***************************************************************************
     * Load header: set up problem domain (including BoxArray)                 *
     *              allocate incflo memory (incflo::AllocateArrays)            *
     *              (by calling MakeNewLevelFromScratch)
     ***************************************************************************/

    std::string File(restart_file + "/Header");

    VisMF::IO_Buffer io_buffer(VisMF::GetIOBufferSize());

    Vector<char> fileCharPtr;
    ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
    std::string fileCharPtrString(fileCharPtr.dataPtr());
    std::istringstream is(fileCharPtrString, std::istringstream::in);

    std::string line, word;

    // Start reading from checkpoint file

    // Title line
    std::getline(is, line);

    // Finest level
    is >> finest_level;
    GotoNextLine(is);

    int nstep;
    amrex::Real cur_time, dt_restart;
    // Step count
    is >> nstep;
    GotoNextLine(is);

    // Current time
    is >> cur_time;
    GotoNextLine(is);

    m_time.set_restart_time(nstep, cur_time);

    // Time step size
    is >> dt_restart;
    GotoNextLine(is);

    is >> m_time.deltaTNm1();
    GotoNextLine(is);

    is >> m_time.deltaTNm2();
    GotoNextLine(is);

    // Low coordinates of domain bounding box
    std::getline(is, line);
    {
        std::istringstream lis(line);
        int i = 0;
        while (lis >> word) {
            prob_lo[i++] = std::stod(word);
        }
    }

    // High coordinates of domain bounding box
    std::getline(is, line);
    {
        std::istringstream lis(line);
        int i = 0;
        while (lis >> word) {
            prob_hi[i++] = std::stod(word);
        }
    }

    // Set up problem domain
    RealBox rb(prob_lo, prob_hi);
    Geometry::ResetDefaultProbDomain(rb);
    for (int lev = 0; lev <= max_level; ++lev) {
        SetGeometry(
            lev, Geometry(
                     Geom(lev).Domain(), rb, Geom(lev).CoordInt(),
                     Geom(lev).isPeriodic()));
    }

    amrex::Vector<amrex::BoxArray> ba_inp(finest_level + 1);
    amrex::Vector<amrex::DistributionMapping> dm_inp(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        // read in level 'lev' BoxArray from Header
        ba_inp[lev].readFrom(is);
        GotoNextLine(is);

        // Create distribution mapping
        dm_inp[lev].define(ba_inp[lev], ParallelDescriptor::NProcs());
        DistributionMapping dm = dm_inp[lev];

        BoxArray ba(ba_inp[lev].simplified());
        ba.maxSize(maxGridSize(lev));
        if (ba != ba_inp[lev]) {
            if ((lev == 0) || refine_grid_layout)
                ChopGrids(lev, ba, ParallelDescriptor::NProcs());
            dm = DistributionMapping{ba, ParallelDescriptor::NProcs()};
        }

        MakeNewLevelFromScratch(lev, m_time.current_time(), ba, dm);
    }

    m_sim.io_manager().read_checkpoint_fields(restart_file, ba_inp, dm_inp);
}
