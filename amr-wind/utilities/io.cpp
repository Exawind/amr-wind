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

    amrex::ParmParse pp("amr");
    amrex::Vector<int> rep{{1, 1, 1}};
    pp.queryarr("replicate", rep);
    IntVect Nrep(rep[0], rep[1], rep[2]);

    bool replicate = (Nrep == IntVect::TheUnitVector()) ? false : true;

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

    if (replicate) {
        for (int d = 0; d < BL_SPACEDIM; d++) {
            prob_lo[d] = Nrep[d] * prob_lo[d];
            prob_hi[d] = Nrep[d] * prob_hi[d];
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

        Box orig_domain(ba_inp[lev].minimalBox());

        if (replicate) {
            amrex::Print() << " OLD BA had " << ba_inp[lev].size() << " GRIDS "
                           << std::endl;
            amrex::Print() << " OLD Domain" << orig_domain << std::endl;
        }

        BoxList bl;
        for (int nb = 0; nb < ba_inp[lev].size(); nb++) {
            for (int k = 0; k < Nrep[2]; k++) {
                for (int j = 0; j < Nrep[1]; j++) {
                    for (int i = 0; i < Nrep[0]; i++) {
                        Box b(ba_inp[lev][nb]);
                        IntVect shift_vec(
                            i * orig_domain.length(0),
                            j * orig_domain.length(1),
                            k * orig_domain.length(2));
                        b.shift(shift_vec);
                        bl.push_back(b);
                    }
                }
            }
        }
        BoxArray ba_rep;
        ba_rep.define(bl);

        if (replicate) {
            amrex::Print() << " NEW BA had " << ba_rep.size() << " GRIDS "
                           << std::endl;
            amrex::Print() << " NEW Domain" << ba_rep.minimalBox() << std::endl;
        }

        // Create distribution mapping
        dm_inp[lev].define(ba_inp[lev], ParallelDescriptor::NProcs());
        DistributionMapping dm = dm_inp[lev];

        BoxArray ba(ba_rep.simplified());
        ba.maxSize(maxGridSize(lev));
        if (refine_grid_layout)
            ChopGrids(lev, ba, ParallelDescriptor::NProcs());
        dm = DistributionMapping{ba, ParallelDescriptor::NProcs()};

        MakeNewLevelFromScratch(lev, m_time.current_time(), ba, dm);
    }

    m_sim.io_manager().read_checkpoint_fields(
        restart_file, ba_inp, dm_inp, Nrep);
}
