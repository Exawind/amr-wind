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

    Real prob_lo[AMREX_SPACEDIM] = {0.0};
    Real prob_hi[AMREX_SPACEDIM] = {0.0};

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

    amrex::Vector<amrex::Real> prob_lo_input(AMREX_SPACEDIM);
    amrex::Vector<amrex::Real> prob_hi_input(AMREX_SPACEDIM);

    {
        amrex::ParmParse pp("geometry");
        pp.getarr("prob_lo", prob_lo_input);
        pp.getarr("prob_hi", prob_hi_input);
    }

    amrex::Vector<int> n_cell_input(AMREX_SPACEDIM);

    {
        amrex::ParmParse pp("amr");
        pp.getarr("n_cell", n_cell_input);
    }

    IntVect rep(1, 1, 1);
    for (int d = 0; d < AMREX_SPACEDIM; d++) {
        AMREX_ALWAYS_ASSERT(prob_hi[d] > prob_lo[d]); // NOLINT

        const amrex::Real domain_ratio =
            (prob_hi_input[d] - prob_lo_input[d]) / (prob_hi[d] - prob_lo[d]);

        rep[d] = static_cast<int>(domain_ratio);

        constexpr amrex::Real domain_eps = 1.0e-6;
        if (std::abs(static_cast<amrex::Real>(rep[d]) - domain_ratio) >
            domain_eps) {
            amrex::Abort(
                "Domain size changed which indicates replication but there is "
                "either some precision issues or a non integer domain size "
                "change was input");
        }
    }

    bool replicate = (rep != IntVect::TheUnitVector());

    if (replicate) {
        amrex::Print() << "replicating restart file: " << rep << std::endl;
    }

    // Set up problem domain
    RealBox rb(prob_lo_input.data(), prob_hi_input.data());
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
    }

    // always use level 0 to check domain size
    constexpr int lev0{0};
    Box orig_domain(ba_inp[lev0].minimalBox());

    if (replicate) {
        amrex::Print() << " OLD BA had " << ba_inp[lev0].size() << " GRIDS "
                       << std::endl;
        amrex::Print() << " OLD Domain" << orig_domain << std::endl;
    }

    for (int lev = 0; lev <= finest_level; ++lev) {
        BoxList bl;
        for (int k = 0; k < rep[2]; k++) {
            for (int j = 0; j < rep[1]; j++) {
                for (int i = 0; i < rep[0]; i++) {
                    for (int nb = 0; nb < ba_inp[lev].size(); nb++) {
                        Box b(ba_inp[lev][nb]);
                        IntVect shift_vec(
                            i * orig_domain.length(0),
                            j * orig_domain.length(1),
                            k * orig_domain.length(2));

                        // equivalent to 2^lev
                        shift_vec *= (1 << lev);

                        b.shift(shift_vec);
                        bl.push_back(b);
                    }
                }
            }
        }
        BoxArray ba_rep;
        ba_rep.define(bl);

        if (replicate && lev == 0) {

            for (int d = 0; d < AMREX_SPACEDIM; d++) {
                auto new_domain = ba_rep.minimalBox();
                const auto* hi_vect = new_domain.hiVect();

                if (hi_vect[d] + 1 != n_cell_input[d]) {
                    amrex::Abort(
                        "input file error, domain size changed which indicated "
                        "replication, but the amr.n_cell is inconsistent with "
                        "that change in domain size, please adjust amr.n_cell, "
                        "or geometry.prob_lo, or geometry.prob_hi");
                }
            }

            amrex::Print() << " NEW BA had " << ba_rep.size() << " GRIDS "
                           << std::endl;
            amrex::Print() << " NEW Domain" << ba_rep.minimalBox() << std::endl;
        }

        // Create distribution mapping
        dm_inp[lev].define(ba_inp[lev], ParallelDescriptor::NProcs());

        BoxArray ba(ba_rep.simplified());
        ba.maxSize(maxGridSize(lev));
        if (refine_grid_layout) {
            ChopGrids(lev, ba, ParallelDescriptor::NProcs());
        }
        DistributionMapping dm =
            DistributionMapping{ba, ParallelDescriptor::NProcs()};

        MakeNewLevelFromScratch(lev, m_time.current_time(), ba, dm);
    }

    m_sim.io_manager().read_checkpoint_fields(
        restart_file, ba_inp, dm_inp, rep);
}
