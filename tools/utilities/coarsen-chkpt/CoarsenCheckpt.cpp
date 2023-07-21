#include "CoarsenCheckpt.H"
#include "AMReX_MultiFabUtil.H"
#include "AMReX_PlotFileUtil.H"

namespace {
void GotoNextLine(std::istream& is)
{
    constexpr std::streamsize bl_ignore_max{100000};
    is.ignore(bl_ignore_max, '\n');
}
} // namespace

namespace amr_wind {
namespace tools {

CoarsenCheckpt::CoarsenCheckpt() : incflo(), m_io_mgr(new IOManager_Mod(sim()))
{}

void CoarsenCheckpt::ErrorEst(
    int lev, amrex::TagBoxArray& tags, amrex::Real, int)
{
    // Tell AMReX to globally refine the mesh
    if (lev < m_orig_ba.size()) {
        tags.setVal(m_orig_ba[lev], amrex::TagBox::SET);
    } else {
        tags.setVal(amrex::TagBox::SET);
    }
}

/*
void CoarsenCheckpt::read_chkpt_file()
{
    BL_PROFILE("refine-chkpt::read_chkpt_file");
    // Initialize I/O manager to enable restart
    auto& io_mgr = sim().io_manager();
    io_mgr.initialize_io();

    // Ensure that the user has provided a valid checkpoint file
    AMREX_ALWAYS_ASSERT(io_mgr.is_restart());

    ReadCheckpointFile();
    if (amrex::ParallelDescriptor::IOProcessor()) {
        amrex::Print() << "Input grid summary: " << std::endl;
        printGridSummary(amrex::OutStream(), 0, finestLevel());
    }

    m_orig_ba.resize(finestLevel() + 1);
    for (int lev = 0; lev < finestLevel() + 1; ++lev) {
        m_orig_ba[lev] = boxArray(lev);
    }
}
*/

void CoarsenCheckpt::read_chkpt_add_baselevel()
{
    BL_PROFILE("amr-wind::incflo::ReadCheckpointFile()");

    const std::string& restart_file = sim().io_manager().restart_file();
    amrex::Print() << "Coarsening checkpoint " << restart_file << std::endl;

    amrex::Real prob_lo[AMREX_SPACEDIM] = {0.0};
    amrex::Real prob_hi[AMREX_SPACEDIM] = {0.0};

    /***************************************************************************
     * Load header: set up problem domain (including BoxArray)                 *
     *              allocate incflo memory (incflo::AllocateArrays)            *
     *              (by calling MakeNewLevelFromScratch)
     ***************************************************************************/

    std::string File(restart_file + "/Header");

    amrex::VisMF::IO_Buffer io_buffer(amrex::VisMF::GetIOBufferSize());

    amrex::Vector<char> fileCharPtr;
    amrex::ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
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

    sim().time().set_restart_time(nstep, cur_time);

    // Time step size
    is >> dt_restart;
    GotoNextLine(is);

    is >> sim().time().deltaTNm1();
    GotoNextLine(is);

    is >> sim().time().deltaTNm2();
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

    amrex::IntVect rep(1, 1, 1);
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

    bool replicate = (rep != amrex::IntVect::TheUnitVector());

    if (replicate) {
        amrex::Print() << "replicating restart file: " << rep << std::endl;
    }

    // Set up problem domain
    amrex::RealBox rb(prob_lo_input.data(), prob_hi_input.data());
    amrex::Geometry::ResetDefaultProbDomain(rb);
    for (int lev = 0; lev <= max_level; ++lev) {
        SetGeometry(
            lev, amrex::Geometry(
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
    amrex::Box orig_domain(ba_inp[lev0].minimalBox());

    // create base level BoxArray

    // create base level from scratch
    // MakeNewLevelFromScratch(0, sim().time().current_time(), ba, dm);

    for (int levnew = 1; levnew <= finest_level + 1; ++levnew) {
        amrex::BoxList bl;
        int levold = levnew - 1;
        for (int k = 0; k < rep[2]; k++) {
            for (int j = 0; j < rep[1]; j++) {
                for (int i = 0; i < rep[0]; i++) {
                    for (int nb = 0; nb < ba_inp[levold].size(); nb++) {
                        amrex::Box b(ba_inp[levold][nb]);
                        amrex::IntVect shift_vec(
                            i * orig_domain.length(0),
                            j * orig_domain.length(1),
                            k * orig_domain.length(2));

                        // equivalent to 2^lev
                        shift_vec *= (1 << levold);

                        b.shift(shift_vec);
                        bl.push_back(b);
                    }
                }
            }
        }
        amrex::BoxArray ba_rep;
        ba_rep.define(bl);

        // Create distribution mapping
        dm_inp[levold].define(
            ba_inp[levold], amrex::ParallelDescriptor::NProcs());

        amrex::BoxArray ba(ba_rep.simplified());
        ba.maxSize(maxGridSize(levnew));

        amrex::DistributionMapping dm =
            amrex::DistributionMapping{ba, amrex::ParallelDescriptor::NProcs()};

        MakeNewLevelFromScratch(levnew, sim().time().current_time(), ba, dm);
    }

    io_manager().read_checkpoint_fields_offset(
        restart_file, ba_inp, dm_inp, rep, 1);
}

void CoarsenCheckpt::coarsen_chkpt_file()
{
    amrex::Print() << "Uniformly refining mesh" << std::endl;
    amrex::Real rstart = amrex::ParallelDescriptor::second();
    regrid(0, sim().time().current_time());
    amrex::Real rend = amrex::ParallelDescriptor::second();
    amrex::Print() << "Refinement time elapsed: " << (rend - rstart)
                   << std::endl;
    if (amrex::ParallelDescriptor::IOProcessor()) {
        amrex::Print() << "Refined grid summary: " << std::endl;
        printGridSummary(amrex::OutStream(), 0, finestLevel());
    }

    // Initialize modified io_manager
    read_chkpt_add_baselevel();
    // average down
}

void CoarsenCheckpt::run_utility()
{

    if ((maxLevel() - finestLevel()) != 1) {
        amrex::Print() << "Checkpoint refinement only supported for one "
                          "refinement level at a time\n"
                       << "  Current finest level = " << finestLevel()
                       << "\n  Requested max level  = " << maxLevel()
                       << "\n  Unable to refine checkpoint file" << std::endl;
    } else {
        coarsen_chkpt_file();
        const int start_level = 1;
        amrex::Print() << "Writing refined levels: " << start_level << " - "
                       << finestLevel() << std::endl;
        sim().io_manager().write_checkpoint_file(start_level);
    }
    // write checkpoint file, only level 0
}

IOManager_Mod::IOManager_Mod(CFDSim& sim) : m_sim(sim), IOManager(sim) {}

void IOManager_Mod::read_checkpoint_fields_offset(
    const std::string& restart_file,
    const amrex::Vector<amrex::BoxArray>& ba_chk,
    const amrex::Vector<amrex::DistributionMapping>& dm_chk,
    const amrex::IntVect& rep,
    const int off)
{
    BL_PROFILE("amr-wind::IOManager::read_checkpoint_fields");

    // Track set of fields that might be missing at this level
    std::set<std::string> missing;
    const std::string level_prefix = "Level_";
    const int nlevels = m_sim.mesh().finestLevel() + 1;

    // always use the level 0 domain
    amrex::Box orig_domain(ba_chk[0].minimalBox());

    for (int lev = 0; lev < nlevels; ++lev) {
        for (auto* fld : m_chk_fields) {
            auto& field = *fld;
            const auto& fab_file = amrex::MultiFabFileFullPrefix(
                lev, restart_file, level_prefix, field.name());

            // Fields might be registered for checkpoint but might not be
            // necessary for actually performing the simulation. Check if the
            // field exists before attempting to read the restart field.
            if (!amrex::VisMF::Exist(fab_file)) {
                missing.insert(field.name());
                continue;
            }

            auto& mfab = field(lev);
            const auto& ba_fab = amrex::convert(ba_chk[lev], mfab.ixType());
            if (mfab.boxArray() == ba_fab &&
                mfab.DistributionMap() == dm_chk[lev]) {
                amrex::VisMF::Read(
                    field(lev),
                    amrex::MultiFabFileFullPrefix(
                        lev, restart_file, level_prefix, field.name()));
            } else {
                amrex::MultiFab tmp(
                    ba_fab, dm_chk[lev], mfab.nComp(), mfab.nGrowVect());
                amrex::VisMF::Read(
                    tmp, amrex::MultiFabFileFullPrefix(
                             lev, restart_file, level_prefix, field.name()));

                for (int k = 0; k < rep[2]; k++) {
                    for (int j = 0; j < rep[1]; j++) {
                        for (int i = 0; i < rep[0]; i++) {

                            amrex::IntVect shift_vec(
                                i * orig_domain.length(0),
                                j * orig_domain.length(1),
                                k * orig_domain.length(2));

                            // equivalent to 2^lev
                            shift_vec *= (1 << lev);

                            tmp.shift(shift_vec);
                            mfab.ParallelCopy(tmp);
                            tmp.shift(-shift_vec);
                        }
                    }
                }

                mfab.setBndry(0.0);
            }
        }
    }

    // If fields were missing, print diagnostic message.
    if (!missing.empty()) {
        amrex::Print() << "\nWARNING: The following fields were missing in the "
                          "restart file for one or more levels. Please check "
                          "your restart file and inputs."
                       << std::endl
                       << "Missing checkpoint fields: " << std::endl;
        for (const auto& ff : missing) {
            amrex::Print() << "  - " << ff << std::endl;
        }
        amrex::Print() << std::endl;
        if (!m_allow_missing_restart_fields) {
            amrex::Abort("Missing fields in restart file.");
        }
    }
}

} // namespace tools
} // namespace amr_wind
