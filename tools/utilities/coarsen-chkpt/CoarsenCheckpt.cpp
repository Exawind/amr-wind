#include "CoarsenCheckpt.H"
#include "amr-wind/utilities/IOManager.H"

namespace {
void GotoNextLine(std::istream& is)
{
    constexpr std::streamsize bl_ignore_max{100000};
    is.ignore(bl_ignore_max, '\n');
}
} // namespace

namespace amr_wind {
namespace tools {

CoarsenCheckpt::CoarsenCheckpt() : incflo() {}

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

void CoarsenCheckpt::create_twolevelgrid_coarsebase() {}

void CoarsenCheckpt::populate_finest_level()
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

    if (replicate) {
        amrex::Print() << " OLD BA had " << ba_inp[lev0].size() << " GRIDS "
                       << std::endl;
        amrex::Print() << " OLD Domain" << orig_domain << std::endl;
    }

    for (int lev = 0; lev <= finest_level; ++lev) {
        amrex::BoxList bl;
        for (int k = 0; k < rep[2]; k++) {
            for (int j = 0; j < rep[1]; j++) {
                for (int i = 0; i < rep[0]; i++) {
                    for (int nb = 0; nb < ba_inp[lev].size(); nb++) {
                        amrex::Box b(ba_inp[lev][nb]);
                        amrex::IntVect shift_vec(
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
        amrex::BoxArray ba_rep;
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
        dm_inp[lev].define(ba_inp[lev], amrex::ParallelDescriptor::NProcs());

        amrex::BoxArray ba(ba_rep.simplified());
        ba.maxSize(maxGridSize(lev));
        if (refine_grid_layout) {
            ChopGrids(lev, ba, amrex::ParallelDescriptor::NProcs());
        }
        amrex::DistributionMapping dm =
            amrex::DistributionMapping{ba, amrex::ParallelDescriptor::NProcs()};

        MakeNewLevelFromScratch(lev, sim().time().current_time(), ba, dm);
    }

    sim().io_manager().read_checkpoint_fields(
        restart_file, ba_inp, dm_inp, rep);
}

void CoarsenCheckpt::read_chkpt_metadata()
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

    if (replicate) {
        amrex::Print() << " OLD BA had " << ba_inp[lev0].size() << " GRIDS "
                       << std::endl;
        amrex::Print() << " OLD Domain" << orig_domain << std::endl;
    }

    for (int lev = 0; lev <= finest_level; ++lev) {
        amrex::BoxList bl;
        for (int k = 0; k < rep[2]; k++) {
            for (int j = 0; j < rep[1]; j++) {
                for (int i = 0; i < rep[0]; i++) {
                    for (int nb = 0; nb < ba_inp[lev].size(); nb++) {
                        amrex::Box b(ba_inp[lev][nb]);
                        amrex::IntVect shift_vec(
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
        amrex::BoxArray ba_rep;
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
        dm_inp[lev].define(ba_inp[lev], amrex::ParallelDescriptor::NProcs());

        amrex::BoxArray ba(ba_rep.simplified());
        ba.maxSize(maxGridSize(lev));
        if (refine_grid_layout) {
            ChopGrids(lev, ba, amrex::ParallelDescriptor::NProcs());
        }
        amrex::DistributionMapping dm =
            amrex::DistributionMapping{ba, amrex::ParallelDescriptor::NProcs()};

        MakeNewLevelFromScratch(lev, sim().time().current_time(), ba, dm);
    }

    sim().io_manager().read_checkpoint_fields(
        restart_file, ba_inp, dm_inp, rep);
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

    read_chkpt_metadata();

    create_twolevelgrid_coarsebase();

    populate_finest_level();

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

} // namespace tools
} // namespace amr_wind
