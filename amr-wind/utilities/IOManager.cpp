#include <chrono>
#include <ctime>

#include "amr-wind/utilities/IOManager.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/utilities/console_io.H"

#include "AMReX_ParmParse.H"
#include "AMReX_PlotFileUtil.H"

namespace amr_wind {
namespace {

/** Generate field names for multi-component fields
 *
 *  Function generates variable names based on the number of components in the
 *  field. If the field is a scalar, it just inserts the variable name. For
 *  vectors it appends "x, y, z" to the components, and for all other fields it
 *  numbers the components start with 0.
 */
void add_var_names(
    amrex::Vector<std::string>& vnames,
    const std::string& fname,
    const int ncomp)
{
    const amrex::Vector<std::string> comp{"x", "y", "z"};

    switch (ncomp) {
    case 1:
        vnames.push_back(fname);
        break;

    case AMREX_SPACEDIM:
        for (int i = 0; i < AMREX_SPACEDIM; ++i)
            vnames.push_back(fname + comp[i]);
        break;

    default:
        for (int i = 0; i < ncomp; ++i)
            vnames.push_back(fname + std::to_string(i));
    }
}
} // namespace

IOManager::IOManager(CFDSim& sim) : m_sim(sim) {}

IOManager::~IOManager() = default;

void IOManager::initialize_io()
{
    amrex::Vector<std::string> out_vars;
    std::set<std::string> outputs;

    amrex::ParmParse pp("io");
    pp.query("output_default_variables", m_output_default_vars);
    pp.query("plot_file", m_plt_prefix);
    pp.query("check_file", m_chk_prefix);
    pp.query("restart_file", m_restart_file);

    // ParmParse requires us to read in a vector
    pp.queryarr("outputs", out_vars);

    // We process the input vector to eliminate duplicates
    for (const auto& name: out_vars)
        outputs.insert(name);

    // If the user hasn't disabled default output variables, then we append them
    // to the list so that we don't end up with any duplicates (in case user
    // also added the variable explicitly in the input file)
    if (m_output_default_vars) {
        for (const auto& name: m_pltvars_default)
            outputs.insert(name);
    }

    amrex::Print() << "Initializing I/O manager" << std::endl;

    // Process output variables information
    auto& repo = m_sim.repo();
    m_plt_num_comp = 0;
    for (const auto& fname: outputs) {
        if (repo.field_exists(fname)) {
            auto& fld = repo.get_field(fname);
            m_plt_num_comp += fld.num_comp();
            m_plt_fields.emplace_back(&fld);
            add_var_names(m_plt_var_names, fld.name(), fld.num_comp());
        } else {
            amrex::Print() << "  Invalid output variable requested: " << fname << std::endl;
        }
    }

    for (const auto& fname: m_chkvars) {
        auto& fld = repo.get_field(fname);
        m_chk_fields.emplace_back(&fld);
    }
}

void IOManager::write_plot_file()
{
    BL_PROFILE("amr-wind::IOManager::write_plot_file");

    amrex::Vector<int> istep(m_sim.mesh().finestLevel() + 1, m_sim.time().time_index());
    auto outfield = m_sim.repo().create_scratch_field(m_plt_num_comp);
    const int nlevels = m_sim.repo().num_active_levels();

    for (int lev=0; lev < nlevels; ++lev) {
        int icomp = 0;
        auto& mf = (*outfield)(lev);

        for (auto* fld: m_plt_fields) {
            amrex::MultiFab::Copy(mf, (*fld)(lev), 0, icomp, fld->num_comp(), 0);
            icomp += fld->num_comp();
        }
    }

    const std::string& plt_filename =
        amrex::Concatenate(m_plt_prefix, m_sim.time().time_index());
    const auto& mesh = m_sim.mesh();
    amrex::Print()
        << "Writing plot file       " << plt_filename << " at time "
        << m_sim.time().new_time() << std::endl;
    amrex::WriteMultiLevelPlotfile(
        plt_filename, nlevels, outfield->vec_const_ptrs(), m_plt_var_names,
        mesh.Geom(), m_sim.time().new_time(), istep, mesh.refRatio());

    write_info_file(plt_filename);
}

void IOManager::write_checkpoint_file()
{
    BL_PROFILE("amr-wind::IOManager::write_checkpoint_file");
    const std::string level_prefix = "Level_";
    const std::string chkname = amrex::Concatenate(m_chk_prefix, m_sim.time().time_index());

    amrex::Print() << "Writing checkpoint file " << chkname << " at time "
                   << m_sim.time().new_time() << std::endl;
    const auto& mesh = m_sim.mesh();
    amrex::PreBuildDirectorHierarchy(chkname, level_prefix, mesh.finestLevel() + 1, true);
    write_header(chkname);
    write_info_file(chkname);

    for (int lev = 0; lev < mesh.finestLevel() + 1; ++lev) {
        for (auto* fld: m_chk_fields) {
            auto& field = *fld;
            amrex::VisMF::Write(
                field(lev),
                amrex::MultiFabFileFullPrefix(
                    lev, chkname, level_prefix, field.name()));
        }
    }
}

void IOManager::read_checkpoint_fields(const std::string& restart_file)
{
    BL_PROFILE("amr-wind::IOManager::read_checkpoint_fields");
    const std::string level_prefix = "Level_";
    const int nlevels = m_sim.mesh().finestLevel() + 1;
    for (int lev = 0; lev < nlevels; ++lev) {
        for (auto* fld: m_chk_fields) {
            auto& field = *fld;
            amrex::VisMF::Read(
                field(lev),
                amrex::MultiFabFileFullPrefix(
                    lev, restart_file, level_prefix, field.name()));
        }
    }
}

void IOManager::write_header(const std::string& chkname)
{
    if (!amrex::ParallelDescriptor::IOProcessor()) return;

    const std::string hdr_name(chkname + "/Header");
    amrex::VisMF::IO_Buffer io_buf(amrex::VisMF::IO_Buffer_Size);

    std::ofstream hdr;
    hdr.rdbuf()->pubsetbuf(io_buf.dataPtr(), io_buf.size());

    hdr.open(hdr_name.c_str(),
             std::ofstream::out | std::ofstream::trunc | std::ofstream::binary);
    if (!hdr.good()) {
        amrex::FileOpenFailed(hdr_name);
    }

    hdr.precision(17);

    const auto& mesh = m_sim.mesh();
    const auto& time = m_sim.time();
    hdr << "Checkpoint version: 1\n"
        << mesh.finestLevel() << "\n"
        << time.time_index() << "\n"
        << time.new_time() << "\n"
        << time.deltaT() << "\n"
        << time.deltaTNm1() << "\n"
        << time.deltaTNm2() << "\n";

    const auto geom = mesh.Geom(0);
    for (int i = 0; i < AMREX_SPACEDIM; ++i) hdr << geom.ProbLo(i) << " ";
    hdr << "\n";
    for (int i = 0; i < AMREX_SPACEDIM; ++i) hdr << geom.ProbHi(i) << " ";
    hdr << "\n";

    for (int lev = 0; lev < mesh.finestLevel() + 1; ++lev) {
        mesh.boxArray(lev).writeOn(hdr);
        hdr << "\n";
    }

    hdr.close();
}

void IOManager::write_info_file(const std::string& path)
{
    if (!amrex::ParallelDescriptor::IOProcessor()) return;

    const std::string dash_line = "\n" + std::string(78, '-') + "\n";
    const std::string fname(path + "/amr_wind_info");
    std::ofstream fh(fname.c_str(), std::ios::out);
    if (!fh.good()) {
        amrex::FileOpenFailed(fname);
    }

    amr_wind::io::print_banner(fh);

    fh << dash_line << "Grid information: " << std::endl;
    const auto& mesh = m_sim.mesh();
    for (int lev=0; lev < m_sim.mesh().finestLevel()+1; ++lev) {
        fh << "  Level: " << lev << "\n"
           << "    num. boxes = " << mesh.boxArray().size() << "\n"
           << "    maximum zones = ";

        for (int dir=0; dir < AMREX_SPACEDIM; ++dir)
            fh << mesh.Geom(lev).Domain().length(dir) << " ";
        fh << "\n";
    }

    fh << dash_line << "Input file parameters: " << std::endl;
    amrex::ParmParse::dumpTable(fh, true);
    fh.close();
}

} // namespace amr_wind
