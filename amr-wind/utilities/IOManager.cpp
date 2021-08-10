#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>
#include <chrono>
#include <ctime>
#include <fstream>
#include <netcdf.h>

#include "amr-wind/utilities/IOManager.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/utilities/console_io.H"
#include "amr-wind/utilities/io_utils.H"
#include "amr-wind/utilities/DerivedQuantity.H"
#include "amr-wind/utilities/DerivedQtyDefs.H"
#include "amr-wind/utilities/ncutils/nc_interface.H"

#include "AMReX_ParmParse.H"
#include "AMReX_PlotFileUtil.H"
#include "AMReX_MultiFabUtil.H"

namespace amr_wind {

IOManager::IOManager(CFDSim& sim)
    : m_sim(sim), m_derived_mgr(new DerivedQtyMgr(m_sim.repo()))
{}

IOManager::~IOManager() = default;

void IOManager::initialize_io()
{
    amrex::Vector<std::string> out_vars;
    amrex::Vector<std::string> out_int_vars;
    amrex::Vector<std::string> out_derived_vars;
    amrex::Vector<std::string> out_skip_vars;
    std::set<std::string> outputs;
    std::set<std::string> skip_outputs;
    std::set<std::string> int_outputs;

    amrex::ParmParse pp("io");
    pp.query("output_default_variables", m_output_default_vars);
    pp.query("plot_file", m_plt_prefix);
    pp.query("check_file", m_chk_prefix);
    pp.query("restart_file", m_restart_file);
    pp.query("allow_missing_restart_fields", m_allow_missing_restart_fields);

    // ParmParse requires us to read in a vector
    pp.queryarr("outputs", out_vars);
    pp.queryarr("int_outputs", out_int_vars);
    pp.queryarr("derived_outputs", out_derived_vars);
    pp.queryarr("skip_outputs", out_skip_vars);

    // We process the input vector to eliminate duplicates
    for (const auto& name : out_vars) {
        outputs.insert(name);
    }
    for (const auto& name : out_skip_vars) {
        skip_outputs.insert(name);
    }
    for (const auto& name : out_int_vars) {
        int_outputs.insert(name);
    }

    // If the user hasn't disabled default output variables, then we append them
    // to the list so that we don't end up with any duplicates (in case user
    // also added the variable explicitly in the input file)
    if (m_output_default_vars) {
        for (const auto& name : m_pltvars_default) {
            if (skip_outputs.find(name) == skip_outputs.end()) {
                outputs.insert(name);
            }
        }

        for (const auto& name : m_int_pltvars_default) {
            int_outputs.insert(name);
        }
    }

    amrex::Print() << "Initializing I/O manager" << std::endl;

    // Process output variables information
    auto& repo = m_sim.repo();
    m_plt_num_comp = 0;
    for (const auto& fname : outputs) {
        if (repo.field_exists(fname)) {
            auto& fld = repo.get_field(fname);
            m_plt_num_comp += fld.num_comp();
            m_plt_fields.emplace_back(&fld);
            ioutils::add_var_names(m_plt_var_names, fld.name(), fld.num_comp());
        } else {
            amrex::Print() << "  Invalid output variable requested: " << fname
                           << std::endl;
        }
    }

    for (const auto& fname : int_outputs) {
        if (repo.int_field_exists(fname)) {
            auto& fld = repo.get_int_field(fname);
            m_plt_num_comp += fld.num_comp();
            m_int_plt_fields.emplace_back(&fld);
            ioutils::add_var_names(m_plt_var_names, fld.name(), fld.num_comp());
        } else {
            amrex::Print() << "  Invalid output variable requested: " << fname
                           << std::endl;
        }
    }

    if (!out_derived_vars.empty()) {
        m_derived_mgr->create(out_derived_vars);
        m_plt_num_comp += m_derived_mgr->num_comp();
        m_derived_mgr->var_names(m_plt_var_names);
    }

    for (const auto& fname : m_chkvars) {
        auto& fld = repo.get_field(fname);
        m_chk_fields.emplace_back(&fld);
    }
}

void IOManager::write_plot_file()
{
    BL_PROFILE("amr-wind::IOManager::write_plot_file");

    amrex::Vector<int> istep(
        m_sim.mesh().finestLevel() + 1, m_sim.time().time_index());
    const int plt_comp = m_plt_num_comp;
    const int start_comp = m_plt_num_comp - m_derived_mgr->num_comp();
    auto outfield = m_sim.repo().create_scratch_field(plt_comp);
    const int nlevels = m_sim.repo().num_active_levels();

    for (int lev = 0; lev < nlevels; ++lev) {
        int icomp = 0;
        auto& mf = (*outfield)(lev);

        for (auto* fld : m_plt_fields) {
            amrex::MultiFab::Copy(
                mf, (*fld)(lev), 0, icomp, fld->num_comp(), 0);
            icomp += fld->num_comp();
        }

        for (auto* fld : m_int_plt_fields) {
            amrex::MultiFab::Copy(
                mf, amrex::ToMultiFab((*fld)(lev)), 0, icomp, fld->num_comp(),
                0);
            icomp += fld->num_comp();
        }
    }

    (*m_derived_mgr)(*outfield, start_comp);

    const std::string& plt_filename =
        amrex::Concatenate(m_plt_prefix, m_sim.time().time_index());
    const auto& mesh = m_sim.mesh();
    amrex::Print() << "Writing plot file       " << plt_filename << " at time "
                   << m_sim.time().new_time() << std::endl;
    amrex::WriteMultiLevelPlotfile(
        plt_filename, nlevels, outfield->vec_const_ptrs(), m_plt_var_names,
        mesh.Geom(), m_sim.time().new_time(), istep, mesh.refRatio());

    write_info_file(plt_filename);
}

void IOManager::write_checkpoint_file(const int start_level)
{
    BL_PROFILE("amr-wind::IOManager::write_checkpoint_file");
    const std::string level_prefix = "Level_";
    const std::string chkname =
        amrex::Concatenate(m_chk_prefix, m_sim.time().time_index());

    amrex::Print() << "Writing checkpoint file " << chkname << " at time "
                   << m_sim.time().new_time() << std::endl;
    const auto& mesh = m_sim.mesh();
    amrex::PreBuildDirectorHierarchy(
        chkname, level_prefix, mesh.finestLevel() + 1 - start_level, true);
    write_header(chkname, start_level);
    write_info_file(chkname);

    for (int lev = start_level; lev < mesh.finestLevel() + 1; ++lev) {
        for (auto* fld : m_chk_fields) {
            auto& field = *fld;
            amrex::VisMF::Write(
                field(lev),
                amrex::MultiFabFileFullPrefix(
                    lev - start_level, chkname, level_prefix, field.name()));
        }
    }
}

void IOManager::read_checkpoint_fields(
    const std::string& restart_file,
    const amrex::Vector<amrex::BoxArray>& ba_chk,
    const amrex::Vector<amrex::DistributionMapping>& dm_chk,
    const amrex::IntVect& rep)
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

void IOManager::write_header(const std::string& chkname, const int start_level)
{
    if (!amrex::ParallelDescriptor::IOProcessor()) {
        return;
    }

    const std::string hdr_name(chkname + "/Header");
    amrex::VisMF::IO_Buffer io_buf(amrex::VisMF::IO_Buffer_Size);

    std::ofstream hdr;
    hdr.rdbuf()->pubsetbuf(io_buf.dataPtr(), io_buf.size());

    hdr.open(
        hdr_name.c_str(),
        std::ofstream::out | std::ofstream::trunc | std::ofstream::binary);
    if (!hdr.good()) {
        amrex::FileOpenFailed(hdr_name);
    }

    hdr.precision(17);

    const auto& mesh = m_sim.mesh();
    const auto& time = m_sim.time();
    hdr << "Checkpoint version: 1\n"
        << mesh.finestLevel() - start_level << "\n"
        << time.time_index() << "\n"
        << time.new_time() << "\n"
        << time.deltaT() << "\n"
        << time.deltaTNm1() << "\n"
        << time.deltaTNm2() << "\n";

    const auto geom = mesh.Geom(0);
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        hdr << geom.ProbLo(i) << " ";
    }
    hdr << "\n";
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        hdr << geom.ProbHi(i) << " ";
    }
    hdr << "\n";

    for (int lev = start_level; lev < mesh.finestLevel() + 1; ++lev) {
        mesh.boxArray(lev).writeOn(hdr);
        hdr << "\n";
    }

    hdr.close();
}

void IOManager::write_info_file(const std::string& path)
{
    if (!amrex::ParallelDescriptor::IOProcessor()) {
        return;
    }

    const std::string dash_line = "\n" + std::string(78, '-') + "\n";
    const std::string fname(path + "/amr_wind_info");
    std::ofstream fh(fname.c_str(), std::ios::out);
    if (!fh.good()) {
        amrex::FileOpenFailed(fname);
    }

    amr_wind::io::print_banner(amrex::ParallelContext::CommunicatorSub(), fh);

    fh << dash_line << "Grid information: " << std::endl;
    const auto& mesh = m_sim.mesh();
    for (int lev = 0; lev < m_sim.mesh().finestLevel() + 1; ++lev) {
        fh << "  Level: " << lev << "\n"
           << "    num. boxes = " << mesh.boxArray().size() << "\n"
           << "    maximum zones = ";

        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            fh << mesh.Geom(lev).Domain().length(dir) << " ";
        }
        fh << "\n";
    }

    fh << dash_line << "Input file parameters: " << std::endl;
    amrex::ParmParse::dumpTable(fh, true);
    fh.close();
}

void IOManager::write_netcdf_file(){

#ifdef AMR_WIND_USE_NETCDF
  BL_PROFILE("amr-wind::IOManager::write_netcdf_file");

  amrex::ParmParse pp("io");

  pp.query("ncf_file", m_netcdf_prefix);

  const std::string& ncf_filename =
      amrex::Concatenate(m_netcdf_prefix, m_sim.time().time_index());

  auto ncf = ncutils::NCFile::create(ncf_filename, NC_CLOBBER | NC_NETCDF4);

  ncf.enter_def_mode();

  amrex::Vector<amrex::Real> factormap; 
  amrex::ParmParse ppMap("ConstantScaling");
  ppMap.queryarr("scaling_factor",factormap,0,AMREX_SPACEDIM);

  // Get normal direction and associated stuff
  const auto& mesh = m_sim.mesh();
  const auto& geom = mesh.Geom()[0];
  amrex::Box const& domain = geom.Domain();
  const auto dlo = amrex::lbound(domain);
  const auto dhi = amrex::ubound(domain);
  const auto* problo = mesh.Geom(0).ProbLo();
  const amrex::Real* dx = mesh.Geom(0).CellSize();
  
  ncf.def_dim("xc", dhi.x - dlo.x + 1);
  ncf.def_dim("yc", dhi.y - dlo.y + 1);
  ncf.def_dim("zc", dhi.z - dlo.z + 1);

  auto xcell = ncf.def_var("xc", NC_DOUBLE, {"xc"});
  auto ycell = ncf.def_var("yc", NC_DOUBLE, {"yc"});
  auto zcell = ncf.def_var("zc", NC_DOUBLE, {"zc"});

  xcell.put_attr("units","m");
  xcell.put_attr("axis", "X");

  ycell.put_attr("units","m");
  ycell.put_attr("axis", "Y");

  zcell.put_attr("units","m");
  zcell.put_attr("axis", "Z");
  
  const std::vector<std::string> dim_3_c{"xc", "yc", "zc"};
  
  auto xvelncf = ncf.def_array("x_velocity", NC_DOUBLE, dim_3_c);
  auto yvelncf = ncf.def_array("y_velocity", NC_DOUBLE, dim_3_c);
  auto zvelncf = ncf.def_array("z_velocity", NC_DOUBLE, dim_3_c);

  ncf.exit_def_mode();
  
  std::vector<double> fill_val_x(dhi.x - dlo.x + 1);
  std::vector<double> fill_val_y(dhi.y - dlo.y + 1);
  std::vector<double> fill_val_z(dhi.z - dlo.z + 1);

  for (int i = dlo.x; i<= dhi.x; i++) {
    fill_val_x[i] = problo[0] + dx[0]*(static_cast<amrex::Real>(i)+0.5)*factormap[0];  // perform mapping to the non-uniform space
  }
  xcell.put(fill_val_x.data());

  for (int j = dlo.y; j<= dhi.y; j++) {
    fill_val_y[j] = problo[1] + dx[1]*(static_cast<amrex::Real>(j)+0.5)*factormap[1];  // perform mapping to the non-uniform space
  }
  ycell.put(fill_val_y.data());
  for (int k = dlo.z; k<= dhi.z; k++) {
    fill_val_z[k] = problo[2] + dx[2]*(static_cast<amrex::Real>(k)+0.5)*factormap[2];  // perform mapping to the non-uniform space
  }
  zcell.put(fill_val_z.data());

  auto& repo = m_sim.repo();
  auto& velocity = repo.get_field("velocity");
  
  if (velocity.is_mesh_mapped()) {
        velocity.to_unmapped_mesh();
  }

  auto velcomp_x = repo.create_scratch_field(1);
  auto velcomp_y = repo.create_scratch_field(1);
  auto velcomp_z = repo.create_scratch_field(1);

  amrex::MultiFab::Copy((*velcomp_x)(0), velocity(0), 0, 0, 1, 0);
  amrex::MultiFab::Copy((*velcomp_y)(0), velocity(0), 1, 0, 1, 0);
  amrex::MultiFab::Copy((*velcomp_z)(0), velocity(0), 2, 0, 1, 0);
  
  amrex::MFItInfo mfi_info;
  for (amrex::MFIter mfi(velocity(0), mfi_info); mfi.isValid();
       ++mfi) {

    amrex::Box const& bx = mfi.validbox();
    auto xvel = (*velcomp_x)(0).const_array(mfi);
    auto yvel = (*velcomp_y)(0).const_array(mfi);
    auto zvel = (*velcomp_z)(0).const_array(mfi);
  
    const auto lo = lbound(bx);
    const auto hi = ubound(bx);

    amrex::Vector<size_t> start = {static_cast<unsigned long>(lo.x), static_cast<unsigned long>(lo.y), static_cast<unsigned long>(lo.z)};
    amrex::Vector<size_t> count = {static_cast<unsigned long>(hi.x-lo.x+1), static_cast<unsigned long>(hi.y-lo.y+1), static_cast<unsigned long>(hi.z-lo.z+1)};
    xvelncf.put(xvel.dataPtr(), start, count);
    yvelncf.put(yvel.dataPtr(), start, count);
    zvelncf.put(zvel.dataPtr(), start, count);
  }
  
  ncf.close();
#endif
}

} // namespace amr_wind
