#include "amr-wind/utilities/sampling/FieldNorms.H"
#include "amr-wind/utilities/io_utils.H"
#include "amr-wind/utilities/ncutils/nc_interface.H"
#include <AMReX_MultiFabUtil.H>
#include "AMReX_ParmParse.H"

namespace amr_wind {
namespace field_norms {

FieldNorms::FieldNorms(CFDSim& sim, const std::string& label)
    : m_sim(sim), m_label(label)
{}

FieldNorms::~FieldNorms() = default;

void FieldNorms::initialize()
{
    BL_PROFILE("amr-wind::FieldNorms::initialize");

    // Labels for the different sampler types
    amrex::Vector<std::string> labels;
    // Fields to be sampled - requested by user
    amrex::Vector<std::string> field_names;

    bool main_fields = false;
    bool all_fields = false;

    {
        amrex::ParmParse pp(m_label);
        pp.queryarr("fields", field_names);
        pp.query("output_frequency", m_out_freq);
        pp.query("output_format", m_out_fmt);
        pp.query("main_equation_fields", main_fields);
        pp.query("all_fields", all_fields);
    }

    if (field_names.size() == 0) {
        main_fields = true;
    }

    if (all_fields) {
        for (const auto& fld : m_sim.repo().fields()) {
            if ((*fld).field_location() == amr_wind::FieldLoc::CELL) {
                field_names.emplace_back((*fld).name());
            }
        }

    } else if (main_fields) {
        field_names.emplace_back("density");
        auto& pde_mgr = m_sim.pde_manager();
        {
            auto& vel = pde_mgr.icns().fields().field;
            field_names.emplace_back(vel.name());
        }
        for (auto& eqn : pde_mgr.scalar_eqns()) {
            auto& fld = eqn->fields().field;
            field_names.emplace_back(fld.name());
        }
    }

    // Process field information
    int ncomp = 0;
    auto& repo = m_sim.repo();
    for (const auto& fname : field_names) {
        if (!repo.field_exists(fname)) {
            amrex::Print()
                << "WARNING: FieldNorms: Non-existent field requested: "
                << fname << std::endl;
            continue;
        }

        auto& fld = repo.get_field(fname);

        if (fld.field_location() != amr_wind::FieldLoc::CELL) {
            amrex::Print() << "WARNING: FieldNorms: Only works for cell based "
                              "fields, skipping: "
                           << fld.name() << std::endl;
            continue;
        }

        ncomp += fld.num_comp();
        m_fields.emplace_back(&fld);
        ioutils::add_var_names(m_var_names, fld.name(), fld.num_comp());
    }

    m_fnorms.resize(m_var_names.size(), 0.0);

    if (m_out_fmt == "netcdf") prepare_netcdf_file();
    if (m_out_fmt == "ascii") prepare_ascii_file();
}

amrex::Real FieldNorms::L2_norm(amr_wind::Field& field, const int comp)
{
    amrex::Real nrm = 0.0;

    AMREX_ASSERT(comp >= 0 && comp < field.num_comp());

    const int finest_level = field.repo().num_active_levels() - 1;
    const auto& geom = field.repo().mesh().Geom();
    const auto& dmap = field.repo().mesh().DistributionMap();
    const auto& grids = field.repo().mesh().boxArray();
    const int nghost = 0;

    for (int lev = 0; lev <= finest_level; lev++) {

        amrex::iMultiFab level_mask;
        if (lev < finest_level) {
            level_mask = amrex::makeFineMask(
                grids[lev], dmap[lev], grids[lev + 1], amrex::IntVect(2), 1, 0);
        } else {
            level_mask.define(grids[lev], dmap[lev], 1, 0, amrex::MFInfo());
            level_mask.setVal(1);
        }

        const amrex::Real cell_vol = geom[lev].CellSize()[0] *
                                     geom[lev].CellSize()[1] *
                                     geom[lev].CellSize()[2];

        nrm += amrex::ReduceSum(
            field(lev), level_mask, nghost,
            [=] AMREX_GPU_HOST_DEVICE(
                amrex::Box const& bx,
                amrex::Array4<amrex::Real const> const& field_arr,
                amrex::Array4<int const> const& mask_arr) -> amrex::Real {
                amrex::Real nrm_fab = 0.0;

                amrex::Loop(bx, [=, &nrm_fab](int i, int j, int k) noexcept {
                    nrm_fab += cell_vol * mask_arr(i, j, k) *
                               field_arr(i, j, k, comp) *
                               field_arr(i, j, k, comp);
                });
                return nrm_fab;
            });
    }

    amrex::ParallelDescriptor::ReduceRealSum(nrm);

    const amrex::Real total_volume = geom[0].ProbDomain().volume();

    return std::sqrt(nrm / total_volume);
}

void FieldNorms::process_field_norms()
{
    int ind = 0;
    for (int field = 0; field < m_fields.size(); ++field) {
        const auto fld = m_fields[field];
        for (int comp = 0; comp < fld->num_comp(); ++comp) {
            m_fnorms[ind++] = L2_norm(*fld, comp);
        }
    }
}

void FieldNorms::post_advance_work()
{
    BL_PROFILE("amr-wind::FieldNorms::post_advance_work");
    const auto& time = m_sim.time();
    const int tidx = time.time_index();
    // Skip processing if it is not an output timestep
    if (!(tidx % m_out_freq == 0)) return;

    process_field_norms();
    process_output();
}

void FieldNorms::process_output()
{
    if (m_out_fmt == "native") {
        impl_write_native();
    } else if (m_out_fmt == "ascii") {
        write_ascii();
    } else if (m_out_fmt == "netcdf") {
        write_netcdf();
    } else {
        amrex::Abort("FieldNorms: Invalid output format encountered");
    }
}

void FieldNorms::impl_write_native()
{
    BL_PROFILE("amr-wind::FieldNorms::write_native");

    const std::string post_dir = "post_processing";
    const std::string sdir =
        amrex::Concatenate(m_label, m_sim.time().time_index());
    amrex::Vector<std::string> int_var_names{"uid", "set_id", "probe_id"};

    amrex::Abort("native field norm not implemented yet");
}

void FieldNorms::prepare_ascii_file()
{
    BL_PROFILE("amr-wind::FieldNorms::prepare_ascii_file");
    amrex::Print()
        << "WARNING: FieldNorms: ASCII output will impact performance"
        << std::endl;

    const std::string post_dir = "post_processing";
    const std::string sname =
        amrex::Concatenate(m_label, m_sim.time().time_index());

    if (!amrex::UtilCreateDirectory(post_dir, 0755)) {
        amrex::CreateDirectoryFailed(post_dir);
    }
    m_out_fname = post_dir + "/" + sname + ".txt";

    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ofstream f;
        f.open(m_out_fname.c_str());
        f << "time_step "
          << "time";
        for (int i = 0; i < m_var_names.size(); ++i) f << ' ' << m_var_names[i];
        f << std::endl;
        f.close();
    }
}

void FieldNorms::write_ascii()
{
    BL_PROFILE("amr-wind::FieldNorms::write_ascii");

    const std::string post_dir = "post_processing";

    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ofstream f;
        f.open(m_out_fname.c_str(), std::ios_base::app);
        f << m_sim.time().time_index() << std::fixed
          << std::setprecision(m_precision) << std::setw(m_width)
          << m_sim.time().new_time();
        for (int i = 0; i < m_fnorms.size(); ++i)
            f << std::setw(m_width) << m_fnorms[i];
        f << std::endl;
        f.close();
    }
}

void FieldNorms::prepare_netcdf_file()
{
#ifdef AMR_WIND_USE_NETCDF

    amrex::Abort("netcdf field norm not implemented yet \n");

    const std::string post_dir = "post_processing";
    const std::string sname =
        amrex::Concatenate(m_label, m_sim.time().time_index());
    if (!amrex::UtilCreateDirectory(post_dir, 0755)) {
        amrex::CreateDirectoryFailed(post_dir);
    }
    m_ncfile_name = post_dir + "/" + sname + ".nc";

    // Only I/O processor handles NetCDF generation
    if (!amrex::ParallelDescriptor::IOProcessor()) return;

    auto ncf = ncutils::NCFile::create(m_ncfile_name, NC_CLOBBER | NC_NETCDF4);
    const std::string nt_name = "num_time_steps";
    const std::string npart_name = "num_points";
    const std::vector<std::string> two_dim{nt_name, npart_name};
    ncf.enter_def_mode();
    ncf.put_attr("title", "AMR-Wind data sampling output");
    ncf.put_attr("version", ioutils::amr_wind_version());
    ncf.put_attr("created_on", ioutils::timestamp());
    ncf.def_dim(nt_name, NC_UNLIMITED);
    ncf.def_dim("ndim", AMREX_SPACEDIM);
    ncf.def_var("time", NC_DOUBLE, {nt_name});
    // Define groups for each sampler
    //    for (const auto& obj : m_samplers) {
    //        auto grp = ncf.def_group(obj->label());
    //
    //        grp.def_dim(npart_name, obj->num_points());
    //        obj->define_netcdf_metadata(grp);
    //        grp.def_var("coordinates", NC_DOUBLE, {npart_name, "ndim"});
    //        for (const auto& vname : m_var_names)
    //            grp.def_var(vname, NC_DOUBLE, two_dim);
    //    }
    ncf.exit_def_mode();
//
//    {
//        const std::vector<size_t> start{0, 0};
//        std::vector<size_t> count{0, AMREX_SPACEDIM};
//        SamplerBase::SampleLocType locs;
//        for (const auto& obj : m_samplers) {
//            auto grp = ncf.group(obj->label());
//            obj->populate_netcdf_metadata(grp);
//            obj->sampling_locations(locs);
//            auto xyz = grp.var("coordinates");
//            count[0] = obj->num_points();
//            xyz.put(&locs[0][0], start, count);
//        }
//    }
#else
    amrex::Abort(
        "NetCDF support was not enabled during build time. Please recompile or "
        "use native format");
#endif
}

void FieldNorms::write_netcdf()
{
#ifdef AMR_WIND_USE_NETCDF
    std::vector<double> buf(m_total_particles * m_var_names.size(), 0.0);
    m_scontainer->populate_buffer(buf);

    if (!amrex::ParallelDescriptor::IOProcessor()) return;
    auto ncf = ncutils::NCFile::open(m_ncfile_name, NC_WRITE);
    const std::string nt_name = "num_time_steps";
    // Index of the next timestep
    const size_t nt = ncf.dim(nt_name).len();
    {
        auto time = m_sim.time().new_time();
        ncf.var("time").put(&time, {nt}, {1});
    }

    std::vector<size_t> start{nt, 0};
    std::vector<size_t> count{1, 0};

    const int nvars = m_var_names.size();
    //    for (int iv = 0; iv < nvars; ++iv) {
    //        start[1] = 0;
    //        count[1] = 0;
    //        int offset = iv * m_scontainer->num_sampling_particles();
    //        for (const auto& obj : m_samplers) {
    //            count[1] = obj->num_points();
    //            auto grp = ncf.group(obj->label());
    //            auto var = grp.var(m_var_names[iv]);
    //            var.put(&buf[offset], start, count);
    //            offset += count[1];
    //        }
    //    }
    ncf.close();
#endif
}

} // namespace field_norms
} // namespace amr_wind
