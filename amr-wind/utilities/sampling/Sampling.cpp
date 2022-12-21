#include <memory>
#include <utility>

#include "amr-wind/utilities/sampling/Sampling.H"
#include "amr-wind/utilities/io_utils.H"
#include "amr-wind/utilities/ncutils/nc_interface.H"

#include "AMReX_ParmParse.H"

namespace amr_wind::sampling {

Sampling::Sampling(CFDSim& sim, std::string label)
    : m_sim(sim), m_label(std::move(label))
{}

Sampling::~Sampling() = default;

void Sampling::initialize()
{
    BL_PROFILE("amr-wind::Sampling::initialize");

    // Labels for the different sampler types
    amrex::Vector<std::string> labels;
    // Fields to be sampled - requested by user
    amrex::Vector<std::string> field_names;

    {
        amrex::ParmParse pp(m_label);
        pp.getarr("labels", labels);
        pp.getarr("fields", field_names);
        pp.query("output_frequency", m_out_freq);
        pp.query("output_format", m_out_fmt);
    }

    // Process field information
    m_ncomp = 0;
    auto& repo = m_sim.repo();
    for (const auto& fname : field_names) {
        if (!repo.field_exists(fname)) {
            amrex::Print()
                << "WARNING: Sampling: Non-existent field requested: " << fname
                << std::endl;
            continue;
        }

        auto& fld = repo.get_field(fname);
        m_ncomp += fld.num_comp();
        m_fields.emplace_back(&fld);
        ioutils::add_var_names(m_var_names, fld.name(), fld.num_comp());
    }

    // Load different probe types, default probe type is line
    int idx = 0;
    m_total_particles = 0;
    for (auto& lbl : labels) {
        const std::string key = m_label + "." + lbl;
        amrex::ParmParse pp1(key);
        std::string stype = "LineSampler";

        pp1.query("type", stype);
        auto obj = SamplerBase::create(stype, m_sim);
        obj->label() = lbl;
        obj->id() = idx++;
        obj->initialize(key);

        m_total_particles += obj->num_points();
        m_samplers.emplace_back(std::move(obj));
    }

    update_container();

    if (m_out_fmt == "netcdf") {
        prepare_netcdf_file();
    }
}

void Sampling::update_container()
{
    BL_PROFILE("amr-wind::Sampling::update_container");

    // Initialize the particle container based on user inputs
    m_scontainer = std::make_unique<SamplingContainer>(m_sim.mesh());
    m_scontainer->setup_container(m_ncomp);
    m_scontainer->initialize_particles(m_samplers);
    // Redistribute particles to appropriate boxes/MPI ranks
    m_scontainer->Redistribute();
    m_scontainer->num_sampling_particles() =
        static_cast<int>(m_total_particles);
}

void Sampling::update_sampling_locations()
{
    BL_PROFILE("amr-wind::Sampling::update_sampling_locations");

    for (const auto& obj : m_samplers) {
        obj->update_sampling_locations();
    }

    update_container();
}

void Sampling::post_advance_work()
{

    BL_PROFILE("amr-wind::Sampling::post_advance_work");
    const auto& time = m_sim.time();
    const int tidx = time.time_index();
    // Skip processing if it is not an output timestep
    if (!(tidx % m_out_freq == 0)) {
        return;
    }

    update_sampling_locations();

    m_scontainer->interpolate_fields(m_fields);

    process_output();
}

void Sampling::post_regrid_actions()
{

    BL_PROFILE("amr-wind::Sampling::post_regrid_actions");
    m_scontainer->Redistribute();
}

void Sampling::process_output()
{
    if (m_out_fmt == "native") {
        impl_write_native();
    } else if (m_out_fmt == "ascii") {
        write_ascii();
    } else if (m_out_fmt == "netcdf") {
        write_netcdf();
    } else {
        amrex::Abort("Sampling: Invalid output format encountered");
    }
}

void Sampling::impl_write_native()
{
    BL_PROFILE("amr-wind::Sampling::write_native");

    const std::string post_dir = "post_processing";
    const std::string name =
        amrex::Concatenate(m_label, m_sim.time().time_index());
    amrex::Vector<std::string> int_var_names{"uid", "set_id", "probe_id"};

    m_scontainer->WritePlotFile(
        post_dir, name, m_var_names, int_var_names,
        [=] AMREX_GPU_HOST_DEVICE(
            const SamplingContainer::SuperParticleType& p) {
            return p.id() > 0;
        });
}

void Sampling::write_ascii()
{
    BL_PROFILE("amr-wind::Sampling::write_ascii");
    amrex::Print() << "WARNING: Sampling: ASCII output will impact performance"
                   << std::endl;

    const std::string post_dir = "post_processing";
    const std::string sname =
        amrex::Concatenate(m_label, m_sim.time().time_index());

    if (!amrex::UtilCreateDirectory(post_dir, 0755)) {
        amrex::CreateDirectoryFailed(post_dir);
    }
    const std::string fname = post_dir + "/" + sname + ".txt";
    m_scontainer->WriteAsciiFile(fname);
}

void Sampling::prepare_netcdf_file()
{
#ifdef AMR_WIND_USE_NETCDF

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
    for (const auto& obj : m_samplers) {
        auto grp = ncf.def_group(obj->label());

        grp.def_dim(npart_name, obj->num_points());
        obj->define_netcdf_metadata(grp);
        grp.def_var("coordinates", NC_DOUBLE, {npart_name, "ndim"});
        for (const auto& vname : m_var_names)
            grp.def_var(vname, NC_DOUBLE, two_dim);
    }
    ncf.exit_def_mode();

    {
        const std::vector<size_t> start{0, 0};
        std::vector<size_t> count{0, AMREX_SPACEDIM};
        SamplerBase::SampleLocType locs;
        for (const auto& obj : m_samplers) {
            auto grp = ncf.group(obj->label());
            obj->populate_netcdf_metadata(grp);
            obj->sampling_locations(locs);
            auto xyz = grp.var("coordinates");
            count[0] = obj->num_points();
            xyz.put(&locs[0][0], start, count);
        }
    }

#else
    amrex::Abort(
        "NetCDF support was not enabled during build time. Please recompile or "
        "use native format");
#endif
}

void Sampling::write_netcdf()
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

    for (const auto& obj : m_samplers) {
        auto grp = ncf.group(obj->label());
        obj->output_netcdf_data(grp, nt);
    }

    std::vector<size_t> start{nt, 0};
    std::vector<size_t> count{1, 0};

    const int nvars = m_var_names.size();
    for (int iv = 0; iv < nvars; ++iv) {
        start[1] = 0;
        count[1] = 0;
        int offset = iv * m_scontainer->num_sampling_particles();
        for (const auto& obj : m_samplers) {
            auto grp = ncf.group(obj->label());
            auto var = grp.var(m_var_names[iv]);
            // Do sampler specific output if needed
            bool do_output = obj->output_netcdf_field(&buf[offset], var);
            // Do generic output if specific output returns true
            if (do_output) {
                count[1] = obj->num_points();
                var.put(&buf[offset], start, count);
                offset += count[1];
            }
        }
    }
    ncf.close();
#endif
}

} // namespace amr_wind::sampling
