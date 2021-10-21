#include "amr-wind/utilities/sampling/IsoSampling.H"
#include "amr-wind/utilities/io_utils.H"
#include "amr-wind/utilities/ncutils/nc_interface.H"

#include "AMReX_ParmParse.H"

namespace amr_wind {
namespace sampling {

IsoSampling::IsoSampling(CFDSim& sim, const std::string& label)
    : m_sim(sim), m_label(label)
{}

IsoSampling::~IsoSampling() = default;

void IsoSampling::initialize()
{
    BL_PROFILE("amr-wind::IsoSampling::initialize");

    // Labels for the different sampler types
    amrex::Vector<std::string> labels;

    {
        amrex::ParmParse pp(m_label);
        pp.getarr("labels", labels);
        pp.query("output_frequency", m_out_freq);
        pp.query("output_format", m_out_fmt);
    }

    // Load different probe types, default probe type is line
    int idx = 0;
    m_total_particles = 0;
    int ncomp = 0;
    auto& repo = m_sim.repo();
    for (auto& lbl : labels) {
        const std::string key = m_label + "." + lbl;
        amrex::ParmParse pp1(key);
        std::string stype = "IsoProbeSampler";

        pp1.query("type", stype);
        auto obj = SamplerBase::create(stype, m_sim);
        obj->label() = lbl;
        obj->id() = idx++;
        obj->initialize(key);

        m_total_particles += obj->num_points();
        m_samplers.emplace_back(std::move(obj));

        // Process field information (one per sampler)
        std::string fname;
        pp1.query("field", fname);
        if (!repo.field_exists(fname)) {
            amrex::Print()
                << "WARNING: IsoSampling: Non-existent field requested: "
                << fname << std::endl;
            continue;
        }
        auto& fld = repo.get_field(fname);
        if (fld.num_comp()!=1) {
            amrex::Abort("IsoSampling: Non-scalar field requested: " + fname);
        }
        ncomp = fld.num_comp(); // Only one data field per particle
        m_fields.emplace_back(&fld);
        ioutils::add_var_names(m_var_names, fld.name(), fld.num_comp());
        SamplerBase::SampleValType fval;
        pp1.query("field_value", fval);
        m_field_values.emplace_back(fval);
    }

    // Initialize the particle container based on user inputs
    m_scontainer.reset(new SamplingContainer(m_sim.mesh()));
    // Real components
    // = current value     (index = 0)
    // + target field value        (1)
    // + current left value        (2)
    // + current right value       (3)
    // + current left location     (4           :3+  spacedim)
    // + current right location    (4+  spacedim:3+2*spacedim)
    // + sample init location      (4+2*spacedim:3+3*spacedim)
    // + sample orientation        (4+3*spacedim:3+4*spacedim)
    // Integer components
    // = flag used in searching algorithm (index=0)
    ncomp += 3 + 4 * AMREX_SPACEDIM;
    const int nicomp = 1;
    // Store number of components
    m_preals = ncomp;
    m_pints = nicomp;
    m_scontainer->setup_container(m_preals,m_pints);
    m_scontainer->initialize_particles(m_samplers,m_field_values);
    // Populate particle real component names
    m_pcomp_names.emplace_back("field"); // will be replaced
    m_pcomp_names.emplace_back("target");
    m_pcomp_names.emplace_back("lval");
    m_pcomp_names.emplace_back("rval");
    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
        m_pcomp_names.emplace_back("lpos"+std::to_string(n+1));
    }
    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
        m_pcomp_names.emplace_back("rpos"+std::to_string(n+1));
    }
    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
        m_pcomp_names.emplace_back("pos0"+std::to_string(n+1));
    }
    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
        m_pcomp_names.emplace_back("ori"+std::to_string(n+1));
    }
    // Redistribute particles to appropriate boxes/MPI ranks
    m_scontainer->Redistribute();
    m_scontainer->num_sampling_particles() = m_total_particles;
    // Set up initial bounds for bisection
    m_scontainer->iso_initbounds(m_fields);

    if (m_out_fmt == "netcdf") prepare_netcdf_file();
}

void IsoSampling::post_advance_work()
{
    BL_PROFILE("amr-wind::IsoSampling::post_advance_work");
    const auto& time = m_sim.time();
    const int tidx = time.time_index();
    // Skip processing if it is not an output timestep
    if (!(tidx % m_out_freq == 0)) return;

    m_scontainer->iso_relocate(m_fields);

    process_output();
}

void IsoSampling::post_regrid_actions()
{

    BL_PROFILE("amr-wind::IsoSampling::post_regrid_actions");
    m_scontainer->Redistribute();
}

void IsoSampling::process_output()
{
    if (m_out_fmt == "native") {
        impl_write_native();
    } else if (m_out_fmt == "ascii") {
        write_ascii();
    } else if (m_out_fmt == "netcdf") {
        write_netcdf();
    } else {
        amrex::Abort("IsoSampling: Invalid output format encountered");
    }
}

void IsoSampling::impl_write_native()
{
    BL_PROFILE("amr-wind::IsoSampling::write_native");

    const std::string post_dir = "post_processing";
    const std::string sdir =
        amrex::Concatenate(m_label, m_sim.time().time_index());
    amrex::Vector<std::string> int_var_names{"uid", "set_id", "probe_id"};

    m_scontainer->WritePlotFile(
        post_dir, sdir, m_var_names, int_var_names,
        [=] AMREX_GPU_HOST_DEVICE(
            const SamplingContainer::SuperParticleType& p) {
            return p.id() > 0;
        });
}

void IsoSampling::write_ascii()
{
    BL_PROFILE("amr-wind::IsoSampling::write_ascii");
    amrex::Print()
        << "WARNING: IsoSampling: ASCII output will impact performance"
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

void IsoSampling::prepare_netcdf_file()
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
    ncf.put_attr("title", "AMR-Wind data isosampling output");
    ncf.put_attr("version", ioutils::amr_wind_version());
    ncf.put_attr("created_on", ioutils::timestamp());
    ncf.def_dim(nt_name, NC_UNLIMITED);
    ncf.def_dim("ndim", AMREX_SPACEDIM);
    ncf.def_var("time", NC_DOUBLE, {nt_name});
    // Define groups for each sampler
    int is = 0;
    for (const auto& obj : m_samplers) {
        auto grp = ncf.def_group(obj->label());

        grp.def_dim(npart_name, obj->num_points());
        obj->define_netcdf_metadata(grp);
        grp.def_var("coordinates", NC_DOUBLE, {npart_name, "ndim"});
        int iv = 0;
        for (const auto& vname : m_pcomp_names) {
            if (iv == 0) {
                grp.def_var(m_var_names[is],NC_DOUBLE, two_dim);
            } else {
                grp.def_var(vname, NC_DOUBLE, two_dim);
            }
            ++iv;
        }
        ++is;
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

void IsoSampling::write_netcdf()
{
#ifdef AMR_WIND_USE_NETCDF
    std::vector<double> buf(m_total_particles * m_preals, 0.0);
    m_scontainer->populate_buffer(buf);
    std::vector<int> ibuf(m_total_particles * m_pints, 0);
    m_scontainer->populate_buffer(ibuf);

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

    const int nvars = m_pcomp_names.size();
    for (int iv = 0; iv < nvars; ++iv) {
        start[1] = 0;
        count[1] = 0;
        int offset = iv * m_scontainer->num_sampling_particles();
        int is = 0;
        for (const auto& obj : m_samplers) {
            count[1] = obj->num_points();
            auto grp = ncf.group(obj->label());
            if (iv == 0) {
                // Use variable name corresponding to current sampler
                auto var = grp.var(m_var_names[is]);
            } else {
                // Use particle component name
                auto var = grp.var(m_pcomp_names[iv]);
            }
            var.put(&buf[offset], start, count);
            offset += count[1];
            ++is;
        }
    }
    ncf.close();
#endif
}

} // namespace sampling
} // namespace amr_wind
