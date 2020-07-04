#include "amr-wind/utilities/sampling/Sampling.H"
#include "amr-wind/utilities/io_utils.H"

#include "AMReX_ParmParse.H"

namespace amr_wind {
namespace sampling {

Sampling::Sampling(CFDSim& sim, const std::string& label)
    : m_sim(sim), m_label(label)
{}

Sampling::~Sampling() = default;

void Sampling::initialize()
{
    // Labels for the different sampler types
    amrex::Vector<std::string> labels;
    amrex::Vector<std::string> field_names;

    {
        amrex::ParmParse pp(m_label);
        pp.getarr("labels", labels);
        pp.getarr("fields", field_names);
        pp.query("output_frequency", m_out_freq);
        pp.query("output_format", m_out_fmt);
    }

    // Process field information
    int ncomp = 0;
    auto& repo = m_sim.repo();
    for (const auto& fname : field_names) {
        if (!repo.field_exists(fname)) {
            amrex::Print()
                << "WARNING: Sampling: Non-existent field requested: " << fname
                << std::endl;
            continue;
        }

        auto& fld = repo.get_field(fname);
        ncomp += fld.num_comp();
        m_fields.emplace_back(&fld);
        ioutils::add_var_names(m_var_names, fld.name(), fld.num_comp());
    }

    int idx = 0;
    for (auto& lbl : labels) {
        const std::string key = m_label + "/" + lbl;
        amrex::ParmParse pp1(key);
        std::string stype = "LineSampler";

        pp1.query("type", stype);
        auto obj = SamplerBase::create(stype, m_sim);
        obj->label() = lbl;
        obj->id() = idx++;
        obj->initialize(key);

        m_samplers.emplace_back(std::move(obj));
    }

    m_scontainer.reset(new SamplingContainer(m_sim.mesh()));
    m_scontainer->setup_container(ncomp);
    m_scontainer->initialize_particles(m_samplers);
    // Redistribute particles to appropriate boxes/MPI ranks
    m_scontainer->Redistribute();
}

void Sampling::post_advance_work()
{
    const auto& time = m_sim.time();
    const int tidx = time.time_index();
    // Skip processing if it is not an output timestep
    if (!(tidx % m_out_freq == 0)) return;

    m_scontainer->interpolate_fields(m_fields);

    process_output();
}

void Sampling::process_output()
{
    if (m_out_fmt == "native") {
        write_native();
    } else if (m_out_fmt == "ascii") {
        amrex::Print() << "WARNING: Sampling: ASCII format will impact performance"
                       << std::endl;
        m_scontainer->WriteAsciiFile(m_label);
    } else if (m_out_fmt == "netcdf") {
        amrex::Abort("Sampling: NetCDF format not yet supported");
    } else {
        amrex::Abort("Sampling: Invalid output format encountered");
    }
}

void Sampling::write_native()
{
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

} // namespace sampling
} // namespace amr_wind
