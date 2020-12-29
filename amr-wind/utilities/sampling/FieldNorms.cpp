#include "amr-wind/utilities/sampling/FieldNorms.H"
#include "amr-wind/utilities/io_utils.H"
#include "amr-wind/utilities/ncutils/nc_interface.H"
#include <AMReX_MultiFabUtil.H>
#include "AMReX_ParmParse.H"
#include "amr-wind/utilities/IOManager.H"

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

    {
        amrex::ParmParse pp(m_label);
        pp.queryarr("fields", field_names);
        pp.query("output_frequency", m_out_freq);
    }

    auto & io_man = m_sim.io_manager();
    int ncomp = 0;
    for(const auto& fld : io_man.plot_fields()){
        ncomp += fld->num_comp();
        ioutils::add_var_names(m_var_names, fld->name(), fld->num_comp());
    }

    m_fnorms.resize(m_var_names.size(), 0.0);

    prepare_ascii_file();

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
    write_ascii();
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

} // namespace field_norms
} // namespace amr_wind
