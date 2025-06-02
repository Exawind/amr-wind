#include "amr-wind/utilities/sampling/FieldNorms.H"
#include "amr-wind/utilities/io_utils.H"
#include "amr-wind/utilities/ncutils/nc_interface.H"
#include <AMReX_MultiFabUtil.H>
#include <utility>
#include "AMReX_ParmParse.H"
#include "amr-wind/utilities/IOManager.H"

namespace amr_wind::field_norms {

FieldNorms::FieldNorms(CFDSim& sim, std::string label)
    : m_sim(sim), m_label(std::move(label))
{}

FieldNorms::~FieldNorms() = default;

void FieldNorms::initialize()
{
    BL_PROFILE("amr-wind::FieldNorms::initialize");

    {
        amrex::ParmParse pp(m_label);
        populate_output_parameters(pp);
        pp.query("mask_redundant_grids", m_use_mask);
    }

    const auto& io_mng = m_sim.io_manager();
    for (const auto& fld : io_mng.plot_fields()) {
        ioutils::add_var_names(m_var_names, fld->name(), fld->num_comp());
    }

    m_fnorms.resize(m_var_names.size(), 0.0);

    prepare_ascii_file();
}

amrex::Real FieldNorms::l2_norm(
    const amr_wind::Field& field, const int comp, const bool use_mask)
{
    amrex::Real nrm = 0.0;

    AMREX_ASSERT(comp >= 0 && comp < field.num_comp());

    const int finest_level = field.repo().num_active_levels() - 1;
    const auto& geom = field.repo().mesh().Geom();
    constexpr int nghost = 0;

    const auto& mesh = field.repo().mesh();

    for (int lev = 0; lev <= finest_level; lev++) {

        const amrex::Real cell_vol = geom[lev].CellSize()[0] *
                                     geom[lev].CellSize()[1] *
                                     geom[lev].CellSize()[2];

        amrex::iMultiFab level_mask;
        if (use_mask) {
            if (lev < finest_level) {
                level_mask = makeFineMask(
                    mesh.boxArray(lev), mesh.DistributionMap(lev),
                    mesh.boxArray(lev + 1), mesh.refRatio(lev), 1, 0);
            } else {
                level_mask.define(
                    mesh.boxArray(lev), mesh.DistributionMap(lev), 1, 0,
                    amrex::MFInfo());
                level_mask.setVal(1);
            }
        } else {
            // Always on
            level_mask.define(
                mesh.boxArray(lev), mesh.DistributionMap(lev), 1, 0,
                amrex::MFInfo());
            level_mask.setVal(1);
        }

        nrm += amrex::ReduceSum(
            field(lev), level_mask, nghost,
            [=] AMREX_GPU_HOST_DEVICE(
                amrex::Box const& bx,
                amrex::Array4<amrex::Real const> const& field_arr,
                amrex::Array4<int const> const& mask_arr) -> amrex::Real {
                amrex::Real nrm_fab = 0.0;

                amrex::Loop(bx, [=, &nrm_fab](int i, int j, int k) noexcept {
                    nrm_fab += cell_vol * field_arr(i, j, k, comp) *
                               field_arr(i, j, k, comp) * mask_arr(i, j, k);
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
    for (const auto& fld : m_sim.io_manager().plot_fields()) {
        for (int comp = 0; comp < fld->num_comp(); ++comp) {
            m_fnorms[ind++] = l2_norm(*fld, comp, m_use_mask);
        }
    }
}

void FieldNorms::output_actions()
{
    BL_PROFILE("amr-wind::FieldNorms::output_actions");

    process_field_norms();
    write_ascii();
}

void FieldNorms::prepare_ascii_file()
{
    BL_PROFILE("amr-wind::FieldNorms::prepare_ascii_file");

    const std::string post_dir = m_sim.io_manager().post_processing_directory();
    const std::string sname =
        amrex::Concatenate(m_label, m_sim.time().time_index());

    m_out_fname = post_dir + "/" + sname + ".txt";

    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ofstream f(m_out_fname.c_str());
        f << "time_step time";
        for (const auto& m_var_name : m_var_names) {
            f << ' ' << m_var_name;
        }
        f << std::endl;
        f.close();
    }
}

void FieldNorms::write_ascii()
{
    BL_PROFILE("amr-wind::FieldNorms::write_ascii");

    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ofstream f(m_out_fname.c_str(), std::ios_base::app);
        f << m_sim.time().time_index() << std::fixed
          << std::setprecision(m_precision) << std::setw(m_width)
          << m_sim.time().new_time();
        for (double m_fnorm : m_fnorms) {
            f << std::setw(m_width) << m_fnorm;
        }
        f << std::endl;
        f.close();
    }
}

} // namespace amr_wind::field_norms
