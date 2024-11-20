#include "amr-wind/utilities/tagging/BoundFieldRefinement.H"
#include "amr-wind/CFDSim.H"

#include "AMReX.H"
#include "AMReX_ParmParse.H"

namespace amr_wind {

BoundFieldRefinement::BoundFieldRefinement(const CFDSim& sim)
    : m_sim(sim)
    , m_mesh(sim.mesh())
    , m_field_error(
          m_sim.mesh().maxLevel() + 1, std::numeric_limits<amrex::Real>::max())
    , m_grad_error(
          m_sim.mesh().maxLevel() + 1, std::numeric_limits<amrex::Real>::max())
{}

void BoundFieldRefinement::initialize(const std::string& key)
{
    amrex::ParmParse pp(key);
    std::string fname;
    pp.query("field_name", fname);
    const auto& repo = m_sim.repo();
    if (repo.field_exists(fname)) {
        m_field = &(m_sim.repo().get_field(fname));
        AMREX_ALWAYS_ASSERT(m_field->num_grow() > amrex::IntVect{0});
    } else if (repo.int_field_exists(fname)) {
        m_int_field = &(m_sim.repo().get_int_field(fname));
        AMREX_ALWAYS_ASSERT(m_int_field->num_grow() > amrex::IntVect{0});
    } else {
        amrex::Abort("BoundFieldRefinement: Cannot find field = " + fname);
    }

    amrex::Vector<amrex::Real> field_err;
    amrex::Vector<amrex::Real> grad_err;
    pp.queryarr("field_error", field_err);
    pp.queryarr("grad_error", grad_err);

    if ((field_err.empty()) && (grad_err.empty())) {
        amrex::Abort(
            "BoundFieldRefinement: Must specify at least one of field_error or "
            "grad_error");
    }

    {
        const int fcount = std::min(
            static_cast<int>(field_err.size()),
            static_cast<int>(m_field_error.size()));
        for (int i = 0; i < fcount; ++i) {
            m_field_error[i] = field_err[i];
        }
        m_max_lev_field = fcount - 1;
    }
    {
        const int fcount = std::min(
            static_cast<int>(grad_err.size()),
            static_cast<int>(m_grad_error.size()));
        for (int i = 0; i < fcount; ++i) {
            m_grad_error[i] = grad_err[i];
        }
        m_max_lev_grad = fcount - 1;
    }
    pp.query("bounds_file", m_boundsfile);
    if (!m_boundsfile.empty()) {
        std::ifstream boundsfile(m_boundsfile, std::ios::in);
        if (!boundsfile.good()) {
            amrex::Abort("Cannot find the bounds file " + m_boundsfile);
        }
        amrex::Real value1, value2, value3, value4, value5;
        while (boundsfile >> value1 >> value2 >> value3 >> value4 >> value5) {
            m_xmin_list.push_back(value1);
            m_xmax_list.push_back(value2);
            m_ymin_list.push_back(value3);
            m_ymax_list.push_back(value4);
            m_level_list.push_back(value5);
        }
    }
}

void BoundFieldRefinement::operator()(
    const int level,
    amrex::TagBoxArray& tags,
    const amrex::Real time,
    const int /*ngrow*/)
{
    const bool tag_grad = level <= m_max_lev_grad;
    if (tag_grad) {
        if (m_field != nullptr) {
            m_field->fillpatch(level, time, (*m_field)(level), 1);
        } else if (m_int_field != nullptr) {
            (*m_int_field)(level).FillBoundary(
                m_sim.repo().mesh().Geom(level).periodicity());
        }
    }

    if (m_field != nullptr) {
        const auto& mfab = (*m_field)(level);
        tag(level, tags, mfab);
    } else if (m_int_field != nullptr) {
        const auto& mfab = (*m_int_field)(level);
        tag(level, tags, mfab);
    }
}

} // namespace amr_wind
