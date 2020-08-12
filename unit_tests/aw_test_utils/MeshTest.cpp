#include "MeshTest.H"
#include "pp_utils.H"

namespace amr_wind_tests {

void MeshTest::create_mesh_instance()
{
    if (!m_mesh) {
        reset_prob_domain();
        m_mesh.reset(new AmrTestMesh());
    }
}

void MeshTest::initialize_mesh()
{
    if (m_need_params) populate_parameters();
    create_mesh_instance();

    m_mesh->initialize_mesh(0.0);
}

void MeshTest::reset_prob_domain()
{
    amrex::Vector<amrex::Real> problo(3), probhi(3);
    amrex::Vector<int> periodic(3);

    amrex::ParmParse pp("geometry");
    pp.getarr("prob_lo",problo,0,AMREX_SPACEDIM);
    pp.getarr("prob_hi",probhi,0,AMREX_SPACEDIM);
    pp.getarr("is_periodic",periodic,0,AMREX_SPACEDIM);

    amrex::RealBox rb(problo.data(), probhi.data());
    amrex::Geometry* gg = amrex::AMReX::top()->getDefaultGeometry();

    if (gg != nullptr) {
        gg->ResetDefaultProbDomain(rb);
        gg->ResetDefaultPeriodicity({periodic[0],periodic[1],periodic[2]});
    }
}

void MeshTest::populate_parameters()
{
    pp_utils::default_time_inputs();
    pp_utils::default_mesh_inputs();
    m_need_params = false;
}

}
