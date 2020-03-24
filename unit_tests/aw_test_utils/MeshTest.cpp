#include "MeshTest.H"
#include "pp_utils.H"

namespace amr_wind_tests {

void MeshTest::create_mesh_instance()
{
    if (!m_mesh) m_mesh.reset(new AmrTestMesh());
}

void MeshTest::initialize_mesh()
{
    if (m_need_params) populate_parameters();
    create_mesh_instance();

    m_mesh->initialize_mesh(0.0);
}

void MeshTest::populate_parameters()
{
    pp_utils::default_time_inputs();
    pp_utils::default_mesh_inputs();
    m_need_params = false;
}

}
