#ifndef MESHTEST_H
#define MESHTEST_H

#include <memory>

#include "aw_test_utils/AmrexTest.H"
#include "aw_test_utils/AmrTestMesh.H"
#include "SimTime.H"

namespace amr_wind_tests {

/** Base class for test fixtures that require a valid AmrMesh instance
 */
class MeshTest : public AmrexTest
{
public:
    MeshTest() = default;

    virtual ~MeshTest() = default;

    /** Populate the ParmParse instance with all the necessary input variables
     *  for this test fixture
     */
    virtual void populate_parameters();

    //! Create the mesh instance
    virtual void create_mesh_instance();

    template<typename T>
    void create_mesh_instance()
    {
        if (!m_mesh) {
            reset_prob_domain();
            m_mesh.reset(new T());
        }
    }

    //! Perform initialization actions
    virtual void initialize_mesh();

    //! Reset the problem domain for this test
    virtual void reset_prob_domain();

    amr_wind::CFDSim& sim()
    {
        return mesh().sim();
    }

    amr_wind::SimTime& time()
    {
        return sim().time();
    }

    AmrTestMesh& mesh() { return *m_mesh; }

    template <typename T>
    T* mesh()
    {
        return dynamic_cast<T*>(m_mesh.get());
    }

protected:
    std::unique_ptr<AmrTestMesh> m_mesh;

    bool m_need_params{true};
};

}

#endif /* MESHTEST_H */
