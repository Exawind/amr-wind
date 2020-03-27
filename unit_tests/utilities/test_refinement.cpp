#include <sstream>

#include "aw_test_utils/AmrexTest.H"
#include "aw_test_utils/MeshTest.H"

#include "AMReX_Box.H"
#include "AMReX_BoxArray.H"
#include "AMReX_BoxList.H"
#include "AMReX_Geometry.H"
#include "AMReX_RealBox.H"
#include "AMReX_Vector.H"

#include "CartBoxRefinement.H"

namespace amr_wind_tests {

//! Custom mesh class to provide error estimator based on refinement criteria
class NestRefineMesh : public AmrTestMesh
{
public:
    amrex::Vector<std::unique_ptr<amr_wind::RefinementCriteria>>&
    refine_criteria_vec()
    {
        return m_refine_crit;
    }

protected:
    virtual void
    ErrorEst(int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow) override
    {
        for (auto& ref: m_refine_crit)
            (*ref)(lev, tags, time, ngrow);
    }

private:
    amrex::Vector<std::unique_ptr<amr_wind::RefinementCriteria>> m_refine_crit;
};

//! Custom test fixture for Cartesian Box refinement
class NestRefineTest : public MeshTest
{};

TEST_F(NestRefineTest, box_refine)
{
    populate_parameters();

    {
        amrex::ParmParse pp("amr");
        amrex::Vector<int> ncell{{16, 128, 16}};

        pp.add("max_level", 1);
        pp.addarr("n_cell", ncell);
    }
    {
        amrex::ParmParse pp("geometry");
        amrex::Vector<amrex::Real> problo {{-20.0, -100.0, 0.0}};
        amrex::Vector<amrex::Real> probhi {{20.0, 100.0, 30.0}};

        pp.addarr("prob_lo", problo);
        pp.addarr("prob_hi", probhi);
    }

    // Create the "input file"
    std::stringstream ss;
    ss << "1 // Number of levels" << std::endl;
    ss << "4 // Number of boxes at this level" << std::endl;
    ss << "-10.0 -75.0 0.0 15.0 -65.0 20.0" << std::endl;
    ss << "-10.0 -55.0 0.0 15.0 -45.0 20.0" << std::endl;
    ss << "-10.0  25.0 0.0 15.0  35.0 20.0" << std::endl;
    ss << "-10.0  65.0 0.0 15.0  75.0 20.0" << std::endl;

    create_mesh_instance<NestRefineMesh>();
    std::unique_ptr<amr_wind::CartBoxRefinement> box_refine(new amr_wind::CartBoxRefinement);
    box_refine->read_inputs(mesh(), ss);

    // Store the target boxarray for future tests
    auto targets = box_refine->boxarray_vec();

    mesh<NestRefineMesh>()->refine_criteria_vec().push_back(std::move(box_refine));
    initialize_mesh();

    auto ba1 = mesh().boxArray(1);
    // NOTE: Box definitions were based on Level 0 to tag cells
    auto ba2 = targets[0].refine(2);
    EXPECT_TRUE(ba1.contains(ba2));
}

}  // amr_wind_tests
