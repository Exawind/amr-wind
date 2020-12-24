#include "aw_test_utils/MeshTest.H"

#include "amr-wind/wind_energy/actuator/ActuatorContainer.H"
#include "amr-wind/wind_energy/actuator/ActuatorModel.H"
#include "amr-wind/wind_energy/actuator/actuator_types.H"
#include "amr-wind/core/gpu_utils.H"
#include "amr-wind/core/vs/vector_space.H"

namespace amr_wind_tests {
namespace {
class ActFlatPlateTest : public MeshTest
{
protected:
    void populate_parameters() override
    {
        MeshTest::populate_parameters();

        {
            amrex::ParmParse pp("amr");
            amrex::Vector<int> ncell{{32, 32, 32}};
            pp.add("max_level", 0);
            pp.add("max_grid_size", 16);
            pp.addarr("n_cell", ncell);
        }
        {
            amrex::ParmParse pp("geometry");
            amrex::Vector<amrex::Real> problo{{0.0, 0.0, 0.0}};
            amrex::Vector<amrex::Real> probhi{{32.0, 32.0, 32.0}};

            pp.addarr("prob_lo", problo);
            pp.addarr("prob_hi", probhi);
        }
    }
};

namespace act = amr_wind::actuator;
namespace vs = amr_wind::vs;

struct FlatPlateData
{
    // Number of point along the wing
    int num_pts{11};
    // Pitch angle (i.e., angle of attack) of wing
    amrex::Real pitch{6.0};
    // User input for gaussian smearing parameter
    vs::Vector eps_inp{2.0, 2.0, 2.0};
    // Normal to the wing (used to compute chord direction)
    vs::Vector normal{0.0, 0.0, 1.0};

    vs::Vector begin{0.0, 0.0, 0.0};
    vs::Vector end{0.0, 0.0, 0.0};
};

struct FlatPlate : public act::ActuatorType
{
    using InfoType = act::ActInfo;
    using GridType = act::ActGrid;
    using MetaType = FlatPlateData;
    using DataType = act::ActDataHolder<FlatPlate>;

    static const std::string identifier() { return "FlatPlate"; }
};

} // namespace

} // namespace amr_wind_tests

namespace amr_wind {
namespace actuator {

namespace ops {

template <>
void read_inputs<::amr_wind_tests::FlatPlate>(
    ::amr_wind_tests::FlatPlate::DataType& data)
{
    auto& info = data.m_info;
    auto& fpobj = data.m_meta;

    fpobj.begin = {16.0, 12.0, 16.0};
    fpobj.end = {16.0, 20.0, 16.0};

    // Create a bounding box extending to max search radius in each direction
    info.bound_box = amrex::RealBox(8.0, 4.0, 8.0, 24.0, 28.0, 24.0);
}

template <>
void init_data_structures<::amr_wind_tests::FlatPlate>(
    ::amr_wind_tests::FlatPlate::DataType& data)
{
    auto& grid = data.m_grid;
    auto& obj = data.m_meta;

    int npts = obj.num_pts;
    grid.resize(npts);

    // Wing span
    auto wspan = obj.end - obj.begin;
    // Compute chord/flow direction as a cross-product
    auto chord = (wspan ^ obj.normal);
    // Set up global to local transformation matrix
    auto tmat = vs::Tensor(chord.unit(), wspan.unit(), obj.normal.unit());

    // Equal spacing along span
    auto dx = (1.0 / static_cast<amrex::Real>(npts - 1)) * wspan;

    for (int i = 0; i < npts; ++i) {
        grid.pos[i] = obj.begin + static_cast<amrex::Real>(i) * dx;
        grid.epsilon[i] = obj.eps_inp;
        grid.orientation[i] = tmat;
    }

    // Initialize remaining data
    grid.force.assign(npts, vs::Vector::zero());
    grid.vel_pos.assign(grid.pos.begin(), grid.pos.end());
    grid.vel.assign(npts, vs::Vector::zero());

    {
        // Check that the wing coordinates were initialized correctly
        auto wing_len = vs::mag(grid.pos.back() - grid.pos.front());
        EXPECT_NEAR(wing_len, 8.0, 1.0e-12);

        // Check that the orientation tensor was initialized correctly
        EXPECT_NEAR(chord.unit() & vs::Vector::ihat(), 1.0, 1.0e-12);
        EXPECT_NEAR(vs::mag_sqr(tmat), 3.0, 1.0e-12);
    }
}

} // namespace ops

template class ::amr_wind::actuator::
    ActModel<::amr_wind_tests::FlatPlate, ::amr_wind::actuator::ActSrcLine>;
} // namespace actuator
} // namespace amr_wind

namespace amr_wind_tests {

TEST_F(ActFlatPlateTest, test_init)
{
    initialize_mesh();

    ::amr_wind::actuator::ActModel<FlatPlate> flat_plate(
        sim(), "flat_plate1", 0);
    flat_plate.pre_init_actions();
    flat_plate.post_init_actions();

    const auto& info = flat_plate.info();
    EXPECT_EQ(info.root_proc, 0);
    EXPECT_EQ(info.procs.size(), amrex::ParallelDescriptor::NProcs());
}

} // namespace amr_wind_tests
