#include "aw_test_utils/MeshTest.H"
#include "test_act_utils.H"

#include "amr-wind/wind_energy/actuator/Actuator.H"
#include "amr-wind/wind_energy/actuator/ActuatorContainer.H"
#include "amr-wind/wind_energy/actuator/ActuatorModel.H"
#include "amr-wind/wind_energy/actuator/ActParser.H"
#include "amr-wind/wind_energy/actuator/actuator_types.H"
#include "amr-wind/wind_energy/actuator/wing/ActuatorWing.H"
#include "amr-wind/wind_energy/actuator/wing/wing_ops.H"
#include "amr-wind/core/gpu_utils.H"
#include "amr-wind/core/vs/vector_space.H"
#include "amr-wind/utilities/trig_ops.H"

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

struct FlatPlate : public act::ActuatorType
{
    using InfoType = act::ActInfo;
    using GridType = act::ActGrid;
    using MetaType = act::WingBaseData;
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
    ::amr_wind_tests::FlatPlate::DataType& data,
    const amr_wind::actuator::utils::ActParser& pp)
{
    amr_wind::actuator::wing::read_inputs(data.m_meta, data.m_info, pp);
}

template <>
void init_data_structures<::amr_wind_tests::FlatPlate>(
    ::amr_wind_tests::FlatPlate::DataType& data)
{
    auto& grid = data.m_grid;
    auto& wdata = data.m_meta;

    amr_wind::actuator::wing::init_data_structures(wdata, grid);

    {
        // Check that the wing coordinates were initialized correctly
        auto wing_len = vs::mag(grid.pos.back() - grid.pos.front());
        EXPECT_NEAR(wing_len, 8.0, 1.0e-12);
        EXPECT_NEAR(
            grid.orientation[0].x().unit() & vs::Vector::ihat(), 1.0, 1.0e-12);
        EXPECT_NEAR(vs::mag_sqr(grid.orientation[0]), 3.0, 1.0e-12);
    }
}

template <>
inline void update_positions<::amr_wind_tests::FlatPlate>(
    ::amr_wind_tests::FlatPlate::DataType&)
{}

template <>
inline void update_velocities<::amr_wind_tests::FlatPlate>(
    ::amr_wind_tests::FlatPlate::DataType& data)
{
    const auto& grid = data.m_grid;
    const auto& pos = grid.vel_pos;
    const auto& vel = grid.vel;
    amrex::Real rerr = 0.0;
    for (int i = 0; i < grid.vel_pos.size(); ++i) {
        const amrex::Real val = pos[i].x() + pos[i].y() + pos[i].z();
        const vs::Vector vgold{val, val, val};
        rerr += vs::mag_sqr(vel[i] - vgold);
    }
    EXPECT_NEAR(rerr, 0.0, 1.0e-12);
}

template <>
inline void compute_forces<::amr_wind_tests::FlatPlate>(
    ::amr_wind_tests::FlatPlate::DataType& data)
{
    auto& grid = data.m_grid;
    const int npts = data.m_meta.num_pts;
    const amrex::Real pitch = ::amr_wind::utils::radians(data.m_meta.pitch);
    for (int ip = 0; ip < npts; ++ip) {
        const auto& tmat = grid.orientation[ip];
        // Effective velocity at the wing control point in local frame
        auto wvel = tmat & grid.vel[ip];
        // Set spanwise component to zero to get a pure 2D velocity
        wvel.y() = 0.0;

        const auto vmag = vs::mag(wvel);
        const auto phi = std::atan2(wvel.z(), wvel.x());
        const auto aoa = pitch - phi;

        // Make up some Cl, Cd values
        const auto cl = amr_wind::utils::two_pi() * aoa;
        const auto cd = cl * std::sin(aoa);

        // Assume unit chord
        const auto qval = 0.5 * vmag * vmag * grid.dx[ip];
        const auto lift = qval * cl;
        const auto drag = qval * cd;
        // Determine unit vector parallel and perpendicular to velocity vector
        const auto drag_dir = (wvel / vmag) & tmat;
        const auto lift_dir = drag_dir ^ tmat.y();

        // Compute force on fluid from this section of wing
        grid.force[ip] = -(lift_dir * lift + drag * drag_dir);
    }
}

} // namespace ops

template class ::amr_wind::actuator::
    ActModel<::amr_wind_tests::FlatPlate, ::amr_wind::actuator::ActSrcLine>;
} // namespace actuator
} // namespace amr_wind

namespace amr_wind_tests {

TEST_F(ActFlatPlateTest, act_model_init)
{
    initialize_mesh();

    {
        amrex::ParmParse pp("Actuator.FlatPlateLine");
        pp.add("num_points", 11);
        pp.addarr("start", amrex::Vector<amrex::Real>{16.0, 12.0, 16.0});
        pp.addarr("end", amrex::Vector<amrex::Real>{16.0, 20.0, 16.0});
        pp.addarr("normal", amrex::Vector<amrex::Real>{0.0, 0.0, 1.0});
        pp.addarr("epsilon", amrex::Vector<amrex::Real>{3.0, 3.0, 3.0});
        pp.add("pitch", 6.0);
    }

    ::amr_wind::actuator::ActModel<FlatPlate> flat_plate(sim(), "F1", 0);
    {
        amr_wind::actuator::utils::ActParser pp(
            "Actuator.FlatPlateLine", "Actuator.F1");
        flat_plate.read_inputs(pp);
    }

    amrex::Vector<int> act_proc_count(amrex::ParallelDescriptor::NProcs(), 0);
    flat_plate.setup_actuator_source(act_proc_count);

    const auto& info = flat_plate.info();
    EXPECT_EQ(info.root_proc, 0);
    EXPECT_EQ(info.procs.size(), amrex::ParallelDescriptor::NProcs());
}

TEST_F(ActFlatPlateTest, actuator_init)
{
    initialize_mesh();
    auto& vel = sim().repo().declare_field("velocity", 3, 3);
    init_field(vel);

    amrex::Vector<std::string> actuators{"T1", "T2"};
    {
        amrex::ParmParse pp("Actuator");
        pp.addarr("labels", actuators);
        pp.add("type", std::string("FlatPlateLine"));
    }
    {
        amrex::ParmParse pp("Actuator.FlatPlateLine");
        pp.add("num_points", 11);
        pp.addarr("normal", amrex::Vector<amrex::Real>{0.0, 0.0, 1.0});
        pp.addarr("epsilon", amrex::Vector<amrex::Real>{1.0, 1.0, 1.0});
        pp.add("pitch", 6.0);
    }
    {
        for (int i = 0; i < 2; ++i) {
            amrex::Real zloc = 8.0 + 16.0 * i;
            amrex::ParmParse pp("Actuator." + actuators[i]);
            pp.addarr("start", amrex::Vector<amrex::Real>{16.0, 12.0, zloc});
            pp.addarr("end", amrex::Vector<amrex::Real>{16.0, 20.0, zloc});
        }
    }

    ::amr_wind::actuator::Actuator act(sim());
    act.pre_init_actions();
    act.post_init_actions();
}

} // namespace amr_wind_tests
