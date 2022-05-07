#include "aw_test_utils/MeshTest.H"
#include "test_act_utils.H"

#include "amr-wind/wind_energy/actuator/Actuator.H"
#include "amr-wind/wind_energy/actuator/ActuatorContainer.H"
#include "amr-wind/wind_energy/actuator/ActuatorModel.H"
#include "amr-wind/wind_energy/actuator/ActParser.H"
#include "amr-wind/wind_energy/actuator/disk/Joukowsky_ops.H"

namespace amr_wind_tests {
namespace {
class ActJoukowskyTest : public MeshTest
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

struct Joukowsky : public act::Joukowsky
{
    static std::string identifier() { return "TestJoukowsky"; }
};

} // namespace

} // namespace amr_wind_tests

namespace amr_wind {
namespace actuator {
namespace ops {

template <>
struct ReadInputsOp<::amr_wind_tests::Joukowsky, ActSrcDisk>
{
    void operator()(
        ::amr_wind_tests::Joukowsky::DataType& data, const utils::ActParser& pp)
    {
        ReadInputsOp<::amr_wind::actuator::Joukowsky, ActSrcDisk> actual_op;
        actual_op(data, pp);
    }
};

template <>
struct InitDataOp<::amr_wind_tests::Joukowsky, ActSrcDisk>
{
    void operator()(::amr_wind_tests::Joukowsky::DataType& /*data*/) {}
};

template <>
struct ComputeForceOp<::amr_wind_tests::Joukowsky, ActSrcDisk>
{
    void operator()(::amr_wind_tests::Joukowsky::DataType& /*data*/) {}
};

template <>
struct ProcessOutputsOp<::amr_wind_tests::Joukowsky, ActSrcDisk>
{
    ProcessOutputsOp<::amr_wind_tests::Joukowsky, ActSrcDisk>(
        ::amr_wind_tests::Joukowsky::DataType& /**/)
    {}
    void operator()(::amr_wind_tests::Joukowsky::DataType& /*data*/) {}
    void read_io_options(const utils::ActParser& /**/) {}
    void prepare_outputs(const std::string& /**/) {}
    void write_outputs(){};
};

} // namespace ops
template class ::amr_wind::actuator::
    ActModel<::amr_wind_tests::Joukowsky, ::amr_wind::actuator::ActSrcDisk>;
} // namespace actuator
} // namespace amr_wind

namespace amr_wind_tests {
class ActPhysicsTest : public ::amr_wind::actuator::Actuator
{
public:
    explicit ActPhysicsTest(::amr_wind::CFDSim& sim)
        : ::amr_wind::actuator::Actuator(sim)
    {}

protected:
    void prepare_outputs() override final {}
};

TEST_F(ActJoukowskyTest, act_model_init)
{
    initialize_mesh();
    sim().repo().declare_field("actuator_src_term", 3, 0);
    amr_wind::actuator::ActuatorContainer::ParticleType::NextID(1U);
    {
        amrex::ParmParse pp("Actuator.TestJoukowskyDisk");
        pp.add("num_points_t", 3);
        pp.add("num_points_r", 2);
    }
    {
        amrex::ParmParse pp("Actuator");
        pp.add("type", std::string("TestJoukowskyDisk"));
        pp.addarr("labels", {"D1"});
    }
    ::amr_wind::actuator::ActModel<Joukowsky, amr_wind::actuator::ActSrcDisk>
        model(sim(), "D1", 0);
    ActPhysicsTest act(sim());
    act.pre_init_actions();
}
} // namespace amr_wind_tests