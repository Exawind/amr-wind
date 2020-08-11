#ifndef MMS_H
#define MMS_H

#include "amr-wind/core/Physics.H"
#include "amr-wind/core/Field.H"
#include "amr-wind/CFDSim.H"

namespace amr_wind {
namespace mms {

using FuncDef = double (*)(double, double, double);

/** Method of manufactured solutions physics
 *  \ingroup physics
 *
 */
class MMS : public Physics::Register<MMS>
{
public:
    static const std::string identifier() { return "MMS"; }

    MMS(const CFDSim& sim);

    virtual ~MMS() = default;

    void initialize_fields(int level, const amrex::Geometry& geom) override;

    void fill_src();

    amrex::Real compute_error(const int comp, const Field& field, FuncDef f);

    void post_init_actions() override;

    void post_regrid_actions() override;

    void pre_advance_work() override {}

    void post_advance_work() override;

private:
    const amr_wind::SimTime& m_time;
    const FieldRepo& m_repo;
    const amrex::AmrCore& m_mesh;
    Field& m_velocity;
    Field& m_density;
    Field& m_mms_vel_source;

    const int m_w = 18;
    const std::string m_output_fname = "mms.log";

    std::unordered_map<std::string, amrex::Real> m_params_map = {
        {"a", 1.05},   {"b", 2.15},    {"c", -3.2},    {"d", 10.1},
        {"beta", 2.2}, {"gamma", 2.4}, {"delta", 2.0}, {"kx", 1.0},
        {"ky", 1.0},   {"kz", 1.0},    {"nu", 1.0}};
};
} // namespace mms
} // namespace amr_wind

#endif /* MMS_H */
