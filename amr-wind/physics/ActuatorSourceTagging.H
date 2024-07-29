#ifndef ACTUATORSOURCETAGGING_H
#define ACTUATORSOURCETAGGING_H

#include "amr-wind/core/Physics.H"
#include "amr-wind/core/Field.H"
#include "amr-wind/core/IntField.H"

namespace amr_wind {

/** Tracer Tagging physics
 *  \ingroup physics
 */
class ActuatorSourceTagging : public Physics::Register<ActuatorSourceTagging>
{
public:
    static const std::string identifier() { return "ActuatorSourceTagging"; }

    explicit ActuatorSourceTagging(CFDSim& sim);

    virtual ~ActuatorSourceTagging() = default;

    void initialize_fields(int level, const amrex::Geometry& geom) override;

    void post_init_actions() override;

    void post_regrid_actions() override {}

    void pre_advance_work() override {}

    void post_advance_work() override;

private:
    CFDSim& m_sim;
    Field* m_act_src{nullptr};
    IntField* m_iblank{nullptr};
    Field* m_tracer{nullptr};
    amrex::Real m_src_threshold{0.1};
    bool m_tag_hole{false};
    bool m_tag_fringe{true};
};

} // namespace amr_wind

#endif /* ACTUATORSOURCETAGGING_H */
