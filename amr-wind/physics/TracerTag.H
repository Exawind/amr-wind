#ifndef TracerTag_H
#define TracerTag_H

#include "amr-wind/core/Physics.H"
#include "amr-wind/core/Field.H"

namespace amr_wind {

/** Tracer Tagging physics
 *  \ingroup physics
 */
class TracerTag : public Physics::Register<TracerTag>
{
public:
    static const std::string identifier() { return "TracerTag"; }

    explicit TracerTag(CFDSim& sim);

    virtual ~TracerTag() = default;
    //! Initialize the temperature and velocity fields for BoussinesqBubble
    //! simulations
    void initialize_fields(int level, const amrex::Geometry& geom) override;

    void post_init_actions() override {}

    void post_regrid_actions() override {}

    void pre_advance_work() override;

    void post_advance_work() override {}

private:
    CFDSim& m_sim;
    Field& m_act_src;
    Field* m_tracer{nullptr};
    amrex::Real m_src_threshold{0.1};
};

} // namespace amr_wind

#endif /* TracerTag_H */
