#ifndef IB_H
#define IB_H

#include "amr-wind/core/Physics.H"

/** \defgroup immersed boundary Immersed boundary module
 *  Immersed boundary source modeling for wind turbines
 *  \ingroup wind
 */

namespace amr_wind {

class Field;

namespace ib {

class ImmersedBoundaryModel;
class IBContainer;

/** Immersed boundary modeling for wind turbines.
 *
 *  \ingroup immersed boundary
 *
 *  This class provides an interface to model wind turbines in a wind farm as
 *  actuator line/disk source terms.
 *
 *  \sa ImmersedBoundaryModel, IBContainer
 */
class IB : public Physics::Register<IB>
{
public:
    static const std::string identifier() { return "IB"; }

    explicit IB(CFDSim&);

    virtual ~IB();

    virtual void initialize_fields(int, const amrex::Geometry&) override {}

    void pre_init_actions() override;

    void post_init_actions() override;

    void post_regrid_actions() override;

    void pre_advance_work() override;

    void post_advance_work() override;

protected:
    //! Total number of actuator components (e.g., turbines) in the flowfield
    int num_ibs() const { return m_ibs.size(); }

    virtual void prepare_outputs();

private:
    void setup_container();

    void update_positions();

    void update_velocities();

    void compute_forces();

    void compute_source_term();

    CFDSim& m_sim;

    Field& m_ib_source;

    std::vector<std::unique_ptr<ImmersedBoundaryModel>> m_ibs;

    std::unique_ptr<IBContainer> m_container;
};

} // namespace ib
} // namespace amr_wind

#endif /* IB_H */
