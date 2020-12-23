#ifndef ACTUATORMODEL_H
#define ACTUATORMODEL_H

#include "amr-wind/core/Factory.H"
#include "amr-wind/wind_energy/actuator/actuator_types.H"
#include "amr-wind/wind_energy/actuator/actuator_ops.H"

namespace amr_wind {

class CFDSim;

namespace actuator {

/** Abstract representation of an actuator source in the flowfield
 *
 *  This class is intended to represent a single object (e.g., turbine, fixed
 *  wing) that exists in the computational domain. A collection of source
 *  objects (e.g., turbines in a wind farm) is modeled using the
 *  amr_wind::Actuator physics class that holds a collection of instances of
 *  ActuatorModel and acts on them.
 */
class ActuatorModel
    : public Factory<ActuatorModel, CFDSim&, const std::string&, const int>
{
public:
    static std::string base_identifier() { return "ActuatorModel"; }

    virtual ~ActuatorModel() = default;

    virtual const std::string& label() const = 0;

    virtual int id() const = 0;

    virtual const ActInfo& info() const = 0;

    virtual void pre_init_actions() = 0;

    virtual void post_init_actions() = 0;
};

/** Concrete implementation of the ActuatorModel for different actuator types
 *
 *  \tparam ActTrait An actuator type trait that defines the behavior of the
 *  different actions (initialization, velocity sampling, force computation, and
 *  momentum source term computations) that are executed by the ActuatorModel
 *  class through the amr_wind::Actuator physics driver class.
 *
 *  \tparam SrcTrait A source type trait that defines how the forcing term is
 *  computed, e.g., if the turbine rotor is represented as a line or a disk. For
 *  most cases, like the fixed wing case, this will just be line, therefore, it
 *  is set as default for convenience.
 */
template <typename ActTrait, typename SrcTrait = ActSrcLine>
class ActModel : public ActuatorModel::Register<ActModel<ActTrait, SrcTrait>>
{
private:
    //! Instance that holds all data related to a particular actuator type
    ActDataHolder<ActTrait> m_data;

public:
    static std::string identifier()
    {
        return ActTrait::identifier() + SrcTrait::identifier();
    }

    /**
     *  \param sim CFD simulation environment instance
     *  \param label Unique string identifier for this instance
     *  \param id Unique integer tag for this instance
     */
    ActModel(CFDSim& sim, const std::string& label, const int id)
        : m_data(sim, label, id)
    {}

    //! Return the unique identifier (name) for this instance
    const std::string& label() const override { return m_data.m_info.label; }

    //! Return the unique tag (integer ID) for this instance
    int id() const override { return m_data.m_info.id; }

    //! Return the core info object for this actuator instance
    const ActInfo& info() const override { return m_data.m_info; }

    void pre_init_actions() override { ops::read_inputs<ActTrait>(m_data); }

    void post_init_actions() override
    {
        ops::init_data_structures<ActTrait>(m_data);
        ops::determine_influenced_procs<ActTrait>(m_data);
    }
};

} // namespace actuator
} // namespace amr_wind

#endif /* ACTUATORMODEL_H */
