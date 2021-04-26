#ifndef IB_H
#define IB_H

#include "amr-wind/core/Physics.H"

/** \defgroup ib immersed boundary source terms
 *  Immersed boundary source modeling for turbine
 *  supporting structures (e.g. tower, nacelle, floater)
 *  \ingroup wind
 */

namespace amr_wind {

class Field;

namespace ib {

class IBModel;
class IBContainer;

/** Bluff bodies modeling for wind turbines.
 *
 *  \ingroup ibm
 *
 *  This class provides an interface to model solid bluf bodies in a wind farm
 *  through a direct forcing source terms.
 *
 *  \sa IBModel, IBContainer
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
    //! Total number of ib components (e.g., towers, nacelles) in the
    //! flowfield
    int num_ib() const { return m_ib.size(); }

    virtual void prepare_outputs();

private:
    void setup_container();

    void update_positions();

    void update_velocities();

    void compute_forces();

    void compute_source_term();

    CFDSim& m_sim;

    Field& m_IB_source;

    std::vector<std::unique_ptr<IBModel>> m_ibmodel;

    std::unique_ptr<IBContainer> m_ibcontainer;
};

} // namespace ib
} // namespace amr_wind

#endif /* IB_H */
