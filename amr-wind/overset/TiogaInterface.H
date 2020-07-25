#ifndef TIOGAINTERFACE_H
#define TIOGAINTERFACE_H

#include "amr-wind/overset/OversetManager.H"

namespace amr_wind {

class IntField;

class TiogaInterface : public OversetManager::Register<TiogaInterface>
{
public:
    static std::string identifier() { return "TIOGA"; }

    TiogaInterface(CFDSim& sim);

    /** Perform one-time initialization actions.
     *
     *  This method is called after the initial mesh and solution fields have
     *  been created.
     */
    virtual void post_init_actions() override;

    /** Perform necessary updates after a regrid.
     */
    virtual void post_regrid_actions() override;

    /** Perform actions before initializing connectivity through overset
     * library.
     *
     *  This method is called before performing connectivity. The implementation
     *  should update all necessary mesh data and metadata that will be used for
     *  overset connectivity.
     */
    virtual void pre_overset_conn_work() override;

    /** Updates necessary after overset connectivity has been updated.
     *
     *  Update dependent fields or data structures as a result of changes in
     *  overset connectivity.
     */
    virtual void post_overset_conn_work() override;

    /** Perform actions necessary to register fields for overset exchange.
     *
     *  This method will be called at least once every timestep and, perhaps,
     *  multiple times within a timestep.
     */
    virtual void register_solution() override;

    /** Additional work after solution fields have been exchanged.
     */
    virtual void update_solution() override;

private:
    //! IBLANK on cell centered fields
    IntField& m_iblank_cell;

    //! IBLANK for nodal fields
    IntField& m_iblank_node;

    //! AMReX mask for linear system solves
    IntField& m_mask_cell;

    IntField& m_mask_node;
};

} // namespace amr_wind

#endif /* TIOGAINTERFACE_H */
