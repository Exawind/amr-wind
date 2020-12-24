#ifndef ACT_OPS_H
#define ACT_OPS_H

#include "amr-wind/wind_energy/actuator/actuator_types.H"
#include "amr-wind/wind_energy/actuator/ActParser.H"
#include "AMReX_Vector.H"

namespace amr_wind {
namespace actuator {
namespace ops {

/** Read user inputs from an input file for a given actuator object
 *
 *  This function is called within Physics::pre_init_actions() and, therefore,
 *  does not have access to mesh data. Actions that require mesh information for
 *  initialization should be deferred until ops::init_data_structures call.
 *
 *  \tparam T An actuator traits type
 *  \param data Data object for the specific actuator instance
 */
template <typename T>
void read_inputs(typename T::DataType&, const utils::ActParser&);

/** Perform one-time initialization of the actuator object
 *
 *  This function is responsible for allocating memory ensuring that the
 *  actuator instance has been initialized and ready to be used by the Actuator
 *  physics instance. This function is called once during
 * Actuator::post_init_actions
 *
 *  \tparam T An actuator traits type
 *  \param data Data object for the specific actuator instance
 */
template <typename T>
void init_data_structures(typename T::DataType&);

/** Determine the list of processes (MPI ranks) that is expected to be
 *  influenced by this actuator source.
 *
 *  The list of processes is used to determine where velocity sampling data is
 *  too be gathered, as well as ensuring that all actuator data is properly
 *  mapped.
 */
template <typename T>
void determine_influenced_procs(typename T::DataType&);

/** During initialization determine the influenced processes and elect a root
 *  process to manage I/O and external solver interfaces for this turbine.
 *
 *  \tparam T An actuator traits type
 *  \param  data Data object for the specific actuator instance
 *
 *  \param act_proc_count A map of the number of turbines managed by each proc
 */
template <typename T>
void determine_root_proc(typename T::DataType&, amrex::Vector<int>&);

/** Update the positions of the actuator nodes
 *
 *  After this call, the interface expects updated position vectors in the
 *  vel_pos array of the actuator grid instance. These positions will be used
 *  the sample the velocity field during the simulation.
 *
 *  \tparam T An actuator traits type
 *  \param  data Data object for the specific actuator instance
 */
template <typename T>
void update_positions(typename T::DataType&);

/** Process updated velocities from CFD flow field
 *
 *  The actuator physics interface will populate the ``vel`` array of the
 *  actuator grid instance before calling this function. The updated velocities
 *  can be further processed by specializing this function.
 *
 *  \tparam T An actuator traits type
 *  \param  data Data object for the specific actuator instance
 */
template <typename T>
void update_velocities(typename T::DataType&);

/** Compute aerodynamic forces at the actuator grid points during a simulation.
 *
 *  \tparam T An actuator traits type
 *  \param  data Data object for the specific actuator instance
 */
template <typename T>
void compute_forces(typename T::DataType&);

} // namespace ops
} // namespace actuator
} // namespace amr_wind

#include "amr-wind/wind_energy/actuator/actuator_opsI.H"

#endif /* ACT_OPS_H */
