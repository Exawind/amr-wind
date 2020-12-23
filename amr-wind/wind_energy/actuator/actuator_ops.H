#ifndef ACT_OPS_H
#define ACT_OPS_H

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
 *  \param Data object for the specific actuator instance
 */
template <typename T>
void read_inputs(typename T::DataType&);

/** Perform one-time initialization of the actuator object
 *
 *  This function is responsible for allocating memory ensuring that the
 *  actuator instance has been initialized and ready to be used by the Actuator
 *  physics instance. This function is called once during
 * Actuator::post_init_actions
 *
 *  \tparam T An actuator traits type
 *  \param Data object for the specific actuator instance
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

} // namespace ops
} // namespace actuator
} // namespace amr_wind

#include "amr-wind/wind_energy/actuator/actuator_opsI.H"

#endif /* ACT_OPS_H */
