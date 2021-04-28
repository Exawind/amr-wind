#ifndef IBOPS_H
#define IBOPS_H

#include "amr-wind/immersed_boundary/IBTypes.H"
#include "amr-wind/immersed_boundary/IBParser.H"
#include "AMReX_Vector.H"

namespace amr_wind {
namespace ib {
namespace ops {

/** Read user inputs from an input file for a given immersed boundary object.
 *
 *  \ingroup immersed boundary
 *
 *  This operator is called within Physics::pre_init_actions() and, therefore,
 *  does not have access to mesh data. Actions that require mesh information for
 *  initialization should be deferred until ops::init_data_structures call.
 */
template <typename GeomTrait, typename = void>
struct ReadInputsOp;

/** Perform one-time initialization of the immersed boundary object.
 *
 *  \ingroup immersed boundary
 *
 *  This operator is responsible for allocating memory ensuring that the
 *  immersed boundary instance has been initialized and ready to be used by the immersed boundary
 *  physics instance. This function is called once during IB::post_init_actions
 *
 */
template <typename GeomTrait, typename = void>
struct InitDataOp;

/** Update the positions of the immersed boundary nodes.
 *
 *  \ingroup immersed boundary
 *
 *  After this call, the interface expects updated position vectors in the
 *  vel_pos array of the immersed boundary grid instance. These positions will be used
 *  the sample the velocity field during the simulation.
 */
template <typename GeomTrait, typename = void>
struct UpdatePosOp;

/** Process updated velocities from CFD flow field.
 *
 *  \ingroup immersed boundary
 *
 *  The immersed boundary physics interface will populate the `vel` array of the
 *  immersed boundary grid instance before calling this function. The updated velocities
 *  can be further processed by specializing this function.
 */
template <typename GeomTrait, typename = void>
struct UpdateVelOp;

/** Compute aerodynamic forces at the immersed boundary grid points during a simulation.
 *
 *  \ingroup immersed boundary
 */
template <typename GeomTrait, typename = void>
struct ComputeForceOp;

/** Compute source term for the momentum equation.
 *
 *  \ingroup immersed boundary
 */
template <typename GeomTrait, typename = void>
class IBSrcOp;

template <typename GeomTrait, typename = void>
struct ProcessOutputsOp;

/** Determine the list of processes (MPI ranks) that is expected to be
 *  influenced by this immersed boundary source.
 *
 *  The list of processes is used to determine where velocity sampling data is
 *  too be gathered, as well as ensuring that all immersed boundary data is properly
 *  mapped.
 */
template <typename T>
void determine_influenced_procs(typename T::DataType&);

/** During initialization determine the influenced processes and elect a root
 *  process to manage I/O and external solver interfaces for this turbine.
 *
 *  \tparam T An immersed boundary traits type
 *  \param  data Data object for the specific immersed boundary instance
 *
 *  \param ib_proc_count A map of the number of turbines managed by each proc
 */
template <typename T>
void determine_root_proc(typename T::DataType&, amrex::Vector<int>&);

} // namespace ops
} // namespace ib
} // namespace amr_wind

#include "amr-wind/immersed_boundary/IBOpsI.H"

#endif /* IBOPS_H */
