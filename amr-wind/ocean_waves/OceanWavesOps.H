#ifndef OCEAN_WAVES_OPS_H
#define OCEAN_WAVES_OPS_H

#include "amr-wind/ocean_waves/OceanWavesTypes.H"
#include "AMReX_Vector.H"
#include "amr-wind/CFDSim.H"

namespace amr_wind::ocean_waves::ops {

/** Read user inputs from an input file for a given immersed boundary object.
 *
 *  \ingroup ocean waves
 *
 *  This operator is called within Physics::pre_init_actions() and, therefore,
 *  does not have access to mesh data. Actions that require mesh information for
 *  initialization should be deferred until ops::init_data_structures call.
 */
template <typename WaveTheoryTrait, typename = void>
struct ReadInputsOp;

/** Perform one-time initialization of the immersed boundary object.
 *
 *  \ingroup ocean waves
 *
 *  This operator is responsible for allocating memory ensuring that the
 *  ocean waves instance has been initialized and ready to be used by the
 * ocean waves physics instance. This function is called once during
 * IB::post_init_actions
 *
 */
template <typename WaveTheoryTrait, typename = void>
struct InitDataOp;

template <typename WaveTheoryTrait, typename = void>
struct UpdateRelaxZonesOp;

template <typename WaveTheoryTrait, typename = void>
struct ApplyRelaxZonesOp;

template <typename WaveTheoryTrait, typename = void>
struct ProcessOutputsOp;

} // namespace amr_wind::ocean_waves::ops

#endif /* OCEAN_WAVES_OPS_H */
