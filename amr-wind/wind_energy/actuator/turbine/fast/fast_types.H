#ifndef FAST_TYPES_H
#define FAST_TYPES_H

#include "amr-wind/wind_energy/actuator/turbine/fast/fast_wrapper.H"
#include "amr-wind/wind_energy/actuator/turbine/external/external_base_types.H"

#include <string>

namespace exw_fast {

/** Representation of a turbine for exchanging data with FAST
 */
struct FastTurbine : public exw_ext_turb::ExternalTurbine
{
    /** Variables exclusive to OpenFAST turbine type
     */

    //! Node cluster type for the chord
    int chord_cluster_type{0};

    //! Does FAST need solution0
    bool is_solution0{true};

    exw_fast::OfInputType to_cfd;
    exw_fast::OfOutputType from_cfd;

#if OPENFAST_VERSION_MAJOR != 4
    exw_fast::SC_DX_InputType to_sc;
    exw_fast::SC_DX_OutputType from_sc;
#endif
};

} // namespace exw_fast

#endif /* FAST_TYPES_H */
