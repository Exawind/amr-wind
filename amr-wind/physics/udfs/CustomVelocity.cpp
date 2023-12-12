#include "amr-wind/physics/udfs/CustomVelocity.H"
#include "amr-wind/core/Field.H"
#include "amr-wind/core/FieldRepo.H"
#include "amr-wind/core/vs/vector.H"
#include "amr-wind/equation_systems/icns/icns.H"

#include "AMReX_ParmParse.H"

namespace amr_wind::udf {

CustomVelocity::CustomVelocity(const Field& fld)
{
    // This is a where the user can set some user defined variables
    // This capability can be activated with the following in the input file:
    // xlo.type = "mass_inflow"
    // xlo.velocity.inflow_type = CustomVelocity
    // CustomVelocity.foo = 1.0
}

} // namespace amr_wind::udf
