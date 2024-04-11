#include "amr-wind/boundary_conditions/scalar_bcs.H"

namespace amr_wind::scalar_bc {
void register_scalar_dirichlet(
    Field& field,
    const amrex::AmrCore& mesh,
    const SimTime& time,
    const std::pair<const std::string, const std::string>& udfs)
{
    const std::string inflow_udf = udfs.first;
    const std::string wall_udf = udfs.second;

    if ((inflow_udf == "ConstDirichlet") && (wall_udf == "ConstDirichlet")) {
        return;
    }

    if (wall_udf != "ConstDirichlet") {
        amrex::Abort(
            "Scalar BC: Only constant dirichlet supported for Wall BC");
    }

    register_inflow_scalar_dirichlet<ConstDirichlet>(
        field, inflow_udf, mesh, time);
}
} // namespace amr_wind::scalar_bc
