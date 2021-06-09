#include "amr-wind/immersed_boundary/bluff_body/bluff_body_ops.H"
#include "amr-wind/core/MultiParser.H"
#include "amr-wind/utilities/ncutils/nc_interface.H"
#include "amr-wind/utilities/io_utils.H"

#include "AMReX_ParmParse.H"

namespace amr_wind {
namespace ib {
namespace bluff_body {

void read_inputs(
    BluffBodyBaseData& wdata, IBInfo&, const ::amr_wind::utils::MultiParser& pp)
{
    pp.query("has_wall_model", wdata.has_wall_model);
    pp.query("is_moving", wdata.is_moving);
    pp.queryarr("vel_bc", wdata.vel_bc);
}

void init_data_structures(BluffBodyBaseData&) {}

void prepare_netcdf_file(
    const std::string& ncfile,
    const BluffBodyBaseData& meta,
    const IBInfo& info)
{
    amrex::ignore_unused(ncfile, meta, info);
}

void write_netcdf(
    const std::string& ncfile,
    const BluffBodyBaseData& meta,
    const IBInfo& info,
    const amrex::Real time)
{
    amrex::ignore_unused(ncfile, meta, info, time);
}

} // namespace bluff_body
} // namespace ib
} // namespace amr_wind
