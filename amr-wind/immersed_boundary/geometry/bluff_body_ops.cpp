#include "amr-wind/immersed_boundary/geometry/bluff_body_ops.H"
#include "amr-wind/immersed_boundary/IBParser.H"
#include "amr-wind/utilities/ncutils/nc_interface.H"
#include "amr-wind/utilities/io_utils.H"

#include "AMReX_ParmParse.H"

namespace amr_wind {
namespace ib {
namespace bluff_body {

void read_inputs(BluffBodyBaseData& wdata, IBInfo&, const utils::IBParser& pp)
{
    pp.get("num_points", wdata.num_pts);
}

void init_data_structures(BluffBodyBaseData&, IBGrid&) {}

void prepare_netcdf_file(
    const std::string& ncfile,
    const BluffBodyBaseData& meta,
    const IBInfo& info,
    const IBGrid& grid)
{
    amrex::ignore_unused(ncfile, meta, info, grid);
}

void write_netcdf(
    const std::string& ncfile,
    const BluffBodyBaseData& meta,
    const IBInfo& info,
    const IBGrid& grid,
    const amrex::Real time)
{
    amrex::ignore_unused(ncfile, meta, info, grid, time);
}

} // namespace bluff_body
} // namespace ib
} // namespace amr_wind
