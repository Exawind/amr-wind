#include "amr-wind/wind_energy/actuator/wing/wing_ops.H"
#include "amr-wind/wind_energy/actuator/ActParser.H"

#include "AMReX_ParmParse.H"

namespace amr_wind {
namespace actuator {
namespace wing {

void read_inputs(WingBaseData& wdata, ActInfo& info, const utils::ActParser& pp)
{
    pp.get("num_points", wdata.num_pts);
    pp.get("start", wdata.start);
    pp.get("end", wdata.end);
    pp.get("epsilon", wdata.eps_inp);
    pp.get("pitch", wdata.pitch);
    pp.query("normal", wdata.normal);

    amrex::Real max_eps =
        *std::max_element(wdata.eps_inp.begin(), wdata.eps_inp.end());
    amrex::Real search_radius = max_eps * 3.0;
    const auto& p1 = wdata.start;
    const auto& p2 = wdata.end;
    // clang-format off
    info.bound_box = amrex::RealBox(
        amrex::min(p1.x(), p2.x()) - search_radius,
        amrex::min(p1.y(), p2.y()) - search_radius,
        amrex::min(p1.z(), p2.z()) - search_radius,
        amrex::max(p1.x(), p2.x()) + search_radius,
        amrex::max(p1.y(), p2.y()) + search_radius,
        amrex::max(p1.z(), p2.z()) + search_radius
    );
    // clang-format on
}

void init_data_structures(WingBaseData& wdata, ActGrid& grid)
{
    int npts = wdata.num_pts;
    grid.resize(npts);

    // Wing span
    auto wspan = wdata.end - wdata.start;
    // Compute chord/flow direction as a cross-product
    auto chord = (wspan ^ wdata.normal);
    // Set up global to local transformation matrix
    auto tmat = vs::Tensor(chord.unit(), wspan.unit(), wdata.normal.unit());

    // Equal spacing along span
    auto dx = (1.0 / static_cast<amrex::Real>(npts - 1)) * wspan;

    for (int i = 0; i < npts; ++i) {
        grid.pos[i] = wdata.start + static_cast<amrex::Real>(i) * dx;
        grid.epsilon[i] = wdata.eps_inp;
        grid.orientation[i] = tmat;
    }

    // Initialize remaining data
    grid.force.assign(npts, vs::Vector::zero());
    grid.vel_pos.assign(grid.pos.begin(), grid.pos.end());
    grid.vel.assign(npts, vs::Vector::zero());

    // Assign length of actuator segments
    grid.dx.assign(npts, vs::mag(dx));
    // The first and last segments have half width
    grid.dx.front() *= 0.5;
    grid.dx.back() *= 0.5;
}

} // namespace wing
} // namespace actuator
} // namespace amr_wind
