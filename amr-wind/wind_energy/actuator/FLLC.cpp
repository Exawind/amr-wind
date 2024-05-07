#include "amr-wind/wind_energy/actuator/FLLC.H"
#include "amr-wind/utilities/linear_interpolation.H"

namespace amr_wind::actuator {

void FLLCInit(
    FLLCData& data, const ComponentView& view, const amrex::Real eps_chord)
{
    const int npts = static_cast<int>(view.pos.size());
    data.different_sizes = view.pos.size() != view.vel_pos.size();
    if (data.different_sizes) {
        data.span_distance_force.resize(view.pos.size());
        data.span_distance_vel.resize(view.vel_pos.size());
        for (size_t i = 0; i < view.pos.size(); ++i) {
            data.span_distance_force[i] = vs::mag(view.pos[i]);
        }
        for (size_t i = 0; i < view.vel_pos.size(); ++i) {
            data.span_distance_vel[i] = vs::mag(view.vel_pos[i]);
        }
    }
    data.dr.resize(npts);
    data.optimal_epsilon.resize(npts);
    data.les_velocity.assign(npts, vs::Vector::zero());
    data.optimal_velocity.assign(npts, vs::Vector::zero());
    data.correction_velocity.assign(npts, vs::Vector::zero());
    data.lift.assign(npts, vs::Vector::zero());
    data.grad_lift.assign(npts, vs::Vector::zero());
    data.force_point_velocity.assign(npts, vs::Vector::zero());

    for (int i = 0; i < npts; ++i) {

        data.optimal_epsilon[i] = view.chord[i] * eps_chord;
    }

    // Central difference
    for (int i = 1; i < npts - 1; ++i) {
        data.dr[i] = vs::mag(view.pos[i + 1] - view.pos[i - 1]) / 2.;
    }
    data.dr[0] = vs::mag(view.pos[1] - view.pos[0]) / 2.;
    data.dr[npts - 1] = vs::mag(view.pos[npts - 1] - view.pos[npts - 2]) / 2.;

    namespace interp = ::amr_wind::interp;

    // Non-uniform radial distribution variables
    size_t idx = 0;
    amrex::Real dr;
    amrex::Real dr1;
    amrex::Real eps1;
    amrex::Real dr2;
    amrex::Real eps2;

    // Non-uniform radial distribution calc
    data.r_n.push_back(
        data.span_distance_vel[0]); // Initialize the first location
    while (data.r_n.back() < data.span_distance_vel.back()) {
        amrex::Real r_ = data.r_n.back(); // The latest radial position
        // Interpolate the value of epsilon to the current radial location
        eps1 = interp::linear(data.span_distance_vel, data.optimal_epsilon, r_);
        dr1 = eps1 / data.eps_dr; // spacing needed to match condition
        // Same calculation as before but for the new location
        eps2 = interp::linear(
            data.span_distance_vel, data.optimal_epsilon, r_ + dr1);
        dr2 = eps2 / data.eps_dr;
        // This will ensure that the spacing will always follow the constraint
        dr = std::min(dr1, dr2);
        assert(dr > 0.);
        data.r_n.push_back(std::min(r_ + dr, data.span_distance_vel.back()));
        idx++;
    }

    // Central difference for non-uniform points
int npts_n = data.r_n.size();
data.dr_n.resize(npts_n);
    for (int i = 1; i < npts_n - 1; ++i) {
        data.dr_n[i] = std::abs(data.r_n[i + 1] - data.r_n[i - 1]) / 2.;
    }
    data.dr_n[0] = std::abs(data.r_n[1] - data.r_n[0]) / 2.;
    data.dr_n[npts_n - 1] = std::abs(data.r_n[npts_n - 1] - data.r_n[npts_n - 2]) / 2.;


    data.initialized = true;
}

void FLLCParse(const utils::ActParser& pp, FLLCData& data)
{
    pp.query("epsilon", data.epsilon);
    pp.query("fllc_relaxation_factor", data.relaxation_factor);
    pp.query("fllc_start_time", data.fllc_start_time);
    std::string typeString = "variable_chord";
    pp.query("fllc_type", typeString);
    data.correction_type = FLLCTypeMap.at(typeString);
    pp.query("non-uniform", data.nonu);
    pp.query("eps_dr", data.eps_dr);

    if (!pp.contains("epsilon") || !pp.contains("epsilon_chord")) {
        amrex::Abort(
            "Actuators using the filtered lifting line correction (FLLC) "
            "require specification 'epsilon' and 'epsilon_chord'");
    }
}

} // namespace amr_wind::actuator
