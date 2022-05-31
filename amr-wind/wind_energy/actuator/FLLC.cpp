#include "amr-wind/wind_energy/actuator/FLLC.H"

namespace amr_wind {
namespace actuator {

void FLLCData::init_data(const ComponentView& view, const amrex::Real eps_chord)
{
    const int npts = view.pos.size();
    dx.resize(npts);
    optimal_epsilon.resize(npts);
    les_velocity.assign(npts, vs::Vector::zero());
    optimal_velocity.assign(npts, vs::Vector::zero());
    correction_velocity.assign(npts, vs::Vector::zero());
    lift.assign(npts, vs::Vector::zero());
    grad_lift.assign(npts, vs::Vector::zero());

    for (int i = 0; i < npts - 1; ++i) {
        dx[i] = vs::mag(view.pos[i + 1] - view.pos[i]);
        optimal_epsilon[i] = view.chord[i] * eps_chord;
    }
    dx[npts - 1] = vs::mag(view.pos[npts - 1] - view.pos[npts - 2]);
    optimal_epsilon[npts - 1] = view.chord[npts - 1] * eps_chord;
}

void FLLCParse(const utils::ActParser& pp, FLLCData& data)
{
    pp.query("epsilon", data.epsilon);
    pp.query("fllc_relaxation_factor", data.relaxation_factor);

    if (!pp.contains("epsilon") || !pp.contains("epsilon_chord")) {
        amrex::Abort(
            "Actuators using the filtered lifting line correction (FLLC) "
            "require specification 'epsilon' and 'epsilon_chord'");
    }
}

} // namespace actuator
} // namespace amr_wind