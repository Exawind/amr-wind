#ifndef FLLC_H
#define FLLC_H

#include "amr-wind/wind_energy/actuator/actuator_types.H"
#include "amr-wind/wind_energy/actuator/FLLC.H"

namespace amr_wind {
namespace actuator {

/** @brief This struct will operate on a blade/wing.
 *  The velocity from the simulation is corrected using the
 *  Filtered Lifting Line Theory Correction from Martinez and Meneveau 2019
 *
 */

struct FLLCOp
{
    void operator()(ComponentView& data, FLLCData& fllc)
    {

        const int npts = fllc.correction_velocity.size();
        auto& du_fllc = fllc.correction_velocity;

        auto& u_les_fllc = fllc.les_velocity;
        auto& u_opt_fllc = fllc.optimal_velocity;
        auto& G_fllc = fllc.lift;
        auto& dG_fllc = fllc.grad_lift;

        /**
         * Step 1
         * Compute the lift force distribution (G)
         * Compute equation 5.3 from Martinez-Tossas and Meneveau 2019
         */
        for (int ip = 0; ip < npts; ++ip) {

            const auto& tmat = data.orientation[ip];
            const auto force = data.force[ip];
            const auto vel = data.vel_rel[ip] & tmat;
            const auto dx = fllc.dx[ip];
            const auto vmag =
                std::max(vs::mag(vel), vs::DTraits<amrex::Real>::eps());
            const auto vmag2 = vmag * vmag;

            const auto fv = force & vel;

            G_fllc[ip] = (force - vel * fv / vmag2) / dx;
        }

        /**
         * Step 2
         * Compute gradient of the lift force distribution (\Delta G)
         * Compute equations 5.4 and 5.5 from Martinez-Tossas and Meneveau 2019
         */
        dG_fllc[0] = G_fllc[0];
        dG_fllc[npts - 1] = -1 * G_fllc[npts - 1];
        for (int ip = 1; ip < npts - 1; ++ip)
            dG_fllc[ip] = 0.5 * (G_fllc[ip + 1] - G_fllc[ip - 1]);

        /**
         * Step 3
         * Compute the induced velocities
         * Compute equations 5.6 and 5.7 from Martinez-Tossas and Meneveau 2019
         */
        for (int ip = 0; ip < npts; ++ip) {

            const auto eps_les = fllc.epsilon;
            const auto eps_opt = fllc.optimal_epsilon[ip];

            for (int jp = 0; jp < npts; ++jp) {

                if (ip == jp) continue;

                const auto dr = vs::mag(data.vel_pos[ip] - data.vel_pos[jp]);
                const auto& vel = data.vel_rel[jp];
                const auto vmag =
                    std::max(vs::mag(vel), vs::DTraits<amrex::Real>::eps());

                auto coefficient =
                    1.0 / (-4.0 * amr_wind::utils::pi() * dr * vmag);
                const auto cLes =
                    1.0 - std::exp(-dr * dr / (eps_les * eps_les));
                const auto cOpt =
                    1.0 - std::exp(-dr * dr / (eps_opt * eps_opt));

                // The sign of the induced velocity depends on which side of the
                // blade we are on
                if (ip < jp) coefficient *= -1;

                u_les_fllc[ip] =
                    u_les_fllc[ip] - dG_fllc[jp] * coefficient * cLes;
                u_opt_fllc[ip] =
                    u_opt_fllc[ip] - dG_fllc[jp] * coefficient * cOpt;
            }

            // Relaxation to compute the induced velocity
            const double f = fllc.relaxation_factor;
            du_fllc[ip] =
                (1. - f) * du_fllc[ip] + f * (u_opt_fllc[ip] - u_les_fllc[ip]);
        }

        /**
         * Step 4
         * Compute the induced velocity difference
         * Compute equation 5.8 from Martinez-Tossas and Meneveau 2019
         */
        for (int ip = 0; ip < npts; ++ip)
            data.vel[ip] = data.vel[ip] + du_fllc[ip];
    }
};

} // namespace actuator
} // namespace amr_wind

#endif /* FLLC_H */
