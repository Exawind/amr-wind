#ifndef FLLC_H
#define FLLC_H

namespace amr_wind {
namespace actuator {


/** @brief This struct will operate on a blade/wing.
 *  The velocity from the simulation is corrected using the
 *  Filtered Lifting Line Theory Correction from Martinez and Meneveau 2019
 *
*/

template <typename ActTrait, typename SrcTrait>
struct FLLCOpBase //<ActTrait, SrcTrait>
{ 

    //! The velocity correction from filtered lifting line theory
    VecList u_les_fllc;
    VecList u_opt_fllc;
    VecList G_fllc;
    VecList dG_fllc;

    FLLCOpBase(typename ActTrait::DataType& data){

        auto& wdata = data.meta();
        const int npts = wdata.num_pts;

         G_fllc.assign(npts, vs::Vector::zero());
         dG_fllc.assign(npts, vs::Vector::zero());
         u_les_fllc.assign(npts, vs::Vector::zero());
         u_opt_fllc.assign(npts, vs::Vector::zero());

    }

   void operator()(typename ActTrait::DataType& data) {

        auto& grid = data.grid();
        auto& wdata = data.meta();
        const int npts = wdata.num_pts;
        auto& du_fllc = wdata.du_fllc;

        /**
        * Step 1
        * Compute the lift force distribution (G)
        * Compute equation 5.3 from Martinez-Tossas and Meneveau 2019
        */
        for (int ip = 0; ip < npts; ++ip) {

            const auto& tmat = grid.orientation[ip];
            auto force = grid.force[ip];
            auto vel =  wdata.vel_rel[ip] & tmat;
            auto dx =  wdata.dx[ip];
            const auto vmag = vs::mag(vel);
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
        dG_fllc[npts-1] = -1 * G_fllc[npts-1];
        for (int ip = 1; ip < npts-1; ++ip) 
            dG_fllc[ip] = 0.5 * (G_fllc[ip+1] - G_fllc[ip - 1]);

        /**
        * Step 3
        * Compute the induced velocities
        * Compute equations 5.6 and 5.7 from Martinez-Tossas and Meneveau 2019
        */
        for (int ip = 0; ip < npts; ++ip) {

            const auto eps_les = wdata.eps_inp[0];
            const auto eps_opt = wdata.epsilon_chord[0] * wdata.chord[ip];

            for (int jp = 0; jp < npts; ++jp) {

                if (ip==jp) continue;

                const auto dr = vs::mag(grid.vel_pos[ip] - grid.vel_pos[jp]);
                const auto& vel =  wdata.vel_rel[jp];
                const auto vmag = vs::mag(vel);

                auto coefficient = 1.0 / (-4.0 * amr_wind::utils::pi() * dr * vmag);
                const auto cLes = 1.0 - std::exp(-dr * dr / (eps_les * eps_les));
                const auto cOpt = 1.0 - std::exp(-dr * dr / (eps_opt * eps_opt));

                // The sign of the induced velocity depends on which side of the blade we are on
                if (ip<jp) coefficient*=-1;

               u_les_fllc[ip] = u_les_fllc[ip] - dG_fllc[jp] * coefficient * cLes;
               u_opt_fllc[ip] = u_opt_fllc[ip] - dG_fllc[jp] * coefficient * cOpt;

            }

            // Relaxation to compute the induced velocity
            const double f = 0.1;
            du_fllc[ip] = (1.-f) * du_fllc[ip] + f * (u_opt_fllc[ip] - u_les_fllc[ip]);

        }

        /**
        * Step 4
        * Compute the induced velocity difference
        * Compute equation 5.8 from Martinez-Tossas and Meneveau 2019
        */
        for (int ip = 0; ip < npts; ++ip) 
            grid.vel[ip] = grid.vel[ip] + du_fllc[ip];


    }
};

template <typename ActTrait, typename SrcTrait>
struct FLLCOpTurbineBase 
{


   FLLCOpTurbineBase(typename ActTrait::DataType& data){}

   void operator()(typename ActTrait::DataType& data) {

    // The cfd grid data
    auto& grid = data.grid();
    // The turbine data
    auto& tdata = data.meta();

    // Introduce FLLC logic
    for (int i=0; i < tdata.num_blades; ++i)
    {

    }
    }
};




}
}

#endif /* FLLC_H */
