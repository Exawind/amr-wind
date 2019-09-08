#include <incflo.H>

#include <cmath>
#include <limits>

using namespace std;

//
// Compute new dt by using the formula derived in
// "A Boundary Condition Capturing Method for Multiphase Incompressible Flow"
// by Kang et al. (JCP).
//
//  dt/2 * ( C+V + sqrt( (C+V)**2 + 4Fx/dx + 4Fy/dy + 4Fz/dz )
//
// where
//
// C = max(|U|)/dx + max(|V|)/dy + max(|W|)/dz    --> Convection
//
// V = 2 * max(eta/ro) * (1/dx^2 + 1/dy^2 +1/dz^2) --> Diffusion
//
// Fx, Fy, Fz = net acceleration due to external forces
//
// WARNING: We use a slightly modified version of C in the implementation below
//
void incflo::ComputeDt(int initialisation)
{
	BL_PROFILE("incflo::ComputeDt");

	// Compute dt for this time step
	Real umax = 0.0;
	Real vmax = 0.0;
	Real wmax = 0.0;
	Real romin = 1.e20;
	Real etamax = 0.0;

    for(int lev = 0; lev <= finest_level; lev++)
    {
        // The functions take the min/max over uncovered cells 
        umax   = amrex::max(umax,   Norm(vel, lev, 0, 0));
        vmax   = amrex::max(vmax,   Norm(vel, lev, 1, 0));
        wmax   = amrex::max(wmax,   Norm(vel, lev, 2, 0));
        romin  = amrex::min(romin,  Norm( ro, lev, 0, 0));
        etamax = amrex::max(etamax, Norm(eta, lev, 0, 0));
    }

    const Real* dx = geom[finest_level].CellSize();
    Real idx = 1.0 / dx[0];
    Real idy = 1.0 / dx[1];
    Real idz = 1.0 / dx[2];

    // Convective term
    Real conv_cfl = std::max(std::max(umax * idx, vmax * idy), wmax * idz);

    // Viscous term
    Real diff_cfl = 2.0 * etamax / romin * (idx * idx + idy * idy + idz * idz);

    // Forcing term
    Real forc_cfl = std::abs(gravity[0] - std::abs(gp0[0])) * idx
                  + std::abs(gravity[1] - std::abs(gp0[1])) * idy
                  + std::abs(gravity[2] - std::abs(gp0[2])) * idz;

    // Combined CFL conditioner
    Real comb_cfl = conv_cfl + diff_cfl + sqrt(pow(conv_cfl + diff_cfl, 2) + 4.0 * forc_cfl);

    // Update dt
    Real dt_new = 2.0 * cfl / comb_cfl;

    // Reduce CFL for initial step
    if(initialisation)
    {
        dt_new *= 0.1;
    }

    // Protect against very small comb_cfl
    // This may happen, for example, when the initial velocity field
    // is zero for an inviscid flow with no external forcing
    Real eps = std::numeric_limits<Real>::epsilon();
    if(comb_cfl <= eps)
    {
        dt_new = 0.5 * dt;
    }

    // Don't let the timestep grow by more than 10% per step.
    if(dt > 0.0 && last_plt != nstep)
    {
        dt_new = amrex::min(dt_new, 1.1 * dt);
    }
    
    // Don't overshoot specified plot times
    if(plot_per > 0.0 && 
            (trunc((cur_time + dt_new + eps) / plot_per) > trunc((cur_time + eps) / plot_per)))
    {
        dt_new = trunc((cur_time + dt_new) / plot_per) * plot_per - cur_time;
    }

    // Don't overshoot the final time if not running to steady state
    if(!steady_state && stop_time > 0.0)
    {
        if(cur_time + dt_new > stop_time)
        {
            dt_new = stop_time - cur_time;
        }
    }

    // Make sure the timestep is not set to zero after a plot_per stop
    if(dt_new < eps)
    {
        dt_new = 0.5 * dt;
    }

    // If using fixed time step, check CFL condition and give warning if not satisfied
	if(fixed_dt > 0.0)
	{
		if(dt_new < fixed_dt)
		{
			amrex::Print() << "WARNING: fixed_dt does not satisfy CFL condition: \n"
						   << "max dt by CFL     : " << dt_new << "\n"
						   << "fixed dt specified: " << fixed_dt << std::endl;
		}
		dt = fixed_dt;
	}
	else
	{
		dt = dt_new;
	}
}
