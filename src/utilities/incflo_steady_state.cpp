#include <incflo.H>

using namespace amrex;

//
// Check if steady state has been reached by verifying that
//
//      max(abs( u^(n+1) - u^(n) )) / dt < tol
//      max(abs( v^(n+1) - v^(n) )) / dt < tol
//      max(abs( w^(n+1) - w^(n) )) / dt < tol
//
//      OR
//
//      sum(abs( u^(n+1) - u^(n) )) / sum(abs( u^(n) )) < tol
//      sum(abs( v^(n+1) - v^(n) )) / sum(abs( v^(n) )) < tol
//      sum(abs( w^(n+1) - w^(n) )) / sum(abs( w^(n) )) < tol
//
bool incflo::SteadyStateReached()
{

    amrex::Abort("xxxxx TODO: SteadyStateReached");
    return false;
#if 0
    BL_PROFILE("incflo::SteadyStateReached()");

    int condition1[finest_level + 1];
    int condition2[finest_level + 1];

    // Make sure velocity is up to date
    incflo_set_velocity_bcs(m_cur_time, vel);

    // Use temporaries to store the difference between current and previous solution
    Vector<std::unique_ptr<MultiFab>> diff_vel;
    diff_vel.resize(finest_level + 1);
    for(int lev = 0; lev <= finest_level; lev++)
    {
#ifdef AMREX_USE_EB
        diff_vel[lev].reset(new MultiFab(grids[lev], dmap[lev], AMREX_SPACEDIM, 0, MFInfo(), *ebfactory[lev]));
#else
        diff_vel[lev].reset(new MultiFab(grids[lev], dmap[lev], AMREX_SPACEDIM, 0, MFInfo()));
#endif
        MultiFab::LinComb(*diff_vel[lev], 1.0, *vel[lev], 0, -1.0, *vel_o[lev], 0, 0, AMREX_SPACEDIM, 0);

        Real max_change = 0.0;
        Real max_relchange = 0.0;
        // Loop over components, only need to check the largest one
        for(int i = 0; i < AMREX_SPACEDIM; i++)
        {
            // max(abs(u^{n+1}-u^n))
            max_change = amrex::max(max_change, Norm(diff_vel, lev, i, 0));

            // sum(abs(u^{n+1}-u^n)) / sum(abs(u^n))
            // TODO: this gives zero often, check for bug
            Real norm1_diff = Norm(diff_vel, lev, i, 1);
            Real norm1_old = Norm(vel_o, lev, i, 1);
            Real relchange = norm1_old > 1.0e-15 ? norm1_diff / norm1_old : 0.0;
            max_relchange = amrex::max(max_relchange, relchange);
        }

        condition1[lev] = (max_change < m_steady_state_tol * dt);
        condition2[lev] = (max_relchange < m_steady_state_tol);

        // Print out info on steady state checks
        if (m_verbose > 0)
        {
            amrex::Print() << "\nSteady state check level " << lev << std::endl; 
            amrex::Print() << "||u-uo||/||uo|| = " << max_relchange
                           << ", du/dt  = " << max_change/dt << std::endl;
        }
    }

    bool reached = true;
    for(int lev = 0; lev <= finest_level; lev++)
    {
        reached = reached && (condition1[lev] || condition2[lev]);
    }

    // Always return negative to first access. This way
    // initial zero velocity field do not test for false positive
    if(m_nstep < 2)
    {
        return false;
    } else {
        return reached;
    }
#endif
}
