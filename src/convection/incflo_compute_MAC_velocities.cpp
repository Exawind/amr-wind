#include <incflo.H>

#ifdef AMREX_USE_EB
#include <AMReX_EB_utils.H>
#endif

using namespace amrex;

//
// Compute the MAC (normal) velocities on face centroids
//
void
incflo::incflo_compute_MAC_velocities( Vector< std::unique_ptr<MultiFab> >& vel_forces_in,
                                       Vector< std::unique_ptr<MultiFab> >& scal_forces_in,
                                       Vector< std::unique_ptr<MultiFab> >& vel_in,
                                       Vector< std::unique_ptr<MultiFab> >& density_in,
                                       Vector< std::unique_ptr<MultiFab> >& tracer_in,
                                       Real time)
{
    BL_PROFILE("incflo::incflo_compute_MAC_velocities");

    // First do FillPatch of {velocity, density, tracer} so we know the ghost cells of
    // these arrays are all filled
    for (int lev = 0; lev <= finest_level; lev++)
    {

#ifdef AMREX_USE_EB
    const FabFactory<FArrayBox>& factory =  *ebfactory[lev];
#else
    const FabFactory<FArrayBox>& factory = FArrayBoxFactory();
#endif

        // State with ghost cells
        MultiFab Sborder_u(grids[lev], dmap[lev], vel_in[lev]->nComp(), nghost, MFInfo(), factory);
        FillPatchVel(lev, time, Sborder_u);

        // Copy each FAB back from Sborder_u into the vel array, complete with filled ghost cells
        MultiFab::Copy (*vel_in[lev], Sborder_u, 0, 0, vel_in[lev]->nComp(), nghost);

        MultiFab Sborder_r(grids[lev], dmap[lev], 1, nghost, MFInfo(), factory);
        FillPatchDensity(lev, time, Sborder_r);
        MultiFab::Copy (*density_in[lev], Sborder_r, 0, 0, 1, nghost);

        if (advect_tracer)
        {
           MultiFab Sborder_s(grids[lev], dmap[lev], ntrac, nghost, MFInfo(), factory);
           FillPatchScalar(lev, time, Sborder_s);
           MultiFab::Copy (*tracer_in[lev], Sborder_s, 0, 0, tracer_in[lev]->nComp(), nghost);
        }
 
        // We need this to avoid FPE
        m_u_mac[lev]->setVal(covered_val);
        m_v_mac[lev]->setVal(covered_val);
        m_w_mac[lev]->setVal(covered_val);
 
        // Predict normal velocity to faces -- note that the {u_mac, v_mac, w_mac}
        //    arrays returned from this call are on face CENTROIDS
        if (use_godunov) 
           incflo_predict_godunov(lev, time, vel_in, vel_forces_in);
        else
           incflo_predict_vels_on_faces(lev, time, vel_in);
    }

    // Do projection on all AMR levels in one shot -- note that the {u_mac, v_mac, w_mac}
    //    arrays returned from this call are in on face CENTROIDS
    apply_MAC_projection (m_u_mac, m_v_mac, m_w_mac, density_in, time);
}
