#include <incflo.H>

using namespace amrex;

void incflo::AllocateArrays(int lev)
{
#ifdef AMREX_USE_EB
    UpdateEBFactory(lev);
#endif

#ifdef AMREX_USE_EB
    const FabFactory<FArrayBox>& factory =  *ebfactory[lev];
#else
    const FabFactory<FArrayBox>& factory = FArrayBoxFactory();
#endif

    // ********************************************************************************
    // Cell-based arrays
    // ********************************************************************************

    // Level Mask = 1 if not covered by a finer patch, else 0 
    // Note we never pass an EBFactory to this one
    level_mask[lev].reset(new iMultiFab(grids[lev], dmap[lev], 1, 0, MFInfo() /*, default factory*/));
    level_mask[lev]->setVal(1);
    
    // Current Density
    density[lev].reset(new MultiFab(grids[lev], dmap[lev], 1, nghost, MFInfo(), factory));
    density[lev]->setVal(0.);

    // Old density
    density_o[lev].reset(new MultiFab(grids[lev], dmap[lev], 1, nghost, MFInfo(), factory));
    density_o[lev]->setVal(0.);

    // Current Tracer; default to 0
    tracer[lev].reset(new MultiFab(grids[lev], dmap[lev], ntrac, nghost, MFInfo(), factory));
    tracer[lev]->setVal(0.);

    // Old tracer; default to 0
    tracer_o[lev].reset(new MultiFab(grids[lev], dmap[lev], ntrac, nghost, MFInfo(), factory));
    tracer_o[lev]->setVal(0.);

    // Current Velocity
    vel[lev].reset(new MultiFab(grids[lev], dmap[lev], AMREX_SPACEDIM, nghost, MFInfo(), factory));
    vel[lev]->setVal(0.);

    // Old velocity
    vel_o[lev].reset(new MultiFab(grids[lev], dmap[lev], AMREX_SPACEDIM, nghost, MFInfo(), factory));
    vel_o[lev]->setVal(0.);

    // Pressure gradient
    gp[lev].reset(new MultiFab(grids[lev], dmap[lev], AMREX_SPACEDIM, nghost, MFInfo(), factory));
    gp[lev]->setVal(0.);

    // Current Viscosity
    eta[lev].reset(new MultiFab(grids[lev], dmap[lev], 1, nghost, MFInfo(), factory));
    eta[lev]->setVal(0.);

    // Old Viscosity
    eta_old[lev].reset(new MultiFab(grids[lev], dmap[lev], 1, nghost, MFInfo(), factory));
    eta_old[lev]->setVal(0.);

    // Strain-rate magnitude
    strainrate[lev].reset(new MultiFab(grids[lev], dmap[lev], 1, nghost, MFInfo(), factory));
    strainrate[lev]->setVal(0.);

    // Vorticity
    vort[lev].reset(new MultiFab(grids[lev], dmap[lev], 1, nghost, MFInfo(), factory));
    vort[lev]->setVal(0.);

    // Drag
    drag[lev].reset(new MultiFab(grids[lev], dmap[lev], 1, nghost, MFInfo(), factory));
    drag[lev]->setVal(0.);

    // Convective terms for velocity
    conv_u[lev].reset(new MultiFab(grids[lev], dmap[lev], AMREX_SPACEDIM, 0, MFInfo(), factory));
    conv_u[lev]->setVal(0.);

    // Convective terms for density 
    conv_r[lev].reset(new MultiFab(grids[lev], dmap[lev], 1, 0, MFInfo(), factory));
    conv_r[lev]->setVal(0.);

    // Convective terms for tracers
    conv_t[lev].reset(new MultiFab(grids[lev], dmap[lev], ntrac, 0, MFInfo(), factory));
    conv_t[lev]->setVal(0.);

    // Old Convective terms for velocity
    conv_u_old[lev].reset(new MultiFab(grids[lev], dmap[lev], AMREX_SPACEDIM, 0, MFInfo(), factory));
    conv_u_old[lev]->setVal(0.);

    // Convective terms for density
    conv_r_old[lev].reset(new MultiFab(grids[lev], dmap[lev], 1, 0, MFInfo(), factory));
    conv_r_old[lev]->setVal(0.);

    // Convective terms for tracers
    conv_t_old[lev].reset(new MultiFab(grids[lev], dmap[lev], ntrac, 0, MFInfo(), factory));
    conv_t_old[lev]->setVal(0.);

    // Divergence of stress tensor terms
    divtau[lev].reset(new MultiFab(grids[lev], dmap[lev], AMREX_SPACEDIM, 0, MFInfo(), factory));
    divtau_old[lev].reset(new MultiFab(grids[lev], dmap[lev], AMREX_SPACEDIM, 0, MFInfo(), factory));
    divtau[lev]->setVal(0.);
    divtau_old[lev]->setVal(0.);

    // Forcing terms on velocity
    vel_forces[lev].reset(new MultiFab(grids[lev], dmap[lev], AMREX_SPACEDIM, 2, MFInfo(), factory));
    vel_forces[lev]->setVal(0.);

    // Forcing terms on density (last component and always zero) and tracers
    scal_forces[lev].reset(new MultiFab(grids[lev], dmap[lev], ntrac+1, 2, MFInfo(), factory));
    scal_forces[lev]->setVal(0.);

    // Scalar diffusion terms
    laps[lev].reset(new MultiFab(grids[lev], dmap[lev], ntrac, 0, MFInfo(), factory));
    laps_old[lev].reset(new MultiFab(grids[lev], dmap[lev], ntrac, 0, MFInfo(), factory));
    laps[lev]->setVal(0.);
    laps_old[lev]->setVal(0.);

    // Slopes in x-direction
    xslopes_u[lev].reset(new MultiFab(grids[lev], dmap[lev], AMREX_SPACEDIM, nghost, MFInfo(), factory));
    xslopes_u[lev]->setVal(0.);
    xslopes_r[lev].reset(new MultiFab(grids[lev], dmap[lev], 1, nghost, MFInfo(), factory));
    xslopes_r[lev]->setVal(0.);
    xslopes_t[lev].reset(new MultiFab(grids[lev], dmap[lev], ntrac, nghost, MFInfo(), factory));
    xslopes_t[lev]->setVal(0.);

    // Slopes in y-direction
    yslopes_u[lev].reset(new MultiFab(grids[lev], dmap[lev], AMREX_SPACEDIM, nghost, MFInfo(), factory));
    yslopes_u[lev]->setVal(0.);
    yslopes_r[lev].reset(new MultiFab(grids[lev], dmap[lev], 1, nghost, MFInfo(), factory));
    yslopes_r[lev]->setVal(0.);
    yslopes_t[lev].reset(new MultiFab(grids[lev], dmap[lev], ntrac, nghost, MFInfo(), factory));
    yslopes_t[lev]->setVal(0.);

    // Slopes in z-direction
    zslopes_u[lev].reset(new MultiFab(grids[lev], dmap[lev], AMREX_SPACEDIM, nghost, MFInfo(), factory));
    zslopes_u[lev]->setVal(0.);
    zslopes_r[lev].reset(new MultiFab(grids[lev], dmap[lev], ntrac, nghost, MFInfo(), factory));
    zslopes_r[lev]->setVal(0.);
    zslopes_t[lev].reset(new MultiFab(grids[lev], dmap[lev], ntrac, nghost, MFInfo(), factory));
    zslopes_t[lev]->setVal(0.);

    // ********************************************************************************
    // Node-based arrays
    // ********************************************************************************

    const BoxArray & nd_grids = amrex::convert(grids[lev], IntVect{1,1,1});

    // Pressure
    p0[lev].reset(new MultiFab(nd_grids, dmap[lev], 1, nghost, MFInfo(), factory));
    p0[lev]->setVal(0.);
    p[lev].reset(new MultiFab(nd_grids, dmap[lev], 1, nghost, MFInfo(), factory));
    p[lev]->setVal(0.);

    // Divergence of velocity field
    divu[lev].reset(new MultiFab(nd_grids, dmap[lev], 1, nghost, MFInfo(), factory));
    divu[lev]->setVal(0.);

    // ********************************************************************************
    // Face-based arrays
    // ********************************************************************************

    // Create a BoxArray on x-faces.
    BoxArray x_edge_ba = grids[lev];
    x_edge_ba.surroundingNodes(0);
    m_u_mac[lev].reset(new MultiFab(x_edge_ba, dmap[lev], 1, nghost, MFInfo(), factory));
    m_u_mac[lev]->setVal(0.);

    // Create a BoxArray on y-faces.
    BoxArray y_edge_ba = grids[lev];
    y_edge_ba.surroundingNodes(1);
    m_v_mac[lev].reset(new MultiFab(y_edge_ba, dmap[lev], 1, nghost, MFInfo(), factory));
    m_v_mac[lev]->setVal(0.);

    // Create a BoxArray on y-faces.
    BoxArray z_edge_ba = grids[lev];
    z_edge_ba.surroundingNodes(2);
    m_w_mac[lev].reset(new MultiFab(z_edge_ba, dmap[lev], 1, nghost, MFInfo(), factory));
    m_w_mac[lev]->setVal(0.);
}

void incflo::RegridArrays(int lev)
{
#ifdef AMREX_USE_EB
    bool need_regrid = UpdateEBFactory(lev);
#else
    bool need_regrid = true;
#endif

#ifdef AMREX_USE_EB
    const FabFactory<FArrayBox>& factory =  *ebfactory[lev];
#else
    const FabFactory<FArrayBox>& factory = FArrayBoxFactory();
#endif

   // exit this function is ebfactory has not been updated because that means
   // that dm and ba haven't changed
   if (!need_regrid)
       return;

   // ********************************************************************************
   // Cell-based arrays
   // ********************************************************************************
   //
   // After calling copy() with dst_ngrow set to ng, we do not need to call
   // FillBoundary().
   //

   // Level Mask = 1 if not covered by a finer patch, else 0
   level_mask[lev].reset(new iMultiFab(grids[lev], dmap[lev], 1, 0, MFInfo() /*, default factory*/));
   level_mask[lev]->setVal(1);
    
   // Density
   std::unique_ptr<MultiFab> density_new(new MultiFab(grids[lev], dmap[lev], 1, nghost, MFInfo(), factory));
   density_new->setVal(0.0);
   density_new->copy(*density[lev], 0, 0, 1, 0, nghost);
   density[lev] = std::move(density_new);

   // Old Density
   std::unique_ptr<MultiFab> density_o_new(new MultiFab(grids[lev], dmap[lev], 1, nghost, MFInfo(), factory));
   density_o_new->setVal(0.0);
   density_o_new->copy(*density_o[lev], 0, 0, 1, 0, nghost);
   density_o[lev] = std::move(density_o_new);

   // Tracer
   std::unique_ptr<MultiFab> tracer_new(new MultiFab(grids[lev], dmap[lev], ntrac, nghost, MFInfo(), factory));
   tracer_new->setVal(0.0);
   tracer_new->copy(*tracer[lev], 0, 0, ntrac, 0, nghost);
   tracer[lev] = std::move(tracer_new);

   // Old Tracer
   std::unique_ptr<MultiFab> tracer_o_new(new MultiFab(grids[lev], dmap[lev], ntrac, nghost, MFInfo(), factory));
   tracer_o_new->setVal(0.0);
   tracer_o_new->copy(*tracer_o[lev], 0, 0, ntrac, 0, nghost);
   tracer_o[lev] = std::move(tracer_o_new);

   // Gas velocity
   std::unique_ptr<MultiFab> vel_new(new MultiFab(grids[lev], dmap[lev], AMREX_SPACEDIM, nghost, MFInfo(), factory));
   vel_new->setVal(0.);
   vel_new->copy(*vel[lev], 0, 0, vel[lev]->nComp(), 0, nghost);
   vel[lev] = std::move(vel_new);

    // Old gas velocity
    std::unique_ptr<MultiFab> vel_o_new(new MultiFab(grids[lev], dmap[lev], AMREX_SPACEDIM, nghost, MFInfo(), factory));
    vel_o_new->setVal(0.);
    vel_o_new->copy(*vel_o[lev], 0, 0, vel_o[lev]->nComp(), 0, nghost);
    vel_o[lev] = std::move(vel_o_new);

    // Pressure gradients
    std::unique_ptr<MultiFab> gp_new(new MultiFab(grids[lev], dmap[lev], AMREX_SPACEDIM, nghost, MFInfo(), factory));
    gp_new->setVal(0.);
    gp_new->copy(*gp[lev], 0, 0, gp[lev]->nComp(), 0, nghost);
    gp[lev] = std::move(gp_new);

    // Apparent viscosity
    std::unique_ptr<MultiFab> eta_new(new MultiFab(grids[lev], dmap[lev], 1, nghost, MFInfo(), factory));
    eta_new->setVal(0.);
    eta_new->copy(*eta[lev], 0, 0, 1, 0, nghost);
    eta[lev] = std::move(eta_new);

    std::unique_ptr<MultiFab> eta_old_new(new MultiFab(grids[lev], dmap[lev], 1, nghost, MFInfo(), factory));
    eta_old_new->setVal(0.);
    eta_old_new->copy(*eta_old[lev], 0, 0, 1, 0, nghost);
    eta_old[lev] = std::move(eta_old_new);

    // ***************************
    // Strain-rate magnitude
    // ***************************
    std::unique_ptr<MultiFab> strainrate_new(new MultiFab(grids[lev], dmap[lev], 1, nghost, MFInfo(), factory));
    strainrate[lev] = std::move(strainrate_new);
    strainrate[lev]->setVal(0.);

    // ***************************
    // Vorticity
    // ***************************
    std::unique_ptr<MultiFab> vort_new(new MultiFab(grids[lev], dmap[lev], 1, nghost, MFInfo(), factory));
    vort[lev] = std::move(vort_new);
    vort[lev]->setVal(0.);

    // ***************************
    // Drag
    // ***************************
    std::unique_ptr<MultiFab> drag_new(new MultiFab(grids[lev], dmap[lev], 1, nghost, MFInfo(), factory));
    drag[lev] = std::move(drag_new);
    drag[lev]->setVal(0.);

    // ***************************
    // Convective terms for velocity
    // ***************************
    std::unique_ptr<MultiFab> conv_u_new(new MultiFab(grids[lev], dmap[lev], conv_u[lev]->nComp(), nghost, MFInfo(), factory));
    conv_u[lev] = std::move(conv_u_new);
    conv_u[lev]->setVal(0.);

    std::unique_ptr<MultiFab> conv_u_old_new(new MultiFab(grids[lev], dmap[lev], conv_u_old[lev]->nComp(), nghost, MFInfo(), factory));
    conv_u_old[lev] = std::move(conv_u_old_new);
    conv_u_old[lev]->setVal(0.);

    // ***************************
    // Convective terms for density
    // ***************************
    std::unique_ptr<MultiFab> conv_r_new(new MultiFab(grids[lev], dmap[lev], conv_r[lev]->nComp(), nghost, MFInfo(), factory));
    conv_r[lev] = std::move(conv_r_new);
    conv_r[lev]->setVal(0.);

    std::unique_ptr<MultiFab> conv_r_old_new(new MultiFab(grids[lev], dmap[lev], conv_r_old[lev]->nComp(), nghost, MFInfo(), factory));
    conv_r_old[lev] = std::move(conv_r_old_new);
    conv_r_old[lev]->setVal(0.);

    // ***************************
    // Convective terms for tracer
    // ***************************
    std::unique_ptr<MultiFab> conv_t_new(new MultiFab(grids[lev], dmap[lev], conv_t[lev]->nComp(), nghost, MFInfo(), factory));
    conv_t[lev] = std::move(conv_t_new);
    conv_t[lev]->setVal(0.);

    std::unique_ptr<MultiFab> conv_t_old_new(new MultiFab(grids[lev], dmap[lev], conv_t_old[lev]->nComp(), nghost, MFInfo(), factory));
    conv_t_old[lev] = std::move(conv_t_old_new);
    conv_t_old[lev]->setVal(0.);

    // ***************************
    // Divergence of stress tensor terms 
    // ***************************
    std::unique_ptr<MultiFab> divtau_new(new MultiFab(grids[lev], dmap[lev], AMREX_SPACEDIM, nghost, MFInfo(), factory));
    divtau[lev] = std::move(divtau_new);
    divtau[lev]->setVal(0.);

    std::unique_ptr<MultiFab> divtau_old_new(new MultiFab(grids[lev], dmap[lev], AMREX_SPACEDIM, nghost, MFInfo(), factory));
    divtau_old[lev] = std::move(divtau_old_new);
    divtau_old[lev]->setVal(0.);

    // ***************************
    // Velocity forcing terms 
    // ***************************
    std::unique_ptr<MultiFab> vel_forces_new(new MultiFab(grids[lev], dmap[lev], AMREX_SPACEDIM, 2, MFInfo(), factory));
    vel_forces[lev] = std::move(vel_forces_new);
    vel_forces[lev]->setVal(0.);

    // ***************************
    // Scalar forcing terms 
    // ***************************
    std::unique_ptr<MultiFab> scal_forces_new(new MultiFab(grids[lev], dmap[lev], scal_forces[lev]->nComp(), 2, MFInfo(), factory));
    scal_forces[lev] = std::move(scal_forces_new);
    scal_forces[lev]->setVal(0.);

    // ***************************
    // Scalar diffusion terms 
    // ***************************
    std::unique_ptr<MultiFab> laps_new(new MultiFab(grids[lev], dmap[lev], ntrac, nghost, MFInfo(), factory));
    laps[lev] = std::move(laps_new);
    laps[lev]->setVal(0.);

    std::unique_ptr<MultiFab> laps_old_new(new MultiFab(grids[lev], dmap[lev], ntrac, nghost, MFInfo(), factory));
    laps_old[lev] = std::move(laps_old_new);
    laps_old[lev]->setVal(0.);

    // ***************************
    // Slopes in x-direction
    // ***************************
    std::unique_ptr<MultiFab> xslopes_u_new(new MultiFab(grids[lev], dmap[lev], AMREX_SPACEDIM, nghost, MFInfo(), factory));
    xslopes_u[lev] = std::move(xslopes_u_new);
    xslopes_u[lev] -> setVal(0.);

    std::unique_ptr<MultiFab> xslopes_r_new(new MultiFab(grids[lev], dmap[lev], 1, nghost, MFInfo(), factory));
    xslopes_r[lev] = std::move(xslopes_r_new);
    xslopes_r[lev] -> setVal(0.);

    std::unique_ptr<MultiFab> xslopes_t_new(new MultiFab(grids[lev], dmap[lev], ntrac, nghost, MFInfo(), factory));
    xslopes_t[lev] = std::move(xslopes_t_new);
    xslopes_t[lev] -> setVal(0.);

    // ***************************
    // Slopes in y-direction
    // ***************************
    std::unique_ptr<MultiFab> yslopes_u_new(new MultiFab(grids[lev], dmap[lev], AMREX_SPACEDIM, nghost, MFInfo(), factory));
    yslopes_u[lev] = std::move(yslopes_u_new);
    yslopes_u[lev] -> setVal(0.);

    std::unique_ptr<MultiFab> yslopes_r_new(new MultiFab(grids[lev], dmap[lev], 1, nghost, MFInfo(), factory));
    yslopes_r[lev] = std::move(yslopes_r_new);
    yslopes_r[lev] -> setVal(0.);

    std::unique_ptr<MultiFab> yslopes_t_new(new MultiFab(grids[lev], dmap[lev], ntrac, nghost, MFInfo(), factory));
    yslopes_t[lev] = std::move(yslopes_t_new);
    yslopes_t[lev] -> setVal(0.);

    // ***************************
    // Slopes in z-direction
    // ***************************
    std::unique_ptr<MultiFab> zslopes_u_new(new MultiFab(grids[lev], dmap[lev], AMREX_SPACEDIM, nghost, MFInfo(), factory));
    zslopes_u[lev] = std::move(zslopes_u_new);
    zslopes_u[lev] -> setVal(0.);

    std::unique_ptr<MultiFab> zslopes_r_new(new MultiFab(grids[lev], dmap[lev], 1, nghost, MFInfo(), factory));
    zslopes_r[lev] = std::move(zslopes_r_new);
    zslopes_r[lev] -> setVal(0.);

    std::unique_ptr<MultiFab> zslopes_t_new(new MultiFab(grids[lev], dmap[lev], ntrac, nghost, MFInfo(), factory));
    zslopes_t[lev] = std::move(zslopes_t_new);
    zslopes_t[lev] -> setVal(0.);

    /****************************************************************************
    * Node-based Arrays                                                        *
    ****************************************************************************/

    // Pressures, projection vars
    const BoxArray & nd_grids = amrex::convert(grids[lev], IntVect{1,1,1});

    std::unique_ptr<MultiFab> p_new(new MultiFab(nd_grids, dmap[lev], 1, nghost, MFInfo(), factory));
    p_new->setVal(0.0);
    p_new->copy(*p[lev],0,0,1,0,nghost);
    p[lev] = std::move(p_new);

    std::unique_ptr<MultiFab> p0_new(new MultiFab(nd_grids, dmap[lev], 1, nghost, MFInfo(), factory));
    p0_new->setVal(0.0);
    p0_new->copy(*p0[lev],0,0,1,0,nghost);
    p0[lev] = std::move(p0_new);

    std::unique_ptr<MultiFab> divu_new(new MultiFab(nd_grids, dmap[lev], 1, nghost, MFInfo(), factory));
    divu[lev] = std::move(divu_new);
    divu[lev]->setVal(0.);

    /****************************************************************************
    * Face-based Arrays                                                        *
    ****************************************************************************/

    BoxArray x_ba = grids[lev];
    x_ba = x_ba.surroundingNodes(0);

    // MAC velocity
    std::unique_ptr<MultiFab> u_mac_new(new MultiFab(x_ba, dmap[lev], 1, nghost, MFInfo(), factory));
    m_u_mac[lev] = std::move(u_mac_new);
    m_u_mac[lev]->setVal(0.0);

    BoxArray y_ba = grids[lev];
    y_ba = y_ba.surroundingNodes(1);

    // MAC velocity
    std::unique_ptr<MultiFab> v_mac_new(new MultiFab(y_ba, dmap[lev], 1, nghost, MFInfo(), factory));
    m_v_mac[lev] = std::move(v_mac_new);
    m_v_mac[lev] -> setVal(0.0);

    BoxArray z_ba = grids[lev];
    z_ba = z_ba.surroundingNodes(2);

    // MAC velocity
    std::unique_ptr<MultiFab> w_mac_new(new MultiFab(z_ba, dmap[lev], 1, nghost, MFInfo(), factory));
    m_w_mac[lev] = std::move(w_mac_new);
    m_w_mac[lev] -> setVal(0.0);
}

// Resize all arrays when instance of incflo class is constructed.
// This is only done at the very start of the simulation. 
void incflo::ResizeArrays()
{
    // Time holders for fillpatch stuff
    t_new.resize(max_level + 1);
    t_old.resize(max_level + 1);

    // Density 
    density.resize(max_level + 1);
    density_o.resize(max_level + 1);

    // Tracer 
    tracer.resize(max_level + 1);
    tracer_o.resize(max_level + 1);

    // Current (vel) and old (vel_o) velocities
    vel.resize(max_level + 1);
    vel_o.resize(max_level + 1);

    // Pressure
    p.resize(max_level + 1);
    p0.resize(max_level + 1);

    // Pressure gradients
    gp.resize(max_level + 1);

    // Derived quantities: viscosity, strainrate, vorticity, div(u)
    eta.resize(max_level + 1);
    eta_old.resize(max_level + 1);
    strainrate.resize(max_level + 1);
    vort.resize(max_level + 1);
    drag.resize(max_level + 1);
    divu.resize(max_level + 1);

    // Convective terms u grad u 
    conv_u.resize(max_level + 1);
    conv_u_old.resize(max_level + 1);
    conv_r.resize(max_level + 1);
    conv_r_old.resize(max_level + 1);
    conv_t.resize(max_level + 1);
    conv_t_old.resize(max_level + 1);

    divtau.resize(max_level + 1);
    divtau_old.resize(max_level + 1);

    vel_forces.resize(max_level + 1);
    scal_forces.resize(max_level + 1);

    laps.resize(max_level + 1);
    laps_old.resize(max_level + 1);

    // MAC velocities used for defining convective term
    m_u_mac.resize(max_level + 1);
    m_v_mac.resize(max_level + 1);
    m_w_mac.resize(max_level + 1);

    // Slopes used for upwinding convective terms
    xslopes_u.resize(max_level + 1);
    yslopes_u.resize(max_level + 1);
    zslopes_u.resize(max_level + 1);
    xslopes_r.resize(max_level + 1);
    yslopes_r.resize(max_level + 1);
    zslopes_r.resize(max_level + 1);
    xslopes_t.resize(max_level + 1);
    yslopes_t.resize(max_level + 1);
    zslopes_t.resize(max_level + 1);

    // BCs
    bc_ilo.resize(max_level + 1);
    bc_ihi.resize(max_level + 1);
    bc_jlo.resize(max_level + 1);
    bc_jhi.resize(max_level + 1);
    bc_klo.resize(max_level + 1);
    bc_khi.resize(max_level + 1);

#ifdef AMREX_USE_EB
    // EB factory
    ebfactory.resize(max_level + 1);
#endif
    
    level_mask.resize(max_level + 1);
}

void incflo::MakeBCArrays()
{
    for(int lev = 0; lev <= max_level; lev++)
    {
        // Define and allocate the integer MultiFab that is the outside adjacent cells of the
        // problem domain.
        Box domainx(geom[lev].Domain());
        domainx.grow(1, nghost);
        domainx.grow(2, nghost);
        Box box_ilo = amrex::adjCellLo(domainx, 0, 1);
        Box box_ihi = amrex::adjCellHi(domainx, 0, 1);

        Box domainy(geom[lev].Domain());
        domainy.grow(0, nghost);
        domainy.grow(2, nghost);
        Box box_jlo = amrex::adjCellLo(domainy, 1, 1);
        Box box_jhi = amrex::adjCellHi(domainy, 1, 1);

        Box domainz(geom[lev].Domain());
        domainz.grow(0, nghost);
        domainz.grow(1, nghost);
        Box box_klo = amrex::adjCellLo(domainz, 2, 1);
        Box box_khi = amrex::adjCellHi(domainz, 2, 1);

        // Note that each of these is a single IArrayBox so every process has a copy of them
        bc_ilo[lev].reset(new IArrayBox(box_ilo, 2));
        bc_ihi[lev].reset(new IArrayBox(box_ihi, 2));
        bc_jlo[lev].reset(new IArrayBox(box_jlo, 2));
        bc_jhi[lev].reset(new IArrayBox(box_jhi, 2));
        bc_klo[lev].reset(new IArrayBox(box_klo, 2));
        bc_khi[lev].reset(new IArrayBox(box_khi, 2));
    }
}

