#include <incflo.H>

void incflo::AllocateArrays(int lev)
{
    incflo_update_ebfactory(lev);

	// ********************************************************************************
	// Cell- or node-based arrays
	// ********************************************************************************

    // Gas density
	ro[lev].reset(new MultiFab(grids[lev], dmap[lev], 1, nghost));
	ro_o[lev].reset(new MultiFab(grids[lev], dmap[lev], 1, nghost));
	ro[lev]->setVal(0.);
	ro_o[lev]->setVal(0.);

    p0[lev].reset(new MultiFab(grids[lev], dmap[lev], 1, nghost));
    p[lev].reset(new MultiFab(grids[lev], dmap[lev], 1, nghost));
    p_o[lev].reset(new MultiFab(grids[lev], dmap[lev], 1, nghost));
    pp[lev].reset(new MultiFab(grids[lev], dmap[lev], 1, nghost));

	p0[lev]->setVal(0.);
	p[lev]->setVal(0.);
	p_o[lev]->setVal(0.);
	pp[lev]->setVal(0.);

	// Presssure gradients
	gp[lev].reset(new MultiFab(grids[lev], dmap[lev], 3, nghost));
	gp0[lev].reset(new MultiFab(grids[lev], dmap[lev], 3, nghost));
	gp[lev]->setVal(0.);
	gp0[lev]->setVal(0.);

	// Molecular viscosity
	eta[lev].reset(new MultiFab(grids[lev], dmap[lev], 1, nghost));
	eta[lev]->setVal(0.);

	// Current velocity
	vel[lev].reset(new MultiFab(grids[lev], dmap[lev], 3, nghost, MFInfo(), *ebfactory[lev]));
	vel[lev]->setVal(0.);

	// Old velocity
	vel_o[lev].reset(new MultiFab(grids[lev], dmap[lev], 3, nghost, MFInfo(), *ebfactory[lev]));
	vel_o[lev]->setVal(0.);
    
	// Strain-rate magnitude
	strainrate[lev].reset(new MultiFab(grids[lev], dmap[lev], 1, nghost));
	strainrate[lev]->setVal(0.);

	// Vorticity
	vort[lev].reset(new MultiFab(grids[lev], dmap[lev], 1, nghost));
	vort[lev]->setVal(0.);

    phi[lev].reset(new MultiFab(grids[lev], dmap[lev], 1, nghost, MFInfo(), *ebfactory[lev]));
    divu[lev].reset(new MultiFab(grids[lev], dmap[lev], 1, nghost, MFInfo(), *ebfactory[lev]));

	phi[lev]->setVal(0.);
	divu[lev]->setVal(0.);

	// Arrays to store the solution and rhs for the diffusion solve
	phi_diff[lev].reset(new MultiFab(grids[lev], dmap[lev], 1, nghost));
	rhs_diff[lev].reset(new MultiFab(grids[lev], dmap[lev], 1, nghost));

	// Slopes in x-direction
	xslopes[lev].reset(new MultiFab(grids[lev], dmap[lev], 3, nghost));
	xslopes[lev]->setVal(0.);

	// Slopes in y-direction
	yslopes[lev].reset(new MultiFab(grids[lev], dmap[lev], 3, nghost));
	yslopes[lev]->setVal(0.);

	// Slopes in z-direction
	zslopes[lev].reset(new MultiFab(grids[lev], dmap[lev], 3, nghost));
	zslopes[lev]->setVal(0.);

	// ********************************************************************************
	// X-face-based arrays
	// ********************************************************************************

    // Create a BoxArray on x-faces.
    BoxArray x_edge_ba = grids[lev];
    x_edge_ba.surroundingNodes(0);
    bcoeff[lev][0].reset(new MultiFab(x_edge_ba, dmap[lev], 1, nghost, 
                         MFInfo(), *ebfactory[lev]));

    // Create a BoxArray on y-faces.
    BoxArray y_edge_ba = grids[lev];
    y_edge_ba.surroundingNodes(1);
    bcoeff[lev][1].reset(new MultiFab(y_edge_ba, dmap[lev], 1, nghost, 
                         MFInfo(), *ebfactory[lev]));

    // Create a BoxArray on y-faces.
    BoxArray z_edge_ba = grids[lev];
    z_edge_ba.surroundingNodes(2);
    bcoeff[lev][2].reset(new MultiFab(z_edge_ba, dmap[lev], 1, nghost, 
                         MFInfo(), *ebfactory[lev]));

	bcoeff[lev][0]->setVal(0.);
	bcoeff[lev][1]->setVal(0.);
	bcoeff[lev][2]->setVal(0.);

	// ****************************************************************

	// Create a BoxArray on x-faces.
	bcoeff_diff[lev][0].reset(new MultiFab(x_edge_ba, dmap[lev], 1, nghost));
	m_u_mac[lev].reset(new MultiFab(x_edge_ba, dmap[lev], 1, nghost, MFInfo(), *ebfactory[lev]));

	// Create a BoxArray on y-faces.
	bcoeff_diff[lev][1].reset(new MultiFab(y_edge_ba, dmap[lev], 1, nghost));
	m_v_mac[lev].reset(new MultiFab(y_edge_ba, dmap[lev], 1, nghost, MFInfo(), *ebfactory[lev]));

	// Create a BoxArray on y-faces.
	bcoeff_diff[lev][2].reset(new MultiFab(z_edge_ba, dmap[lev], 1, nghost));
	m_w_mac[lev].reset(new MultiFab(z_edge_ba, dmap[lev], 1, nghost, MFInfo(), *ebfactory[lev]));

	bcoeff_diff[lev][0]->setVal(0.);
	bcoeff_diff[lev][1]->setVal(0.);
	bcoeff_diff[lev][2]->setVal(0.);

	m_u_mac[lev]->setVal(0.);
	m_v_mac[lev]->setVal(0.);
	m_w_mac[lev]->setVal(0.);
}

void incflo::RegridArrays(int lev)
{
    incflo_update_ebfactory(lev);

	// ********************************************************************************
	// Cell-based arrays
	// ********************************************************************************
    //
    // After calling copy() with dst_ngrow set to ng, we do not need to call
    // FillBoundary().
    // 
    //

	// Gas density
	int ng = ro[lev]->nGrow();
	std::unique_ptr<MultiFab> ro_new(new MultiFab(grids[lev], dmap[lev], 1, ng));
	ro_new->setVal(0.0);
	ro_new->copy(*ro[lev], 0, 0, 1, 0, ng);
	ro[lev] = std::move(ro_new);

	// Old gas density
	ng = ro_o[lev]->nGrow();
	std::unique_ptr<MultiFab> ro_o_new(new MultiFab(grids[lev], dmap[lev], 1, ng));
	ro_o_new->setVal(0.0);
	ro_o_new->copy(*ro_o[lev], 0, 0, 1, ng, ng);
	ro_o[lev] = std::move(ro_o_new);

    ng = p[lev]->nGrow();
    std::unique_ptr<MultiFab> p_new(new MultiFab(grids[lev], dmap[lev], 1, ng));
    p_new->setVal(0.);
    p_new->copy(*p[lev], 0, 0, 1, 0, ng);
    p[lev] = std::move(p_new);

    ng = p_o[lev]->nGrow();
    std::unique_ptr<MultiFab> p_o_new(new MultiFab(grids[lev], dmap[lev], 1, ng));
    p_o_new->setVal(0.);
    p_o_new->copy(*p_o[lev], 0, 0, 1, 0, ng);
    p_o[lev] = std::move(p_o_new);

    ng = p0[lev]->nGrow();
    std::unique_ptr<MultiFab> p0_new(new MultiFab(grids[lev], dmap[lev], 1, ng));
    p0_new->setVal(0.);
    p0_new->copy(*p0[lev], 0, 0, 1, 0, ng);
    p0[lev] = std::move(p0_new);

    ng = pp[lev]->nGrow();
    std::unique_ptr<MultiFab> pp_new(new MultiFab(grids[lev], dmap[lev], 1, ng));
    pp_new->setVal(0.);
    pp_new->copy(*pp[lev], 0, 0, 1, 0, ng);
    pp[lev] = std::move(pp_new);

    ng = phi[lev]->nGrow();
    std::unique_ptr<MultiFab> phi_new(new MultiFab(grids[lev], dmap[lev], 1, ng, MFInfo(), *ebfactory[lev]));
    phi[lev] = std::move(phi_new);
    phi[lev]->setVal(0.);

    ng = divu[lev]->nGrow();
    std::unique_ptr<MultiFab> divu_new(new MultiFab(grids[lev], dmap[lev], 1, ng, MFInfo(), *ebfactory[lev]));
    divu[lev] = std::move(divu_new);
    divu[lev]->setVal(0.);

    // Cell-centered pressure uses face-based coefficients
    BoxArray x_ba = grids[lev];
    x_ba = x_ba.surroundingNodes(0);
    std::unique_ptr<MultiFab> bc0_new(new MultiFab(x_ba, dmap[lev], 1, nghost, MFInfo(), *ebfactory[lev]));
    bcoeff[lev][0] = std::move(bc0_new);
    bcoeff[lev][0]->setVal(0.0);

    BoxArray y_ba = grids[lev];
    y_ba = y_ba.surroundingNodes(1);
    std::unique_ptr<MultiFab> bc1_new(new MultiFab(y_ba, dmap[lev], 1, nghost, MFInfo(), *ebfactory[lev]));
    bcoeff[lev][1] = std::move(bc1_new);
    bcoeff[lev][1]->setVal(0.0);

    BoxArray z_ba = grids[lev];
    z_ba = z_ba.surroundingNodes(2);
    std::unique_ptr<MultiFab> bc2_new(new MultiFab(z_ba, dmap[lev], 1, nghost, MFInfo(), *ebfactory[lev]));
    bcoeff[lev][2] = std::move(bc2_new);
    bcoeff[lev][2]->setVal(0.0);

	// Molecular viscosity
	ng = eta[lev]->nGrow();
	std::unique_ptr<MultiFab> eta_new(new MultiFab(grids[lev], dmap[lev], 1, ng));
	eta_new->setVal(0.);
	eta_new->copy(*eta[lev], 0, 0, 1, 0, ng);
	eta[lev] = std::move(eta_new);

	// Gas velocity
	ng = vel[lev]->nGrow();
	std::unique_ptr<MultiFab> vel_new(new MultiFab(grids[lev], dmap[lev], vel[lev]->nComp(), ng, MFInfo(), *ebfactory[lev]));
	vel_new->setVal(0.);
	vel_new->copy(*vel[lev], 0, 0, vel[lev]->nComp(), 0, ng);
	vel[lev] = std::move(vel_new);

	// Old gas velocity
	ng = vel_o[lev]->nGrow();
    std::unique_ptr<MultiFab> vel_o_new(new MultiFab(grids[lev], dmap[lev], 
                                        vel_o[lev]->nComp(), ng, MFInfo(), *ebfactory[lev]));
	vel_o_new->setVal(0.);
	vel_o_new->copy(*vel_o[lev], 0, 0, vel_o[lev]->nComp(), 0, ng);
	vel_o[lev] = std::move(vel_o_new);

	// Pressure gradients
	ng = gp[lev]->nGrow();
	std::unique_ptr<MultiFab> gp_new(new MultiFab(grids[lev], dmap[lev], 3, ng));
    gp_new->setVal(0.);
	gp_new->copy(*gp[lev], 0, 0, 1, 0, ng);
	gp[lev] = std::move(gp_new);

	// Pressure gradients
	ng = gp0[lev]->nGrow();
	std::unique_ptr<MultiFab> gp0_new(new MultiFab(grids[lev], dmap[lev], 3, ng));
    gp0_new->setVal(0.);
	gp0_new->copy(*gp0[lev], 0, 0, 1, 0, ng);
	gp0[lev] = std::move(gp0_new);

	// Strain-rate magnitude
	ng = strainrate[lev]->nGrow();
	std::unique_ptr<MultiFab> strainrate_new(new MultiFab(grids[lev], dmap[lev], 1, ng));
	strainrate[lev] = std::move(strainrate_new);
	strainrate[lev]->setVal(0.);

	// Vorticity
	ng = vort[lev]->nGrow();
	std::unique_ptr<MultiFab> vort_new(new MultiFab(grids[lev], dmap[lev], 1, ng));
	vort[lev] = std::move(vort_new);
	vort[lev]->setVal(0.);

	/****************************************************************************
    * Face-based Arrays                                                        *
    ****************************************************************************/

	// MAC velocity
	std::unique_ptr<MultiFab> u_mac_new(new MultiFab(x_ba, dmap[lev], 1, nghost, MFInfo(), *ebfactory[lev]));
	m_u_mac[lev] = std::move(u_mac_new);
	m_u_mac[lev]->setVal(0.0);

	// Diffusion coefficient on x-faces
	bcoeff_diff[lev][0] = std::move(bc0_new);
	bcoeff_diff[lev][0]->setVal(0.0);

    // MAC velocity
    std::unique_ptr<MultiFab> v_mac_new(new MultiFab(y_ba, dmap[lev], 1, nghost, MFInfo(), *ebfactory[lev]));
    m_v_mac[lev] = std::move(v_mac_new);
    m_v_mac[lev] -> setVal(0.0);

    // Diffusion coefficient on y-faces
    bcoeff_diff[lev][1] = std::move(bc1_new);
    bcoeff_diff[lev][1] -> setVal(0.0);

    // MAC velocity
    std::unique_ptr<MultiFab> w_mac_new(new MultiFab(z_ba, dmap[lev], 1, nghost, MFInfo(), *ebfactory[lev]));
    m_w_mac[lev] = std::move(w_mac_new);
    m_w_mac[lev] -> setVal(0.0);

    // Diffusion coefficient on z-faces
    bcoeff[lev][2] = std::move(bc2_new);
    bcoeff[lev][2] -> setVal(0.0);

	//****************************************************************************

    // Slopes in x-direction
    ng = xslopes[lev] -> nGrow();
    std::unique_ptr<MultiFab> xslopes_new(new MultiFab(grids[lev], dmap[lev], xslopes[lev]->nComp(), nghost));
    xslopes[lev] = std::move(xslopes_new);
    xslopes[lev] -> setVal(0.);

    // Slopes in y-direction
    ng = yslopes[lev] -> nGrow();
    std::unique_ptr<MultiFab> yslopes_new(new MultiFab(grids[lev], dmap[lev], yslopes[lev]->nComp(), nghost));
    yslopes[lev] = std::move(yslopes_new);
    yslopes[lev] -> setVal(0.);

    // Slopes in z-direction
    ng = zslopes[lev] -> nGrow();
    std::unique_ptr<MultiFab> zslopes_new(new MultiFab(grids[lev], dmap[lev], zslopes[lev]->nComp(), nghost));
    zslopes[lev] = std::move(zslopes_new);
    zslopes[lev] -> setVal(0.);
    
	/****************************************************************************
    * Nodal Arrays                                                             *
    ****************************************************************************/

	// ********************************************************************************
	// Make sure we fill the ghost cells as appropriate -- this is copied from init_fluid
	// ********************************************************************************

	fill_mf_bc(lev, *ro[lev]);
	fill_mf_bc(lev, *ro_o[lev]);

	fill_mf_bc(lev, *eta[lev]);

    fill_mf_bc(lev, *p[lev]);
    fill_mf_bc(lev, *p_o[lev]);
}

void incflo::ResizeArrays()
{
	int nlevs_max = maxLevel() + 1;

	p.resize(nlevs_max);
	p_o.resize(nlevs_max);

	p0.resize(nlevs_max);
	pp.resize(nlevs_max);

	ro.resize(nlevs_max);
	ro_o.resize(nlevs_max);

	phi.resize(nlevs_max);
	divu.resize(nlevs_max);

	// RHS and solution arrays for diffusive solve
	rhs_diff.resize(nlevs_max);
	phi_diff.resize(nlevs_max);

	// Current (vel) and old (vel_o) velocities
	vel.resize(nlevs_max);
	vel_o.resize(nlevs_max);

	// Pressure gradients
	gp.resize(nlevs_max);
	gp0.resize(nlevs_max);

	eta.resize(nlevs_max);

    strainrate.resize(nlevs_max);
	vort.resize(nlevs_max);

	// MAC velocities used for defining convective term
	m_u_mac.resize(nlevs_max);
	m_v_mac.resize(nlevs_max);
	m_w_mac.resize(nlevs_max);

	xslopes.resize(nlevs_max);
	yslopes.resize(nlevs_max);
	zslopes.resize(nlevs_max);

	bcoeff.resize(nlevs_max);
	for(int i = 0; i < nlevs_max; ++i)
	{
		bcoeff[i].resize(3);
	}

	bcoeff_diff.resize(nlevs_max);
	for(int i = 0; i < nlevs_max; ++i)
	{
		bcoeff_diff[i].resize(3);
	}

	fluid_cost.resize(nlevs_max);

	// EB factory
	ebfactory.resize(nlevs_max);
}

