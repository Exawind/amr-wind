#include <incflo_level.H>

void incflo_level::AllocateArrays(int lev)
{
	// ********************************************************************************
	// Cell- or node-based arrays
	// ********************************************************************************

	// Gas density
	ro_g[lev].reset(new MultiFab(grids[lev], dmap[lev], 1, nghost));
	ro_go[lev].reset(new MultiFab(grids[lev], dmap[lev], 1, nghost));
	ro_g[lev]->setVal(0.);
	ro_go[lev]->setVal(0.);

	if(nodal_pressure)
	{
		const BoxArray& nd_grids = amrex::convert(grids[lev], IntVect{1, 1, 1});

		p0_g[lev].reset(new MultiFab(nd_grids, dmap[lev], 1, 0));
		p_g[lev].reset(new MultiFab(nd_grids, dmap[lev], 1, 0));
		p_go[lev].reset(new MultiFab(nd_grids, dmap[lev], 1, 0));
		pp_g[lev].reset(new MultiFab(nd_grids, dmap[lev], 1, 0));
	}
	else
	{

		p0_g[lev].reset(new MultiFab(grids[lev], dmap[lev], 1, nghost));
		p_g[lev].reset(new MultiFab(grids[lev], dmap[lev], 1, nghost));
		p_go[lev].reset(new MultiFab(grids[lev], dmap[lev], 1, nghost));
		pp_g[lev].reset(new MultiFab(grids[lev], dmap[lev], 1, nghost));
	}

	p0_g[lev]->setVal(0.);
	p_g[lev]->setVal(0.);
	p_go[lev]->setVal(0.);
	pp_g[lev]->setVal(0.);

	// Presssure gradients
	gp[lev].reset(new MultiFab(grids[lev], dmap[lev], 3, nghost));
	gp0[lev].reset(new MultiFab(grids[lev], dmap[lev], 3, nghost));
	gp[lev]->setVal(0.);
	gp0[lev]->setVal(0.);

	// Molecular viscosity
	mu_g[lev].reset(new MultiFab(grids[lev], dmap[lev], 1, nghost));
	mu_g[lev]->setVal(0.);

	// Coefficient of grad(div(u)) in viscous terms
	lambda_g[lev].reset(new MultiFab(grids[lev], dmap[lev], 1, nghost));
	lambda_g[lev]->setVal(0.);

	// Current velocity
	vel_g[lev].reset(new MultiFab(grids[lev], dmap[lev], 3, nghost, MFInfo(), *ebfactory[lev]));
	vel_g[lev]->setVal(0.);

	// Old velocity
	vel_go[lev].reset(new MultiFab(grids[lev], dmap[lev], 3, nghost, MFInfo(), *ebfactory[lev]));
	vel_go[lev]->setVal(0.);
}

void incflo_level::AllocateTempArrays(int lev)
{
	// Div(u)
	trD_g[lev].reset(new MultiFab(grids[lev], dmap[lev], 1, nghost));
	trD_g[lev]->setVal(0.);

	// Vorticity
	vort[lev].reset(new MultiFab(grids[lev], dmap[lev], 1, nghost));
	vort[lev]->setVal(0.);

	if(nodal_pressure)
	{
		const BoxArray& nd_grids = amrex::convert(grids[lev], IntVect{1, 1, 1});

		phi[lev].reset(new MultiFab(nd_grids, dmap[lev], 1, 0));
		diveu[lev].reset(new MultiFab(nd_grids, dmap[lev], 1, 0));
	}
	else
	{

		phi[lev].reset(new MultiFab(grids[lev], dmap[lev], 1, nghost, MFInfo(), *ebfactory[lev]));
		diveu[lev].reset(new MultiFab(grids[lev], dmap[lev], 1, nghost, MFInfo(), *ebfactory[lev]));
	}

	phi[lev]->setVal(0.);
	diveu[lev]->setVal(0.);

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

	// When the pressure is on nodes, bcoeff is at cell centers
	if(nodal_pressure)
	{
		bcoeff[lev][0].reset(new MultiFab(grids[lev], dmap[lev], 1, nghost));
		bcoeff[lev][1].reset(new MultiFab(grids[lev], dmap[lev], 1, nghost));
		bcoeff[lev][2].reset(new MultiFab(grids[lev], dmap[lev], 1, nghost));
	}
	else
	{

		// Create a BoxArray on x-faces.
		BoxArray x_edge_ba = grids[lev];
		x_edge_ba.surroundingNodes(0);
		bcoeff[lev]
			  [0].reset(new MultiFab(x_edge_ba, dmap[lev], 1, nghost, MFInfo(), *ebfactory[lev]));

		// Create a BoxArray on y-faces.
		BoxArray y_edge_ba = grids[lev];
		y_edge_ba.surroundingNodes(1);
		bcoeff[lev]
			  [1].reset(new MultiFab(y_edge_ba, dmap[lev], 1, nghost, MFInfo(), *ebfactory[lev]));

		// Create a BoxArray on y-faces.
		BoxArray z_edge_ba = grids[lev];
		z_edge_ba.surroundingNodes(2);
		bcoeff[lev]
			  [2].reset(new MultiFab(z_edge_ba, dmap[lev], 1, nghost, MFInfo(), *ebfactory[lev]));
	}

	bcoeff[lev][0]->setVal(0.);
	bcoeff[lev][1]->setVal(0.);
	bcoeff[lev][2]->setVal(0.);

	// ****************************************************************

	// Create a BoxArray on x-faces.
	BoxArray x_edge_ba = grids[lev];
	x_edge_ba.surroundingNodes(0);
	bcoeff_diff[lev][0].reset(new MultiFab(x_edge_ba, dmap[lev], 1, nghost));
	m_u_mac[lev].reset(new MultiFab(x_edge_ba, dmap[lev], 1, nghost, MFInfo(), *ebfactory[lev]));

	// Create a BoxArray on y-faces.
	BoxArray y_edge_ba = grids[lev];
	y_edge_ba.surroundingNodes(1);
	bcoeff_diff[lev][1].reset(new MultiFab(y_edge_ba, dmap[lev], 1, nghost));
	m_v_mac[lev].reset(new MultiFab(y_edge_ba, dmap[lev], 1, nghost, MFInfo(), *ebfactory[lev]));

	// Create a BoxArray on y-faces.
	BoxArray z_edge_ba = grids[lev];
	z_edge_ba.surroundingNodes(2);
	bcoeff_diff[lev][2].reset(new MultiFab(z_edge_ba, dmap[lev], 1, nghost));
	m_w_mac[lev].reset(new MultiFab(z_edge_ba, dmap[lev], 1, nghost, MFInfo(), *ebfactory[lev]));

	bcoeff_diff[lev][0]->setVal(0.);
	bcoeff_diff[lev][1]->setVal(0.);
	bcoeff_diff[lev][2]->setVal(0.);

	m_u_mac[lev]->setVal(0.);
	m_v_mac[lev]->setVal(0.);
	m_w_mac[lev]->setVal(0.);
}

void incflo_level::RegridArrays(int lev, BoxArray& new_grids, DistributionMapping& new_dmap)
{
	// ********************************************************************************
	// Cell-based arrays
	// ********************************************************************************

	// Gas density
	int ng = ro_g[lev]->nGrow();
	std::unique_ptr<MultiFab> ro_g_new(new MultiFab(new_grids, new_dmap, 1, ro_g[lev]->nGrow()));
	ro_g_new->copy(*ro_g[lev], 0, 0, 1, ng, ng);
	ro_g_new->FillBoundary(geom[lev].periodicity());
	ro_g[lev] = std::move(ro_g_new);

	// Old gas density
	ng = ro_go[lev]->nGrow();
	std::unique_ptr<MultiFab> ro_go_new(new MultiFab(new_grids, new_dmap, 1, ro_go[lev]->nGrow()));
	ro_go_new->copy(*ro_go[lev], 0, 0, 1, ng, ng);
	ro_go_new->FillBoundary(geom[lev].periodicity());
	ro_go[lev] = std::move(ro_go_new);

	if(nodal_pressure)
	{
		const BoxArray& nd_grids = amrex::convert(new_grids, IntVect{1, 1, 1});

		ng = p_g[lev]->nGrow();
		std::unique_ptr<MultiFab> p_g_new(new MultiFab(nd_grids, new_dmap, 1, p_g[lev]->nGrow()));
		p_g_new->copy(*p_g[lev], 0, 0, 1, ng, ng);
		p_g_new->FillBoundary(geom[lev].periodicity());
		p_g[lev] = std::move(p_g_new);

		ng = p_go[lev]->nGrow();
		std::unique_ptr<MultiFab> p_go_new(new MultiFab(nd_grids, new_dmap, 1, p_go[lev]->nGrow()));
		p_go_new->copy(*p_go[lev], 0, 0, 1, ng, ng);
		p_go_new->FillBoundary(geom[lev].periodicity());
		p_go[lev] = std::move(p_go_new);

		ng = p0_g[lev]->nGrow();
		std::unique_ptr<MultiFab> p0_g_new(new MultiFab(nd_grids, new_dmap, 1, ng));
		p0_g_new->copy(*p0_g[lev], 0, 0, 1, ng, ng);
		p0_g_new->FillBoundary(p0_periodicity);
		p0_g[lev] = std::move(p0_g_new);

		ng = pp_g[lev]->nGrow();
		std::unique_ptr<MultiFab> pp_g_new(new MultiFab(nd_grids, new_dmap, 1, pp_g[lev]->nGrow()));
		pp_g_new->copy(*pp_g[lev], 0, 0, 1, ng, ng);
		pp_g_new->FillBoundary(geom[lev].periodicity());
		pp_g[lev] = std::move(pp_g_new);

		std::unique_ptr<MultiFab> diveu_new(
			new MultiFab(nd_grids, new_dmap, 1, diveu[lev]->nGrow()));
		diveu[lev] = std::move(diveu_new);
		diveu[lev]->setVal(0.);

		std::unique_ptr<MultiFab> phi_new(new MultiFab(nd_grids, new_dmap, 1, phi[lev]->nGrow()));
		phi[lev] = std::move(phi_new);
		phi[lev]->setVal(0.);

		std::unique_ptr<MultiFab> bc0_new(
			new MultiFab(new_grids, new_dmap, 1, bcoeff[lev][0]->nGrow()));
		bcoeff[lev][0] = std::move(bc0_new);
		bcoeff[lev][0]->setVal(0.);

		std::unique_ptr<MultiFab> bc1_new(
			new MultiFab(new_grids, new_dmap, 1, bcoeff[lev][1]->nGrow()));
		bcoeff[lev][1] = std::move(bc1_new);
		bcoeff[lev][1]->setVal(0.);

		std::unique_ptr<MultiFab> bc2_new(
			new MultiFab(new_grids, new_dmap, 1, bcoeff[lev][2]->nGrow()));
		bcoeff[lev][2] = std::move(bc2_new);
		bcoeff[lev][2]->setVal(0.);
	}
	else
	{

		ng = p_g[lev]->nGrow();
		std::unique_ptr<MultiFab> p_g_new(new MultiFab(new_grids, new_dmap, 1, p_g[lev]->nGrow()));
		p_g_new->copy(*p_g[lev], 0, 0, 1, ng, ng);
		p_g_new->FillBoundary(geom[lev].periodicity());
		p_g[lev] = std::move(p_g_new);

		ng = p_go[lev]->nGrow();
		std::unique_ptr<MultiFab> p_go_new(
			new MultiFab(new_grids, new_dmap, 1, p_go[lev]->nGrow()));
		p_go_new->copy(*p_go[lev], 0, 0, 1, ng, ng);
		p_go_new->FillBoundary(geom[lev].periodicity());
		p_go[lev] = std::move(p_go_new);

		ng = p0_g[lev]->nGrow();
		std::unique_ptr<MultiFab> p0_g_new(
			new MultiFab(new_grids, new_dmap, 1, p0_g[lev]->nGrow()));
		p0_g_new->copy(*p0_g[lev], 0, 0, 1, ng, ng);
		p0_g_new->FillBoundary(p0_periodicity);
		p0_g[lev] = std::move(p0_g_new);

		ng = pp_g[lev]->nGrow();
		std::unique_ptr<MultiFab> pp_g_new(
			new MultiFab(new_grids, new_dmap, 1, pp_g[lev]->nGrow()));
		pp_g_new->copy(*pp_g[lev], 0, 0, 1, ng, ng);
		pp_g_new->FillBoundary(geom[lev].periodicity());
		pp_g[lev] = std::move(pp_g_new);

		std::unique_ptr<MultiFab> phi_new(new MultiFab(new_grids, new_dmap, 1, phi[lev]->nGrow()));
		phi[lev] = std::move(phi_new);
		phi[lev]->setVal(0.);

		std::unique_ptr<MultiFab> diveu_new(
			new MultiFab(new_grids, new_dmap, 1, diveu[lev]->nGrow()));
		diveu[lev] = std::move(diveu_new);
		diveu[lev]->setVal(0.);

		// Cell-centered pressure uses face-based coefficients
		BoxArray x_ba = new_grids;
		x_ba = x_ba.surroundingNodes(0);
		std::unique_ptr<MultiFab> bc0_new(new MultiFab(x_ba, new_dmap, 1, nghost));
		bcoeff[lev][0] = std::move(bc0_new);
		bcoeff[lev][0]->setVal(0.0);

		BoxArray y_ba = new_grids;
		y_ba = y_ba.surroundingNodes(1);
		std::unique_ptr<MultiFab> bc1_new(new MultiFab(y_ba, new_dmap, 1, nghost));
		bcoeff[lev][1] = std::move(bc1_new);
		bcoeff[lev][1]->setVal(0.0);

		BoxArray z_ba = new_grids;
		z_ba = z_ba.surroundingNodes(2);
		std::unique_ptr<MultiFab> bc2_new(new MultiFab(z_ba, new_dmap, 1, nghost));
		bcoeff[lev][2] = std::move(bc2_new);
		bcoeff[lev][2]->setVal(0.0);
	}

	// Molecular viscosity
	ng = mu_g[lev]->nGrow();
	std::unique_ptr<MultiFab> mu_g_new(new MultiFab(new_grids, new_dmap, 1, mu_g[lev]->nGrow()));
	mu_g_new->copy(*mu_g[lev], 0, 0, 1, ng, ng);
	mu_g_new->FillBoundary(geom[lev].periodicity());
	mu_g[lev] = std::move(mu_g_new);

	// Lambda
	ng = lambda_g[lev]->nGrow();
	std::unique_ptr<MultiFab> lambda_g_new(
		new MultiFab(new_grids, new_dmap, 1, lambda_g[lev]->nGrow()));
	lambda_g_new->copy(*lambda_g[lev], 0, 0, 1, ng, ng);
	lambda_g_new->FillBoundary(geom[lev].periodicity());
	lambda_g[lev] = std::move(lambda_g_new);

	// Gas velocity
	ng = vel_g[lev]->nGrow();
	std::unique_ptr<MultiFab> vel_g_new(new MultiFab(new_grids, new_dmap, vel_g[lev]->nComp(), ng));
	vel_g_new->copy(*vel_g[lev], 0, 0, vel_g[lev]->nComp(), ng, ng);
	vel_g_new->FillBoundary(geom[lev].periodicity());
	vel_g[lev] = std::move(vel_g_new);

	// Old gas velocity
	ng = vel_go[lev]->nGrow();
	std::unique_ptr<MultiFab> vel_go_new(
		new MultiFab(new_grids, new_dmap, vel_go[lev]->nComp(), ng));
	vel_go_new->copy(*vel_go[lev], 0, 0, vel_go[lev]->nComp(), ng, ng);
	vel_go_new->FillBoundary(geom[lev].periodicity());
	vel_go[lev] = std::move(vel_go_new);

	// Pressure gradients
	ng = gp[lev]->nGrow();
	std::unique_ptr<MultiFab> gp_new(new MultiFab(new_grids, new_dmap, 1, gp[lev]->nGrow()));
	gp_new->copy(*gp[lev], 0, 0, 1, ng, ng);
	gp_new->FillBoundary(geom[lev].periodicity());
	gp[lev] = std::move(gp_new);

	// Pressure gradients
	ng = gp0[lev]->nGrow();
	std::unique_ptr<MultiFab> gp0_new(new MultiFab(new_grids, new_dmap, 1, gp0[lev]->nGrow()));
	gp0_new->copy(*gp0[lev], 0, 0, 1, ng, ng);
	gp0_new->FillBoundary(geom[lev].periodicity());
	gp0[lev] = std::move(gp0_new);

	// Trace(D)
	ng = trD_g[lev]->nGrow();
	std::unique_ptr<MultiFab> trD_g_new(new MultiFab(new_grids, new_dmap, 1, trD_g[lev]->nGrow()));
	trD_g[lev] = std::move(trD_g_new);
	trD_g[lev]->setVal(0.);

	// Vorticity
	ng = vort[lev]->nGrow();
	std::unique_ptr<MultiFab> vort_new(new MultiFab(new_grids, new_dmap, 1, vort[lev]->nGrow()));
	vort[lev] = std::move(vort_new);
	vort[lev]->setVal(0.);

	/****************************************************************************
    * Face-based Arrays                                                        *
    ****************************************************************************/

	BoxArray x_ba = new_grids;
	x_ba = x_ba.surroundingNodes(0);

	// MAC velocity
	std::unique_ptr<MultiFab> u_mac_new(new MultiFab(x_ba, new_dmap, 1, nghost));
	m_u_mac[lev] = std::move(u_mac_new);
	m_u_mac[lev]->setVal(0.0);

	// Diffusion coefficient on x-faces
	std::unique_ptr<MultiFab> bc0_new(new MultiFab(x_ba, new_dmap, 1, nghost));
	bcoeff_diff[lev][0] = std::move(bc0_new);
	bcoeff_diff[lev][0]->setVal(0.0);

	//****************************************************************************

	BoxArray y_ba = new_grids;
	y_ba = y_ba.surroundingNodes(1);

	// MAC velocity
	std::unique_ptr<MultiFab> v_mac_new(new MultiFab(y_ba, new_dmap, 1, nghost));
	m_v_mac[lev] = std::move(v_mac_new);
	m_v_mac[lev]->setVal(0.0);

	// Diffusion coefficient on y-faces
	std::unique_ptr<MultiFab> bc1_new(new MultiFab(y_ba, new_dmap, 1, nghost));
	bcoeff_diff[lev][1] = std::move(bc1_new);
	bcoeff_diff[lev][1]->setVal(0.0);

	//****************************************************************************

	BoxArray z_ba = new_grids;
	z_ba = z_ba.surroundingNodes(2);

	// MAC velocity
	std::unique_ptr<MultiFab> w_mac_new(new MultiFab(z_ba, new_dmap, 1, nghost));
	m_w_mac[lev] = std::move(w_mac_new);
	m_w_mac[lev]->setVal(0.0);

	// Diffusion coefficient on z-faces
	std::unique_ptr<MultiFab> bc2_new(new MultiFab(z_ba, new_dmap, 1, nghost));
	bcoeff[lev][2] = std::move(bc2_new);
	bcoeff[lev][2]->setVal(0.0);

	/****************************************************************************
    * Nodal Arrays                                                             *
    ****************************************************************************/

	// ********************************************************************************
	// Make sure we fill the ghost cells as appropriate -- this is copied from init_fluid
	// ********************************************************************************

	fill_mf_bc(lev, *ro_g[lev]);
	fill_mf_bc(lev, *ro_go[lev]);

	fill_mf_bc(lev, *mu_g[lev]);
	fill_mf_bc(lev, *lambda_g[lev]);

	if(!nodal_pressure)
	{
		fill_mf_bc(lev, *p_g[lev]);
		fill_mf_bc(lev, *p_go[lev]);
	}
}
