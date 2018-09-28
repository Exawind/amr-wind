#include <incflo_level.H>
#include <mac_F.H>
#include <convection_F.H>

//
// Compute acc using the vel passed in
//
void incflo_level::incflo_compute_ugradu_predictor(int lev,
												   MultiFab& conv,
												   Vector<std::unique_ptr<MultiFab>>& vel)
{
	BL_PROFILE("incflo_level::incflo_compute_ugradu");
	Box domain(geom[lev].Domain());

	incflo_compute_velocity_at_faces(lev, vel);

	mac_projection->apply_projection(m_u_mac, m_v_mac, m_w_mac, ro);

    // Get EB geometric info
    Array< const MultiCutFab*,AMREX_SPACEDIM> areafrac;
    Array< const MultiCutFab*,AMREX_SPACEDIM> facecent;
    const amrex::MultiFab*                    volfrac;
    const amrex::MultiCutFab*                 bndrycent;

    areafrac  =   ebfactory[lev] -> getAreaFrac();
    facecent  =   ebfactory[lev] -> getFaceCent();
    volfrac   = &(ebfactory[lev] -> getVolFrac());
    bndrycent = &(ebfactory[lev] -> getBndryCent());
       
#ifdef _OPENMP
#pragma omp parallel
#endif
	for(MFIter mfi(*vel[lev], true); mfi.isValid(); ++mfi)
	{
		// Tilebox
		Box bx = mfi.tilebox();

		// this is to check efficiently if this tile contains any eb stuff
		const EBFArrayBox& vel_fab = dynamic_cast<EBFArrayBox const&>((*vel[lev])[mfi]);
		const EBCellFlagFab& flags = vel_fab.getEBCellFlagFab();

		if(flags.getType(amrex::grow(bx, 0)) == FabType::covered)
		{
			// If tile is completely covered by EB geometry, set slopes
			// value to some very large number so we know if
			// we accidentaly use these covered slopes later in calculations
			conv.setVal(1.2345e300, bx, 0, 3);
		}
		else
		{
			// No cut cells in tile + nghost-cell witdh halo -> use non-eb routine
			if(flags.getType(amrex::grow(bx, nghost)) == FabType::regular)
			{
				compute_ugradu(BL_TO_FORTRAN_BOX(bx),
							   BL_TO_FORTRAN_ANYD(conv[mfi]),
							   BL_TO_FORTRAN_ANYD((*vel[lev])[mfi]),
							   BL_TO_FORTRAN_ANYD((*m_u_mac[lev])[mfi]),
							   BL_TO_FORTRAN_ANYD((*m_v_mac[lev])[mfi]),
							   BL_TO_FORTRAN_ANYD((*m_w_mac[lev])[mfi]),
							   (*xslopes[lev])[mfi].dataPtr(),
							   (*yslopes[lev])[mfi].dataPtr(),
							   BL_TO_FORTRAN_ANYD((*zslopes[lev])[mfi]),
							   domain.loVect(),
							   domain.hiVect(),
							   bc_ilo.dataPtr(),
							   bc_ihi.dataPtr(),
							   bc_jlo.dataPtr(),
							   bc_jhi.dataPtr(),
							   bc_klo.dataPtr(),
							   bc_khi.dataPtr(),
							   geom[lev].CellSize(),
							   &nghost,
							   &ugradu_type);
			}
			else
			{
				compute_ugradu_eb(BL_TO_FORTRAN_BOX(bx),
								  BL_TO_FORTRAN_ANYD(conv[mfi]),
								  BL_TO_FORTRAN_ANYD((*vel[lev])[mfi]),
								  BL_TO_FORTRAN_ANYD((*m_u_mac[lev])[mfi]),
								  BL_TO_FORTRAN_ANYD((*m_v_mac[lev])[mfi]),
								  BL_TO_FORTRAN_ANYD((*m_w_mac[lev])[mfi]),
								  BL_TO_FORTRAN_ANYD((*areafrac[0])[mfi]),
								  BL_TO_FORTRAN_ANYD((*areafrac[1])[mfi]),
								  BL_TO_FORTRAN_ANYD((*areafrac[2])[mfi]),
								  BL_TO_FORTRAN_ANYD((*facecent[0])[mfi]),
								  BL_TO_FORTRAN_ANYD((*facecent[1])[mfi]),
								  BL_TO_FORTRAN_ANYD((*facecent[2])[mfi]),
								  BL_TO_FORTRAN_ANYD(flags),
								  BL_TO_FORTRAN_ANYD((*volfrac)[mfi]),
                                  BL_TO_FORTRAN_ANYD((*bndrycent)[mfi]),
								  (*xslopes[lev])[mfi].dataPtr(),
								  (*yslopes[lev])[mfi].dataPtr(),
								  BL_TO_FORTRAN_ANYD((*zslopes[lev])[mfi]),
								  domain.loVect(),
								  domain.hiVect(),
								  bc_ilo.dataPtr(),
								  bc_ihi.dataPtr(),
								  bc_jlo.dataPtr(),
								  bc_jhi.dataPtr(),
								  bc_klo.dataPtr(),
								  bc_khi.dataPtr(),
								  geom[lev].CellSize(),
								  &nghost);
			}
		}
	}
}

//
// Compute acc using the vel passed in
//
void incflo_level::incflo_compute_ugradu_corrector(int lev,
												   MultiFab& conv,
												   Vector<std::unique_ptr<MultiFab>>& vel)
{
	BL_PROFILE("incflo_level::incflo_compute_ugradu");
	Box domain(geom[lev].Domain());

	// First compute the slopes
	incflo_compute_velocity_slopes(lev, vel);
	incflo_compute_velocity_at_faces(lev, vel);

	mac_projection->apply_projection(m_u_mac, m_v_mac, m_w_mac, ro);
	// Get EB geometric info
	Array<const MultiCutFab*, AMREX_SPACEDIM> areafrac;
	Array<const MultiCutFab*, AMREX_SPACEDIM> facecent;
	const amrex::MultiFab* volfrac;
    const amrex::MultiCutFab* bndrycent;

	areafrac = ebfactory[lev]->getAreaFrac();
	facecent = ebfactory[lev]->getFaceCent();
	volfrac = &(ebfactory[lev]->getVolFrac());
    bndrycent = &(ebfactory[lev]->getBndryCent());

#ifdef _OPENMP
#pragma omp parallel
#endif
	for(MFIter mfi(*vel[lev], true); mfi.isValid(); ++mfi)
	{
		// Tilebox
		Box bx = mfi.tilebox();

		// this is to check efficiently if this tile contains any eb stuff
		const EBFArrayBox& vel_fab = dynamic_cast<EBFArrayBox const&>((*vel[lev])[mfi]);
		const EBCellFlagFab& flags = vel_fab.getEBCellFlagFab();

		if(flags.getType(amrex::grow(bx, 0)) == FabType::covered)
		{
			// If tile is completely covered by EB geometry, set slopes
			// value to some very large number so we know if
			// we accidentaly use these covered slopes later in calculations
			conv.setVal(1.2345e300, bx, 0, 3);
		}
		else
		{
			// No cut cells in tile + nghost-cell witdh halo -> use non-eb routine
			if(flags.getType(amrex::grow(bx, nghost)) == FabType::regular)
			{
				compute_ugradu(BL_TO_FORTRAN_BOX(bx),
							   BL_TO_FORTRAN_ANYD(conv[mfi]),
							   BL_TO_FORTRAN_ANYD((*vel[lev])[mfi]),
							   BL_TO_FORTRAN_ANYD((*m_u_mac[lev])[mfi]),
							   BL_TO_FORTRAN_ANYD((*m_v_mac[lev])[mfi]),
							   BL_TO_FORTRAN_ANYD((*m_w_mac[lev])[mfi]),
							   (*xslopes[lev])[mfi].dataPtr(),
							   (*yslopes[lev])[mfi].dataPtr(),
							   BL_TO_FORTRAN_ANYD((*zslopes[lev])[mfi]),
							   domain.loVect(),
							   domain.hiVect(),
							   bc_ilo.dataPtr(),
							   bc_ihi.dataPtr(),
							   bc_jlo.dataPtr(),
							   bc_jhi.dataPtr(),
							   bc_klo.dataPtr(),
							   bc_khi.dataPtr(),
							   geom[lev].CellSize(),
							   &nghost,
							   &ugradu_type);
			}
			else
			{
				compute_ugradu_eb(BL_TO_FORTRAN_BOX(bx),
								  BL_TO_FORTRAN_ANYD(conv[mfi]),
								  BL_TO_FORTRAN_ANYD((*vel[lev])[mfi]),
								  BL_TO_FORTRAN_ANYD((*m_u_mac[lev])[mfi]),
								  BL_TO_FORTRAN_ANYD((*m_v_mac[lev])[mfi]),
								  BL_TO_FORTRAN_ANYD((*m_w_mac[lev])[mfi]),
								  BL_TO_FORTRAN_ANYD((*areafrac[0])[mfi]),
								  BL_TO_FORTRAN_ANYD((*areafrac[1])[mfi]),
								  BL_TO_FORTRAN_ANYD((*areafrac[2])[mfi]),
								  BL_TO_FORTRAN_ANYD((*facecent[0])[mfi]),
								  BL_TO_FORTRAN_ANYD((*facecent[1])[mfi]),
								  BL_TO_FORTRAN_ANYD((*facecent[2])[mfi]),
								  BL_TO_FORTRAN_ANYD(flags),
								  BL_TO_FORTRAN_ANYD((*volfrac)[mfi]),
                                  BL_TO_FORTRAN_ANYD((*bndrycent)[mfi]),
								  (*xslopes[lev])[mfi].dataPtr(),
								  (*yslopes[lev])[mfi].dataPtr(),
								  BL_TO_FORTRAN_ANYD((*zslopes[lev])[mfi]),
								  domain.loVect(),
								  domain.hiVect(),
								  bc_ilo.dataPtr(),
								  bc_ihi.dataPtr(),
								  bc_jlo.dataPtr(),
								  bc_jhi.dataPtr(),
								  bc_klo.dataPtr(),
								  bc_khi.dataPtr(),
								  geom[lev].CellSize(),
								  &nghost);
			}
		}
	}
}

//
// Compute the slopes of each velocity component in the
// three directions.
//
void incflo_level::incflo_compute_velocity_slopes(int lev, Vector<std::unique_ptr<MultiFab>>& vel)
{
	BL_PROFILE("incflo_level::incflo_compute_velocity_slopes");

	Box domain(geom[lev].Domain());

#ifdef _OPENMP
#pragma omp parallel
#endif
	for(MFIter mfi(*vel[lev], true); mfi.isValid(); ++mfi)
	{
		// Tilebox
		Box bx = mfi.tilebox();

		// this is to check efficiently if this tile contains any eb stuff
		const EBFArrayBox& vel_fab = dynamic_cast<EBFArrayBox const&>((*vel[lev])[mfi]);
		const EBCellFlagFab& flags = vel_fab.getEBCellFlagFab();

		if(flags.getType(amrex::grow(bx, 0)) == FabType::covered)
		{
			// If tile is completely covered by EB geometry, set slopes
			// value to some very large number so we know if
			// we accidentaly use these covered slopes later in calculations
			xslopes[lev]->setVal(1.2345e300, bx, 0, 3);
			yslopes[lev]->setVal(1.2345e300, bx, 0, 3);
			zslopes[lev]->setVal(1.2345e300, bx, 0, 3);
		}
		else
		{
			// No cut cells in tile + 1-cell witdh halo -> use non-eb routine
			if(flags.getType(amrex::grow(bx, 1)) == FabType::regular)
			{
				compute_slopes(BL_TO_FORTRAN_BOX(bx),
							   BL_TO_FORTRAN_ANYD((*vel[lev])[mfi]),
							   (*xslopes[lev])[mfi].dataPtr(),
							   (*yslopes[lev])[mfi].dataPtr(),
							   BL_TO_FORTRAN_ANYD((*zslopes[lev])[mfi]),
							   domain.loVect(),
							   domain.hiVect(),
							   bc_ilo.dataPtr(),
							   bc_ihi.dataPtr(),
							   bc_jlo.dataPtr(),
							   bc_jhi.dataPtr(),
							   bc_klo.dataPtr(),
							   bc_khi.dataPtr(),
							   &nghost);
			}
			else
			{
				compute_slopes_eb(BL_TO_FORTRAN_BOX(bx),
								  BL_TO_FORTRAN_ANYD((*vel[lev])[mfi]),
								  (*xslopes[lev])[mfi].dataPtr(),
								  (*yslopes[lev])[mfi].dataPtr(),
								  BL_TO_FORTRAN_ANYD((*zslopes[lev])[mfi]),
								  BL_TO_FORTRAN_ANYD(flags),
								  domain.loVect(),
								  domain.hiVect(),
								  bc_ilo.dataPtr(),
								  bc_ihi.dataPtr(),
								  bc_jlo.dataPtr(),
								  bc_jhi.dataPtr(),
								  bc_klo.dataPtr(),
								  bc_khi.dataPtr(),
								  &nghost);
			}
		}
	}

	xslopes[lev]->FillBoundary(geom[lev].periodicity());
	yslopes[lev]->FillBoundary(geom[lev].periodicity());
	zslopes[lev]->FillBoundary(geom[lev].periodicity());
}

void incflo_level::incflo_compute_velocity_at_faces(int lev, Vector<std::unique_ptr<MultiFab>>& vel)
{
	BL_PROFILE("incflo_level::incflo_compute_velocity_at_faces");
	Box domain(geom[lev].Domain());

	// First compute the slopes
	incflo_compute_velocity_slopes(lev, vel);

	// Get EB geometric info
	Array<const MultiCutFab*, AMREX_SPACEDIM> areafrac;
	Array<const MultiCutFab*, AMREX_SPACEDIM> facecent;

	areafrac = ebfactory[lev]->getAreaFrac();
	facecent = ebfactory[lev]->getFaceCent();

// Then compute velocity at faces
#ifdef _OPENMP
#pragma omp parallel
#endif
	for(MFIter mfi(*vel[lev], true); mfi.isValid(); ++mfi)
	{
		// Tilebox
		Box bx = mfi.tilebox();
		Box ubx = mfi.tilebox(e_x);
		Box vbx = mfi.tilebox(e_y);
		Box wbx = mfi.tilebox(e_z);

		// this is to check efficiently if this tile contains any eb stuff
		const EBFArrayBox& vel_fab = dynamic_cast<EBFArrayBox const&>((*vel[lev])[mfi]);
		const EBCellFlagFab& flags = vel_fab.getEBCellFlagFab();

		if(flags.getType(amrex::grow(bx, 0)) == FabType::covered)
		{
			m_u_mac[lev]->setVal(1.2345e300, ubx, 0, 1);
			m_v_mac[lev]->setVal(1.2345e300, vbx, 0, 1);
			m_w_mac[lev]->setVal(1.2345e300, wbx, 0, 1);
		}
		else
		{
			// No cut cells in tile + 1-cell witdh halo -> use non-eb routine
			if(flags.getType(amrex::grow(bx, 1)) == FabType::regular)
			{
				compute_velocity_at_faces(BL_TO_FORTRAN_BOX(bx),
										  BL_TO_FORTRAN_ANYD((*m_u_mac[lev])[mfi]),
										  BL_TO_FORTRAN_ANYD((*m_v_mac[lev])[mfi]),
										  BL_TO_FORTRAN_ANYD((*m_w_mac[lev])[mfi]),
										  BL_TO_FORTRAN_ANYD((*vel[lev])[mfi]),
										  BL_TO_FORTRAN_ANYD((*xslopes[lev])[mfi]),
										  (*yslopes[lev])[mfi].dataPtr(),
										  (*zslopes[lev])[mfi].dataPtr(),
										  bc_ilo.dataPtr(),
										  bc_ihi.dataPtr(),
										  bc_jlo.dataPtr(),
										  bc_jhi.dataPtr(),
										  bc_klo.dataPtr(),
										  bc_khi.dataPtr(),
										  &nghost,
										  domain.loVect(),
										  domain.hiVect());
			}
			else
			{
				compute_velocity_at_x_faces_eb(BL_TO_FORTRAN_BOX(ubx),
											   BL_TO_FORTRAN_ANYD((*m_u_mac[lev])[mfi]),
											   BL_TO_FORTRAN_ANYD((*vel[lev])[mfi]),
											   BL_TO_FORTRAN_ANYD((*xslopes[lev])[mfi]),
											   BL_TO_FORTRAN_ANYD((*areafrac[0])[mfi]),
											   BL_TO_FORTRAN_ANYD((*facecent[0])[mfi]),
											   BL_TO_FORTRAN_ANYD(flags),
											   bc_ilo.dataPtr(),
											   bc_ihi.dataPtr(),
											   &nghost,
											   domain.loVect(),
											   domain.hiVect());

				compute_velocity_at_y_faces_eb(BL_TO_FORTRAN_BOX(vbx),
											   BL_TO_FORTRAN_ANYD((*m_v_mac[lev])[mfi]),
											   BL_TO_FORTRAN_ANYD((*vel[lev])[mfi]),
											   BL_TO_FORTRAN_ANYD((*yslopes[lev])[mfi]),
											   BL_TO_FORTRAN_ANYD((*areafrac[1])[mfi]),
											   BL_TO_FORTRAN_ANYD((*facecent[1])[mfi]),
											   BL_TO_FORTRAN_ANYD(flags),
											   bc_jlo.dataPtr(),
											   bc_jhi.dataPtr(),
											   &nghost,
											   domain.loVect(),
											   domain.hiVect());

				compute_velocity_at_z_faces_eb(BL_TO_FORTRAN_BOX(wbx),
											   BL_TO_FORTRAN_ANYD((*m_w_mac[lev])[mfi]),
											   BL_TO_FORTRAN_ANYD((*vel[lev])[mfi]),
											   BL_TO_FORTRAN_ANYD((*zslopes[lev])[mfi]),
											   BL_TO_FORTRAN_ANYD((*areafrac[2])[mfi]),
											   BL_TO_FORTRAN_ANYD((*facecent[2])[mfi]),
											   BL_TO_FORTRAN_ANYD(flags),
											   bc_klo.dataPtr(),
											   bc_khi.dataPtr(),
											   &nghost,
											   domain.loVect(),
											   domain.hiVect());
			}
		}
	}
}
