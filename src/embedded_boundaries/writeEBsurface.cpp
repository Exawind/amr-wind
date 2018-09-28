#include <algorithm>
#include <eb_F.H>
#include <incflo_level.H>

void incflo_level::WriteEBSurface(int lev)
{
	if(Geom(0).isAllPeriodic())
		return;

	const Real* dx = Geom(lev).CellSize();

	BoxArray ba = grids[lev];

	// This creates the associated Distribution Mapping
	// DistributionMapping dm(ba, ParallelDescriptor::NProcs());

	MultiFab mf_ba(ba, dmap[lev], 1, 0, MFInfo(), *ebfactory[lev]);

	// // // Deliberately didn't time this loop.
	for(MFIter mfi(mf_ba); mfi.isValid(); ++mfi)
	{

		const auto& sfab = dynamic_cast<EBFArrayBox const&>((mf_ba)[mfi]);
		const auto& my_flag = sfab.getEBCellFlagFab();

		const Box& bx = mfi.validbox();

		if(my_flag.getType(bx) == FabType::covered or my_flag.getType(bx) == FabType::regular)
			continue;

		std::array<const MultiCutFab*, AMREX_SPACEDIM> areafrac;
		const MultiCutFab* bndrycent;

		areafrac = ebfactory[lev]->getAreaFrac();
		bndrycent = &(ebfactory[lev]->getBndryCent());

		incflo_eb_to_polygon(dx,
							 bx.loVect(),
							 bx.hiVect(),
							 my_flag.dataPtr(),
							 my_flag.loVect(),
							 my_flag.hiVect(),
							 (*bndrycent)[mfi].dataPtr(),
							 (*bndrycent)[mfi].loVect(),
							 (*bndrycent)[mfi].hiVect(),
							 (*areafrac[0])[mfi].dataPtr(),
							 (*areafrac[0])[mfi].loVect(),
							 (*areafrac[0])[mfi].hiVect(),
							 (*areafrac[1])[mfi].dataPtr(),
							 (*areafrac[1])[mfi].loVect(),
							 (*areafrac[1])[mfi].hiVect(),
							 (*areafrac[2])[mfi].dataPtr(),
							 (*areafrac[2])[mfi].loVect(),
							 (*areafrac[2])[mfi].hiVect());
	}

	int cpu = ParallelDescriptor::MyProc();
	int nProcs = ParallelDescriptor::NProcs();

	incflo_write_eb_vtp(&cpu);

	if(ParallelDescriptor::IOProcessor())
		incflo_write_pvtp(&nProcs);

	// // // Deliberately didn't time this loop.
	for(MFIter mfi(mf_ba); mfi.isValid(); ++mfi)
	{

		const auto& sfab = dynamic_cast<EBFArrayBox const&>((mf_ba)[mfi]);
		const auto& my_flag = sfab.getEBCellFlagFab();

		const Box& bx = mfi.validbox();

		if(my_flag.getType(bx) == FabType::covered or my_flag.getType(bx) == FabType::regular)
			continue;

		incflo_eb_grid_coverage(&cpu,
								dx,
								bx.loVect(),
								bx.hiVect(),
								my_flag.dataPtr(),
								my_flag.loVect(),
								my_flag.hiVect());
	}
}
