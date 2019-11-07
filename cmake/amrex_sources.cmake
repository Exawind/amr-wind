function(get_amrex_sources)
   #set(AMREX_SOURCE_DIR "${CMAKE_SOURCE_DIR}/submods/amrex/Src/Amr")
   #add_sources(GlobalSourceList
   #  ${AMREX_SOURCE_DIR}/AMReX_Amr.cpp
   #  ${AMREX_SOURCE_DIR}/AMReX_AmrLevel.cpp
   #  ${AMREX_SOURCE_DIR}/AMReX_Derive.cpp
   #  ${AMREX_SOURCE_DIR}/AMReX_StateData.cpp
   #  ${AMREX_SOURCE_DIR}/AMReX_StateDescriptor.cpp
   #  ${AMREX_SOURCE_DIR}/AMReX_AuxBoundaryData.cpp
   #  ${AMREX_SOURCE_DIR}/AMReX_Extrapolater.cpp
   #  ${AMREX_SOURCE_DIR}/AMReX_extrapolater_${AMREX_DIM}d.f90
   #)
   set(AMREX_SOURCE_DIR "${CMAKE_SOURCE_DIR}/submods/amrex/Src/AmrCore")
   add_sources(GlobalSourceList
     ${AMREX_SOURCE_DIR}/AMReX_AmrCore.cpp
     ${AMREX_SOURCE_DIR}/AMReX_Cluster.cpp
     ${AMREX_SOURCE_DIR}/AMReX_ErrorList.cpp 
     ${AMREX_SOURCE_DIR}/AMReX_FillPatchUtil.cpp
     ${AMREX_SOURCE_DIR}/AMReX_FluxRegister.cpp
     ${AMREX_SOURCE_DIR}/AMReX_Interpolater.cpp
     ${AMREX_SOURCE_DIR}/AMReX_TagBox.cpp
     ${AMREX_SOURCE_DIR}/AMReX_AmrMesh.cpp
     ${AMREX_SOURCE_DIR}/AMReX_FLUXREG_nd.F90
     ${AMREX_SOURCE_DIR}/AMReX_INTERP_${AMREX_DIM}D.F90
     ${AMREX_SOURCE_DIR}/AMReX_FillPatchUtil_${AMREX_DIM}d.F90
   )
   set(AMREX_SOURCE_DIR "${CMAKE_SOURCE_DIR}/submods/amrex/Src/Base")
   add_sources(GlobalSourceList
     ${AMREX_SOURCE_DIR}/AMReX.cpp
     ${AMREX_SOURCE_DIR}/AMReX_error_fi.cpp
     ${AMREX_SOURCE_DIR}/AMReX.cpp
     ${AMREX_SOURCE_DIR}/AMReX_error_fi.cpp
     ${AMREX_SOURCE_DIR}/AMReX_ParmParse.cpp
     ${AMREX_SOURCE_DIR}/AMReX_parmparse_fi.cpp
     ${AMREX_SOURCE_DIR}/AMReX_Utility.cpp
     ${AMREX_SOURCE_DIR}/AMReX_Random.cpp
     ${AMREX_SOURCE_DIR}/AMReX_DistributionMapping.cpp
     ${AMREX_SOURCE_DIR}/AMReX_ParallelDescriptor.cpp
     ${AMREX_SOURCE_DIR}/AMReX_ForkJoin.cpp
     ${AMREX_SOURCE_DIR}/AMReX_ParallelContext.cpp
     ${AMREX_SOURCE_DIR}/AMReX_VisMF.cpp 
     ${AMREX_SOURCE_DIR}/AMReX_Arena.cpp
     ${AMREX_SOURCE_DIR}/AMReX_BArena.cpp
     ${AMREX_SOURCE_DIR}/AMReX_CArena.cpp
     ${AMREX_SOURCE_DIR}/AMReX_DArena.cpp
     ${AMREX_SOURCE_DIR}/AMReX_EArena.cpp
     ${AMREX_SOURCE_DIR}/AMReX_NFiles.cpp  
     ${AMREX_SOURCE_DIR}/AMReX_parstream.cpp
     ${AMREX_SOURCE_DIR}/AMReX_FabConv.cpp  
     ${AMREX_SOURCE_DIR}/AMReX_FPC.cpp
     ${AMREX_SOURCE_DIR}/AMReX_VectorIO.cpp
     ${AMREX_SOURCE_DIR}/AMReX_IntConv.cpp
     ${AMREX_SOURCE_DIR}/AMReX_Box.cpp
     ${AMREX_SOURCE_DIR}/AMReX_BoxIterator.cpp
     ${AMREX_SOURCE_DIR}/AMReX_IntVect.cpp
     ${AMREX_SOURCE_DIR}/AMReX_IndexType.cpp
     ${AMREX_SOURCE_DIR}/AMReX_Orientation.cpp
     ${AMREX_SOURCE_DIR}/AMReX_Periodicity.cpp
     ${AMREX_SOURCE_DIR}/AMReX_RealBox.cpp
     ${AMREX_SOURCE_DIR}/AMReX_RealVect.cpp
     ${AMREX_SOURCE_DIR}/AMReX_BoxList.cpp 
     ${AMREX_SOURCE_DIR}/AMReX_BoxArray.cpp
     ${AMREX_SOURCE_DIR}/AMReX_BoxDomain.cpp
     ${AMREX_SOURCE_DIR}/AMReX_FArrayBox.cpp
     ${AMREX_SOURCE_DIR}/AMReX_IArrayBox.cpp 
     ${AMREX_SOURCE_DIR}/AMReX_BaseFab.cpp
     ${AMREX_SOURCE_DIR}/AMReX_MultiFab.cpp 
     ${AMREX_SOURCE_DIR}/AMReX_MFCopyDescriptor.cpp
     ${AMREX_SOURCE_DIR}/AMReX_iMultiFab.cpp
     ${AMREX_SOURCE_DIR}/AMReX_FabArrayBase.cpp
     ${AMREX_SOURCE_DIR}/AMReX_MFIter.cpp
     ${AMREX_SOURCE_DIR}/AMReX_CoordSys.cpp 
     ${AMREX_SOURCE_DIR}/AMReX_Geometry.cpp 
     ${AMREX_SOURCE_DIR}/AMReX_MultiFabUtil.cpp
     ${AMREX_SOURCE_DIR}/AMReX_BCRec.cpp
     ${AMREX_SOURCE_DIR}/AMReX_PhysBCFunct.cpp
     ${AMREX_SOURCE_DIR}/AMReX_BCUtil.cpp
     ${AMREX_SOURCE_DIR}/AMReX_PlotFileUtil.cpp
     ${AMREX_SOURCE_DIR}/AMReX_PlotFileDataImpl.cpp
     ${AMREX_SOURCE_DIR}/AMReX_BLutil_F.F90
     ${AMREX_SOURCE_DIR}/AMReX_BLProfiler_F.F90
     ${AMREX_SOURCE_DIR}/AMReX_FILCC_${AMREX_DIM}D.F90
     ${AMREX_SOURCE_DIR}/AMReX_filcc_mod.F90
     ${AMREX_SOURCE_DIR}/AMReX_omp_mod.F90
     ${AMREX_SOURCE_DIR}/AMReX_acc_mod.F90
     ${AMREX_SOURCE_DIR}/AMReX_fort_mod.F90
     ${AMREX_SOURCE_DIR}/AMReX_constants_mod.f90
     ${AMREX_SOURCE_DIR}/AMReX_error_mod.F90
     ${AMREX_SOURCE_DIR}/AMReX_parmparse_mod.F90
     ${AMREX_SOURCE_DIR}/AMReX_string_mod.F90
     ${AMREX_SOURCE_DIR}/AMReX_bc_types_mod.F90
     ${AMREX_SOURCE_DIR}/AMReX_ParallelDescriptor_F.F90
     ${AMREX_SOURCE_DIR}/AMReX_io_mod.F90
     ${AMREX_SOURCE_DIR}/AMReX_GpuControl.cpp
     ${AMREX_SOURCE_DIR}/AMReX_GpuLaunch.cpp
     ${AMREX_SOURCE_DIR}/AMReX_GpuDevice.cpp
     ${AMREX_SOURCE_DIR}/AMReX_GpuUtility.cpp
     ${AMREX_SOURCE_DIR}/AMReX_GpuAsyncArray.cpp
     ${AMREX_SOURCE_DIR}/AMReX_GpuElixir.cpp
     ${AMREX_SOURCE_DIR}/AMReX_CudaAllocators.cpp
     ${AMREX_SOURCE_DIR}/AMReX_Machine.cpp
     ${AMREX_SOURCE_DIR}/AMReX_MemPool.cpp
     ${AMREX_SOURCE_DIR}/AMReX_mempool_mod.F90
     ${AMREX_SOURCE_DIR}/AMReX_BLProfiler.cpp
     ${AMREX_SOURCE_DIR}/AMReX_BLBackTrace.cpp
   )
   #set(AMREX_SOURCE_DIR "${CMAKE_SOURCE_DIR}/submods/amrex/Src/F_Interfaces/AmrCore")
   #add_sources(GlobalSourceList
   #   ${AMREX_SOURCE_DIR}/AMReX_FAmrCore.cpp
   #   ${AMREX_SOURCE_DIR}/AMReX_amr_mod.F90
   #   ${AMREX_SOURCE_DIR}/AMReX_amrcore_fi.cpp
   #   ${AMREX_SOURCE_DIR}/AMReX_amrcore_mod.F90
   #   ${AMREX_SOURCE_DIR}/AMReX_fillpatch_fi.cpp
   #   ${AMREX_SOURCE_DIR}/AMReX_fillpatch_mod.F90
   #   ${AMREX_SOURCE_DIR}/AMReX_fluxregister_fi.cpp
   #   ${AMREX_SOURCE_DIR}/AMReX_fluxregister_mod.F90
   #   ${AMREX_SOURCE_DIR}/AMReX_interpolater_mod.F90
   #   ${AMREX_SOURCE_DIR}/AMReX_tagbox_fi.cpp
   #   ${AMREX_SOURCE_DIR}/AMReX_tagbox_mod.F90
   #)
   #set(AMREX_SOURCE_DIR "${CMAKE_SOURCE_DIR}/submods/amrex/Src/F_Interfaces/Base")
   #add_sources(GlobalSourceList
   #   ${AMREX_SOURCE_DIR}/AMReX_FPhysBC.cpp
   #   ${AMREX_SOURCE_DIR}/AMReX_base_mod.F90
   #   ${AMREX_SOURCE_DIR}/AMReX_box_fi.cpp
   #   ${AMREX_SOURCE_DIR}/AMReX_box_mod.F90
   #   ${AMREX_SOURCE_DIR}/AMReX_boxarray_fi.cpp
   #   ${AMREX_SOURCE_DIR}/AMReX_boxarray_mod.F90
   #   ${AMREX_SOURCE_DIR}/AMReX_distromap_fi.cpp
   #   ${AMREX_SOURCE_DIR}/AMReX_distromap_mod.F90
   #   ${AMREX_SOURCE_DIR}/AMReX_fab_mod.F90
   #   ${AMREX_SOURCE_DIR}/AMReX_geometry_fi.cpp
   #   ${AMREX_SOURCE_DIR}/AMReX_geometry_mod.F90
   #   ${AMREX_SOURCE_DIR}/AMReX_init_fi.cpp
   #   ${AMREX_SOURCE_DIR}/AMReX_init_mod.F90
   #   ${AMREX_SOURCE_DIR}/AMReX_multifab_fi.cpp
   #   ${AMREX_SOURCE_DIR}/AMReX_multifab_mod.F90
   #   ${AMREX_SOURCE_DIR}/AMReX_multifabutil_fi.cpp
   #   ${AMREX_SOURCE_DIR}/AMReX_multifabutil_mod.F90
   #   ${AMREX_SOURCE_DIR}/AMReX_parallel_mod.F90
   #   ${AMREX_SOURCE_DIR}/AMReX_physbc_fi.cpp
   #   ${AMREX_SOURCE_DIR}/AMReX_physbc_mod.F90
   #   ${AMREX_SOURCE_DIR}/AMReX_plotfile_fi.cpp
   #   ${AMREX_SOURCE_DIR}/AMReX_plotfile_mod.F90
   #   ${AMREX_SOURCE_DIR}/AMReX_vismf_fi.cpp
   #)
   if(AMREX_ENABLE_EB)
     set(AMREX_SOURCE_DIR "${CMAKE_SOURCE_DIR}/submods/amrex/Src/EB")
     add_sources(GlobalSourceList
       ${AMREX_SOURCE_DIR}/AMReX_EBAmrUtil.cpp
       ${AMREX_SOURCE_DIR}/AMReX_EBDataCollection.cpp
       ${AMREX_SOURCE_DIR}/AMReX_EBFArrayBox.cpp
       ${AMREX_SOURCE_DIR}/AMReX_EBInterpolater.cpp
       ${AMREX_SOURCE_DIR}/AMReX_EBCellFlag.cpp
       ${AMREX_SOURCE_DIR}/AMReX_EBFabFactory.cpp
       ${AMREX_SOURCE_DIR}/AMReX_EBFluxRegister.cpp  
       ${AMREX_SOURCE_DIR}/AMReX_EBMultiFabUtil.cpp
       ${AMREX_SOURCE_DIR}/AMReX_MultiCutFab.cpp
       ${AMREX_SOURCE_DIR}/AMReX_EB_levelset.cpp
       ${AMREX_SOURCE_DIR}/AMReX_EB_utils.cpp
       ${AMREX_SOURCE_DIR}/AMReX_EB_LSCoreBase.cpp 
       ${AMREX_SOURCE_DIR}/AMReX_EBFluxRegister_nd.F90
       ${AMREX_SOURCE_DIR}/AMReX_ebcellflag_mod.F90   
       ${AMREX_SOURCE_DIR}/AMReX_compute_normals.F90
       ${AMREX_SOURCE_DIR}/AMReX_eb_to_pvd.F90
       ${AMREX_SOURCE_DIR}/AMReX_EB_geometry.F90
       ${AMREX_SOURCE_DIR}/AMReX_EB_levelset_F.F90
       ${AMREX_SOURCE_DIR}/AMReX_EB_Tagging.F90
       ${AMREX_SOURCE_DIR}/AMReX_EB_bc_fill_nd.F90
       ${AMREX_SOURCE_DIR}/AMReX_distFcnElement.cpp
       ${AMREX_SOURCE_DIR}/AMReX_EB2.cpp
       ${AMREX_SOURCE_DIR}/AMReX_EB2_Level.cpp
       ${AMREX_SOURCE_DIR}/AMReX_EB2_MultiGFab.cpp
       ${AMREX_SOURCE_DIR}/AMReX_EB2_${AMREX_DIM}D_C.cpp
       ${AMREX_SOURCE_DIR}/AMReX_algoim.cpp
       ${AMREX_SOURCE_DIR}/AMReX_WriteEBSurface.cpp
     )
   endif()
   set(AMREX_SOURCE_DIR "${CMAKE_SOURCE_DIR}/submods/amrex/Src/Boundary")
   add_sources(GlobalSourceList
      ${AMREX_SOURCE_DIR}/AMReX_BndryData.cpp
      ${AMREX_SOURCE_DIR}/AMReX_BndryRegister.cpp
      ${AMREX_SOURCE_DIR}/AMReX_FabSet.cpp
      ${AMREX_SOURCE_DIR}/AMReX_InterpBndryData.cpp
      ${AMREX_SOURCE_DIR}/AMReX_LO_UTIL.F90
      ${AMREX_SOURCE_DIR}/AMReX_MacBndry.cpp
      ${AMREX_SOURCE_DIR}/AMReX_Mask.cpp
      ${AMREX_SOURCE_DIR}/AMReX_MultiMask.cpp
      ${AMREX_SOURCE_DIR}/AMReX_YAFluxRegister.cpp
      ${AMREX_SOURCE_DIR}/AMReX_lo_bctypes_mod.F90
   )
   set(AMREX_SOURCE_DIR "${CMAKE_SOURCE_DIR}/submods/amrex/Src/LinearSolvers/MLMG")
   add_sources(GlobalSourceList
      ${AMREX_SOURCE_DIR}/AMReX_MLMG.cpp
      ${AMREX_SOURCE_DIR}/AMReX_MLMGBndry.cpp
      ${AMREX_SOURCE_DIR}/AMReX_MLLinOp.cpp
      ${AMREX_SOURCE_DIR}/AMReX_MLLinOp_nd.F90
      ${AMREX_SOURCE_DIR}/AMReX_MLCellLinOp.cpp
      ${AMREX_SOURCE_DIR}/AMReX_MLNodeLinOp.cpp
      ${AMREX_SOURCE_DIR}/AMReX_MLCellABecLap.cpp
      ${AMREX_SOURCE_DIR}/AMReX_MLCGSolver.cpp
      ${AMREX_SOURCE_DIR}/AMReX_MLABecLaplacian.cpp
      ${AMREX_SOURCE_DIR}/AMReX_MLALaplacian.cpp
      ${AMREX_SOURCE_DIR}/AMReX_MLPoisson.cpp
      ${AMREX_SOURCE_DIR}/AMReX_MLNodeLaplacian.cpp
      ${AMREX_SOURCE_DIR}/AMReX_MLNodeLap_${AMREX_DIM}d.F90
      ${AMREX_SOURCE_DIR}/AMReX_MLTensorOp.cpp
   )
   if(AMREX_ENABLE_EB)
      add_sources(GlobalSourceList
         ${AMREX_SOURCE_DIR}/AMReX_MLEBABecLap.cpp
         ${AMREX_SOURCE_DIR}/AMReX_MLEBTensorOp.cpp
      )
   endif()
   set(AMREX_SOURCE_DIR "${CMAKE_SOURCE_DIR}/submods/amrex/Src/LinearSolvers/Projections")
   add_sources(GlobalSourceList
      ${AMREX_SOURCE_DIR}/AMReX_MacProjector.cpp
      ${AMREX_SOURCE_DIR}/AMReX_NodalProjector.cpp
   )
endfunction(get_amrex_sources)
