#include <incflo.H>

using namespace amrex;

void
incflo::incflo_correct_small_cells (Vector<MultiFab*      > const& vel_in,
                                    Vector<MultiFab const*> const& u_mac,
                                    Vector<MultiFab const*> const& v_mac,
                                    Vector<MultiFab const*> const& w_mac)

{
#ifdef AMREX_USE_EB
    BL_PROFILE("incflo::incflo_correct_small_cells");

    for (int lev = 0; lev <= finest_level; lev++)
    {
       for (MFIter mfi(*vel_in[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
       {
          // Tilebox
          const Box bx = mfi.tilebox();

          EBCellFlagFab const& flags = EBFactory(lev).getMultiEBCellFlagFab()[mfi];

          // Face-centered velocity components
          const auto& umac_fab = (u_mac[lev])->array(mfi);
          const auto& vmac_fab = (v_mac[lev])->array(mfi);
          const auto& wmac_fab = (w_mac[lev])->array(mfi);

          if (flags.getType(amrex::grow(bx,0)) == FabType::covered )
          {
            // do nothing
          }
  
          // No cut cells in this FAB
          else if (flags.getType(amrex::grow(bx,1)) == FabType::regular )
          {
            // do nothing
          }

          // Cut cells in this FAB
          else
          {
             // Face-centered areas
             const auto& apx_fab   = EBFactory(lev).getAreaFrac()[0]->const_array(mfi);
             const auto& apy_fab   = EBFactory(lev).getAreaFrac()[1]->const_array(mfi);
             const auto& apz_fab   = EBFactory(lev).getAreaFrac()[2]->const_array(mfi);
             const auto& vfrac_fab = EBFactory(lev).getVolFrac().const_array(mfi);

             const auto& ccvel_fab = vel_in[lev]->array(mfi);

             // This FAB has cut cells -- we define the centroid value in terms of the MAC velocities onfaces
             amrex::ParallelFor(bx,
               [vfrac_fab,apx_fab,apy_fab,apz_fab,ccvel_fab,umac_fab,vmac_fab,wmac_fab]
               AMREX_GPU_DEVICE (int i, int j, int k) noexcept
             {
                 if (vfrac_fab(i,j,k) > 0.0 && vfrac_fab(i,j,k) < 1.e-4)
                 {
                    Real u_avg = (apx_fab(i,j,k) * umac_fab(i,j,k) + apx_fab(i+1,j,k) * umac_fab(i+1,j,k)) 
                               / (apx_fab(i,j,k) + apx_fab(i+1,j,k));
                    Real v_avg = (apy_fab(i,j,k) * vmac_fab(i,j,k) + apy_fab(i,j+1,k) * vmac_fab(i,j+1,k)) 
                               / (apy_fab(i,j,k) + apy_fab(i,j+1,k));
                    Real w_avg = (apz_fab(i,j,k) * wmac_fab(i,j,k) + apz_fab(i,j,k+1) * wmac_fab(i,j,k+1)) 
                               / (apz_fab(i,j,k) + apz_fab(i,j,k+1));

                    ccvel_fab(i,j,k,0) = u_avg;
                    ccvel_fab(i,j,k,1) = v_avg;
                    ccvel_fab(i,j,k,2) = w_avg;

                 }
             });
          }
       }
    }
#endif
}
