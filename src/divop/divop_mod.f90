module divop_mod

   use amrex_fort_module,       only: ar => amrex_real
   use iso_c_binding ,          only: c_int
   use param,                   only: zero, half, one
   use amrex_ebcellflag_module, only: is_covered_cell, is_single_valued_cell, &
        &                             get_neighbor_cells

   implicit none

contains

   ! 
   ! Compute the divergence operator for EB geometries
   !
   ! OUTPUTS:
   !
   !   div    the cell-centered divergence
   !
   ! INPUTS:
   !
   !   fx       the fluxes at the x-face CENTER (not centroid!)
   !   fy       the fluxes at the y-face CENTER (not centroid!)
   !   fz       the fluxes at the z-face CENTER (not centroid!)
   !
   ! WARNING: fx, fy, fz HAS to be filled with at least 3 GHOST nodes
   !          
   ! 
   ! 
   subroutine compute_divop ( lo, hi, &
        div, dlo, dhi, &
        fx,  ulo, uhi, &
        fy,  vlo, vhi, &
        fz,  wlo, whi, &
        afrac_x, axlo, axhi, &
        afrac_y, aylo, ayhi, &
        afrac_z, azlo, azhi, &
        cent_x,  cxlo, cxhi, &
        cent_y,  cylo, cyhi, &
        cent_z,  czlo, czhi, &
        flags,    flo,  fhi, &
        vfrac,   vflo, vfhi, &
        dx, ng ) bind(C)

      use eb_interpolation_mod,    only: interpolate_to_face_centroid

      ! Tile bounds (cell centered)
      integer(c_int),  intent(in   ) :: lo(3),  hi(3)

      ! Array Bounds
      integer(c_int),  intent(in   ) :: dlo(3), dhi(3)
      integer(c_int),  intent(in   ) :: ulo(3), uhi(3)
      integer(c_int),  intent(in   ) :: vlo(3), vhi(3)
      integer(c_int),  intent(in   ) :: wlo(3), whi(3)
      integer(c_int),  intent(in   ) :: axlo(3), axhi(3)
      integer(c_int),  intent(in   ) :: aylo(3), ayhi(3)
      integer(c_int),  intent(in   ) :: azlo(3), azhi(3)
      integer(c_int),  intent(in   ) :: cxlo(3), cxhi(3)
      integer(c_int),  intent(in   ) :: cylo(3), cyhi(3)
      integer(c_int),  intent(in   ) :: czlo(3), czhi(3)
      integer(c_int),  intent(in   ) ::  flo(3),  fhi(3)
      integer(c_int),  intent(in   ) :: vflo(3), vfhi(3)
      integer(c_int),  intent(in   ) :: ng

      ! Grid
      real(ar),        intent(in   ) :: dx(3)

      ! Arrays
      real(ar),        intent(in   ) ::                       &
           & fx(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3),3), &
           & fy(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3),3), &
           & fz(wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3),3), &
           & afrac_x(axlo(1):axhi(1),axlo(2):axhi(2),axlo(3):axhi(3)), &
           & afrac_y(aylo(1):ayhi(1),aylo(2):ayhi(2),aylo(3):ayhi(3)), &
           & afrac_z(azlo(1):azhi(1),azlo(2):azhi(2),azlo(3):azhi(3)), &
           & cent_x(cxlo(1):cxhi(1),cxlo(2):cxhi(2),cxlo(3):cxhi(3),2),&
           & cent_y(cylo(1):cyhi(1),cylo(2):cyhi(2),cylo(3):cyhi(3),2),&
           & cent_z(czlo(1):czhi(1),czlo(2):czhi(2),czlo(3):czhi(3),2),&
           & vfrac(vflo(1):vfhi(1),vflo(2):vfhi(2),vflo(3):vfhi(3))           

      real(ar),        intent(inout) ::                           &
           & div(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),3)

      ! BC types
      integer(c_int), intent(in   ) ::  &
           & flags(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3))

      
      ! Conservative div and EB stuff
      real(ar)  ::    &
           &  divc(lo(1)-2:hi(1)+2,lo(2)-2:hi(2)+2,lo(3)-2:hi(3)+2), &
           & optmp(lo(1)-2:hi(1)+2,lo(2)-2:hi(2)+2,lo(3)-2:hi(3)+2), &
           &  delm(lo(1)-2:hi(1)+2,lo(2)-2:hi(2)+2,lo(3)-2:hi(3)+2)


      ! Local variables
      real(ar)          :: &
           & fx_c(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3),3), &
           & fy_c(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3),3), &
           & fz_c(wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3),3)
           
      integer(c_int)                 :: i, j, k, n, nbr(-1:1,-1:1,-1:1)
      integer(c_int)                 :: lo_grow(3), hi_grow(3)
      real(ar)                       :: idx, idy, idz

      idx = one / dx(1)
      idy = one / dx(2)
      idz = one / dx(3)

      ! Check number of ghost cells
      if (ng < 4) then
         write(*,*) "ERROR: EB divop requires at least 3 ghost cells"
         stop
      end if

      !
      ! Fist we interpolate the fluxes at the face centroid
      ! We need to have the fluxes defined on the faces of all the cells in the tile
      ! plus two layers of ghost cells (cfr computation of divc below). 
      !
      lo_grow = lo - 2 
      
      ! X direction
      hi_grow = hi + [ 3, 2, 2 ]
      call interpolate_to_face_centroid( lo_grow, hi_grow, fx_c, fx, ulo, uhi, 3, &
           afrac_x, axlo, axhi, cent_x, cxlo, cxhi, flags, flo, fhi, 1 )
      
      ! Y direction
      hi_grow = hi + [ 2, 3, 2 ]
      call interpolate_to_face_centroid( lo_grow, hi_grow, fy_c, fy, vlo, vhi, 3, &
           afrac_y, aylo, ayhi, cent_y, cylo, cyhi, flags, flo, fhi, 2 )

      ! Z direction
      hi_grow = hi + [ 2, 2, 3 ]      
      call interpolate_to_face_centroid( lo_grow, hi_grow, fz_c, fz, wlo, whi, 3, &
           afrac_z, azlo, azhi, cent_z, czlo, czhi, flags, flo, fhi, 3 )
      

      !
      ! The we use the EB algorithmm to compute the divergence at cell centers
      ! 
      ncomp_loop: do n = 1, 3

         !
         ! Step 1: compute conservative divergence on stencil (lo-2,hi+2)
         !
         block
            real(ar)   :: fxp, fxm, fyp, fym, fzp, fzm

            do k = lo(3)-2, hi(3)+2
               do j = lo(2)-2, hi(2)+2
                  do i = lo(1)-2, hi(1)+2
                     if (is_covered_cell(flags(i,j,k))) then
                        divc(i,j,k) = huge(one)
                     else 
                        fxp = fx_c(i+1,j,k,n) * afrac_x(i+1,j,k)
                        fxm = fx_c(i  ,j,k,n) * afrac_x(i  ,j,k)

                        fyp = fy_c(i,j+1,k,n) * afrac_y(i,j+1,k)
                        fym = fy_c(i,j  ,k,n) * afrac_y(i,j  ,k)

                        fzp = fz_c(i,j,k+1,n) * afrac_z(i,j,k+1)
                        fzm = fz_c(i,j,k  ,n) * afrac_z(i,j,k  )

                        divc(i,j,k) = ( ( fxp - fxm ) * idx + &
                             &          ( fyp - fym ) * idy + &
                             &          ( fzp - fzm ) * idz ) &
                             &        / vfrac(i,j,k) 
                     end if
                  end do
               end do
            end do
         end block

         !
         ! Step 2: compute delta M ( mass gain or loss ) on (lo-1,hi+1)
         !
         optmp = zero
         block
            integer   :: ii, jj, kk
            real(ar)  :: divnc, vtot

            do k = lo(3)-1, hi(3)+1
               do j = lo(2)-1, hi(2)+1
                  do i = lo(1)-1, hi(1)+1
                     if (is_single_valued_cell(flags(i,j,k))) then
                        divnc = zero
                        vtot  = zero
                        call get_neighbor_cells(flags(i,j,k),nbr)
                        do kk = -1, 1
                           do jj = -1, 1
                              do ii = -1, 1
                                 ! Check if we have to include also cell (i,j,k) itself
                                 if ( ( ii /= 0 .or. jj /= 0 .or. kk /= 0) &
                                      .and. (nbr(ii,jj,kk)==1) ) then
                                    vtot    = vtot  + vfrac(i+ii,j+jj,k+kk) 
                                    divnc   = divnc + vfrac(i+ii,j+jj,k+kk)*divc(i+ii,j+jj,k+kk)
                                 end if
                              end do
                           end do
                        end do
                        divnc   = divnc / vtot
                        optmp(i,j,k) = (one - vfrac(i,j,k)) * ( divnc-divc(i,j,k) )
                        delm(i,j,k) = -vfrac(i,j,k) * optmp(i,j,k)
                     else
                        delm(i,j,k) = zero
                     end if
                  end do
               end do
            end do

         end block

         !
         ! Step 3: redistribute excess/loss of mass
         !
         block 
            real(ar) :: wtot
            integer  :: ii, jj, kk

            do k = lo(3)-1, hi(3)+1
               do j = lo(2)-1, hi(2)+1
                  do i = lo(1)-1, hi(1)+1
                     if (is_single_valued_cell(flags(i,j,k))) then
                        wtot = zero
                        call get_neighbor_cells(flags(i,j,k),nbr)
                        do kk = -1, 1
                           do jj = -1, 1
                              do ii = -1, 1
                                 ! Check if we have to include also cell (i,j,k) itself
                                 if ( ( ii /= 0 .or. jj /= 0 .or. kk /= 0) &
                                      .and. (nbr(ii,jj,kk)==1) ) then
                                    wtot = wtot + vfrac(i+ii,j+jj,k+kk)
                                 end if
                              end do
                           end do
                        end do

                        wtot = one/wtot

                        do kk = -1, 1
                           do jj = -1, 1
                              do ii = -1, 1
                                 ! Check if we have to include also cell (i,j,k) itself
                                 if ( ( ii /= 0 .or. jj /= 0 .or. kk /= 0) &
                                      .and. (nbr(ii,jj,kk)==1) ) then
                                    optmp(i+ii,j+jj,k+kk) = optmp(i+ii,j+jj,k+kk) + delm(i,j,k) * wtot
                                 end if
                              end do
                           end do
                        end do
                     end if
                  end do
               end do
            end do
         end block


         ! ****************************************************
         ! Return the negative 
         ! ****************************************************
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  div(i,j,k,n) = divc(i,j,k) + optmp(i,j,k)
               end do
            end do
         end do

      end do ncomp_loop

   end subroutine compute_divop


end module divop_mod
