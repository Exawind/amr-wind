module divop_mod

   use amrex_fort_module,       only: ar => amrex_real
   use iso_c_binding ,          only: c_int
   use param,                   only: zero, half, one, my_huge
   use amrex_error_module,      only: amrex_abort
   use amrex_ebcellflag_module, only: is_covered_cell, is_single_valued_cell, &
        &                             get_neighbor_cells

   implicit none

contains

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
                             vel,vllo,vlhi, &
                             fx, fxlo, fxhi, &
                             fy, fylo, fyhi, &
                             fz, fzlo, fzhi, &
                             afrac_x, axlo, axhi, &
                             afrac_y, aylo, ayhi, &
                             afrac_z, azlo, azhi, &
                             cent_x,  cxlo, cxhi, &
                             cent_y,  cylo, cyhi, &
                             cent_z,  czlo, czhi, &
                             flags,    flo,  fhi, &
                             vfrac,   vflo, vfhi, &
                             bcent,    blo,  bhi, &
                             domlo, domhi,        &
                             dx, ng, eta,          &
                             do_explicit_diffusion ) bind(C)

      use bc
      use eb_interpolation_mod, only: interp_to_face_centroid
      use eb_wallflux_mod,      only: compute_diff_wallflux

      ! Tile bounds (cell centered)
      integer(c_int),  intent(in   ) :: lo(3),  hi(3)

      ! Array Bounds
      integer(c_int),  intent(in   ) :: dlo(3), dhi(3)
      integer(c_int),  intent(in   ) :: vllo(3), vlhi(3)
      integer(c_int),  intent(in   ) :: fxlo(3), fxhi(3)
      integer(c_int),  intent(in   ) :: fylo(3), fyhi(3)
      integer(c_int),  intent(in   ) :: fzlo(3), fzhi(3)
      integer(c_int),  intent(in   ) :: axlo(3), axhi(3)
      integer(c_int),  intent(in   ) :: aylo(3), ayhi(3)
      integer(c_int),  intent(in   ) :: azlo(3), azhi(3)
      integer(c_int),  intent(in   ) :: cxlo(3), cxhi(3)
      integer(c_int),  intent(in   ) :: cylo(3), cyhi(3)
      integer(c_int),  intent(in   ) :: czlo(3), czhi(3)
      integer(c_int),  intent(in   ) ::  flo(3),  fhi(3)
      integer(c_int),  intent(in   ) :: vflo(3), vfhi(3)
      integer(c_int),  intent(in   ) ::  blo(3),  bhi(3)

      integer(c_int),  intent(in   ) :: domlo(3), domhi(3)

      ! Number of ghost cells
      integer(c_int),  intent(in   ) :: ng

      ! Grid
      real(ar),        intent(in   ) :: dx(3)

      ! Arrays
      real(ar),        intent(in   ) ::                       &
           & fx(fxlo(1):fxhi(1),fxlo(2):fxhi(2),fxlo(3):fxhi(3),3), &
           & fy(fylo(1):fyhi(1),fylo(2):fyhi(2),fylo(3):fyhi(3),3), &
           & fz(fzlo(1):fzhi(1),fzlo(2):fzhi(2),fzlo(3):fzhi(3),3), &
           & afrac_x(axlo(1):axhi(1),axlo(2):axhi(2),axlo(3):axhi(3)), &
           & afrac_y(aylo(1):ayhi(1),aylo(2):ayhi(2),aylo(3):ayhi(3)), &
           & afrac_z(azlo(1):azhi(1),azlo(2):azhi(2),azlo(3):azhi(3)), &
           & cent_x(cxlo(1):cxhi(1),cxlo(2):cxhi(2),cxlo(3):cxhi(3),2),&
           & cent_y(cylo(1):cyhi(1),cylo(2):cyhi(2),cylo(3):cyhi(3),2),&
           & cent_z(czlo(1):czhi(1),czlo(2):czhi(2),czlo(3):czhi(3),2),&
           & vfrac(vflo(1):vfhi(1),vflo(2):vfhi(2),vflo(3):vfhi(3)),   &
           &   vel(vllo(1):vlhi(1),vllo(2):vlhi(2),vllo(3):vlhi(3),3), &
           & bcent(blo(1):bhi(1),blo(2):bhi(2),blo(3):bhi(3),3)

      ! Optional arrays (only for viscous calculations)
      real(ar),        intent(in   ), optional  ::                &
           &     eta(vflo(1):vfhi(1),vflo(2):vfhi(2),vflo(3):vfhi(3))

      real(ar),        intent(inout) ::                           &
           & div(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),3)

      ! BC types
      integer(c_int), intent(in   ) ::  &
           & flags(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3))

      ! If true  then we include all the diffusive terms in this explicit result
      ! If false then we include all only the off-diagonal terms here -- we do this
      !     by computing the full tensor then subtracting the diagonal terms
      integer(c_int),  intent(in   ), optional :: do_explicit_diffusion

      ! Conservative div and EB stuff
      real(ar)  ::    &
           &  divc(lo(1)-2:hi(1)+2,lo(2)-2:hi(2)+2,lo(3)-2:hi(3)+2), &
           & optmp(lo(1)-2:hi(1)+2,lo(2)-2:hi(2)+2,lo(3)-2:hi(3)+2), &
           &  delm(lo(1)-2:hi(1)+2,lo(2)-2:hi(2)+2,lo(3)-2:hi(3)+2), &
           &  mask(lo(1)-2:hi(1)+2,lo(2)-2:hi(2)+2,lo(3)-2:hi(3)+2)

      real(ar), allocatable          :: divdiff_w(:,:)
      integer(c_int)                 :: i, j, k, n, nbr(-1:1,-1:1,-1:1)
      integer(c_int)                 :: nwalls
      real(ar)                       :: idx, idy, idz
      logical                        :: is_viscous

      idx = one / dx(1)
      idy = one / dx(2)
      idz = one / dx(3)

      ! Check number of ghost cells
      if (ng < 5) call amrex_abort( "compute_divop(): ng must be >= 5")

      ! Check if we are computing divergence for viscous term
      if ( present(eta) ) then
         is_viscous = .true.
      else
         is_viscous = .false.
      end if

      if ( abs(dx(1) - dx(2)) > epsilon(0.0_ar) .or.&
           abs(dx(1) - dx(3)) > epsilon(0.0_ar) .or.&
           abs(dx(3) - dx(2)) > epsilon(0.0_ar) ) then
         call amrex_abort("compute_divop(): grid spacing must be uniform")
      end if

      !
      ! Allocate arrays to host viscous wall fluxes
      !
      nwalls = 0
      if (is_viscous) then
         do k = lo(3)-2, hi(3)+2
            do j = lo(2)-2, hi(2)+2
               do i = lo(1)-2, hi(1)+2
                  if (is_single_valued_cell(flags(i,j,k))) nwalls = nwalls + 1
               end do
            end do
         end do
         allocate( divdiff_w(3,nwalls) )
         divdiff_w = zero
      end if

      !
      ! Array "mask" is used to sever the link to ghost cells when the BCs are not periodic
      ! It is set to 1 when the a cell can be used in computations, 0 otherwise
      !
      do k = lo(3)-2, hi(3)+2
         do j = lo(2)-2, hi(2)+2
            do i = lo(1)-2, hi(1)+2
               if ( ( .not. cyclic_x .and. (i < domlo(1) .or. i > domhi(1)) ) .or. &
                    ( .not. cyclic_y .and. (j < domlo(2) .or. j > domhi(2)) ) .or. &
                    ( .not. cyclic_z .and. (k < domlo(3) .or. k > domhi(3)) ) ) then
                  mask(i,j,k) = zero
               else
                  mask(i,j,k) = one
               end if
            end do
         end do
      end do

      !
      ! We use the EB algorithmm to compute the divergence at cell centers
      !
      ncomp_loop: do n = 1, 3

         !
         ! Step 1: compute conservative divergence on stencil (lo-2,hi+2)
         !
         compute_divc: block
            real(ar) :: fxp, fxm, fyp, fym, fzp, fzm
            integer  :: iwall

            iwall = 0

            divc = zero

            do k = lo(3)-2, hi(3)+2
               do j = lo(2)-2, hi(2)+2
                  do i = lo(1)-2, hi(1)+2

                     if (is_covered_cell(flags(i,j,k))) then

                        divc(i,j,k) = my_huge

                     else if (is_single_valued_cell(flags(i,j,k))) then

                        call get_neighbor_cells( flags(i,j,k), nbr )

                        ! interp_to_face_centroid returns the proper flux multiplied
                        ! by the face area
                        fxp = interp_to_face_centroid( i+1, j, k, 1, fx, fxlo, n,  &
                             & afrac_x, axlo, cent_x, cxlo, nbr )

                        fxm = interp_to_face_centroid( i  , j, k, 1, fx, fxlo, n,  &
                             & afrac_x, axlo, cent_x, cxlo, nbr )

                        fyp = interp_to_face_centroid( i, j+1, k, 2, fy, fylo, n,  &
                             & afrac_y, aylo, cent_y, cylo, nbr )

                        fym = interp_to_face_centroid( i, j, k, 2, fy, fylo, n,  &
                             & afrac_y, aylo, cent_y, cylo, nbr )

                        fzp = interp_to_face_centroid( i, j, k+1, 3, fz, fzlo, n,  &
                             & afrac_z, azlo, cent_z, czlo, nbr )

                        fzm = interp_to_face_centroid( i, j, k, 3, fz, fzlo, n,  &
                             & afrac_z, azlo, cent_z, czlo, nbr )

                        divc(i,j,k) = ( ( fxp - fxm ) * idx + &
                             &          ( fyp - fym ) * idy + &
                             &          ( fzp - fzm ) * idz ) / vfrac(i,j,k)

                        ! Add viscous wall fluxes (compute three components only
                        ! during the first pass, i.e. for n=1
                        iwall = iwall + 1
                        if (is_viscous) then
                           if (n==1) then
                              call compute_diff_wallflux(divdiff_w(:,iwall), &
                                                         dx, i, j, k, &
                                                         vel, vllo, vlhi, &
                                                         eta, vflo, vfhi, &
                                                         bcent, blo, bhi, &
                                                         afrac_x, axlo, axhi, &
                                                         afrac_y, aylo, ayhi, &
                                                         afrac_z, azlo, azhi, &
                                                         do_explicit_diffusion)
                           end if
                           divc(i,j,k) = divc(i,j,k) - divdiff_w(n,iwall) / (dx(n) * vfrac(i,j,k))
                        end if

                     else

                        divc(i,j,k) = ( fx(i+1,j  ,k  ,n) - fx(i,j,k,n) ) * idx  &
                             &      + ( fy(i  ,j+1,k  ,n) - fy(i,j,k,n) ) * idy  &
                             &      + ( fz(i  ,j  ,k+1,n) - fz(i,j,k,n) ) * idz

                     end if

                  end do
               end do
            end do
         end block compute_divc

         !
         ! Step 2: compute delta M ( mass gain or loss ) on (lo-1,hi+1)
         !
         optmp = zero
         block
            integer   :: ii, jj, kk
            real(ar)  :: divnc, vtot
            real(ar)  :: vfrac_mask

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
                                    vfrac_mask  = vfrac(i+ii,j+jj,k+kk) * mask(i+ii,j+jj,k+kk)
                                    vtot        = vtot  + vfrac_mask
                                    divnc       = divnc + vfrac_mask * divc(i+ii,j+jj,k+kk)
                                 end if
                              end do
                           end do
                        end do
                        divnc = divnc / vtot
                        optmp(i,j,k) = (one - vfrac(i,j,k)) * ( divnc - divc(i,j,k) )
                        delm(i,j,k) = - vfrac(i,j,k) * optmp(i,j,k)
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
                                    wtot = wtot + vfrac(i+ii,j+jj,k+kk) * mask(i+ii,j+jj,k+kk)
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
                                    optmp(i+ii,j+jj,k+kk) = optmp(i+ii,j+jj,k+kk) &
                                                          + delm(i,j,k) * wtot &
                                                          * mask(i+ii,j+jj,k+kk)
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
         ! Resume the correct sign, AKA return the negative
         ! ****************************************************
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  div(i,j,k,n) = divc(i,j,k) + optmp(i,j,k)
               end do
            end do
         end do

      end do ncomp_loop

      !
      ! Delete working arrays
      ! 
      if (is_viscous) deallocate(divdiff_w)

   end subroutine compute_divop

end module divop_mod
