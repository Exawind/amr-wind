module convection_eb_mod
   
   use amrex_fort_module,       only: ar => amrex_real
   use iso_c_binding ,          only: c_int
   use param,                   only: zero, half, one
   use bc,                      only: minf_, nsw_, fsw_, psw_, pinf_, pout_
   use eb_interpolation_mod,    only: interpolate_to_face_centroid
   use amrex_ebcellflag_module, only: is_covered_cell, is_single_valued_cell, &
        &                             get_neighbor_cells

   implicit none
   private
   
contains


   !$************************************************
   !$ WARNING:
   !$
   !$ In compute_convective_fluxes_<dir>, we loop over
   !$ a tile + some layers of ghost nodes.
   !$ Make sure this approach is not problematic with
   !$ OpenMP
   !$    
   !$************************************************


   !$************************************************
   !$ WARNING 2:
   !$
   !$ For now the convective term is only div(u_mac . u)
   !$ The term -u*grad(u_mac) is omitted
   !$    
   !$************************************************  
    
   subroutine compute_convective_fluxes_x( lo, hi, &
        fx, fxlo, fxhi, &
        vel, vello, velhi, &
        u, ulo, uhi, &
        v, vlo, vhi, &
        w, wlo, whi, &
        xslopes, slo, shi, &
        areafrac, alo, ahi, &
        cent, clo, chi, &
        flags, flo, fhi, &
        domlo, domhi, &
        bc_ilo, bc_ihi, ng ) bind(C)

      ! Tile bounds (cell-centered)
      integer(c_int),  intent(in   ) :: lo(3),  hi(3)

      ! Array Bounds
      integer(c_int),  intent(in   ) :: slo(3),   shi(3)
      integer(c_int),  intent(in   ) :: fxlo(3),  fxhi(3)
      integer(c_int),  intent(in   ) :: vello(3), velhi(3)
      integer(c_int),  intent(in   ) :: ulo(3),   uhi(3)
      integer(c_int),  intent(in   ) :: vlo(3),   vhi(3)
      integer(c_int),  intent(in   ) :: wlo(3),   whi(3)
      integer(c_int),  intent(in   ) :: alo(3),   ahi(3)
      integer(c_int),  intent(in   ) :: clo(3),   chi(3)
      integer(c_int),  intent(in   ) :: flo(3),   fhi(3)
      integer(c_int),  intent(in   ) :: domlo(3), domhi(3), ng

      ! Arrays
      real(ar),        intent(in   ) ::                            &
           & vel(vello(1):velhi(1),vello(2):velhi(2),vello(3):velhi(3),3)    , &
           & xslopes(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),3), & 
           & u(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3)), &
           & v(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3)), &
           & w(wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3)), &
           & areafrac(alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3)),       &
           & cent(clo(1):chi(1),clo(2):chi(2),clo(3):chi(3),2)
           

      real(ar),        intent(  out) ::                           &
           & fx(fxlo(1):fxhi(1),fxlo(2):fxhi(2),fxlo(3):fxhi(3),3)

      ! BC types
      integer(c_int), intent(in   ) ::  &
           & bc_ilo(domlo(2)-ng:domhi(2)+ng,domlo(3)-ng:domhi(3)+ng,2), &
           & bc_ihi(domlo(2)-ng:domhi(2)+ng,domlo(3)-ng:domhi(3)+ng,2), &
           & flags(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3))

      ! Temporary arrays: CC velocity at faces
      real(ar) :: fx_tmp(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3),3)

      ! Local variables
      integer(c_int)                 :: i, j, k, n, lo_grow(3), hi_grow(3)
      real(ar)                       :: upls, umns, u_face
      integer, parameter             :: bc_list(6) = [MINF_, NSW_, FSW_, PSW_, PINF_, POUT_]

      ! 
      ! First compute the fluxes at the x-face center
      ! Do this on ALL x-faces on the tile, i.e. INCLUDE as many ghost faces as
      ! possible      
      !
      do n = 1, 3
         do k = lo(3)-ng, hi(3)+ng
            do j = lo(2)-ng, hi(2)+ng
               do i = lo(1)-ng+1, hi(1)+ng

                  if ( areafrac(i,j,k) > zero ) then
                     if ( i >= domlo(1) .and. any(bc_ilo(j,k,1) == bc_list) ) then
                        u_face =  vel(domlo(1)-1,j,k,n)
                     else if ( i >= domhi(1)+1 .and. any(bc_ihi(j,k,1) == bc_list ) ) then
                        u_face =  vel(domhi(1)+1,j,k,n)
                     else
                        upls  = vel(i  ,j,k,n) - half * xslopes(i  ,j,k,n)
                        umns  = vel(i-1,j,k,n) + half * xslopes(i-1,j,k,n)

                        u_face = upwind( umns, upls, u(i,j,k) )
                     end if                     
                  else
!                    u_face = huge(one)
                     u_face = 1.d200
                  end if

                  fx_tmp(i,j,k,n) = u(i,j,k) * u_face

!                 if ( k.eq.2.and.j.eq.5) print *,'MAKING FX ',i,n,fx_tmp(i,j,k,n),u(i,j,k),u_face
!                 if ( k.eq.2.and.j.eq.5) print *,'AREA FRAC ',areafrac(i,j,k)
!                 if ( k.eq.2.and.j.eq.5) print *,'  USING SLOPES ', vel(i,j,k,n), xslopes(i,j,k,n)
               end do
            end do
         end do
      end do

      ! 
      ! Interpolate at face centroid. Interpolated value now lives in fx
      ! Interpolation must be performed in the ghost nodes as well, except for the
      ! outermost halo layer in y and z (because of interpolation stencil limitations)
      lo_grow = lo + [-ng+1, -ng+1, -ng+1 ]
      hi_grow = hi + [ ng  ,  ng-1,  ng-1 ]
      
      call interpolate_to_face_centroid( lo_grow, hi_grow, fx, fx_tmp, ulo, uhi, 3, &
           areafrac, alo, ahi, cent, clo, chi, flags, flo, fhi, 1 )
      
   end subroutine compute_convective_fluxes_x

   subroutine compute_convective_fluxes_y( lo, hi, &
        fy, fylo, fyhi, &
        vel, vello, velhi, &
        u, ulo, uhi, &
        v, vlo, vhi, &
        w, wlo, whi, &
        yslopes, slo, shi, &
        areafrac, alo, ahi, &
        cent, clo, chi, &
        flags, flo, fhi, &
        domlo, domhi, &
        bc_jlo, bc_jhi, ng ) bind(C)

      ! Tile bounds (cell centered)
      integer(c_int),  intent(in   ) :: lo(3),  hi(3)

      ! Array Bounds
      integer(c_int),  intent(in   ) :: slo(3),   shi(3)
      integer(c_int),  intent(in   ) :: fylo(3),  fyhi(3)
      integer(c_int),  intent(in   ) :: vello(3), velhi(3)
      integer(c_int),  intent(in   ) :: ulo(3),   uhi(3)
      integer(c_int),  intent(in   ) :: vlo(3),   vhi(3)
      integer(c_int),  intent(in   ) :: wlo(3),   whi(3)
      integer(c_int),  intent(in   ) :: alo(3),   ahi(3)
      integer(c_int),  intent(in   ) :: clo(3),   chi(3)
      integer(c_int),  intent(in   ) :: flo(3),   fhi(3)
      integer(c_int),  intent(in   ) :: domlo(3), domhi(3), ng

      ! Arrays
      real(ar),        intent(in   ) ::                            &
           & vel(vello(1):velhi(1),vello(2):velhi(2),vello(3):velhi(3),3)    , &
           & yslopes(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),3), & 
           & u(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3)), &
           & v(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3)), &
           & w(wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3)), &
           & areafrac(alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3)),       &
           & cent(clo(1):chi(1),clo(2):chi(2),clo(3):chi(3),2)
           

      real(ar),        intent(  out) ::                           &
           & fy(fylo(1):fyhi(1),fylo(2):fyhi(2),fylo(3):fyhi(3),3)

      ! BC types
      integer(c_int), intent(in   ) ::  &
           & bc_jlo(domlo(1)-ng:domhi(1)+ng,domlo(3)-ng:domhi(3)+ng,2), &
           & bc_jhi(domlo(1)-ng:domhi(1)+ng,domlo(3)-ng:domhi(3)+ng,2), &
           & flags(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3))

      ! Temporary arrays: CC velocity at faces
      real(ar) :: fy_tmp(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3),3)

      ! Local variables
      integer(c_int)                 :: i, j, k, n, lo_grow(3), hi_grow(3)
      real(ar)                       :: vpls, vmns, v_face
      integer, parameter             :: bc_list(6) = [MINF_, NSW_, FSW_, PSW_, PINF_, POUT_]

      ! 
      ! First compute the fluxes at the y-face center
      ! Do this on ALL y-faces on the tile, i.e. INCLUDE as many ghost faces as
      ! possible      
      !
      do n = 1, 3
         do k = lo(3)-ng,     hi(3)+ng
            do j = lo(2)-ng+1,   hi(2)+ng
               do i = lo(1)-ng,    hi(1)+ng
                  if ( areafrac(i,j,k) > zero ) then
                     if ( j <= domlo(2) .and. any(bc_jlo(i,k,1) == bc_list) ) then
                        v_face =  vel(i,domlo(2)-1,k,n)
                     else if ( j >= domhi(2)+1 .and. any(bc_jhi(i,k,1) == bc_list ) ) then
                        v_face =  vel(i,domhi(2)+1,k,n)
                     else
                        vpls  = vel(i,j  ,k,n) - half * yslopes(i,j  ,k,n)
                        vmns  = vel(i,j-1,k,n) + half * yslopes(i,j-1,k,n)

                        v_face = upwind( vmns, vpls, v(i,j,k) )
                     end if
                  else
!                    v_face = huge(one)
                     v_face = 1.d200
                  end if
                  fy_tmp(i,j,k,n) = v(i,j,k) * v_face
               end do
            end do
         end do
      end do

      ! 
      ! Interpolate at face centroid. Interpolated value are returned in fy
      ! Interpolation must be performed in the ghost nodes as well, except for the
      ! outermost halo layer in x and z (because of interpolation stencil limitations)
      lo_grow = lo + [-ng+1, -ng+1, -ng+1 ]
      hi_grow = hi + [ ng-1,  ng  ,  ng-1 ]
       
      call interpolate_to_face_centroid( lo_grow, hi_grow, fy, fy_tmp, vlo, vhi, 3, &
           areafrac, alo, ahi, cent, clo, chi, flags, flo, fhi, 2 )
      
   end subroutine compute_convective_fluxes_y

   subroutine compute_convective_fluxes_z( lo, hi, &
        fz, fzlo, fzhi, &
        vel, vello, velhi, &
        u, ulo, uhi, &
        v, vlo, vhi, &
        w, wlo, whi, &
        zslopes, slo, shi, &
        areafrac, alo, ahi, &
        cent, clo, chi, &
        flags, flo, fhi, &
        domlo, domhi, &
        bc_klo, bc_khi, ng ) bind(C)

      ! Tile bounds (cell centered)
      integer(c_int),  intent(in   ) :: lo(3),  hi(3)

      ! Array Bounds
      integer(c_int),  intent(in   ) :: slo(3),   shi(3)
      integer(c_int),  intent(in   ) :: fzlo(3),  fzhi(3)
      integer(c_int),  intent(in   ) :: vello(3), velhi(3)
      integer(c_int),  intent(in   ) :: ulo(3),   uhi(3)
      integer(c_int),  intent(in   ) :: vlo(3),   vhi(3)
      integer(c_int),  intent(in   ) :: wlo(3),   whi(3)
      integer(c_int),  intent(in   ) :: alo(3),   ahi(3)
      integer(c_int),  intent(in   ) :: clo(3),   chi(3)
      integer(c_int),  intent(in   ) :: flo(3),   fhi(3)
      integer(c_int),  intent(in   ) :: domlo(3), domhi(3), ng

      ! Arrays
      real(ar),        intent(in   ) ::                            &
           & vel(vello(1):velhi(1),vello(2):velhi(2),vello(3):velhi(3),3)    , &
           & zslopes(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),3), & 
           & u(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3)), &
           & v(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3)), &
           & w(wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3)), &
           & areafrac(alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3)),       &
           & cent(clo(1):chi(1),clo(2):chi(2),clo(3):chi(3),2)
           

      real(ar),        intent(  out) ::                           &
           & fz(fzlo(1):fzhi(1),fzlo(2):fzhi(2),fzlo(3):fzhi(3),3)

      ! BC types
      integer(c_int), intent(in   ) ::  &
           & bc_klo(domlo(1)-ng:domhi(1)+ng,domlo(2)-ng:domhi(2)+ng,2), &
           & bc_khi(domlo(1)-ng:domhi(1)+ng,domlo(2)-ng:domhi(2)+ng,2), &
           & flags(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3))

      ! Temporary arrays: CC velocity at faces
      real(ar) :: fz_tmp(wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3),3)

      ! Local variables
      integer(c_int)                 :: i, j, k, n, lo_grow(3), hi_grow(3)
      real(ar)                       :: wpls, wmns, w_face
      integer, parameter             :: bc_list(6) = [MINF_, NSW_, FSW_, PSW_, PINF_, POUT_]

      ! 
      ! First compute the fluxes at the z-face center
      ! Do this on ALL z-faces on the tile, i.e. INCLUDE as many ghost faces as
      ! possible      
      !
      do n = 1, 3
         do k = lo(3)-ng+1, hi(3)+ng
            do j = lo(2)-ng, hi(2)+ng
               do i = lo(1)-ng, hi(1)+ng
                  if ( areafrac(i,j,k) > zero ) then
                     if ( k <= domlo(3) .and. any(bc_klo(i,j,1) == bc_list) ) then
                        w_face =  vel(i,j,domlo(3)-1,n)
                     else if ( k >= domhi(3)+1 .and. any(bc_khi(i,j,1) == bc_list ) ) then
                        w_face =  vel(i,j,domhi(3)+1,n)
                     else
                        wpls  = vel(i,j,k  ,n) - half * zslopes(i,j,k  ,n)
                        wmns  = vel(i,j,k-1,n) + half * zslopes(i,j,k-1,n)

                        w_face = upwind( wmns, wpls, w(i,j,k) )
                     end if
                  else
!                    w_face = huge(one)
                     w_face = 1.d200
                  end if
                  fz_tmp(i,j,k,n) = w(i,j,k) * w_face
               end do
            end do
         end do
      end do

      ! 
      ! Interpolate at face centroid. Interpolated value now lives in fz
      ! Interpolation must be performed in the ghost nodes as well, except for the
      ! outermost halo layer in x and y (because of interpolation stencil limitations)
      ! 
      lo_grow = lo + [-ng+1, -ng+1, -ng+1 ]
      hi_grow = hi + [ ng-1,  ng-1,  ng   ]
      
      call interpolate_to_face_centroid( lo_grow, hi_grow, fz, fz_tmp, wlo, whi, 3, &
           areafrac, alo, ahi, cent, clo, chi, flags, flo, fhi, 3 )
     
   end subroutine compute_convective_fluxes_z
   
   subroutine compute_ugradu_eb ( lo, hi, &
        ugradu, glo, ghi, &
        vel, vello, velhi, &
        u, ulo, uhi, &
        v, vlo, vhi, &
        w, wlo, whi, &
        afrac_x, axlo, axhi, &
        afrac_y, aylo, ayhi, &
        afrac_z, azlo, azhi, &
        cent_x,  cxlo, cxhi, &
        cent_y,  cylo, cyhi, &
        cent_z,  czlo, czhi, &
        flags,    flo,  fhi, &
        vfrac,   vflo, vfhi, &
        xslopes, yslopes, zslopes, slo, shi, &
        domlo, domhi, &
        bc_ilo, bc_ihi, &
        bc_jlo, bc_jhi, &
        bc_klo, bc_khi, dx, ng ) bind(C)


      ! Tile bounds
      integer(c_int),  intent(in   ) :: lo(3),  hi(3)

      ! Array Bounds
      integer(c_int),  intent(in   ) :: slo(3), shi(3)
      integer(c_int),  intent(in   ) :: glo(3), ghi(3)
      integer(c_int),  intent(in   ) :: vello(3), velhi(3)
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
      integer(c_int),  intent(in   ) :: domlo(3), domhi(3), ng

      ! Grid
      real(ar),        intent(in   ) :: dx(3)

      ! Arrays
      real(ar),        intent(in   ) ::                            &
           & vel(vello(1):velhi(1),vello(2):velhi(2),vello(3):velhi(3),3)    , &
           & xslopes(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),3), &
           & yslopes(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),3), &
           & zslopes(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),3), &
           & u(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3)), &
           & v(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3)), &
           & w(wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3)), &
           & afrac_x(axlo(1):axhi(1),axlo(2):axhi(2),axlo(3):axhi(3)),       &
           & afrac_y(aylo(1):ayhi(1),aylo(2):ayhi(2),aylo(3):ayhi(3)),       &
           & afrac_z(azlo(1):azhi(1),azlo(2):azhi(2),azlo(3):azhi(3)), &
           & cent_x(cxlo(1):cxhi(1),cxlo(2):cxhi(2),cxlo(3):cxhi(3),2),&
           & cent_y(cylo(1):cyhi(1),cylo(2):cyhi(2),cylo(3):cyhi(3),2),&
           & cent_z(czlo(1):czhi(1),czlo(2):czhi(2),czlo(3):czhi(3),2),&
           & vfrac(vflo(1):vfhi(1),vflo(2):vfhi(2),vflo(3):vfhi(3))           
      
      real(ar),        intent(  out) ::                           &
           & ugradu(glo(1):ghi(1),glo(2):ghi(2),glo(3):ghi(3),3)

      ! BC types
      integer(c_int), intent(in   ) ::  &
           & bc_ilo(domlo(2)-ng:domhi(2)+ng,domlo(3)-ng:domhi(3)+ng,2), &
           & bc_ihi(domlo(2)-ng:domhi(2)+ng,domlo(3)-ng:domhi(3)+ng,2), &
           & bc_jlo(domlo(1)-ng:domhi(1)+ng,domlo(3)-ng:domhi(3)+ng,2), &
           & bc_jhi(domlo(1)-ng:domhi(1)+ng,domlo(3)-ng:domhi(3)+ng,2), &
           & bc_klo(domlo(1)-ng:domhi(1)+ng,domlo(2)-ng:domhi(2)+ng,2), &
           & bc_khi(domlo(1)-ng:domhi(1)+ng,domlo(2)-ng:domhi(2)+ng,2), &
           & flags(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3))
      
      ! Temporary arrays
      ! Fluxes (this array can probably be shrinked to tile size)
      real(ar)  ::    &
           & fx(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3),3), &
           & fy(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3),3), &
           & fz(wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3),3)
      ! Conservative div and EB stuff
      real(ar)  ::    &
           &  divc(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng), &
           & optmp(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng), &
           &  delm(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng)


      ! Local variables
      integer(c_int)                 :: i, j, k, n, nbr(-1:1,-1:1,-1:1)
      real(ar)                       :: idx, idy, idz

      idx = one / dx(1)
      idy = one / dx(2)
      idz = one / dx(3)

      ! Check number of ghost cells
      if (ng < 3) then
         write(*,*) "ERROR: EB convection term requires at least 3 ghost cells"
         stop
      end if

      ! First we compute the fluxes
      call compute_convective_fluxes_x( lo, hi, fx, ulo, uhi, vel, vello, velhi, &
           u, ulo, uhi, v, vlo, vhi, w, wlo, whi,                  &
           xslopes, slo, shi, afrac_x, axlo, axhi, cent_x, cxlo, cxhi,               &
           flags, flo, fhi, domlo, domhi, bc_ilo, bc_ihi, ng )
      
      call compute_convective_fluxes_y( lo, hi, fy, vlo, vhi, vel, vello, velhi, &
           u, ulo, uhi, v, vlo, vhi, w, wlo, whi,                  &
           yslopes, slo, shi, afrac_y, aylo, ayhi, cent_y, cylo, cyhi,               &
           flags, flo, fhi, domlo, domhi, bc_jlo, bc_jhi, ng )
      
      call compute_convective_fluxes_z( lo, hi, fz, wlo, whi, vel, vello, velhi, &
           u, ulo, uhi, v, vlo, vhi, w, wlo, whi,                  &
           zslopes, slo, shi, afrac_z, azlo, azhi, cent_z, czlo, czhi,               &
           flags, flo, fhi, domlo, domhi, bc_klo, bc_khi, ng )      


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
                        fxp = fx(i+1,j,k,n) * afrac_x(i+1,j,k)
                        fxm = fx(i  ,j,k,n) * afrac_x(i  ,j,k)
                  
                        fyp = fy(i,j+1,k,n) * afrac_y(i,j+1,k)
                        fym = fy(i,j  ,k,n) * afrac_y(i,j  ,k)
                        
                        fzp = fz(i,j,k+1,n) * afrac_z(i,j,k+1)
                        fzm = fz(i,j,k  ,n) * afrac_z(i,j,k  )
                        
                        divc(i,j,k) = ( ( fxp - fxm ) * idx + ( fyp - fym ) * idy + &
                             &          ( fzp - fzm ) * idz ) / vfrac(i,j,k) 
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
            real(ar)  :: epvfrac 

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
                                    epvfrac = vfrac(i+ii,j+jj,k+kk) 
                                    vtot    = vtot  + epvfrac
                                    divnc   = divnc + epvfrac*divc(i+ii,j+jj,k+kk)
                                 end if
                              end do
                           end do
                        end do
                        divnc   = divnc / vtot
                        epvfrac = vfrac(i,j,k)
                        optmp(i,j,k) = (one - epvfrac) * ( divnc-divc(i,j,k) )
                        delm(i,j,k) = -epvfrac * optmp(i,j,k)
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
                  ugradu(i,j,k,n) = -( divc(i,j,k) + optmp(i,j,k) )        
               end do
            end do
         end do
         
      end do ncomp_loop
      
   end subroutine compute_ugradu_eb



   ! Upwind non-normal velocity
   function upwind ( velmns, velpls, uedge ) result (ev)

      ! Small value to protect against tiny velocities used in upwinding
      real(ar),        parameter     :: small_vel = 1.0d-10

      real(ar), intent(in) :: velmns, velpls, uedge
      real(ar)             :: ev

      if ( abs(uedge) .lt. small_vel) then
         ev = half * ( velpls + velmns )
      else 
         ev = merge ( velmns, velpls, uedge >= zero ) 
      end if

   end function upwind
   
   
end module convection_eb_mod
