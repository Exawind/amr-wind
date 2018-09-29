!
!
!  This module contains the subroutines used directly in advance()
!
!
module advance_mod

   use amrex_fort_module, only: ar => amrex_real
   use iso_c_binding ,    only: c_int
   use param,             only: zero, half, one

   implicit none
   private

   ! Define here the unit vectors
   ! This is used to shift index  based on how the variable is staggered
   ! Check e_x, e_y and e_z in incflo_level.H
   integer(c_int), parameter :: e_i(3,3) = reshape ( [1,0,0,0,1,0,0,0,1], [3,3] )

contains

   !
   ! Compute new dt by using the formula derived in
   ! "A Boundary Condition Capturing Method for Multiphase Incompressible Flow"
   ! by Kang et al. (JCP).
   !
   !  dt/2 * ( C+V + sqrt( (C+V)**2 + 4Fx/dx + 4Fy/dy + 4Fz/dz )
   !
   ! where
   !
   ! C = max(|U|)/dx + max(|V|)/dy + max(|W|)/dz    --> Convection
   !
   ! V = 2 * max(mu/ro) * (1/dx^2 + 1/dy^2 +1/dz^2) --> Diffusion
   !
   ! Fx, Fy, Fz = net acceleration due to external forces
   !
   ! WARNING: We use a slightly modified version of C in the implementation below
   !
   subroutine compute_new_dt ( umax, vmax, wmax, romin, mumax, gradp0max, &
        dx, cfl, steady_state, time, stop_time, dt ) &
        & bind(C)

      use constant, only: gravity

      integer(c_int), intent(in   ) :: steady_state
      real(ar),       intent(in   ) :: umax, vmax, wmax
      real(ar),       intent(in   ) :: mumax, romin
      real(ar),       intent(in   ) :: dx(3), cfl
      real(ar),       intent(in   ) :: gradp0max(3)
      real(ar),       intent(in   ) :: time, stop_time
      real(ar),       intent(inout) :: dt
      real(ar)                      :: old_dt
      real(ar)                      :: c_cfl, v_cfl, f_cfl
      real(ar)                      :: odx, ody, odz
      real(ar)                      :: tmp
      real(ar),       parameter     :: two = 2.0_ar, four = two*two
      real(ar),       parameter     :: eps = epsilon (zero)

      odx    = one / dx(1)
      ody    = one / dx(2)
      odz    = one / dx(3)
      c_cfl  = zero
      v_cfl  = zero
      f_cfl  = zero
      old_dt = dt

      ! Convection
      c_cfl = max ( umax*odx, vmax*ody, wmax*odz )

      ! Viscous
      v_cfl = two * ( mumax / romin ) * ( odx**2 + ody**2 + odz**2 )

      ! Gravity and/or gradient of p0
      f_cfl = abs(gravity(1)-gradp0max(1)) * odx + &
              abs(gravity(2)-gradp0max(2)) * ody + &
              abs(gravity(3)-gradp0max(3)) * odz

      ! Put all together
      tmp = (c_cfl + v_cfl)  + sqrt ( (c_cfl + v_cfl)**2 + four * f_cfl )
      dt  = cfl * two / tmp

      ! Protect against tmp very small
      ! This may happen, for example, when the initial velocity field
      ! is zero for an inviscid flow with no external forcing
      if ( tmp <= eps ) then
         dt = .5 * old_dt
      end if

      ! Don't let the timestep grow by more than 1% per step.
      if (old_dt > 0.0) &
         dt = min ( dt, 1.01*old_dt )

      ! Don't overshoot the final time if not running to steady state
      if (steady_state .eq. 0 .and. stop_time .ge. 0.) then
         if (time+dt .gt. stop_time) &
              dt = stop_time - time
      end if

   end subroutine compute_new_dt

   subroutine compute_gradp0_max ( lo, hi, p0, slo, shi, gp0_max, dx, nodal_pressure) &
        bind (C)

      ! Loop bounds
      integer(c_int),  intent(in   ) :: lo(3),  hi(3)

      ! Array bounds
      integer(c_int),  intent(in   ) :: slo(3), shi(3)

      ! Grid and time spacing
      real(ar),        intent(in   ) :: dx(3)
      real(ar),        intent(  out) :: gp0_max(3)

      ! Arrays
      real(ar),        intent(in   ) ::                       &
           p0(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      integer(c_int),  intent(in   ) :: nodal_pressure

      ! Local variables
      integer(c_int) :: i, j, k
      real(ar)       :: odx, ody, odz

      odx = one / dx(1)
      ody = one / dx(2)
      odz = one / dx(3)

      gp0_max(:) = zero

      if ( nodal_pressure == 1 ) then

         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)

                  gp0_max(1) = max( gp0_max(1), abs( &
                       p0(i+1,j,k) + p0(i+1,j+1,k) + p0(i+1,j,k+1) + p0(i+1,j+1,k+1) &
                       - p0(i  ,j,k) - p0(i  ,j+1,k) - p0(i  ,j,k+1) - p0(i  ,j+1,k+1) ) )

                  gp0_max(2) = max( gp0_max(2), abs( &
                       p0(i,j+1,k) + p0(i+1,j+1,k) + p0(i,j+1,k+1) + p0(i+1,j+1,k+1) &
                       - p0(i,j  ,k) - p0(i+1,j  ,k) - p0(i,j  ,k+1) - p0(i+1,j  ,k+1) ) )

                  gp0_max(3) = max( gp0_max(3), abs( &
                       p0(i,j,k+1) + p0(i+1,j,k+1) + p0(i,j+1,k+1) + p0(i+1,j+1,k+1) &
                       - p0(i,j,k  ) - p0(i+1,j,k  ) - p0(i,j+1,k  ) - p0(i+1,j+1,k  ) ) )


               end do
            end do
         end do

         gp0_max(1) = gp0_max(1) * odx * 0.25d0
         gp0_max(2) = gp0_max(2) * ody * 0.25d0
         gp0_max(3) = gp0_max(3) * odz * 0.25d0

      else

         ! The box is cell-centered
         do k = lo(3), hi(3)+1
            do j = lo(2), hi(2)+1
               do i = lo(1), hi(1)+1

                  gp0_max(1) = max( gp0_max(1), abs(p0(i,j,k) - p0(i-1,j,k)))
                  gp0_max(2) = max( gp0_max(2), abs(p0(i,j,k) - p0(i,j-1,k)))
                  gp0_max(3) = max( gp0_max(3), abs(p0(i,j,k) - p0(i,j,k-1)))

               end do
            end do
         end do

         gp0_max(1) = gp0_max(1) * odx
         gp0_max(2) = gp0_max(2) * ody
         gp0_max(3) = gp0_max(3) * odz

      end if

   end subroutine compute_gradp0_max

   !
   ! Add forcing (acceleration) terms to velocity
   !
   subroutine add_forcing ( lo, hi, vel, ulo, uhi, &
        & ro, slo, shi, domlo, domhi, dx, dt )  bind(C)

      use constant, only: gravity

      ! Loop bounds
      integer(c_int), intent(in   ) ::  lo(3), hi(3)

      ! Array bounds
      integer(c_int), intent(in   ) :: slo(3), shi(3)
      integer(c_int), intent(in   ) :: ulo(3), uhi(3)

      ! Time step width
      real(ar),       intent(in   ) :: dt

      ! Domain bounds
      integer(c_int), intent(in   ) :: domlo(3), domhi(3)

      ! Grid
      real(ar),       intent(in   ) :: dx(3)

      ! Arrays
      real(ar),       intent(in   ) :: &
           ro(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(ar),       intent(inout) :: &
           vel(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3),3)

      ! Local variables
      integer(c_int)                :: i, j, k , n

      do n = 1, 3
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)

                  vel(i,j,k,n) = vel(i,j,k,n) + dt * gravity(n) 

               end do
            end do
         end do
      end do

   end subroutine add_forcing

   !
   ! Compute the cell-centered divergence of  u
   !
   subroutine compute_diveucc ( lo, hi, diveu, slo, shi, vel, ulo, uhi, dx) &
        bind(C)

      ! Loop bounds
      integer(c_int), intent(in   ) ::  lo(3), hi(3)

      ! Arrays bounds
      integer(c_int), intent(in   ) :: slo(3),shi(3)
      integer(c_int), intent(in   ) :: ulo(3),uhi(3)

      ! Grid
      real(ar),       intent(in   ) :: dx(3)

      ! Array
      real(ar),       intent(  out) :: &
           diveu(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(ar),       intent(in   ) :: &
           vel(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3),3)

      ! Local variables
      integer  :: i, j, k
      real(ar) :: odx, ody, odz
      real(ar) :: eu_e, eu_w, ev_n, ev_s, ew_t, ew_b

      odx = one / dx(1)
      ody = one / dx(2)
      odz = one / dx(3)

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)

               ! Face values
               eu_e = half * (  vel(i+1,j,k,1) +  vel(i,j,k,1) )
               eu_w = half * (  vel(i-1,j,k,1) +  vel(i,j,k,1) )

               ev_n = half * (  vel(i,j+1,k,2) +  vel(i,j,k,2) )
               ev_s = half * (  vel(i,j-1,k,2) +  vel(i,j,k,2) )

               ew_t = half * (  vel(i,j,k+1,3) +  vel(i,j,k,3) )
               ew_b = half * (  vel(i,j,k-1,3) +  vel(i,j,k,3) )

               ! Divergence
               diveu(i,j,k) = (eu_e - eu_w) * odx + (ev_n - ev_s) * ody + &
                    &         (ew_t - ew_b) * odz
            end do
         end do
      end do
      stop

   end subroutine compute_diveucc

   subroutine compute_diveund ( lo, hi, diveu, slo, shi, vec, ulo, uhi, dx) &
        bind(C)

      ! Loop bounds
      integer(c_int), intent(in   ) ::  lo(3), hi(3)

      ! Arrays bounds
      integer(c_int), intent(in   ) :: slo(3),shi(3)
      integer(c_int), intent(in   ) :: ulo(3),uhi(3)

      ! Grid
      real(ar),       intent(in   ) :: dx(3)

      ! Array
      real(ar),       intent(  out) :: &
           diveu(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(ar),       intent(in   ) :: &
           vec(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3),3)

      ! Local variables
      integer  :: i, j, k
      real(ar) :: odx, ody, odz
      real(ar) :: eu_x, eu_y, eu_z

      odx = one / dx(1)
      ody = one / dx(2)
      odz = one / dx(3)

      ! Note these are nodal indices
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)

               ! Divergence
               eu_x = ( vec(i  ,j  ,k  ,1) + vec(i  ,j-1,k  ,1) &
                    +vec(i  ,j  ,k-1,1) + vec(i  ,j-1,k-1,1) &
                    -vec(i-1,j  ,k  ,1) - vec(i-1,j-1,k  ,1) &
                    -vec(i-1,j  ,k-1,1) - vec(i-1,j-1,k-1,1) )

               eu_y = ( vec(i  ,j  ,k  ,2) + vec(i-1,j  ,k  ,2) &
                    +vec(i  ,j  ,k-1,2) + vec(i-1,j  ,k-1,2) &
                    -vec(i  ,j-1,k  ,2) - vec(i-1,j-1,k  ,2) &
                    -vec(i  ,j-1,k-1,2) - vec(i-1,j-1,k-1,2) )

               eu_z = ( vec(i  ,j  ,k  ,3) + vec(i-1,j  ,k  ,3) &
                    +vec(i  ,j-1,k  ,3) + vec(i-1,j-1,k  ,3) &
                    -vec(i  ,j  ,k-1,3) - vec(i-1,j  ,k-1,3) &
                    -vec(i  ,j-1,k-1,3) - vec(i-1,j-1,k-1,3) )

               diveu(i,j,k) = 0.25d0 * (eu_x*odx + eu_y*ody + eu_z*odz)

            end do
         end do
      end do

   end subroutine compute_diveund

   !
   ! Average to faces in chosen direction  -- note we only average the "idir"th 
   !    component of cc onto the idir'th face
   ! 
   subroutine average_cc_to_fc ( lo, hi, fx, fxlo, fxhi, fy, fylo, fyhi, fz, fzlo, fzhi,  &
                                 cc, slo, shi) bind(C) 

      ! Loop bounds (assumed face centered!)
      integer(c_int), intent(in   ) ::  lo(3), hi(3)

      ! Arrays bounds
      integer(c_int), intent(in   ) :: fxlo(3),fxhi(3)
      integer(c_int), intent(in   ) :: fylo(3),fyhi(3)
      integer(c_int), intent(in   ) :: fzlo(3),fzhi(3)
      integer(c_int), intent(in   ) :: slo(3),shi(3)

      ! Array
      real(ar),       intent(inout) :: &
           fx(fxlo(1):fxhi(1),fxlo(2):fxhi(2),fxlo(3):fxhi(3)), &
           fy(fylo(1):fyhi(1),fylo(2):fyhi(2),fylo(3):fyhi(3)), &
           fz(fzlo(1):fzhi(1),fzlo(2):fzhi(2),fzlo(3):fzhi(3))

      real(ar),       intent(in   ) :: &
           cc(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),3)
      
      ! Local variables
      integer  :: i, j, k

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)+1
               fx(i,j,k) = half * ( cc(i-1,j,k,1) + cc(i,j,k,1) )  
            end do
         end do
      end do

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)+1
            do i = lo(1), hi(1)
               fy(i,j,k) = half * ( cc(i,j-1,k,2) + cc(i,j,k,2) )  
            end do
         end do
      end do

      do k = lo(3), hi(3)+1
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               fz(i,j,k) = half * ( cc(i,j,k-1,3) + cc(i,j,k,3) )  
            end do
         end do
      end do

   end subroutine average_cc_to_fc

end module advance_mod
