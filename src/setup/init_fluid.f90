module init_fluid_module

      use amrex_fort_module, only: rt => amrex_real
      use iso_c_binding,     only: c_int

contains
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: init_fluid                                              !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   subroutine init_fluid(slo, shi, lo, hi, &
                         domlo, domhi, p, vel, density, tracer, eta, & 
                         dx, dy, dz, xlength, ylength, zlength, probtype) &
      bind(C, name="init_fluid")

      use constant, only: ro_0, mu
      use constant, only: ic_u, ic_v, ic_w, ntrac

      implicit none

! Dummy arguments .....................................................//
      integer(c_int), intent(in   ) ::  lo(3),  hi(3)
      integer(c_int), intent(in   ) :: slo(3), shi(3)
      integer(c_int), intent(in   ) :: domlo(3),domhi(3)

      real(rt), intent(inout) ::       p(slo(1):shi(1), slo(2):shi(2), slo(3):shi(3)  )
      real(rt), intent(inout) ::     vel(slo(1):shi(1), slo(2):shi(2), slo(3):shi(3),3)
      real(rt), intent(inout) :: density(slo(1):shi(1), slo(2):shi(2), slo(3):shi(3))
      real(rt), intent(inout) ::  tracer(slo(1):shi(1), slo(2):shi(2), slo(3):shi(3),ntrac)
      real(rt), intent(inout) ::     eta(slo(1):shi(1), slo(2):shi(2), slo(3):shi(3)  )

      real(rt), intent(in   ) :: dx, dy, dz
      real(rt), intent(in   ) :: xlength, ylength, zlength
      integer,  intent(in   ) :: probtype

      ! Set the initial fluid density and viscosity
      vel(:,:,:,1) = ic_u
      vel(:,:,:,2) = ic_v
      vel(:,:,:,3) = ic_w

      density = ro_0
      tracer  = 0.0
      eta     = mu
      
      if (probtype ==  1) call taylor_green(lo, hi, vel, slo, shi, dx, dy, dz, domlo)
      if (probtype ==  2) call double_shear_layer(lo, hi, vel, slo, shi, dx, dy, dz, domlo)
      if (probtype ==  3) call taylor_green3d(lo, hi, vel, slo, shi, dx, dy, dz, domlo)
      if (probtype == 31) call plane_poiseuille(lo, hi, vel, tracer, slo, shi, dx, dy, dz, domlo, domhi, probtype)
      if (probtype == 32) call plane_poiseuille(lo, hi, vel, tracer, slo, shi, dx, dy, dz, domlo, domhi, probtype)
      if (probtype == 33) call plane_poiseuille(lo, hi, vel, tracer, slo, shi, dx, dy, dz, domlo, domhi, probtype)
      if (probtype ==  4) call couette(lo, hi, vel, slo, shi, dx, dy, dz, domlo, domhi)
      if (probtype == 11) call tuscan(lo, hi, vel, density, tracer, slo, shi, dx, dy, dz, domlo, domhi)

   end subroutine init_fluid

   subroutine taylor_green(lo, hi, vel, slo, shi, dx, dy, dz, domlo)

      use constant, only: zero, half, one

      implicit none

      ! Array bounds
      integer(c_int), intent(in   ) :: lo(3), hi(3)
      integer(c_int), intent(in   ) :: domlo(3)
      integer(c_int), intent(in   ) :: slo(3), shi(3)

      ! Arrays
      real(rt),       intent(inout) :: vel(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),3)

      real(rt),       intent(in   ) :: dx, dy, dz


      ! Local variables
      integer(c_int)                  :: i, j, k
      real(rt)                        :: x, y
      real(rt)                        :: twopi = 8.0_rt * atan(one)

      ! x-direction
      do j = lo(2), hi(2)
         y =  ( real(j,rt) + half ) * dy
         do i = lo(1), hi(1)
            x =  ( real(i,rt) + half ) * dx
            do k = lo(3), hi(3)
               vel(i,j,k,1) =  sin(twopi * x) * cos(twopi * y)
               vel(i,j,k,2) = -cos(twopi * x) * sin(twopi * y)
               vel(i,j,k,3) = zero
            end do
         end do
      end do

   end subroutine taylor_green

    subroutine taylor_green3d(lo, hi, vel, slo, shi, dx, dy, dz, domlo)

       use constant, only: zero, half, one

       implicit none

       ! Array bounds
       integer(c_int), intent(in   ) :: lo(3), hi(3)
       integer(c_int), intent(in   ) :: domlo(3)
       integer(c_int), intent(in   ) :: slo(3), shi(3)

       ! Arrays
       real(rt),       intent(inout) :: vel(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),3)

       real(rt),       intent(in   ) :: dx, dy, dz

       ! Local variables
       integer(c_int)                  :: i, j, k
       real(rt)                        :: x, y, z
       real(rt)                        :: twopi = 8.0_rt * atan(one)

       do j = lo(2), hi(2)
          y =  ( real(j,rt) + half ) * dy
          do i = lo(1), hi(1)
             x =  ( real(i,rt) + half ) * dx
             do k = lo(3), hi(3)
                z =  ( real(k,rt) + half ) * dz
                vel(i,j,k,1) =  sin(twopi * x) * cos(twopi * y) * cos(twopi * z)
                vel(i,j,k,2) = -cos(twopi * x) * sin(twopi * y) * cos(twopi * z)
                vel(i,j,k,3) = zero
             end do
          end do
       end do

    end subroutine taylor_green3d

   subroutine double_shear_layer(lo, hi, vel, slo, shi, dx, dy, dz, domlo)

      use constant, only: zero, half, one

      implicit none

      ! Array bounds
      integer(c_int), intent(in   ) :: lo(3), hi(3)
      integer(c_int), intent(in   ) :: domlo(3)
      integer(c_int), intent(in   ) :: slo(3), shi(3)

      ! Arrays
      real(rt),       intent(inout) :: vel(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),3)

      real(rt),       intent(in   ) :: dx, dy, dz

      ! Local variables
      integer(c_int)                :: i, j, k, plane
      real(rt)                      :: x, y, z
      real(rt)                      :: twopi = 8.0_rt * atan(one)

      plane = 1

      select case ( plane )

      case (1)  ! x-y plane

         ! x-direction
         do j = lo(2), hi(2)
            y =  ( real(j,rt) + half ) * dy
            do i = lo(1), hi(1)
               x =  ( real(i,rt) + half ) * dx
               do k = lo(3), hi(3)
                  vel(i,j,k,1) = tanh ( 30.0_rt * (0.25_rt - abs ( y - 0.5_rt ) ) )
                  vel(i,j,k,2) = 0.05_rt * sin ( twopi * x )
                  vel(i,j,k,3) = zero
               end do
            end do
         end do

      case (2)  ! x-z plane

         ! x-direction
         do k = lo(3), hi(3)
            z =  ( real(k,rt) + half ) * dz
            do i = lo(1), hi(1)
               x =  ( real(i,rt) + half ) * dx
               do j = lo(2), hi(2)
                  vel(i,j,k,1) = tanh ( 30.0_rt * (0.25_rt - abs ( z - 0.5_rt ) ) )
                  vel(i,j,k,2) = zero
                  vel(i,j,k,3) = 0.05_rt * sin ( twopi * x )
               end do
            end do
         end do

      case (3)  ! y-z plane

         ! x-direction
         do k = lo(3), hi(3)
            z =  ( real(k,rt) + half ) * dz
            do j = lo(2), hi(2)
               y =  ( real(j,rt) + half ) * dy
               do i = lo(1), hi(1)
                  vel(i,j,k,1) = zero
                  vel(i,j,k,2) = tanh ( 30.0_rt * (0.25_rt - abs ( z - 0.5_rt ) ) )
                  vel(i,j,k,3) = 0.05_rt * sin ( twopi * y )
               end do
            end do
         end do

      end select

   end subroutine double_shear_layer

   subroutine plane_poiseuille(lo, hi, vel, tracer, slo, shi, dx, dy, dz, domlo, domhi, probtype)

      use constant, only: ic_u, ic_v, ic_w, ntrac

      implicit none

      integer(c_int),   intent(in   ) ::    lo(3),    hi(3)
      integer(c_int),   intent(in   ) :: domlo(3), domhi(3)
      integer(c_int),   intent(in   ) ::   slo(3),   shi(3)
      integer(c_int),   intent(in   ) :: probtype

      real(rt),         intent(in   ) :: dx, dy, dz
      real(rt),         intent(inout) ::    vel(slo(1):shi(1), slo(2):shi(2), slo(3):shi(3), 3)
      real(rt),         intent(inout) :: tracer(slo(1):shi(1), slo(2):shi(2), slo(3):shi(3), ntrac)

      ! Local variables
      integer(c_int)                  :: i, j, k
      integer(c_int)                  :: num_cells
      real(rt)                        :: x, y, z

      if (probtype .eq. 31) then

         num_cells = domhi(2) - domlo(2) + 1

!        if (ntrac .ne. 3) then
!           print *,"ntrac here is ", ntrac
!           print *,"probtype = 31 is supposed to have 3 tracers!"
!           stop
!        end if

         do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               y =  (real(j,rt) + 0.5d0) / dble(num_cells)
               x =  (real(i,rt) + 0.5d0) / dble(num_cells)
               vel(i,j,k,1) = 6.0d0 * ic_u * y * (1.0d0 - y)
               vel(i,j,k,2) = 0.0d0
               vel(i,j,k,3) = 0.0d0
               tracer(i,j,k,1:ntrac) = 0.0d0
               if (ntrac .gt. 0 .and. i .le. domhi(1)/8) &
                  tracer(i,j,k,1) = 1.d0
               if (ntrac .gt. 1 .and. i .le. domhi(1)/2) &
                  tracer(i,j,k,2) = 2.d0
               if (ntrac .gt. 2 .and. i .le. 3*domhi(1)/4) &
                  tracer(i,j,k,3) = 3.d0
            end do
         end do
         end do

      else if (probtype .eq. 32) then

         num_cells = domhi(3) - domlo(3) + 1

         do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               z =  (real(k,rt) + 0.5d0) / dble(num_cells)
               vel(i,j,k,1) = 0.0d0
               vel(i,j,k,2) = 6.0d0 * ic_v * z * (1.0d0 - z)
               vel(i,j,k,3) = 0.0d0
               tracer(i,j,k,1:ntrac) = 0.0d0
               if (ntrac .gt. 0 .and. j .le. domhi(2)/8) &
                  tracer(i,j,k,1) = 1.d0
               if (ntrac .gt. 1 .and. j .le. domhi(2)/2) &
                  tracer(i,j,k,2) = 2.d0
               if (ntrac .gt. 2 .and. j .le. 3*domhi(2)/4) &
                  tracer(i,j,k,3) = 3.d0
            end do
         end do
         end do

      else if (probtype .eq. 33) then

         num_cells = domhi(1) - domlo(1) + 1

         do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               x =  (real(i,rt) + 0.5d0) / dble(num_cells)
               vel(i,j,k,1) = 0.0d0
               vel(i,j,k,2) = 0.0d0
               vel(i,j,k,3) = 6.0d0 * ic_w * x * (1.0d0 - x)
               tracer(i,j,k,1:ntrac) = 0.0d0
               if (ntrac .gt. 0 .and. k .le. domhi(3)/8) &
                  tracer(i,j,k,1) = 1.d0
               if (ntrac .gt. 1 .and. k .le. domhi(3)/2) &
                  tracer(i,j,k,2) = 2.d0
               if (ntrac .gt. 2 .and. k .le. 3*domhi(3)/4) &
                  tracer(i,j,k,3) = 3.d0
            end do
         end do
         end do

      else 
          print *,'Dont know this probtype in plane_poiseuille: ', probtype
          stop
      end if

   end subroutine plane_poiseuille

   subroutine couette(lo, hi, vel, slo, shi, dx, dy, dz, domlo, domhi)

      use constant,          only: zero, half, one
      use constant,          only: ic_u

      implicit none

      integer(c_int),   intent(in   ) ::    lo(3),    hi(3)
      integer(c_int),   intent(in   ) :: domlo(3), domhi(3)
      integer(c_int),   intent(in   ) ::   slo(3),   shi(3)

      real(rt),         intent(inout) :: vel(slo(1):shi(1), slo(2):shi(2), slo(3):shi(3), 3)

      real(rt),         intent(in   ) :: dx, dy, dz

      ! Local variables
      integer(c_int)                  :: i, j, k
      integer(c_int)                  :: num_cells_y
      real(rt)                        :: y

      num_cells_y = domhi(2) - domlo(2) + 1

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            y =  (real(j,rt) + half) / num_cells_y
            do i = lo(1), hi(1)
               vel(i,j,k,1) = ic_u * (y - half)
               vel(i,j,k,2) = zero
               vel(i,j,k,3) = zero
            end do
         end do
      end do

   end subroutine couette

   subroutine tuscan(lo, hi, vel, density, tracer, slo, shi, dx, dy, dz, domlo, domhi)

      use constant,          only: zero, half, one

      implicit none

      integer(c_int),   intent(in   ) ::    lo(3),    hi(3)
      integer(c_int),   intent(in   ) :: domlo(3), domhi(3)
      integer(c_int),   intent(in   ) ::   slo(3),   shi(3)

      real(rt),         intent(inout) ::     vel(slo(1):shi(1), slo(2):shi(2), slo(3):shi(3), 3)
      real(rt),         intent(inout) :: density(slo(1):shi(1), slo(2):shi(2), slo(3):shi(3))
      real(rt),         intent(inout) ::  tracer(slo(1):shi(1), slo(2):shi(2), slo(3):shi(3))

      real(rt),         intent(in   ) :: dx, dy, dz

      ! Local variables
      integer(c_int)                  :: i, j, k
      integer(c_int)                  :: num_cells
      real(rt)                        :: pert

      num_cells = domhi(3) - domlo(3) + 1

      pert = 0.01d0

      ! Start the flow at rest
      vel = zero

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               if (k .le. num_cells/2) then
                  density(i,j,k) = one
                  tracer (i,j,k) = zero
               else
                  density(i,j,k) = one
                  tracer (i,j,k) = zero + pert
               endif 
            end do
         end do
      end do

   end subroutine tuscan

end module init_fluid_module

