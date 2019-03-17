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
                         domlo, domhi, ro, p, vel, eta, & 
                         dx, dy, dz, xlength, ylength, zlength, probtype) &
      bind(C, name="init_fluid")

      use constant, only: ro_0, mu
      use constant, only: ic_u, ic_v, ic_w

      implicit none

! Dummy arguments .....................................................//
      integer(c_int), intent(in   ) ::  lo(3),  hi(3)
      integer(c_int), intent(in   ) :: slo(3), shi(3)
      integer(c_int), intent(in   ) :: domlo(3),domhi(3)

      real(rt), intent(inout) ::  ro(slo(1):shi(1), slo(2):shi(2), slo(3):shi(3)  )
      real(rt), intent(inout) ::   p(slo(1):shi(1), slo(2):shi(2), slo(3):shi(3)  )
      real(rt), intent(inout) :: vel(slo(1):shi(1), slo(2):shi(2), slo(3):shi(3),3)
      real(rt), intent(inout) :: eta(slo(1):shi(1), slo(2):shi(2), slo(3):shi(3)  )

      real(rt), intent(in   ) :: dx, dy, dz
      real(rt), intent(in   ) :: xlength, ylength, zlength
      integer,  intent(in   ) :: probtype

      ! Set the initial fluid density and viscosity
      ro  = ro_0
      eta = mu
      vel(:,:,:,1) = ic_u
      vel(:,:,:,2) = ic_v
      vel(:,:,:,3) = ic_w
      
      if(probtype == 1) call taylor_green(lo, hi, vel, slo, shi, dx, dy, dz, domlo)
      if(probtype == 2) call double_shear_layer(lo, hi, vel, slo, shi, dx, dy, dz, domlo)
      if(probtype == 3) call plane_poiseuille(lo, hi, vel, slo, shi, dx, dy, dz, domlo, domhi)

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
      real(rt)                        :: x, y, z
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

   subroutine plane_poiseuille(lo, hi, vel, slo, shi, dx, dy, dz, domlo, domhi)

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

      do j = lo(2), hi(2)
         y =  (real(j,rt) + half) / num_cells_y
         do i = lo(1), hi(1)
            do k = lo(3), hi(3)
               vel(i,j,k,1) = 6.0 * ic_u * y * (one - y)
               vel(i,j,k,2) = zero
               vel(i,j,k,3) = zero
            end do
         end do
      end do

   end subroutine plane_poiseuille

end module init_fluid_module

