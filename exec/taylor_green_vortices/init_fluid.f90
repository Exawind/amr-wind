module init_fluid_module
contains
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: init_fluid                                              !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   subroutine init_fluid(slo, shi, lo, hi, &
                         domlo, domhi, ro, p, vel, &
                         eta, dx, dy, dz, xlength, ylength, zlength) &
      bind(C, name="init_fluid")

      use amrex_fort_module, only : rt => amrex_real
      use iso_c_binding , only: c_int

      use constant      , only: ro_0, mu

      implicit none

! Dummy arguments .....................................................//
      integer(c_int), intent(in   ) ::  lo(3),  hi(3)
      integer(c_int), intent(in   ) :: slo(3), shi(3)
      integer(c_int), intent(in   ) :: domlo(3),domhi(3)

      real(rt), intent(inout) :: ro&
                                 (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(rt), intent(inout) :: p&
                                 (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(rt), intent(inout) :: vel&
                                 (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),3)

      real(rt), intent(inout) :: eta&
                                 (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(rt), intent(in   ) :: dx, dy, dz
      real(rt), intent(in   ) :: xlength, ylength, zlength

      call init_taylor_green ( lo, hi, vel, slo, shi, dx, dy, dz, domlo)

      ! Set the initial fluid density and viscosity
      ro  = ro_0
      eta = mu

   end subroutine init_fluid

   subroutine init_taylor_green(lo, hi, vel, slo, shi, dx, dy, dz, domlo)

      use amrex_fort_module, only: ar => amrex_real
      use iso_c_binding ,    only: c_int
      use param,             only: zero, half, one

      implicit none

      ! Array bounds
      integer(c_int),   intent(in   ) :: slo(3), shi(3)

      ! Tile bounds
      integer(c_int),   intent(in   ) ::  lo(3),  hi(3)

      ! Grid and domain lower bound
      integer(c_int),   intent(in   ) :: domlo(3)
      real(ar),         intent(in   ) :: dx, dy, dz

      ! Arrays
      real(ar),         intent(inout) ::                   &
           & vel(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),3)

      ! Local variables
      integer(c_int)                  :: i, j, k
      real(ar)                        :: x, y, z
      real(ar)                        :: twopi = 8.0_ar * atan(one)

      ! x-direction
      do j = lo(2), hi(2)
         y =  ( real(j,ar) + half ) * dy
         do i = lo(1), hi(1)
            x =  ( real(i,ar) + half ) * dx
            do k = lo(3), hi(3)
               vel(i,j,k,1) =  sin(twopi * x) * cos(twopi * y)
               vel(i,j,k,2) = -cos(twopi * x) * sin(twopi * y)
               vel(i,j,k,3) = zero
            end do
         end do
      end do

   end subroutine init_taylor_green

end module init_fluid_module
