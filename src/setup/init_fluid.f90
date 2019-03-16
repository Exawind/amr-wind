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
                         domlo, domhi, ro, p, vel, &
                         eta, dx, dy, dz, xlength, ylength, zlength) &
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

      ! Set the initial fluid density and viscosity
      ro  = ro_0
      eta = mu
      vel(:,:,:,1) = ic_u
      vel(:,:,:,2) = ic_v
      vel(:,:,:,3) = ic_w

   end subroutine init_fluid

end module init_fluid_module
