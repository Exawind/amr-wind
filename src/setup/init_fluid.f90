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

      ! Set user specified initial conditions (IC)
      call set_ic(slo, shi, domlo, domhi, dx, dy, dz, vel)

      ! Set the initial fluid density and viscosity
      ro  = ro_0
      eta = mu

   end subroutine init_fluid

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: init_fluid_restart                                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   subroutine init_fluid_restart(slo, shi, lo, hi, eta) &
      bind(C, name="init_fluid_restart")

      use amrex_fort_module, only : rt => amrex_real
      use iso_c_binding , only: c_int
      use constant, only: mu

      implicit none

! Dummy arguments .....................................................//
      integer(c_int), intent(in   ) ::  lo(3),  hi(3)
      integer(c_int), intent(in   ) :: slo(3), shi(3)

      real(rt), intent(inout) :: eta&
                                 (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      eta = mu

   end subroutine init_fluid_restart

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: SET_IC                                                  !
!  Author: M. Syamlal                                 Date: 21-JAN-92  !
!                                                                      !
!  Purpose: This module sets all the initial conditions.               !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   subroutine set_ic(slo, shi, domlo, domhi, dx, dy, dz, vel)

      use ic, only: dim_ic, ic_defined
      use ic, only: ic_p, ic_u, ic_v, ic_w
      use ic, only: ic_x_e, ic_y_n, ic_z_t
      use ic, only: ic_x_w, ic_y_s, ic_z_b
      use param, only: undefined, is_defined

      use amrex_fort_module, only : rt => amrex_real
      use iso_c_binding , only: c_int

      use calc_cell_module, only: calc_cell_ic

      implicit none

      integer(c_int), intent(in   ) :: slo(3), shi(3)
      integer(c_int), intent(in   ) :: domlo(3),domhi(3)
      real(rt), intent(in   ) :: dx, dy, dz

      real(rt), intent(inout) ::  vel&
                                 (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),3)

!-----------------------------------------------
! Local variables
!-----------------------------------------------
! indices
      integer :: istart, iend
      integer :: jstart, jend
      integer :: kstart, kend
      ! local index for initial condition
      integer :: icv

      ! Temporary variables for storing IC values
      real(rt) :: pgx, ugx, vgx, wgx

      integer :: i_w, j_s, k_b
      integer :: i_e, j_n, k_t

!  Set the initial conditions.
      do icv = 1, dim_ic
         if (ic_defined(icv)) then

            call calc_cell_ic(dx, dy, dz, &
                              ic_x_w(icv), ic_y_s(icv), ic_z_b(icv), &
                              ic_x_e(icv), ic_y_n(icv), ic_z_t(icv), &
                              i_w, i_e, j_s, j_n, k_b, k_t)

            pgx = ic_p(icv)
            ugx = ic_u(icv)
            vgx = ic_v(icv)
            wgx = ic_w(icv)

            if (is_defined(ugx)) then
               istart = max(slo(1), i_w)
               jstart = max(slo(2), j_s)
               kstart = max(slo(3), k_b)
               iend   = min(shi(1), i_e)
               jend   = min(shi(2), j_n)
               kend   = min(shi(3), k_t)
               vel(istart:iend,jstart:jend,kstart:kend,1) = ugx
               if (slo(1).lt.domlo(1) .and. domlo(1) == istart) &
                  vel(slo(1):istart-1,jstart:jend,kstart:kend,1) = ugx
               if (shi(1).gt.domhi(1) .and. domhi(1) == iend  ) &
                  vel(iend+1:shi(1)  ,jstart:jend,kstart:kend,1) = ugx
            end if

            if (is_defined(vgx)) then
               istart = max(slo(1), i_w)
               jstart = max(slo(2), j_s)
               kstart = max(slo(3), k_b)
               iend   = min(shi(1), i_e)
               jend   = min(shi(2), j_n)
               kend   = min(shi(3), k_t)
               vel(istart:iend,jstart:jend,kstart:kend,2) = vgx
               if (slo(2).lt.domlo(2) .and. domlo(2) == jstart) &
                  vel(istart:iend,slo(2):jstart-1,kstart:kend,2) = vgx
               if (shi(2).gt.domhi(2) .and. domhi(2) == jend  ) &
                  vel(istart:iend,jend+1:shi(2)  ,kstart:kend,2) = vgx
            end if

            if (is_defined(wgx)) then
               istart = max(slo(1), i_w)
               jstart = max(slo(2), j_s)
               kstart = max(slo(3), k_b)
               iend   = min(shi(1), i_e)
               jend   = min(shi(2), j_n)
               kend   = min(shi(3), k_t)
               vel(istart:iend,jstart:jend,kstart:kend,3) = wgx
               if (slo(3).lt.domlo(3) .and. domlo(3) == kstart) &
                  vel(istart:iend,jstart:jend,slo(3):kstart-1,3) = wgx
               if (shi(3).gt.domhi(3) .and. domhi(3) == kend  ) &
                  vel(istart:iend,jstart:jend,kend+1:shi(3)  ,3) = wgx
            end if

         endif
      enddo

   end subroutine set_ic

end module init_fluid_module
