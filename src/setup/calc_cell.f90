MODULE CALC_CELL_MODULE

   use amrex_fort_module, only : rt => amrex_real
   use iso_c_binding , only: c_int

   private

   public :: calc_cell_ic
   public :: calc_cell_bc_wall, calc_cell_bc_flow

contains
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine name: CALC_CELL                                          !
!  Purpose: calculate the i, j or k cell index for the corresponding   !
!     x y or z reactor location. the index returned depends on which   !
!     half of the i, j or k cell that the x, y, or z position          !
!     intersects                                                       !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   pure integer function calc_cell(location, dx)

      use param, only : half

      implicit none

      ! the x, y or z location
      real(rt), intent(in) :: location

      ! the cell lengths along the corresponding axis (dx, dy or dz)
      real(rt), intent(in) :: dx

      calc_cell = floor(location/dx + half) - 1

   end function calc_cell

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: calc_cell_ic                                            !
!  Purpose: calculate the i, j or k cell index for IC regions.         !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   subroutine calc_cell_ic(dx, dy, dz, x_w, y_s, z_b, x_e, y_n, z_t, &
                           i_w, i_e, j_s, j_n, k_b, k_t)

      use param, only: equal

      implicit none

      real(rt), intent(in   ) :: dx, dy, dz
      real(rt), intent(in   ) :: x_w, y_s, z_b, x_e, y_n, z_t

      integer,      intent(  out) :: i_w, j_s, k_b, i_e, j_n, k_t

      i_w = calc_cell(x_w, dx) + 1
      i_e = calc_cell(x_e, dx)

      j_s = calc_cell(y_s, dy) + 1
      j_n = calc_cell(y_n, dy)

      k_b = calc_cell(z_b, dz) + 1
      k_t = calc_cell(z_t, dz)

   end subroutine calc_cell_ic

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: calc_cell_bc_wall                                       !
!  Purpose: calculate the i, j or k cell index for wall BCs.           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   subroutine calc_cell_bc_wall(domlo, domhi, xlength, ylength, &
                                zlength, dx, dy, dz, x_w, y_s, z_b, x_e, y_n, z_t, &
                                i_w, i_e, j_s, j_n, k_b, k_t)

      use param, only: zero, equal

      implicit none

      integer,      intent(in   ) :: domlo(3), domhi(3)
      real(rt), intent(in   ) :: dx, dy, dz
      real(rt), intent(in   ) :: xlength, ylength, zlength
      real(rt), intent(in   ) :: x_w, y_s, z_b, x_e, y_n, z_t

      integer,      intent(  out) :: i_w, j_s, k_b, i_e, j_n, k_t

      i_w = calc_cell (x_w, dx) + 1
      i_e = calc_cell (x_e, dx)
      if(equal(x_w, x_e)) then
         if(equal(x_w,0.0d0)) then
            i_w = domlo(1)-1
            i_e = domlo(1)-1
         elseif(equal(x_w,xlength)) then
            i_w = domhi(1)+1
            i_e = domhi(1)+1
         endif
      endif

      j_s = calc_cell (y_s, dy) + 1
      j_n = calc_cell (y_n, dy)
      if(equal(y_s, y_n)) then
         if(equal(y_s,zero)) then
            j_s = domlo(2)-1
            j_n = domlo(2)-1
         else if (equal(y_s,ylength)) then
            j_s = domhi(2)+1
            j_n = domhi(2)+1
         endif
      endif

      k_b = calc_cell (z_b, dz) + 1
      k_t = calc_cell (z_t, dz)
      if(equal(z_b, z_t)) then
         if(equal(z_b,zero)) then
            k_b = domlo(3)-1
            k_t = domlo(3)-1
         elseif(equal(z_b,zlength)) then
            k_b = domhi(3)+1
            k_t = domhi(3)+1
         endif
      endif

   end subroutine calc_cell_bc_wall

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: calc_cell_bc_wall                                       !
!  Purpose: calculate the i, j or k cell index for wall BCs.           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   subroutine calc_cell_bc_flow(xlength, ylength, zlength, &
                                dx, dy, dz, x_w, y_s, z_b, x_e, y_n, z_t, &
                                i_w, i_e, j_s, j_n, k_b, k_t)

      use param, only: equal

      implicit none

      real(rt), intent(in   ) :: dx, dy, dz
      real(rt), intent(in   ) :: xlength, ylength, zlength
      real(rt), intent(in   ) :: x_w, y_s, z_b, x_e, y_n, z_t

      integer,      intent(  out) :: i_w, j_s, k_b, i_e, j_n, k_t

      i_w = calc_cell (x_w, dx)
      i_e = calc_cell (x_e, dx)
      if (.not.equal(x_w, x_e)) then
         i_w = i_w + 1
      else if(equal(x_w,xlength)) then
         i_w = i_w + 1
         i_e = i_w
      endif

      j_s = calc_cell (y_s, dy)
      j_n = calc_cell (y_n, dy)
      if(.not.equal(y_s, y_n)) then
         j_s = j_s + 1
      else if(equal(y_s,ylength)) then
         j_s = j_s + 1
         j_n = j_s
      endif

      k_b = calc_cell (z_b, dz)
      k_t = calc_cell (z_t, dz)
      if(.not.equal(z_b, z_t)) then
         k_b = k_b + 1
      else if(equal(z_b,zlength)) then
         k_b = k_b + 1
         k_t = k_b
      endif

   end subroutine calc_cell_bc_flow

END MODULE CALC_CELL_MODULE
