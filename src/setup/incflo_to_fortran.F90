module incflo_to_fortran_module
! _________________________________________________________________

   use amrex_fort_module, only : rt => amrex_real
   use iso_c_binding , only: c_int, c_char

   implicit none

contains

!**************************************************************************!
! Take constants from incflo and make them available in the Fortran module !
! "constant", so that they can be accessed from f90 functions.             !
!**************************************************************************!

   subroutine fortran_get_data(is_cyclic_in, &
                               delp_in, gravity_in, ro_0_in, mu_in, &
                               ic_u_in, ic_v_in, ic_w_in, ic_p_in, &
                               n_in, tau_0_in, papa_reg_in, eta_0_in, &
                               fluid_model_name, fluid_model_namelength) &
                           bind(C, name="fortran_get_data")

      use bc, only: cyclic_x, cyclic_y, cyclic_z
      use constant

      implicit none

      logical,                intent(in) :: is_cyclic_in(3)
      real(rt),               intent(in) :: delp_in(3)
      real(rt),               intent(in) :: gravity_in(3)
      real(rt),               intent(in) :: ro_0_in, mu_in
      real(rt),               intent(in) :: ic_u_in, ic_v_in, ic_w_in, ic_p_in
      real(rt),               intent(in) :: n_in, tau_0_in, papa_reg_in, eta_0_in
      character(kind=c_char), intent(in) :: fluid_model_name(*)
      integer(c_int),         intent(in), value :: fluid_model_namelength

      ! Local 
      integer :: i

      cyclic_x = is_cyclic_in(1)
      cyclic_y = is_cyclic_in(2)
      cyclic_z = is_cyclic_in(3)
      delp(:) = delp_in(:)
      gravity(:) = gravity_in(:)
      ro_0 = ro_0_in
      mu = mu_in
      ic_u = ic_u_in
      ic_v = ic_v_in
      ic_w = ic_w_in
      ic_p = ic_p_in
      n = n_in
      tau_0 = tau_0_in
      papa_reg = papa_reg_in
      eta_0 = eta_0_in
      
      allocate(character(fluid_model_namelength) :: fluid_model)
      forall(i = 1:fluid_model_namelength) fluid_model(i:i) = fluid_model_name(i)

   end subroutine fortran_get_data

end module incflo_to_fortran_module
