module rheology_module

   use amrex_fort_module, only : rt => amrex_real
   use constant, only: mu_0, tau_0, papa_reg
   use iso_c_binding , only: c_int
   use param, only: zero, half, one, two

   implicit none

contains

   !
   ! Compute the viscosity distribution
   !
   subroutine compute_viscosity(lo, hi, mu, mulo, muhi, strainrate, slo, shi, dx) &
         bind(C, name="compute_viscosity")
      
      integer(c_int), intent(in   ) :: mulo(3), muhi(3)
      integer(c_int), intent(in   ) :: slo(3), shi(3)
      integer(c_int), intent(in   ) ::  lo(3), hi(3)

      real(rt),   intent(in   ) :: dx(3)
      real(rt),   intent(  out) :: &
         mu(mulo(1):muhi(1),mulo(2):muhi(2),mulo(3):muhi(3))

      real(rt), intent(in   ) :: &
         strainrate(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(rt) :: nu

      ! Local variables
      !-----------------------------------------------
      integer      :: i, j, k

      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)

               ! Newtonian part
               mu(i,j,k) = mu_0

               ! Non-zero yield stress
               if (tau_0 .gt. 0.0d0) then 

                  nu = strainrate(i,j,k) / papa_reg

                  ! Bingham model
                  if (nu .lt. 1.0e-14) then 
                     ! Avoid overflow
                     mu(i,j,k) = mu(i,j,k) + tau_0 * & 
                        (one - half * nu + nu**2 / 6.0d0 - nu**3 / 24.0d0) / papa_reg
                  else
                     mu(i,j,k) = mu(i,j,k) + tau_0 * &
                        (one - exp(-nu)) / strainrate(i,j,k)
                  end if

               end if

            end do
         end do
      end do

   end subroutine compute_viscosity

end module rheology_module
