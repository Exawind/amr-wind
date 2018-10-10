module rheology_module

   use amrex_fort_module, only : rt => amrex_real
   use iso_c_binding , only: c_int
   use param, only: zero, half, one, two

   implicit none

contains

   !
   ! Compute the viscosity distribution of a 
   ! Power-law fluid: 
   !
   ! eta = mu dot(gamma)^(n-1)
   !
   subroutine powerlaw_viscosity(lo, hi, eta, elo, ehi, sr, slo, shi) &
         bind(C, name="powerlaw_viscosity")

      use constant, only: mu, n
      
      integer(c_int), intent(in   ) :: elo(3), ehi(3)
      integer(c_int), intent(in   ) :: slo(3), shi(3)
      integer(c_int), intent(in   ) ::  lo(3), hi(3)

      real(rt),   intent(  out) :: &
         eta(elo(1):ehi(1),elo(2):ehi(2),elo(3):ehi(3))

      real(rt), intent(in   ) :: &
         sr(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      ! Local variables
      !-----------------------------------------------
      integer      :: i, j, k

      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               
               eta(i,j,k) = mu * sr(i,j,k)**(n - one)

            end do
         end do
      end do

   end subroutine powerlaw_viscosity

   !
   ! Compute the viscosity distribution of a 
   ! Papanastasiou-regularised Bingham fluid: 
   !
   ! eta = (mu dot(gamma) + tau_0) (1 - exp(-dot(gamma) / eps)) / dot(gamma)
   !
   subroutine bingham_viscosity(lo, hi, eta, elo, ehi, sr, slo, shi) &
         bind(C, name="bingham_viscosity")

      use constant, only: mu, tau_0, papa_reg
      
      integer(c_int), intent(in   ) :: elo(3), ehi(3)
      integer(c_int), intent(in   ) :: slo(3), shi(3)
      integer(c_int), intent(in   ) ::  lo(3), hi(3)

      real(rt),   intent(  out) :: &
         eta(elo(1):ehi(1),elo(2):ehi(2),elo(3):ehi(3))

      real(rt), intent(in   ) :: &
         sr(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      ! Local variables
      !-----------------------------------------------
      integer      :: i, j, k
      real(rt) :: nu

      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)

               nu = sr(i,j,k) / papa_reg
               eta(i,j,k) = (mu * sr(i,j,k) + tau_0) * expterm(nu) / papa_reg

            end do
         end do
      end do

   end subroutine bingham_viscosity

   !
   ! Compute the viscosity distribution of a 
   ! Papanastasiou-regularised Herschel-Bulkley fluid: 
   !
   ! eta = (mu dot(gamma)^n + tau_0) (1 - exp(-dot(gamma) / eps)) / dot(gamma)
   !
   subroutine hb_viscosity(lo, hi, eta, elo, ehi, sr, slo, shi) &
         bind(C, name="hb_viscosity")

      use constant, only: mu, n, tau_0, papa_reg
      
      integer(c_int), intent(in   ) :: elo(3), ehi(3)
      integer(c_int), intent(in   ) :: slo(3), shi(3)
      integer(c_int), intent(in   ) ::  lo(3), hi(3)

      real(rt),   intent(  out) :: &
         eta(elo(1):ehi(1),elo(2):ehi(2),elo(3):ehi(3))

      real(rt), intent(in   ) :: &
         sr(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      ! Local variables
      !-----------------------------------------------
      integer  :: i, j, k
      real(rt) :: nu

      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)

               nu = sr(i,j,k) / papa_reg
               eta(i,j,k) = (mu * sr(i,j,k)**n + tau_0) * expterm(nu) / papa_reg

            end do
         end do
      end do

   end subroutine hb_viscosity

   !
   ! Compute the viscosity distribution of a 
   ! de Souza Mendes - Dutra fluid: 
   !
   ! eta = (mu dot(gamma)^n + tau_0) (1 - exp(-eta_0 dot(gamma) / tau_0)) / dot(gamma)
   !
   subroutine smd_viscosity(lo, hi, eta, elo, ehi, sr, slo, shi) &
         bind(C, name="smd_viscosity")

      use constant, only: mu, n, tau_0, eta_0
      
      integer(c_int), intent(in   ) :: elo(3), ehi(3)
      integer(c_int), intent(in   ) :: slo(3), shi(3)
      integer(c_int), intent(in   ) ::  lo(3), hi(3)

      real(rt),   intent(  out) :: &
         eta(elo(1):ehi(1),elo(2):ehi(2),elo(3):ehi(3))

      real(rt), intent(in   ) :: &
         sr(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      ! Local variables
      !-----------------------------------------------
      integer  :: i, j, k
      real(rt) :: nu

      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)

               nu = eta_0 * sr(i,j,k) / tau_0
               eta(i,j,k) = (mu * sr(i,j,k)**n + tau_0) * expterm(nu) * eta_0 / tau_0

            end do
         end do
      end do

   end subroutine smd_viscosity

   real(rt) function expterm(nu)
      real(rt), intent(in) :: nu
      ! Avoid overflow 
      if (nu .lt. 1.0e-14) then 
         expterm = one - half * nu + nu**2 / 6.0d0 - nu**3 / 24.0d0
      else
         expterm = (one - exp(-nu)) / nu
      end if
   end function expterm

end module rheology_module
