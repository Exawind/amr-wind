module rheology_module

   use amrex_error_module, only: amrex_abort
   use amrex_fort_module, only : rt => amrex_real
   use constant, only: half, one
   use constant, only: mu, n, tau_0, eta_0, papa_reg, fluid_model

   implicit none

contains

   real(rt) function viscosity(sr) & 
         bind(C, name="viscosity")

      real(rt), value :: sr
      real(rt)        :: nu

      if (fluid_model == "newtonian") then 

         ! Viscosity is constant
         viscosity = mu

      else if (fluid_model == "powerlaw") then 

         ! Power-law fluid: 
         !
         ! eta = mu dot(gamma)^(n-1)
         viscosity = mu * sr**(n - one)

      else if (fluid_model == "bingham") then 

         ! Papanastasiou-regularised Bingham fluid: 
         !
         ! eta = mu + tau_0 (1 - exp(-dot(gamma) / eps)) / dot(gamma)
         nu = sr / papa_reg
         viscosity = mu + tau_0 * expterm(nu) / papa_reg

      else if (fluid_model == "hb") then 

         ! Papanastasiou-regularised Herschel-Bulkley fluid: 
         !
         ! eta = (mu dot(gamma)^n + tau_0) (1 - exp(-dot(gamma) / eps)) / dot(gamma)
         
         nu = sr / papa_reg
         viscosity = (mu * sr**n + tau_0) * expterm(nu) / papa_reg

      else if (fluid_model == "smd") then

         ! de Souza Mendes - Dutra fluid: 
         !
         ! eta = (mu dot(gamma)^n + tau_0) (1 - exp(-eta_0 dot(gamma) / tau_0)) / dot(gamma)
         
         nu = eta_0 * sr / tau_0
         viscosity = (mu * sr**n + tau_0) * expterm(nu) * eta_0 / tau_0

      else

         ! This should have been caught earlier, but doesn't hurt to double check
         call amrex_abort("Unknown fluid_model! Choose either newtonian, powerlaw, bingham, hb, smd")

      end if

   end function viscosity


   ! 
   ! Compute the exponential term:
   !
   !  ( 1 - exp(-nu) ) / nu ,
   !
   ! making sure to avoid overflow for small nu by using the exponential Taylor series
   !
   real(rt) function expterm(nu)
      real(rt), intent(in) :: nu
      ! Avoid overflow 
      if (nu .lt. 1.0e-9) then 
         expterm = one - half * nu + nu**2 / 6.0d0 - nu**3 / 24.0d0
      else
         expterm = (one - exp(-nu)) / nu
      end if
   end function expterm

end module rheology_module
