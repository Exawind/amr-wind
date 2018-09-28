!
!
!  This module contains the subroutines to perform some of the steps of the
!  projection method.
!
!
module projection_mod

   use amrex_fort_module, only: ar => amrex_real
   use iso_c_binding ,    only: c_int
   use param,             only: zero, half, one

   implicit none
   private

   ! Define here the unit vectors
   ! This is used to shift index  based on how the variable is staggered
   ! Check e_x, e_y and e_z in incflo_level.H
   integer(c_int), parameter :: e_i(3,3) = reshape ( [1,0,0,0,1,0,0,0,1], [3,3] )

contains

   !
   ! Compute the coefficients of the PPE, i.e. 1 / ro,
   ! at the faces of the pressure cells along the "dir"-axis.
   !
   subroutine compute_bcoeff_cc ( lo, hi, bcoeff, blo, bhi, &
        ro, slo, shi, dir )  bind(C)

      ! Loop bounds
      integer(c_int), intent(in   ) ::  lo(3), hi(3)

      ! Array bounds
      integer(c_int), intent(in   ) :: slo(3),shi(3)
      integer(c_int), intent(in   ) :: blo(3),bhi(3)

      ! Direction
      integer(c_int), intent(in   ) :: dir

      ! Arrays
      real(ar),       intent(in   ) :: &
           ro(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(ar),       intent(  out) :: &
           bcoeff(blo(1):bhi(1),blo(2):bhi(2),blo(3):bhi(3))

      integer      :: i, j, k, i0, j0, k0

      i0 = e_i(dir,1)
      j0 = e_i(dir,2)
      k0 = e_i(dir,3)

      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               bcoeff(i,j,k) = half * ( one/ro(i,j,k) + one/ro(i-i0,j-j0,k-k0) )
            end do
         end do
      end do

   end subroutine compute_bcoeff_cc

   subroutine compute_bcoeff_nd ( lo, hi, bcoeff, blo, bhi, &
        ro, slo, shi, dir )  bind(C)

      ! Loop bounds
      integer(c_int), intent(in   ) ::  lo(3), hi(3)

      ! Array bounds
      integer(c_int), intent(in   ) :: slo(3),shi(3)
      integer(c_int), intent(in   ) :: blo(3),bhi(3)

      ! Direction
      integer(c_int), intent(in   ) :: dir

      ! Arrays
      real(ar),       intent(in   ) :: &
           ro(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(ar),       intent(  out) :: &
           bcoeff(blo(1):bhi(1),blo(2):bhi(2),blo(3):bhi(3))

      integer      :: i, j, k

      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               bcoeff(i,j,k) =  one / ro(i,j,k)
            end do
         end do
      end do

   end subroutine compute_bcoeff_nd

   !
   ! Set the boundary condition for Pressure Poisson Equation (PPE)
   !
   ! MLMG expects the BC type to be the uniform on each domain wall.
   ! Since incflo allows for BC patches on each wall, we first check that
   ! the user-provided BCs are uniform, and then return a single BC type for
   ! each domain wall.
   !
   subroutine set_ppe_bc ( bc_lo, bc_hi, domlo, domhi, ng, bct_ilo, bct_ihi, &
        & bct_jlo, bct_jhi, bct_klo, bct_khi)  bind(C)

      use amrex_lo_bctypes_module
      use bc

      ! Array of global BC types
      integer(c_int), intent(  out) :: bc_lo(3), bc_hi(3)

      ! Domain bounds
      integer(c_int), intent(in   ) :: domlo(3), domhi(3), ng

      ! Arrays of point-by-point BC types
      integer(c_int), intent(in   )  ::                                 &
           & bct_ilo(domlo(2)-ng:domhi(2)+ng,domlo(3)-ng:domhi(3)+ng,2), &
           & bct_ihi(domlo(2)-ng:domhi(2)+ng,domlo(3)-ng:domhi(3)+ng,2), &
           & bct_jlo(domlo(1)-ng:domhi(1)+ng,domlo(3)-ng:domhi(3)+ng,2), &
           & bct_jhi(domlo(1)-ng:domhi(1)+ng,domlo(3)-ng:domhi(3)+ng,2), &
           & bct_klo(domlo(1)-ng:domhi(1)+ng,domlo(2)-ng:domhi(2)+ng,2), &
           & bct_khi(domlo(1)-ng:domhi(1)+ng,domlo(2)-ng:domhi(2)+ng,2)

      ! Local variables
      integer(c_int)                :: bc_face

      !
      ! By default, all the BCs are Neumann
      !
      bc_lo    = amrex_lo_neumann
      bc_hi    = amrex_lo_neumann

      !
      ! BC -- X direction
      !
      if ( cyclic_x ) then
         bc_lo(1) = amrex_lo_periodic
         bc_hi(1) = amrex_lo_periodic
      else

         ! X at domlo(1)
         bc_face = get_bc_face(bct_ilo,ng)
         if ( (bc_face == pinf_) .or. (bc_face == pout_) ) then
            bc_lo(1) = amrex_lo_dirichlet
         end if

         ! X at domhi(1)
         bc_face = get_bc_face(bct_ihi,ng)
         if ( (bc_face == pinf_) .or. (bc_face == pout_) ) then
            bc_hi(1) = amrex_lo_dirichlet
         end if

      end if


      !
      ! BC -- Y direction
      !
      if ( cyclic_y ) then
         bc_lo(2) = amrex_lo_periodic
         bc_hi(2) = amrex_lo_periodic
      else

         ! Y at domlo(2)
         bc_face = get_bc_face(bct_jlo,ng)
         if ( (bc_face == pinf_) .or. (bc_face == pout_) ) then
            bc_lo(2) = amrex_lo_dirichlet
         end if

         ! Y at domhi(2)
         bc_face = get_bc_face(bct_jhi,ng)
         if ( (bc_face == pinf_) .or. (bc_face == pout_) ) then
            bc_hi(2) = amrex_lo_dirichlet
         end if

      end if

      !
      ! BC -- Z direction
      !
      if ( cyclic_z ) then
         bc_lo(3) = amrex_lo_periodic
         bc_hi(3) = amrex_lo_periodic
      else

         ! Z at domlo(3)
         bc_face = get_bc_face(bct_klo,ng)
         if ( (bc_face == pinf_) .or. (bc_face == pout_) ) then
            bc_lo(3) = amrex_lo_dirichlet
         end if

         ! Z at domhi(3)
         bc_face = get_bc_face(bct_khi,ng)
         if ( (bc_face == pinf_) .or. (bc_face == pout_) ) then
            bc_hi(3) = amrex_lo_dirichlet
         end if

      end if

   contains

      !
      ! Test whether the BC type is the same everywhere on
      ! the face. If BC is uniform on face, it returns its value
      !
      function get_bc_face (bct_array,nghost) result (bc_face)
         integer(c_int), intent(in   ) :: bct_array(:,:,:)
         integer(c_int), intent(in   ) :: nghost
         integer                       :: bc_face
         integer                       :: is, ie, js, je

         ! Do not consider the edges: they may cause problems
         is = nghost+1
         ie = size(bct_array,1) - nghost
         js = nghost+1
         je = size(bct_array,2) - nghost

         bc_face = bct_array(is,js,1)

         if ( .not. all (bct_array(is:ie,js:je,1) == bc_face) ) then
            stop "BC type must be uniform on each face of the domain"
         end if

      end function get_bc_face

   end subroutine set_ppe_bc

end module projection_mod
