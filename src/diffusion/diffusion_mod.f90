!
!
!  This module contains the subroutines to compute the three components
!  of the diffusion term div(tau) where
!
!      tau = mu ( grad(u) + grad(u)^T )
!
!  Author: Michele Rosso
!
!  Date: October 12, 2017
!
!
module diffusion_mod

   use amrex_fort_module, only: rt => amrex_real
   use iso_c_binding ,    only: c_int
   use param,             only: zero, half, one, two
   use amrex_mempool_module, only: amrex_allocate, amrex_deallocate

   implicit none
   private

   real(rt), parameter :: q4 = one / ( two * two )

   public compute_divtau
   public fill_vel_diff_bc

   ! Define here the unit vectors
   ! This is used to shift index  based on how the variable is staggered
   ! Check e_x, e_y and e_z in incflo_level.H
   integer(c_int), parameter :: e_i(3,3) = reshape ( [1,0,0,0,1,0,0,0,1], [3,3] )

contains
   !
   ! Computes  d(txx)/dx + d(txy)/dy + d(txz)/dz
   !
   !  txx = 2 * mu * du/dx
   !  txy =  mu * ( du/dy + dv/dx )
   !  txz =  mu * ( du/dz + dw/dx )
   !
   subroutine compute_divtau (lo, hi, &
        divtau, dlo, dhi, &
        vel_in, vinlo, vinhi, &
        mu, lambda, ro, &
        slo, shi, &
        domlo, domhi, &
        bc_ilo_type, bc_ihi_type, &
        bc_jlo_type, bc_jhi_type, &
        bc_klo_type, bc_khi_type, &
        dx, ng, &
        do_explicit_diffusion) bind(C)


    ! Loops bounds (cell-centered)
    integer(c_int),  intent(in   ) :: lo(3),  hi(3)

    ! Number of ghost cells
    integer(c_int),  intent(in   ) :: ng

    ! If true  then we include all the diffusive terms in this explicit result
    ! If false then we include all only the off-diagonal terms here -- we do this
    !     by computing the full tensor then subtracting the diagonal terms
    integer(c_int),  intent(in   ) :: do_explicit_diffusion

    ! Array bounds
    integer(c_int),  intent(in   ) :: vinlo(3), vinhi(3)
    integer(c_int),  intent(in   ) :: slo(3), shi(3)
    integer(c_int),  intent(in   ) :: dlo(3), dhi(3)
    integer(c_int),  intent(in   ) :: domlo(3), domhi(3)

    ! Grid
    real(rt),        intent(in   ) :: dx(3)

    ! Arrays
    real(rt),        intent(in   ) ::                           &
         & vel_in(vinlo(1):vinhi(1),vinlo(2):vinhi(2),vinlo(3):vinhi(3),3), &
         & ro(      slo(1):  shi(1),  slo(2):  shi(2),  slo(3):  shi(3)),  &
         & mu(      slo(1):  shi(1),  slo(2):  shi(2),  slo(3):  shi(3)),  &
         & lambda(  slo(1):  shi(1),  slo(2):  shi(2),  slo(3):  shi(3))

    real(rt),        intent(inout) ::                        &
         & divtau(  dlo(1):  dhi(1),  dlo(2):  dhi(2),  dlo(3):  dhi(3),3)

    ! BC types
    integer(c_int), intent(in   ) ::  &
         & bc_ilo_type(domlo(2)-ng:domhi(2)+ng,domlo(3)-ng:domhi(3)+ng,2), &
         & bc_ihi_type(domlo(2)-ng:domhi(2)+ng,domlo(3)-ng:domhi(3)+ng,2), &
         & bc_jlo_type(domlo(1)-ng:domhi(1)+ng,domlo(3)-ng:domhi(3)+ng,2), &
         & bc_jhi_type(domlo(1)-ng:domhi(1)+ng,domlo(3)-ng:domhi(3)+ng,2), &
         & bc_klo_type(domlo(1)-ng:domhi(1)+ng,domlo(2)-ng:domhi(2)+ng,2), &
         & bc_khi_type(domlo(1)-ng:domhi(1)+ng,domlo(2)-ng:domhi(2)+ng,2)


    ! Temporary array just to handle bc's
    integer(c_int) :: vlo(3), vhi(3)
    real(rt), dimension(:,:,:,:), pointer, contiguous :: vel

    ! Temporaty array to handle div(u) at the nodes
    real(rt), dimension(:,:,:  ), pointer, contiguous :: divu

    integer(c_int)                 :: i, j, k, n
    real(rt)                       :: idx, idy, idz
    real(rt)                       :: du, dv, dw

    real(rt)  :: txx, tyy, tzz
    real(rt)  :: mu_e, mu_w, mu_n, mu_s, mu_t, mu_b
    real(rt)  :: lambda_e, lambda_w, lambda_n, lambda_s, lambda_t, lambda_b
    real(rt)  :: divu_e, divu_w, divu_n, divu_s, divu_t, divu_b
    real(rt)  :: txx_e, txx_w, txy_n, txy_s, txz_t, txz_b
    real(rt)  :: txy_e, txy_w, tyy_n, tyy_s, tyz_t, tyz_b
    real(rt)  :: txz_e, txz_w, tyz_n, tyz_s, tzz_t, tzz_b

    idx = one / dx(1)
    idy = one / dx(2)
    idz = one / dx(3)

    vlo = lo - ng
    vhi = hi + ng
    call amrex_allocate( vel, vlo(1), vhi(1)  , vlo(2), vhi(2)  , vlo(3), vhi(3)  , 1, 3)
    call amrex_allocate(divu, vlo(1), vhi(1)+1, vlo(2), vhi(2)+1, vlo(3), vhi(3)+1)

    ! Put values into ghost cells so we can easy take derivatives
    call fill_vel_diff_bc(vel_in, vinlo, vinhi, vel, lo, hi, domlo, domhi, ng, &
                          bc_ilo_type, bc_ihi_type, &
                          bc_jlo_type, bc_jhi_type, &
                          bc_klo_type, bc_khi_type)

    !
    ! Compute div(u) at the nodes
    !
    do k = lo(3), hi(3)+1
       do j = lo(2), hi(2)+1
          do i = lo(1), hi(1)+1

             ! Divergence
             du = (   vel(i  ,j  ,k  ,1) + vel(i  ,j-1,k  ,1) &
                  &  +vel(i  ,j  ,k-1,1) + vel(i  ,j-1,k-1,1) &
                  &  -vel(i-1,j  ,k  ,1) - vel(i-1,j-1,k  ,1) &
                  &  -vel(i-1,j  ,k-1,1) - vel(i-1,j-1,k-1,1) )

             dv = (   vel(i  ,j  ,k  ,2) + vel(i-1,j  ,k  ,2) &
                  &  +vel(i  ,j  ,k-1,2) + vel(i-1,j  ,k-1,2) &
                  &  -vel(i  ,j-1,k  ,2) - vel(i-1,j-1,k  ,2) &
                  &  -vel(i  ,j-1,k-1,2) - vel(i-1,j-1,k-1,2) )

             dw = (   vel(i  ,j  ,k  ,3) + vel(i-1,j  ,k  ,3) &
                  &  +vel(i  ,j-1,k  ,3) + vel(i-1,j-1,k  ,3) &
                  &  -vel(i  ,j  ,k-1,3) - vel(i-1,j  ,k-1,3) &
                  &  -vel(i  ,j-1,k-1,3) - vel(i-1,j-1,k-1,3) )

             divu(i,j,k) = ( du*idx + dv*idy + dw*idz ) * q4

          end do
       end do
    end do



    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             mu_w = half * (mu(i,j,k) + mu(i-1,j,k))
             mu_e = half * (mu(i,j,k) + mu(i+1,j,k))
             mu_s = half * (mu(i,j,k) + mu(i,j-1,k))
             mu_n = half * (mu(i,j,k) + mu(i,j+1,k))
             mu_b = half * (mu(i,j,k) + mu(i,j,k-1))
             mu_t = half * (mu(i,j,k) + mu(i,j,k+1))

             lambda_w = half * (lambda(i,j,k) + lambda(i-1,j,k))
             lambda_e = half * (lambda(i,j,k) + lambda(i+1,j,k))
             lambda_s = half * (lambda(i,j,k) + lambda(i,j-1,k))
             lambda_n = half * (lambda(i,j,k) + lambda(i,j+1,k))
             lambda_b = half * (lambda(i,j,k) + lambda(i,j,k-1))
             lambda_t = half * (lambda(i,j,k) + lambda(i,j,k+1))

             !*************************************
             !         div(tau)_x
             !*************************************

             ! X
             txx_e = two * mu_e * ( vel(i+1,j,k,1) - vel(i  ,j,k,1) ) * idx
             txx_w = two * mu_w * ( vel(i  ,j,k,1) - vel(i-1,j,k,1) ) * idx

             ! Y north
             du = vel(i,j+1,k,1) - vel(i,j,k,1)

             dv = (  vel(i+1,j,k,2) + vel(i+1,j+1,k,2) + vel(i  ,j,k,2) + vel(i  ,j+1,k,2) &
                  &- vel(i  ,j,k,2) - vel(i  ,j+1,k,2) - vel(i-1,j,k,2) - vel(i-1,j+1,k,2) ) * q4

             txy_n = mu_n * ( du*idy + dv*idx )

             ! Y south
             du = vel(i,j,k,1) - vel(i,j-1,k,1)

             dv = (  vel(i+1,j-1,k,2) + vel(i+1,j,k,2) + vel(i  ,j-1,k,2) + vel(i  ,j,k,2) &
                  &- vel(i  ,j-1,k,2) - vel(i  ,j,k,2) - vel(i-1,j-1,k,2) - vel(i-1,j,k,2) ) * q4

             txy_s = mu_s * ( du*idy + dv*idx )

             ! Z top
             du = vel(i,j,k+1,1) - vel(i,j,k,1)

             dw = (  vel(i+1,j,k,3) + vel(i+1,j,k+1,3) + vel(i  ,j,k,3) + vel(i  ,j,k+1,3) &
                  &- vel(i  ,j,k,3) - vel(i  ,j,k+1,3) - vel(i-1,j,k,3) - vel(i-1,j,k+1,3) ) * q4

             txz_t = mu_t * ( du*idz + dw*idx )

             ! Z bottom
             du = vel(i,j,k,1) - vel(i,j,k-1,1)

             dw = (  vel(i+1,j,k-1,3) + vel(i+1,j,k,3) + vel(i  ,j,k-1,3) + vel(i  ,j,k,3) &
                  &- vel(i  ,j,k-1,3) - vel(i  ,j,k,3) - vel(i-1,j,k-1,3) - vel(i-1,j,k,3) ) * q4

             txz_b = mu_b * ( du*idz + dw*idx )


             ! Div term
             divu_e = lambda_e * ( divu(i+1,j,k) + divu(i+1,j+1,k) + divu(i+1,j,k+1) + divu(i+1,j+1,k+1)) * q4
             divu_w = lambda_w * ( divu(i  ,j,k) + divu(i  ,j+1,k) + divu(i  ,j,k+1) + divu(i  ,j+1,k+1)) * q4


             ! Assemble
             divtau(i,j,k,1) = ( txx_e - txx_w ) * idx  + &
                  &            ( txy_n - txy_s ) * idy  + &
                  &            ( txz_t - txz_b ) * idz  + &
                  &            ( divu_e - divu_w ) * idx

             !*************************************
             !         div(tau)_y
             !*************************************
             ! X east
             du = (   vel(i+1,j  ,k,1) + vel(i+1,j+1,k,1) + vel(i,j  ,k,1) + vel(i,j+1,k,1) &
                  & - vel(i+1,j-1,k,1) - vel(i+1,j  ,k,1) - vel(i,j-1,k,1) - vel(i,j  ,k,1) ) * q4

             dv = vel(i+1,j,k,2) - vel(i,j,k,2)

             txy_e = mu_e * ( du*idy + dv*idx )

             ! X west
             du = (   vel(i,j  ,k,1) + vel(i,j+1,k,1) + vel(i-1,j  ,k,1) + vel(i-1,j+1,k,1) &
                  & - vel(i,j-1,k,1) - vel(i,j  ,k,1) - vel(i-1,j-1,k,1) - vel(i-1,j  ,k,1) ) * q4

             dv = vel(i,j,k,2) - vel(i-1,j,k,2)

             txy_w = mu_w * ( du*idy + dv*idx )

             ! Y
             tyy_n = two * mu_n * ( vel(i,j+1,k,2) - vel(i,j  ,k,2) ) * idy
             tyy_s = two * mu_s * ( vel(i,j  ,k,2) - vel(i,j-1,k,2) ) * idy

             ! Z top
             dv = vel(i,j,k+1,2) - vel(i,j,k,2)

             dw = (   vel(i,j  ,k+1,3) + vel(i,j+1,k+1,3) + vel(i,j  ,k,3) + vel(i,j+1,k,3) &
                  & - vel(i,j-1,k+1,3) - vel(i,j  ,k+1,3) - vel(i,j-1,k,3) - vel(i,j  ,k,3) ) * q4

             tyz_t = mu_t * ( dv*idz + dw*idy )

             ! Z bottom
             dv = vel(i,j,k,2) - vel(i,j,k-1,2)

             dw = (   vel(i,j  ,k,3) + vel(i,j+1,k,3) + vel(i,j  ,k-1,3) + vel(i,j+1,k-1,3) &
                  & - vel(i,j-1,k,3) - vel(i,j  ,k,3) - vel(i,j-1,k-1,3) - vel(i,j  ,k-1,3) ) * q4

             tyz_b = mu_b * ( dv*idz + dw*idy )


             ! Div term
             divu_n = lambda_n * ( divu(i,j+1,k) + divu(i,j+1,k+1) + divu(i+1,j+1,k+1) + divu(i+1,j+1,k) ) * q4
             divu_s = lambda_s * ( divu(i,j  ,k) + divu(i,j  ,k+1) + divu(i+1,j  ,k+1) + divu(i+1,j  ,k) ) * q4

             ! Assemble
             divtau(i,j,k,2) = ( txy_e - txy_w ) * idx  + &
                  &            ( tyy_n - tyy_s ) * idy  + &
                  &            ( tyz_t - tyz_b ) * idz  + &
                  &            ( divu_n - divu_s ) * idy


             !*************************************
             !         div(tau)_z
             !*************************************

             ! X east
             dw = vel(i+1,j,k,3) - vel(i,j,k,3)

             du = (   vel(i+1,j,k  ,1) + vel(i+1,j,k+1,1) + vel(i,j,k  ,1) + vel(i,j,k+1,1) &
                  & - vel(i+1,j,k-1,1) - vel(i+1,j,k  ,1) - vel(i,j,k-1,1) - vel(i,j,k  ,1) ) * q4

             txz_e = mu_e * ( du*idz + dw*idx )

             ! X west
             dw = vel(i,j,k,3) - vel(i-1,j,k,3)

             du = (   vel(i,j,k  ,1) + vel(i,j,k+1,1) + vel(i-1,j,k  ,1) + vel(i-1,j,k+1,1) &
                  & - vel(i,j,k-1,1) - vel(i,j,k  ,1) - vel(i-1,j,k-1,1) - vel(i-1,j,k  ,1) ) * q4

             txz_w = mu_w * ( du*idz + dw*idx )

             ! Y north
             dw = vel(i,j+1,k,3) - vel(i,j,k,3)

             dv = (   vel(i,j,k+1,2) + vel(i,j+1,k+1,2) + vel(i,j,k  ,2) + vel(i,j+1,k  ,2) &
                  & - vel(i,j,k  ,2) - vel(i,j+1,k  ,2) - vel(i,j,k-1,2) - vel(i,j+1,k-1,2) ) * q4

             tyz_n = mu_n * ( dv*idz + dw*idy )

             ! Y south
             dw = vel(i,j,k,3) - vel(i,j-1,k,3)

             dv = (   vel(i,j-1,k+1,2) + vel(i,j,k+1,2) + vel(i,j-1,k  ,2) + vel(i,j,k  ,2) &
                  & - vel(i,j-1,k  ,2) - vel(i,j,k  ,2) - vel(i,j-1,k-1,2) - vel(i,j,k-1,2) ) * q4

             tyz_s = mu_s * ( dv*idz + dw*idy )

             ! Z
             tzz_t = two * mu_t * ( vel(i,j,k+1,3) - vel(i,j,k  ,3) ) * idz
             tzz_b = two * mu_b * ( vel(i,j,k  ,3) - vel(i,j,k-1,3) ) * idz

             ! Div term
             divu_t = lambda_t * ( divu(i,j,k+1) + divu(i+1,j,k+1) + divu(i+1,j+1,k+1) + divu(i,j+1,k+1) ) * q4
             divu_b = lambda_b * ( divu(i,j,k  ) + divu(i+1,j,k  ) + divu(i+1,j+1,k  ) + divu(i,j+1,k  ) ) * q4


             ! Assemble
             divtau(i,j,k,3) = ( txz_e - txz_w ) * idx  + &
                  &            ( tyz_n - tyz_s ) * idy  + &
                  &            ( tzz_t - tzz_b ) * idz  + &
                  &            ( divu_t - divu_b ) * idz

             if (do_explicit_diffusion .eq. 0) then
                !
                ! Subtract diagonal terms of stress tensor, to be obtained through implicit solve
                ! Note that the variable names are misleading, but we want to avoid declaring new ones
                !
                do n = 1, 3
                   txx = ( mu_e * ( vel(i+1,j,k,n) - vel(i  ,j,k,n) ) &
                          -mu_w * ( vel(i  ,j,k,n) - vel(i-1,j,k,n) ) ) * idx * idx
                   tyy = ( mu_n * ( vel(i,j+1,k,n) - vel(i,j  ,k,n) ) &
                          -mu_s * ( vel(i,j  ,k,n) - vel(i,j-1,k,n) ) ) * idy * idy
                   tzz = ( mu_t * ( vel(i,j,k+1,n) - vel(i,j,k  ,n) ) &
                          -mu_b * ( vel(i,j,k  ,n) - vel(i,j,k-1,n) ) ) * idz * idz
                   divtau(i,j,k,n) = divtau(i,j,k,n) - (txx + tyy + tzz)
                end do
             end if

             !*************************************
             !         div(tau)/ro
             !*************************************
             divtau(i,j,k,:) = divtau(i,j,k,:) / ro(i,j,k)

          end do
       end do
    end do

    call amrex_deallocate(vel)
    call amrex_deallocate(divu)

 end subroutine compute_divtau

   !
   ! Compute the coefficients for the diffusion solve
   ! at the faces of the cells along the "dir"-axis.
   !
   subroutine compute_bcoeff_diff ( lo, hi, bcoeff, blo, bhi, &
        mu, slo, shi, dir )  bind(C)

      ! Loop bounds
      integer(c_int), intent(in   ) ::  lo(3), hi(3)

      ! Array bounds
      integer(c_int), intent(in   ) :: slo(3),shi(3)
      integer(c_int), intent(in   ) :: blo(3),bhi(3)

      ! Direction
      integer(c_int), intent(in   ) :: dir

      ! Arrays
      real(rt),       intent(in   ) :: &
           mu(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(rt),       intent(  out) :: &
           bcoeff(blo(1):bhi(1),blo(2):bhi(2),blo(3):bhi(3))

      integer      :: i, j, k, i0, j0, k0

      i0 = e_i(dir,1)
      j0 = e_i(dir,2)
      k0 = e_i(dir,3)

      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               bcoeff(i,j,k) = half * (mu(i,j,k) + mu(i-i0,j-j0,k-k0))
            end do
         end do
      end do

   end subroutine compute_bcoeff_diff

   !
   ! Set the boundary condition for diffusion solve
   !
   ! MLMG expects the BC type to be the uniform on each domain wall.
   ! Since incflo allows for BC patches on each wall, we first check that
   ! the user-provided BCs are uniform, and then return a single BC type for
   ! each domain wall.
   !
   subroutine set_diff_bc ( bc_lo, bc_hi, domlo, domhi, ng, bct_ilo, bct_ihi, &
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
      ! By default, all the BCs are Dirichlet
      !
      bc_lo    = amrex_lo_dirichlet
      bc_hi    = amrex_lo_dirichlet

      !
      ! BC -- X direction
      !
      if ( cyclic_x ) then
         bc_lo(1) = amrex_lo_periodic
         bc_hi(1) = amrex_lo_periodic
      else

         ! X at domlo(1)
         bc_face = get_bc_face(bct_ilo, ng)
         if ( (bc_face == pinf_) .or. (bc_face == pout_) .or. (bc_face == fsw_) ) then
            bc_lo(1) = amrex_lo_neumann
         end if

         ! X at domhi(1)
         bc_face = get_bc_face(bct_ihi, ng)
         if ( (bc_face == pinf_) .or. (bc_face == pout_) .or. (bc_face == fsw_) ) then
            bc_hi(1) = amrex_lo_neumann
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
         bc_face = get_bc_face(bct_jlo, ng)
         if ( (bc_face == pinf_) .or. (bc_face == pout_) .or. (bc_face == fsw_) ) then
            bc_lo(2) = amrex_lo_neumann
         end if

         ! Y at domhi(2)
         bc_face = get_bc_face(bct_jhi, ng)
         if ( (bc_face == pinf_) .or. (bc_face == pout_) .or. (bc_face == fsw_) ) then
            bc_hi(2) = amrex_lo_neumann
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
         bc_face = get_bc_face(bct_klo, ng)
         if ( (bc_face == pinf_) .or. (bc_face == pout_) .or. (bc_face == fsw_) ) then
            bc_lo(3) = amrex_lo_neumann
         end if

         ! Z at domhi(3)
         bc_face = get_bc_face(bct_khi, ng)
         if ( (bc_face == pinf_) .or. (bc_face == pout_) .or. (bc_face == fsw_) ) then
            bc_hi(3) = amrex_lo_neumann
         end if

      end if

   contains

      !
      ! Test whether the BC type is the same everywhere on
      ! the face. If BC is uniform on face, it returns its value
      !
      function get_bc_face (bct_array, nghost) result (bc_face)
         integer(c_int), intent(in   ) :: bct_array(:,:,:)
         integer       , intent(in   ) :: nghost
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

    end subroutine set_diff_bc

!-----------------------------------------------------------------------!
!-----------------------------------------------------------------------!
!-----------------------------------------------------------------------!
!-----------------------------------------------------------------------!

   subroutine fill_vel_diff_bc(vel_in, vinlo, vinhi, vel, lo, hi, domlo, domhi, ng, &
                       bc_ilo_type, bc_ihi_type, &
                       bc_jlo_type, bc_jhi_type, &
                       bc_klo_type, bc_khi_type)


   use bc, only: minf_, nsw_, fsw_, psw_

   integer,  intent(in   ) ::   vinlo(3),   vinhi(3)
   integer,  intent(in   ) ::    lo(3),    hi(3)
     integer,  intent(in   ) :: domlo(3), domhi(3)

     real(rt), intent(in   ) :: vel_in(vinlo(1):vinhi(1),vinlo(2):vinhi(2),vinlo(3):vinhi(3),3)
     real(rt), intent(  out) ::    vel(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng,3)

     ! BC types
     integer(c_int), intent(in   ) :: ng,  &
          & bc_ilo_type(domlo(2)-ng:domhi(2)+ng,domlo(3)-ng:domhi(3)+ng,2), &
          & bc_ihi_type(domlo(2)-ng:domhi(2)+ng,domlo(3)-ng:domhi(3)+ng,2), &
          & bc_jlo_type(domlo(1)-ng:domhi(1)+ng,domlo(3)-ng:domhi(3)+ng,2), &
          & bc_jhi_type(domlo(1)-ng:domhi(1)+ng,domlo(3)-ng:domhi(3)+ng,2), &
          & bc_klo_type(domlo(1)-ng:domhi(1)+ng,domlo(2)-ng:domhi(2)+ng,2), &
          & bc_khi_type(domlo(1)-ng:domhi(1)+ng,domlo(2)-ng:domhi(2)+ng,2)


     integer :: i,j,k,n

      do k = lo(3)-ng, hi(3)+ng
         do j = lo(2)-ng, hi(2)+ng
            do i = lo(1)-ng, hi(1)+ng
               vel(i,j,k,:) = vel_in(i,j,k,:)
            end do
         end do
      end do

      if ( lo(1) == domlo(1) ) then
         i = lo(1)
         do n = 1, 3
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)

                  if ( ( bc_ilo_type(j,k,1) == MINF_ ) .or. &
                       ( bc_ilo_type(j,k,1) == NSW_ )  .or. &
                       ( bc_ilo_type(j,k,1) == FSW_ )  .or. &
                       ( bc_ilo_type(j,k,1) == PSW_ )  ) then

                     vel(:lo(1)-1,j,k,n) = 2.d0*vel_in(lo(1)-1,j,k,n) - vel_in(lo(1),j,k,n)

                  end if
               end do
            end do
         end do
      end if

      if ( hi(1) == domhi(1) ) then

         i = hi(1)

         do n = 1, 3
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)

                  if ( ( bc_ihi_type(j,k,1) == MINF_ ) .or. &
                       ( bc_ihi_type(j,k,1) == NSW_ )  .or. &
                       ( bc_ihi_type(j,k,1) == FSW_ )  .or. &
                       ( bc_ihi_type(j,k,1) == PSW_ )  ) then

                     vel(hi(1)+1:,j,k,n) = 2.d0*vel_in(hi(1)+1,j,k,n) - vel_in(hi(1),j,k,n)

                  end if
               end do
            end do
         end do
      end if

      if ( lo(2) == domlo(2) ) then

         j = lo(2)

         do n = 1, 3
            do k = lo(3), hi(3)
               do i = lo(1)-1, hi(1)+1

                  if ( ( bc_jlo_type(i,k,1) == MINF_ ) .or. &
                       ( bc_jlo_type(i,k,1) == NSW_ )  .or. &
                       ( bc_jlo_type(i,k,1) == FSW_ )  .or. &
                       ( bc_jlo_type(i,k,1) == PSW_ )  ) then

                     
                     vel(i,:lo(2)-1,k,n) = 2.d0*vel_in(i,lo(2)-1,k,n) - vel_in(i,lo(2),k,n)

                  end if
               end do
            end do
         end do
      end if

      if ( hi(2) == domhi(2) ) then

         j = hi(2)

         do n = 1, 3
            do k = lo(3), hi(3)
               do i = lo(1)-1, hi(1)+1

                  if ( ( bc_jhi_type(i,k,1) == MINF_ ) .or. &
                       ( bc_jhi_type(i,k,1) == NSW_ )  .or. &
                       ( bc_jhi_type(i,k,1) == FSW_ )  .or. &
                       ( bc_jhi_type(i,k,1) == PSW_ )  ) then

                     vel(i,hi(2)+1:,k,n) = 2.d0*vel_in(i,hi(2)+1,k,n) - vel_in(i,hi(2),k,n)

                  end if
               end do
            end do
         end do
      end if


      if ( lo(3) == domlo(3) ) then

         k = lo(3)

         do n = 1, 3
            do j = lo(2)-1, hi(2)+1
               do i = lo(1)-1, hi(1)+1

                  if ( ( bc_klo_type(i,j,1) == MINF_ ) .or. &
                       ( bc_klo_type(i,j,1) == NSW_ )  .or. &
                       ( bc_klo_type(i,j,1) == FSW_ )  .or. &
                       ( bc_klo_type(i,j,1) == PSW_ )  ) then

                     vel(i,j,:lo(3)-1,n) = 2.d0*vel_in(i,j,lo(3)-1,n) - vel_in(i,j,lo(3),n)

                  end if
               end do
            end do
         end do
      end if

      if ( hi(3) == domhi(3) ) then

         k = hi(3)

         do n = 1, 3
            do j = lo(2)-1, hi(2)+1
               do i = lo(1)-1, hi(1)+1

                  if ( ( bc_khi_type(i,j,1) == MINF_ ) .or. &
                       ( bc_khi_type(i,j,1) == NSW_ )  .or. &
                       ( bc_khi_type(i,j,1) == FSW_ )  .or. &
                       ( bc_khi_type(i,j,1) == PSW_ )  ) then

                     vel(i,j,hi(3)+1:,n) = 2.d0*vel_in(i,j,hi(3)+1,n) - vel_in(i,j,hi(3),n)

                  end if
               end do
            end do
         end do
      end if



      !
      ! WHAT'S THE POINT OF THE CODE BELOW??????
      ! 

      ! ! Revisit these
      ! if ( lo(1) == domlo(1) ) then
      !    i = lo(1)
      !    do n = 1, 3
      !       do k = lo(3)-1, hi(3)+1
      !          do j = lo(2)-1, hi(2)+1

      !             if ( ( bc_ilo_type(j,k,1) == MINF_ ) .or. &
      !                  ( bc_ilo_type(j,k,1) == NSW_ )  .or. &
      !                  ( bc_ilo_type(j,k,1) == FSW_ )  .or. &
      !                  ( bc_ilo_type(j,k,1) == PSW_ )  ) then

      !                vel(lo(1)-1,j,k,n) = 2.d0*vel_in(lo(1)-1,j,k,n) - vel_in(lo(1),j,k,n)

      !             end if
      !          end do
      !       end do
      !    end do
      ! end if

      ! ! Revisit these
      ! if ( hi(1) == domhi(1) ) then

      !    i = hi(1)

      !    do n = 1, 3
      !       do k = lo(3)-1, hi(3)+1
      !          do j = lo(2)-1, hi(2)+1


      !             if ( ( bc_ihi_type(j,k,1) == MINF_ ) .or. &
      !                  ( bc_ihi_type(j,k,1) == NSW_ )  .or. &
      !                  ( bc_ihi_type(j,k,1) == FSW_ )  .or. &
      !                  ( bc_ihi_type(j,k,1) == PSW_ )  ) then

      !                vel(hi(1)+1,j,k,n) = 2.d0*vel_in(hi(1)+1,j,k,n) - vel_in(hi(1),j,k,n)

      !             end if
      !          end do
      !       end do
      !    end do
      ! end if

      ! ! Revisit these
      ! if ( lo(2) == domlo(2) ) then

      !    j = lo(2)

      !    do n = 1, 3
      !       do k = lo(3)-1, hi(3)+1
      !          do i = lo(1)-1, hi(1)+1

      !             if ( ( bc_jlo_type(i,k,1) == MINF_ ) .or. &
      !                  ( bc_jlo_type(i,k,1) == NSW_ )  .or. &
      !                  ( bc_jlo_type(i,k,1) == FSW_ )  .or. &
      !                  ( bc_jlo_type(i,k,1) == PSW_ )  ) then

      !                vel(i,lo(2)-1,k,n) = 2.d0*vel_in(i,lo(2)-1,k,n) - vel_in(i,lo(2),k,n)

      !             end if
      !          end do
      !       end do
      !    end do
      ! end if

      ! ! Revisit these
      ! if ( hi(2) == domhi(2) ) then

      !    j = hi(2)

      !    do n = 1, 3
      !       do k = lo(3)-1, hi(3)+1
      !          do i = lo(1)-1, hi(1)+1

      !             if ( ( bc_jhi_type(i,k,1) == MINF_ ) .or. &
      !                  ( bc_jhi_type(i,k,1) == NSW_ )  .or. &
      !                  ( bc_jhi_type(i,k,1) == FSW_ )  .or. &
      !                  ( bc_jhi_type(i,k,1) == PSW_ )  ) then

      !                vel(i,hi(2)+1,k,n) = 2.d0*vel_in(i,hi(2)+1,k,n) - vel_in(i,hi(2),k,n)

      !             end if
      !          end do
      !       end do
      !    end do
      ! end if


    end subroutine fill_vel_diff_bc

end module diffusion_mod
