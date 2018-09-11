!                                                                      !
!  Subroutine: set_p0                                                  !
!  Author: M. Syamlal                                 Date: 21-JAN-92  !
!                                                                      !
!  Purpose: Set the pressure field inside the bed assuming gravity     !
!           is acting in the negative y-direction.                     !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      subroutine set_p0(lo, hi, domlo, domhi, &
                        p0_g, slo, shi, &
                        dx, dy, dz, xlength, ylength, zlength, delp_dir, &
                        bct_ilo, bct_ihi, bct_jlo, bct_jhi, &
                        bct_klo, bct_khi, ng, nodal_pressure) &
                 bind(C, name="set_p0")

      use bc       , only: delp_x, delp_y, delp_z
      use bc       , only: dim_bc, bc_type, bc_p_g, bc_defined
      use bc       , only: pinf_, pout_, minf_
      use constant , only: gravity
      use fld_const, only: ro_g0
      use ic       , only: ic_p_g, ic_defined
      use scales   , only: scale_pressure

      use amrex_fort_module, only : ar => amrex_real
      use iso_c_binding , only: c_int
      use param   , only: zero, undefined
      use param   , only: is_defined, is_undefined
      use param, only: dim_ic

      implicit none

      integer, intent(in) :: slo(3), shi(3), lo(3), hi(3)
      integer, intent(in) :: domlo(3), domhi(3), ng, nodal_pressure

      real(ar), intent(inout) :: p0_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(ar), intent(in) :: dx, dy, dz
      real(ar), intent(in) :: xlength, ylength, zlength
      integer     , intent(in) :: delp_dir

      integer(c_int), intent(in   ) :: &
           bct_ilo(domlo(2)-ng:domhi(2)+ng,domlo(3)-ng:domhi(3)+ng,2), &
           bct_ihi(domlo(2)-ng:domhi(2)+ng,domlo(3)-ng:domhi(3)+ng,2), &
           bct_jlo(domlo(1)-ng:domhi(1)+ng,domlo(3)-ng:domhi(3)+ng,2), &
           bct_jhi(domlo(1)-ng:domhi(1)+ng,domlo(3)-ng:domhi(3)+ng,2), &
           bct_klo(domlo(1)-ng:domhi(1)+ng,domlo(2)-ng:domhi(2)+ng,2), &
           bct_khi(domlo(1)-ng:domhi(1)+ng,domlo(2)-ng:domhi(2)+ng,2)

!-----------------------------------------------
! Local variables
!-----------------------------------------------
      ! indices
      integer :: i, j, k, icv, bcv
      integer :: nlft, nbot, ndwn, nrgt, ntop, nup

      ! Gas pressure at the axial location j
      real(ar) :: pj

      ! Average pressure drop per unit length
      real(ar) :: dpodx, dpody, dpodz

! ---------------------------------------------------------------->>>

      !  Make sure that ic_p_g is set if using delp pressure conditions
      do icv = 1, dim_ic
         if (ic_defined(icv)) then
            if ( delp_dir .ge. 0 ) then
               if (.not. is_defined(ic_p_g(icv))) then
                  print *,'MUST DEFINE ic_p_g if using the DELP pressure condition'
                  stop
               end if
               pj = ic_p_g(icv)
            else
               if (is_undefined(ic_p_g(icv))) goto 60
               if (gravity(1).ne.0.d0 .or. gravity(2).ne.0.d0 .or. gravity(3).ne.0.d0) goto 60
               p0_g(:,:,:) = ic_p_g(icv)
            end if
         end if
      end do

! ---------------------------------------------------------------->>>

      ! Here the pressure in each cell is determined from a specified pressure
      ! drop across the domain length. This section requires that the pressure
      ! is already defined in all initial condition regions (otherwise this
      ! section would be skipped)
      block

         !  This hack allows to set the IC pressure  at L-dx/2 for both
         !  nodal and CC pressure -> reference value for pressure, AKA IC_P_G,
         !  is set at the last cell center location.
         real(ar) :: offset 
         
         if ( nodal_pressure == 1 ) then
            offset = - 0.5_ar
         else
            offset = 0.0_ar
         end if

         if (abs(delp_x) > epsilon(zero)) then
            dpodx = delp_x/xlength
            pj = pj - dpodx*dx*(hi(1)-domhi(1)+2 + offset )
            do i = shi(1), slo(1), -1
               pj = pj + dpodx*dx
               p0_g(i,slo(2):shi(2),slo(3):shi(3)) = scale_pressure(pj)
            enddo
         endif

         if (abs(delp_y) > epsilon(zero)) then
            dpody = delp_y/ylength
            pj = pj - dpody*dy*(hi(2)-domhi(2)+2 + offset )
            do j = shi(2), slo(2), -1
               pj = pj + dpody*dy
               p0_g(slo(1):shi(1),j,slo(3):shi(3)) = scale_pressure(pj)
            enddo
         endif

         if (abs(delp_z) > epsilon(zero)) then
            dpodz = delp_z/zlength
            pj = pj - dpodz*dz*(hi(3)-domhi(3)+2 + offset )
            do k = shi(3), slo(3), -1
               pj = pj + dpodz*dz
               p0_g(slo(1):shi(1),slo(2):shi(2),k) = scale_pressure(pj)
            end do
         endif

      end block
         
      GOTO 100   ! pressure in all intial condition region cells was defined

! ----------------------------------------------------------------<<<

   60 CONTINUE   ! pressure in an initial condition region cell was undefined

! ---------------------------------------------------------------->>>

      ! Search for an outflow boundary condition where pressure is specified
      pj = undefined
      do icv = 1, dim_bc
         if (bc_defined(icv)) then
            if(bc_type(icv)=='P_OUTFLOW' .or. bc_type(icv)=='PO') &
               pj = bc_p_g(icv)
         endif
      enddo

      ! Either a PO was not specified and/or a PO was specified but not the
      ! pressure at the outlet
      if (is_undefined(pj)) then
         p0_g = zero
         goto 100
      endif

! ----------------------------------------------------------------<<<

      ! Set an approximate pressure field assuming that the pressure drop
      ! balances the weight of the bed, if the initial pressure-field is not
      ! specified

      if (abs(gravity(1)) > epsilon(0.0d0)) then

         ! Find the average weight per unit area over an x-z slice
         dpodx = -gravity(1)*ro_g0

         if (gravity(1) <= 0.0d0) then
            do i = domhi(1)+1, domlo(1), -1
               if (i <= shi(1) .and. i >= slo(1)) &
                  p0_g(i,:,:) = scale_pressure(pj)
               pj = pj + dpodx*dx
            enddo
         else
            do i = domlo(1), domhi(1)+1
               if (i <= shi(1) .and. i >= slo(1)) &
                  p0_g(i,:,:) = scale_pressure(pj)
               pj = pj - dpodx*dx
            enddo
         endif

      else if (abs(gravity(2)) > epsilon(0.0d0)) then

         dpody = -gravity(2)*ro_g0

         if (gravity(2) <= 0.0d0) then
            do j = domhi(2)+1, domlo(2), -1
               if (j <= shi(2) .and. j >= slo(2)) &
                  p0_g(:,j,:) = scale_pressure(pj)
               pj = pj + dpody*dy
            enddo
         else
            do j = domlo(2),domhi(2)+1
               if (j <= shi(2) .and. j >= slo(2)) &
                  p0_g(:,j,:) = scale_pressure(pj)
               pj = pj - dpody*dy
            enddo
         endif

      else if (abs(gravity(3)) > epsilon(0.0d0)) then

         dpodz = -gravity(3)*ro_g0

         if(gravity(3) <= 0.0d0) then
            do k = domhi(3)+1, domlo(3), -1
               if (k <= shi(3) .and. k >= slo(3)) &
                  p0_g(:,:,k) = scale_pressure(pj)
               pj = pj + dpodz*dz
            enddo
         else
            do k = domlo(3),domhi(3)+1
               if (k <= shi(3) .and. k >= slo(3)) &
                  p0_g(:,:,k) = scale_pressure(pj)
               pj = pj - dpodz*dz
            enddo
         endif
      endif

! ----------------------------------------------------------------<<<

  100 continue

! ---------------------------------------------------------------->>>

      block
         integer offset

         if (nodal_pressure == 1) then
            offset = 1
         else
            offset = 0
         end if
         
         nlft = max(0,domlo(1)-slo(1)+offset)
         nbot = max(0,domlo(2)-slo(2)+offset)
         ndwn = max(0,domlo(3)-slo(3)+offset)
         
         nrgt = max(0,shi(1)-domhi(1))
         ntop = max(0,shi(2)-domhi(2))
         nup  = max(0,shi(3)-domhi(3))
      end block

      if (nlft .gt. 0) then
         do k=slo(3),shi(3)
            do j=slo(2),shi(2)

               select case ( bct_ilo(j,k,1) )

               case (pinf_, pout_)

                   bcv = bct_ilo(j,k,2)
                   if (nodal_pressure .eq. 1) then
                       p0_g(slo(1):domlo(1)  ,j,k) = scale_pressure(bc_p_g(bcv))
                   else
                       p0_g(slo(1):domlo(1)-1,j,k) = scale_pressure(bc_p_g(bcv))
                   endif

               case (minf_)

                   if (nodal_pressure .eq. 0) &
                       p0_g(slo(1):domlo(1)-1,j,k) = &
                           2.d0 * p0_g(domlo(1),j,k) - p0_g(domlo(1)+1,j,k)

               end select
            end do
         end do
      endif

      if (nrgt .gt. 0) then
         do k=slo(3),shi(3)
            do j=slo(2),shi(2)

               select case ( bct_ihi(j,k,1) )

               case (pinf_, pout_)

                   bcv = bct_ihi(j,k,2)
                   p0_g(domhi(1)+1:shi(1),j,k) = scale_pressure(bc_p_g(bcv))

               case (minf_)

                   if (nodal_pressure .eq. 0) &
                       p0_g(domhi(1)+1:shi(1),j,k) = &
                           2.d0 * p0_g(domhi(1),j,k) - p0_g(domhi(1)-1,j,k)

               end select
            end do
         end do
      endif

      if (nbot .gt. 0) then
         do k=slo(3),shi(3)
            do i=slo(1),shi(1)


               select case ( bct_jlo(i,k,1) )

               case (pinf_, pout_)

                   bcv = bct_jlo(i,k,2)
                   if (nodal_pressure .eq. 1) then
                       p0_g(i,slo(2):domlo(2)  ,k) = scale_pressure(bc_p_g(bcv))
                   else
                       p0_g(i,slo(2):domlo(2)-1,k) = scale_pressure(bc_p_g(bcv))
                   endif

               case (minf_)

                   if (nodal_pressure .eq. 0) &
                       p0_g(i,slo(2):domlo(2)-1,k) = &
                           2.d0 * p0_g(i,domlo(2),k) - p0_g(i,domlo(2)+1,k)

               end select
            end do
         end do
      endif

      if (ntop .gt. 0) then
         do k = slo(3),shi(3)
            do i = slo(1),shi(1)

               select case ( bct_jhi(i,k,1) )

               case (pinf_, pout_)

                  bcv = bct_jhi(i,k,2)
                  p0_g(i,domhi(2)+1:shi(2),k) = scale_pressure(bc_p_g(bcv))

               case (minf_)

                   if (nodal_pressure .eq. 0) &
                       p0_g(i,domhi(2)+1:shi(2),k) = &
                           2.d0 * p0_g(i,domhi(2),k) - p0_g(i,domhi(2)-1,k)

               end select
            end do
         end do
      endif

      if (ndwn .gt. 0) then
         do j=slo(2),shi(2)
            do i=slo(1),shi(1)

               select case ( bct_klo(i,j,1) )

               case (pinf_, pout_)

                  bcv = bct_klo(i,j,2)
                  if (nodal_pressure .eq. 1) then
                      p0_g(i,j,slo(3):domlo(3)  ) = scale_pressure(bc_p_g(bcv))
                  else
                      p0_g(i,j,slo(3):domlo(3)-1) = scale_pressure(bc_p_g(bcv))
                  endif

               case (minf_)

                   if (nodal_pressure .eq. 0) &
                       p0_g(i,j,slo(3):domlo(3)-1) = &
                           2.d0 * p0_g(i,j,domlo(3)) - p0_g(i,j,domlo(3)+1)

               end select
            end do
         end do
      endif

      if (nup .gt. 0) then
         do j=slo(2),shi(2)
            do i=slo(1),shi(1)

              select case ( bct_khi(i,j,1) )

               case (pinf_, pout_)

                   bcv = bct_khi(i,j,2)
                   p0_g(i,j,domhi(3)+1:shi(3)) = scale_pressure(bc_p_g(bcv))

               case (minf_)

                   if (nodal_pressure .eq. 0) &
                       p0_g(i,j,domhi(3)+1:shi(3)) = &
                           2.d0 * p0_g(i,j,domhi(3)) - p0_g(i,j,domhi(3)-1)

               end select
            end do
         end do
      endif

   end subroutine set_p0
