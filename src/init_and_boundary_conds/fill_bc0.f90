! ::: -----------------------------------------------------------
! ::: This routine is intended to be a generic fill function
! ::: for cell centered data.  It knows how to exrapolate,
! ::: and reflect data and can be used to suppliment problem
! ::: specific fill functions (ie. EXT_DIR).
! :::
! ::: INPUTS/OUTPUTS:
! ::: q        <=  array to fill
! ::: DIMS(q)   => index extent of q array
! ::: domlo,hi  => index extent of problem domain
! ::: dx        => cell spacing
! ::: xlo       => physical location of lower left hand
! :::            corner of q array
! ::: bc	=> array of boundary flags bc(SPACEDIM,lo:hi)
! :::
! ::: NOTE: corner data not used in computing soln but must have
! :::       reasonable values for arithmetic to live
! ::: -----------------------------------------------------------
subroutine fill_bc0(s, slo, shi, &
      bc_ilo_type, bc_ihi_type, bc_jlo_type, bc_jhi_type, &
      bc_klo_type, bc_khi_type, domlo, domhi, ng) &
   bind(C, name="fill_bc0")


! Global modules
!--------------------------------------------------------------------//
      use iso_c_binding , only: c_int
      use amrex_fort_module, only : rt => amrex_real
      use bc, only: NSW_, FSW_, PSW_

      implicit none


! Dummy arguments
!``````````````````````````````````````````````````````````````````````
      integer(c_int), intent(in   ) :: slo(3),shi(3)
      integer(c_int), intent(in   ) :: domlo(3),domhi(3),ng

      real(rt), intent(inout) ::  s&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      integer(c_int), intent(in   ) :: bc_ilo_type&
         (domlo(2)-ng:domhi(2)+ng,domlo(3)-ng:domhi(3)+ng,2)
      integer(c_int), intent(in   ) :: bc_ihi_type&
         (domlo(2)-ng:domhi(2)+ng,domlo(3)-ng:domhi(3)+ng,2)
      integer(c_int), intent(in   ) :: bc_jlo_type&
         (domlo(1)-ng:domhi(1)+ng,domlo(3)-ng:domhi(3)+ng,2)
      integer(c_int), intent(in   ) :: bc_jhi_type&
         (domlo(1)-ng:domhi(1)+ng,domlo(3)-ng:domhi(3)+ng,2)
      integer(c_int), intent(in   ) :: bc_klo_type&
         (domlo(1)-ng:domhi(1)+ng,domlo(2)-ng:domhi(2)+ng,2)
      integer(c_int), intent(in   ) :: bc_khi_type&
         (domlo(1)-ng:domhi(1)+ng,domlo(2)-ng:domhi(2)+ng,2)


! Local variables
!``````````````````````````````````````````````````````````````````````
      integer    nlft, nrgt, nbot, ntop, nup, ndwn
      integer    ilo, ihi, jlo, jhi, klo, khi
      integer    i, j, k
!......................................................................

      nlft = max(0,domlo(1)-slo(1))
      nbot = max(0,domlo(2)-slo(2))
      ndwn = max(0,domlo(3)-slo(3))

      nrgt = max(0,shi(1)-domhi(1))
      ntop = max(0,shi(2)-domhi(2))
      nup  = max(0,shi(3)-domhi(3))

      if (nlft .gt. 0) then
         ilo = domlo(1)
         do i = 1, nlft
            do k=slo(3),shi(3)
               do j=slo(2),shi(2)
               if(bc_ilo_type(j,k,1) == NSW_ .or. &
                  bc_ilo_type(j,k,1) == FSW_ .or. &
                  bc_ilo_type(j,k,1) == PSW_) then
                     s(ilo-i,j,k) = s(ilo,j,k)
                  endif
               end do
            end do
         end do
      endif

      if (nrgt .gt. 0) then
         ihi = domhi(1)
         do i = 1, nrgt
            do k=slo(3),shi(3)
               do j=slo(2),shi(2)
                  if(bc_ihi_type(j,k,1) == NSW_ .or. &
                     bc_ihi_type(j,k,1) == FSW_ .or. &
                     bc_ihi_type(j,k,1) == PSW_) then
                     s(ihi+i,j,k) = s(ihi,j,k)
                  endif
               end do
            end do
         end do
      endif

      if (nbot .gt. 0) then
         jlo = domlo(2)
         do j = 1, nbot
            do k=slo(3),shi(3)
               do i=slo(1),shi(1)
               if(bc_jlo_type(i,k,1) == NSW_ .or. &
                  bc_jlo_type(i,k,1) == FSW_ .or. &
                  bc_jlo_type(i,k,1) == PSW_) then
                     s(i,jlo-j,k) = s(i,jlo,k)
                  endif
               end do
            end do
         end do
      endif

      if (ntop .gt. 0) then
         jhi = domhi(2)
         do j = 1, ntop
            do k=slo(3),shi(3)
               do i=slo(1),shi(1)
               if(bc_jhi_type(i,k,1) == NSW_ .or. &
                  bc_jhi_type(i,k,1) == FSW_ .or. &
                  bc_jhi_type(i,k,1) == PSW_) then
                     s(i,jhi+j,k) = s(i,jhi,k)
                  endif
               end do
            end do
         end do
      endif
      
      if (ndwn .gt. 0) then
         klo = domlo(3)
         do k = 1, ndwn
            do j=slo(2),shi(2)
               do i=slo(1),shi(1)
               if(bc_klo_type(i,j,1) == NSW_ .or. &
                  bc_klo_type(i,j,1) == FSW_ .or. &
                  bc_klo_type(i,j,1) == PSW_) then
                     s(i,j,klo-k) = s(i,j,klo)
                  endif
               end do
            end do
         end do
      endif

      if (nup .gt. 0) then
         khi = domhi(3)
         do k = 1, nup
            do j=slo(2),shi(2)
               do i=slo(1),shi(1)
               if(bc_khi_type(i,j,1) == NSW_ .or. &
                  bc_khi_type(i,j,1) == FSW_ .or. &
                  bc_khi_type(i,j,1) == PSW_) then
                     s(i,j,khi+k) = s(i,j,khi)
                  endif
               end do
            end do
         end do
      endif

! ------------------------------------------------------
! fill edges -------------------------------------------
! ------------------------------------------------------

      do i = 1, nlft
         do j = 1, nbot
            do k=slo(3)+ndwn,shi(3)-nup
               if ( (bc_ilo_type(j,k,1) == NSW_ .or. &
                     bc_ilo_type(j,k,1) == FSW_ .or. &
                     bc_ilo_type(j,k,1) == PSW_) .and. &
                    (bc_jlo_type(i,k,1) == NSW_ .or. &
                     bc_jlo_type(i,k,1) == FSW_ .or. &
                     bc_jlo_type(i,k,1) == PSW_) ) then
                  s(domlo(1)-i,domlo(2)-j,k) = s(domlo(1),domlo(2),k)
               else if (bc_ilo_type(j,k,1) == NSW_ .or. &
                        bc_ilo_type(j,k,1) == FSW_ .or. &
                        bc_ilo_type(j,k,1) == PSW_) then
                  s(domlo(1)-i,domlo(2)-j,k) = s(domlo(1),domlo(2)-j,k)
               else if (bc_jlo_type(i,k,1) == NSW_ .or. &
                        bc_jlo_type(i,k,1) == FSW_ .or. &
                        bc_jlo_type(i,k,1) == PSW_) then
                  s(domlo(1)-i,domlo(2)-j,k) = s(domlo(1)-i,domlo(2),k)
               end if
            end do
         end do
      end do

      do i = 1, nlft
         do j = 1, ntop
            do k=slo(3)+ndwn,shi(3)-nup
               if ( (bc_ilo_type(j,k,1) == NSW_ .or. &
                     bc_ilo_type(j,k,1) == FSW_ .or. &
                     bc_ilo_type(j,k,1) == PSW_) .and. &
                    (bc_jhi_type(i,k,1) == NSW_ .or. &
                     bc_jhi_type(i,k,1) == FSW_ .or. &
                     bc_jhi_type(i,k,1) == PSW_) ) then
                  s(domlo(1)-i,domhi(2)+j,k) = s(domlo(1),domhi(2),k)
               else if (bc_ilo_type(j,k,1) == NSW_ .or. &
                        bc_ilo_type(j,k,1) == FSW_ .or. &
                        bc_ilo_type(j,k,1) == PSW_) then
                  s(domlo(1)-i,domhi(2)+j,k) = s(domlo(1),domhi(2)+j,k)
               else if (bc_jhi_type(i,k,1) == NSW_ .or. &
                        bc_jhi_type(i,k,1) == FSW_ .or. &
                        bc_jhi_type(i,k,1) == PSW_) then
                  s(domlo(1)-i,domhi(2)+j,k) = s(domlo(1)-i,domhi(2),k)
               end if
            end do
         end do
      end do

      do i = 1, nlft
         do k = 1, ndwn
            do j=slo(2)+nbot,shi(2)-ntop
               if ( (bc_ilo_type(j,k,1) == NSW_ .or. &
                     bc_ilo_type(j,k,1) == FSW_ .or. &
                     bc_ilo_type(j,k,1) == PSW_) .and. &
                    (bc_klo_type(i,j,1) == NSW_ .or. &
                     bc_klo_type(i,j,1) == FSW_ .or. &
                     bc_klo_type(i,j,1) == PSW_) ) then
                  s(domlo(1)-i,j,domlo(3)-k) = s(domlo(1),j,domlo(3))
               else if (bc_ilo_type(j,k,1) == NSW_ .or. &
                        bc_ilo_type(j,k,1) == FSW_ .or. &
                        bc_ilo_type(j,k,1) == PSW_) then
                  s(domlo(1)-i,j,domlo(3)-k) = s(domlo(1),j,domlo(3)-k)
               else if (bc_klo_type(i,j,1) == NSW_ .or. &
                        bc_klo_type(i,j,1) == FSW_ .or. &
                        bc_klo_type(i,j,1) == PSW_) then
                  s(domlo(1)-i,j,domlo(3)-k) = s(domlo(1)-i,j,domlo(3))
               end if
            end do
         end do
      end do

      do i = 1, nlft
         do k = 1, nup
            do j=slo(2)+nbot,shi(2)-ntop
               if ( (bc_ilo_type(j,k,1) == NSW_ .or. &
                     bc_ilo_type(j,k,1) == FSW_ .or. &
                     bc_ilo_type(j,k,1) == PSW_) .and. &
                    (bc_khi_type(i,k,1) == NSW_ .or. &
                     bc_khi_type(i,k,1) == FSW_ .or. &
                     bc_khi_type(i,k,1) == PSW_) ) then
                  s(domlo(1)-i,j,domhi(3)+k) = s(domlo(1),j,domhi(3))
               else if (bc_ilo_type(j,k,1) == NSW_ .or. &
                        bc_ilo_type(j,k,1) == FSW_ .or. &
                        bc_ilo_type(j,k,1) == PSW_) then
                  s(domlo(1)-i,j,domhi(3)+k) = s(domlo(1),j,domhi(3)+k)
               else if (bc_khi_type(i,k,1) == NSW_ .or. &
                        bc_khi_type(i,k,1) == FSW_ .or. &
                        bc_khi_type(i,k,1) == PSW_) then
                  s(domlo(1)-i,j,domhi(3)+k) = s(domlo(1)-i,j,domhi(3))
               end if
            end do
         end do
      end do

      do i = 1, nrgt
         do j = 1, nbot
            do k=slo(3)+ndwn,shi(3)-nup
               if ( (bc_ihi_type(j,k,1) == NSW_ .or. &
                     bc_ihi_type(j,k,1) == FSW_ .or. &
                     bc_ihi_type(j,k,1) == PSW_) .and. &
                    (bc_jlo_type(i,k,1) == NSW_ .or. &
                     bc_jlo_type(i,k,1) == FSW_ .or. &
                     bc_jlo_type(i,k,1) == PSW_) ) then
                  s(domhi(1)+i,domlo(2)-j,k) = s(domhi(1),domlo(2),k)
               else if (bc_ihi_type(j,k,1) == NSW_ .or. &
                        bc_ihi_type(j,k,1) == FSW_ .or. &
                        bc_ihi_type(j,k,1) == PSW_) then
                  s(domhi(1)+i,domlo(2)-j,k) = s(domhi(1),domlo(2)-j,k)
               else if (bc_jlo_type(i,k,1) == NSW_ .or. &
                        bc_jlo_type(i,k,1) == FSW_ .or. &
                        bc_jlo_type(i,k,1) == PSW_) then
                  s(domhi(1)+i,domlo(2)-j,k) = s(domhi(1)+i,domlo(2),k)
               end if
            end do
         end do
      end do

      do i = 1, nrgt
         do j = 1, ntop
            do k=slo(3)+ndwn,shi(3)-nup
               if ( (bc_ihi_type(j,k,1) == NSW_ .or. &
                     bc_ihi_type(j,k,1) == FSW_ .or. &
                     bc_ihi_type(j,k,1) == PSW_) .and. &
                    (bc_jhi_type(i,k,1) == NSW_ .or. &
                     bc_jhi_type(i,k,1) == FSW_ .or. &
                     bc_jhi_type(i,k,1) == PSW_) ) then
                  s(domhi(1)+i,domhi(2)+j,k) = s(domhi(1),domhi(2),k)
               else if (bc_ihi_type(j,k,1) == NSW_ .or. &
                        bc_ihi_type(j,k,1) == FSW_ .or. &
                        bc_ihi_type(j,k,1) == PSW_) then
                  s(domhi(1)+i,domhi(2)+j,k) = s(domhi(1),domhi(2)+j,k)
               else if (bc_jhi_type(i,k,1) == NSW_ .or. &
                        bc_jhi_type(i,k,1) == FSW_ .or. &
                        bc_jhi_type(i,k,1) == PSW_) then
                  s(domhi(1)+i,domhi(2)+j,k) = s(domhi(1)+i,domhi(2),k)
               end if
            end do
         end do
      end do

      do i = 1, nrgt
         do k = 1, ndwn
            do j=slo(2)+nbot,shi(2)-ntop
               if ( (bc_ihi_type(j,k,1) == NSW_ .or. &
                     bc_ihi_type(j,k,1) == FSW_ .or. &
                     bc_ihi_type(j,k,1) == PSW_) .and. &
                    (bc_klo_type(i,j,1) == NSW_ .or. &
                     bc_klo_type(i,j,1) == FSW_ .or. &
                     bc_klo_type(i,j,1) == PSW_) ) then
                  s(domhi(1)+i,j,domlo(3)-k) = s(domhi(1),j,domlo(3))
               else if (bc_ihi_type(j,k,1) == NSW_ .or. &
                        bc_ihi_type(j,k,1) == FSW_ .or. &
                        bc_ihi_type(j,k,1) == PSW_) then
                  s(domhi(1)+i,j,domlo(3)-k) = s(domhi(1),j,domlo(3)-k)
               else if (bc_klo_type(i,j,1) == NSW_ .or. &
                        bc_klo_type(i,j,1) == FSW_ .or. &
                        bc_klo_type(i,j,1) == PSW_) then
                  s(domhi(1)+i,j,domlo(3)-k) = s(domhi(1)+i,j,domlo(3))
               end if
            end do
         end do
      end do

      do i = 1, nrgt
         do k = 1, nup
            do j=slo(2)+nbot,shi(2)-ntop
               if ( (bc_ihi_type(j,k,1) == NSW_ .or. &
                     bc_ihi_type(j,k,1) == FSW_ .or. &
                     bc_ihi_type(j,k,1) == PSW_) .and. &
                    (bc_khi_type(i,j,1) == NSW_ .or. &
                     bc_khi_type(i,j,1) == FSW_ .or. &
                     bc_khi_type(i,j,1) == PSW_) ) then
                  s(domhi(1)+i,j,domhi(3)+k) = s(domhi(1),j,domhi(3))
               else if (bc_ihi_type(j,k,1) == NSW_ .or. &
                        bc_ihi_type(j,k,1) == FSW_ .or. &
                        bc_ihi_type(j,k,1) == PSW_) then
                  s(domhi(1)+i,j,domhi(3)+k) = s(domhi(1),j,domhi(3)+k)
               else if (bc_khi_type(i,j,1) == NSW_ .or. &
                        bc_khi_type(i,j,1) == FSW_ .or. &
                        bc_khi_type(i,j,1) == PSW_) then
                  s(domhi(1)+i,j,domhi(3)+k) = s(domhi(1)+i,j,domhi(3))
               end if
            end do
         end do
      end do

      do j = 1, nbot
         do k = 1, ndwn
            do i=slo(1)+nlft,shi(1)-nrgt
               if ( (bc_klo_type(i,j,1) == NSW_ .or. &
                     bc_klo_type(i,j,1) == FSW_ .or. &
                     bc_klo_type(i,j,1) == PSW_) .and. &
                    (bc_jlo_type(i,k,1) == NSW_ .or. &
                     bc_jlo_type(i,k,1) == FSW_ .or. &
                     bc_jlo_type(i,k,1) == PSW_) ) then
                  s(i,domlo(2)-j,domlo(3)-k) = s(i,domlo(2),domlo(3))
               else if (bc_klo_type(i,j,1) == NSW_ .or. &
                        bc_klo_type(i,j,1) == FSW_ .or. &
                        bc_klo_type(i,j,1) == PSW_) then
                  s(i,domlo(2)-j,domlo(3)-k) = s(i,domlo(2)-j,domlo(3))
               else if (bc_jlo_type(i,k,1) == NSW_ .or. &
                        bc_jlo_type(i,k,1) == FSW_ .or. &
                        bc_jlo_type(i,k,1) == PSW_) then
                  s(i,domlo(2)-j,domlo(3)-k) = s(i,domlo(2),domlo(3)-k)
               end if
            end do
         end do
      end do

      do j = 1, ntop
         do k = 1, ndwn
            do i=slo(1)+nlft,shi(1)-nrgt
               if ( (bc_klo_type(i,j,1) == NSW_ .or. &
                     bc_klo_type(i,j,1) == FSW_ .or. &
                     bc_klo_type(i,j,1) == PSW_) .and. &
                    (bc_jhi_type(i,k,1) == NSW_ .or. &
                     bc_jhi_type(i,k,1) == FSW_ .or. &
                     bc_jhi_type(i,k,1) == PSW_) ) then
                  s(i,domhi(2)+j,domlo(3)-k) = s(i,domhi(2),domlo(3))
               else if (bc_klo_type(i,j,1) == NSW_ .or. &
                        bc_klo_type(i,j,1) == FSW_ .or. &
                        bc_klo_type(i,j,1) == PSW_) then
                  s(i,domhi(2)+j,domlo(3)-k) = s(i,domhi(2)+j,domlo(3))
               else if (bc_jhi_type(i,k,1) == NSW_ .or. &
                        bc_jhi_type(i,k,1) == FSW_ .or. &
                        bc_jhi_type(i,k,1) == PSW_) then
                  s(i,domhi(2)+j,domlo(3)-k) = s(i,domhi(2),domlo(3)-k)
               end if
            end do
         end do
      end do

      do j = 1, nbot
         do k = 1, nup
            do i=slo(1)+nlft,shi(1)-nrgt
               if ( (bc_khi_type(i,j,1) == NSW_ .or. &
                     bc_khi_type(i,j,1) == FSW_ .or. &
                     bc_khi_type(i,j,1) == PSW_) .and. &
                    (bc_jlo_type(i,k,1) == NSW_ .or. &
                     bc_jlo_type(i,k,1) == FSW_ .or. &
                     bc_jlo_type(i,k,1) == PSW_) ) then
                  s(i,domlo(2)-j,domhi(3)+k) = s(i,domlo(2),domhi(3))
               else if (bc_khi_type(i,j,1) == NSW_ .or. &
                        bc_khi_type(i,j,1) == FSW_ .or. &
                        bc_khi_type(i,j,1) == PSW_) then
                  s(i,domlo(2)-j,domhi(3)+k) = s(i,domlo(2)-j,domhi(3))
               else if (bc_jlo_type(i,k,1) == NSW_ .or. &
                        bc_jlo_type(i,k,1) == FSW_ .or. &
                        bc_jlo_type(i,k,1) == PSW_) then
                  s(i,domlo(2)-j,domhi(3)+k) = s(i,domlo(2),domhi(3)+k)
               end if
            end do
         end do
      end do

      do j = 1, ntop
         do k = 1, nup
            do i=slo(1)+nlft,shi(1)-nrgt
               if ( (bc_khi_type(i,j,1) == NSW_ .or. &
                     bc_khi_type(i,j,1) == FSW_ .or. &
                     bc_khi_type(i,j,1) == PSW_) .and. &
                    (bc_jhi_type(i,k,1) == NSW_ .or. &
                     bc_jhi_type(i,k,1) == FSW_ .or. &
                     bc_jhi_type(i,k,1) == PSW_) ) then
                  s(i,domhi(2)+j,domhi(3)+k) = s(i,domhi(2),domhi(3))
               else if (bc_khi_type(i,j,1) == NSW_ .or. &
                        bc_khi_type(i,j,1) == FSW_ .or. &
                        bc_khi_type(i,j,1) == PSW_) then
                  s(i,domhi(2)+j,domhi(3)+k) = s(i,domhi(2)+j,domhi(3))
               else if (bc_jhi_type(i,k,1) == NSW_ .or. &
                        bc_jhi_type(i,k,1) == FSW_ .or. &
                        bc_jhi_type(i,k,1) == PSW_) then
                  s(i,domhi(2)+j,domhi(3)+k) = s(i,domhi(2),domhi(3)+k)
               end if
            end do
         end do
      end do

      return
   end subroutine fill_bc0
