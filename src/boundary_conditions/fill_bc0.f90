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
! ::: bc        => array of boundary flags bc(SPACEDIM,lo:hi)
! :::
! ::: NOTE: corner data not used in computing soln but must have
! :::       reasonable values for arithmetic to live
! ::: -----------------------------------------------------------
subroutine fill_bc0(s, slo, shi, &
                    bc_ilo_type, bc_ihi_type, bc_jlo_type, bc_jhi_type, &
                    bc_klo_type, bc_khi_type, domlo, domhi, ng) &
   bind(C, name="fill_bc0")

   use iso_c_binding , only: c_int
   use amrex_fort_module, only : rt => amrex_real
   use bc, only: NSW_, FSW_, PSW_

   implicit none

   ! Array bounds
   integer(c_int), intent(in   ) :: slo(3),shi(3)

   ! Domain bounds
   integer(c_int), intent(in   ) :: domlo(3),domhi(3)

   ! Number of ghost cells
   integer(c_int), intent(in   ) ::   ng

   ! Arrays
   real(rt), intent(inout) :: s(slo(1):shi(1), slo(2):shi(2), slo(3):shi(3))

   integer(c_int), intent(in) :: bc_ilo_type(domlo(2)-ng:domhi(2)+ng, domlo(3)-ng:domhi(3)+ng, 2)
   integer(c_int), intent(in) :: bc_ihi_type(domlo(2)-ng:domhi(2)+ng, domlo(3)-ng:domhi(3)+ng, 2)
   integer(c_int), intent(in) :: bc_jlo_type(domlo(1)-ng:domhi(1)+ng, domlo(3)-ng:domhi(3)+ng, 2)
   integer(c_int), intent(in) :: bc_jhi_type(domlo(1)-ng:domhi(1)+ng, domlo(3)-ng:domhi(3)+ng, 2)
   integer(c_int), intent(in) :: bc_klo_type(domlo(1)-ng:domhi(1)+ng, domlo(2)-ng:domhi(2)+ng, 2)
   integer(c_int), intent(in) :: bc_khi_type(domlo(1)-ng:domhi(1)+ng, domlo(2)-ng:domhi(2)+ng, 2)



   ! Local variables
   integer :: nlft, nrgt, nbot, ntop, nup, ndwn
   integer :: ilo, ihi, jlo, jhi, klo, khi
   integer :: i, j, k

   ! These are the BCS for which we need to extrapolate
   integer, parameter :: valid_bcs(3) = [NSW_, FSW_, PSW_]

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
               if(any(bc_ilo_type(j,k,1) == valid_bcs)) s(ilo-i,j,k) = s(ilo,j,k)
            end do
         end do
      end do
   endif

   if (nrgt .gt. 0) then
      ihi = domhi(1)
      do i = 1, nrgt
         do k=slo(3),shi(3)
            do j=slo(2),shi(2)
               if(any(bc_ihi_type(j,k,1) == valid_bcs)) s(ihi+i,j,k) = s(ihi,j,k)
            end do
         end do
      end do
   endif

   if (nbot .gt. 0) then
      jlo = domlo(2)
      do j = 1, nbot
         do k=slo(3),shi(3)
            do i=slo(1),shi(1)
               if(any(bc_jlo_type(i,k,1) == valid_bcs)) s(i,jlo-j,k) = s(i,jlo,k)
            end do
         end do
      end do
   endif

   if (ntop .gt. 0) then
      jhi = domhi(2)
      do j = 1, ntop
         do k=slo(3),shi(3)
            do i=slo(1),shi(1)
               if(any(bc_jhi_type(i,k,1) == valid_bcs)) s(i,jhi+j,k) = s(i,jhi,k)
            end do
         end do
      end do
   endif

   if (ndwn .gt. 0) then
      klo = domlo(3)
      do k = 1, ndwn
         do j=slo(2),shi(2)
            do i=slo(1),shi(1)
               if(any(bc_klo_type(i,j,1) == valid_bcs)) s(i,j,klo-k) = s(i,j,klo)
            end do
         end do
      end do
   endif

   if (nup .gt. 0) then
      khi = domhi(3)
      do k = 1, nup
         do j=slo(2),shi(2)
            do i=slo(1),shi(1)
               if(any(bc_khi_type(i,j,1) == valid_bcs)) s(i,j,khi+k) = s(i,j,khi)
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
            if ( any( bc_ilo_type(j,k,1) == valid_bcs )   .and. &
             &   any( bc_jlo_type(i,k,1) == valid_bcs ) ) then
               s(domlo(1)-i,domlo(2)-j,k) = s(domlo(1),domlo(2),k)
            else if ( any( bc_ilo_type(j,k,1) == valid_bcs ) ) then 
               s(domlo(1)-i,domlo(2)-j,k) = s(domlo(1),domlo(2)-j,k)
            else if ( any( bc_jlo_type(i,k,1) == valid_bcs ) ) then
               s(domlo(1)-i,domlo(2)-j,k) = s(domlo(1)-i,domlo(2),k)
            end if
         end do
      end do
   end do

   do i = 1, nlft
      do j = 1, ntop
         do k=slo(3)+ndwn,shi(3)-nup
            if ( any( bc_ilo_type(j,k,1) == valid_bcs )   .and. &
             &   any( bc_jhi_type(i,k,1) == valid_bcs ) ) then
               s(domlo(1)-i,domhi(2)+j,k) = s(domlo(1),domhi(2),k)
            else if ( any( bc_ilo_type(j,k,1) == valid_bcs ) ) then
               s(domlo(1)-i,domhi(2)+j,k) = s(domlo(1),domhi(2)+j,k)
            else if ( any( bc_jhi_type(i,k,1) == valid_bcs ) ) then
               s(domlo(1)-i,domhi(2)+j,k) = s(domlo(1)-i,domhi(2),k)
            end if
         end do
      end do
   end do

   do i = 1, nlft
      do k = 1, ndwn
         do j=slo(2)+nbot,shi(2)-ntop
            if ( any( bc_ilo_type(j,k,1) == valid_bcs )   .and. &
             &   any( bc_klo_type(i,j,1) == valid_bcs ) ) then 
               s(domlo(1)-i,j,domlo(3)-k) = s(domlo(1),j,domlo(3))
            else if ( any( bc_ilo_type(j,k,1) == valid_bcs ) ) then
               s(domlo(1)-i,j,domlo(3)-k) = s(domlo(1),j,domlo(3)-k)
            else if ( any( bc_klo_type(i,j,1) == valid_bcs ) ) then 
               s(domlo(1)-i,j,domlo(3)-k) = s(domlo(1)-i,j,domlo(3))
            end if
         end do
      end do
   end do

   do i = 1, nlft
      do k = 1, nup
         do j=slo(2)+nbot,shi(2)-ntop
            if ( any( bc_ilo_type(j,k,1) == valid_bcs )   .and. & 
             &   any( bc_khi_type(i,k,1) == valid_bcs ) ) then 
               s(domlo(1)-i,j,domhi(3)+k) = s(domlo(1),j,domhi(3))
            else if ( any( bc_ilo_type(j,k,1) == valid_bcs ) ) then 
               s(domlo(1)-i,j,domhi(3)+k) = s(domlo(1),j,domhi(3)+k)
            else if ( any( bc_khi_type(i,k,1) == valid_bcs ) ) then 
               s(domlo(1)-i,j,domhi(3)+k) = s(domlo(1)-i,j,domhi(3))
            end if
         end do
      end do
   end do

   do i = 1, nrgt
      do j = 1, nbot
         do k=slo(3)+ndwn,shi(3)-nup
            if ( any( bc_ihi_type(j,k,1) == valid_bcs )   .and. &
             &   any( bc_jlo_type(i,k,1) == valid_bcs ) ) then
               s(domhi(1)+i,domlo(2)-j,k) = s(domhi(1),domlo(2),k)
            else if ( any( bc_ihi_type(j,k,1) == valid_bcs ) ) then
               s(domhi(1)+i,domlo(2)-j,k) = s(domhi(1),domlo(2)-j,k)
            else if ( any( bc_jlo_type(i,k,1) == valid_bcs ) ) then
               s(domhi(1)+i,domlo(2)-j,k) = s(domhi(1)+i,domlo(2),k)
            end if
         end do
      end do
   end do

   do i = 1, nrgt
      do j = 1, ntop
         do k=slo(3)+ndwn,shi(3)-nup
            if ( any( bc_ihi_type(j,k,1) == valid_bcs )   .and. &
             &   any( bc_jhi_type(i,k,1) == valid_bcs ) ) then
               s(domhi(1)+i,domhi(2)+j,k) = s(domhi(1),domhi(2),k)
            else if ( any( bc_ihi_type(j,k,1) == valid_bcs ) ) then 
               s(domhi(1)+i,domhi(2)+j,k) = s(domhi(1),domhi(2)+j,k)
            else if ( any( bc_jhi_type(i,k,1) == valid_bcs ) ) then 
               s(domhi(1)+i,domhi(2)+j,k) = s(domhi(1)+i,domhi(2),k)
            end if
         end do
      end do
   end do

   do i = 1, nrgt
      do k = 1, ndwn
         do j=slo(2)+nbot,shi(2)-ntop
            if ( any( bc_ihi_type(j,k,1) == valid_bcs )   .and. &
             &   any( bc_klo_type(i,j,1) == valid_bcs ) ) then 
               s(domhi(1)+i,j,domlo(3)-k) = s(domhi(1),j,domlo(3))
            else if ( any( bc_ihi_type(j,k,1) == valid_bcs ) ) then
               s(domhi(1)+i,j,domlo(3)-k) = s(domhi(1),j,domlo(3)-k)
            else if ( any( bc_klo_type(i,j,1) == valid_bcs ) ) then 
               s(domhi(1)+i,j,domlo(3)-k) = s(domhi(1)+i,j,domlo(3))
            end if
         end do
      end do
   end do

   do i = 1, nrgt
      do k = 1, nup
         do j=slo(2)+nbot,shi(2)-ntop
            if ( any( bc_ihi_type(j,k,1) == valid_bcs )   .and. &
             &   any( bc_khi_type(i,j,1) == valid_bcs ) ) then 
               s(domhi(1)+i,j,domhi(3)+k) = s(domhi(1),j,domhi(3))
            else if ( any( bc_ihi_type(j,k,1) == valid_bcs ) ) then
               s(domhi(1)+i,j,domhi(3)+k) = s(domhi(1),j,domhi(3)+k)
            else if ( any( bc_khi_type(i,j,1) == valid_bcs ) ) then 
               s(domhi(1)+i,j,domhi(3)+k) = s(domhi(1)+i,j,domhi(3))
            end if
         end do
      end do
   end do

   do j = 1, nbot
      do k = 1, ndwn
         do i=slo(1)+nlft,shi(1)-nrgt
            if ( any( bc_klo_type(i,j,1) == valid_bcs )   .and. &
             &   any( bc_jlo_type(i,k,1) == valid_bcs ) ) then
               s(i,domlo(2)-j,domlo(3)-k) = s(i,domlo(2),domlo(3))
            else if ( any( bc_klo_type(i,j,1) == valid_bcs ) ) then
               s(i,domlo(2)-j,domlo(3)-k) = s(i,domlo(2)-j,domlo(3))
            else if ( any( bc_jlo_type(i,k,1) == valid_bcs ) ) then
               s(i,domlo(2)-j,domlo(3)-k) = s(i,domlo(2),domlo(3)-k)
            end if
         end do
      end do
   end do

   do j = 1, ntop
      do k = 1, ndwn
         do i=slo(1)+nlft,shi(1)-nrgt
            if ( any( bc_klo_type(i,j,1) == valid_bcs )   .and. &
             &   any( bc_jhi_type(i,k,1) == valid_bcs ) ) then
               s(i,domhi(2)+j,domlo(3)-k) = s(i,domhi(2),domlo(3))
            else if ( any( bc_klo_type(i,j,1) == valid_bcs ) ) then
               s(i,domhi(2)+j,domlo(3)-k) = s(i,domhi(2)+j,domlo(3))
            else if ( any( bc_jhi_type(i,k,1) == valid_bcs ) ) then 
               s(i,domhi(2)+j,domlo(3)-k) = s(i,domhi(2),domlo(3)-k)
            end if
         end do
      end do
   end do

   do j = 1, nbot
      do k = 1, nup
         do i=slo(1)+nlft,shi(1)-nrgt
            if ( any( bc_khi_type(i,j,1) == valid_bcs )   .and. &
             &   any( bc_jlo_type(i,k,1) == valid_bcs ) ) then 
               s(i,domlo(2)-j,domhi(3)+k) = s(i,domlo(2),domhi(3))
            else if ( any( bc_khi_type(i,j,1) == valid_bcs ) ) then 
               s(i,domlo(2)-j,domhi(3)+k) = s(i,domlo(2)-j,domhi(3))
            else if ( any( bc_jlo_type(i,k,1) == valid_bcs ) ) then 
               s(i,domlo(2)-j,domhi(3)+k) = s(i,domlo(2),domhi(3)+k)
            end if
         end do
      end do
   end do

   do j = 1, ntop
      do k = 1, nup
         do i=slo(1)+nlft,shi(1)-nrgt
            if ( any( bc_khi_type(i,j,1) == valid_bcs )   .and. &
               & any( bc_jhi_type(i,k,1) == valid_bcs ) ) then 
               s(i,domhi(2)+j,domhi(3)+k) = s(i,domhi(2),domhi(3))
            else if ( any( bc_khi_type(i,j,1) == valid_bcs ) ) then 
               s(i,domhi(2)+j,domhi(3)+k) = s(i,domhi(2)+j,domhi(3))
            else if ( any( bc_jhi_type(i,k,1) == valid_bcs ) ) then 
               s(i,domhi(2)+j,domhi(3)+k) = s(i,domhi(2),domhi(3)+k)
            end if
         end do
      end do
   end do

   return
end subroutine fill_bc0
