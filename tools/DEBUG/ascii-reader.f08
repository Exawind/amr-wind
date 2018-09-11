!!----------------------------------------------------------------------!
!!                                                                      !
!! Subprogram to read AMReX generated particle ASCII files. Below is an !
!! example program that would use this routine.                         !
!!                                                                      !
!!----------------------------------------------------------------------!
!!program example
!!
!!  use data
!!
!!  implicit none
!!
!!  integer :: lc, ios
!!  integer :: pc, rc, ic
!!
!!  double precision, allocatable :: rdata(:,:)
!!  integer,          allocatable :: idata(:,:)
!!
!!  integer :: id, min_id, max_id
!!  integer, allocatable :: ids(:)
!!
!!  call open_data(pc,rc, ic)
!!
!!  allocate(rdata(pc,rc))
!!  allocate(idata(pc,ic))
!!
!!  call read_data(pc, rc, ic, rdata, idata)
!!
!!end program example
!!
!!----------------------------------------------------------------------!
!!
module ascii_reader

  integer, parameter :: fdata = 1000
  logical, parameter :: debug = .false.

  private
  public :: open_ascii_file, read_ascii_file

contains

!----------------------------------------------------------------------!
!                                                                      !
!                                                                      !
!                                                                      !
!----------------------------------------------------------------------!
subroutine open_ascii_file(fbase, pc, rc, ic)

  implicit none

  integer :: ios, tc
  character(len=*), intent(in   ) :: fbase
  integer,          intent(  out) :: pc, rc, ic

  if(len_trim(fbase) == 0) then
     write(*,*) 'No file name given.'
     stop 1000
  endif

  open(unit=fdata, file=trim(fbase), iostat = ios)
  if(ios /= 0)then
     write(*,*) 'Error opening AMReX particle ASCII file.'
     stop 1001
  endif

  read(fdata,*) pc  ! Number of particles
  read(fdata,*) rc  ! Number of reals
  read(fdata,*) ic  ! Number of integers
  read(fdata,*) tc  ! Not sure
  read(fdata,*) tc  ! Not sure

  rc = rc + 3 ! Position not included in count
  ic = ic + 2 ! ID and CPU not included in count

  if(rc /= 17 .or. ic /= 4) then
     write(*,*) 'Unknown file format!'
     write(*,*) '    Real count:', rc, '  Expected 17.'
     write(*,*) ' Integer count:', ic, '  Expected  4.'
  endif

  if(debug) then
     write(*,*) 'Particle count:', pc
     write(*,*) '    Real count:', rc
     write(*,*) ' Integer count:', ic
  endif

end subroutine open_ascii_file


!----------------------------------------------------------------------!
!                                                                      !
!                                                                      !
!                                                                      !
!----------------------------------------------------------------------!
subroutine read_ascii_file(pc, rc, ic, rdata, idata)

  implicit none

  integer,          intent(in   ) :: pc, rc, ic

  double precision, intent(  out), allocatable :: rdata(:,:)
  integer         , intent(  out), allocatable :: idata(:,:)

  character(len=2048) :: buff
  integer :: lc, is, ie, p
  integer :: ios

  call allocate_rdata(pc,rc,rdata)
  call allocate_idata(pc,ic,idata)

  do p=1,pc
     read(fdata,"(A)",iostat=ios)buff

     is=1
     ie = scan(buff, " ") -1

     do lc=1, rc
        read(buff(is:ie),*,iostat=ios) rdata(p,lc)
        is = ie + 2
        ie = is + scan(buff(is:), " ") - 2
     enddo
     do lc=1, ic
        read(buff(is:ie),*,iostat=ios) idata(p,lc)
        is = ie + 2
        ie = is + scan(buff(is:), " ") - 2
     enddo
  enddo
  close(fdata)
end subroutine read_ascii_file


!----------------------------------------------------------------------!
!                                                                      !
!                                                                      !
!                                                                      !
!----------------------------------------------------------------------!
subroutine allocate_rdata(pc, rc, rdata)

  implicit none

  integer,          intent(in   ) :: pc, rc

  double precision, intent(  out), allocatable :: rdata(:,:)

  if(allocated(rdata)) then
     if(size(rdata,1) /= pc .or. size(rdata,2) /= rc) then
        write(*,*) 'Resizing rdata.'
        deallocate(rdata)
        allocate(rdata(pc,rc))
     endif
  else
     allocate(rdata(pc,rc))
  endif

end subroutine allocate_rdata

!----------------------------------------------------------------------!
!                                                                      !
!                                                                      !
!                                                                      !
!----------------------------------------------------------------------!
subroutine allocate_idata(pc, ic, idata)

  implicit none

  integer, intent(in   ) :: pc, ic

  integer, intent(  out), allocatable :: idata(:,:)

  if(allocated(idata)) then
     if(size(idata,1) /= pc .or. size(idata,2) /= ic) then
        write(*,*) 'Resizing idata.'
        deallocate(idata)
        allocate(idata(pc,ic))
     endif
  else
     allocate(idata(pc,ic))
  endif

end subroutine allocate_idata

end module ascii_reader
