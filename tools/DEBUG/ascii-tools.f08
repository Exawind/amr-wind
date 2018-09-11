!----------------------------------------------------------------------!
!                                                                      !
!                                                                      !
!                                                                      !
!----------------------------------------------------------------------!
program main

  use ascii_reader

  implicit none

  integer :: pc, rc, ic

  character(len=512) :: fbase=' '

  double precision, allocatable :: rdata(:,:)
  integer,          allocatable :: idata(:,:)

  integer :: sort = 10

  logical :: group_ids   = .false.
  logical :: get_max_vel = .false.

  call read_inputs()

  call open_ascii_file(fbase, pc, rc, ic)

  call read_ascii_file(pc, rc, ic, rdata, idata)

  if(group_ids) call group_pids(pc, rc, ic, rdata, idata)
  if(get_max_vel) call max_pvel(pc, rc, ic, rdata, idata)

contains

!-----------------------------------------------------------------------!
!                                                                       !
!                                                                       !
!                                                                       !
!-----------------------------------------------------------------------!
subroutine read_inputs ()

  implicit none

  integer        :: length
  character(512) :: val1, val2

  integer :: lc, nargs, ios, cbrt

  nargs = command_argument_count()

  lc = 1
  do while(lc <= nargs)

     call get_command_argument(lc, val1, length); lc = lc+1

     select case (trim (val1))
     case ( "--file", "-f" );
        call get_command_argument (lc, val2, length  ); lc = lc+1
        fbase = trim (val2)

     case ( "--ids" );
        group_ids = .true.

     case ( "--max-vel" );
        get_max_vel = .true.

     case ( "--sort", "-s" );
        call get_command_argument (lc, val2, length  ); lc = lc+1
        read(val2,*) sort

     case default
        write(*,*) "Option "//trim(val1)//" not recognized"
        stop 2299
     end select

  end do

end subroutine read_inputs

!-----------------------------------------------------------------------!
!                                                                       !
!                                                                       !
!                                                                       !
!-----------------------------------------------------------------------!
  subroutine group_pids(pc, rc, ic, rdata, idata)

    implicit none

    integer,          intent(in   ) :: pc, rc, ic
    double precision, intent(in   ) :: rdata(pc,rc)
    integer,          intent(in   ) :: idata(pc,ic)

    integer :: lc, id, min_id, max_id
    integer, allocatable :: ids(:)

    max_id = -1
    min_id = pc

    do lc=1, pc
       id = idata(lc,1)
       max_id = max(max_id, id)
       min_id = min(min_id, id)
    enddo

    allocate(ids(min_id:max_id))
    ids = 0

    do lc=1,pc
       id = idata(lc,1)
       ids(id) = ids(id) + 1
    enddo

    do lc=min_id,max_id
       write(*,*) lc, ids(lc)
    enddo

    deallocate(ids)

  end subroutine group_pids


!-----------------------------------------------------------------------!
!                                                                       !
!                                                                       !
!                                                                       !
!-----------------------------------------------------------------------!
  subroutine max_pvel(pc, rc, ic, rdata, idata)

    implicit none

    integer,          intent(in   ) :: pc, rc, ic
    double precision, intent(in   ) :: rdata(pc,rc)
    integer,          intent(in   ) :: idata(pc,ic)

    integer :: lc1, lc2
    double precision :: vel(3), vel_mag

    integer :: id(sort), cpu(sort)
    double precision :: max_vel(sort), rad(sort)
    double precision :: min_vel

    integer :: tid, tcpu
    double precision :: tvel, trad




    max_vel = -1.0
    min_vel =  1.0

    do lc1=1, pc
       vel = rdata(lc1,9:11)
       vel_mag =  sqrt(dot_product(vel,vel))
       if(vel_mag > max_vel(sort)) then

          max_vel(sort) = vel_mag
          id(sort)      = idata(lc1,1)
          cpu(sort)     = idata(lc1,2)
          rad(sort)     = rdata(lc1,4)

          lc2 = sort-1
          do while(max_vel(lc2) < max_vel(lc2+1) .and. lc2 > 0)

             tvel = max_vel(lc2+1)
             tid  = id(lc2+1)
             tcpu = cpu(lc2+1)
             trad = rad(lc2+1)

             max_vel(lc2+1) = max_vel(lc2)
             id(lc2+1)      = id(lc2)
             cpu(lc2+1)     = cpu(lc2)
             rad(lc2+1)     = rad(lc2)

             max_vel(lc2)   = tvel
             id(lc2)        = tid
             cpu(lc2)       = tcpu
             rad(lc2)       = trad

             lc2 = lc2-1
          enddo
       endif
       min_vel = min(vel_mag, min_vel)

    enddo

    ! write(*,*)'>>>>',min_vel, max_vel(1)
    do lc1=1, sort
       write(*,*) id(lc1), cpu(lc1), rad(lc1), max_vel(lc1)
    enddo
  end subroutine max_pvel

end program main
