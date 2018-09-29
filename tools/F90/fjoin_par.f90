! Program to compare particles data written in ASCII file  by AMReX
!
! Usage:
!
!      fjoin_par  --file STRING --int INTEGER --id INTEGER  \
!         --var INTEGER --dt REAL --verbose
!
!  --file    The base name used to generate ASCII particle output
!            files. This is the same as `amr.par_ascii_file` in
!            the input deck. <par>
!
!  --start
!  --end
!
!  --id      ID of particle for data extraction.
!
!  --dt      This is the simulation dt separating the ascii output
!            files. When this is passed, the output tries to reconstruct
!            the simulation output time.
!
!  --var     Index of particle property for extraction. This input
!            may be passed multiple times.
!                1 - position-x
!                2 - position-y
!                3 - position-z
!                4 - radius
!                5 - volume
!                6 - mass
!                7 - density
!                8 - oneOverI
!                9 - velocity-x
!               10 - velocity-y
!               11 - velocity-z
!               12 - omega-x
!               13 - omega-y
!               14 - omega-z
!               15 - drag-x
!               16 - drag-y
!               17 - drag-z
!
!              100 - kinetic energy
!

program fjoin_par

   implicit none

   integer, parameter :: dp = selected_real_kind (2*precision (1.0))
   integer, parameter :: nr_min = 3  ! Minimum number of reals
   integer, parameter :: ni_min = 2  ! Minimum number of ints
   integer, parameter :: max_char_len = 5000
   integer, parameter :: max_var = 10

   character(:), allocatable  :: fbase
   integer :: id  = -1
   integer :: interval =  1
   integer :: istart =    0
   integer :: iend   =    0
   integer :: iformat =  10
   integer :: var(10) =  -1
   real (dp) :: dt =   -1.0

   logical :: verbose = .false.

   ! This is the data type holding the particle infos
   type particle_t
      real (dp),  allocatable :: rdata (:)
      integer,    allocatable :: idata (:)
   end type particle_t

   type (particle_t), allocatable :: particles(:)
   integer, parameter :: funit = 1000
   integer :: tmpi
   integer :: nr, ni, np, nro, nio, npo
   integer :: lc1, lc2, lc3, lc4
   character(len=5) :: clc2
   logical :: fexist
   integer :: var_count = 0
   character(len=1000) :: output
   real (dp) :: val(10)

   character(len=32) ::  io_format
   integer :: fcount = 0
   integer :: err = 0

   call read_inputs ()
   call check_inputs ()
   if(verbose) call print_inputs ()

   fcount = 0
   lc2=0
   do
      write(clc2,"(i5.5)") istart + lc2
      inquire(file = trim(fbase)//clc2, exist=fexist)
      if(fexist) then

         if(verbose) write(*,*)'inspecting --> ',trim(fbase)//clc2
         open ( unit=funit, file=trim(fbase)//clc2 )
         read (funit,*) np
         read (funit,*) nr
         read (funit,*) ni
         close ( funit )

         if(lc2 > 0) then
            if(npo /= np) err=1
            if(nro /= nr) err=1
            if(nio /= ni) err=1
         endif
         if(err == 0 ) fcount = fcount + 1
      endif
      if(istart + lc2 >= iend .or. err /=0) exit

      npo = np
      nro = nr
      nio = ni

      lc2 = lc2 + 1
   enddo

   if(fcount == 0)then
      write(*,*) 'No files detected.'
      stop 2
   endif

   if(verbose) then
      write(*,*) 'Number of files:    ',fcount
      write(*,*) 'Number of particles:',np
      write(*,*) 'Number of reals:    ',nr
      write(*,*) 'Number of integers: ',ni
   endif

   call alloc_particles

   io_format=''
   if(iformat < 10) then
      write(io_format,"('(f24.',I1,',1x)')") iformat
   else
      write(io_format,"('(f24.',I2,',1x)')") iformat
   endif

   lc1=0
   lc2=0
   do
      write(clc2,"(i5.5)") istart + lc2
      inquire(file = trim(fbase)//clc2, exist=fexist)
      if(fexist) then

         if(verbose) write(*,*)'reading    --> ',trim(fbase)//clc2
         open ( unit=funit, file=trim(fbase)//clc2 )
         call read_particle_data
         close ( funit )
         lc1 = lc1 + 1

         output = ''
         if(dt /= -1.0) write(*,"(3x,f15.6)",advance='no') dt*lc2

         do lc4=1,var_count
            call write_var(var(lc4))
         enddo

      endif
      if(istart + lc2 >= iend) exit
      lc2 = lc2 + 1
   enddo

contains

!-----------------------------------------------------------------------!
!                                                                       !
!                                                                       !
!                                                                       !
!-----------------------------------------------------------------------!
   subroutine alloc_particles

      integer :: lc1

      allocate( particles(np))

      do lc1=1,np
         allocate ( particles(lc1) % rdata ( nr + nr_min ) )
         allocate ( particles(lc1) % idata ( ni + ni_min ) )
      enddo

   end subroutine alloc_particles

!-----------------------------------------------------------------------!
!                                                                       !
!                                                                       !
!                                                                       !
!-----------------------------------------------------------------------!
   subroutine read_particle_data

      implicit none

      character(max_char_len)           :: record
      integer                           :: i, is, ie, ios
      integer :: lc1
      logical :: lverbose

      lverbose = verbose .and. .false.

      read (funit,*)
      read (funit,*)
      read (funit,*)
      read (funit,*)
      read (funit,*)

      do lc1=1,np

         read (funit,' (A)', iostat = ios) record
         call check (  (ios==0), "cannot read input file " )

         is = 1
         ie = scan (record, " ") - 1
         do i = 1, nr + nr_min
            read ( record (is:ie), * ) particles(lc1) % rdata (i)
            if(lverbose) write(*,*) i, particles(lc1) % rdata (i)
            is = ie + 2 ! 2 = 1 space + 1 next character
            ie = is + scan ( record (is:), " " ) - 2
         end do
         do i = 1, ni + ni_min
            read ( record (is:ie), * ) particles(lc1) % idata (i)
            if(lverbose) write(*,*) i, particles(lc1) % idata (i)
            is = ie + 2 ! 2 = 1 space + 1 next character
            ie = is + scan ( record (is:), " " ) - 2
         end do
      enddo

   end subroutine read_particle_data

!-----------------------------------------------------------------------!
!                                                                       !
!                                                                       !
!                                                                       !
!-----------------------------------------------------------------------!
   subroutine check ( condition, msg )

      logical,      intent (in) :: condition
      character (*), intent (in) :: msg

      if  ( .not. condition ) then
         write (*,*)
         write (*,*) "ERROR: "//msg
         write (*,*) "STOP"
         error stop 1
      end if

   end subroutine check

!-----------------------------------------------------------------------!
!                                                                       !
!                                                                       !
!                                                                       !
!-----------------------------------------------------------------------!
   subroutine write_var(llc)

      implicit none

      integer, intent(in) :: llc
      integer :: lc
      real(dp) :: tmp

      if(llc == 100) then
         write(*,trim(io_format),advance='no') &
            clean_value(calc_granular_temperature())
      else
         if(id == -1) then
            do lc=1,np
               write(*,trim(io_format),advance='no') &
                  clean_value(particles(lc) % rdata(llc))
            enddo
         else
            write(*,trim(io_format),advance='no') &
               clean_value(particles(id) % rdata(llc))
         endif
      endif
      write(*,*)' '

   end subroutine write_var

   real(dp) function clean_value(vv)
      real(dp),  intent(in) :: vv
      clean_value = vv
      if(abs(clean_value) < epsilon(0.0d0)) clean_value = 0.0d0
   end function clean_value

!-----------------------------------------------------------------------!
!                                                                       !
!                                                                       !
!                                                                       !
!-----------------------------------------------------------------------!
   real(dp) function calc_granular_temperature()

      implicit none

      integer :: lc
      real(dp) :: gtmp, tvel

      calc_granular_temperature = 0.0d0
      if(np <= 0) return

      gtmp=0.0d0
      do lc=1,np
         gtmp = gtmp + dot_product( &
                particles(lc) % rdata(9:11), &
                particles(lc) % rdata(9:11) )
      enddo
      calc_granular_temperature = gtmp/(3.0d0*dble(np))

   end function calc_granular_temperature

!-----------------------------------------------------------------------!
!                                                                       !
!                                                                       !
!                                                                       !
!-----------------------------------------------------------------------!
   subroutine read_inputs ()

      integer           :: i, length
      character (500)    :: val1, val2

      do i = 1, command_argument_count (), 2
         call get_command_argument ( i, val1, length  )
         call check ( length <= max_char_len, &
                     "Argument length exceeds max char length" )
         call get_command_argument ( i+1, val2, length  )
         call check ( length <= max_char_len, &
                     "Command length exceeds max char length" )
         call set_inputs ( val1, val2 )
      end do

   end subroutine read_inputs

!-----------------------------------------------------------------------!
!                                                                       !
!                                                                       !
!                                                                       !
!-----------------------------------------------------------------------!
   subroutine set_inputs ( arg_name, arg_value )

      character (*), intent (in) :: arg_name, arg_value

      select case (trim (arg_name))
      case ( "--file", "-f" ); fbase = trim (arg_value)
      case ( "--id"    ); read ( arg_value,*) id
      case ( "--int"   ); read ( arg_value,*) interval
      case ( "--dt"    ); read ( arg_value,*) dt
      case ( "--format"); read ( arg_value,*) iformat
      case ( "--start" ); read ( arg_value,*) istart
      case ( "--end"   ); read ( arg_value,*) iend
      case ( "--verbose" ); verbose = .true.
      case ( "--var" )
         var_count = var_count + 1
         call check ( var_count < max_var, &
                     "Too many variables specified" )
         read ( arg_value,*) var(var_count)
      case default
         call check ( .false., &
                     "Option "//trim (arg_name)//" not recognized")
      end select

   end subroutine set_inputs

!-----------------------------------------------------------------------!
!                                                                       !
!                                                                       !
!                                                                       !
!-----------------------------------------------------------------------!
   subroutine check_inputs

      integer :: lc

      ! Check that file 1 has been provided and exists
      call check ( allocated (fbase), &
                  "Base file name has not been provided" )

      ! Check that nreals and nints are non-negative
      call check ( id > 0 .or. id == -1, "Particle ID must be positive")
      call check (var_count > 0, "No variables specified")
      do lc=1,var_count
         call check (var(lc) > 0, "Variable index must be positive")
      enddo

   end subroutine check_inputs

!-----------------------------------------------------------------------!
!                                                                       !
!                                                                       !
!                                                                       !
!-----------------------------------------------------------------------!
   subroutine print_inputs ()

      integer :: lc

      write (*,"(/,A/)")  repeat ("<",36) // " fjoin_par " // repeat (">",36)
      write (*,"(3X,A)")      "fbase    = "//fbase
      write (*,"(3X,A,I0)")   "ID       = ", id
      if(dt /= -1.0) write (*,"(3X,A,es12.6/)")  "dt       = ", dt
      do lc=1,var_count
         write (*,"(3X,A,I0/)")  "Variable = ", var(lc)
      enddo

   end subroutine print_inputs

end program fjoin_par
