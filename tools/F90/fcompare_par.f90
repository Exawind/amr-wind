! Program to compare particles data written in ASCII file
! by AMReX
! Usage:
!      fcompare_par  --file1 <file1> --file2 <file2> --nreals <NR> --nints <NI>
!
program fcompare_par

   implicit none

   integer,        parameter  :: dp = selected_real_kind (2*precision (1.0))
   integer,        parameter  :: nr_min = 3  ! Minimum number of reals  (position)
   integer,        parameter  :: ni_min = 2  ! Minimum number of ints   (cpu and id)
   integer,        parameter  :: max_char_len = 5000
   real(dp),       parameter  :: tol = 10.0 ** ( - precision(1.0_dp) )

   ! This is the data type holding the particle infos
   type particle_t
      real (dp),  allocatable :: rdata (:)
      integer,    allocatable :: idata (:)
   end type particle_t

   integer                         :: nr=0, ni=0, np, nfails
   integer                         :: fh1, fh2, p, tmp
   integer                         :: comparison_failed
   integer,           allocatable  :: fails(:)
   logical                         :: brief = .false.
   type (particle_t)               :: single_particle
   type (particle_t), allocatable  :: particles1 (:), particles2 (:)
   character(:),      allocatable  :: file1, file2

   call read_inputs ()
   call check_inputs ()
   call print_inputs ()

   open ( newunit = fh1, file = file1 )
   open ( newunit = fh2, file = file2 )

   ! First element of both files is number of particles
   read (fh1,*) np
   read (fh2,*) tmp

   call check ( ( np == tmp ), "number of particles in file1 and file2 differs" )

   ! Allocate types
   allocate ( particles1 (np), particles2 (np), fails (np) )
   call alloc_particle ( single_particle ) ! this is used as tmp during reads

   ! Read the data into the types
   do p = 1, np
      call read_particle_data ( particles1, fh1 )
      call read_particle_data ( particles2, fh2 )
   end do

   ! Do the actual comparison
   nfails = 0
   comparison_failed = 0
   fails  = 0

   do p = 1, np
      if ( all( ( abs (particles1(p) % rdata - particles2(p) % rdata) ) <= tol ) .and. &
          all( ( particles1(p) % idata == particles2(p) % idata ) ) ) cycle
      nfails = nfails + 1
      fails (nfails) = p
      write (*,'(2X,A,I0)') "Comparison failed for particle ID ", p
   end do

   close ( fh1 )
   close ( fh2 )

   if (  nfails > 0  ) then
      if ( .not. brief ) &
         call print_diff ( particles1, particles2, fails(1:nfails) )
      call check ( .false. , " file1 and file2 differ" )
   end if

contains

   subroutine read_inputs ()

      integer           :: i, length
      character (500)    :: val1, val2

      do i = 1, command_argument_count (), 2
         call get_command_argument ( i, val1, length  )
         call check ( length <= max_char_len, &
              & "Argument length exceeds max char length" )
         call get_command_argument ( i+1, val2, length  )
         call check ( length <= max_char_len, &
              & "Command length exceeds max char length" )
         call set_inputs ( val1, val2 )
      end do

   end subroutine read_inputs

   subroutine set_inputs ( arg_name, arg_value )

      character (*), intent (in) :: arg_name, arg_value

      select case (trim (arg_name))
      case ( "--file1", "-f1" )
         file1 = trim (arg_value)
      case ( "--file2", "-f2" )
         file2 = trim (arg_value)
      case ( "--nreals", "-r" )
         read ( arg_value,' (I8)' ) nr
      case ( "--nints", "-i" )
         read ( arg_value,' (I8)' ) ni
      case ( "--brief" )
         brief = .true.
      case default
         call check ( .false., "Option "//trim (arg_name)//" not recognized")
      end select

   end subroutine set_inputs

   subroutine print_inputs ()

      write (*,"(/,A/)")  repeat ("<",36) // " fcompare_par " // repeat (">",36)
      write (*,"(3X,A)")      "file 1 = "//file1
      write (*,"(3X,A)")      "file 2 = "//file2
      write (*,"(3X,A,I0)")   "nreals = ", nr
      write (*,"(3X,A,I0/)")  "nints  = ", ni

   end subroutine print_inputs

   subroutine check_inputs

      logical :: ex

      ! Check that file 1 has been provided and exists
      call check ( allocated (file1), "file1 has not been provided" )
      inquire ( file = file1, exist = ex )
      call check ( ex, "file1 = "//file1//" does not exist")

      ! Check that file 2 has been provided and exists
      call check ( allocated (file2), "file2 has not been provided" )
      inquire ( file = file2, exist = ex )
      call check ( ex, "file2 = "//file2//" does not exist")

      ! Check that nreals and nints are non-negative
      call check ( ni >= 0, "nints must be non-negative")
      call check ( nr >= 0, "nreals must be non-negative")

   end subroutine check_inputs

   subroutine alloc_particle ( a_particle )

      type (particle_t), intent (out) :: a_particle

      allocate ( a_particle % rdata ( nr + nr_min ) )
      allocate ( a_particle % idata ( ni + ni_min ) )

   end subroutine alloc_particle

   subroutine read_particle_data ( particles, fh )

      type (particle_t), intent (inout) :: particles (:)
      integer,           intent (in)    :: fh
      character(max_char_len)           :: record
      integer                           :: is, ie, ios
      integer                           :: n_r, n_i

      read (fh,' (A)', iostat = ios) record
      call check (  (ios==0), "cannot read input file " )

      n_r = 0
      n_i = 0
      is  = 0
      ie  = 0

      do

         if ( ( n_i == ni + ni_min ) .and. (n_r == nr + nr_min ) ) exit

         is = find_next_start ( record, ie + 1 )
         ie = find_next_end ( record, is )

         if ( n_r < ( nr + nr_min ) ) then
            n_r = n_r + 1
            read ( record (is:ie), * ) single_particle % rdata (n_r)
         else
            n_i = n_i + 1
            read ( record (is:ie), * ) single_particle % idata (n_i)
         end if

      end do

      call alloc_particle ( particles ( single_particle % idata (1) ) )
      particles( single_particle % idata (1) ) = single_particle

   end subroutine read_particle_data

   function find_next_start ( record, istart )  result ( next_start )

      character(*), intent(in) :: record
      integer,      intent(in) :: istart
      integer                  :: next_start, i, ascii

      i = istart

      do
         ascii = iachar ( record(i:i) )
         if ( ( ascii >= iachar ("0") ) .and. ( ascii <= iachar("9") ) .or. &
             ( ascii == iachar ("-") ) ) then
            next_start = i
            exit
         else
            i = i + 1
         end if
      end do

   end function find_next_start

   function find_next_end ( record, istart ) result ( next_end )

      character(*), intent(in) :: record
      integer,      intent(in) :: istart
      integer                  :: next_end
      integer                  :: ascii, i
      character(1), parameter  :: HT = achar(9)

      i = istart

      do
         ascii = iachar ( record(i:i) )
         if ( ( ascii == iachar (" ") ) .or. ( ascii == iachar(HT) ) ) then
            next_end = i - 1
            exit
         else
            i = i + 1
         end if
      end do

   end function find_next_end

   subroutine print_diff ( p1, p2, ids )

      type(particle_t), intent(in) :: p1(:), p2(:)
      integer,          intent(in) :: ids(:)
      integer                      :: i, n, p
      real(dp)                     :: diff, diff_pc

      do n = 1, size(ids)

         p = ids(n)

         write(*,'(/,A,I0,A)')   repeat ("+",32) // " Diff for particle ",p, &
              & " "//repeat ("+",32)

         write(*,'(/,4(A19),/)') "File1", "File2", "Diff", "  % Rel Diff"

         do i = 1, nr_min + nr
            diff    = p1(p) % rdata(i) - p2(p) % rdata(i)
            diff_pc = diff / ( p1(p) % rdata(i) + epsilon ( p1(P) % rdata(i) ) )
            diff_pc = diff_pc * 100.0_dp
            write(*,'(4(es20.6))') p1(p) % rdata(i), p2(p) % rdata(i), diff, &
                 & diff_pc
         end do

         do i = 1, ni_min + ni
            diff    = p1(p) % idata(i) - p2(p) % idata(i)
            diff_pc = diff / ( p1(p) % idata(i) + epsilon ( 0.0_dp ) )
            diff_pc = diff_pc * 100.0_dp
            write(*,'(3(I20),es20.6)') p1(p) % idata(i), p2(p) % idata(i), &
                 & int(diff), diff_pc
         end do

         write(*,*)

      end do

   end subroutine print_diff

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

end program fcompare_par
