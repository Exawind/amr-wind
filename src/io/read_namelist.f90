MODULE read_namelist_module

   integer, private :: argc = 0
   character(len=80), private :: argv(32)

   contains
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!     Module name: READ_NAMELIST(POST)                                 !
!     Author: P. Nicoletti                            Date: 25-NOV-91  !
!                                                                      !
!     Purpose: Read in the NAMELIST variables                          !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE READ_NAMELIST()

      use bc
      use ic, only: ic_pack_type, ic_p
      use ic, only: ic_u, ic_v, ic_w
      use ic, only: ic_x_e, ic_x_w, ic_y_n, ic_y_s, ic_z_b, ic_z_t
      use utilities, only: blank_line, seek_comment
      use utilities, only: make_upper_case, replace_tab
      use param, only: undefined

      use remove_comment_module, only: remove_comment
      use remove_comment_module, only: remove_par_blanks

      implicit none

! Local Variables:
!------------------------------------------------------------------------//
      integer, parameter :: unit_dat = 51
! LINE_STRING(1:MAXCOL) has valid input data
      integer, PARAMETER :: MAXCOL = 80
! Holds one line in the input file
      CHARACTER(LEN=512) :: LINE_STRING
! Length of noncomment string
      integer :: LINE_LEN
! Line number
      integer :: LINE_NO
! Indicate whether to do a namelist read on the line
      logical :: READ_FLAG
! Logical to check if file exits.
      logical :: lEXISTS
! Error flag
      logical :: ERROR

      CHARACTER(len=256) :: STRING
      integer :: IOS

! Flags restricting what data from the incflo.dat to process

      READ_FLAG = .TRUE.
      LINE_NO = 0

      ! Open the incflo.dat file. Report errors if the file is not located or
      ! there is difficulties opening it.
      inquire(file='incflo.dat',exist=lEXISTS)
      IF(.NOT.lEXISTS) THEN
         WRITE(*,1000)
         stop 20010

 1000 FORMAT(2/,1X,70('*')/' From: READ_NAMELIST',/' Error 1000: ',    &
         'The input data file, incflo.dat, is missing. Aborting.',/1x,   &
         70('*'),2/)

      ELSE
         OPEN(UNIT=UNIT_DAT, FILE='incflo.dat', STATUS='OLD', IOSTAT=IOS)
         IF(IOS /= 0) THEN
            WRITE (*,1100)
            stop 20011
         ENDIF
      ENDIF

 1100 FORMAT(//1X,70('*')/1x,'From: READ_NAMELIST',/1x,'Error 1100: ', &
         'Line ',A,' in file incflo.dat is too long. Input lines should', &
         /1x,'not pass column ',A,'.',2/3x,A,2/1x,&
         'Please correct input deck.',/1X,70('*'),2/)


     ! Loop through the incflo.dat file and process the input data.
      READ_LP: DO
         READ (UNIT_DAT,"(A)",IOSTAT=IOS) LINE_STRING
         IF(IOS < 0) EXIT READ_LP

         LINE_NO = LINE_NO + 1

         LINE_LEN = SEEK_COMMENT(LINE_STRING,LEN(LINE_STRING)) - 1
         call remove_comment(LINE_STRING, LINE_LEN+1, LEN(LINE_STRING))

         IF(LINE_LEN <= 0) CYCLE READ_LP           ! comment line
         IF(BLANK_LINE(LINE_STRING)) CYCLE READ_LP ! blank line
  
         CALL SET_KEYWORD(ERROR)

      ENDDO READ_LP

      CLOSE(UNIT=UNIT_DAT)

      RETURN

      CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: SET_KEYWORD(ERROR)                                       !
! Author: P. Nicoletti                                Date: 25-NOV-91  !
!                                                                      !
! Purpose: Process LINE_STRING for incflo keyword data.                !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE SET_KEYWORD(ERROR)

      IMPLICIT NONE

      logical, intent(OUT) ::ERROR

! External namelist files:
!---------------------------------------------------------------------//
      include 'geometry.inc'
      include 'initial_conditions.inc'
      include 'boundary_conditions.inc'

      ERROR = .FALSE.

      CALL MAKE_UPPER_CASE (LINE_STRING, LINE_LEN)
      CALL REPLACE_TAB (LINE_STRING, LINE_LEN)
      call remove_par_blanks(LINE_STRING)

! Write the current line to a scratch file
! and read the scratch file in NAMELIST format
      IF(.NOT.READ_FLAG) RETURN

! Geometry and discretization keywords
      STRING=''; STRING = '&GEOMETRY '//&
         trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
      READ(STRING, NML=GEOMETRY, IOSTAT=IOS)
      IF(IOS == 0)  RETURN


! Initial condtion keywords
      STRING=''; STRING = '&INITIAL_CONDITIONS '//&
         trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
      READ(STRING, NML=INITIAL_CONDITIONS, IOSTAT=IOS)
      IF(IOS == 0)  RETURN


! Boundary condition keywords
      STRING=''; STRING = '&BOUNDARY_CONDITIONS '//&
         trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
      READ(STRING, NML=BOUNDARY_CONDITIONS, IOSTAT=IOS)
      IF(IOS == 0)  RETURN

       ERROR = .TRUE.

      RETURN
      END SUBROUTINE SET_KEYWORD

END SUBROUTINE READ_NAMELIST


subroutine add_argument(fname, nlen) &
   bind(c,name='incflo_add_argument')

   use iso_c_binding, only: c_int, c_float, c_char
   implicit none
   integer(c_int), intent(in) :: nlen
   character(kind=c_char), intent(in) :: fname(nlen)
   integer :: lc, bnd

   argc = argc + 1
   if(argc > 32) then
      write(*,*) 'from add_argument:'
      write(*,*) 'too many arguments!!!!'
      stop 986532
   endif

! Copy over the string, character-by-character.
   bnd = min(nlen,80)
   do lc=1, bnd
      argv(argc)(lc:lc) = fname(lc)
   enddo
! Clear out the remaining string.
   do lc=bnd+1,80
      argv(argc)(lc:lc) = ' '
   enddo

end subroutine add_argument


END MODULE read_namelist_module
