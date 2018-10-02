MODULE utilities

   use amrex_fort_module, only : rt => amrex_real
   use iso_c_binding , only: c_int

   IMPLICIT NONE

CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Function name: SEEK_COMMENT (LINE_MAXCOL)                           !
!  Author: P.Nicoletti                                Date: 25-NOV-91  !
!                                                                      !
!  Purpose: Returns the index to where a comment character was found   !
!  in the input data line.  Equals MAXCOL + 1 if no-comment characters !
!  in the line.                                                        !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   PURE integer FUNCTION SEEK_COMMENT (LINE, MAXCOL)

      IMPLICIT NONE

! Dummy Arguments
!---------------------------------------------------------------------//
! Input data line
      CHARACTER(len=*), intent(IN) :: LINE
! Maximum column of input data line to search
      integer, intent(IN) :: MAXCOL

! Local Variables
!---------------------------------------------------------------------//
! The number of designated comment characters
      integer, PARAMETER :: DIM_COMMENT = 2
! The comment characters
      CHARACTER, PARAMETER :: COMMENT_CHAR(DIM_COMMENT) = (/'#', '!'/)
! Loop indicies
      integer :: L, L2
!.......................................................................!

      DO L = 1, MAXCOL
         DO L2 = 1, DIM_COMMENT
            IF (LINE(L:L) == COMMENT_CHAR(L2)) THEN
               SEEK_COMMENT = L
               RETURN
            ENDIF
         END DO
      END DO
      SEEK_COMMENT = MAXCOL + 1
!
      RETURN
   END FUNCTION SEEK_COMMENT

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Function: BLANK_LINE                                                !
! Author: P. Nicoletti                                Date: 25-NOV-91  !
!                                                                      !
! Purpose: Return .TRUE. if a line contains no input or only spaces.   !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   PURE logical FUNCTION BLANK_LINE (line)

      IMPLICIT NONE

! Dummy Arguments
!---------------------------------------------------------------------//
      CHARACTER, intent(IN) :: LINE*(*)

! Local Variables
!---------------------------------------------------------------------//
      integer :: L
!.......................................................................!

      BLANK_LINE = .FALSE.
      DO L=1, len(line)
         IF(line(L:L)/=' ' .and. line(L:L)/='    ')RETURN
      ENDDO

      BLANK_LINE = .TRUE.
      RETURN
   END FUNCTION BLANK_LINE

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: MAKE_UPPER_CASE (LINE_STRING,MAXCOL)                   C
!  Author: P.Nicoletti                                Date: 26-NOV-91  C
!                                                                      C
!  Purpose: change lowercase characters to uppercase in input line     C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
   SUBROUTINE MAKE_UPPER_CASE(LINE_STRING, MAXCOL)

      IMPLICIT NONE

! Dummy Arguments
!---------------------------------------------------------------------//
! Input line to change to uppercase
      CHARACTER(len=*), intent(inout) :: LINE_STRING
! Number of characters to look at in LINE_STRING
      integer, intent(IN) :: MAXCOL

! Local Variables
!---------------------------------------------------------------------//
! ICHAR value for UPPERCASE A, lowercase a, lowercase z
      integer, PARAMETER :: A_UP = ICHAR('A')
      integer, PARAMETER :: A_LO = ICHAR('a')
      integer, PARAMETER :: Z_LO = ICHAR('z')
! ICHAR differnce between lower and uppercase letters
      integer, PARAMETER :: A_DIFF = A_LO - A_UP
! Holds ICHAR value of current character
      integer :: INT_C
! loop index
      integer :: L
!.......................................................................!

      DO L = 1, MAXCOL
         INT_C = ICHAR(LINE_STRING(L:L))
         IF (A_LO<=INT_C .AND. INT_C<=Z_LO) THEN
            INT_C = INT_C - A_DIFF
            LINE_STRING(L:L) = CHAR(INT_C)
         ENDIF
      END DO
      RETURN
   END SUBROUTINE MAKE_UPPER_CASE

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: REPLACE_TAB (LINE_STRING,MAXCOL)                       !
!  Author: M. Syamlal                                 Date: 10-JUL-03  !
!                                                                      !
!  Purpose: replace tab characters with space                          !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   SUBROUTINE REPLACE_TAB(LINE_STRING, MAXCOL)

      IMPLICIT NONE

! Dummy Arguments
!---------------------------------------------------------------------//
! Input line to change to uppercase
      CHARACTER(len=*), intent(inout) :: LINE_STRING
! Number of characters to look at in LINE_STRING
      integer, intent(IN) :: MAXCOL

! Local Variables
!---------------------------------------------------------------------//
      CHARACTER, PARAMETER :: TAB = CHAR(9)
      CHARACTER, PARAMETER :: CRET = CHAR(13)
! Loop index
      integer :: L
!.......................................................................!

      DO L = 1, MAXCOL
         if(LINE_STRING(L:L) .eq. TAB) LINE_STRING(L:L) = ' '
         if(LINE_STRING(L:L) .eq. CRET) LINE_STRING(L:L) = ' '
      END DO
      RETURN
   END SUBROUTINE REPLACE_TAB

END MODULE utilities
