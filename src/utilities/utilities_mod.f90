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

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine Name: REMOVE_COMMENT (LINE, LSTART, MAXCOL)              !
!  Author: P.Nicoletti                                Date: -unknown-  !
!  Reviewer: J.Musser                                 Date: 19-Sept-13 !
!                                                                      !
!  Purpose: Remove comments                                            !
!                                                                      !
!           Example:  IN: > LINE ::  "MW_g( 3 ) =  32.0 ! Oxygen"      !
!                     OUT > LINE ::  "MW_g( 3 ) =  32.0         "      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   SUBROUTINE REMOVE_COMMENT(LINE, LSTART, MAXCOL)

      IMPLICIT NONE

! Passed Variables: Dummy argument format required by ODEPACK.
!---------------------------------------------------------------------//
! Input data line
      CHARACTER(len=*), intent(inout) :: LINE
!Start of comments
      integer, intent(IN) :: LSTART
! Maximum column of input data line to search
      integer, intent(IN) :: MAXCOL

! Local Variables:
!---------------------------------------------------------------------//
! Loop index
      integer :: L

      DO L = LSTART, MAXCOL
         LINE(L:L) = ' '
      END DO

      RETURN
   END SUBROUTINE REMOVE_COMMENT

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine Name: REMOVE_PAR_BLANKS (LINE)                           !
!  Author: J.Musser                                   Date: 19-Spt-13  !
!                                                                      !
!  Purpose: Remove blanks within parentheses. This addition was need   !
!           to resolve portability issues on ALCF machines.            !
!                                                                      !
!           Example:  IN: > LINE :: "MW_g( 3 ) = 32.0"                 !
!                     OUT > LINE :: "MW_g(3)   = 32.0"                 !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   SUBROUTINE REMOVE_PAR_BLANKS(LINE)

      IMPLICIT NONE

! Passed Variables: Dummy argument format required by ODEPACK.
!---------------------------------------------------------------------//
! Input data line
      CHARACTER(len=*), intent(inout) :: LINE

! Local Variables:
!---------------------------------------------------------------------//
! Loop index
      integer :: L
! Search positions
      integer :: POS, lP, rP
! Flag for space replacement
      logical :: searchRight
! Internal search flag
      logical :: replace
! Debug flag
      logical, parameter :: verbose = .FALSE.

! Exit if the string is empty.
      IF(len_trim(LINE) == 0) return

! Get the position of the first left parentheses.
      lP = index(LINE,"(")

! Initialize the loop flag.
      searchRight = (lP /= 0)
      DO WHILE(searchRight)
! Find the position of the first right parentheses.
         rP = lP + index(LINE(lP:),")")
! Check if there are any blank spaces:
         IF(index(LINE(lP:rP-1)," ") /= 0) THEN

            IF(verbose) WRITE(*,"(3X,'Removing spaces: ')")
            IF(verbose) WRITE(*,"(5X,'Before: ',A)") trim(LINE)

! Initialize the loop flag and sub-string position.
            replace = .TRUE.
            POS = lP+1
            DO WHILE(replace)
! If a blank is located, slide all entries to the left one position and
! add a blank to the end of the sub-string.
               IF(LINE(POS:POS) == " ") THEN
                  DO L=POS,rP-2
                     LINE(L:L) = LINE(L+1:L+1)
                  ENDDO
                  LINE(rP-1:rP-1) = " "
               ELSE
! If there the character is not a space, increment character index.
                  POS = POS + 1
               ENDIF
! Exit if all that remains in the sub-string are empty spaces.
               replace = (len_trim(LINE(POS:rP-1)) /= 0)
            ENDDO
            IF(verbose) WRITE(*,"(5X,'After:  ',A)") trim(LINE)
         ENDIF
! Check if there is another set of parentheses.
         lP = rP + index(LINE(rP+1:),"(")
! Exit if no addition parentheses pair is found.
         searchRight = (lP.NE.rP)
      ENDDO

      return
   END SUBROUTINE REMOVE_PAR_BLANKS

END MODULE utilities
