MODULE REMOVE_COMMENT_MODULE
CONTAINS
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
END MODULE REMOVE_COMMENT_MODULE
