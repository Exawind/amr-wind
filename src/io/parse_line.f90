module parse_line_module


      use param, only: one


   contains

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: PARSE_LINE(LINE, LMAX,           READ_FLAG)            C
!  Author: M. Syamlal                                 Date: 27-JUN-97  C
!                                                                      C
!  Purpose: Parse input line                                           C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      subroutine PARSE_LINE(LINE, LMAX, READ_FLAG)

      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------

! Line to be parsed.
      character(len=*), intent(IN) :: LINE

! Length of of LINE.
      integer, intent(IN) :: LMAX

! Indicate whether to do a namelist read on the line. A namelist read
! is still preformed when an arithmetic operation is found.
      logical, intent(OUT) :: READ_FLAG

! Start and end locations for the search parameters.
      integer LSTART, LEND

! The string is empty. No need to parse.
      IF (LMAX == 0) THEN
         READ_FLAG = .FALSE.
         RETURN
      ENDIF

! Check to see if the input line contains '@('. If this string is found,
! then the line contains information to parsed. This string indicates
! that one of two actions need to occur"
! 1) there is an expression to evalute; @(6.0/2.0) = 3.0
      LSTART = INDEX(LINE,'@(')

! If the returned index is not zero, the input line contain the string.
      IF (LSTART /= 0) THEN

! Look for the ending parenthesis. If none exists, flag the error and
! exit incflo.
         LEND = LSTART - 1 + INDEX(LINE(LSTART:LMAX),')')
         IF (LEND <= LSTART) THEN
            WRITE (*, 1000) 0,LINE(LSTART:LMAX)
            stop 20015
         ENDIF
      ENDIF ! IF (LSTART /= 0) THEN
!
      LSTART = INDEX(LINE,'@(')             !Arithmetic processing ?
      IF (LSTART /= 0) CALL PARSE_ARITH (LINE, LMAX)
      READ_FLAG = .TRUE.
!
      RETURN

 1000 FORMAT(//1X,70('*')/' (PE ',I6,'): From: PARSE_LINE',/&
         ' Message: An evaluation statement "@(" was found in ',&
         'the input line,',/' but no ending parenthesis was located:',/&
         ' INPUT: ',A,/1X,70('*')//)

      END subroutine PARSE_LINE


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: PARSE_ARITH(LINE, LMAX)                                C
!  Author: M. Syamlal                                 Date: 10-AUG-92  C
!                                                                      C
!  Purpose: Complete arithmetic operations and expand the line         C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      subroutine PARSE_ARITH(LINE, LMAX)

      use amrex_fort_module, only : rt => amrex_real
      use iso_c_binding , only: c_int

      use utilities, ONLY: seek_end

      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!C
!                      The part of LINE containing input
      integer LMAX
!                      Input line with arithmetic operations.  Out put
!                      line with completed arithmetic statements.
!
      character LINE*(*)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!
!                      Value of pi
      real(rt) PI
!
!     Cumulative value and sub value
      real(rt) VALUE, SUB_VALUE
!
!                      Start and end locations for the arithmetic operation
      integer          LSTART, LEND
!
!                      Length of arithmetic operation string
      integer          LENGTH
!
!                      22 - LENGTH
      integer          LDIF
!
!                      Locations in SUB_STR, and LINE
      integer          LSUB, L
!
!                      Operator symbol (Legal values: *, /)
      character(LEN=1) :: OPERATION
!
!                      Substring taken from LINE
      character(LEN=80) :: SUB_STR
!
!-----------------------------------------------
!
!
      PI = 4.0D0*ATAN(ONE)
!
!  Search for arithmetic operation
!
   10 CONTINUE
      LMAX = SEEK_END(LINE,LEN(LINE))
!
      LSTART = INDEX(LINE,'@(')
!
      IF (LSTART == 0) RETURN
!
      LEND = LSTART - 1 + INDEX(LINE(LSTART:LMAX),')')
      IF (LEND <= LSTART) THEN
         WRITE (*, 1000) 0,LINE(LSTART:LMAX)
         stop 20016
      ENDIF
!
!    Do the arithmetic
!
      VALUE = ONE
      OPERATION = '*'
      LSUB = 1
      DO L = LSTART + 2, LEND
         IF (LINE(L:L)=='*' .OR. LINE(L:L)=='/' .OR. LINE(L:L)==')') THEN
            IF (LSUB == 1) THEN
               WRITE (*, 1015) 0,LINE(LSTART:LEND)
               stop 20017
            ENDIF
            IF (SUB_STR(1:LSUB-1) == 'PI') THEN
               SUB_VALUE = PI
            ELSE
               READ (SUB_STR(1:LSUB-1), *, ERR=900) SUB_VALUE
            ENDIF
            IF (OPERATION == '*') THEN
               VALUE = VALUE*SUB_VALUE
            ELSE IF (OPERATION == '/') THEN
               VALUE = VALUE/SUB_VALUE
            ENDIF
            LSUB = 1
            OPERATION = LINE(L:L)
         ELSE IF (LINE(L:L) == ' ') THEN
         ELSE
            SUB_STR(LSUB:LSUB) = LINE(L:L)
            LSUB = LSUB + 1
         ENDIF
      END DO
      LENGTH = LEND - LSTART + 1
      IF (LENGTH > 22) THEN
         DO L = LSTART + 22, LEND
            LINE(L:L) = ' '
         END DO
      ELSE IF (LENGTH < 22) THEN
         LMAX = SEEK_END(LINE,LEN(LINE))
         LDIF = 22 - LENGTH
         IF (LMAX + LDIF > LEN(LINE)) THEN
            WRITE (*, 1020) 0,LINE(1:80)
            stop 20018
         ENDIF
         DO L = LMAX, LEND + 1, -1
            LINE(L+LDIF:L+LDIF) = LINE(L:L)
         END DO
      ENDIF
!
!  Transfer the value to LINE
!
      WRITE (SUB_STR, '(G22.15)') VALUE
      L = LSTART
      DO LSUB = 1, 22
         LINE(L:L) = SUB_STR(LSUB:LSUB)
         L = L + 1
      END DO
      GO TO 10
!
  900 CONTINUE
      WRITE (*, 1010) 0, SUB_STR(1:LSUB-1)
      stop 20019
 1000 FORMAT(/1X,70('*')//'(PE ',I6,'): From: PARSE_ARITH',/&
         ' Message: No ending ) found in the input line: ',/9X,A,/1X,70('*')/)
 1010 FORMAT(/1X,70('*')//'(PE ',I6,'): From: PARSE_ARITH',/&
         ' Message: Error reading the input string: ',/9X,A,/1X,70('*')/)
 1015 FORMAT(/1X,70('*')//'(PE ',I6,'): From: PARSE_ARITH',/&
         ' Message: Invalid operator in the input string: ',/9X,A,/1X,70('*')/)
 1020 FORMAT(/1X,70('*')//'(PE ',I6,'): From: PARSE_ARITH',/&
         ' Message: Too many arithmetic operations in the line: ',/1X,A,/1X,70(&
         '*')/)
      end subroutine PARSE_ARITH

   end module parse_line_module
