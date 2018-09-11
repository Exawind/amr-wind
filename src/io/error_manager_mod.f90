!----------------------------------------------------------------------!
! Module: ERROR_MANAGER                                                !
!                                                                      !
! Purpose: Unify error message handeling.                              !
!                                                                      !
!----------------------------------------------------------------------!
      MODULE ERROR_MANAGER

      use amrex_fort_module, only : rt => amrex_real
      use iso_c_binding , only: c_int



      implicit none

! Interface
!---------------------------------------------------------------------//
      interface iVal
         module procedure iVal_int
         module procedure iVal_dbl
         module procedure iVal_log
      end interface

! Maximum number of lines a message can have before a flush is needed.
      integer, PARAMETER :: LINE_COUNT  = 32
! Maximum number of characters per line.
      integer, PARAMETER :: LINE_LENGTH = 256

! Character string for storing the error message.
      CHARACTER(LEN=LINE_LENGTH), DIMENSION(LINE_COUNT) :: ERR_MSG

! Depth that the current call tree can go.
      integer, PARAMETER, PRIVATE :: MAX_CALL_DEPTH = 16
! Current call depth.
      integer, PRIVATE :: CALL_DEPTH

! The name of the calling routine. Set by calling: INIT_ERR_MSG
      CHARACTER(LEN=128), DIMENSION(MAX_CALL_DEPTH), PRIVATE :: CALLERS

! Error Flag.
      integer :: IER_EM

    contains

!``````````````````````````````````````````````````````````````````````!
! Subroutine: INIT_ERROR_MANAGER                                       !
!                                                                      !
! Purpose: Initialize the error manager. This routine also opens the   !
! .LOG file(s) based on user input settings.                           !
!......................................................................!
      SUBROUTINE INIT_ERROR_MANAGER

! Global Variables:
!---------------------------------------------------------------------//

! Undefined character string.
      use param, only: UNDEFINED_C

      implicit none

! Local Variables:
!---------------------------------------------------------------------//
! Log file name.
      CHARACTER(len=255) :: LOGFILE
      CHARACTER(len=255) :: FILE_NAME
! First non-blank character in run_name.
      integer :: NB
! Integer error flag
      integer :: IER(0:0)

! Initizilae the error flags.
      IER = 0
      IER_EM = 0
! Initialize the call tree depth.
      CALL_DEPTH = 0
! Clear the error message storage container.
      ERR_MSG = ''
! Clear the caller routine information.
      CALLERS = ''

! Verify that the .LOG file was successfully opened. Otherwise, flag the
! error and abort.
      IF(sum(IER) /= 0) THEN
         WRITE(*,1001) trim(FILE_NAME)
        stop 20003
      ENDIF

      RETURN

 1001 FORMAT(2/,1X,70('*')/' From: INIT_ERROR_MANAGER',/               &
         ' Error 1001: Failed to open log file: ',A,/' Aborting run.'/,&
         1x,70('*'),2/)

      END SUBROUTINE INIT_ERROR_MANAGER


!``````````````````````````````````````````````````````````````````````!
! Subroutine: INIT_ERR_MSG                                             !
!                                                                      !
! Purpose: Initialize the error manager for the local routine. This    !
! call is needed to set the caller routines name for error messages.   !
!......................................................................!
      SUBROUTINE INIT_ERR_MSG(CALLER)

! Rank ID of process

      implicit none

      CHARACTER(LEN=*), intent(IN) :: CALLER

! Verify that the maximum call dept will not be exceeded.  If so, flag
! the error and exit.
      IF(CALL_DEPTH + 1 > MAX_CALL_DEPTH) THEN
         WRITE(*,1000) CALL_DEPTH
         CALL SHOW_CALL_TREE
         stop 20004
      ELSE
! Store the caller routines name.
         CALL_DEPTH = CALL_DEPTH + 1
         CALLERS(CALL_DEPTH) = trim(CALLER)
      ENDIF

! Clear out the error manager.
      ERR_MSG=''

      RETURN

 1000 FORMAT(/1X,70('*')/' From: ERROR_MANAGER --> INIT_ERR_MSG',/     &
         ' Error 1000: Invalid ERROR_MANAGER usage. The maximum call', &
         ' depth ',/' was exceeded. The calls to INIT_ERR_MSG should', &
         ' have corresponding',/' calls to FINL_ERR_MSG. The current', &
         ' CALL tree depth is: ',I4)

      END SUBROUTINE INIT_ERR_MSG


!``````````````````````````````````````````````````````````````````````!
! Subroutine: FINL_ERR_MSG                                             !
!                                                                      !
! Purpose: Finalize the error manager. The call is needed to clear out !
! old information and unset the lock.                                  !
!......................................................................!
      SUBROUTINE FINL_ERR_MSG

! Rank ID of process

      implicit none

! Single line.
      CHARACTER(LEN=LINE_LENGTH) :: LINE
! Line length with trailing space removed.
      integer :: LENGTH
! Line Counter
      integer :: LC
! Number of non-empty lines.
      integer :: COUNT

! The current calling routine.
      CHARACTER(LEN=128) :: CALLER

! Verify that at the INIT routine was called.
      IF(CALL_DEPTH < 1) THEN
         WRITE(*,1000)
         stop 20005
      ELSE
! Store the current caller, clear the array position, and decrement
! the counter.
         CALLER = CALLERS(CALL_DEPTH)
         CALLERS(CALL_DEPTH) = ''
         CALL_DEPTH = CALL_DEPTH - 1
      ENDIF

! Verify that the error message container is empty.
      COUNT = 0
      DO LC = 1, LINE_COUNT
         LINE = ERR_MSG(LC)
         LENGTH = len_trim(LINE)
         IF(0 < LENGTH .AND. LENGTH < 256 ) COUNT = COUNT + 1
      ENDDO

! If the error message container is not empty, report the error, dump
! the error message and abort incflo.
      IF(COUNT /= 0) THEN
         WRITE(*,1001) trim(CALLER)
! Write out the error message container contents.
         DO LC = 1, LINE_COUNT
            LINE = ERR_MSG(LC)
            LENGTH = len_trim(LINE)
            IF(0 < LENGTH .AND. LENGTH < 256 ) THEN
               WRITE(*,1002)LC, LENGTH, trim(LINE)
            ENDIF
         ENDDO
         WRITE(*,1003)
         stop 20006
      ENDIF

! This shouldn't be needed, but it doesn't hurt.
      ERR_MSG = ''

      RETURN

 1000 FORMAT(/1X,70('*')/' From: ERROR_MANAGER --> FINL_ERR_MSG',/     &
         ' Error 1000: Ivalid ERROR_MANAGER usage. A call to FINL_ERR',&
         '_MSG was',/' made while the call tree is empty. This can',   &
         ' occur if a call to',/' FINL_ERR_MSG was made without a',    &
         ' corresponding call to INIT_ERR_MSG.',/' Aborting incflo.'/  &
         1x,70('*'),2/)

 1001 FORMAT(/1X,70('*')/' From: ERROR_MANAGER --> FINL_ERR_MSG',/     &
         ' Error 1001: Error container ERR_MSG not empty.',/           &
         ' CALLERS: ',A,2/' Contents:')

 1002 FORMAT(' LC ',I2.2,': LEN: ',I3.3,1x,A)

 1003 FORMAT(/,1x,'Aborting incflo.',1x,70('*'),2/)

      END SUBROUTINE FINL_ERR_MSG



!``````````````````````````````````````````````````````````````````````!
!                                                                      !
!......................................................................!
      SUBROUTINE FLUSH_ERR_MSG(DEBUG, HEADER, FOOTER, ABORT, LOG, &
         CALL_TREE)

! Dummy Arguments:
!---------------------------------------------------------------------//
! Debug flag.
      logical, INTENT(IN), OPTIONAL :: DEBUG
! Flag to suppress the message header.
      logical, INTENT(IN), OPTIONAL :: HEADER
! Flag to suppress the message footer.
      logical, INTENT(IN), OPTIONAL :: FOOTER
! Flag to abort execution
      logical, INTENT(IN), OPTIONAL :: ABORT
! Flag to force (or override) writing data to the log file.
      logical, INTENT(IN), OPTIONAL :: LOG
! Provide the call tree in error message.
      logical, INTENT(IN), OPTIONAL :: CALL_TREE

! Local Variables:
!---------------------------------------------------------------------//
! Single line.
      CHARACTER(LEN=LINE_LENGTH) :: LINE
! Line length with trailing space removed.
      integer :: LENGTH
! Index of last line in the message.
      integer :: LAST_LINE
! Line Counter
      integer :: LC
! Local debug flag.
      logical :: D_FLAG
! Local flag to suppress writing the header.
      logical :: H_FLAG
! Local flag to suppress writing the footer.
      logical :: F_FLAG
! Local abort flag.
      logical :: A_FLAG
! Local call tree flag.
      logical :: CT_FLAG

! The current calling routine.
      CHARACTER(LEN=128) :: CALLER

! Set the abort flag. Continue running by default.
      IF(PRESENT(ABORT))THEN
         A_FLAG = ABORT
      ELSE
         A_FLAG = .FALSE.
      ENDIF

! Set the local debug flag. Suppress debugging messages by default.
      IF(PRESENT(DEBUG)) THEN
         D_FLAG = DEBUG
      ELSE
         D_FLAG = .FALSE.
      ENDIF

! Set the header flag. Write the header by default.
      IF(PRESENT(HEADER)) THEN
         H_FLAG = HEADER
      ELSE
         H_FLAG = .TRUE.
      ENDIF

! Set the footer flag. Write the footer by default.
      IF(PRESENT(FOOTER))THEN
         F_FLAG = FOOTER
      ELSE
         F_FLAG = .TRUE.
      ENDIF


! Set the call tree flag. Suppress the call tree by default.
      IF(PRESENT(CALL_TREE)) THEN
         CT_FLAG = CALL_TREE
      ELSE
         CT_FLAG = .FALSE.
      ENDIF

! Write out header infomration.
      IF(H_FLAG) THEN
! Set the current caller.
         CALLER = CALLERS(CALL_DEPTH)
         IF(D_FLAG) THEN
            WRITE(*,2000) trim(CALLER)
         ELSE
            WRITE(*,1000) trim(CALLER)
         ENDIF
      ENDIF

! Find the end of the message.
      LAST_LINE = 0
      DO LC = 1, LINE_COUNT
         LINE = ERR_MSG(LC)
         LENGTH = len_trim(LINE)
         IF(0 < LENGTH .AND. LENGTH < 256 ) LAST_LINE = LC
      ENDDO

! Write the message body.
      IF(D_FLAG)THEN
         DO LC = 1, LINE_COUNT
            LINE = ERR_MSG(LC)
            LENGTH = len_trim(LINE)
            IF(LENGTH == 0) THEN
               WRITE(*,2001) LC, LENGTH, "EMPTY."
            ELSEIF(LENGTH >=  LINE_LENGTH)THEN
               WRITE(*,2001) LC, LENGTH, "OVERFLOW."
            ELSE
               WRITE(*,2001) LC, LENGTH, trim(LINE)
            ENDIF
         ENDDO
      ELSE
         DO LC = 1, LAST_LINE
            LINE = ERR_MSG(LC)
            LENGTH = len_trim(LINE)
            IF(0 < LENGTH .AND. LENGTH < 256 ) THEN
               WRITE(*,1001) trim(LINE)
            ELSE
               WRITE(*,"('  ')")
            ENDIF
         ENDDO
         IF(LAST_LINE == 0) THEN
            WRITE(*,"('  ')")
         ENDIF
      ENDIF

! Print footer.
      IF(F_FLAG) THEN
         IF(D_FLAG) THEN
            WRITE(*, 2002)
         ELSE
            WRITE(*, 1002)
         ENDIF
      ENDIF

! Clear the message array.
      ERR_MSG=''


! Abort the run if specified.
      IF(A_FLAG) stop 20007

      RETURN

 1000 FORMAT(2/,1x,70('*'),/' From: ',A)
 1001 FORMAT(1x,A)
 1002 FORMAT(1x,70('*'))

 2000 FORMAT(2/,'--- HEADER ---> ',70('*'),/'--- HEADER ---> From: ',A)
 2001 FORMAT('LC ',I2.2,': LEN: ',I3.3,1x,A)
 2002 FORMAT('--- FOOTER --->',1x,70('*'))

      END SUBROUTINE FLUSH_ERR_MSG


!``````````````````````````````````````````````````````````````````````!
!                                                                      !
!......................................................................!
      SUBROUTINE SHOW_CALL_TREE(HEADER, FOOTER)


! Dummy Arguments:
!---------------------------------------------------------------------//
! Flag to suppress the message header.
      logical, INTENT(IN), OPTIONAL :: HEADER
! Flag to suppress the message footer.
      logical, INTENT(IN), OPTIONAL :: FOOTER

! Local Variables:
!---------------------------------------------------------------------//
! Local flag to suppress writing the header.
      logical :: H_FLAG
! Local flag to suppress writing the footer.
      logical :: F_FLAG
! Generic loop counters.
      integer ::  LC, SL

! Set the header flag. Write the header by default.
      H_FLAG = merge(HEADER, .TRUE., PRESENT(HEADER))
! Set the footer flag. Write the footer by default.
      F_FLAG = merge(FOOTER, .TRUE., PRESENT(FOOTER))

! Header
      IF(H_FLAG) THEN
         WRITE(*,1000)
      ENDIF

! Call Tree
      DO LC=1,MAX_CALL_DEPTH
         DO SL=1,LC
            WRITE(*,1001,ADVANCE='NO')
         ENDDO
         WRITE(*,1002,ADVANCE='YES') CALLERS(LC)
      ENDDO

! Footer.
      IF(F_FLAG) THEN
         WRITE(*,1003)
      ENDIF

      RETURN

 1000 FORMAT(2/,1x,70('*'),' CALL TREE INFORMATION')
 1001 FORMAT(' ')
 1002 FORMAT('> ',A)
 1003 FORMAT(/1x,70('*'))

      END SUBROUTINE SHOW_CALL_TREE



!``````````````````````````````````````````````````````````````````````!
!                                                                      !
!......................................................................!
      CHARACTER(len=32) FUNCTION iVar(VAR, i1, i2, i3)

      CHARACTER(len=*), intent(in) :: VAR

      integer,  intent(in) :: i1
      integer, OPTIONAL, intent(in) :: i2
      integer, OPTIONAL, intent(in) :: i3

      CHARACTER(len=16) :: iASc
      CHARACTER(len=64) :: tVAR

      iASc=''; WRITE(iASc,*)i1
      tVar=''; WRITE(tVar,"(A,'(',A)") &
         trim(adjustl(VAR)), trim(adjustl(iASc))


      IF(PRESENT(i2))THEN
         iASc=''; WRITE(iASc,*)i2
         WRITE(tVar,"(A,',',A)") trim(tVar), trim(adjustl(iASc))
      ENDIF

      IF(PRESENT(i3))THEN
         iASc=''; WRITE(iASc,*)i3
         WRITE(tVar,"(A,',',A)") trim(tVar), trim(adjustl(iASc))
      ENDIF

      WRITE(tVar,"(A,')')") trim(tVar)

      iVar = trim(adjustl(tVar))

      RETURN
      END FUNCTION iVar


!``````````````````````````````````````````````````````````````````````!
!                                                                      !
!......................................................................!
      CHARACTER(len=32) FUNCTION iVal_int(VAL)
      integer, intent(in) :: VAL

      CHARACTER(len=32) :: iASc

      WRITE(iASc,*) VAL
      iVal_int = trim(adjustl(iASc))

      END FUNCTION iVal_int


!``````````````````````````````````````````````````````````````````````!
!                                                                      !
!......................................................................!
      CHARACTER(len=32) FUNCTION iVal_dbl(VAL)

      real(rt), intent(in) :: VAL

      CHARACTER(len=32) :: dASc

      IF(abs(VAL) < 1.0d-2 .AND. abs(VAL) < 1.0d2) THEN
         WRITE(dASc,"(F18.4)") VAL
      ELSE
         WRITE(dASc,"(G18.4)") VAL
      ENDIF

      iVal_dbl = trim(adjustl(dASc))

      END FUNCTION iVal_dbl


!``````````````````````````````````````````````````````````````````````!
!                                                                      !
!......................................................................!
      CHARACTER(len=32) FUNCTION iVal_log(VAL)
      logical, intent(in) :: VAL

      IF(VAL) THEN
         iVal_log = ".TRUE."
      ELSE
         iVal_log = ".FALSE."
      ENDIF

      RETURN
      END FUNCTION iVal_log


      END MODULE ERROR_MANAGER
