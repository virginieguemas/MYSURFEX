!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
      SUBROUTINE INIT_WRITE_TXT (DGU, &
                                 HREC,OWFL)
!     ######################
!
!!****  *INIT_WRITE_TXT_n* Initialize array name to be written and associated
!!                         unit number
!!
!!    PURPOSE
!!    -------
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!      A. LEMONSU     *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!
!
!*       0.   DECLARATIONS
!             ------------
!
!
!
USE MODD_DIAG_SURF_ATM_n, ONLY : DIAG_SURF_ATM_t
!
USE MODD_IO_SURF_TXT,ONLY:NMASK, NFULL, CMASK
USE MODD_WRITE_TXT,  ONLY:NUNIT0, NVAR, CVAR, CVARN, JPVAR, NIND
!
USE MODD_WRITE_SURF_ATM, ONLY : LFIRST_WRITE, NCPT_WRITE
!
USE MODI_ABOR1_SFX
USE MODI_TEST_RECORD_LEN
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!
TYPE(DIAG_SURF_ATM_t), INTENT(INOUT) :: DGU
!
 CHARACTER(LEN=12),   INTENT(IN)     :: HREC    
LOGICAL,             INTENT(INOUT)  :: OWFL
INTEGER                             :: IP, IVAR, IFIELD, JFIELD
LOGICAL                             :: LMATCH
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('INIT_WRITE_TXT',0,ZHOOK_HANDLE)
!
IVAR=NUNIT0
DO IP=1, JPVAR
  IF (HREC==CVAR(IP)) THEN
    IVAR=NVAR(IP)
    EXIT
  ELSEIF(HREC==CVARN(IP)) THEN
    IVAR=-1
    EXIT
  ENDIF
ENDDO
!
!
IF (IVAR.LT.0) THEN
!
  OWFL=.FALSE.
!
ELSEIF (IVAR.NE.NUNIT0) THEN
!
  OWFL=.TRUE.
!
ELSE
!
  IF (CVAR(1).NE.'                ') IVAR=MAXVAL(NVAR(:))
!
!
  IF (.NOT.DGU%LSELECT) THEN
!
    IF ( (HREC(5:7)/='_OC'                          ) .AND.  & 
         (HREC(4:6)/='_OC'                          )   ) THEN  

      IVAR = IVAR+1
      IF (IVAR-NUNIT0>JPVAR) THEN
        CALL ABOR1_SFX('TOO MANY FIELDS TO BE WRITTEN IN THE "TEXTE" TYPE TIMESERIES')
      END IF
      CVAR(IVAR-NUNIT0) = HREC
      NVAR(IVAR-NUNIT0) = IVAR
      OPEN(UNIT=IVAR,FILE=TRIM(HREC)//'.TXT',FORM='FORMATTED')
      OWFL=.TRUE.
   
    ELSE
      IP = 1
      DO WHILE (CVARN(IP).NE.'                ') 
        IP=IP+1
      ENDDO
      CVARN(IP) = HREC
      OWFL=.FALSE.
    ENDIF
!
  ELSE
!        
    IFIELD=0
    DO JFIELD=1,SIZE(DGU%CSELECT)
      IF (DGU%CSELECT(JFIELD)== '            ') EXIT
      IFIELD=IFIELD+1
    ENDDO
  
    CALL TEST_RECORD_LEN("ASCII ",HREC,DGU%LSELECT,DGU%CSELECT,LMATCH)
    NCPT_WRITE=NCPT_WRITE-1

    IF (.NOT. LMATCH ) THEN

      IVAR = IVAR+1
      IF (IVAR-NUNIT0>JPVAR) THEN
        CALL ABOR1_SFX('TOO MANY FIELDS TO BE WRITTEN IN THE "TEXTE" TYPE TIMESERIES')
      END IF
      CVAR(IVAR-NUNIT0) = HREC
      NVAR(IVAR-NUNIT0) = IVAR
      OPEN(UNIT=IVAR,FILE=TRIM(HREC)//'.TXT',FORM='FORMATTED')
      OWFL=.TRUE.

    ELSE
      OWFL=.FALSE.
    ENDIF

  ENDIF
ENDIF

NIND=IVAR
IF (LHOOK) CALL DR_HOOK('INIT_WRITE_TXT',1,ZHOOK_HANDLE)
!
!------------------------------------------------------------------------------
!
END SUBROUTINE INIT_WRITE_TXT
