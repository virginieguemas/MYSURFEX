!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
      SUBROUTINE WRITESURF_SEAFLUX_n (DGU, U, &
                                       O, OR, S, &
                                      HPROGRAM)
!     ########################################
!
!!****  *WRITE_SEAFLUX_n* - writes SEAFLUX fields
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!	V. Masson   *Meteo France*	
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/2003 
!!      Modified    01/2014 : S. Senesi : handle seaice scheme
!!      S. Belamari 03/2014   Include sea surface salinity XSSS
!!      R. Séférian 01/2015 : introduce interactive ocean surface albedo
!!      S. Senesi   08/2015 : fix units in some HCOMMENTs
!!      A. Voldoire 09/2016 : Switch to tile the fluxes calculation over sea and seaice
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
!
!
!
!
USE MODD_DIAG_SURF_ATM_n, ONLY : DIAG_SURF_ATM_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
!
USE MODD_OCEAN_n, ONLY : OCEAN_t
USE MODD_OCEAN_REL_n, ONLY : OCEAN_REL_t
USE MODD_SEAFLUX_n, ONLY : SEAFLUX_t
!
USE MODN_SFX_OASIS,      ONLY : LSEAICE_2FLX
!
USE MODI_WRITE_SURF
USE MODI_WRITESURF_OCEAN_n
USE MODI_WRITESURF_SEAICE_N
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
!
!
TYPE(DIAG_SURF_ATM_t), INTENT(INOUT) :: DGU
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
!
TYPE(OCEAN_t), INTENT(INOUT) :: O
TYPE(OCEAN_REL_t), INTENT(INOUT) :: OR
TYPE(SEAFLUX_t), INTENT(INOUT) :: S
!
 CHARACTER(LEN=6),  INTENT(IN)  :: HPROGRAM ! program calling

!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER           :: JMTH, INMTH
CHARACTER(LEN=2 ) :: YMTH
!
INTEGER           :: IRESP          ! IRESP  : return-code if a problem appears
CHARACTER(LEN=12) :: YRECFM         ! Name of the article to be read
CHARACTER(LEN=100):: YCOMMENT       ! Comment string
!
LOGICAL           :: GHANDLE_SIC
!
REAL(KIND=JPRB)   :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!
IF (LHOOK) CALL DR_HOOK('WRITESURF_SEAFLUX_N',0,ZHOOK_HANDLE)
!
CALL WRITESURF_OCEAN_n(DGU, U, &
                       O, OR, &
                       HPROGRAM)
!
!*       2.     Sea-ice prognostic fields:
!               --------------------------
!
!* flag to tell if Sea Ice model is used
!
GHANDLE_SIC = S%LHANDLE_SIC
!
!coupled configuration case : HANDLE_SIC must be false in the restar file
IF (LSEAICE_2FLX) GHANDLE_SIC = .FALSE.
!
YCOMMENT='flag to handle sea ice cover'
CALL WRITE_SURF(DGU, U, &
                HPROGRAM,'HANDLE_SIC',GHANDLE_SIC,IRESP,YCOMMENT)
!
IF (GHANDLE_SIC) CALL WRITESURF_SEAICE_n(DGU, U, S, &
                                         HPROGRAM)
!
!
!*       3.     Prognostic fields:
!               -----------------
!
!* water temperature
!
IF(S%LINTERPOL_SST)THEN
!
  INMTH=SIZE(S%XSST_MTH,2)
!
  DO JMTH=1,INMTH
     WRITE(YMTH,'(I2)') (JMTH-1)
     YRECFM='SST_MTH'//ADJUSTL(YMTH(:LEN_TRIM(YMTH)))
     YCOMMENT='SST at month t'//ADJUSTL(YMTH(:LEN_TRIM(YMTH)))//' (K)'
     CALL WRITE_SURF(DGU, U, &
                HPROGRAM,YRECFM,S%XSST_MTH(:,JMTH),IRESP,HCOMMENT=YCOMMENT)
  ENDDO
!
ENDIF
!
YRECFM='SST'
YCOMMENT='SST (K)'
CALL WRITE_SURF(DGU, U, &
                HPROGRAM,YRECFM,S%XSST(:),IRESP,HCOMMENT=YCOMMENT)  
!
!-------------------------------------------------------------------------------
!
!*       4.     Semi-prognostic fields:
!               ----------------------
!
!* roughness length
!
YRECFM='Z0SEA'
YCOMMENT='Z0SEA (m)'
CALL WRITE_SURF(DGU, U, &
                HPROGRAM,YRECFM,S%XZ0(:),IRESP,HCOMMENT=YCOMMENT)
!
!
!* sea surface salinity
!
IF(S%LINTERPOL_SSS)THEN
   !
   INMTH=SIZE(S%XSSS_MTH,2)
   !
   DO JMTH=1,INMTH
      WRITE(YMTH,'(I2)') (JMTH-1)
      YRECFM='SSS_MTH'//ADJUSTL(YMTH(:LEN_TRIM(YMTH)))
      YCOMMENT='Sea Surface Salinity at month t'//ADJUSTL(YMTH(:LEN_TRIM(YMTH)))//' (psu)'
      CALL WRITE_SURF(DGU, U, &
                HPROGRAM,YRECFM,S%XSSS_MTH(:,JMTH),IRESP,HCOMMENT=YCOMMENT)
   ENDDO
!
ENDIF
!
YRECFM='SSS'
YCOMMENT='Sea Surface Salinity (psu)'
CALL WRITE_SURF(DGU, U, &
                HPROGRAM,YRECFM,S%XSSS(:),IRESP,HCOMMENT=YCOMMENT)  
!
!
!* ocean surface albedo (direct and diffuse fraction)
!
IF(S%CSEA_ALB=='RS14')THEN
!
  YRECFM='OSA_DIR'
  YCOMMENT='direct ocean surface albedo (-)'
  CALL WRITE_SURF(DGU, U, &
                HPROGRAM,YRECFM,S%XDIR_ALB_SEA(:),IRESP,HCOMMENT=YCOMMENT)
!
  YRECFM='OSA_SCA'
  YCOMMENT='diffuse ocean surface albedo (-)'
  CALL WRITE_SURF(DGU, U, &
                HPROGRAM,YRECFM,S%XSCA_ALB_SEA(:),IRESP,HCOMMENT=YCOMMENT)
!
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       5.  Time
!            ----
!
YRECFM='DTCUR'
YCOMMENT='s'
 CALL WRITE_SURF(DGU, U, &
                HPROGRAM,YRECFM,S%TTIME,IRESP,HCOMMENT=YCOMMENT)
IF (LHOOK) CALL DR_HOOK('WRITESURF_SEAFLUX_N',1,ZHOOK_HANDLE)
!
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE WRITESURF_SEAFLUX_n
