!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #############################################################
      SUBROUTINE INIT_SEAFLUX_n (DTCO, DGU, UG, U, SM, &
                                 HPROGRAM,HINIT,                            &
                                  KI,KSV,KSW,                                &
                                  HSV,PCO2,PRHOA,                            &
                                  PZENITH,PAZIM,PSW_BANDS,PDIR_ALB,PSCA_ALB, &
                                  PEMIS,PTSRAD,PTSURF,                       &
                                  KYEAR, KMONTH,KDAY, PTIME,                 &
                                  HATMFILE,HATMFILETYPE,                     &
                                  HTEST                                      )  
!     #############################################################
!
!!****  *INIT_SEAFLUX_n* - routine to initialize SEAFLUX
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
!!      V. Masson   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/2003
!!      Modified    01/2006 : sea flux parameterization.
!!                  01/2008 : coupling with 1D ocean
!!      B. Decharme 08/2009 : specific treatment for sea/ice in the Earth System Model 
!!      B. Decharme 07/2011 : read pgd+prep 
!!      B. Decharme 04/2013 : new coupling variables
!!      S. Senesi   01/2014 : introduce sea-ice model 
!!      S. Belamari 03/2014 : add NZ0 (to choose PZ0SEA formulation)
!!      R. Séférian 01/2015 : introduce interactive ocean surface albedo
!!      A. Voldoire 09/2016 : Switch to tile the fluxes calculation over sea and seaice
!!      R. Séférian    11/16 : Implement carbon cycle coupling (Earth system model)
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
USE MODD_SURFEX_n, ONLY : SEAFLUX_MODEL_t
!
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_t
USE MODD_DIAG_SURF_ATM_n, ONLY : DIAG_SURF_ATM_t
USE MODD_SURF_ATM_GRID_n, ONLY : SURF_ATM_GRID_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
!
USE MODN_SFX_OASIS,      ONLY : LSEAICE_2FLX
USE MODD_SFX_OASIS,      ONLY : LCPL_SEA, LCPL_SEAICE, LCPL_SEACARB
!
USE MODD_READ_NAMELIST,  ONLY : LNAM_READ
USE MODD_CSTS,           ONLY : XTTS
USE MODD_SNOW_PAR,       ONLY : XZ0HSN
USE MODD_SURF_PAR,       ONLY : XUNDEF, NUNDEF
USE MODD_CHS_AEROSOL,    ONLY: LVARSIGI, LVARSIGJ
USE MODD_DST_SURF,       ONLY: LVARSIG_DST, NDSTMDE, NDST_MDEBEG, LRGFIX_DST
USE MODD_SLT_SURF,       ONLY: LVARSIG_SLT, NSLTMDE, NSLT_MDEBEG, LRGFIX_SLT
!
USE MODI_INIT_IO_SURF_n
USE MODI_DEFAULT_CH_DEP
!
USE MODI_DEFAULT_SEAFLUX
USE MODI_DEFAULT_DIAG_SEAFLUX
USE MODI_READ_DEFAULT_SEAFLUX_n
USE MODI_READ_SEAFLUX_CONF_n
USE MODI_READ_SEAFLUX_n
!
USE MODI_READ_OCEAN_n
!
USE MODI_DEFAULT_SEAICE
USE MODI_READ_SEAICE_n
!
USE MODI_READ_PGD_SEAFLUX_n
USE MODI_DIAG_SEAFLUX_INIT_n
USE MODI_END_IO_SURF_n
USE MODI_GET_LUOUT
USE MODI_READ_SURF
USE MODI_READ_SEAFLUX_DATE
USE MODI_READ_NAM_PREP_SEAFLUX_n
USE MODI_INIT_CHEMICAL_n
USE MODI_PREP_CTRL_SEAFLUX
USE MODI_UPDATE_RAD_SEA
USE MODI_READ_SEAFLUX_SBL_n
USE MODI_ABOR1_SFX
!
USE MODI_SET_SURFEX_FILEIN
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
!
TYPE(DATA_COVER_t), INTENT(INOUT) :: DTCO
TYPE(DIAG_SURF_ATM_t), INTENT(INOUT) :: DGU
TYPE(SURF_ATM_GRID_t), INTENT(INOUT) :: UG
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
TYPE(SEAFLUX_MODEL_t), INTENT(INOUT) :: SM
!
 CHARACTER(LEN=6),                 INTENT(IN)  :: HPROGRAM  ! program calling surf. schemes
 CHARACTER(LEN=3),                 INTENT(IN)  :: HINIT     ! choice of fields to initialize
INTEGER,                          INTENT(IN)  :: KI        ! number of points
INTEGER,                          INTENT(IN)  :: KSV       ! number of scalars
INTEGER,                          INTENT(IN)  :: KSW       ! number of short-wave spectral bands
 CHARACTER(LEN=6), DIMENSION(KSV), INTENT(IN)  :: HSV       ! name of all scalar variables
REAL,             DIMENSION(KI),  INTENT(IN)  :: PCO2      ! CO2 concentration (kg/m3)
REAL,             DIMENSION(KI),  INTENT(IN)  :: PRHOA     ! air density
REAL,             DIMENSION(KI),  INTENT(IN)  :: PZENITH   ! solar zenithal angle
REAL,             DIMENSION(KI),  INTENT(IN)  :: PAZIM     ! solar azimuthal angle (rad from N, clock)
REAL,             DIMENSION(KSW), INTENT(IN)  :: PSW_BANDS ! middle wavelength of each band
REAL,             DIMENSION(KI,KSW),INTENT(OUT) :: PDIR_ALB  ! direct albedo for each band
REAL,             DIMENSION(KI,KSW),INTENT(OUT) :: PSCA_ALB  ! diffuse albedo for each band
REAL,             DIMENSION(KI),  INTENT(OUT) :: PEMIS     ! emissivity
REAL,             DIMENSION(KI),  INTENT(OUT) :: PTSRAD    ! radiative temperature
REAL,             DIMENSION(KI),  INTENT(OUT) :: PTSURF    ! surface effective temperature         (K)
INTEGER,                          INTENT(IN)  :: KYEAR     ! current year (UTC)
INTEGER,                          INTENT(IN)  :: KMONTH    ! current month (UTC)
INTEGER,                          INTENT(IN)  :: KDAY      ! current day (UTC)
REAL,                             INTENT(IN)  :: PTIME     ! current time since
                                                           !  midnight (UTC, s)
!
 CHARACTER(LEN=28),                INTENT(IN)  :: HATMFILE    ! atmospheric file name
 CHARACTER(LEN=6),                 INTENT(IN)  :: HATMFILETYPE! atmospheric file type
 CHARACTER(LEN=2),                 INTENT(IN)  :: HTEST       ! must be equal to 'OK'
!
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER           :: ILU    ! sizes of SEAFLUX arrays
INTEGER           :: ILUOUT ! unit of output listing file
INTEGER           :: IRESP  ! return code
LOGICAL           :: GSIC
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!         Initialisation for IO
!
IF (LHOOK) CALL DR_HOOK('INIT_SEAFLUX_N',0,ZHOOK_HANDLE)
 CALL GET_LUOUT(HPROGRAM,ILUOUT)
!
IF (HTEST/='OK') THEN
  CALL ABOR1_SFX('INIT_SEAFLUXN: FATAL ERROR DURING ARGUMENT TRANSFER')
END IF
!
!
!         Others litlle things
!
PDIR_ALB = XUNDEF
PSCA_ALB = XUNDEF
PEMIS    = XUNDEF
PTSRAD   = XUNDEF
PTSURF   = XUNDEF
!
SM%O%LMERCATOR = .FALSE.
SM%O%LCURRENT  = .FALSE.
!
IF (LNAM_READ) THEN
 !
 !*       0.     Defaults
 !               --------
 !
 !        0.1. Hard defaults
 !      
 
 CALL DEFAULT_SEAFLUX(SM%S%XTSTEP,SM%S%XOUT_TSTEP,SM%S%CSEA_ALB,SM%S%CSEA_FLUX,         &
                      SM%S%CSEA_SFCO2, SM%S%LPWG,                                       &
                      SM%S%LPRECIP,SM%S%LPWEBB,SM%S%NZ0,SM%S%NGRVWAVES,SM%O%LPROGSST,   &
                      SM%O%NTIME_COUPLING,SM%O%XOCEAN_TSTEP,SM%S%XICHCE,SM%S%CINTERPOL_SST,&
                      SM%S%CINTERPOL_SSS                            )
 CALL DEFAULT_SEAICE(HPROGRAM,                                   &
                     SM%S%CINTERPOL_SIC,SM%S%CINTERPOL_SIT, SM%S%XFREEZING_SST, &
                     SM%S%XSEAICE_TSTEP, SM%S%XSIC_EFOLDING_TIME,          &
                     SM%S%XSIT_EFOLDING_TIME, SM%S%XCD_ICE_CST, SM%S%XSI_FLX_DRV)     
 !                     
 CALL DEFAULT_CH_DEP(SM%CHS%CCH_DRY_DEP) 
 !            
 CALL DEFAULT_DIAG_SEAFLUX(SM%DGS%N2M,SM%DGS%LSURF_BUDGET,SM%DGS%L2M_MIN_ZS,&
                        SM%DGS%LRAD_BUDGET,SM%DGS%LCOEF,SM%DGS%LSURF_VARS,&
                           SM%DGO%LDIAG_OCEAN,SM%DGSI%LDIAG_SEAICE,SM%DGS%LSURF_BUDGETC,&
                          SM%DGS%LRESET_BUDGETC,SM%DGS%XDIAG_TSTEP )  

ENDIF
!
!
!        0.2. Defaults from file header
!    
 CALL READ_DEFAULT_SEAFLUX_n(SM%CHS, SM%DGO, SM%DGS, SM%DGSI, SM%O, SM%S, &
                             HPROGRAM)
!
!*       1.1    Reading of configuration:
!               -------------------------
!
 CALL READ_SEAFLUX_CONF_n(SM%CHS, SM%DGO, SM%DGS, SM%DGSI, SM%O, SM%S, &
                          HPROGRAM)
!
SM%S%LINTERPOL_SST=.FALSE.
SM%S%LINTERPOL_SSS=.FALSE.
SM%S%LINTERPOL_SIC=.FALSE.
SM%S%LINTERPOL_SIT=.FALSE.
IF(LCPL_SEA)THEN 
! No STT / SSS interpolation in Earth System Model
  SM%S%CINTERPOL_SST='NONE  '
  SM%S%CINTERPOL_SSS='NONE  '
  SM%S%CINTERPOL_SIC='NONE  '
  SM%S%CINTERPOL_SIT='NONE  '
ELSE
   IF(TRIM(SM%S%CINTERPOL_SST)/='NONE'.AND.TRIM(SM%S%CINTERPOL_SST)/='READAY')THEN
      SM%S%LINTERPOL_SST=.TRUE.
   ENDIF
   IF(TRIM(SM%S%CINTERPOL_SSS)/='NONE'.AND.TRIM(SM%S%CINTERPOL_SSS)/='READAY')THEN
      SM%S%LINTERPOL_SSS=.TRUE.
   ENDIF
   IF(TRIM(SM%S%CINTERPOL_SIC)/='NONE'.AND.TRIM(SM%S%CINTERPOL_SIC)/='READAY')THEN
      SM%S%LINTERPOL_SIC=.TRUE.
   ENDIF
   IF(TRIM(SM%S%CINTERPOL_SIT)/='NONE'.AND.TRIM(SM%S%CINTERPOL_SIT)/='READAY')THEN
      SM%S%LINTERPOL_SIT=.TRUE.
   ENDIF
ENDIF
!
!*       1.     Cover fields and grid:
!               ---------------------
!* date
!
SELECT CASE (HINIT)
!
  CASE ('PGD')
!
    SM%S%TTIME%TDATE%YEAR = NUNDEF
    SM%S%TTIME%TDATE%MONTH= NUNDEF
    SM%S%TTIME%TDATE%DAY  = NUNDEF
    SM%S%TTIME%TIME       = XUNDEF
!
  CASE ('PRE')
!
    CALL PREP_CTRL_SEAFLUX(SM%DGS%N2M,SM%DGS%LSURF_BUDGET,SM%DGS%L2M_MIN_ZS,&
                                SM%DGS%LRAD_BUDGET,SM%DGS%LCOEF,SM%DGS%LSURF_VARS,&
                             SM%DGO%LDIAG_OCEAN,SM%DGSI%LDIAG_SEAICE,ILUOUT,SM%DGS%LSURF_BUDGETC ) 
    IF (LNAM_READ) CALL READ_NAM_PREP_SEAFLUX_n(HPROGRAM)      
    CALL READ_SEAFLUX_DATE(SM%O, &
                           HPROGRAM,HINIT,ILUOUT,HATMFILE,HATMFILETYPE,KYEAR,KMONTH,KDAY,PTIME,SM%S%TTIME)
!
  CASE DEFAULT
!
 CALL INIT_IO_SURF_n(DTCO, DGU, U, &
                        HPROGRAM,'SEA   ','SEAFLX','READ ')
    CALL READ_SURF(&
                   HPROGRAM,'DTCUR',SM%S%TTIME,IRESP)
    CALL END_IO_SURF_n(HPROGRAM)
!
END SELECT
!
!-----------------------------------------------------------------------------------------------------
! READ PGD FILE
!-----------------------------------------------------------------------------------------------------
!
!         Initialisation for IO
!
 CALL SET_SURFEX_FILEIN(HPROGRAM,'PGD ') ! change input file name to pgd name
 CALL INIT_IO_SURF_n(DTCO, DGU, U, &
                        HPROGRAM,'SEA   ','SEAFLX','READ ')
!
!         Reading of the fields
!
 CALL READ_PGD_SEAFLUX_n(DTCO, SM%DTS, SM%SG, SM%S, U, &
                         HPROGRAM)
!
 CALL END_IO_SURF_n(HPROGRAM)
 CALL SET_SURFEX_FILEIN(HPROGRAM,'PREP') ! restore input file name
!-------------------------------------------------------------------------------
!
!* if only physiographic fields are to be initialized, stop here.
!
IF (HINIT/='ALL' .AND. HINIT/='SOD') THEN
  IF (LHOOK) CALL DR_HOOK('INIT_SEAFLUX_N',1,ZHOOK_HANDLE)
  RETURN
END IF
!
!-------------------------------------------------------------------------------
!
!         Initialisation for IO
!
 CALL INIT_IO_SURF_n(DTCO, DGU, U, &
                        HPROGRAM,'SEA   ','SEAFLX','READ ')
!
!*       2.     Prognostic fields:
!               ----------------
!
IF(SM%S%LINTERPOL_SST.OR.SM%S%LINTERPOL_SSS.OR.SM%S%LINTERPOL_SIC.OR.SM%S%LINTERPOL_SIT)THEN
!  Initialize current Month for SST interpolation
   SM%S%TZTIME%TDATE%YEAR  = SM%S%TTIME%TDATE%YEAR
   SM%S%TZTIME%TDATE%MONTH = SM%S%TTIME%TDATE%MONTH
   SM%S%TZTIME%TDATE%DAY   = SM%S%TTIME%TDATE%DAY
   SM%S%TZTIME%TIME        = SM%S%TTIME%TIME        
ENDIF
!
 CALL READ_SEAFLUX_n(DTCO, SM%SG, SM%S, U, &
                     HPROGRAM,ILUOUT)
!
IF (HINIT/='ALL') THEN
  CALL END_IO_SURF_n(HPROGRAM)
  IF (LHOOK) CALL DR_HOOK('INIT_SEAFLUX_N',1,ZHOOK_HANDLE)
  RETURN
END IF
!-------------------------------------------------------------------------------
!
!*       2.1    Ocean fields:
!               -------------
!
 CALL READ_OCEAN_n(DTCO, SM%O, SM%OR, U, &
                   HPROGRAM)
!
!-------------------------------------------------------------------------------
!
ILU = SIZE(SM%S%XCOVER,1)
!
ALLOCATE(SM%S%XSST_INI    (ILU))
SM%S%XSST_INI(:) = SM%S%XSST(:)
!
ALLOCATE(SM%S%XZ0H(ILU))
WHERE (SM%S%XSST(:)>=XTTS)
  SM%S%XZ0H(:) = SM%S%XZ0(:)
ELSEWHERE
  SM%S%XZ0H(:) = XZ0HSN
ENDWHERE
!
!* ocean surface albedo (direct and diffuse fraction)
!
ALLOCATE(SM%S%XDIR_ALB (ILU))
ALLOCATE(SM%S%XSCA_ALB (ILU))
!
SM%S%XDIR_ALB(:) = XUNDEF
SM%S%XSCA_ALB(:) = XUNDEF
!
!* sea-ice cover (default = 0)
!
ALLOCATE(SM%S%XSIC(ILU))
SM%S%XSIC(:)=0.0
!
!-------------------------------------------------------------------------------
!
!*       3.     Specific fields when using earth system model or sea-ice scheme
!               (Sea current and Sea-ice temperature)
!               -----------------------------------------------------------------
!
IF(LCPL_SEA)THEN       
! 
  ALLOCATE(SM%S%XUMER   (ILU))
  ALLOCATE(SM%S%XVMER   (ILU))
!
  SM%S%XUMER   (:)=XUNDEF
  SM%S%XVMER   (:)=XUNDEF
!
ELSE
! 
  ALLOCATE(SM%S%XUMER   (0))
  ALLOCATE(SM%S%XVMER   (0))
!
ENDIF
!
IF(LCPL_SEAICE.OR.SM%S%LHANDLE_SIC)THEN       
  ALLOCATE(SM%S%XTICE   (ILU))
  ALLOCATE(SM%S%XICE_ALB(ILU))
  SM%S%XTICE   (:)=XUNDEF
  SM%S%XICE_ALB(:)=XUNDEF
ELSE
  ALLOCATE(SM%S%XTICE   (0))
  ALLOCATE(SM%S%XICE_ALB(0))
ENDIF
!
IF(LCPL_SEACARB .AND. LCPL_SEA)THEN       
  ALLOCATE(SM%S%XSFCO2(ILU))
  SM%S%XSFCO2    (:)=XUNDEF
ELSE
  ALLOCATE(SM%S%XSFCO2(0))
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       4.     Seaice prognostic variables and forcings :
!
IF (SM%S%LHANDLE_SIC) THEN
  CALL READ_SEAICE_n(&
                   SM%SG, SM%S, &
                   HPROGRAM,ILU,ILUOUT)
ENDIF
!
! Activation of double flux param for ESM
!
IF (LCPL_SEAICE.AND.LSEAICE_2FLX) SM%S%LHANDLE_SIC=.TRUE.
!
!-------------------------------------------------------------------------------
!
!*       5.     Albedo, emissivity and temperature fields on the mix (open sea + sea ice)
!               -----------------------------------------------------------------
!
ALLOCATE(SM%S%XEMIS(ILU))
SM%S%XEMIS(:) = 0.0
!
IF (.NOT.LCPL_SEA) THEN
!
  CALL UPDATE_RAD_SEA(SM%S%CSEA_ALB,SM%S%XSST,PZENITH,XTTS,SM%S%XEMIS,   &
                      SM%S%XDIR_ALB_SEA,SM%S%XSCA_ALB_SEA,SM%S%XDIR_ALB, &
                      SM%S%XSCA_ALB,PDIR_ALB,PSCA_ALB,PEMIS,PTSRAD,      &
                      SM%S%LHANDLE_SIC,SM%S%XTICE,SM%S%XSIC,SM%S%XICE_ALB)  
!
  IF (SM%S%LHANDLE_SIC) THEN
     PTSURF(:) = SM%S%XSST(:) * ( 1.0 - SM%S%XSIC(:)) + SM%S%XTICE(:) * SM%S%XSIC(:)
  ELSE
     PTSURF(:) = SM%S%XSST(:)
  ENDIF
!
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       6.     SBL air fields:
!               --------------
!
 CALL READ_SEAFLUX_SBL_n(DTCO, SM%S, SM%SSB, U, &
                         HPROGRAM)
!
!-------------------------------------------------------------------------------
!
!*       7.     Chemistry /dust
!               ---------
!
 CALL INIT_CHEMICAL_n(ILUOUT, KSV, HSV, SM%CHS%SVS,     &
                     SM%CHS%CCH_NAMES, SM%CHS%CAER_NAMES,     &
                     HDSTNAMES=SM%CHS%CDSTNAMES, HSLTNAMES=SM%CHS%CSLTNAMES        )
!
!* deposition scheme
!
IF (SM%CHS%SVS%NBEQ>0 .AND. SM%CHS%CCH_DRY_DEP=='WES89') THEN
  ALLOCATE(SM%CHS%XDEP(ILU,SM%CHS%SVS%NBEQ))
ELSE
  ALLOCATE(SM%CHS%XDEP(0,0))
END IF
!
!-------------------------------------------------------------------------------
!
!*       8.     diagnostics initialization
!               --------------------------
!
GSIC=(SM%S%LHANDLE_SIC.AND.(SM%S%CSEAICE_SCHEME /= 'NONE  '))
!
IF(HINIT=='ALL'.AND.(DGU%LDIAG_CMIP.OR.GSIC))THEN
!
  SM%DGS%N2M           = 2
  SM%DGS%LSURF_BUDGET  = .TRUE.
!
  IF(DGU%LDIAG_CMIP)THEN
    SM%DGS%LRAD_BUDGET   = .TRUE.
    SM%DGS%LSURF_BUDGETC = .FALSE.
  ENDIF
!
  IF(LCPL_SEAICE)THEN
    SM%DGSI%LDIAG_SEAICE = .FALSE.
  ELSEIF(SM%S%LHANDLE_SIC)THEN
    SM%DGSI%LDIAG_SEAICE = .TRUE.
  ENDIF
!
ENDIF
!
IF(HINIT=='ALL'.AND.(.NOT.SM%S%LHANDLE_SIC).AND.(.NOT.LCPL_SEAICE))THEN
  SM%DGSI%LDIAG_SEAICE=.FALSE.
ENDIF
!
IF(LCPL_SEA.AND.SM%DGS%N2M<1)THEN
   CALL ABOR1_SFX('INIT_SEAFLUX_n: N2M must be set >0 in case of LCPL_SEA')
ENDIF
!
 CALL DIAG_SEAFLUX_INIT_n(&
                         SM%DGO, SM%DGS, SM%DGSI, DGU, SM%S, &
                         HPROGRAM,ILU,KSW)
!
!-------------------------------------------------------------------------------
!
!         End of IO
!
 CALL END_IO_SURF_n(HPROGRAM)
IF (LHOOK) CALL DR_HOOK('INIT_SEAFLUX_N',1,ZHOOK_HANDLE)
!
!
END SUBROUTINE INIT_SEAFLUX_n
