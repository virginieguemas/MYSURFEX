!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
    SUBROUTINE COARE30_FLUX (S, &
                             PZ0SEA,PTA,PEXNA,PRHOA,PSST,PEXNS,PQA,  &
            PVMOD,PZREF,PUREF,PPS,PQSAT,PSFTH,PSFTQ,PUSTAR,PCD,PCDN,PCH,PCE,PRI,&
            PRESA,PRAIN,PZ0HSEA)  
!     #######################################################################
!
!
!!****  *COARE25_FLUX*  
!!
!!    PURPOSE
!!    -------
!      Calculate the surface fluxes of heat, moisture, and momentum over
!      sea surface with COARE3.0 bulk algorithm.
!     
!!**  METHOD
!!    ------
!      transfer coefficients were obtained using a dataset which combined COARE
!      data with those from three other ETL field experiments, and reanalysis of
!      the HEXMAX data (DeCosmos et al. 1996). 
!      ITERMAX=3 
!      Take account of the surface gravity waves on the velocity roughness and 
!      hence the momentum transfer coefficient
!        NGRVWAVES=0 no gravity waves action (Charnock) !default value
!        NGRVWAVES=1 wave age parameterization of Oost et al. 2002
!        NGRVWAVES=2 model of Taylor and Yelland 2001
!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!!      
!!    REFERENCE
!!    ---------
!!      Fairall et al (2003), J. of Climate, vol. 16, 571-591
!!      Fairall et al (1996), JGR, 3747-3764
!!      Gosnell et al (1995), JGR, 437-442
!!      Fairall et al (1996), JGR, 1295-1308
!!      
!!    AUTHOR
!!    ------
!!     C. Lebeaupin  *Météo-France* (adapted from C. Fairall's code)
!!
!!    MODIFICATIONS
!!    -------------
!!      Original     1/06/2006
!!      B. Decharme    06/2009 limitation of Ri
!!      B. Decharme    09/2012 Bug in Ri (temperature in Celsius instead of Kelvin) 
!!                             and limitation of Ri in surface_ri.F90
!!      B. Decharme    06/2013 bug in z0 (output) computation 
!!      J.Escobar      06/2013  for REAL4/8 add EPSILON management
!!      C. Lebeaupin   03/2014 bug if PTA=PSST and PEXNA=PEXNS: set a minimum value
!!                             add abort if no convergence
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
!
USE MODD_SEAFLUX_n, ONLY : SEAFLUX_t
!
USE MODD_CSTS,       ONLY : XKARMAN, XG, XSTEFAN, XRD, XRV, XPI, &
                            XLVTT, XCL, XCPD, XCPV, XRHOLW, XTT, &
                            XP00
USE MODD_SURF_ATM,   ONLY : XVZ0CM
!
USE MODD_SURF_PAR,   ONLY : XUNDEF, XSURF_EPSILON
USE MODD_WATER_PAR
!
USE MODI_SURFACE_RI
USE MODI_WIND_THRESHOLD
USE MODE_COARE30_PSI
!
USE MODE_THERMOS
!
!
USE MODI_ABOR1_SFX
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
!
!
TYPE(SEAFLUX_t), INTENT(INOUT) :: S
!
REAL, DIMENSION(:), INTENT(IN)       :: PTA   ! air temperature at PZREF atm. level (K)
REAL, DIMENSION(:), INTENT(IN)       :: PQA   ! air humidity at PZREF atm. level (kg/kg)
REAL, DIMENSION(:), INTENT(IN)       :: PEXNA ! Exner function at PZREF atm. level
REAL, DIMENSION(:), INTENT(IN)       :: PRHOA ! air density at PZREF atm. level
REAL, DIMENSION(:), INTENT(IN)       :: PVMOD ! module of wind at PUREF atm. wind level (m/s)
REAL, DIMENSION(:), INTENT(IN)       :: PZREF ! atm. level for temp. and humidity (m)
REAL, DIMENSION(:), INTENT(IN)       :: PUREF ! atm. level for wind (m)
REAL, DIMENSION(:), INTENT(IN)       :: PSST  ! Sea Surface Temperature (K)
REAL, DIMENSION(:), INTENT(IN)       :: PEXNS ! Exner function at sea surface
REAL, DIMENSION(:), INTENT(IN)       :: PPS   ! air pressure at sea surface (Pa)
REAL, DIMENSION(:), INTENT(IN)       :: PRAIN ! precipitation rate (kg/s/m2)
!
REAL, DIMENSION(:), INTENT(INOUT)    :: PZ0SEA! aerodynamic roughness length over the ocean - what for ?
!                                             ! We already have a better estimate with ZO
!                                                                                 
!  surface fluxes : latent heat, sensible heat, friction fluxes
REAL, DIMENSION(:), INTENT(OUT)      :: PSFTH ! upward sensible heat flux (W/m2)
REAL, DIMENSION(:), INTENT(OUT)      :: PSFTQ ! upward water flux (kg/m2/s)
REAL, DIMENSION(:), INTENT(OUT)      :: PUSTAR! friction velocity (m/s)
!
! diagnostics
REAL, DIMENSION(:), INTENT(OUT)      :: PQSAT ! sea surface specific humidity
REAL, DIMENSION(:), INTENT(OUT)      :: PCD   ! momentum transfer coefficient at PUREF
REAL, DIMENSION(:), INTENT(OUT)      :: PCDN  ! neutral momentum transfer coefficient at 10m
REAL, DIMENSION(:), INTENT(OUT)      :: PCH   ! heat transfer coefficient at PZREF
REAL, DIMENSION(:), INTENT(OUT)      :: PCE   ! moisture transfer coefficient at PZREF
REAL, DIMENSION(:), INTENT(OUT)      :: PRI   ! bulk Richardson number
REAL, DIMENSION(:), INTENT(OUT)      :: PRESA ! aerodynamical resistance
REAL, DIMENSION(:), INTENT(OUT)      :: PZ0HSEA ! heat roughness length - what for ?
!                                             ! We already have a better estimate with ZOT
!
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PTA))      :: ZVMOD    ! wind intensity at PUREF atm level (m/s)
REAL, DIMENSION(SIZE(PTA))      :: ZPA      ! pressure at PZREF atm. level (Pa) 
REAL, DIMENSION(SIZE(PTA))      :: ZTA      ! temperature at PZREF atm. level (Kelvin)
REAL, DIMENSION(SIZE(PTA))      :: ZQASAT   ! specific humidity at saturation at PZREF (kg/kg)
!
REAL, DIMENSION(SIZE(PTA))      :: ZO       ! aerodynamic roughness (m)
REAL, DIMENSION(SIZE(PTA))      :: ZWG      ! gustiness correction (m/s)
!
REAL, DIMENSION(SIZE(PTA))      :: ZDU,ZDT,ZDQ,ZDUWG !differences between meterological parameters at
!                                              ! PUREF/PZREF and ocean
!
REAL, DIMENSION(SIZE(PTA))      :: ZUSR        ! velocity scaling parameter "ustar" (m/s) = friction velocity
REAL, DIMENSION(SIZE(PTA))      :: ZTSR        ! temperature sacling parameter "thetastar" (degC)
REAL, DIMENSION(SIZE(PTA))      :: ZQSR        ! humidity scaling parameter "qstar" (kg/kg)
!
REAL, DIMENSION(SIZE(PTA))      :: ZU10,ZT10   ! neutral wind speed and temperature at 10m  
REAL, DIMENSION(SIZE(PTA))      :: ZVISA       ! kinematic viscosity of dry air
REAL, DIMENSION(SIZE(PTA))      :: ZO10,ZOT10  ! first guess for roughness lengths from neutral profiles
REAL, DIMENSION(SIZE(PTA))      :: ZCD         ! first guess of neutral momentum transfer coef. at PUREF
REAL, DIMENSION(SIZE(PTA))      :: ZCT         ! first guess of neutral heat transfer coef. at PZREF
REAL, DIMENSION(SIZE(PTA))      :: ZCC         ! C constant as in equation (10) of Grachev and Fairall (1997, JAM)
REAL, DIMENSION(SIZE(PTA))      :: ZCD10       ! first guess of neutral momentum transfer coef. at 10m
REAL, DIMENSION(SIZE(PTA))      :: ZCT10       ! first guess of neutral thermal component of heat transfer coef. at 10m 
                                               ! (see Fairall et al 1996, section 2.1, equation (5) & (9))
REAL, DIMENSION(SIZE(PTA))      :: ZRIBU,ZRIBCU ! Richardson bulk number Rib, saturation Richardson number 
REAL, DIMENSION(SIZE(PTA))      :: ZETU        ! Estimate of Monin-Obukhov stability parameter zeta = z/L
                                               ! with z=PUREF from Richardson bulk number (for initial guess)
REAL, DIMENSION(SIZE(PTA))      :: ZL10        ! Monin-Obukhov stability length L estimated from Richardson
                                               ! bulk number (initial guess)
!
REAL, DIMENSION(SIZE(PTA))      :: ZCHARN                      ! Charnock parameter depends on neutral wind module
REAL, DIMENSION(SIZE(PTA))      :: ZTWAVE,ZHWAVE,ZCWAVE,ZLWAVE ! parameters of gravity wave models
!
REAL, DIMENSION(SIZE(PTA))      :: ZZL,ZZTL!,ZZQL    ! Monin-Obukhov stability parameter zeta = z/L for U,T and Q
                                                     ! inside the iterative loop
REAL, DIMENSION(SIZE(PTA))      :: ZRR               ! roughness Reynolds number inside the iterative loop
REAL, DIMENSION(SIZE(PTA))      :: ZOT,ZOQ           ! scalar roughness lengths inside the iterative loop
REAL, DIMENSION(SIZE(PTA))      :: ZPUZ,ZPTZ,ZPQZ    ! stability correction PSI for U,T,Q inside the iterative loop
!
REAL, DIMENSION(SIZE(PTA))      :: ZBF               ! g/T * Q0v (to be used for gustiness correction)
!
REAL, DIMENSION(SIZE(PTA))      :: ZTAU       ! upward momentum flux (W/m2)
REAL, DIMENSION(SIZE(PTA))      :: ZHF        ! upward sensible heat flux (W/m2)
REAL, DIMENSION(SIZE(PTA))      :: ZEF        ! upward latent heat flux (W/m2)
REAL, DIMENSION(SIZE(PTA))      :: ZWBAR      ! Webb correction to latent heat flux but not used here after
REAL, DIMENSION(SIZE(PTA))      :: ZTAUR      ! momentum flux due to rain (W/m2)
REAL, DIMENSION(SIZE(PTA))      :: ZRF        ! sensible heat flux due to rain (W/m2)
REAL, DIMENSION(SIZE(PTA))      :: ZCHN,ZCEN  ! neutral transfer coefficients for heat and vapor at 10m
!
REAL, DIMENSION(SIZE(PTA))      :: ZLV        ! latent heat of vaporisation at surface temperature
!
REAL, DIMENSION(SIZE(PTA))      :: ZTAC       ! atmospheric temperature in Celsius
REAL, DIMENSION(SIZE(PTA))      :: ZDQSDT     ! Clausius-Clapeyron relation for dqsat/dT
REAL, DIMENSION(SIZE(PTA))      :: ZDTMP,ZDWAT! heat and water vapour diffusivity in air
REAL, DIMENSION(SIZE(PTA))      :: ZALFAC     ! alpha factor in formula for sensible heat flux due to rain
REAL, DIMENSION(SIZE(PTA))      :: ZXLR       ! latent heat of vaporisation at rain temperature
REAL, DIMENSION(SIZE(PTA))      :: ZCPLW      ! specific heat for water at rain temperature
!
REAL, DIMENSION(SIZE(PTA))      :: ZUSTAR2    ! final square of friction velocity
!
REAL, DIMENSION(SIZE(PTA))      :: ZDIRCOSZW  ! orography slope cosine (=1 on water!)
REAL, DIMENSION(SIZE(PTA))      :: ZAC        ! aerodynamical conductance
!
!
INTEGER, DIMENSION(SIZE(PTA))   :: ITERMAX             ! maximum number of iterations
!
REAL    :: ZRVSRDM1,ZRDSRV,ZR2 ! thermodynamic constants
REAL    :: ZBETAGUST           ! gustiness factor
REAL    :: ZZBL                ! atm. boundary layer depth (m)
REAL    :: ZVISW               ! m2/s kinematic viscosity of water
REAL    :: ZS                  ! reference height for atmospheric profiles (10m)
REAL    :: ZCH10               ! first guess of heat transfer coef. at 10m
!
INTEGER :: J, JLOOP    !loop indice
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!       1.     Initializations
!              ---------------
!
!       1.1   Constants and parameters
!
IF (LHOOK) CALL DR_HOOK('COARE30_FLUX',0,ZHOOK_HANDLE)
!
ZRVSRDM1  = XRV/XRD-1. ! 0.607766
ZRDSRV    = XRD/XRV    ! 0.62198
ZR2       = 1.-ZRDSRV  ! pas utilisé dans cette routine
ZBETAGUST = 1.25       ! value based on TOGA-COARE experiment
ZZBL      = 600.       ! Set a default value for boundary layer depth
ZS        = 10.        ! Standard heigth = 10m
ZCH10     = 0.00115    ! Heat transfer coefficient at 10m
!
ZVISW     = 1.E-6
!
!       1.2   Array initialization by undefined values
!
PSFTH (:)=XUNDEF
PSFTQ (:)=XUNDEF
PUSTAR(:)=XUNDEF
!
PCD(:) = XUNDEF
PCDN(:) = XUNDEF
PCH(:) = XUNDEF
PCE(:) =XUNDEF
PRI(:) = XUNDEF
!
PRESA(:)=XUNDEF
!
!       1.3     Temperature 
!
! Set a non-zero value for the temperature gradient
!
WHERE((PTA(:)*PEXNS(:)/PEXNA(:)-PSST(:))==0.) 
      ZTA(:)=PTA(:)-1E-3
ELSEWHERE
      ZTA(:)=PTA(:)      
ENDWHERE

!       1.4     Wind and humidity 
!
! Sea surface specific humidity 
!
PQSAT(:)=QSAT_SEAWATER(PSST(:),PPS(:))         
!              
! Set a minimum value to wind 
!
ZVMOD(:) = WIND_THRESHOLD(PVMOD(:),PUREF(:))
!
! Specific humidity at saturation at the atm. level 
!
ZPA(:) = XP00* (PEXNA(:)**(XCPD/XRD))
ZQASAT(:) = QSAT(ZTA(:),ZPA(:)) 
!
!-------------------------------------------------------------------------------
!
!        2. INITIAL GUESS FOR THE ITERATIVE METHOD
!        -----------------------------------------
!
!        2.1    A few parameters
!
ZO(:)  = 0.0001              ! first guess for aerodynamic roughness length
ZWG(:) = 0.                  ! first guess for gustiness correction
IF (S%LPWG) ZWG(:) = 0.5
!
DO J=1,SIZE(PTA)
  !
  !      2.2    Atmospheric profiles
  !    
  ZDU(J) = ZVMOD(J)          ! wind speed difference with surface (=0) (m/s)
  ZDT(J) = -(ZTA(J)/PEXNA(J)) + (PSST(J)/PEXNS(J)) ! potential temperature difference
  ZDQ(J) = PQSAT(J)-PQA(J)                         ! specific humidity difference
  ! defined once and for all
  !
  ZDUWG(J) = SQRT(ZDU(J)**2+ZWG(J)**2)     ! wind speed difference including gustiness correction
  ! gustiness correction will depend on turbulent heat flux and change at each iteration
  !
  !      2.3   First guess of neutral profiles and neutral transfer coefficients
  !
  ZU10(J)  = ZDUWG(J)*LOG(ZS/ZO(J))/LOG(PUREF(J)/ZO(J)) ! first guess of neutral wind speed at 10m
  ZUSR(J)  = 0.035*ZU10(J)                              ! first guess of ustar
  ZVISA(J) = 1.326E-5*(1.+6.542E-3*(ZTA(J)-XTT)+&
             8.301E-6*(ZTA(J)-XTT)**2-4.84E-9*(ZTA(J)-XTT)**3) !Andrea (1989) CRREL Rep. 89-11
  ! 
  ZCHARN(J) = MAX(0.011,MIN(0.018,0.011+0.007*(ZU10(J)-10.)/8.)) ! Charnock parameter varies between 0.011
  ! and 0.018 as a function of the neutral wind speed at 10m according to Fairall et al (2003), see section 3c
  ZO10(J) = ZCHARN(J)*ZUSR(J)*ZUSR(J)/XG+0.11*ZVISA(J)/ZUSR(J) ! second guess for aerodynamic roughness length
                                                               ! based on Charnock model and neutral profile
  ZCD(J)  = (XKARMAN/LOG(PUREF(J)/ZO10(J)))**2   ! neutral momentum transfer coefficient at PUREF height
  ZCD10(J)= (XKARMAN/LOG(ZS/ZO10(J)))**2         ! neutral momentum transfer coefficient at 10m
  ZCT10(J)= ZCH10/SQRT(ZCD10(J))                 ! neutral thermal component of heat transfer coef. at 10m 
                                                 ! (see Fairall et al 1996, section 2.1, equation (5) & (9))
  ZOT10(J)= ZS/EXP(XKARMAN/ZCT10(J))             ! first guess of scalar roughness length estimated from first
                                                 ! guess of neutral profile
  ZCT(J) = XKARMAN/LOG(PZREF(J)/ZOT10(J))        ! thermal component of neutral heat transfer coefficient at 
                                                 ! PZREF height ((see Fairall et al 1996, sec. 2.1, eq. (5) & (9))
  !
  !      2.4   Estimate of zeta=z/L from bulk Richardson number to speed up subsequent iterative loop
  !
  ZCC(J) = XKARMAN*ZCT(J)/ZCD(J)                 ! C constant as in equation (10) of Grachev and Fairall (1997, JAM)
                                                 ! used in their zeta = f (Rib) model
  !
  ZRIBCU(J) = -PUREF(J)/(ZZBL*0.003*ZBETAGUST**3) !saturation Richardson Bulk number Ribc according to equation
                                                 ! (24) in Grachev and Fairall (1997, JAM).
                                                 ! 
  ! Warning : Grachev and Fairall (1997) assume temperature, humidity and winds to be taken at the same
  ! level in ZCC, ZRIBCU, ZRIBU, ZETU, which is not necessarily the case in this routine (PUREF, PZREF).
  ! The bulk Richardson number should contain an additional PUREF/PZREF factor but the adaptation of
  ! Grachev and Fairall (1997) model to account for two different levels for U and T/Q should also
  ! include an additional PZREF/PUREF in ZCC and PUREF/PZREF in RIBCU which would cancel out.
  ZRIBU(J)  = -XG*PUREF(J)*(ZDT(J)+ZRVSRDM1*ZTA(J)*ZDQ(J))/&
               (ZTA(J)*ZDU(J)**2)              ! Richardson Bulk number Rib approximation according to 
                                                 ! Grachev and Fairall (1997, JAM) equation (11)
  !
  IF (ZRIBU(J)<0.) THEN
    ZETU(J) = ZCC(J)*ZRIBU(J)/(1.+ZRIBU(J)/ZRIBCU(J))    ! Equation (25) in Grachev and Fairall (1997, JAM)
                                                         ! for unstable cases
  ELSE
    ZETU(J) = ZCC(J)*ZRIBU(J)/(1.+27./9.*ZRIBU(J)/ZCC(J))! Equation for stable cases - literature ?
  ENDIF
  !
  ZL10(J) = PUREF(J)/ZETU(J)                     ! Monin-Obukhov stability length L
  !
ENDDO
  !
  !      2.5   First guess of u*,theta*,q* scaling parameters accounting for stability effects
  !
ZUSR(:) = ZDUWG(:)*XKARMAN/(LOG(PUREF(:)/ZO10(:))-PSIFCTU(PUREF(:)/ZL10(:)))
ZTSR(:) = -ZDT(:)*XKARMAN/(LOG(PZREF(:)/ZOT10(:))-PSIFCTT(PZREF(:)/ZL10(:)))
ZQSR(:) = -ZDQ(:)*XKARMAN/(LOG(PZREF(:)/ZOT10(:))-PSIFCTT(PZREF(:)/ZL10(:)))
ZO(:) = ZO10(:) ! second guess of aerodynamic roughness length
!
ZZL(:) = 0.0   ! Why initializing this to 0 when ZL10 is available as well as PUREF and PZREF ?
!
DO J=1,SIZE(PTA)
  !
  IF (ZETU(J)>50.) THEN
    ITERMAX(J) = 1
  ELSE
    ITERMAX(J) = 3 ! maximum number of iterations, only 3 since estimation of zeta from Rib
                   ! speeds up convergence (Fairall et al, 2003, section 3)
  ENDIF
  !
  !
ENDDO
!
!-------------------------------------------------------------------------------
!
  !      3.  ITERATIVE LOOP TO COMPUTE USR, TSR, QSR
  !      -------------------------------------------
  !
DO JLOOP=1,MAXVAL(ITERMAX) ! begin iterative loop
  !
  DO J=1,SIZE(PTA)
    !
    IF (JLOOP.GT.ITERMAX(J)) CYCLE
  !
  !      3.1   Aerodynamic rougness length
  !
    ZU10(J)  = ZUSR(J)/XKARMAN*LOG(ZS/ZO(J)) ! neutral wind speed at 10m required for
                                             ! all roughness length models

    IF (S%NGRVWAVES==0) THEN
      ZCHARN(J) = MAX(0.011,MIN(0.018,0.011+0.007*(ZU10(J)-10.)/8.))
            ! Smith (1988) adapted by Fairall et al (2003) with a varying ZCHARN
            ! according to neutral wind speed (see section 3c)
      ZO(J) = ZCHARN(J)*ZUSR(J)*ZUSR(J)/XG + 0.11*ZVISA(J)/ZUSR(J)
    ELSE
            ! Parameters for gravity wave models which depend on neutral wind speed
      ZHWAVE(J) = 0.018*ZU10(J)*ZU10(J)*(1.+0.015*ZU10(J))
                                     ! hs used in equation (23) in Fairall et al (2003)
                                     ! Where does this formula comes from ?
      ZTWAVE(J) = 0.729*ZU10(J)      ! Tp in equation (27) from Fairall et al (2003)
                                     ! neutral wind should be used instead of ZVMOD
      ZCWAVE(J) = XG*ZTWAVE(J)/(2.*XPI) ! Cp in equation (26) from Fairall et al (2003)
      ZLWAVE(J) = ZTWAVE(J)*ZCWAVE(J)   ! Lp in equation (26) from Fairall et al (2003)

      IF (S%NGRVWAVES==1) THEN
            ! OOst et al (2002) accounting for gravity waves
            ! Equation (25b) in Fairall et al (2003)
        ZO(J) = (50./(2.*XPI))*ZLWAVE(J)*(ZUSR(J)/ZCWAVE(J))**4.5 &
              + 0.11*ZVISA(J)/ZUSR(J)
      ELSE IF (S%NGRVWAVES==2) THEN
            ! Taylor and Yelland (2001) accounting for gravity waves
            ! Equation (25a) in Fairall et al (2003)
        ZO(J) = 1200.*ZHWAVE(J)*(ZHWAVE(J)/ZLWAVE(J))**4.5 &
              + 0.11*ZVISA(J)/ZUSR(J)
      ENDIF
    ENDIF
    !
    !    3.2   Scalar rougness lengths
    !
    ZRR(J) = ZO(J)*ZUSR(J)/ZVISA(J)            ! Roughness Reynolds number
    ZOQ(J) = MIN(1.1E-4 , 5.5E-5/ZRR(J)**0.6)  ! Equation (28) in Fairall et al (2003)
    ZOT(J) = ZOQ(J)                            ! Idem
    !
    !    3.3   Stability effects : zeta and psi
    !
    ZZL(J) = XKARMAN * XG * PUREF(J) * &       ! zeta = z/L for U (at PUREF)
              ( ZTSR(J)*(1.+ZRVSRDM1*PQA(J)) + ZRVSRDM1*ZTA(J)*ZQSR(J) ) / &
              ( ZTA(J)*ZUSR(J)*ZUSR(J)*(1.+ZRVSRDM1*PQA(J)) )  
    ! L uses here beta=g/Tv instead of beta=g/thetav and thetavstar =
    ! thetastar*(1+0.61*q) + 0.61*T*qstar instead of thetavstar =
    ! thetastar*(1+0.61*q) + 0.61*theta*qstar but Fairall et al (2003) makes
    ! even more approximation (beta=g/T, thetavstar = thetastar+0.61*T*qstar)
    ZZTL(J)= ZZL(J)*PZREF(J)/PUREF(J)          ! zeta = z/L for T & Q (at PZREF)
!    ZZQL(J)=ZZL(J)*PZREF(J)/PUREF(J)  ! for Q
  ENDDO
  !
  ZPUZ(:) = PSIFCTU(ZZL(:))    ! Stability correction psi for wind profile
  ZPTZ(:) = PSIFCTT(ZZTL(:))   ! Stability correction psi for T and Q profiles
  !
  DO J=1,SIZE(PTA)
    !
    ! ZPQZ(J)=PSIFCTT(ZZQL(J))    
    ZPQZ(J) = ZPTZ(J)
    !
    !    3.4   Updated estimates of ustar, thetastar, qstar scaling parameters
    !
    ZUSR(J) = ZDUWG(J)*XKARMAN/(LOG(PUREF(J)/ZO(J)) -ZPUZ(J))
    ZTSR(J) = -ZDT(J) *XKARMAN/(LOG(PZREF(J)/ZOT(J))-ZPTZ(J))
    ZQSR(J) = -ZDQ(J) *XKARMAN/(LOG(PZREF(J)/ZOQ(J))-ZPQZ(J))
    !
    !    3.5   Updated estimate of gustiness correction (ZWG)
    !
    IF(S%LPWG) THEN
      ZBF(J) = -XG/ZTA(J)*ZUSR(J)*(ZTSR(J)+ZRVSRDM1*ZTA(J)*ZQSR(J))
      IF (ZBF(J)>0.) THEN
        ZWG(J) = ZBETAGUST*(ZBF(J)*ZZBL)**(1./3.)
      ELSE
        ZWG(J) = 0.2
      ENDIF
    ENDIF  
    ZDUWG(J) = SQRT(ZVMOD(J)**2 + ZWG(J)**2)
    !
  ENDDO
  !
ENDDO
!
!-------------------------------------------------------------------------------
!
!        4.  COMPUTE transfer coefficients PCD, PCH, ZCE and SURFACE FLUXES
!        ------------------------------------------------------------------
!
ZTAU(:) = XUNDEF
ZHF(:)  = XUNDEF
ZEF(:)  = XUNDEF
!
ZWBAR(:) = 0.
ZTAUR(:) = 0.
ZRF(:)   = 0.
!
DO J=1,SIZE(PTA)
  !
  !
  !      4.1 Transfert coefficients PCD at PUREF and PCH and PCE at PZREF
  !
  PCD(J) = (ZUSR(J)/ZDUWG(J))**2.
  PCH(J) = ZUSR(J)*ZTSR(J)/(ZDUWG(J)*(-ZDT(J)))
  PCE(J) = ZUSR(J)*ZQSR(J)/(ZDUWG(J)*(PQA(J)-PQSAT(J)))
  !
  !      4.2 Neutral transfer coefficients at 10m
  !          (10m is set in section 1.1, l. 220)
  !
  PCDN(J) = (XKARMAN/LOG(ZS/ZO(J)))**2.
  ZCHN(J) = (XKARMAN/LOG(ZS/ZO(J)))*(XKARMAN/LOG(ZS/ZOT(J)))
  ZCEN(J) = (XKARMAN/LOG(ZS/ZO(J)))*(XKARMAN/LOG(ZS/ZOQ(J)))
  !
  !      4.3 Surface fluxes
  !
  ZLV(J) = XLVTT + (XCPV-XCL)*(PSST(J)-XTT) ! Bolton (1980) linear approximation
                                            ! for latent heat of vaporisation
  !
  IF (ABS(PCDN(J))>1.E-2) THEN   !!!! secure COARE3.0 CODE 
    write(*,*) 'pb PCDN in COARE30: ',PCDN(J)
    write(*,*) 'point: ',J,"/",SIZE(PTA)
    write(*,*) 'roughness: ', ZO(J)
    write(*,*) 'ustar: ',ZUSR(J)
    write(*,*) 'wind: ',ZDUWG(J)
    CALL ABOR1_SFX('COARE30: PCDN too large -> no convergence')
  ELSE
    ZTSR(J) = -ZTSR(J)
    ZQSR(J) = -ZQSR(J)
    ZTAU(J) = -PRHOA(J)*ZUSR(J)*ZUSR(J)*ZVMOD(J)/ZDUWG(J) ! upward momentum flux
    ! Correction by U/S because ustar^2 = Cd*S*U while ustart^2 = Cd*S^2
    ! where S is the wind corrected for gustiness (eqs 9 & 18 in Fairall et al, 1996)
    ZHF(J)  =  PRHOA(J)*XCPD*ZUSR(J)*ZTSR(J)              ! upward sensible heat flux
    ZEF(J)  =  PRHOA(J)*ZLV(J)*ZUSR(J)*ZQSR(J)            ! upward latent heat flux
    !    
    !    4.4 Contribution of precipitation to surface  fluxes
    !
    IF (S%LPRECIP) THEN
      ! 
      !  4.4.a Sensible heat flux
      !
      ZTAC(J)  = ZTA(J)-XTT                    ! atmospheric temperature in Celsius degrees
      !
      ZXLR(J)  = XLVTT + (XCPV-XCL)* ZTAC(J)   ! latent heat of vaporisation at rain temperature
      ZDQSDT(J)= ZQASAT(J) * ZXLR(J) / (XRV*ZTA(J)**2)                  ! Clausius-Clapeyron relation
      ZDTMP(J) = (1.0 + 3.309e-3*ZTAC(J) -1.44e-6*ZTAC(J)*ZTAC(J)) * &  ! heat diffusivity in air
                  0.02411 / (PRHOA(J)*XCPD)
      !
      ZDWAT(J) = 2.11e-5 * (XP00/ZPA(J)) * (ZTA(J)/XTT)**1.94           ! water vapour diffusivity from eq (13.3)
      !                                                                 ! of Pruppacher and Klett (1978)      
      ZALFAC(J)= 1.0 / (1.0 + &                                         ! Equation (11) in Gosnell et al (1995)
                   ZRDSRV*ZDQSDT(J)*ZXLR(J)*ZDWAT(J)/(ZDTMP(J)*XCPD))   ! Clausius-Clapeyron wet-bulb factor (dimensionless)
      ! ZRDSRV = Rd/Rv might be there to account for humidity effect but I do not understand where it comes from.
      ! On limitation in Gosnell et al (1995) would rather be the use of Cpd. We could replace by
      ! Cpd * (1-q) + Cpv * q. I would remove the ZRDSRV.
      ZCPLW(J) = 4224.8482 + ZTAC(J) * &
                              ( -4.707 + ZTAC(J) * &
                                (0.08499 + ZTAC(J) * &
                                  (1.2826e-3 + ZTAC(J) * &
                                    (4.7884e-5 - 2.0027e-6* ZTAC(J))))) ! specific heat of water at rain temperature
      !       
      ZRF(J)   = PRAIN(J) * ZCPLW(J) * ZALFAC(J) * &                    ! Equation (23) in Fairall et al (1996)
                   (PSST(J) - ZTA(J) + (PQSAT(J)-PQA(J))*ZXLR(J)/XCPD ) ! Equation (12) in Gosnell et al (1995)
      ! Warning : the exact formula in Gosnell et al (1995) and Fairall et al (1996) gives -ZRF and not ZRF
      ! because it is the heat transferred to the ocean by the rain. Here the convention is upward heat flux.
      !
      !  4.4.b Momentum flux
      !
      ZTAUR(J)=-0.85*(PRAIN(J) *ZVMOD(J)) ! Equation (8) in Caldwell and Elliott (1971) discussed in last paragraph of
                                          ! section 2.5 p 3752 in Fairall et al (1996) - water density is implicit
      ! Warning : ZVMOD needs to be taken at 10m height for this formula to be valid and the formula gives -ZTAUR
      !
    ENDIF
    !
    !    4.5 Webb correction to latent heat flux
    ! 
    ZWBAR(J)=- (1./ZRDSRV)*ZUSR(J)*ZQSR(J) + (1.0+(1./ZRDSRV)*PQA(J))* &   ! Equation (21) in Fairall et al (1996)
               (- ZUSR(J)*ZTSR(J)/ZTA(J))         
    !
    !    4.6 Friction velocity which contains correction due to rain
    !
    ZUSTAR2(J)= - (ZTAU(J) + ZTAUR(J)) / PRHOA(J)
    PUSTAR(J) =  SQRT(ZUSTAR2(J))
    ! Warning : the rain effect ZTAUR acts on the ocean but not on the atmosphere
    !
    !    4.7 Total surface fluxes
    !           
    PSFTH (J) =  ZHF(J) + ZRF(J)   ! Upward sensible heat flux
    PSFTQ (J) =  ZEF(J) / ZLV(J)   ! Upward water flux
    ! 
  ENDIF
ENDDO                      
!-------------------------------------------------------------------------------
!
!        5.  FINAL STEP : TOTAL SURFACE FLUXES AND DERIVED DIAGNOSTICS
!        -------------------------------------------------------------
!
!        5.1    Richardson number
!             
!
ZDIRCOSZW(:) = 1.
 CALL SURFACE_RI(PSST,PQSAT,PEXNS,PEXNA,ZTA,ZQASAT,&
                PZREF,PUREF,ZDIRCOSZW,PVMOD,PRI   )  
!
!        5.2     Aerodynamical conductance and resistance
!             
ZAC(:) = PCH(:)*ZVMOD(:)
PRESA(:) = 1. / MAX(ZAC(:),XSURF_EPSILON)
!
!        5.3 Z0 and Z0H over sea
!
PZ0SEA(:) =  ZCHARN(:) * ZUSTAR2(:) / XG + XVZ0CM * PCD(:) / PCDN(:)
! Why ? We already have Z0 et Z0T which are more robust estimates.
PZ0HSEA(:) = PZ0SEA(:)
!
IF (LHOOK) CALL DR_HOOK('COARE30_FLUX',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE COARE30_FLUX
