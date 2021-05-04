!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
      SUBROUTINE ICE_SEA_FLUX(PZ0ICE,                                       &
                              PTA, PEXNA, PRHOA, PTICE, PEXNS, PQA, PRR, PRS, &
                              PVMOD, PZREF, PUREF,                            &
                              PPS, PQSAT,                                     &
                              PSFTH, PSFTQ, PUSTAR,                           &
                              PCD, PCDN, PCH, PRI, PRESA, PZ0HICE             )  
!     #######################################################################
!
!
!!****  *ICE_SEA_FLUX*  
!!
!!    PURPOSE
!!    -------
!      Calculate the surface fluxes of heat, moisture, and momentum over
!       sea ice. adapted from WATER_FLUX  
!     
!!**  METHOD
!!    ------
!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!!    XCD_ICE_CST, from MODD_SEAFLUX
!!
!!      
!!    REFERENCE
!!    ---------
!!      
!!    AUTHOR
!!    ------
!!      S. Belair           * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original      01/09/95 
!!      (J.Stein)     16/11/95  use PUSLOPE and Theta to compute Ri
!!      (P.Lacarrere) 19/03/96  bug in the ZTHVI and ZTHVIS computations
!!      (J.Stein)     27/03/96  use only H and LE in the soil scheme
!!      (P.Jabouille) 12/11/96  bug in the Z0 computation
!!      (V.Masson)    01/02/00  detection of sea ice
!!      (P. Tulet)    01/10/03  aerodynamical resistance output
!!      (P. LeMoigne) 29/03/04  bug in the heat flux computation
!!      (P. LeMoigne) 09/02/06  Z0H as output
!!      (B. Decharme)    06/09  limitation of Ri
!!      (B. Decharme)    09/09  limitation of Ri in surface_ri.F90
!!      (S.Senesi)       01/14  use XCD_ICE_CST (if /= 0) as value for for Cd, Cdn and Ch
!!      (V.Guemas)       05/21  Output scalar and aerodynamic roughness consistent
!!                              with XCD_ICE_CST when /=0
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_CSTS,       ONLY : XG, XCPD, XKARMAN
USE MODD_SURF_PAR,   ONLY : XUNDEF
USE MODD_SURF_ATM,   ONLY : LDRAG_COEF_ARP, LRRGUST_ARP, XRRSCALE, &
                            XRRGAMMA, XUTILGUST     
USE MODD_SNOW_PAR,   ONLY : XZ0SN, XZ0HSN
USE MODN_SEAFLUX_n,  ONLY : XCD_ICE_CST
!
USE MODI_SURFACE_RI
USE MODI_SURFACE_AERO_COND
USE MODI_SURFACE_CD
USE MODI_SURFACE_CDCH_1DARP
USE MODI_WIND_THRESHOLD
!
USE MODE_THERMOS
!
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
REAL, DIMENSION(:), INTENT(IN)       :: PTA   ! air temperature at PZREF atm. level (K)
REAL, DIMENSION(:), INTENT(IN)       :: PQA   ! air humidity at PZREF atm. level (kg/kg)
REAL, DIMENSION(:), INTENT(IN)       :: PEXNA ! Exner function at PZREF atm. level
REAL, DIMENSION(:), INTENT(IN)       :: PRHOA ! air density at PZREF atm. level (kg/m3)
REAL, DIMENSION(:), INTENT(IN)       :: PVMOD ! module of wind at PUREF atm. wind level (m/s)
REAL, DIMENSION(:), INTENT(IN)       :: PZREF ! atm. level for temp. and humidity (m)
REAL, DIMENSION(:), INTENT(IN)       :: PUREF ! atm. level for wind (m)
REAL, DIMENSION(:), INTENT(IN)       :: PTICE ! Sea ice Surface Temperature (K)
REAL, DIMENSION(:), INTENT(IN)       :: PEXNS ! Exner function at ice surface
REAL, DIMENSION(:), INTENT(IN)       :: PPS   ! air pressure at ice and sea surface (Pa)
REAL, DIMENSION(:), INTENT(IN)       :: PRR   ! rain rate (kg/s/m2)
REAL, DIMENSION(:), INTENT(IN)       :: PRS   ! snow rate (kg/s/m2)
!
REAL, DIMENSION(:), INTENT(INOUT)    :: PZ0ICE! aerodynamical roughness length over sea ice
!                                         
!                                         
!  surface fluxes : latent heat, sensible heat, momentum fluxes
REAL, DIMENSION(:), INTENT(OUT)      :: PSFTH ! upward heat flux  (W/m2)
REAL, DIMENSION(:), INTENT(OUT)      :: PSFTQ ! upward water flux (kg/m2/s)
REAL, DIMENSION(:), INTENT(OUT)      :: PUSTAR! friction velocity (m/s)
!
!  diagnostics
REAL, DIMENSION(:), INTENT(OUT)      :: PQSAT ! near ice saturation specific humidity
REAL, DIMENSION(:), INTENT(OUT)      :: PCD   ! momentum transfer coefficient at PUREF
REAL, DIMENSION(:), INTENT(OUT)      :: PCDN  ! neutral momentum transfer coefficient at PUREF
REAL, DIMENSION(:), INTENT(OUT)      :: PCH   ! heat transfer coefficient at PZREF
REAL, DIMENSION(:), INTENT(OUT)      :: PRI   ! bulk Richardson number
REAL, DIMENSION(:), INTENT(OUT)      :: PRESA   ! aerodynamical resistance
REAL, DIMENSION(:), INTENT(OUT)      :: PZ0HICE ! scalar roughness length over sea ice
!
!
!*      0.2    declarations of local variables
!
!
REAL, DIMENSION(SIZE(PTA)) :: ZVMOD     ! wind modulus
REAL, DIMENSION(SIZE(PTA)) :: ZUSTAR2   ! square of friction velocity
REAL, DIMENSION(SIZE(PTA)) :: ZAC       ! Aerodynamical conductance
REAL, DIMENSION(SIZE(PTA)) :: ZRA       ! Aerodynamical resistance
REAL, DIMENSION(SIZE(PTA)) :: ZDIRCOSZW ! orography slope cosine (=1 on water!)
REAL, DIMENSION(SIZE(PTA)) :: ZFP       ! working variable
REAL, DIMENSION(SIZE(PTA)) :: ZRRCOR    ! correction od CD, CH, CDN due to moist-gustiness
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!
!-------------------------------------------------------------------------------
!
!       1.     Initializations
!              ---------------
!
IF (LHOOK) CALL DR_HOOK('ICE_SEA_FLUX',0,ZHOOK_HANDLE)
ZDIRCOSZW=1.
!
PRI(:) = XUNDEF
PCH(:) = XUNDEF
PCD(:) = XUNDEF
PCDN(:) = XUNDEF
!
PSFTH (:)=XUNDEF
PSFTQ (:)=XUNDEF
PUSTAR(:)=XUNDEF
PRESA(:)=XUNDEF
!
!
!       1.1    Saturation specific humidity above sea ice
!              -------------------------------------------
!
PQSAT(:) = QSAT(PTICE(:),PPS(:))
!
!       1.2    Minimum wind value depending on height
!              ---------------------------------------
ZVMOD(:)=WIND_THRESHOLD(PVMOD(:),PUREF(:))
!
!       1.3    Calculate the bulk Richardson number
!              -------------------------------------
CALL SURFACE_RI(PTICE,PQSAT,PEXNS,PEXNA,PTA,PQA, &
                     PZREF, PUREF, ZDIRCOSZW,PVMOD,PRI)
!
!-------------------------------------------------------------------------------
!
!       2.     Transfer coefficients for momentum, heat and moisture
!              -----------------------------------------------------
!
!
IF ( XCD_ICE_CST == 0.0 ) THEN ! Namelist parameter to allow backward compatibility
  !
  PZ0HICE(:) = XZ0HSN ! scalar roughness length which can be selected in
                      ! NAM_SURF_SNOW_CSTS. Default : 1E-4
  PZ0ICE (:) = XZ0SN  ! aerodynamic roughness length which can be selected in
                      ! NAM_SURF_SNOW_CSTS. Default : 1E-3
  !
  IF (LDRAG_COEF_ARP) THEN
     ! Computation of transfer coefficients for momentum, heat and moisture
     ! following the formulation which was in ARPEGE/ALADIN
     CALL SURFACE_CDCH_1DARP(PZREF, PZ0ICE, PZ0HICE , ZVMOD, PTA, PTICE, &
                             PQA, PQSAT, PCD, PCDN, PCH                 )  
!
     ZRA(:) = 1. / ( PCH(:) * ZVMOD(:) )
!
  ELSE

     ! Computation of transfer coefficient for momentum following formulations
     ! for land surfaces
     CALL SURFACE_CD(PRI, PZREF, PUREF, PZ0ICE, PZ0HICE, PCD, PCDN)
     ! Computation of transfer coefficient for heat and moisture following
     ! formulations for land surfaces
     CALL SURFACE_AERO_COND(PRI, PZREF, PUREF, ZVMOD, PZ0ICE, PZ0HICE, ZAC, ZRA, PCH)
!
  ENDIF
!
ELSE
!
! Using variable transfer coefficients is not appropriate on seaice 
! with simple bulk functions.
! A constant value (e.g. 1.5.e-3 ) is preferable, and used except if the  
! user request backward compatibility by setting XCD_ICE_CST to 0 (DEFAULT). 
!
   PCD (:)=XCD_ICE_CST
   PCDN(:)=XCD_ICE_CST
   PCH (:)=XCD_ICE_CST
   ZRA (:)=1./(PCH(:)*ZVMOD(:))

   PZ0ICE (:) = PUREF(:)/ EXP( SQRT( XKARMAN**2/XCD_ICE_CST ))
   ! Simplified version of 
   !PZ0ICE (:) = PUREF(:)/ EXP( SQRT( XKARMAN**2/PCDN(:) ))
   !
   PZ0HICE (:) = PZREF(:) / EXP(XKARMAN/SQRT(XCD_ICE_CST)) 
   ! Simplified version of :
   !PZ0HICE (:) = PZREF(:) / EXP( XKARMAN*SQRT(PCDN(:))/PCH(:)) 
              ! In the same way as CDN is assimilated to CD, CHN is assilated to CH
!
ENDIF
!
ZUSTAR2(:) = PCD(:)*ZVMOD(:)*ZVMOD(:)
!
PRESA(:) = ZRA(:)
!
IF (LRRGUST_ARP) THEN ! Correction of CD, CH, CDN due to moist gustiness
                      ! References ?
  ZFP(:)=MAX(0.0,PRR(:)+PRS(:))  ! Total precipitation
  ZRRCOR(:)=SQRT(1.0+((((ZFP(:)/(ZFP(:)+XRRSCALE))**XRRGAMMA)*XUTILGUST)**2) &
      /(PCD(:)*ZVMOD(:)**2))  

  PCD  = PCD*ZRRCOR
  PCH  = PCH*ZRRCOR
  PCDN = PCDN*ZRRCOR
ENDIF
!
!-------------------------------------------------------------------------------
!
!       3.     Upward fluxes of heat and moisture, downward momentum flux
!              -----------------------------------------------------------
!
PSFTH (:) =  XCPD * PRHOA(:) * PCH(:) * ZVMOD(:) * ( PTICE(:) -PTA(:) * PEXNS(:) / PEXNA(:) ) / PEXNS(:)
! Using Heat transfer coefficient CH for vapor transfer coefficient CE !
PSFTQ (:) =  PRHOA(:) * PCH(:) * ZVMOD(:) * ( PQSAT(:)-PQA(:) )
! No Lv because this is the moisture flux and not the latent heat flux
PUSTAR(:) = SQRT(ZUSTAR2(:))
!
IF (LHOOK) CALL DR_HOOK('ICE_SEA_FLUX',1,ZHOOK_HANDLE)
!
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE ICE_SEA_FLUX
