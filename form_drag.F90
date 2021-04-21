!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
     SUBROUTINE FORM_DRAG (SM, PZ0SEA, PZ0ICE, PUREF, PLMOI, PLMOO, &
                           PRATIO, PCDF, PCHF)
!     #######################################################################
!
!
!!****  *FORM_DRAG*
!!
!!    PURPOSE
!!    -------
!      The momentum transfer coefficient over a mixture of ice and sea can be 
!      decomposed into : CD = XSIC * CDICE + (1-XSIC) * CDOCEAN + CDFORMDRAG
!      where CDICE and CDOCEAN are the skin drag coefficients, for ice
!      and ocean respectively, related to the surface characteristics
!      and CDFORMDRAG is the form drag induced by the ice topography or
!      the vertical faces formed by the alternating leads, floes and
!      melt ponds.
!
!      FORM_DRAG estimates the contribution of form drag to the
!      momentum, heat and moisture transfer coefficients.
!
!!**  METHOD
!!    ------
!      This routine applies the equations found in Lupkes and Gryanik
!      (2015, JGR)
!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!
!     Andreas, EL, Horst TW, Grachev AA, Persson POG, Fairall CW,
!       Guest PS, Jordan RE (2010) Parametrizing turbulent exchange over 
!       summer sea ice and the marginal ice zone, Q. J. R. Meteorol. Soc.,
!       138, 927–943.
!     Lüpkes C, Gryanik VM, Hartmann J, Andreas EL (2012), A parametrization, 
!       based on sea ice morphology, of the neutral atmospheric drag coefficients 
!       for weather prediction and climate models, J. Geophys. Res., 117, D13112, 
!       doi:10.1029/2012JD017630.
!     Lüpkes, C and Gryanik, VM (2015) A stability-dependent parametrization
!       of transfer coefficients for momentum and heat over polar sea ice to be 
!       used in climate models. J. Geophys. Res. Atm., 120, 552–581
!!
!!    AUTHOR
!!    ------
!!     V. Guemas  *Météo-France* 
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       04/2021
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_SURFEX_n, ONLY  : SEAFLUX_MODEL_t
USE MODN_SEAFLUX_n, ONLY : XCE, XBETAFORM  
                         ! XCE = Effective coefficient of resistance of an 
                         ! individual floe or pond edge
                         ! XBETAFORM = beta factor determining the shape of ZDF
USE MODE_COARE30_PSI     ! Stability correction functions psi as a function of zeta
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
TYPE(SEAFLUX_MODEL_t), INTENT(IN) :: SM    ! contains sea ice concentration
!
REAL, DIMENSION(:), INTENT(IN)   :: PZ0SEA ! skin drag coefficient over water
REAL, DIMENSION(:), INTENT(IN)   :: PZ0ICE ! skin drag coefficient over ice
REAL, DIMENSION(:), INTENT(IN)   :: PUREF  ! atmospheric height for wind speed
REAL, DIMENSION(:), INTENT(IN)   :: PLMOI  ! Monin-Obukhov length above ice
REAL, DIMENSION(:), INTENT(IN)   :: PLMOO  ! Monin-Obukhov length above ocean
REAL, DIMENSION(:), INTENT(IN)   :: PRATIO ! ratio between aerodynamic and scalar
                                           ! and aerodynamic roughness lengths
REAL, DIMENSION(:), INTENT(OUT)  :: PCDF   ! output form drag coefficient 
REAL, DIMENSION(:), INTENT(OUT)  :: PCHF   ! Equivalent form drag coefficient for heat
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PZ0SEA))    :: ZSIC                 ! sea ice concentration
REAL, DIMENSION(SIZE(PZ0SEA))    :: ZSC                  ! sheltering of atmospheric flow
                                 ! between obstacles
REAL, DIMENSION(SIZE(PZ0SEA))    :: ZDF                  ! characteristic length of floes
REAL, DIMENSION(SIZE(PZ0SEA))    :: ZHF                  ! ice freeboard
REAL, DIMENSION(SIZE(PZ0SEA))    :: ZPSIU                ! psi stability correction for
                                 ! momentum
REAL, DIMENSION(SIZE(PZ0SEA))    :: ZFI, ZFO             ! Louis correction of neutral 
                                 ! drag coefficients over ice and ocean
REAL, DIMENSION(SIZE(PZ0SEA))    :: ZCDFK                ! Part of neutral form drag 
REAL, DIMENSION(SIZE(PZ0SEA))    :: ZCDNFI, ZCDNFO       ! Neutral form drag coefficients
                                 ! contributions from flow influenced by ice and ocean
!
REAL :: ZHFC      ! mean ice freeboard (m)
REAL :: ZDMIN     ! minimum floe length
REAL :: ZE        !
!
!
!-------------------------------------------------------------------------------
! 
!       1.     Initialization
!              --------------
!
ZSIC = SM%S%XSIC
ZHFC      = 0.41  ! Mean value of REFLEX Fram strait data (Lupkes et al, 2012)
ZDMIN     = 8
ZE        = EXP(1.0)
!
!-------------------------------------------------------------------------------
! 
!       3.     Setting a few ice parameters
!              ----------------------------
!
ZHF (:) = ZHFC ! A constant value allows to recover Andreas et al (2010) unified
             ! quadratic form according to Lupkes et al al (2015) and this is 
             ! the value advised by Lupkes and Gryanik (2015) when the ice freeboard
             ! is not a model parameter. See their equation (40).
             ! When GELATO is activated, we can replace by GELATO value.

ZDF (:) = ZDMIN/(1-ZSIC(:))**XBETAFORM ! This simplified form for floe lengths
             ! would compensate for a sheltering function set to 1 that we need
             ! to recover Andreas et al (2010) unified quadratic form for CDF 
             ! valid both in the inner Arctic and in the marginal ice zone.
             ! This form is the one advised by Lupkes and Gryanik (2015),
             ! see their equation (41).

ZSC(:)  = 1  ! Lupkes et al (2012) advises to use the same sheltering function
             ! in the inner Arctic and the marginal ice zone. This sheltering 
             ! function accounts for the sheltering effect of upstream floes
             ! and ridges.
             ! A constant sheltering function set to 1 is advised by Lupkes and
             ! Gryanik (2015), see their equation (43).
!
!-------------------------------------------------------------------------------
! 
!       4.     Neutral form drag
!              -----------------
!
ZCDFK(:) = XCE/2 * ZSC(:)**2 * ZHF(:)/ZDF(:) * ZSIC(:)   
         ! Part of equation (21) in Lupkes and Gryanik (2015)
         ! This equation is the marginal ice zone one.
         ! The inner Arctic one should use the melt pond depth and melt pond
         ! diameters instead of floe freeboard and diameters as well as (1-ZSIC)
         ! instead of ZSIC. But the chosen parameters above for ZHF and ZDF
         ! allow to recover a unified form for inner Arctic and MIZ.
ZCDNFI(:) = ZCDFK(:) * ( LOG(ZHF(:) / (ZE*PZ0ICE(:))) / LOG(PUREF/PZ0ICE(:)) )**2
         ! Equation (21) in Lupkes and Gryanik (2015) for ice
ZCDNFO(:) = ZCDFK(:) * ( LOG(ZHF(:) / (ZE*PZ0SEA(:))) / LOG(PUREF/PZ0SEA(:)) )**2
         ! Equation (21) in Lupkes and Gryanik (2015) for water
!
!-------------------------------------------------------------------------------
! 
!       5.     Stability correction for momentum
!              ---------------------------------
!
ZPSIUO (:) = PSIFCTU(PUREF(:)/PLMOI(:))
!        ! Stability correction above ice
ZPSIUI (:) = PSIFCTU(PUREF(:)/PLMOO(:))
         ! Stability correction above ocean
!
ZFI (:)   = (1 - ZPSIU(:)/LOG(PUREF(:)/PZ0ICE(:)))**(-2)
         ! Equation (A1) in Lupkes and Gryanik (2015) for ice
ZFO (:)   = (1 - ZPSIU(:)/LOG(PUREF(:)/PZ0SEA(:)))**(-2) 
         ! Equation (A1) in Lupkes and Gryanik (2015) for ocean
!
!-------------------------------------------------------------------------------
! 
!       6.     Stability-dependant form drag
!              -----------------------------
!
PCDF(:) = ZCDFO(:) * ZFO(:) * (1-ZSIC(:)) + ZCDFI(:) * ZFI(:) * ZSIC(:) 
         ! Equation (38) in Lupkes and Gryanik (2015)
!
END SUBROUTINE FORM_DRAG
