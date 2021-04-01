!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
MODULE MODE_COARE30_PSI
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
USE MODD_CSTS ,ONLY : XPI
USE MODN_SEAFLUX_n, ONLY : CPSISTAB, CPSIUNSTAB, XGAMMA
!
 INTERFACE PSIFCTU
  MODULE PROCEDURE PSIFUNCTU
 END INTERFACE
 INTERFACE PSIFCTT
  MODULE PROCEDURE PSIFUNCTT
 END INTERFACE
!
 CONTAINS
!
!---------------------------------------------------------------------------------------
!
!#######################################################################################
FUNCTION PSIFUNCTU(PZL) RESULT(PSIFCTU)
!#######################################################################################
!
!****  *PSIFUNCTU*
!
!       PURPOSE
!       -------
!       To evaluate the psi stability function which corrects the Monin-Obukhov wind
!       speed profile from the Monin-Obukhov stability parameter zeta = z/L.
!
!       EXTERNAL
!       --------
!
!       IMPLICIT ARGUMENTS
!       ------------------
!
!       REFERENCE
!       ---------
!
!       Original formulations :
!         Liu et al (1979) : Liu, W. T., K. B. Katsaros, and J. A. Businger (1979)
!          Bulk parameterization of air-sea exchanges of heat and water vapor including
!          the molecular constraints at the interface. J. Atm. Sci., 36, 1722--1735.
!         Dyer and Hicks (1970) : Dyer, A. J., and B. B. Hicks (1970) Flux-gradient relationship
!          in the constant flux layer. Quart. J. Roy. Meteor. Soc., 96, 715--721.
!
!       Formulations included in 2021 :
!         Beljaars and Holtslag (1991) : Beljaars, A. C. M. and Holtslag, A. A. M.
!          (1991) Flux Parameterization over Land Surfaces for Atmospheric Models.
!          Journal of Applied Meteorology (1988-2005), 30, 327--341.
!          URL:http://www.jstor.org/stable/26186639
!         Dyer (1974) : Dyer, A. J (1974) A review of flux-profile relationships.
!          Boundary--LayerMeteorology, 7, 363--372.
!         Fairall et al (1996) : Fairall, C. W., Bradley, E. F., Rogers, D. P. ,
!          Edson, J. B. and Young, G. S. (1996) Bulk parameterization of air-sea
!          fluxes for Tropical Ocean-Global Atmosphere Coupled-Ocean Atmosphere
!          Response Experiment. Journal of Geophysical Research: Oceans, 101,
!          3747--3764.URL:https://doi.org/10.1029/95JC03205
!         Grachev et al (2000) : Grachev, A. A., Fairall, C. W. and Bradley, E. F.
!          (2000) Convective Profile Constants Revisited. Boundary-Layer Meteorology,
!          494–515.
!         Grachev et al (2007) : Grachev, A. A., Andreas, E. L., Fairall, C. W.,
!          Guest, P. S. and Persson, P. O. G. (2007) SHEBA flux–profile relationships
!          in the stable atmospheric boundary layer. Boundary-Layer Meteorology, 124,
!          315–333. URL:https://doi.org/10.1007/s10546-007-9177-6
!         Holtslag and de Bruin (1988) : Holtslag, A. A. M. and De Bruin, H. A. R.
!          (1988) Applied Modeling of the Nighttime Surface Energy Balance over Land.
!          Journal of Applied Meteorology, 27, 689--704.
!         Lettau (1979) : Lettau, H. H. (1979) Wind and temperature profile prediction
!          for diabatic surface layers including strong inversion cases. Boundary-Layer
!          Meteorology, 17, 443--464. URL:https://doi.org/10.1007/BF00118610.
!         Paulson (1970) : Paulson, C. A. (1970) The Mathematical Representation of Wind
!          Speed and Temperature Profiles in the Unstable Atmospheric Surface Layer.
!          Journal of Applied Meteorology, 9, 857--861.
!
!
!       AUTHOR
!       ------
!
!         C. Lebeaupin  *Météo-France*
!
!       MODIFICATIONS
!       -------------
!
!         Original     06/2006
!         V. Guemas    04/2021  Different options for psi functions
!                               The original structure in two functions was kept but
!                               it might make sense to merge them to avoid code
!                               redundancies
!-------------------------------------------------------------------------------
IMPLICIT NONE
!
!        0.  Declarations
!        =================
!
!        0.1 Declaration of arguments
!
REAL, DIMENSION(:), INTENT(IN)    :: PZL       ! Monin-Obukhov stability parameter
                                               ! zeta = z/L
REAL, DIMENSION(SIZE(PZL))        :: PSIFCTU   ! Output value of psi function
!
!        0.2 declaration of local variables
!
REAL, DIMENSION(SIZE(PZL)) :: ZY,ZX,ZC,ZPSIC,ZPSIK,ZF
INTEGER :: JJ
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!        1. Initialization
!        ==================
!
IF (LHOOK) CALL DR_HOOK('MODE_COARE30_PSI:PSIFUNCTU',0,ZHOOK_HANDLE)

PSIFCTU(:) = 0.
!
!        2. PSI estimation
!        ==================
!
SELECT CASE (CPSISTAB)

  CASE ('LOGLIN')             ! Log-linear expression for psi as given in
    DO JJ=1,SIZE(PZL)         ! Dyer and Hicks (1970) and Dyer (1974) review
      IF(PZL(JJ)>0.) THEN
        PSIFCTU(JJ) = -XGAMMA*PZL(JJ)
      ENDIF
    ENDDO

  CASE ('LETTAU79')           ! Lettau (1979)
    DO JJ=1,SIZE(PZL)
      IF(PZL(JJ)>0.) THEN
        ZX(JJ) = 1+4.5*PZL(JJ)
        PSIFCTU(JJ) = LOG((ZX(JJ)**(0.5)+1)/2.) + 2*LOG((ZX(JJ)**(0.25)+1)/2) &
                -2*ATAN(ZX(JJ)**(0.25)) +XPI/2 +4./3.*(1-ZX(JJ)**(0.75))
      ENDIF
    ENDDO

  CASE ('HOLTSLAG-BRUIN')     ! Holtslag and de Bruin (1988)
    DO JJ=1,SIZE(PZL)
      IF(PZL(JJ)>0.) THEN
        PSIFCTU(JJ) = -3.75/0.35 -0.7*PZL(JJ) + 3.75/0.35*EXP(-0.35*PZL(JJ)) &
                -0.75*PZL(JJ)*EXP(-0.35*PZL(JJ))
      ENDIF
    ENDDO

  CASE ('BELJAARS-HOLTSLAG')  ! Beljaars and Holtslag (1991)
    DO JJ=1,SIZE(PZL)
      IF(PZL(JJ)>0.) THEN
        PSIFCTU(JJ) = -3.335/0.35 -PZL(JJ) +3.335/0.35*EXP(-0.35*PZL(JJ)) &
                -0.667*PZL(JJ)*EXP(-0.35*PZL(JJ))
      ENDIF
    ENDDO

  CASE ('GRACHEV07')          ! Grachev (2007) fit obtained from SHEBA Arctic
    DO JJ=1,SIZE(PZL)         ! campaign (1997-1998)
      IF(PZL(JJ)>0.) THEN
        ZX(JJ) = (PZL(JJ)+1)**(1./3.)
        PSIFCTU(JJ) = -19.5*(ZX(JJ)-1) +3.25*0.3**(1./3.)*(2*LOG((ZX(JJ)+0.3**(1./3.))&
                /(1+0.3**(1./3.))) - LOG((ZX(JJ)**2-0.3**(1./3.)*ZX(JJ)+0.3**(2./.3))/&
                (1-0.3**(1./3.)+0.3**(2./3.))) +2*SQRT(3.)*(ATAN((2*ZX(JJ)-0.3**(1./3.))/&
                (SQRT(3.)*0.3**(1./3.)))-ATAN((2-0.3**(1./3.))/(SQRT(3.)*0.3**(1./3.)))))
      ENDIF
    ENDDO

  CASE DEFAULT                ! Maintained for backward compatibility
    DO JJ=1,SIZE(PZL)         ! Formulation initially used with COARE3.0
      IF(PZL(JJ)>0.) THEN     ! in SURFEX since 2005
        ZC(JJ)=MIN(50.,0.35*PZL(JJ))
        PSIFCTU(JJ)=-((1.+1.*PZL(JJ))**1. + 0.6667*(PZL(JJ)-14.28)/EXP(ZC(JJ)) + 8.525)
      ENDIF
    ENDDO

END SELECT

SELECT CASE (CPSIUNSTAB)

  CASE ('BUSINGER-DYER')    ! Derived independently by Joost Businger(1966) and Arch
    DO JJ=1,SIZE(PZL)       ! Dyer in the mod-1960 (Businger, 1988) - formulation
      IF(PZL(JJ)<0.) THEN   ! taken in Paulson (1970) review
        ZX(JJ) = (1 -16*PZL(JJ))**(-0.25)
        PSIFCTU(JJ) = 2*LOG((1+ZX(JJ)**(-1))/2) +LOG((1+ZX(JJ)**(-2))/2) &
                -2*ATAN(ZX(JJ)**(-1)) + XPI/2
      ENDIF
    ENDDO

  CASE ('FAIRALL96')        ! Fairall et al (1996)
    DO JJ=1,SIZE(PZL)
      IF(PZL(JJ)<0.) THEN
        ZX(JJ) = (1 -16*PZL(JJ))**(-0.25)          ! Businger-Dyer
        ZPSIK(JJ) = 2*LOG((1+ZX(JJ)**(-1))/2) +LOG((1+ZX(JJ)**(-2))/2) &
                -2*ATAN(ZX(JJ)**(-1)) + XPI/2
        !
        ZY(JJ) = (1 -12.87*PZL(JJ))**(1/3)         ! Convective
        ZPSIC(JJ)= 1.5 * LOG((ZY(JJ)*ZY(JJ)+ZY(JJ)+1.)/3.) &
                 - (3.0**0.5) * atan((2.0*ZY(JJ)+1.0)/(3.0**0.5)) &
                 + 4.0        * atan(1.0)/(3.0**0.5)
        !
        ZF(JJ)   = PZL(JJ) * PZL(JJ) / (1.0 + PZL(JJ)*PZL(JJ))
        PSIFCTU(JJ)=(1.-ZF(JJ)) * ZPSIK(JJ) + ZF(JJ) * ZPSIC(JJ)
      ENDIF
    ENDDO

  CASE ('GRACHEV00')       ! Grachev et al (2000) used in COARE3.0
    DO JJ=1,SIZE(PZL)
      IF(PZL(JJ)<0.) THEN
        ZX(JJ) = (1 -16*PZL(JJ))**(-0.25)          ! Businger-Dyer
        ZPSIK(JJ) = 2*LOG((1+ZX(JJ)**(-1))/2) +LOG((1+ZX(JJ)**(-2))/2) &
                -2*ATAN(ZX(JJ)**(-1)) + XPI/2
        !
        ZY(JJ) = (1 -10.15*PZL(JJ))**(1/3)         ! Convective
        ZPSIC(JJ)= 1.5 * LOG((ZY(JJ)*ZY(JJ)+ZY(JJ)+1.)/3.) &
                 - (3.0**0.5) * atan((2.0*ZY(JJ)+1.0)/(3.0**0.5)) &
                 + 4.0        * atan(1.0)/(3.0**0.5)
        !
        ZF(JJ)   = PZL(JJ) * PZL(JJ) / (1.0 + PZL(JJ)*PZL(JJ))
        PSIFCTU(JJ)=(1.-ZF(JJ)) * ZPSIK(JJ) + ZF(JJ) * ZPSIC(JJ)
      ENDIF
    ENDDO

  CASE DEFAULT                ! Maintained for backward compatibility
    DO JJ=1,SIZE(PZL)         ! Formulation initially used with COARE3.0
      IF(PZL(JJ)<0.) THEN     ! in SURFEX since 2005
        ZX(JJ)   = (1.0 - 15. * PZL(JJ))**0.25         ! Kansas unstable
        ZPSIK(JJ)= 2.0 * LOG((1.0+ZX(JJ)       )/2.0) &
                 +       LOG((1.0+ZX(JJ)*ZX(JJ))/2.0) &
                 - 2.0 * atan(ZX(JJ)) &
                 + 2.0 * atan(1.0)
        !
        ZY(JJ)   = (1.0 - 10.15 * PZL(JJ))**0.3333     ! Convective
        ZPSIC(JJ)= 1.5 * LOG((ZY(JJ)*ZY(JJ)+ZY(JJ)+1.)/3.) &
                 - (3.0**0.5) * atan((2.0*ZY(JJ)+1.0)/(3.0**0.5)) &
                 + 4.0        * atan(1.0)/(3.0**0.5)
        !
        ZF(JJ)   =PZL(JJ) * PZL(JJ) / (1.0+PZL(JJ)*PZL(JJ))
        !
        PSIFCTU(JJ)=(1.-ZF(JJ)) * ZPSIK(JJ) + ZF(JJ) * ZPSIC(JJ)
      ENDIF
    ENDDO

END SELECT

IF (LHOOK) CALL DR_HOOK('MODE_COARE30_PSI:PSIFUNCTU',1,ZHOOK_HANDLE)

END FUNCTION PSIFUNCTU
!---------------------------------------------------------------------------------------
!
!#######################################################################################
FUNCTION PSIFUNCTT(PZL) RESULT(PSIFCTT)
!#######################################################################################
!
!****  *PSIFUNCTT*
!
!       PURPOSE
!       -------
!       To evaluate the psi stability function which corrects the Monin-Obukhov temperature
!       and humidity profiles from the Monin-Obukhov stability parameter zeta = z/L.
!
!       EXTERNAL
!       --------
!
!       IMPLICIT ARGUMENTS
!       ------------------
!
!       REFERENCE
!       ---------
!
!       Original formulations :
!         Liu et al (1979) : Liu, W. T., K. B. Katsaros, and J. A. Businger (1979)
!          Bulk parameterization of air-sea exchanges of heat and water vapor including
!          the molecular constraints at the interface. J. Atm. Sci., 36, 1722--1735.
!         Dyer and Hicks (1970) : Dyer, A. J., and B. B. Hicks (1970) Flux-gradient relationship
!          in the constant flux layer. Quart. J. Roy. Meteor. Soc., 96, 715--721.
!
!       Formulations included in 2021 :
!         Beljaars and Holtslag (1991) : Beljaars, A. C. M. and Holtslag, A. A. M.
!          (1991) Flux Parameterization over Land Surfaces for Atmospheric Models.
!          Journal of Applied Meteorology (1988-2005), 30, 327--341.
!          URL:http://www.jstor.org/stable/26186639
!         Dyer (1974) : Dyer, A. J (1974) A review of flux-profile relationships.
!          Boundary--LayerMeteorology, 7, 363--372.
!         Fairall et al (1996) : Fairall, C. W., Bradley, E. F., Rogers, D. P. ,
!          Edson, J. B. and Young, G. S. (1996) Bulk parameterization of air-sea
!          fluxes for Tropical Ocean-Global Atmosphere Coupled-Ocean Atmosphere
!          Response Experiment. Journal of Geophysical Research: Oceans, 101,
!          3747--3764.URL:https://doi.org/10.1029/95JC03205
!         Grachev et al (2000) : Grachev, A. A., Fairall, C. W. and Bradley, E. F.
!          (2000) Convective Profile Constants Revisited. Boundary-Layer Meteorology,
!          494–515.
!         Grachev et al (2007) : Grachev, A. A., Andreas, E. L., Fairall, C. W.,
!          Guest, P. S. and Persson, P. O. G. (2007) SHEBA flux–profile relationships
!          in the stable atmospheric boundary layer. Boundary-Layer Meteorology, 124,
!          315–333. URL:https://doi.org/10.1007/s10546-007-9177-6
!         Holtslag and de Bruin (1988) : Holtslag, A. A. M. and De Bruin, H. A. R.
!          (1988) Applied Modeling of the Nighttime Surface Energy Balance over Land.
!          Journal of Applied Meteorology, 27, 689--704.
!         Lettau (1979) : Lettau, H. H. (1979) Wind and temperature profile prediction
!          for diabatic surface layers including strong inversion cases. Boundary-Layer
!          Meteorology, 17, 443--464. URL:https://doi.org/10.1007/BF00118610.
!         Paulson (1970) : Paulson, C. A. (1970) The Mathematical Representation of Wind
!          Speed and Temperature Profiles in the Unstable Atmospheric Surface Layer.
!          Journal of Applied Meteorology, 9, 857--861.
!
!
!       AUTHOR
!       ------
!
!         C. Lebeaupin  *Météo-France*
!
!       MODIFICATIONS
!       -------------
!
!         Original     06/2006
!         V. Guemas    04/2021  Different options for psi functions
!                               The original structure in two functions was kept but
!                               it might make sense to merge them to avoid code
!-------------------------------------------------------------------------------
IMPLICIT NONE
!
!        0.  Declarations
!        =================
!
!        0.1 Declaration of arguments
!
REAL, DIMENSION(:), INTENT(IN)    :: PZL       ! Monin-Obukhov stability parameter
                                               ! zeta = z/L
REAL, DIMENSION(SIZE(PZL))        :: PSIFCTT   ! Output value of psi function
!
!        0.2 declaration of local variables
!
REAL, DIMENSION(SIZE(PZL)) :: ZY,ZX,ZC,ZPSIC,ZPSIK,ZF
INTEGER :: JJ
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!        1. Initialization
!        ==================
!
IF (LHOOK) CALL DR_HOOK('MODE_COARE30_PSI:PSIFUNCTT',0,ZHOOK_HANDLE)
!
PSIFCTT(:) = 0.
!
!        2. PSI estimation
!        ==================
!
SELECT CASE (CPSISTAB)

  CASE ('LOGLIN')             ! Log-linear expression for psi as given in
    DO JJ=1,SIZE(PZL)         ! Dyer and Hicks (1970) and Dyer (1974) review
      IF(PZL(JJ)>0.) THEN
        PSIFCTT(JJ) = -XGAMMA*PZL(JJ)
      ENDIF
    ENDDO

  CASE ('LETTAU79')           ! Lettau (1979)
    DO JJ=1,SIZE(PZL)
      IF(PZL(JJ)>0.) THEN
        ZX(JJ) = 1+4.5*PZL(JJ)
        PSIFCTT(JJ) = 2*LOG((ZX(JJ)**(0.25)+1)/2.) -2*(ZX(JJ)**(1.5)/3. &
                +ZX(JJ)**0.5 -4./3.)
      ENDIF
    ENDDO

  CASE ('HOLTSLAG-BRUIN')     ! Holtslag and de Bruin (1988)
    DO JJ=1,SIZE(PZL)
      IF(PZL(JJ)>0.) THEN
        PSIFCTT(JJ) = -3.75/0.35 -0.7*PZL(JJ) + 3.75/0.35*EXP(-0.35*PZL(JJ)) &
                -0.75*PZL(JJ)*EXP(-0.35*PZL(JJ))
      ENDIF
    ENDDO

  CASE ('BELJAARS-HOLTSLAG')  ! Beljaars and Holtslag (1991)
    DO JJ=1,SIZE(PZL)         ! used in COARE3.0
      IF(PZL(JJ)>0.) THEN
        PSIFCTT(JJ) = -2.985/0.35 -(1+2./3.*PZL(JJ))**(3./2.) &
                +3.335/0.35*EXP(-0.35*PZL(JJ)) -0.667*PZL(JJ)*EXP(-0.35*PZL(JJ))
      ENDIF
    ENDDO

  CASE ('GRACHEV07')          ! Grachev (2007) fit obtained from SHEBA Arctic
    DO JJ=1,SIZE(PZL)         ! campaign (1997-1998)
      IF(PZL(JJ)>0.) THEN
        ZX(JJ) = (PZL(JJ)+1)**(1./3.)
        PSIFCTT(JJ) = -2.5*LOG(1+3*PZL(JJ)+PZL(JJ)**2) + 5/(2*SQRT(5.))*&
                (LOG((2*PZL(JJ)+3-SQRT(5.))/(2*PZL(JJ)+3+SQRT(5.))) &
                -LOG((3-SQRT(5.))/(3+SQRT(5.))) )
      ENDIF
    ENDDO

  CASE DEFAULT                ! Maintained for backward compatibility
    DO JJ=1,SIZE(PZL)         ! Formulation initially used with COARE3.0
      IF(PZL(JJ)>0.) THEN     ! in SURFEX since 2005
        ZC(JJ)=MIN(50.,0.35*PZL(JJ))
        PSIFCTT(JJ)=-((1.+2.*PZL(JJ)/3.)**1.5 + 0.6667*(PZL(JJ)-14.28)/EXP(ZC(JJ)) + 8.525)
      ENDIF
    ENDDO

END SELECT

SELECT CASE (CPSIUNSTAB)

  CASE ('BUSINGER-DYER')    ! Derived independently by Joost Businger(1966) and Arch
    DO JJ=1,SIZE(PZL)       ! Dyer in the mod-1960 (Businger, 1988) - formulation
      IF(PZL(JJ)<0.) THEN   ! taken in Paulson (1970) review
        ZX(JJ) = (1 -16*PZL(JJ))**(-0.5)
        PSIFCTT(JJ) = 2*LOG((1+ZX(JJ)**(-1))/2)
      ENDIF
    ENDDO

  CASE ('FAIRALL96')        ! Fairall et al (1996)
    DO JJ=1,SIZE(PZL)
      IF(PZL(JJ)<0.) THEN
        ZX(JJ) = (1 -16*PZL(JJ))**(-0.5)           ! Businger-Dyer
        ZPSIK(JJ) = 2*LOG((1+ZX(JJ)**(-1))/2)
        !
        ZY(JJ) = (1 -12.87*PZL(JJ))**(1/3)         ! Convective
        ZPSIC(JJ)= 1.5 * LOG((ZY(JJ)*ZY(JJ)+ZY(JJ)+1.)/3.) &
                 - (3.0**0.5) * atan((2.0*ZY(JJ)+1.0)/(3.0**0.5)) &
                 + 4.0        * atan(1.0)/(3.0**0.5)
        !
        ZF(JJ)   = PZL(JJ) * PZL(JJ) / (1.0 + PZL(JJ)*PZL(JJ))
        PSIFCTT(JJ)=(1.-ZF(JJ)) * ZPSIK(JJ) + ZF(JJ) * ZPSIC(JJ)
      ENDIF
    ENDDO

  CASE ('GRACHEV00')       ! Grachev et al (2000) used in COARE3.0
    DO JJ=1,SIZE(PZL)
      IF(PZL(JJ)<0.) THEN
        ZX(JJ) = (1 -16*PZL(JJ))**(-0.5)           ! Businger-Dyer
        ZPSIK(JJ) = 2*LOG((1+ZX(JJ)**(-1))/2)
        !
        ZY(JJ) = (1 -34.15*PZL(JJ))**(1/3)         ! Convective
        ZPSIC(JJ)= 1.5 * LOG((ZY(JJ)*ZY(JJ)+ZY(JJ)+1.)/3.) &
                 - (3.0**0.5) * atan((2.0*ZY(JJ)+1.0)/(3.0**0.5)) &
                 + 4.0        * atan(1.0)/(3.0**0.5)
        !
        ZF(JJ)   = PZL(JJ) * PZL(JJ) / (1.0 + PZL(JJ)*PZL(JJ))
        PSIFCTT(JJ)=(1.-ZF(JJ)) * ZPSIK(JJ) + ZF(JJ) * ZPSIC(JJ)
      ENDIF
    ENDDO

  CASE DEFAULT                ! Maintained for backward compatibility
    DO JJ=1,SIZE(PZL)         ! Formulation initially used with COARE3.0
      IF(PZL(JJ)<0.) THEN     ! in SURFEX since 2005
        ZX(JJ)   = (1. - 15. * PZL(JJ))**.5         ! Kansas unstable
        ZPSIK(JJ)= 2.0 * LOG((1.0+ZX(JJ)       )/2.0)
        !
        ZY(JJ)   = (1.0 - 34.15 * PZL(JJ))**0.3333  ! Convective
        ZPSIC(JJ)= 1.5 * LOG((ZY(JJ)*ZY(JJ)+ZY(JJ)+1.0)/3.) &
                 - (3.0**0.5) * atan((2.0*ZY(JJ)+1.0)/(3.0**0.5)) &
                 + 4.0        * atan(1.0)/(3.0**0.5)
        !
        ZF(JJ)   = PZL(JJ) * PZL(JJ) / (1.0+PZL(JJ)*PZL(JJ))
        !
        PSIFCTT(JJ)= (1.-ZF(JJ)) * ZPSIK(JJ) + ZF(JJ) * ZPSIC(JJ)
      ENDIF
    ENDDO

END SELECT

IF (LHOOK) CALL DR_HOOK('MODE_COARE30_PSI:PSIFUNCTT',1,ZHOOK_HANDLE)

END FUNCTION PSIFUNCTT
!
END MODULE MODE_COARE30_PSI
