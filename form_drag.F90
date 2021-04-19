!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
     SUBROUTINE FORM_DRAG ( )
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
IMPLICIT NONE


END SUBROUTINE FORM_DRAG
