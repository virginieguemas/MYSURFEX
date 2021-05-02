!depfile:form_drag.F90
MODULE MODI_FORM_DRAG
INTERFACE
     SUBROUTINE FORM_DRAG (PSIC, PZ0SEA, PZ0ICE, PUREF, PZREF, PLMOI, PLMOO, &
                           PRATIOI, PRATIOO, PCDF, PCHF)
REAL, DIMENSION(:), INTENT(IN)   :: PSIC    ! sea ice concentration
REAL, DIMENSION(:), INTENT(IN)   :: PZ0SEA  ! skin drag coefficient over water
REAL, DIMENSION(:), INTENT(IN)   :: PZ0ICE  ! skin drag coefficient over ice
REAL, DIMENSION(:), INTENT(IN)   :: PUREF   ! atmospheric height for wind speed
REAL, DIMENSION(:), INTENT(IN)   :: PZREF   ! atmospheric height for temperature and moisture
REAL, DIMENSION(:), INTENT(IN)   :: PLMOI   ! Monin-Obukhov length above ice
REAL, DIMENSION(:), INTENT(IN)   :: PLMOO   ! Monin-Obukhov length above ocean
REAL, DIMENSION(:), INTENT(IN)   :: PRATIOI ! ratio between scalar and aerodynamic roughness
REAL, DIMENSION(:), INTENT(IN)   :: PRATIOO ! ratio between scalar and aerodynamic roughness
REAL, DIMENSION(:), INTENT(OUT)  :: PCDF    ! output form drag coefficient
REAL, DIMENSION(:), INTENT(OUT)  :: PCHF    ! Equivalent form drag coefficient for heat
END SUBROUTINE FORM_DRAG

END INTERFACE
END MODULE MODI_FORM_DRAG

