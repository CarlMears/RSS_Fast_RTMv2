module goff_gratch 

    contains
        SUBROUTINE goff_gratch_vap(T,RH,P,  P_V,RHO_V)
        IMPLICIT NONE
        !
        !    
        !     CALCULATES WATER VAPOR PRESSURE P_V AND DENSITY RHO_V
        !     FROM RELATIVE HUMIDITY RH
        !     S. CRUZ POL, C. RUF AND S. KEIHM, RADIO SCIENCE 33 (5),1319 (1998)
        ! 
        !    WATER VAPOR PRESSURE: P_V = P_S * RH
        !    WATER VAPOR DENSITY   RHO_V = F_W *(P_V * EPS) / (R_D * T)
        !
            REAL(4) :: T,RH,P,P_V,RHO_V,P_S
            REAL(4) :: F_W,XI1,XI2,XI3  
            REAL(4)  :: P_STANDARD = 1013.246  ! HPA
            REAL(4)  :: TS = 373.14     ! K  (WATER BOILING)
            REAL(4)  :: T0 = 273.14  ! K 
            REAL(4)  :: EPS = 1./1.607795  ! M(H2O)/M(AIR)
            REAL(4)  :: R_D = 287.05 ! J /(K*KG)  GAS CONSTANT FOR DRY AIR
        !
            F_W = 1.0 + 1.E-4 * (5.92854 + 3.740346E-2 * P + 1.971198E-4 * (T-T0) * (800-P)  +  6.045511E-6 * P * (T-T0)**2 ) 
            ! DEVIATION FROM IDEAL GAS
        !
            XI1 = -7.90298*(TS/T - 1) + 5.02808*LOG10(TS/T)
            XI2 = -1.3816E-7 * 10**(11.344*(1.-T/TS) -1 )
            XI3 = 8.1328E-3 * (10**(-3.49149*(TS/T - 1)) - 1)
        !
            P_S = P_STANDARD * 10**(XI1 + XI2 + XI3)
            P_V = P_S * RH* 0.01                            !  MBAR
            RHO_V = ((F_W * P_V * EPS) / (R_D * T)) * 1.E5  !  [G/M**3]
        !
        RETURN
        END SUBROUTINE goff_gratch_vap 
end module goff_gratch