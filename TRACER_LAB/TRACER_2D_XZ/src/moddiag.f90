!################   CALCUL DES RHS   #################!
!#                                                   #!
!# auteur : Voitus F.                                #!
!# sujet  : Discrétisation temporelle                #!
!#                                                   #!
!#####################################################!

!=====================================================!
!====================  MODULE DIAG ===================!
!=====================================================!

MODULE MOD_DIAG

  USE MOD_SHARE,     ONLY : NLEV,NXPT,RG,RR,RCP,RPSUR,RP00, &
                          & ZERO,ONE,TWO,HALF,RDX,RDZ,RDZTOP
  
CONTAINS

  !*******************************************************************************!
  !*******************************************************************************!
  !*******************************************************************************!
  !*****************************  LISTE DES ROUTINES *****************************!
  !*******************************************************************************!
  !*******************************************************************************!
  !*******************************************************************************!
  SUBROUTINE DIAG_PIS_MDOT(PM_DOT,PPIS_DOT,PPIS,PM,PDIV, &
                         & PETADOTH,PDELA,PDELB,PRDETA)
    
    IMPLICIT NONE

    REAL(8), DIMENSION(NLEV),   INTENT(OUT) :: PM_DOT
    REAL(8),                    INTENT(OUT) :: PPIS_DOT
    REAL(8), DIMENSION(NLEV),   INTENT(IN)  :: PM,PDIV
    REAL(8), DIMENSION(NLEV),   INTENT(IN)  :: PDELA,PDELB
    REAL(8), DIMENSION(0:NLEV), INTENT(IN)  :: PETADOTH
    REAL(8), DIMENSION(NLEV),   INTENT(IN)  :: PRDETA
    REAL(8),                    INTENT(IN)  :: PPIS
    
    INTEGER(8)                              :: J
    REAL(8), DIMENSION(NLEV)                :: ZDELPI

    DO J=1,NLEV
       ZDELPI(J) = PDELA(J) + PDELB(J)*PPIS
    END DO
    
    CALL OP_N(PPIS_DOT,-PDIV,ZDELPI)

    DO J=1,NLEV
       PM_DOT(J) = -PM(J)*( PDIV(J)  +   &
                 & PRDETA(J)*(PETADOTH(J)-PETADOTH(J-1)) )
    END DO
    
  END SUBROUTINE DIAG_PIS_MDOT
  !*******************************************************************************!
  !*******************************************************************************!
  !*******************************************************************************!
  !*******************************************************************************!
  !*******************************************************************************!
  SUBROUTINE DIAG_EULERIAN_ORO(PADV_PIS_REF,PU,PDPISDX_REF,PDELB)
  
    IMPLICIT NONE
  
    REAL(8),                   INTENT(OUT)   :: PADV_PIS_REF
    REAL(8), DIMENSION(NLEV),  INTENT(IN)    :: PU,PDELB
    REAL(8),                   INTENT(IN)    :: PDPISDX_REF
  
    INTEGER                                  :: J
    REAL(8), DIMENSION(NLEV)                 :: ZADV_PP

    DO J=1,NLEV
       ZADV_PP(J) = -PU(J)*PDPISDX_REF
    END DO
  
    CALL OP_N(PADV_PIS_REF,ZADV_PP,PDELB)
    
  END SUBROUTINE DIAG_EULERIAN_ORO
  !***********************************************************************************!
  !***********************************************************************************!
  !***********************************************************************************!
  !***********************************************************************************!
  SUBROUTINE DEF_ETADOT(PETADOT,PETADOTH,PU,PDIV,PPIS,PDPISDX, &
                      & PDELA,PDELB,PRDETA,PAH,PBH,PRDETAH)

    IMPLICIT NONE

    REAL(8), DIMENSION(NLEV),   INTENT(OUT) :: PETADOT
    REAL(8), DIMENSION(0:NLEV), INTENT(OUT) :: PETADOTH
    REAL(8), DIMENSION(NLEV),   INTENT(IN)  :: PU,PDIV
    REAL(8), DIMENSION(NLEV),   INTENT(IN)  :: PDELA,PDELB,PRDETA
    REAL(8), DIMENSION(0:NLEV), INTENT(IN)  :: PAH,PBH,PRDETAH
    REAL(8),                    INTENT(IN)  :: PPIS,PDPISDX
    
    INTEGER(8)                              :: J
    REAL(8), DIMENSION(NLEV)                :: ZPIF
    REAL(8), DIMENSION(0:NLEV)              :: ZMETADOT,ZSDIV_MV

    !* vertical momentum flux at interfaces
    !------------------------------------------------------------!
    ZSDIV_MV(0)    = ZERO 
    DO J = 1, NLEV
       ZSDIV_MV(J) = ZSDIV_MV(J-1)                       &
                   & + ((PDELA(J)+PDELB(J)*PPIS)*PDIV(J) &
                   & + PDELB(J)*PU(J)*PDPISDX)
    END DO
    DO J = 1, NLEV-1
       ZMETADOT(J) = PBH(J)*ZSDIV_MV(NLEV)-ZSDIV_MV(J)
    END DO
    
    ! Material boundary conditions
    ZMETADOT(0)    = ZERO
    ZMETADOT(NLEV) = ZERO
    
    !------------------------------------------------------------!
    DO J=1,NLEV
       PETADOT(J)  = ((ZMETADOT(J)+ZMETADOT(J-1))/TWO)   &
                   & /(PRDETA(J)*(PDELA(J)+PDELB(J)*PPIS))
    END DO

    !* vertical velocity at interfaces
    CALL DEF_P(ZPIF,PPIS,PAH,PBH,1)

    DO J = 1,NLEV-1
       PETADOTH(J) = ZMETADOT(J)/(PRDETAH(J)*(ZPIF(J+1)-ZPIF(J)))
    END DO   
    PETADOTH(0)    = ZERO
    PETADOTH(NLEV) = ZERO
        
  END SUBROUTINE DEF_ETADOT
  !*******************************************************************************!
  !*******************************************************************************!
  !*******************************************************************************!
  !*******************************************************************************!
   SUBROUTINE DEF_M(PM,PPIS,PDELA,PDELB,PRDETA)

    IMPLICIT NONE

    REAL(8), DIMENSION(NLEV),   INTENT(OUT) :: PM
    REAL(8), DIMENSION(NLEV),   INTENT(IN)  :: PRDETA
    REAL(8), DIMENSION(NLEV),   INTENT(IN)  :: PDELA,PDELB
    REAL(8),                    INTENT(IN)  :: PPIS
    
    INTEGER(8)                              :: J

    DO J=1,NLEV
       PM(J) = (PDELA(J)+PDELB(J)*PPIS)*PRDETA(J)
    END DO

  END SUBROUTINE DEF_M  
  !*******************************************************************************!
  !*******************************************************************************!
  !*******************************************************************************!
  !*******************************************************************************!
  SUBROUTINE DEF_Z(KTOP,PZ,PT,PPIS,PAH,PBH,PZS)

   IMPLICIT NONE

   INTEGER,                       INTENT(IN)   :: KTOP
   REAL(8), DIMENSION(KTOP:NLEV), INTENT(OUT)  :: PZ
   REAL(8), DIMENSION(NLEV),      INTENT(IN)   :: PT
   REAL(8), DIMENSION(0:NLEV),    INTENT(IN)   :: PAH,PBH
   REAL(8),                       INTENT(IN)   :: PPIS,PZS
      
   INTEGER(8)                                  :: J
   REAL(8)                                     :: ZFULL,ZZ
   REAL(8), DIMENSION(NLEV)                    :: ZALFA,ZDELTA
   
   ZFULL= ONE
   IF (KTOP == 0) ZFULL = ZERO

   ZALFA(1)       = ONE
   ZDELTA(1)      = ONE+(RCP/RR)
   DO J=2,NLEV
      ZDELTA(J)   = ((PAH(J)-PAH(J-1))+PPIS*(PBH(J)-PBH(J-1)))/&
                  & DSQRT((PAH(J-1)+PPIS*PBH(J-1))*(PAH(J)+PPIS*PBH(J)))
      ZALFA(J)    = ONE-(DSQRT((PAH(J-1)+PPIS*PBH(J-1))*(PAH(J)+PPIS*PBH(J)))/&
                  & (PAH(J)+PPIS*PBH(J)))
   END DO
   
   ZZ             = PZS
   PZ(NLEV)       = ZZ+ZFULL*(RR/RG)*PT(NLEV)*ZALFA(NLEV)     
   DO J = 0,NLEV-2
      ZZ          = ZZ+(RR/RG)*PT(NLEV-J)*ZDELTA(NLEV-J)
      PZ(NLEV-J-1)= ZZ+ZFULL*(RR/RG)*PT(NLEV-J-1)*ZALFA(NLEV-J-1)
   END DO
   ZZ             = ZZ + (RR/RG)*PT(1)*ZDELTA(1) 

   IF (KTOP == 0) PZ(KTOP) = ZZ      

  END SUBROUTINE DEF_Z
  !************************************************************************************!
  !************************************************************************************!
  !************************************************************************************!
  SUBROUTINE DEF_P(PPI,PPIS,PAH,PBH,KTOP)

    IMPLICIT NONE

    REAL(8), DIMENSION(KTOP:NLEV), INTENT(OUT)  :: PPI
    REAL(8), DIMENSION(0:NLEV),    INTENT(IN)   :: PAH,PBH
    REAL(8),                       INTENT(IN)   :: PPIS
    INTEGER,                       INTENT(IN)   :: KTOP
    
    INTEGER                                     :: J

    IF (KTOP == 0) THEN
       DO J=0,NLEV
          PPI(J)=PAH(J)+PBH(J)*PPIS
       END DO   
    ELSE IF (KTOP == 1) THEN   
       PPI(1) = ((PAH(1)-PAH(0))+PPIS*(PBH(1)-PBH(0))) &
              & /(ONE+(RCP/RR))
       DO J=2,NLEV
          PPI(J) = DSQRT((PAH(J)+PPIS*PBH(J))*(PAH(J-1)&
                 & +PPIS*PBH(J-1)))   
       END DO          
    END IF
    
  END SUBROUTINE DEF_P
  !*************************************************************************************!
  !*************************************************************************************!
  !*************************************************************************************!
  SUBROUTINE DEF_DZ(PDZ,PDZMIN,PDZMAX,PZLIM,KILEV,LREGULAR)
    
    IMPLICIT NONE

    REAL(8), DIMENSION(NLEV),      INTENT(OUT) :: PDZ
    REAL(8),                       INTENT(IN)  :: PDZMIN,PDZMAX,PZLIM
    INTEGER(8),                    INTENT(IN)  :: KILEV
    LOGICAL,                       INTENT(IN)  :: LREGULAR
    
    INTEGER                                    :: L
    REAL(8)                                    :: ZRI,ZR0,ZL                   
    
    IF (LREGULAR) THEN
       PDZ(2:NLEV) = RDZ
       PDZ(1)      = RDZTOP
    ELSE
       ZRI = PZLIM*PDZMAX+(ONE-PZLIM)*PDZMIN
       ZR0 = (REAL(NLEV-2,8)*ZRI-REAL(KILEV-1,8)*PDZMAX) &
           & /REAL(NLEV-1-KILEV,8)
       
       DO L=1,KILEV
          ZL            = REAL(L-1,8)/REAL(KILEV-1,8)
          PDZ(NLEV-L+1) = ZR0 + ZL*(ZRI-ZR0)             &
                        & + ((ONE-ZL)**2)*(PDZMIN-ZR0)
       END DO
       DO L=KILEV+1,NLEV-1
          ZL            = REAL(L-1,8)/REAL(NLEV-2,8)
          PDZ(NLEV-L+1) = ZR0 + ZL*(PDZMAX-ZR0) 
       END DO
       PDZ(1) = RDZTOP
    END IF   
    
  END SUBROUTINE DEF_DZ
  !*******************************************************************************!
  !*******************************************************************************!
  !*******************************************************************************!
  !*******************************************************************************!
  !*******************************************************************************!
  SUBROUTINE OP_N(POUT,PIN,PDELPI)

    IMPLICIT NONE

    REAL(8),                  INTENT(OUT) :: POUT
    REAL(8), DIMENSION(NLEV), INTENT(IN)  :: PIN,PDELPI
    INTEGER                               :: J

    POUT = PDELPI(1)*PIN(1)
    DO J = 2, NLEV
       POUT = POUT + PDELPI(J)*PIN(J)
    END DO

  END SUBROUTINE OP_N
  !*******************************************************************************!
  !*******************************************************************************!
  !*******************************************************************************!
  !*******************************************************************************!
  SUBROUTINE GP_NORMS(PX,PNORMX,CNORM_TYPE)
   
   IMPLICIT NONE

   REAL(8), DIMENSION(NLEV,NXPT),  INTENT(IN)    :: PX   
   REAL(8),                        INTENT(OUT)   :: PNORMX
   CHARACTER(LEN=2),               INTENT(IN)    :: CNORM_TYPE
   INTEGER(8)                                    :: I,J
   REAL(8)                                       :: ZNORM, ZGPTOT
   
   ZGPTOT = REAL(NXPT,8)*REAL(NLEV,8)

   IF (CNORM_TYPE == 'L2') THEN
     ZNORM = ZERO
     DO I=1,NXPT
        DO J =1,NLEV
        ZNORM = ZNORM + PX(J,I)**2
        END DO
     END DO   
     PNORMX = DSQRT(ZNORM/ZGPTOT)
   ELSE IF (CNORM_TYPE == 'L1') THEN
     ZNORM = ZERO
     DO I=1,NXPT
        DO J =1,NLEV
        ZNORM = ZNORM + ABS(PX(J,I))
        END DO
     END DO
     PNORMX = ZNORM/ZGPTOT
   ELSE 
     PNORMX = MAXVAL(ABS(PX))
   END IF

  END SUBROUTINE GP_NORMS
 !**************************************************************************************!
 !**************************************************************************************!
 !**************************************************************************************!
 
END MODULE MOD_DIAG

!=====================================================!
!=====================================================!
!=====================================================!
