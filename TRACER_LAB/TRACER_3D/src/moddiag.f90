!################   CALCUL DES RHS   #################!
!#                                                   #!
!# auteur : C.Colavolpe & Voitus F.                  #!
!# sujet  : Discrétisation temporelle                #!
!#                                                   #!
!#####################################################!

!=====================================================!
!====================  MODULE DIAG ===================!
!=====================================================!

MODULE MOD_DIAGNOSTIC

  USE MOD_PARAMETER, ONLY : NLEV,GXPT,GYPT,RKAPPA,RG,RR,RCP,RCV,RPSUR,RP00,&
                          & ZERO,ONE,TWO,HALF,RSLAG,NVERINT
  
CONTAINS

 !*****************************************************!
 !***************  LISTE DES ROUTINES  ****************!
 !*****************************************************!
  
 !*******************************************************************************!
 !*******************************************************************************!
 SUBROUTINE DIAG_UF2UH(PUH,PVH,PUF,PVF,PEPS)

    IMPLICIT NONE

    REAL(8), DIMENSION(0:NLEV), INTENT(OUT) :: PUH,PVH
    REAL(8), DIMENSION(NLEV),   INTENT(IN)  :: PUF,PVF
    REAL(8), DIMENSION(0:NLEV), INTENT(IN)  :: PEPS

    INTEGER(8)                              :: J

    PUH(0) = (ONE-PEPS(0))*PUF(1)
    PVH(0) = (ONE-PEPS(0))*PVF(1) 
    DO J = 1, NLEV-1
       PUH(J) = (ONE-PEPS(J))*PUF(J+1) + PEPS(J)*PUF(J)
       PVH(J) = (ONE-PEPS(J))*PVF(J+1) + PEPS(J)*PVF(J)
    ENDDO
    PUH(NLEV) = PEPS(NLEV)*PUF(NLEV)
    PVH(NLEV) = PEPS(NLEV)*PVF(NLEV)

 END SUBROUTINE DIAG_UF2UH
 !*******************************************************************************!
 !*******************************************************************************!
 !*******************************************************************************!
 !*******************************************************************************!
 SUBROUTINE DIAG_MASS_METRIC(PMF,PMH,PDELTAPIF,PIF,PDETAF,PIH,PDETAH)

    IMPLICIT NONE

    REAL(8), DIMENSION(NLEV),      INTENT(OUT) :: PMF
    REAL(8), DIMENSION(0:NLEV),    INTENT(OUT) :: PMH
    REAL(8), DIMENSION(NLEV),      INTENT(IN)  :: PIF,PDELTAPIF,PDETAF
    REAL(8), DIMENSION(0:NLEV),    INTENT(IN)  :: PIH,PDETAH

    INTEGER(8)                                 :: J

    DO J=1,NLEV
       PMF(J) = PDELTAPIF(J)/PDETAF(J)
    END DO
    
    PMH(0)    = (PIF(1)-PIH(0))/PDETAH(0)
    PMH(NLEV) = (PIH(NLEV)-PIF(NLEV))/PDETAH(NLEV)
    
    DO J=1,NLEV-1
       PMH(J) = (PIF(J+1)-PIF(J))/PDETAH(J)
    END DO   

 END SUBROUTINE DIAG_MASS_METRIC
 !*******************************************************************************!
 !*******************************************************************************!
 !*******************************************************************************!
 SUBROUTINE DIAG_VERTICAL_METRICS(PPIH,PPIF,PDELTP,PDELTA,PALPHA,PAH,PBH,PPIS)

    IMPLICIT NONE

    REAL(8), DIMENSION(NLEV),   INTENT(OUT)  :: PPIF,PDELTP,PALPHA
    REAL(8), DIMENSION(NLEV),   INTENT(OUT)  :: PDELTA
    REAL(8), DIMENSION(0:NLEV), INTENT(OUT)  :: PPIH
    REAL(8), DIMENSION(0:NLEV), INTENT(IN)   :: PAH,PBH
    REAL(8)                   , INTENT(IN)   :: PPIS
   
    INTEGER                                  :: J
    
    PPIH(0) = PAH(0)+PBH(0)*ZPIS
    
    DO J=1,NLEV
       PPIH(J)    = PAH(J)+PBH(J)*ZPIS
       PDELTP(J)  = (PAH(J)-PAH(J-1))+(PBH(J)-PBH(J-1))*PPIS
    END DO

    PALPHA(1)    = ONE
    PDELTA(1)    = ONE + ONE/RKAPPA
    PPIF(1)      = PDELTP(1)/PDELTA(1)
    
    DO J = 2,NLEV
       PPIF(J)   = DSQRT(PPIH(J)*PPIH(J-1))
       PALPHA(J) = ONE-(PPIF(J)/PPIH(J))
       PDELTA(J) = PDELTP(J)/PPIF(J)
    END DO
   

 END SUBROUTINE DIAG_VERTICAL_METRICS
 !*******************************************************************************!
 !*******************************************************************************!
 !*******************************************************************************!
 !*******************************************************************************!
 SUBROUTINE DIAG_PISDOT(PN_MDIV,PU,PV,PDIV,PPIS,PPIS_REF, &
             & PDPISDX_REF,PDPISDY_REF,PDELTB,PAH)
    
    IMPLICIT NONE

    REAL(8),                    INTENT(OUT) :: PN_MDIV
    REAL(8),                    INTENT(IN)  :: PPIS,PPIS_REF
    REAL(8),                    INTENT(IN)  :: PDPISDX_REF,PDPISDY_REF
    REAL(8), DIMENSION(NLEV),   INTENT(IN)  :: PU,PV,PDIV
    REAL(8), DIMENSION(0:NLEV), INTENT(IN)  :: PAH
    REAL(8), DIMENSION(NLEV),   INTENT(IN)  :: PDELTB
    
    INTEGER(8)                              :: J
    REAL(8)                                 :: ZPIS_ADV
    REAL(8), DIMENSION(NLEV)                :: ZDELTAPI

    ZPIS_ADV = ZERO
    DO J = 1, NLEV
       ZPIS_ADV    = ZPIS_ADV          &
                   & + PDELTB(J)*(PU(J)*PDPISDX+PV(J)*PDPISDY)
       ZDELTAPI(J) = (PAH(J)-PAH(J-1)) &
                   & + PDELTB(J)*(RSLAG*PPIS+(ONE-RSLAG)*PPIS_REF)
    END DO

    CALL OP_N(PN_MDIV,PDIV,ONE,ZDELTAPI)

    PN_MDIV = PN_MDIV + ZPIS_ADV

 END SUBROUTINE DIAG_PIDOT
 !************************************************************************************!
 !************************************************************************************!
 !************************************************************************************!
 SUBROUTINE DIAG_METADOT(PMETADOT,PU,PV,PDIV,PDPISDX,PDPISDY,PDELTAPI,PBH,PDELTB)

    IMPLICIT NONE

    REAL(8), DIMENSION(0:NLEV), INTENT(OUT) :: PMETADOT
    REAL(8), DIMENSION(NLEV),   INTENT(IN)  :: PU,PV,PDIV
    REAL(8), DIMENSION(NLEV),   INTENT(IN)  :: PDELTAPI,PDELTB
    REAL(8), DIMENSION(0:NLEV), INTENT(IN)  :: PBH
    REAL(8),                    INTENT(IN)  :: PDPISDX,PDPISDY

    INTEGER(8)                              :: J
    REAL(8), DIMENSION(NLEV)                :: ZSDIV_MV

    ZSDIV_MV(1)    = PDELTAPI(1)*PDIV(1) &
                   & + PDELTB(1)*(PU(1)*PDPISDX+PV(1)*PDPISDY)
    DO J = 2, NLEV
       ZSDIV_MV(J) = ZSDIV_MV(J-1) + PDELTAPI(J)*PDIV(J) & 
                   & + PDELTB(J)*(PU(J)*PDPISDX+PV(J)*PDPISDY)
    END DO
    DO J = 1, NLEV-1
       PMETADOT(J) = PBH(J)*ZSDIV_MV(NLEV) - ZSDIV_MV(J)
    END DO
    PMETADOT(0)    = ZERO
    PMETADOT(NLEV) = ZERO

 END SUBROUTINE DIAG_METADOT
 !*********************************************************************************!
 !*********************************************************************************!
 SUBROUTINE DIAG_ETADOT(PETADOTF,PETADOTH,PMF,PMH,PU,PV,PDIV,PDPISDX,PDPISDY, &
                         & PDELTAPI,PDELTB,PBH)

    IMPLICIT NONE

    REAL(8), DIMENSION(NLEV),   INTENT(OUT) :: PETADOTF
    REAL(8), DIMENSION(0:NLEV), INTENT(OUT) :: PETADOTH 
    REAL(8), DIMENSION(NLEV),   INTENT(IN)  :: PU,PV,PDIV
    REAL(8), DIMENSION(NLEV),   INTENT(IN)  :: PDELTAPI
    REAL(8), DIMENSION(NLEV),   INTENT(IN)  :: PDELTB
    REAL(8), DIMENSION(NLEV),   INTENT(IN)  :: PMF
    REAL(8), DIMENSION(0:NLEV), INTENT(IN)  :: PMH,PBH
    REAL(8),                    INTENT(IN)  :: PDPISDX,PDPISDY
    
    INTEGER(8)                              :: J
    REAL(8), DIMENSION(0:NLEV)              :: ZMETADOT
    
    CALL DIAG_METADOT(ZMETADOT,PU,PV,PDIV,PDPISDX,PDPISDY,PDELTAPI,PBH,PDELTB)

    PETADOTH(0)    = ZERO
    DO J = 1, NLEV-1
       PETADOTF(J) = HALF*(ZMETADOT(J)+ZMETADOT(J-1))/PMF(J)
       PETADOTH(J) = ZMETADOT(J)/PMH(J)
    END DO   
    PETADOTF(NLEV) = HALF*(ZMETADOT(NLEV-1))/PMF(NLEV)
    PETADOTH(NLEV) = ZERO
    
 END SUBROUTINE DIAG_ETADOT
 !*******************************************************************************!
 !*******************************************************************************!
 !*******************************************************************************!
 !*******************************************************************************!
 !*******************************************************************************!
 SUBROUTINE DIAG_PHI(KTOP,PHI,PT,PALFA,PDELTA,PZS)

   IMPLICIT NONE

   INTEGER(8),                    INTENT(IN)   :: KTOP
   REAL(8), DIMENSION(KTOP:NLEV), INTENT(OUT)  :: PHI
   REAL(8), DIMENSION(NLEV),      INTENT(IN)   :: PT
   REAL(8), DIMENSION(NLEV),      INTENT(IN)   :: PALFA,PDELTA
   REAL(8),                       INTENT(IN)   :: PZS
      
   INTEGER(8)                                  :: J
   REAL(8)                                     :: ZFULL,ZPHI
   
   ZFULL     = REAL(KTOP,8)

   ZPHI      = RG*PZS
   PHI(NLEV) = ZPHI + ZFULL*RR*PT(NLEV)*PALFA(NLEV)     

   DO J = 0, NLEV-2
      ZPHI          = ZPHI + RR*PT(NLEV-J)*PDELTA(NLEV-J)
      PHI(NLEV-J-1) = ZPHI + ZFULL*RR*PT(NLEV-J-1)*PALFA(NLEV-J-1)
   END DO
   
   ZPHI = ZPHI + RR*PT(1)*PDELTA(1) 
   IF (KTOP==0) PHI(KTOP) = ZPHI      

 END SUBROUTINE DIAG_PHI 
 !*******************************************************************************!
 !*******************************************************************************!
 !*******************************************************************************! 
 !*******************************************************************************!
 SUBROUTINE DIAG_ALPHADELTA(PPIF,PALPHA,PDELTA,PDELTAPI,PPIS,PAH,PBH)

    IMPLICIT NONE

    REAL(8), DIMENSION(NLEV),   INTENT(OUT)  :: PPIF,PDELTAPI
    REAL(8), DIMENSION(NLEV),   INTENT(OUT)  :: PALPHA,PDELTA
    REAL(8), DIMENSION(0:NLEV), INTENT(IN)   :: PAH,PBH
    REAL(8),                    INTENT(IN)   :: PPIS
    
    INTEGER(8)                               :: J
    REAL(8), DIMENSION(0:NLEV)               :: ZPIH

    DO J=0,NLEV
       ZPIH(J) = PAH(J)+ PBH(J)*PPIS
    END DO

    PDELTAPI(1)    = ZPIH(1)-ZPIH(0)
    PALPHA(1)      = ONE
    PDELTA(1)      = ONE + (ONE/RKAPPA)
    PPIF(1)        = PDELTAPI(1)/PDELTA(1)
       
    DO J = 2,NLEV
       PDELTAPI(J) = ZPIH(J)-ZPIH(J-1)
       PPIF(J)     = DSQRT(ZPIH(J)*ZPIH(J-1))
       PALPHA(J)   = ONE-(PPIF(J)/ZPIH(J))
       PDELTA(J)   = PDELTAPI(J)/PPIF(J)
    END DO
          
 END SUBROUTINE DIAG_ALPHADELTA
 !*******************************************************************************!
 !*******************************************************************************!
 !*******************************************************************************!
 !*******************************************************************************!
 !*******************************************************************************!
 SUBROUTINE OP_N(POUT,PIN,PPIS,PDELTAPI)

    IMPLICIT NONE

    REAL(8),                  INTENT(OUT) :: POUT
    REAL(8),                  INTENT(IN)  :: PPIS
    REAL(8), DIMENSION(NLEV), INTENT(IN)  :: PIN,PDELTAPI
    INTEGER                               :: J

    POUT = (PDELTAPI(1)/PPIS)*PIN(1)
    DO J = 2, NLEV
       POUT = POUT + (PDELTAPI(J)/PPIS)*PIN(J)
    END DO

 END SUBROUTINE OP_N
 !*******************************************************************************!
 !*******************************************************************************!
 !*******************************************************************************!
 !*******************************************************************************!
 !*******************************************************************************!
 SUBROUTINE SP_NORMS(P,PNORM,CNORM_TYPE)
   
   IMPLICIT NONE

   REAL(8), DIMENSION(NLEV,GXPT,GYPT),  INTENT(IN)    :: P   
   REAL(8),                             INTENT(OUT)   :: PNORM
   CHARACTER(LEN=2),                    INTENT(IN)    :: CNORM_TYPE
   INTEGER(8)                                         :: NX,NY,JL
   REAL(8)                                            :: ZNORM, ZGPTOT
   REAL(8), DIMENSION(NLEV)                           :: ZMAXLEV

   ZGPTOT = REAL(GYPT,8)*REAL(GXPT,8)*REAL(NLEV,8)

   IF (CNORM_TYPE == 'L2') THEN
     ZNORM = 0.d0
     DO JL=1,NLEV 
       DO NX=1,GXPT
         DO NY =1,GYPT
            ZNORM = ZNORM + P(JL,NX,NY)**2
         END DO
       END DO
     END DO
     PNORM = DSQRT(ZNORM/ZGPTOT)
   ELSE IF (CNORM_TYPE == 'L1') THEN
     ZNORM = 0.d0
     DO JL =1,NLEV
       DO NX=1,GXPT
         DO NY =1,GYPT
            ZNORM = ZNORM + ABS(P(JL,NX,NY))
         END DO
       END DO  
     END DO
     PNORM = (ZNORM/ZGPTOT)
   ELSE
     DO JL =1,NLEV
        ZMAXLEV(JL) = MAXVAL(ABS(P(JL,:,:)))
     END DO   
     PNORM = MAXVAL(ABS(ZMAXLEV))   
   END IF

  END SUBROUTINE SP_NORMS
  !*******************************************************************************!
  !*******************************************************************************!
  !*******************************************************************************!

  
  
END MODULE MOD_DIAGNOSTIC

!=====================================================!
!=====================================================!
!=====================================================!
