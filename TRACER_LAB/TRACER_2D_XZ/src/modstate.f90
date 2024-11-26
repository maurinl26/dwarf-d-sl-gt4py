!#################################  MODELE JOUET 2D  ##################################!
!#                                                                                    #!
!# auteurs : F.Voitus                                                                 #!
!# sujet   : Module servant à faire l'ensemble des calculs dynamiques                 #!
!#                                                                                    #!
!######################################################################################!

!=====================================================================!
!=========================  MODULE PRINCIPAL  ========================!
!=====================================================================!

MODULE MOD_STATE

  USE MOD_SHARE, ONLY : NXPT,NLEV,RR,RG,RPI,ONE,TWO,RDT,ZERO,HALF
  USE MOD_SETUP, ONLY : GEOMETRY,RHSVAR,STATEVAR,BACKGROUND
 
CONTAINS

 !*******************************************************************************!
 !*******************************************************************************!
 !*******************************************************************************!
 !*******************************************************************************! 
  SUBROUTINE INITIAL_STATE(ST0,ST9,STB,GEO)

    USE MOD_SETUP,  ONLY : SET_DOMAIN_GEOMETRY, SET_BACKGROUND_STATE, &
                         & SET_INITIAL_PERTURBATION
    
    IMPLICIT NONE

    TYPE(STATEVAR),         INTENT(INOUT) :: ST0
    TYPE(STATEVAR),         INTENT(OUT)   :: ST9
    TYPE(BACKGROUND),       INTENT(INOUT) :: STB
    TYPE(GEOMETRY),         INTENT(INOUT) :: GEO
    

    
    !* define domain geometry
    CALL SET_DOMAIN_GEOMETRY(GEO)
    !* define background
    CALL SET_BACKGROUND_STATE(STB,GEO)
    !* define initial state
    CALL SET_INITIAL_PERTURBATION(ST0,STB,GEO)
    !* swapping X(t) --> X(t-dt)
    CALL SWAPP_VAR(ST9,ST0)
    
   
  END SUBROUTINE INITIAL_STATE
  !****************************************************************************!
  !****************************************************************************!
  !****************************************************************************!
  !****************************************************************************!
  SUBROUTINE RHS_MDL(RHS,ST,GEO)

    USE MOD_SHARE, ONLY : RDELTR,LSLAG
    USE MOD_DIAG,  ONLY : DIAG_EULERIAN_ORO,DIAG_PIS_MDOT,OP_N
    
    IMPLICIT NONE

    TYPE(RHSVAR),                  INTENT(OUT)    :: RHS
    TYPE(STATEVAR),                INTENT(IN)     :: ST
    TYPE(GEOMETRY),                INTENT(IN)     :: GEO
    
    INTEGER(8)                                    :: I,J
    REAL(8)                                       :: ZPIS_DOT,ZADV_PIS_REF 
    REAL(8), DIMENSION(NLEV)                      :: ZM_DOT
    
    ! Non-linear model tendencies   
    IF (LSLAG) THEN
       DO I = 1, NXPT

          CALL DIAG_PIS_MDOT(ZM_DOT,ZPIS_DOT,ST%PIS(NLEV,I), &
               & ST%M(:,I),ST%DIV(:,I),ST%ETADOTH(:,I),      &
               & GEO%DELA,GEO%DELB,GEO%RDETA)
          CALL DIAG_EULERIAN_ORO(ZADV_PIS_REF,ST%U(:,I),     &
               & GEO%DPISDX_REF(NLEV,I),GEO%DELB)
          
          DO J = 1, NLEV
             RHS%M(J,I) = ZM_DOT(J)
             RHS%Q(J,I) = ZERO   
          END DO
          RHS%PIS(I)    = ZPIS_DOT + RDELTR*ZADV_PIS_REF
          
       END DO
    ELSE
       DO I = 1,NXPT
          
          CALL OP_N(ZPIS_DOT,-ST%DIV(:,I),GEO%DELA)
       
          DO J = 1, NLEV
             RHS%M(J,I) = ZERO
             RHS%Q(J,I) = ZERO   
          END DO
          RHS%PIS(I)    = ZPIS_DOT
          
       END DO
    END IF   
   
  END SUBROUTINE RHS_MDL
  !****************************************************************************!
  !****************************************************************************!
  !****************************************************************************!
  !****************************************************************************!
  SUBROUTINE RHS_EXP(RHS_A,ST,PM0,RHS_M,GEO)

    USE MOD_SHARE, ONLY : LSLAG,RDELTR
    
    IMPLICIT NONE

    TYPE(RHSVAR),     INTENT(INOUT) :: RHS_A
    TYPE(STATEVAR),   INTENT(IN)    :: ST
    REAL(8),          INTENT(IN)    :: PM0
    TYPE(RHSVAR),     INTENT(IN)    :: RHS_M
    TYPE(GEOMETRY),   INTENT(IN)    :: GEO

    INTEGER(8)                      :: I,J
     
    RHS_A%Q(:,:)    = ST%Q(:,:)

    IF (LSLAG) THEN
       RHS_A%M(:,:) = ST%M(:,:)      + PM0*RHS_M%M(:,:)
       RHS_A%PIS(:) = ST%PIS(NLEV,:) + PM0*RHS_M%PIS(:) &
                  & -RDELTR*GEO%PIS_REF(NLEV,:)
    ELSE
       RHS_A%M(:,:) = ST%M(:,:)
       RHS_A%PIS(:) = ST%PIS(NLEV,:) + PM0*RHS_M%PIS(:)
    END IF
    
  END SUBROUTINE RHS_EXP
  !****************************************************************************!
  !****************************************************************************!
  !****************************************************************************!
  !****************************************************************************!
  !****************************************************************************!
  !****************************************************************************!
  SUBROUTINE RHS_IMP(ST,ST1,STB,RHS,GEO)

    USE MOD_SHARE, ONLY : LSLAG,LADV_PIS,RDELTR
    USE MOD_DIAG,  ONLY : OP_N, DEF_M
    USE MOD_SETUP, ONLY : SET_WIND_STATE
    
    IMPLICIT NONE

    TYPE(STATEVAR),                INTENT(INOUT) :: ST
    TYPE(STATEVAR),                INTENT(IN)    :: ST1
    TYPE(BACKGROUND),              INTENT(IN)    :: STB
    TYPE(RHSVAR),                  INTENT(IN)    :: RHS
    TYPE(GEOMETRY),                INTENT(IN)    :: GEO
    
    INTEGER(8)                                   :: I,J
    REAL(8), DIMENSION(NXPT)                     :: ZDIVA,ZDIVB
    REAL(8), DIMENSION(NLEV,NXPT)                :: ZDIVM
    
    !* Compute divergence terms
    DO I=1,NXPT
       CALL OP_N(ZDIVB(I),ST1%DIV(:,I),GEO%DELB)
       DO J=1,NLEV
          ZDIVM(J,I)     = ST1%DIV(J,I)+GEO%RDETA(J)*(         &
                         & ST1%ETADOTH(J,I)-ST1%ETADOTH(J-1,I))
       END DO
       CALL OP_N(ZDIVA(I),ST1%DIV(:,I),GEO%DELA)
    END DO

    !* Time-stepping scheme
    IF (LSLAG) THEN
       DO I = 1,NXPT
          ST%PIS(NLEV,I) = (RHS%PIS(I)-(RDT/TWO)*ZDIVA(I)      & 
                         & +  RDELTR*GEO%PIS_REF(NLEV,I))      &
                         & /(ONE+(RDT/TWO)*ZDIVB(I))
          DO J=1,NLEV
             ST%M(J,I)   = RHS%M(J,I)/(ONE+(RDT/TWO)*ZDIVM(J,I))
          END DO   
       END DO
    ELSE
       DO I = 1,NXPT
          ST%PIS(NLEV,I) = RHS%PIS(I)-(RDT/TWO)*ZDIVA(I)
          ST%M(:,I)      = RHS%M(:,I)
       END DO   
    END IF

    !* Transformation pis <---> m= d\pi/d\eta
    IF (.NOT.LADV_PIS) THEN
       DO I=1,NXPT
          ST%PIS(NLEV,I)    = ZERO
          DO J=1,NLEV
             ST%PIS(NLEV,I) = ST%PIS(NLEV,I)                  &
                            & + ST%M(J,I)*GEO%DETA(J)
          END DO   
       END DO
    ELSE
       DO I=1,NXPT
          CALL DEF_M(ST%M(:,I),ST%PIS(NLEV,I),GEO%DELA,       & 
               & GEO%DELB,GEO%RDETA)
       END DO   
    END IF

    ST%Q(:,:) = RHS%Q(:,:)

    !* computation of wind at t+dt
    CALL SET_WIND_STATE(ST,STB,GEO)

    
  END SUBROUTINE RHS_IMP
  !******************************************************************************!
  !******************************************************************************!
  !******************************************************************************!
  !******************************************************************************!
  SUBROUTINE SWAPP_VAR(STOUT,STIN,LWIND)

    IMPLICIT NONE

    TYPE(STATEVAR),           INTENT(OUT) :: STOUT
    TYPE(STATEVAR),           INTENT(IN)  :: STIN
    LOGICAL,        OPTIONAL, INTENT(IN)  :: LWIND
    
    INTEGER(8)                  :: I
    LOGICAL                     :: LL_WIND
    
    IF (PRESENT(LWIND)) THEN
       LL_WIND = LWIND
    ELSE
       LL_WIND = .FALSE.
    END IF  

    DO I = 1, NXPT
       STOUT%U(:,I)         = STIN%U(:,I)
       STOUT%ETADOT(:,I)    = STIN%ETADOT(:,I)
       STOUT%ETADOTH(:,I)   = STIN%ETADOTH(:,I)
       STOUT%DIV(:,I)       = STIN%DIV(:,I)
    END DO
    
    IF (.NOT.LL_WIND) THEN 
       DO I = 1, NXPT
          STOUT%Q(:,I)         = STIN%Q(:,I)
          STOUT%PIS(NLEV,I)    = STIN%PIS(NLEV,I)
          STOUT%M(:,I)         = STIN%M(:,I)
       END DO
    END IF   

  END SUBROUTINE SWAPP_VAR
!*******************************************************************************!
!*******************************************************************************!
!*******************************************************************************!
!*******************************************************************************!
!###############################################################################!
!###############################################################################!  
!*******************************************************************************!
!*******************************************************************************!
!*******************************************************************************!  
!*******************************************************************************!
!###############################################################################!
!###############################################################################!

END MODULE MOD_STATE
!===============================================================================!
!===============================================================================!
