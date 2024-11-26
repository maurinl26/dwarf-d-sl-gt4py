!##################   PARAMETRES   ###################!
!#                                                   #!
!# auteur : F. Voitus                                #!
!# sujet  : Définit les parametres du modèle         #!
!#                                                   #!
!#####################################################!

!=====================================================!
!================  MODULE PARAMETRES  ================!
!=====================================================!

MODULE MOD_SETUP


  USE MOD_PARAM, ONLY :  GXPT,GYPT,GEXT,ZERO,ONE,TWO,HALF,RDX,RDY,RLX,RLY,RDT, & 
                      &  RG,RCP,RR,RKAPPA,RPI,RU00,RV00,RQ00,RBOYD,RAX,RAY,RP, &
                      &  RY0,RX0,RXS,RYS,RXC,RYC,RXA,RYA,STREAM_FUNCTION,LSPEC,&
                      &  NPERIOX,NPERIOY,LPERIO,LALIAZ,RX1,RY1
                          
  IMPLICIT NONE
  
  ! type prognotic variable
  TYPE PROGVAR
     REAL(8), DIMENSION(:,:), ALLOCATABLE   :: U,V
     REAL(8), DIMENSION(:,:), ALLOCATABLE   :: M
     REAL(8), DIMENSION(:,:), ALLOCATABLE   :: Q     
     REAL(8), DIMENSION(:,:), ALLOCATABLE   :: DUDX, DVDY
  END TYPE PROGVAR
  
  !* type variable rhs 
  TYPE RHSVAR
     REAL(8), DIMENSION(:,:), ALLOCATABLE :: M
     REAL(8), DIMENSION(:,:), ALLOCATABLE :: Q
  END TYPE RHSVAR

  !* type domaine
  TYPE DOMAIN
     REAL(8), DIMENSION(GXPT)       :: X,XFD
     REAL(8), DIMENSION(GYPT)       :: Y,YFD
     REAL(8), DIMENSION(GXPT,GYPT)  :: RELAX
  END TYPE DOMAIN

  !* type smilag
  TYPE SMILAG
     REAL(8),    DIMENSION(:,:),   ALLOCATABLE  :: XDOT,   YDOT
     REAL(8),    DIMENSION(:,:),   ALLOCATABLE  :: XTRAJ,  YTRAJ
     INTEGER(8), DIMENSION(:,:,:), ALLOCATABLE  :: NXLAG,  NYLAG
     REAL(8),    DIMENSION(:,:,:), ALLOCATABLE  :: XWEI,   YWEI
     REAL(8),    DIMENSION(:,:),   ALLOCATABLE  :: DXDOT_DX, DYDOT_DX, ALPHA_DX
     REAL(8),    DIMENSION(:,:),   ALLOCATABLE  :: DXDOT_DY, DYDOT_DY, ALPHA_DY
  END TYPE SMILAG

  !* type mpdata
  TYPE MPDATA
     REAL(8),    DIMENSION(:,:,:), ALLOCATABLE  :: M
     REAL(8),    DIMENSION(:,:),   ALLOCATABLE  :: XDOT
     REAL(8),    DIMENSION(:,:),   ALLOCATABLE  :: YDOT
  END TYPE MPDATA

  
CONTAINS

  !*****************************************************!
  !***************  LISTE DES ROUTINES  ****************!
  !*****************************************************!

  !***************************************************************************************
  !***************************************************************************************
  !***************************************************************************************
  SUBROUTINE SET_DOMAIN_GEOMETRY(DOM)

    USE MOD_PARAM, ONLY : RBOYD,NLBC,LPERIO
    
    IMPLICIT NONE

    TYPE(DOMAIN),               INTENT(INOUT) :: DOM
    
    INTEGER(8)                                :: I,J,L
    REAL(8)                                   :: ZZ,ZX,ZY,ZC
    REAL(8)                                   :: ZRELAX,ZPERIO

    
    !*************************************************************************************!    
    ! Définition de la géométrie du problème 
    DO I = 1,GXPT
       DOM%X(I) = REAL(I-1,8)*RDX     
    END DO
    DO J = 1,GYPT 
       DOM%Y(J) = REAL(J-1,8)*RDY      
    END DO

    !*************************************************************************************!
    !.... Relaxation Coupling 2D XY function
    ZPERIO = ONE
    IF (LPERIO) ZPERIO = ZERO
    DOM%RELAX(NLBC+1:GXPT-GEXT-NLBC,NLBC+1:GYPT-GEXT-NLBC)=ZERO
   
    DO J=1,NLBC
       DO I=1,NLBC
          ZC     = DSQRT(((REAL(NLBC-I,8)/REAL(NLBC-1,8))**2) &
                 & + ((REAL(NLBC-J,8)/REAL(NLBC-1,8))**2) )
          ZRELAX = ZPERIO*((RP+ONE)*(ZC**RP)-RP*(ZC**(RP+ONE)))
          DOM%RELAX(GXPT-GEXT-I,GYPT-GEXT-J)       = ZRELAX    ! NE corner
          DOM%RELAX(I,J)                           = ZRELAX    ! SW corner
          DOM%RELAX(I,GYPT-GEXT-J)                 = ZRELAX    ! NW corner
          DOM%RELAX(GXPT-GEXT-I,J)                 = ZRELAX    ! SE corner 
       END DO
       ZZ      = REAL(NLBC-J,8)/REAL(NLBC-1,8)
       ZRELAX  = ZPERIO*((RP+ONE)*(ZZ**RP)-RP*(ZZ**(RP+ONE)))
       DOM%RELAX(NLBC+1:GXPT-GEXT-NLBC,J)            = ZRELAX
       DOM%RELAX(NLBC+1:GXPT-GEXT-NLBC,GYPT-GEXT-J)  = ZRELAX
       DOM%RELAX(J,NLBC+1:GYPT-GEXT-NLBC)            = ZRELAX
       DOM%RELAX(GXPT-GEXT-J,NLBC+1:GYPT-GEXT-NLBC)  = ZRELAX
    END DO
       
    DOM%RELAX(:,GYPT-GEXT:GYPT)                    = ZPERIO
    DOM%RELAX(GXPT-GEXT:GXPT,1:GYPT-GEXT)          = ZPERIO

    !*************************************************************************************!

    DOM%XFD(:) = HALF
    DOM%YFD(:) = HALF

    IF (.NOT.LPERIO) THEN
       DOM%XFD(1)    = ONE
       DOM%YFD(1)    = ONE
       DOM%XFD(GXPT) = ONE
       DOM%YFD(GYPT) = ONE
    END IF
    
  END SUBROUTINE SET_DOMAIN_GEOMETRY
  !***************************************************************************************!
  !***************************************************************************************!
  !***************************************************************************************!
  !***************************************************************************************!
  !***************************************************************************************!
  SUBROUTINE SET_INITIAL_PERTURBATION(STG,DOM,KCASE)
    
    USE MOD_PARAM, ONLY : STREAM_FUNCTION
    
    IMPLICIT NONE

    TYPE(PROGVAR),             INTENT(INOUT) :: STG
    TYPE(DOMAIN),              INTENT(IN)    :: DOM
    INTEGER,                   INTENT(IN)    :: KCASE
        
    INTEGER                                  :: NX,NY
    REAL(8)                                  :: ZR,ZR2,ZR1,ZR0

    
    IF (KCASE == 1) THEN
       DO NX = 1,GXPT
          DO NY = 1,GYPT    
             ZR = MIN(DSQRT((((DOM%X(NX)-RX0)/RLX)**2)        &
                & +(((DOM%Y(NY)-RY0)/RLY)**2)),ONE)
             STG%Q(NX,NY) = RQ00+(((ONE+DCOS(RPI*ZR1))/TWO)**2)
          END DO   
       END DO
    ELSE IF (KCASE == 2) THEN
       DO NX=1,GXPT
          DO NY=1,GYPT
             ZR0 = 1.5d0
             ZR1 = MAX(ABS(((DOM%X(NX)-RX0)/RLX))   &
                 & ,ABS(((DOM%Y(NY)-RY0)/RLY)))
             ZR  = ZR1/ZR0
             IF (ZR.LE.ONE) THEN
                STG%Q(NX,NY) = ONE
             ELSE
                STG%Q(NX,NY) = ZERO
             END IF   
          END DO   
       END DO
    ELSE IF (KCASE == 3) THEN
       DO NX = 1,GXPT
          DO NY = 1,GYPT    
             ZR0 = (((DOM%X(NX)-RX0))**2)+(((DOM%Y(NY)-RY0))**2)
             ZR1 = (((DOM%X(NX)-RX1))**2)+(((DOM%Y(NY)-RY1))**2)
             STG%Q(NX,NY) = RQ00*DEXP(-ZR0/RXA)+RQ00*DEXP(-ZR1/RYA)
          END DO   
       END DO
    END IF
    
    CALL SET_VELOCITY(STG,DOM,0_8)
        
  END SUBROUTINE SET_INITIAL_PERTURBATION
  !********************************************************************!
  !********************************************************************!
  !********************************************************************!
  !********************************************************************!
  SUBROUTINE SET_VELOCITY(ST,DOM,KSTEP)

    USE MOD_FFT,       ONLY : FFTDER_DXDY
    
    IMPLICIT NONE

    TYPE(PROGVAR),    INTENT(INOUT) :: ST
    TYPE(DOMAIN),     INTENT(IN)    :: DOM
    INTEGER(8),       INTENT(IN)    :: KSTEP

    INTEGER(8)                      :: NX,NY
    REAL(8)                         :: ZX,ZY,ZT
    REAL(8)                         :: ZPSI(GXPT,GYPT)
    REAL(8)                         :: ZDPSIDX(GXPT,GYPT),ZDPSIDY(GXPT,GYPT)

    ZT=REAL(KSTEP,8)*RDT   

    DO NX=1,GXPT
       DO NY=1,GYPT
          ZX=REAL(NX-1,8)*RDX
          ZY=REAL(NX-1,8)*RDY
          ZPSI(NX,NY) = STREAM_FUNCTION(ZX,ZY,ZT)
       END DO
    END DO   

    IF (LSPEC .AND. LPERIO) THEN
       CALL FFTDER_DXDY(ZPSI,ZDPSIDX,ZDPSIDY,LALIAZ)
       ST%U(:,:) =  ZDPSIDY(:,:)
       ST%V(:,:) = -ZDPSIDX(:,:)
       CALL FFTDER_DXDY(ST%U,ST%DUDX,ZDPSIDY,LALIAZ)
       CALL FFTDER_DXDY(ST%V,ZDPSIDX,ST%DVDY,LALIAZ)
    ELSE   
       DO NX=1,GXPT
          DO NY=1,GYPT
             ST%U(NX,NY)    =  DOM%YFD(NY)*(ZPSI(NX,NPERIOY(NY+1)) &
                            & -ZPSI(NX,NPERIOY(NY-1)))/RDY
             ST%V(NX,NY)    = -DOM%XFD(NX)*(ZPSI(NPERIOX(NX+1),NY) &
                            & -ZPSI(NPERIOX(NX-1),NY))/RDX
          END DO
       END DO
       DO NX=1,GXPT
          DO NY=1,GYPT
             ST%DVDY(NX,NY) =  DOM%YFD(NY)*(ST%V(NX,NPERIOY(NY+1)) &
                            & -ST%V(NX,NPERIOY(NY-1)))/RDY
             ST%DUDX(NX,NY) =  DOM%XFD(NX)*(ST%U(NPERIOX(NX+1),NY) &
                            & -ST%U(NPERIOX(NX-1),NY))/RDX
          END DO
       END DO
    END IF
    
  END SUBROUTINE SET_VELOCITY
  !*******************************************************************************!
  !*******************************************************************************!
  !*******************************************************************************!
  !*******************************************************************************!
  !*******************************************************************************!
  !*******************************************************************************!
  !*******************************************************************************! 
  !*******************************************************************************!
  SUBROUTINE SUALLOC_PART1(X1,X0,X9)

    IMPLICIT NONE

    TYPE(PROGVAR),        INTENT(INOUT) :: X1,X0,X9 
    
    ALLOCATE(X0%U(GXPT,GYPT),X0%V(GXPT,GYPT))
    ALLOCATE(X0%DUDX(GXPT,GYPT),X0%DVDY(GXPT,GYPT))
    ALLOCATE(X0%Q(GXPT,GYPT),X0%M(GXPT,GYPT))

    ALLOCATE(X1%U(GXPT,GYPT),X1%V(GXPT,GYPT))
    ALLOCATE(X1%DUDX(GXPT,GYPT),X1%DVDY(GXPT,GYPT))
    ALLOCATE(X1%Q(GXPT,GYPT),X1%M(GXPT,GYPT))

    ALLOCATE(X9%U(GXPT,GYPT),X9%V(GXPT,GYPT))
    ALLOCATE(X9%DUDX(GXPT,GYPT),X9%DVDY(GXPT,GYPT))
    ALLOCATE(X9%Q(GXPT,GYPT),X9%M(GXPT,GYPT))
    
  END SUBROUTINE SUALLOC_PART1
  !*******************************************************************************!
  !*******************************************************************************! 
  !*******************************************************************************!
  SUBROUTINE SUALLOC_PART2(YSL,YMP,RHSA,RHSB)

    IMPLICIT NONE

    TYPE(SMILAG),                  INTENT(INOUT) :: YSL
    TYPE(MPDATA),                  INTENT(INOUT) :: YMP
    TYPE(RHSVAR),                  INTENT(INOUT) :: RHSA,RHSB 
    
    ALLOCATE(RHSA%Q(GXPT,GYPT),RHSA%M(GXPT,GYPT))
    ALLOCATE(RHSB%Q(GXPT,GYPT),RHSB%M(GXPT,GYPT))

    ALLOCATE(YSL%XDOT(GXPT,GYPT),YSL%YDOT(GXPT,GYPT))
    ALLOCATE(YSL%XTRAJ(GXPT,GYPT),YSL%YTRAJ(GXPT,GYPT))
    ALLOCATE(YSL%NXLAG(6,GXPT,GYPT),YSL%NYLAG(6,GXPT,GYPT))
    ALLOCATE(YSL%XWEI(6,GXPT,GYPT),YSL%YWEI(6,GXPT,GYPT))
    ALLOCATE(YSL%DXDOT_DX(GXPT,GYPT),YSL%DYDOT_DX(GXPT,GYPT))
    ALLOCATE(YSL%DXDOT_DY(GXPT,GYPT),YSL%DYDOT_DY(GXPT,GYPT))
    ALLOCATE(YSL%ALPHA_DX(GXPT,GYPT),YSL%ALPHA_DY(GXPT,GYPT))

    ALLOCATE(YMP%XDOT(GXPT,GYPT),YMP%YDOT(GXPT,GYPT))
    ALLOCATE(YMP%M(0:1,GXPT,GYPT))
    
    
  END SUBROUTINE SUALLOC_PART2
  !**********************************************************************************
  !**********************************************************************************
  !**********************************************************************************
  !**********************************************************************************
  !**********************************************************************************
  !**********************************************************************************
  SUBROUTINE FREEMEM(X0,X9,X1,YSL,YMP,RHSA,RHSB)

    IMPLICIT NONE

    TYPE(PROGVAR),                 INTENT(INOUT) :: X1,X0,X9
    TYPE(RHSVAR),                  INTENT(INOUT) :: RHSA,RHSB
    TYPE(SMILAG),                  INTENT(INOUT) :: YSL
    TYPE(MPDATA),                  INTENT(INOUT) :: YMP

    DEALLOCATE(X0%U,X0%V,X0%Q,X0%M,X0%DUDX,X0%DVDY)
    DEALLOCATE(X9%U,X9%V,X9%Q,X9%M,X9%DUDX,X9%DVDY)
    DEALLOCATE(X1%U,X1%V,X1%Q,X1%M,X1%DUDX,X1%DVDY)
    DEALLOCATE(RHSA%M,RHSA%Q,RHSB%M,RHSB%Q)
    
    DEALLOCATE(YSL%XDOT,YSL%YDOT)
    DEALLOCATE(YSL%XDOT,YSL%YDOT)
    DEALLOCATE(YSL%XTRAJ,YSL%YTRAJ)
    DEALLOCATE(YSL%NXLAG,YSL%NYLAG)
    DEALLOCATE(YSL%XWEI,YSL%YWEI)
    DEALLOCATE(YSL%DXDOT_DX,YSL%DYDOT_DX)
    DEALLOCATE(YSL%DXDOT_DY,YSL%DYDOT_DY)
    DEALLOCATE(YSL%ALPHA_DX,YSL%ALPHA_DY)
    DEALLOCATE(YMP%XDOT,YMP%YDOT,YMP%M)
    
  END SUBROUTINE FREEMEM
  !**********************************************************************************
  !**********************************************************************************
  !**********************************************************************************

END MODULE MOD_SETUP

!=====================================================!
!=====================================================!
!=====================================================!
