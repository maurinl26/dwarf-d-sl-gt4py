!#################################  MODELE JOUET 2D  ##################################!
!#                                                                                    #!
!# auteurs : F.Voitus                                                                 #!
!# sujet   : Module servant à faire l'ensemble des calculs dynamiques                 #!
!#                                                                                    #!
!######################################################################################!

!=====================================================================!
!=========================  MODULE PRINCIPAL  ========================!
!=====================================================================!

MODULE MOD_SLAG
  
  USE MOD_SETUP, ONLY : SMILAG,DOMAIN,RHSVAR
  USE MOD_PARAM, ONLY : GXPT,GYPT,RPI,ONE,TWO,ZERO,HALF,RQ,RDT,RDX,RDY,      &
                      & NITMP,NPERIOX,NPERIOY,LPERIO,LPRINTLEV,LNESC,LCOMAD, &
                      & LSLTRAJ_PREC,LSETTLS,NCOMAD_OPT,C_SLTRAJ,C_SLINTP
 
CONTAINS

  !********************************************************************************************!
  !********************************************************************************************!
  !********************************************************************************************!
  !*******************************  LISTE DES ROUTINES ****************************************!
  !********************************************************************************************!  
  !********************************************************************************************!
  !********************************************************************************************!
  SUBROUTINE SMILAG_TRANSPORT_SCHEME(RHSG_DEP,RHSG_ARR,SLSTRUC,   &
             & PXDOT9,PYDOT9,PXDOT0,PYDOT0,PDXDOT0_DX,PDYDOT0_DY, &
             & PXDOT1,PYDOT1,PDXDOT1_DX,PDYDOT1_DY,KSTEP,CD_ADV)
    
  IMPLICIT NONE

  TYPE(SMILAG),                        INTENT(INOUT)  :: SLSTRUC
  TYPE(RHSVAR),                        INTENT(INOUT)  :: RHSG_ARR
  TYPE(RHSVAR),                        INTENT(INOUT)  :: RHSG_DEP
  REAL(8),                             INTENT(IN)     :: PXDOT9(GXPT,GYPT),PYDOT9(GXPT,GYPT)
  REAL(8),                             INTENT(IN)     :: PXDOT0(GXPT,GYPT),PDXDOT0_DX(GXPT,GYPT)
  REAL(8),                             INTENT(IN)     :: PYDOT0(GXPT,GYPT),PDYDOT0_DY(GXPT,GYPT)
  REAL(8),                             INTENT(IN)     :: PXDOT1(GXPT,GYPT),PDXDOT1_DX(GXPT,GYPT)
  REAL(8),                             INTENT(IN)     :: PYDOT1(GXPT,GYPT),PDYDOT1_DY(GXPT,GYPT)
  INTEGER(8),                          INTENT(IN)     :: KSTEP
  CHARACTER(LEN=1),                    INTENT(IN)     :: CD_ADV

  INTEGER(8)                                          :: I,J,K
  REAL(8)                                             :: ZLIN,ZWIND_HOR_MAX
  REAL(8)                                             :: XDOT_F(GXPT,GYPT),YDOT_F(GXPT,GYPT)

  ZLIN = ONE
  IF (CD_ADV == 'LI') ZLIN =ZERO
  
  IF (LNESC) THEN   
     XDOT_F(:,:)    = ZLIN*PXDOT0(:,:)
     YDOT_F(:,:)    = ZLIN*PYDOT0(:,:)
  ELSE IF (LSETTLS) THEN
     XDOT_F(:,:)    = ZLIN*(TWO*PXDOT0(:,:)-PXDOT9(:,:))
     YDOT_F(:,:)    = ZLIN*(TWO*PYDOT0(:,:)-PYDOT9(:,:))
  ELSE
     XDOT_F(:,:)    = ZLIN*PXDOT1(:,:)
     YDOT_F(:,:)    = ZLIN*PYDOT1(:,:)
  ENDIF
 
  SLSTRUC%XDOT(:,:) = ZLIN*PXDOT0(:,:)
  SLSTRUC%YDOT(:,:) = ZLIN*PYDOT0(:,:)
  
  !-----------------------------------------------------------------!
  ! STABILITY CHECK-POINT MAXIMUM OF HORIZONTAL WIND VELOCITY       !
  ZWIND_HOR_MAX = MAX(MAXVAL(SLSTRUC%XDOT),MAXVAL(SLSTRUC%YDOT))    !
  IF (LPRINTLEV) THEN                                               !
     WRITE(*,*) 'NSTEP : ',KSTEP,                         &         !
              & '--> MAX.WIND_HOR : ',ZWIND_HOR_MAX                 !
  END IF                                                            !
  IF (ZWIND_HOR_MAX .GE. 300.d0) THEN                               !
     WRITE(*,*) 'WIND TOO STRONG EXPLOSION AT TIMESTEP = ',KSTEP    !
     STOP                                                           !
  END IF                                                            !
  !-----------------------------------------------------------------!

  SLSTRUC%DXDOT_DX(:,:) = ZERO
  SLSTRUC%DXDOT_DY(:,:) = ZERO
  SLSTRUC%DYDOT_DX(:,:) = ZERO
  SLSTRUC%DYDOT_DY(:,:) = ZERO
        
  !---------------------------------------------------------------------!
  ! COMPUTE DEFORMATIONAL TERMS
  IF (LSLTRAJ_PREC) THEN
     CALL COMPUTE_DEFORMATIONAL(SLSTRUC%DXDOT_DX,SLSTRUC%DXDOT_DY,      &
          & SLSTRUC%DYDOT_DX,SLSTRUC%DYDOT_DY,SLSTRUC%XDOT,SLSTRUC%YDOT)
  END IF
     
  !---------------------------------------------------------------------!
  ! RESEARCH OF DEPARTURE POINTS AND INTERPOLATION WEIGHTS 
  CALL LARCINA(SLSTRUC,XDOT_F,YDOT_F)

  IF (LCOMAD) THEN
     !Compute comad stretching coefficients
     CALL COMPUTE_COMAD_STRETCHING(SLSTRUC%ALPHA_DX,SLSTRUC%ALPHA_DY,   &
          & SLSTRUC%XDOT,SLSTRUC%YDOT,PDXDOT0_DX,PDYDOT0_DY,            &
          & XDOT_F,YDOT_F,PDXDOT1_DX,PDYDOT1_DY,                        &
          & SLSTRUC%NXLAG,SLSTRUC%XWEI,SLSTRUC%NYLAG,SLSTRUC%YWEI)
     !Compute new SL weights for Sources
     CALL ELASCAW_2D(SLSTRUC%NXLAG,SLSTRUC%NYLAG,SLSTRUC%XWEI,          &
          & SLSTRUC%YWEI,SLSTRUC%XTRAJ,SLSTRUC%YTRAJ,                   &
          & PALPHA_DX=SLSTRUC%ALPHA_DX,PALPHA_DY=SLSTRUC%ALPHA_DY,      &
          & LDCOMAD=.TRUE.)
  END IF
  !---------------------------------------------------------------------!
  ! SL INTERPOLATIONS AT DEPARTURE POINTS                  
  IF (CD_ADV == 'LI') THEN 
     RHSG_DEP%Q(:,:) = RHSG_ARR%Q(:,:)
     RHSG_DEP%M(:,:) = RHSG_ARR%M(:,:)
  ELSE 
     CALL LARCINB(RHSG_DEP%Q,RHSG_ARR%Q,RHSG_DEP%M,RHSG_ARR%M,          &
          & SLSTRUC%XWEI,SLSTRUC%NXLAG,SLSTRUC%YWEI,SLSTRUC%NYLAG)
  END IF
     
  END SUBROUTINE SMILAG_TRANSPORT_SCHEME
  !################################################################################################
  !################################################################################################
  !################################################################################################
  !################################################################################################  
  SUBROUTINE LARCINA(SLSTRUC,XDOT_F,YDOT_F)  

    IMPLICIT NONE
               
    TYPE(SMILAG),                       INTENT(INOUT) :: SLSTRUC
    REAL(8), DIMENSION(GXPT,GYPT),      INTENT(IN)    :: XDOT_F,       YDOT_F

    REAL(8), DIMENSION(GXPT,GYPT)                     :: XDOT_DEP,     YDOT_DEP
    REAL(8), DIMENSION(GXPT,GYPT)                     :: DXDOT_DX_DEP, DXDOT_DY_DEP
    REAL(8), DIMENSION(GXPT,GYPT)                     :: DYDOT_DX_DEP, DYDOT_DY_DEP
    INTEGER(8)                                        :: JITER

    !* starting point for SL departure point 
    XDOT_DEP(:,:)     = SLSTRUC%XDOT(:,:)
    YDOT_DEP(:,:)     = SLSTRUC%YDOT(:,:)
    
    DXDOT_DX_DEP(:,:) = SLSTRUC%DXDOT_DX(:,:)
    DXDOT_DY_DEP(:,:) = SLSTRUC%DXDOT_DY(:,:)
    DYDOT_DX_DEP(:,:) = SLSTRUC%DYDOT_DX(:,:)  
    DYDOT_DY_DEP(:,:) = SLSTRUC%DYDOT_DY(:,:)
    
    DO JITER=1,NITMP
       !Compute SL trajectories
       CALL ELARCHE(SLSTRUC%XTRAJ,SLSTRUC%YTRAJ,                      & 
            & XDOT_F,YDOT_F,XDOT_DEP,YDOT_DEP,                        &
            & PDXDOT_DX=DXDOT_DX_DEP,PDXDOT_DY=DXDOT_DY_DEP,          &
            & PDYDOT_DX=DYDOT_DX_DEP,PDYDOT_DY=DYDOT_DY_DEP)
       !Compute SL weights
       CALL ELASCAW_2D(SLSTRUC%NXLAG,SLSTRUC%NYLAG,SLSTRUC%XWEI,      &
            & SLSTRUC%YWEI,SLSTRUC%XTRAJ,SLSTRUC%YTRAJ)
       !Interpolations to Departure points
       CALL ELARMES(XDOT_DEP,YDOT_DEP,DXDOT_DX_DEP,DXDOT_DY_DEP,      &
            & DYDOT_DX_DEP,DYDOT_DY_DEP,SLSTRUC%XDOT,SLSTRUC%YDOT,    &
            & SLSTRUC%DXDOT_DX,SLSTRUC%DXDOT_DY,SLSTRUC%DYDOT_DX,     &
            & SLSTRUC%DYDOT_DY,SLSTRUC%NXLAG,SLSTRUC%NYLAG,           &
            & SLSTRUC%XWEI,SLSTRUC%YWEI)
    END DO

  END SUBROUTINE LARCINA
  !*******************************************************************************!
  !*******************************************************************************!
  !*******************************************************************************!
  !*******************************************************************************!
  !*******************************************************************************!
  !*******************************************************************************!
  !*******************************************************************************!
  !*******************************************************************************!
  !*******************************************************************************!
  SUBROUTINE ELARCHE(PXTRAJ,PYTRAJ,PXDOT_ARR,PYDOT_ARR,PXDOT_DEP,PYDOT_DEP, &
                   & PDXDOT_DX,PDXDOT_DY,PDYDOT_DX,PDYDOT_DY)
    
    IMPLICIT NONE

    REAL(8),DIMENSION(GXPT,GYPT), INTENT(OUT) :: PXTRAJ,PYTRAJ      
    REAL(8),DIMENSION(GXPT,GYPT), INTENT(IN)  :: PXDOT_ARR,PYDOT_ARR
    REAL(8),DIMENSION(GXPT,GYPT), INTENT(IN)  :: PXDOT_DEP,PYDOT_DEP
    REAL(8), OPTIONAL,            INTENT(IN)  :: PDXDOT_DX(GXPT,GYPT),PDYDOT_DY(GXPT,GYPT)
    REAL(8), OPTIONAL,            INTENT(IN)  :: PDXDOT_DY(GXPT,GYPT),PDYDOT_DX(GXPT,GYPT)

    INTEGER(8)                                :: II,  JJ,  KK
    REAL(8)                                   :: ZTRAJX  ,ZTRAJY  ,ZDET   ,ZDET0  


    ZDET0 = ABS(100*EPSILON(ONE))

    ! COMPUTE SL TRAJECTORIES AT DEPARTURE POINT
    IF (LSLTRAJ_PREC) THEN
       DO II=1,GXPT
          DO JJ=1,GYPT
             ZDET             = ONE+(RDT/TWO)*(PDXDOT_DX(II,JJ)+PDYDOT_DY(II,JJ))   &
                              & + ((RDT/TWO)**2)*(PDXDOT_DX(II,JJ)*PDYDOT_DY(II,JJ) &
                              & - PDXDOT_DY(II,JJ)*PDYDOT_DX(II,JJ))
             IF (ABS(ZDET) > ZDET0) THEN
                ZTRAJX        = REAL(II-1,8)*RDX                                    &
                              & - (RDT/TWO)*( PXDOT_ARR(II,JJ)+PXDOT_DEP(II,JJ) )   &
                              & + (RDT/TWO)*( PDXDOT_DX(II,JJ)*PXTRAJ(II,JJ)        &
                              & + PDXDOT_DY(II,JJ)*PYTRAJ(II,JJ) ) 
                ZTRAJY        = REAL(JJ-1,8)*RDY                                    &
                              & - (RDT/TWO)*(PYDOT_ARR(II,JJ)+PYDOT_DEP(II,JJ))     &
                              & + (RDT/TWO)*( PDYDOT_DY(II,JJ)*PYTRAJ(II,JJ)        &
                              & + PDYDOT_DX(II,JJ)*PXTRAJ(II,JJ) )
            
                PXTRAJ(II,JJ) = ( (ONE+(RDT/TWO)*PDYDOT_DY(II,JJ))*ZTRAJX           &
                              & - (RDT/TWO)*PDXDOT_DY(II,JJ)*ZTRAJY )/ZDET
                PYTRAJ(II,JJ) = ( (ONE+(RDT/TWO)*PDXDOT_DX(II,JJ))*ZTRAJY           &
                              & - (RDT/TWO)*PDYDOT_DX(II,JJ)*ZTRAJX )/ZDET
             ELSE
                PXTRAJ(II,JJ) = REAL(II-1,8)*RDX                                    &
                              & -(RDT/TWO)*(PXDOT_ARR(II,JJ)+PXDOT_DEP(II,JJ))
                PYTRAJ(II,JJ) = REAL(JJ-1,8)*RDY                                    &
                              & -(RDT/TWO)*(PYDOT_ARR(II,JJ)+PYDOT_DEP(II,JJ))
             END IF   
          ENDDO
       ENDDO
    ELSE
       DO II=1,GXPT
          DO JJ=1,GYPT   
             PXTRAJ(II,JJ)    = REAL(II-1,8)*RDX                                    &
                              & -(RDT/TWO)*(PXDOT_ARR(II,JJ)+PXDOT_DEP(II,JJ))
             PYTRAJ(II,JJ)    = REAL(JJ-1,8)*RDY                                    &
                              & -(RDT/TWO)*(PYDOT_ARR(II,JJ)+PYDOT_DEP(II,JJ))
          ENDDO
       ENDDO
    END IF
    
END SUBROUTINE ELARCHE
!****************************************************************************************!
!****************************************************************************************!
!****************************************************************************************!
!****************************************************************************************!
!****************************************************************************************!
!****************************************************************************************!
!****************************************************************************************!
SUBROUTINE ELASCAW_2D(KXLAG,KYLAG,PXWEI,PYWEI,PXTRAJ,PYTRAJ,&
           & PALPHA_DX,PALPHA_DY,LDCOMAD)

IMPLICIT NONE

REAL(8),   DIMENSION(6,GXPT,GYPT),          INTENT(OUT) :: PXWEI,PYWEI
INTEGER(8),DIMENSION(6,GXPT,GYPT),          INTENT(OUT) :: KXLAG,KYLAG
REAL(8),   DIMENSION(GXPT,GYPT),            INTENT(IN)  :: PXTRAJ,PYTRAJ
REAL(8),   OPTIONAL, DIMENSION(GXPT,GYPT),  INTENT(IN)  :: PALPHA_DX,PALPHA_DY
LOGICAL,   OPTIONAL,                        INTENT(IN)  :: LDCOMAD

INTEGER(8)                     :: II       ,JJ       ,LL       ,KK
INTEGER(8)                     :: IP       ,JP             
REAL(8)                        :: ZXP1     ,ZX0      ,ZXM1     ,ZXM2     ,ZXD 
REAL(8)                        :: ZYP1     ,ZY0      ,ZYM1     ,ZYM2     ,ZYD      
REAL(8)                        :: ZEPS     ,ZALX     ,ZALY     ,ZWXC     ,ZWYC
LOGICAL                        :: LL_COMAD

ZEPS = ABS(EPSILON(ONE))


IF (PRESENT(LDCOMAD)) THEN
   LL_COMAD = LDCOMAD
ELSE
   LL_COMAD =.FALSE.
END IF   

!****************************************!
! COMPUTE WEIGHTS FOR SL INTERPOLATIONS *!
!****************************************!

IF (LL_COMAD) THEN
  DO II=1,GXPT   
     DO JJ=1,GYPT
     
           !**************
           ! X-DIRECTION !
           !**************
           ZXD   = PXTRAJ(II,JJ) 
           ZALX  = REAL(II-1,8)-(ZXD/RDX)
           IF (ABS(ZALX).LT.ZEPS) ZALX = ZERO

           IP = FLOOR(ZALX)

           ZXM1=REAL(II-IP-1-1,8)*RDX
           ZWXC=((ZXD-ZXM1)/RDX)+(((ZXD-ZXM1)/RDX)-HALF)*(PALPHA_DX(II,JJ)-ONE)
           
           ! cubic interpolation between II-IP-2, II-IP-1, II-IP and II-IP+1
           DO LL=1,4
              KXLAG(LL,II,JJ) = NPERIOX(II-IP+2-LL) 
           END DO  
           PXWEI(1,II,JJ)=((ZWXC**2)-ONE)*(ZWXC/6.d0)
           PXWEI(2,II,JJ)=ZWXC+((ZWXC**2)/TWO)*(ONE-ZWXC)
           PXWEI(3,II,JJ)=ONE-(ZWXC**2)+(ZWXC/TWO)*((ZWXC**2)-ONE)
           PXWEI(4,II,JJ)=ONE-PXWEI(1,II,JJ)-PXWEI(2,II,JJ)-PXWEI(3,II,JJ)

           ! Linear interpolating weights between II-IP-1, II-IP
           DO LL=5,6
            KXLAG(LL,II,JJ) = NPERIOX(II-IP+5-LL) 
           END DO
           PXWEI(5,II,JJ)=ZWXC
           PXWEI(6,II,JJ)=ONE-ZWXC
           !---------------------------------------------------------------------------------------
           
           !**************
           ! Y-DIRECTION !
           !**************
           ZYD   = PYTRAJ(II,JJ) 
           ZALY  = REAL(JJ-1,8)-(ZYD/RDY)
           IF (ABS(ZALY).LT.ZEPS) ZALY = ZERO
           
           JP = FLOOR(ZALY)

           ZYM1=REAL(JJ-JP-1-1,8)*RDY
           ZWYC=((ZYD-ZYM1)/RDX)+(((ZYD-ZYM1)/RDX)-HALF)*(PALPHA_DY(II,JJ)-ONE)
           
           ! cubic interpolation between JJ-JP-2, JJ-JP-1, JJ-JP and JJ-JP+1
           DO LL=1,4
              KYLAG(LL,II,JJ) = NPERIOY(JJ-JP+2-LL) 
           END DO
           PYWEI(1,II,JJ)=((ZWYC**2)-ONE)*(ZWYC/6.d0)
           PYWEI(2,II,JJ)=ZWYC+((ZWYC**2)/TWO)*(ONE-ZWYC)
           PYWEI(3,II,JJ)=ONE-(ZWYC**2)+(ZWYC/TWO)*((ZWYC**2)-ONE)
           PYWEI(4,II,JJ)=ONE-PYWEI(1,II,JJ)-PYWEI(2,II,JJ)-PYWEI(3,II,JJ)
           
           ! Linear interpolating weights between JJ-JP-1, JJ-JP 
           DO LL=5,6
              KYLAG(LL,II,JJ) = NPERIOY(JJ-JP+5-LL) 
           END DO
           PYWEI(5,II,JJ)=ZWYC
           PYWEI(6,II,JJ)=ONE-ZWYC
           !---------------------------------------------------------------------------------------
           
     END DO
  END DO
ELSE   !********** NOT COMAD ******!
  DO II=1,GXPT   
     DO JJ=1,GYPT
     
           !**************
           ! X-DIRECTION !
           !**************
           ZXD   = PXTRAJ(II,JJ) 
           ZALX  = REAL(II-1,8)-(ZXD/RDX)
           IF (ABS(ZALX).LT.ZEPS) ZALX = ZERO
         
           IP = FLOOR(ZALX)
           
           ZXP1=REAL(II-IP+1-1,8)*RDX
           ZX0 =REAL(II-IP+0-1,8)*RDX
           ZXM1=REAL(II-IP-1-1,8)*RDX
           ZXM2=REAL(II-IP-2-1,8)*RDX

           ! cubic interpolation between II-IP-2, II-IP-1, II-IP and II-IP+1
           DO LL=1,4
              KXLAG(LL,II,JJ) = NPERIOX(II-IP+2-LL) 
           END DO
           PXWEI(1,II,JJ)=((ZXD-ZX0)/(ZXP1-ZX0))*((ZXD-ZXM1)/(ZXP1-ZXM1))*((ZXD-ZXM2)/(ZXP1-ZXM2))
           PXWEI(2,II,JJ)=((ZXD-ZXP1)/(ZX0-ZXP1))*((ZXD-ZXM1)/(ZX0-ZXM1))*((ZXD-ZXM2)/(ZX0-ZXM2))
           PXWEI(3,II,JJ)=((ZXD-ZXP1)/(ZXM1-ZXP1))*((ZXD-ZX0)/(ZXM1-ZX0))*((ZXD-ZXM2)/(ZXM1-ZXM2))
           PXWEI(4,II,JJ)=((ZXD-ZXP1)/(ZXM2-ZXP1))*((ZXD-ZX0)/(ZXM2-ZX0))*((ZXD-ZXM1)/(ZXM2-ZXM1))
           ! Linear interpolating weights between II-IP-1, II-IP
           DO LL=5,6
            KXLAG(LL,II,JJ) = NPERIOX(II-IP+5-LL) 
           END DO
           PXWEI(5,II,JJ)=(ZXD-ZXM1)/(ZX0-ZXM1)
           PXWEI(6,II,JJ)=(ZXD-ZX0)/(ZXM1-ZX0)
           !---------------------------------------------------------------------------------------
           
           !**************
           ! Y-DIRECTION !
           !**************

           ZYD  = PYTRAJ(II,JJ) 
           ZALY = REAL(JJ-1,8)-(ZYD/RDY)
           IF (ABS(ZALY).LT.ZEPS) ZALY = ZERO
           
           JP = FLOOR(ZALY)
           
           ZYP1=REAL(JJ-JP+1-1,8)*RDY
           ZY0 =REAL(JJ-JP+0-1,8)*RDY
           ZYM1=REAL(JJ-JP-1-1,8)*RDY
           ZYM2=REAL(JJ-JP-2-1,8)*RDY
           
           ! cubic interpolation between JJ-JP-2, JJ-JP-1, JJ-JP and JJ-JP+1
           DO LL=1,4
              KYLAG(LL,II,JJ) = NPERIOY(JJ-JP+2-LL) 
           END DO
           PYWEI(1,II,JJ)=((ZYD-ZY0)/(ZYP1-ZY0))*((ZYD-ZYM1)/(ZYP1-ZYM1))*((ZYD-ZYM2)/(ZYP1-ZYM2))
           PYWEI(2,II,JJ)=((ZYD-ZYP1)/(ZY0-ZYP1))*((ZYD-ZYM1)/(ZY0-ZYM1))*((ZYD-ZYM2)/(ZY0-ZYM2))
           PYWEI(3,II,JJ)=((ZYD-ZYP1)/(ZYM1-ZYP1))*((ZYD-ZY0)/(ZYM1-ZY0))*((ZYD-ZYM2)/(ZYM1-ZYM2))
           PYWEI(4,II,JJ)=((ZYD-ZYP1)/(ZYM2-ZYP1))*((ZYD-ZY0)/(ZYM2-ZY0))*((ZYD-ZYM1)/(ZYM2-ZYM1))
           ! Linear interpolating weights between JJ-JP-1, JJ-JP
           DO LL=5,6
              KYLAG(LL,II,JJ) = NPERIOY(JJ-JP+5-LL) 
           END DO
           PYWEI(5,II,JJ)=(ZYD-ZYM1)/(ZY0-ZYM1)
           PYWEI(6,II,JJ)=(ZYD-ZY0)/(ZYM1-ZY0)
           !---------------------------------------------------------------------------------------
     END DO
  END DO
END IF
  
END SUBROUTINE ELASCAW_2D
!***********************************************************************************************!
!***********************************************************************************************!
!***********************************************************************************************!
!***********************************************************************************************!
SUBROUTINE ELARMES(PXDOT_DEP,PYDOT_DEP,PDXDOT_DX_DEP,PDXDOT_DY_DEP,PDYDOT_DX_DEP,PDYDOT_DY_DEP, &
                 & PXDOT_ARR,PYDOT_ARR,PDXDOT_DX_ARR,PDXDOT_DY_ARR,PDYDOT_DX_ARR,PDYDOT_DY_ARR, &
                 & KXLAG,KYLAG,PXWEI,PYWEI)

IMPLICIT NONE

REAL(8),DIMENSION(GXPT,GYPT),      INTENT(OUT) :: PXDOT_DEP     ,PYDOT_DEP
REAL(8),DIMENSION(GXPT,GYPT),      INTENT(OUT) :: PDXDOT_DX_DEP ,PDXDOT_DY_DEP
REAL(8),DIMENSION(GXPT,GYPT),      INTENT(OUT) :: PDYDOT_DX_DEP ,PDYDOT_DY_DEP
REAL(8),DIMENSION(GXPT,GYPT),      INTENT(IN)  :: PXDOT_ARR     ,PYDOT_ARR
REAL(8),DIMENSION(GXPT,GYPT),      INTENT(OUT) :: PDXDOT_DX_ARR ,PDXDOT_DY_ARR
REAL(8),DIMENSION(GXPT,GYPT),      INTENT(OUT) :: PDYDOT_DX_ARR ,PDYDOT_DY_ARR
REAL(8),DIMENSION(6,GXPT,GYPT),    INTENT(IN)  :: PXWEI         ,PYWEI    
INTEGER(8),DIMENSION(6,GXPT,GYPT), INTENT(IN)  :: KXLAG         ,KYLAG     

INTEGER(8)                    :: KK  ,II  ,JJ
INTEGER(8)                    :: K   ,I   ,J
INTEGER                       :: ISTA   ,IEND

IF (C_SLTRAJ == "C") THEN
   ISTA = 1
   IEND = 4
ELSE IF (C_SLTRAJ == "L") THEN
   ISTA = 5
   IEND = 6
ELSE
   WRITE(*,*) "UNKNOWN TRAJECTORY INTERPOLATION ORDER : ",C_SLTRAJ
   STOP
END IF

!* interpolated advecting velocities
DO II=1,GXPT
   DO JJ=1,GYPT
      PXDOT_DEP(II,JJ) = ZERO
      PYDOT_DEP(II,JJ) = ZERO
      DO I=ISTA,IEND
         DO J=ISTA,IEND
            PXDOT_DEP(II,JJ) = PXDOT_DEP(II,JJ)               &
                 & + PXWEI(I,II,JJ)*PYWEI(J,II,JJ)*(          &
                 &   PXDOT_ARR(KXLAG(I,II,JJ),KYLAG(J,II,JJ)) )
            PYDOT_DEP(II,JJ) = PYDOT_DEP(II,JJ)               & 
                 & + PXWEI(I,II,JJ)*PYWEI(J,II,JJ)*(          &
                 &   PYDOT_ARR(KXLAG(I,II,JJ),KYLAG(J,II,JJ)) )
         ENDDO
      ENDDO
   ENDDO
ENDDO

!* set deformational term to zero
PDXDOT_DX_DEP(:,:) = ZERO
PDXDOT_DY_DEP(:,:) = ZERO
PDYDOT_DX_DEP(:,:) = ZERO
PDYDOT_DY_DEP(:,:) = ZERO

IF (LSLTRAJ_PREC) THEN
  DO II=1,GXPT
     DO JJ=1,GYPT
        DO I=ISTA,IEND
           DO J=ISTA,IEND
              PDXDOT_DX_DEP(II,JJ) = PDXDOT_DX_DEP(II,JJ)           &
                   & + PXWEI(I,II,JJ)*PYWEI(J,II,JJ)*(              &
                   &   PDXDOT_DX_ARR(KXLAG(I,II,JJ),KYLAG(J,II,JJ)) )
              PDXDOT_DY_DEP(II,JJ) = PDXDOT_DY_DEP(II,JJ)           &
                   & + PXWEI(I,II,JJ)*PYWEI(J,II,JJ)*(              &
                   &   PDXDOT_DY_ARR(KXLAG(I,II,JJ),KYLAG(J,II,JJ)) )
              PDYDOT_DX_DEP(II,JJ) = PDYDOT_DX_DEP(II,JJ)           &
                   & + PXWEI(I,II,JJ)*PYWEI(J,II,JJ)*(              &
                   &   PDYDOT_DX_ARR(KXLAG(I,II,JJ),KYLAG(J,II,JJ)) )
              PDYDOT_DY_DEP(II,JJ) = PDYDOT_DY_DEP(II,JJ)           &
                   & + PXWEI(I,II,JJ)*PYWEI(J,II,JJ)*(              &
                   &   PDYDOT_DY_ARR(KXLAG(I,II,JJ),KYLAG(J,II,JJ)) )
           ENDDO
        ENDDO
     ENDDO
  ENDDO
END IF
  
END SUBROUTINE ELARMES
!**********************************************************************************************!
!**********************************************************************************************!
!**********************************************************************************************!
!**********************************************************************************************!
!**********************************************************************************************!
!**********************************************************************************************!
!**********************************************************************************************!
!**********************************************************************************************!
SUBROUTINE LARCINB(RQ_DEP,RQ_ARR,RM_DEP,RM_ARR,PXWEI,KXLAG,PYWEI,KYLAG)
  
IMPLICIT NONE

REAL(8),          INTENT(OUT) :: RQ_DEP(GXPT,GYPT),RM_DEP(GXPT,GYPT)
REAL(8),          INTENT(IN)  :: RQ_ARR(GXPT,GYPT),RM_ARR(GXPT,GYPT)
REAL(8),          INTENT(IN)  :: PXWEI(6,GXPT,GYPT),PYWEI(6,GXPT,GYPT)
INTEGER(8),       INTENT(IN)  :: KXLAG(6,GXPT,GYPT),KYLAG(6,GXPT,GYPT)

INTEGER(8)                    :: II ,JJ   ,I  ,J
INTEGER                       :: ISTA   ,IEND

IF (C_SLINTP == "C") THEN
   ISTA = 1
   IEND = 4
ELSE IF (C_SLINTP == "L") THEN
   ISTA = 5
   IEND = 6
ELSE
   WRITE(*,*) "UNKNOWN TRAJECTORY INTERPOLATION ORDER : ",C_SLINTP
   STOP
END IF

DO II=1,GXPT
   DO JJ=1,GYPT
      RM_DEP(II,JJ) = ZERO
      RQ_DEP(II,JJ) = ZERO
      DO I=ISTA,IEND
         DO J=ISTA,IEND
            RM_DEP(II,JJ) = RM_DEP(II,JJ)               &
              & + PXWEI(I,II,JJ)*PYWEI(J,II,JJ)*(       &
              &   RM_ARR(KXLAG(I,II,JJ),KYLAG(J,II,JJ)) )
            RQ_DEP(II,JJ) = RQ_DEP(II,JJ)               &
              & + PXWEI(I,II,JJ)*PYWEI(J,II,JJ)*(       &
              &   RQ_ARR(KXLAG(I,II,JJ),KYLAG(J,II,JJ)) )
         ENDDO
      ENDDO         
   ENDDO
ENDDO

END SUBROUTINE LARCINB 
!*************************************************************************************!
!*************************************************************************************!
!*************************************************************************************!
!*************************************************************************************!
!*************************************************************************************!
!*************************************************************************************!
!*************************************************************************************!
SUBROUTINE COMPUTE_DEFORMATIONAL(PDXDOT_DX,PDXDOT_DY,PDYDOT_DX,PDYDOT_DY,PXDOT,PYDOT)
    
IMPLICIT NONE
  

REAL(8),           INTENT(IN)    :: PXDOT(GXPT,GYPT),      PYDOT(GXPT,GYPT)
REAL(8),           INTENT(OUT)   :: PDXDOT_DX(GXPT,GYPT),  PDYDOT_DY(GXPT,GYPT)
REAL(8),           INTENT(OUT)   :: PDXDOT_DY(GXPT,GYPT),  PDYDOT_DX(GXPT,GYPT)

INTEGER(8)                       :: II         ,JJ
REAL(8)                          :: ZVP(GXPT,GYPT)
REAL(8)                          :: ZDISCR,  ZEPS,  ZLIPCHITZ, ZTAU


ZEPS = ABS(EPSILON(ONE))
ZTAU = ONE-ZEPS

!----------------------------------
! COMPUTE DEFORMATIONAL TERM
!----------------------------------
DO II=1,GXPT
   DO JJ=1,GYPT
      PDXDOT_DX(II,JJ) = 0.5d0*( PXDOT(NPERIOX(II+1),JJ)   &
                       & - PXDOT(NPERIOX(II-1),JJ) )/RDX
      PDYDOT_DX(II,JJ) = 0.5d0*( PYDOT(NPERIOX(II+1),JJ)   &
                       & - PYDOT(NPERIOX(II-1),JJ) )/RDX
      PDXDOT_DY(II,JJ) = 0.5d0*( PXDOT(II,NPERIOY(JJ+1))   &
                       & - PXDOT(II,NPERIOY(JJ-1)) )/RDY
      PDYDOT_DY(II,JJ) = 0.5d0*( PYDOT(II,NPERIOY(JJ+1))   &
                       & - PYDOT(II,NPERIOY(JJ-1)) )/RDY
   ENDDO
ENDDO   

!----------------------------------
! COMPUTE LIPCHITZ NUMBERS
!----------------------------------
DO II=1,GXPT
   DO JJ=1,GYPT
      ZDISCR        = ((PDXDOT_DX(II,JJ)-PDYDOT_DY(II,JJ))**2) &
                    & +4.d0*PDXDOT_DY(II,JJ)*PDYDOT_DX(II,JJ)
      IF (ZDISCR .GE. ZEPS) THEN
         ZVP(II,JJ) = 0.25d0*RDT*MAX(ABS((PDXDOT_DX(II,JJ)     &
                    & +PDYDOT_DY(II,JJ))+SQRT(ZDISCR)),        &
                    & ABS((PDXDOT_DX(II,JJ)+PDYDOT_DY(II,JJ))  &
                    & -SQRT(ZDISCR)))
      ELSE
         ZVP(II,JJ) = 0.25d0*RDT*SQRT(((PDXDOT_DX(II,JJ)       &
                    & +PDYDOT_DY(II,JJ))**2)+ABS(ZDISCR))
      END IF
   ENDDO
ENDDO 

!***************************************************************************!
IF (LPRINTLEV) THEN                                                         !
   IF (MAXVAL(ZVP) .GE. ZTAU) THEN                                          !
      WRITE(*,*) 'VIOLATION OF THE LIPCHITZ CRITERIUM AT POINTS : '         !       
      !DO II=1,GXPT                                                         !
      !   DO JJ=1,GYPT                                                      !
      !      ZLIPCHITZ=ZVP(II,JJ)                                           !
      !      IF (ZLIPCHITZ .GE. ZTAU) THEN                                  !
      !         WRITE (*,*) '(I,J) = ','(',II,',',JJ,')',' .... ',ZLIPCHITZ !
      !      END IF                                                         !
      !   END DO                                                            !
      !END DO                                                               !
   END IF                                                                   !             
END IF                                                                      !
!***************************************************************************!


END SUBROUTINE COMPUTE_DEFORMATIONAL
!*************************************************************************************************!
!*************************************************************************************************!
!*************************************************************************************************!
!*************************************************************************************************!
!*************************************************************************************************!
SUBROUTINE COMPUTE_COMAD_STRETCHING(PALPHA_DX,PALPHA_DY,PXDOT,PYDOT,PDXDOT_DX,PDYDOT_DY, &
           & PXDOT_F,PYDOT_F,PDXDOT_DX_F,PDYDOT_DY_F,KXLAG,PXWEI,KYLAG,PYWEI)
    
IMPLICIT NONE

REAL(8),           INTENT(IN)    :: PXDOT(GXPT,GYPT),       PYDOT(GXPT,GYPT)
REAL(8),           INTENT(IN)    :: PXDOT_F(GXPT,GYPT),     PYDOT_F(GXPT,GYPT)
REAL(8),           INTENT(IN)    :: PDXDOT_DX(GXPT,GYPT),   PDYDOT_DY(GXPT,GYPT)
REAL(8),           INTENT(IN)    :: PDXDOT_DX_F(GXPT,GYPT), PDYDOT_DY_F(GXPT,GYPT)
REAL(8),           INTENT(IN)    :: PXWEI(6,GXPT,GYPT),     PYWEI(6,GXPT,GYPT)
INTEGER(8),        INTENT(IN)    :: KXLAG(6,GXPT,GYPT),     KYLAG(6,GXPT,GYPT)
REAL(8),           INTENT(OUT)   :: PALPHA_DX(GXPT,GYPT),   PALPHA_DY(GXPT,GYPT)

INTEGER(8)                       :: II,  JJ,  KK,  LL
INTEGER                          :: ISTA,     IEND
REAL(8)                          :: ZDIVY,    ZEPS,      ZDIVX,   ZSDT
REAL(8)                          :: ZDXDOT_DX(GXPT,GYPT),   ZDYDOT_DY(GXPT,GYPT)
REAL(8)                          :: ZDXDOT_DX_F(GXPT,GYPT), ZDYDOT_DY_F(GXPT,GYPT)
REAL(8)                          :: ZDXDOT_DX_O(GXPT,GYPT), ZDYDOT_DY_O(GXPT,GYPT)

ZEPS = ABS(1000*EPSILON(ONE))
ZSDT = (TWO/RDT)*(ONE-ZEPS)


IF (NCOMAD_OPT == 0) THEN
   DO II=1,GXPT
      DO JJ=1,GYPT
         PALPHA_DX(II,JJ) = MAX(0.0001d0,MIN((ONE+RDT*PDXDOT_DX(II,JJ)),ONE))
         PALPHA_DY(II,JJ) = MAX(0.0001d0,MIN((ONE+RDT*PDYDOT_DY(II,JJ)),ONE))
      ENDDO     
   ENDDO
ELSE
   IF (NCOMAD_OPT == 1) THEN
      ZDXDOT_DX(:,:)   = PDXDOT_DX(:,:)
      ZDYDOT_DY(:,:)   = PDYDOT_DY(:,:)
      ZDXDOT_DX_F(:,:) = PDXDOT_DX_F(:,:)
      ZDYDOT_DY_F(:,:) = PDYDOT_DY_F(:,:)
   ELSE IF (NCOMAD_OPT == 2) THEN   
     !----------------------------------
     ! COMPUTE DEFORMATIONAL TERMS
     DO II=1,GXPT
        DO JJ=1,GYPT
           ZDXDOT_DX(II,JJ)   = 0.5d0*( PXDOT(NPERIOX(II+1),JJ)    &
                              & - PXDOT(NPERIOX(II-1),JJ) )/RDX
           ZDYDOT_DY(II,JJ)   =  0.5d0*( PYDOT(II,NPERIOY(JJ+1))   &
                              & - PYDOT(II,NPERIOY(JJ-1)) )/RDY
           ZDXDOT_DX_F(II,JJ) = 0.5d0*( PXDOT_F(NPERIOX(II+1),JJ)  &
                              & - PXDOT_F(NPERIOX(II-1),JJ) )/RDX
           ZDYDOT_DY_F(II,JJ) =  0.5d0*( PYDOT_F(II,NPERIOY(JJ+1)) &
                              & - PYDOT_F(II,NPERIOY(JJ-1)) )/RDY
        ENDDO
     ENDDO
   END IF  
   !----------------------------------
   ! INTERPOLATED DEFORMATIONAL
   IF (C_SLTRAJ == "C") THEN
      ISTA = 1
      IEND = 4
   ELSE IF (C_SLTRAJ == "L") THEN
      ISTA = 5
      IEND = 6
   END IF

   DO II=1,GYPT
      DO JJ=1,GYPT
         ZDXDOT_DX_O(II,JJ) = ZERO
         ZDYDOT_DY_O(II,JJ) = ZERO
         DO KK=ISTA,IEND
            DO LL=ISTA,IEND
               ZDXDOT_DX_O(II,JJ) = ZDXDOT_DX_O(II,JJ)                          &
                                  & + PXWEI(LL,II,JJ)*PYWEI(KK,II,JJ)*(         &
                                  &   ZDXDOT_DX(KXLAG(LL,II,JJ),KYLAG(KK,II,JJ)))
               ZDYDOT_DY_O(II,JJ) = ZDYDOT_DY_O(II,JJ)                          &
                                  & + PXWEI(LL,II,JJ)*PYWEI(KK,II,JJ)*(         &
                                  &   ZDYDOT_DY(KXLAG(LL,II,JJ),KYLAG(KK,II,JJ)))
            ENDDO
         ENDDO
      ENDDO
   ENDDO
   !----------------------------------
   ! COMPUTE STRETCHING COEFFICIENTS
   DO II=1,GXPT
      DO JJ=1,GYPT
         ZDIVX = ZDXDOT_DX_O(II,JJ)
         IF (ZDIVX .LE. ZSDT) THEN
            PALPHA_DX(II,JJ) = MAX(0.0001d0,MIN((ONE+(RDT/TWO)*ZDXDOT_DX_F(II,JJ)) &
                             & /(ONE-(RDT/TWO)*ZDXDOT_DX_O(II,JJ)),ONE))
         ELSE
            PALPHA_DX(II,JJ) = ONE
         END IF    
         ZDIVY = ZDYDOT_DY_O(II,JJ)
         IF (ZDIVY .LE. ZSDT) THEN
            PALPHA_DY(II,JJ) = MAX(0.0001d0,MIN((ONE+(RDT/TWO)*ZDYDOT_DY_F(II,JJ)) &
                             & /(ONE-(RDT/TWO)*ZDYDOT_DY_O(II,JJ)),ONE))
         ELSE
            PALPHA_DY(II,JJ) = ONE
         END IF
      ENDDO     
   ENDDO       
END IF

END SUBROUTINE COMPUTE_COMAD_STRETCHING
!*************************************************************************************************!
!*************************************************************************************************!
!*************************************************************************************************!
!*************************************************************************************************!
!*************************************************************************************************!
!*************************************************************************************************!
!*************************************************************************************************!
SUBROUTINE ILMC_FILTER(PDEP,PARR,KXLAG,KYLAG,PM)
  
IMPLICIT NONE

REAL(8),           INTENT(INOUT) :: PDEP(GXPT,GYPT)
REAL(8),           INTENT(IN)    :: PARR(GXPT,GYPT)
INTEGER(8),        INTENT(IN)    :: KXLAG(6,GXPT,GYPT),KYLAG(6,GXPT,GYPT)
REAL(8), OPTIONAL, INTENT(IN)    :: PM(GXPT,GYPT)

INTEGER(8)                       :: I,J
REAL(8)                          :: ZM(48)
REAL(8)                          :: ZPM(GXPT,GYPT)
REAL(8)                          :: ZMAX_LOC(GXPT,GYPT),ZMIN_LOC(GXPT,GYPT)


IF (PRESENT(PM)) THEN
   ZPM(:,:) = PM(:,:)
ELSE
   ZPM(:,:) = ONE
END IF   

DO I=1,GXPT
   DO J=1,GYPT
      ZMAX_LOC(I,J) =  MAX(                  &
          & PARR(KXLAG(2,I,J),KYLAG(2,I,J)), &
          & PARR(KXLAG(2,I,J),KYLAG(3,I,J)), &
          & PARR(KXLAG(3,I,J),KYLAG(2,I,J)), &
          & PARR(KXLAG(3,I,J),KYLAG(3,I,J))  )
      ZMIN_LOC(I,J) =  MIN(                  &
          & PARR(KXLAG(2,I,J),KYLAG(2,I,J)), &
          & PARR(KXLAG(2,I,J),KYLAG(3,I,J)), &
          & PARR(KXLAG(3,I,J),KYLAG(2,I,J)), &
          & PARR(KXLAG(3,I,J),KYLAG(3,I,J))  )
   ENDDO
ENDDO

!_____._____._____._____._____._____._____.
!     !     !     !     !     !     !     !
! 42  !  43 !  44 !  45 !  46 !  47 !  48 ! 
!_____!_____!_____!_____!_____!_____!_____!
!     !     !     !     !     !     !     !
! 40  !  20 !  21 !  22 !  23 !  24 !  41 !
!_____!_____!_____!_____!_____!_____!_____!
!     !     !     !     !     !     !     !
! 38  !  18 !  6  !  7  !  8  !  19 !  39 !     
!_____!_____!_____!_____!_____!_____!_____!
!     !     !     !     !     !     !     !
! 36  !  16 !  4  !  X  !  5  !  17 !  37 !
!_____!_____!_____!_____!_____!_____!_____!
!     !     !     !     !     !     !     !
! 34  !  14 !  1  !  2  !  3  !  15 !  35 !
!_____!_____!_____!_____!_____!_____!_____!
!     !     !     !     !     !     !     !
! 32  !  9  !  10 !  11 !  12 !  13 !  33 !
!_____!_____!_____!_____!_____!_____!_____!
!     !     !     !     !     !     !     !
! 25  !  26 !  27 !  28 !  29 !  30 !  31 !
!_____!_____!_____!_____!_____!_____!_____!

DO I=1,GXPT
   DO J=1,GYPT
      !----------------------------!
      !  DEALING WITH OVERSHOOTS   !
      !----------------------------!
      IF (PDEP(I,J) > ZMAX_LOC(I,J)) THEN
         ZM(0) = (PDEP(I,J)-ZMAX_LOC(I,J))*ZPM(I,J)
         ZM(1) = MAX(ZMAX_LOC(NPERIOX(I-1),NPERIOY(J-1))-PDEP(NPERIOX(I-1),NPERIOY(J-1)),ZERO)*( &
               & ZPM(NPERIOX(I-1),NPERIOY(J-1)) )
         ZM(2) = MAX(ZMAX_LOC(I,NPERIOY(J-1))-PDEP(I,NPERIOY(J-1)),ZERO)*(                       &
               & ZPM(I,NPERIOY(J-1)))
         ZM(3) = MAX(ZMAX_LOC(NPERIOX(I+1),NPERIOY(J-1))-PDEP(NPERIOX(I+1),NPERIOY(J-1)),ZERO)*( &
               & ZPM(NPERIOX(I+1),NPERIOY(J-1)) )
         ZM(4) = MAX(ZMAX_LOC(NPERIOX(I-1),J)-PDEP(NPERIOX(I-1),J),ZERO)*(                       &
               & ZPM(NPERIOX(I-1),J))
         ZM(5) = MAX(ZMAX_LOC(NPERIOX(I+1),J)-PDEP(NPERIOX(I+1),J),ZERO)*(                       &
               & ZPM(NPERIOX(I+1),J))
         ZM(6) = MAX(ZMAX_LOC(NPERIOX(I-1),NPERIOY(J+1))-PDEP(NPERIOX(I-1),NPERIOY(J+1)),ZERO)*( &
               & ZPM(NPERIOX(I-1),NPERIOY(J+1))
         ZM(7) = MAX(ZMAX_LOC(I,NPERIOY(J+1))-PDEP(I,NPERIOY(J+1)),ZERO)*(                       &
               & ZPM(I,NPERIOY(J+1)))
         ZM(8) = MAX(ZMAX_LOC(NPERIOX(I+1),NPERIOY(J+1))-PDEP(NPERIOX(I+1),NPERIOY(J+1)),ZERO)*( &
               & ZPM(NPERIOX(I+1),NPERIOY(J+1)) )
         ZM_A  = ZERO
         DO K=1,8
            ZM_A = ZM_A + ZM(K)
         END DO
         IF (ZM_A > ZM(0)) THEN
            PDEP(I,J) = ZMAX_LOC(I,J)
            PDEP(NPERIOX(I-1),NPERIOY(J-1)) = PDEP(NPERIOX(I-1),NPERIOY(J-1)) &
                                            & + ZM(0)*(ZM(1)/ZM_A)/ZPM(NPERIOX(I-1),NPERIOY(J-1))
            PDEP(I,NPERIOY(J-1))            = PDEP(I,NPERIOY(J-1))            &
                                            & + ZM(0)*(ZM(2)/ZM_A)/ZPM(I,NPERIOY(J-1))
            PDEP(NPERIOX(I+1),NPERIOY(J-1)) = PDEP(NPERIOX(I+1),NPERIOY(J-1)) &
                                            & + ZM(0)*(ZM(3)/ZM_A)/ZPM(NPERIOX(I+1),NPERIOY(J-1))
            PDEP(NPERIOX(I-1),J)            = PDEP(NPERIOX(I-1),J)            &
                                            & + ZM(0)*(ZM(4)/ZM_A)/ZPM(NPERIOX(I-1),J) 
            PDEP(NPERIOX(I+1),J)            = PDEP(NPERIOX(I+1),J)            &
                                            & + ZM(0)*(ZM(5)/ZM_A)/ZPM(NPERIOX(I+1),J)
            PDEP(NPERIOX(I-1),NPERIOY(J+1)) = PDEP(NPERIOX(I-1),NPERIOY(J+1)) &
                                            & + ZM(0)*(ZM(6)/ZM_A)/ZPM(NPERIOX(I-1),NPERIOY(J+1))
            PDEP(I,NPERIOY(J+1))            = PDEP(I,NPERIOY(J+1))            &
                                            & + ZM(0)*(ZM(7)/ZM_A)/ZPM(I,NPERIOY(J+1))
            PDEP(NPERIOX(I+1),NPERIOY(J+1)) = PDEP(NPERIOX(I+1),NPERIOY(J+1)) &
                                            & + ZM(0)*(ZM(8)/ZM_A)/ZPM(NPERIOX(I+1),NPERIOY(J+1))
         ELSE
            ZM(9)  = MAX(ZMAX_LOC(NPERIOX(I-2),NPERIOY(J-2))-PDEP(NPERIOX(I-2),NPERIOY(J-2)),ZERO)*( &
                   & ZPM(NPERIOX(I-2),NPERIOY(J-2)))
            ZM(10) = MAX(ZMAX_LOC(NPERIOX(I-1),NPERIOY(J-2))-PDEP(NPERIOX(I-1),NPERIOY(J-2)),ZERO)*( &
                   & ZPM(NPERIOX(I-1),NPERIOY(J-2)))
            ZM(11) = MAX(ZMAX_LOC(I,NPERIOY(J-2))-PDEP(I,NPERIOY(J-2)),ZERO)*(                       &
                   & ZPM(I,NPERIOY(J-2))
            ZM(12) = MAX(ZMAX_LOC(NPERIOX(I+1),NPERIOY(J-2))-PDEP(NPERIOX(I+1),NPERIOY(J-2)),ZERO)*( &
                   & ZPM(NPERIOX(I+1),NPERIOY(J-2)))
            ZM(13) = MAX(ZMAX_LOC(NPERIOX(I+2),NPERIOY(J-2))-PDEP(NPERIOX(I+2),NPERIOY(J-2)),ZERO)*( &
                   & ZPM(NPERIOX(I+2),NPERIOY(J-2)))
            ZM(14) = MAX(ZMAX_LOC(NPERIOX(I-2),NPERIOY(J-1))-PDEP(NPERIOX(I-2),NPERIOY(J-1)),ZERO)*( &
                   & ZPM(NPERIOX(I-2),NPERIOY(J-1)))
            ZM(15) = MAX(ZMAX_LOC(NPERIOX(I+2),NPERIOY(J-1))-PDEP(NPERIOX(I+2),NPERIOY(J-1)),ZERO)*( &
                   & ZPM(NPERIOX(I+2),NPERIOY(J-1)))
            ZM(16) = MAX(ZMAX_LOC(NPERIOX(I-2),J)-PDEP(NPERIOX(I-2),J),ZERO)*(                       &
                   & ZPM(NPERIOX(I-2),J))
            ZM(17) = MAX(ZMAX_LOC(NPERIOX(I+2),J)-PDEP(NPERIOX(I+2),J),ZERO)*(                       &
                   & ZPM(NPERIOX(I+2),J))
            ZM(18) = MAX(ZMAX_LOC(NPERIOX(I-2),NPERIOY(J+1))-PDEP(NPERIOX(I-2),NPERIOY(J+1)),ZERO)*( &
                   & ZPM(NPERIOX(I-2),NPERIOY(J+1)))
            ZM(19) = MAX(ZMAX_LOC(NPERIOX(I+2),NPERIOY(J+1))-PDEP(NPERIOX(I+2),NPERIOY(J+1)),ZERO)*( &
                   & ZPM(NPERIOX(I+2),NPERIOY(J+1)))
            ZM(20) = MAX(ZMAX_LOC(NPERIOX(I-2),NPERIOY(J+2))-PDEP(NPERIOX(I-2),NPERIOY(J+2)),ZERO)*( &
                   & ZPM(NPERIOX(I-2),NPERIOY(J+2)))
            ZM(21) = MAX(ZMAX_LOC(NPERIOX(I-1),NPERIOY(J+2))-PDEP(NPERIOX(I-1),NPERIOY(J+2)),ZERO)*( &
                   & ZPM(NPERIOX(I-1),NPERIOY(J+2)))
            ZM(22) = MAX(ZMAX_LOC(I,NPERIOY(J+2))-PDEP(I,NPERIOY(J+2)),ZERO)*(                       &
                   & ZPM(I,NPERIOY(J+2)))
            ZM(23) = MAX(ZMAX_LOC(NPERIOX(I+1),NPERIOY(J+2))-PDEP(NPERIOX(I+1),NPERIOY(J+2)),ZERO)*( &
                   & ZPM(NPERIOX(I+1),NPERIOY(J+2)))
            ZM(24) = MAX(ZMAX_LOC(NPERIOX(I+2),NPERIOY(J+2))-PDEP(NPERIOX(I+2),NPERIOY(J+2)),ZERO)*( &
                   & ZPM(NPERIOX(I+2),NPERIOY(J+2)))
            DO K=9,24
               ZM_A = ZM_A + ZM(K)
            END DO
            IF (ZM_A > ZM(0)) THEN
               PDEP(I,J) = ZMAX_LOC(I,J)
               PDEP(NPERIOX(I-1),NPERIOY(J-1)) = PDEP(NPERIOX(I-1),NPERIOY(J-1)) &
                                               & + ZM(0)*(ZM(1)/ZM_A)/ZPM(NPERIOX(I-1),NPERIOY(J-1))
               PDEP(I,NPERIOY(J-1))            = PDEP(I,NPERIOY(J-1))            &
                                               & + ZM(0)*(ZM(2)/ZM_A)/ZPM(I,NPERIOY(J-1))
               PDEP(NPERIOX(I+1),NPERIOY(J-1)) = PDEP(NPERIOX(I+1),NPERIOY(J-1)) &
                                               & + ZM(0)*(ZM(3)/ZM_A)/ZPM(NPERIOX(I+1),NPERIOY(J-1))
               PDEP(NPERIOX(I-1),J)            = PDEP(NPERIOX(I-1),J)            &
                                               & + ZM(0)*(ZM(4)/ZM_A)/ZPM(NPERIOX(I-1),J) 
               PDEP(NPERIOX(I+1),J)            = PDEP(NPERIOX(I+1),J)            &
                                               & + ZM(0)*(ZM(5)/ZM_A)/ZPM(NPERIOX(I+1),J)
               PDEP(NPERIOX(I-1),NPERIOY(J+1)) = PDEP(NPERIOX(I-1),NPERIOY(J+1)) &
                                               & + ZM(0)*(ZM(6)/ZM_A)/ZPM(NPERIOX(I-1),NPERIOY(J+1))
               PDEP(I,NPERIOY(J+1))            = PDEP(I,NPERIOY(J+1))            &          
                                               & + ZM(0)*(ZM(7)/ZM_A)/ZPM(I,NPERIOY(J+1))
               PDEP(NPERIOX(I+1),NPERIOY(J+1)) = PDEP(NPERIOX(I+1),NPERIOY(J+1)) &
                                               & + ZM(0)*(ZM(8)/ZM_A)/ZPM(NPERIOX(I+1),NPERIOY(J+1))
               PDEP(NPERIOX(I-2),NPERIOY(J-2)) = PDEP(NPERIOX(I-2),NPERIOY(J-2)) &
                                               & + ZM(0)*(ZM(9)/ZM_A)/ZPM(NPERIOX(I-2),NPERIOY(J-2))
               PDEP(NPERIOX(I-1),NPERIOY(J-2)) = PDEP(NPERIOX(I-1),NPERIOY(J-2)) &
                                               & + ZM(0)*(ZM(10)/ZM_A)/ZPM(NPERIOX(I-1),NPERIOY(J-2))
               PDEP(I,NPERIOY(J-2))            = PDEP(I,NPERIOY(J-2))            &
                                               & + ZM(0)*(ZM(11)/ZM_A)/ZPM(I,NPERIOY(J-2))
               PDEP(NPERIOX(I+1),NPERIOY(J-2)) = PDEP(NPERIOX(I+1),NPERIOY(J-2)) &
                                               & + ZM(0)*(ZM(12)/ZM_A)/ZPM(NPERIOX(I+1),NPERIOY(J-2))
               PDEP(NPERIOX(I+2),NPERIOY(J-2)) = PDEP(NPERIOX(I+2),NPERIOY(J-2)) &
                                               & + ZM(0)*(ZM(13)/ZM_A)/ZPM(NPERIOX(I+2),NPERIOY(J-2))
               PDEP(NPERIOX(I-2),NPERIOY(J-1)) = PDEP(NPERIOX(I-2),NPERIO(J-1))  &
                                               & + ZM(0)*(ZM(14)/ZM_A)/ZPM(NPERIOX(I-2),NPERIO(J-1))
               PDEP(NPERIOX(I+2),NPERIOY(J-1)) = PDEP(NPERIOX(I+2),NPERIOY(J-1)) &
                                               & + ZM(0)*(ZM(15)/ZM_A)/ZPM(NPERIOX(I+2),NPERIOY(J-1))
               PDEP(NPERIOX(I-2),J)            = PDEP(NPERIOX(I-2),J)            &
                                               & + ZM(0)*(ZM(16)/ZM_A)/ZPM(NPERIOX(I-2),J)
               PDEP(NPERIOX(I+2),J)            = PDEP(NPERIOX(I+2),J)            &
                                               & + ZM(0)*(ZM(17)/ZM_A)/ZPM(NPERIOX(I+2),J)
               PDEP(NPERIOX(I-2),NPERIOY(J+1)) = PDEP(NPERIOX(I-2),NPERIOY(J+1)) &
                                               & + ZM(0)*(ZM(18)/ZM_A)/ZPM(NPERIOX(I-2),NPERIOY(J+1))
               PDEP(NPERIOX(I+2),NPERIOY(J+1)) = PDEP(NPERIOX(I+2),NPERIOY(J+1)) &
                                               & + ZM(0)*(ZM(19)/ZM_A)/ZPM(NPERIOX(I+2),NPERIOY(J+1))
               PDEP(NPERIOX(I-2),NPERIOY(J+2)) = PDEP(NPERIOX(I-2),NPERIOY(J+2)) &
                                               & + ZM(0)*(ZM(20)/ZM_A)/ZPM(NPERIOX(I-2),NPERIOY(J+2))
               PDEP(NPERIOX(I-1),NPERIOY(J+2)) = PDEP(NPERIOX(I-1),NPERIOY(J+2)) &
                                               & + ZM(0)*(ZM(21)/ZM_A)/ZPM(NPERIOX(I-1),NPERIOY(J+2))
               PDEP(I,NPERIOY(J+2))            = PDEP(I,NPERIOY(J+2))            &
                                               & + ZM(0)*(ZM(22)/ZM_A)/PM(I,NPERIOY(J+2))
               PDEP(NPERIOX(I+1),NPERIOY(J+2)) = PDEP(NPERIOX(I+1),NPERIOY(J+2)) &
                                               & + ZM(0)*(ZM(23)/ZM_A)/ZPM(NPERIOX(I+1),NPERIOY(J+2))
               PDEP(NPERIOX(I+2),NPERIOY(J+2)) = PDEP(NPERIOX(I+2),NPERIOY(J+2)) &
                                               & + ZM(0)*(ZM(24)/ZM_A)/ZPM(NPERIOX(I+2),NPERIOY(J+2))
            ELSE
               ZM(25) = MAX(ZMAX_LOC(NPERIOX(I-3),NPERIOY(J-3))-PDEP(NPERIOX(I-3),NPERIOY(J-3)),ZERO)*( &
                      & ZPM(NPERIOX(I-3),NPERIOY(J-3)) )
               ZM(26) = MAX(ZMAX_LOC(NPERIOX(I-2),NPERIOY(J-3))-PDEP(NPERIOX(I-2),NPERIOY(J-3)),ZERO)*( &
                      & ZPM(NPERIOX(I-2),NPERIOY(J-3)))
               ZM(27) = MAX(ZMAX_LOC(NPERIOX(I-1),NPERIOY(J-3))-PDEP(NPERIOX(I-1),NPERIOY(J-3)),ZERO)*( &
                      & ZPM(NPERIOX(I-1),NPERIOY(J-3)) )
               ZM(28) = MAX(ZMAX_LOC(NPERIOX(I),NPERIOY(J-3))-PDEP(NPERIOX(I),NPERIOY(J-3)),ZERO)*(     &
                      & ZPM(NPERIOX(I),NPERIOY(J-3)))
               ZM(29) = MAX(ZMAX_LOC(NPERIOX(I+1),NPERIOY(J-3))-PDEP(NPERIOX(I+1),NPERIOY(J-3)),ZERO)*( &
                      & ZPM(NPERIOX(I+1),NPERIOY(J-3)))
               ZM(30) = MAX(ZMAX_LOC(NPERIOX(I+2),NPERIOY(J-3))-PDEP(NPERIOX(I+2),NPERIOY(J-3)),ZERO)*( &
                      & ZPM(NPERIOX(I+2),NPERIOY(J-3))
               ZM(31) = MAX(ZMAX_LOC(NPERIOX(I+3),NPERIOY(J-3))-PDEP(NPERIOX(I+3),NPERIOY(J-3)),ZERO)*( &
                      & ZPM(NPERIOX(I+3),NPERIOY(J-3)))
               ZM(32) = MAX(ZMAX_LOC(NPERIOX(I-3),NPERIOY(J-2))-PDEP(NPERIOX(I-3),NPERIOY(J-2)),ZERO)*( &
                      & ZPM(NPERIOX(I-3),NPERIOY(J-2)) )
               ZM(33) = MAX(ZMAX_LOC(NPERIOX(I+3),NPERIOY(J-2))-PDEP(NPERIOX(I+3),NPERIOY(J-2)),ZERO)*( &
                      & ZPM(NPERIOX(I+3),NPERIOY(J-2)))
               ZM(34) = MAX(ZMAX_LOC(NPERIOX(I-3),NPERIOY(J-1))-PDEP(NPERIOX(I-3),NPERIOY(J-1)),ZERO)*( &
                      & ZPM(NPERIOX(I-3),NPERIOY(J-1)))
               ZM(35) = MAX(ZMAX_LOC(NPERIOX(I+3),NPERIOY(J-1))-PDEP(NPERIOX(I+3),NPERIOY(J-1)),ZERO)*( &
                      & ZPM(NPERIOX(I+3),NPERIOY(J-1)))
               ZM(36) = MAX(ZMAX_LOC(NPERIOX(I-3),NPERIOY(J))-PDEP(NPERIOX(I-3),NPERIOY(J)),ZERO)*(     &
                      & ZPM(NPERIOX(I-3),NPERIOY(J)))
               ZM(37) = MAX(ZMAX_LOC(NPERIOX(I+3),NPERIOY(J))-PDEP(NPERIOX(I+3),NPERIOY(J)),ZERO)*(     &
                      & ZPM(NPERIOX(I+3),NPERIOY(J)))
               ZM(38) = MAX(ZMAX_LOC(NPERIOX(I-3),NPERIOY(J+1))-PDEP(NPERIOX(I-3),NPERIOY(J+1)),ZERO)*( &
                      & ZPM(NPERIOX(I-3),NPERIOY(J+1)))
               ZM(39) = MAX(ZMAX_LOC(NPERIOX(I+3),NPERIOY(J+1))-PDEP(NPERIOX(I+3),NPERIOY(J+1)),ZERO)*( &
                      & ZPM(NPERIOX(I+3),NPERIOY(J+1)))
               ZM(40) = MAX(ZMAX_LOC(NPERIOX(I-3),NPERIOY(J+2))-PDEP(NPERIOX(I-3),NPERIOY(J+2)),ZERO)*( &
                      & ZPM(NPERIOX(I-3),NPERIOY(J+2)))
               ZM(41) = MAX(ZMAX_LOC(NPERIOX(I+3),NPERIOY(J+2))-PDEP(NPERIOX(I+3),NPERIOY(J+2)),ZERO)*( &
                      & ZPM(NPERIOX(I+3),NPERIOY(J+2)))
               ZM(42) = MAX(ZMAX_LOC(NPERIOX(I-3),NPERIOY(J+3))-PDEP(NPERIOX(I-3),NPERIOY(J+3)),ZERO)*( &
                      & ZPM(NPERIOX(I-3),NPERIOY(J+3)))
               ZM(43) = MAX(ZMAX_LOC(NPERIOX(I-2),NPERIOY(J+3))-PDEP(NPERIOX(I-2),NPERIOY(J+3)),ZERO)*( &
                      & ZPM(NPERIOX(I-2),NPERIOY(J+3)))
               ZM(44) = MAX(ZMAX_LOC(NPERIOX(I-1),NPERIOY(J+3))-PDEP(NPERIOX(I-1),NPERIOY(J+3)),ZERO)*( &
                      & ZPM(NPERIOX(I-1),NPERIOY(J+2)))
               ZM(45) = MAX(ZMAX_LOC(NPERIOX(I),NPERIOY(J+3))-PDEP(NPERIOX(I),NPERIOY(J+3)),ZERO)*(     &
                      & ZPM(NPERIOX(I),NPERIOY(J+3)))
               ZM(46) = MAX(ZMAX_LOC(NPERIOX(I+1),NPERIOY(J+3))-PDEP(NPERIOX(I+1),NPERIOY(J+3)),ZERO)*( &
                      & ZPM(NPERIOX(I+1),NPERIOY(J+3)))
               ZM(47) = MAX(ZMAX_LOC(NPERIOX(I+2),NPERIOY(J+3))-PDEP(NPERIOX(I+2),NPERIOY(J+3),ZERO)*(  &
                      & ZPM(NPERIOX(I+2),NPERIOY(J+3)))
               ZM(48) = MAX(ZMAX_LOC(NPERIOX(I+3),NPERIOY(J+3))-PDEP(NPERIOX(I+3),NPERIOY(J+3)),ZERO)*( &
                      & ZPM(NPERIOX(I+3),NPERIOY(J+3)))
               DO K=25,48
                  ZM_A = ZM_A + ZM(K)
               END DO
               IF (ZM_A > ZM(0)) THEN
                  PDEP(I,J) = ZMAX_LOC(I,J)
                  PDEP(NPERIOX(I-1),NPERIOY(J-1)) = PDEP(NPERIOX(I-1),NPERIOY(J-1)) &
                                                  & + ZM(0)*(ZM(1)/ZM_A)/ZPM(NPERIOX(I-1),NPERIOY(J-1))
                  PDEP(I,NPERIOY(J-1))            = PDEP(I,NPERIOY(J-1))            &
                                                  & + ZM(0)*(ZM(2)/ZM_A)/ZPM(I,NPERIOY(J-1))
                  PDEP(NPERIOX(I+1),NPERIOY(J-1)) = PDEP(NPERIOX(I+1),NPERIOY(J-1)) &
                                                  & + ZM(0)*(ZM(3)/ZM_A)/ZPM(NPERIOX(I+1),NPERIOY(J-1))
                  PDEP(NPERIOX(I-1),J)            = PDEP(NPERIOX(I-1),J)            &
                                                  & + ZM(0)*(ZM(4)/ZM_A)/ZPM(NPERIOX(I-1),J) 
                  PDEP(NPERIOX(I+1),J)            = PDEP(NPERIOX(I+1),J)            &
                                                  & + ZM(0)*(ZM(5)/ZM_A)/ZPM(NPERIOX(I+1),J)
                  PDEP(NPERIOX(I-1),NPERIOY(J+1)) = PDEP(NPERIOX(I-1),NPERIOY(J+1)) &
                                                  & + ZM(0)*(ZM(6)/ZM_A)/ZPM(NPERIOX(I-1),NPERIOY(J+1))
                  PDEP(I,NPERIOY(J+1))            = PDEP(I,NPERIOY(J+1))            &
                                                  & + ZM(0)*(ZM(7)/ZM_A)/ZPM(I,NPERIOY(J+1))
                  PDEP(NPERIOX(I+1),NPERIOY(J+1)) = PDEP(NPERIOX(I+1),NPERIOY(J+1)) &
                                                  & + ZM(0)*(ZM(8)/ZM_A)/ZPM(NPERIOX(I+1),NPERIOY(J+1))
                  PDEP(NPERIOX(I-2),NPERIOY(J-2)) = PDEP(NPERIOX(I-2),NPERIOY(J-2)) &
                                                  & + ZM(0)*(ZM(9)/ZM_A)/ZPM(NPERIOX(I-2),NPERIOY(J-2))
                  PDEP(NPERIOX(I-1),NPERIOY(J-2)) = PDEP(NPERIOX(I-1),NPERIOY(J-2)) &
                                                  & + ZM(0)*(ZM(10)/ZM_A)/ZPM(NPERIOX(I-1),NPERIOY(J-2))
                  PDEP(I,NPERIOY(J-2))            = PDEP(I,NPERIOY(J-2))            &
                                                  & + ZM(0)*(ZM(11)/ZM_A)/ZPM(I,NPERIOY(J-2))
                  PDEP(NPERIOX(I+1),NPERIOY(J-2)) = PDEP(NPERIOX(I+1),NPERIOY(J-2)) &
                                                  & + ZM(0)*(ZM(12)/ZM_A)/ZPM(NPERIOX(I+1),NPERIOY(J-2))
                  PDEP(NPERIOX(I+2),NPERIOY(J-2)) = PDEP(NPERIOX(I+2),NPERIOY(J-2)) &
                                                  & + ZM(0)*(ZM(13)/ZM_A)/ZPM(NPERIOX(I+2),NPERIOY(J-2))
                  PDEP(NPERIOX(I-2),NPERIOY(J-1)) = PDEP(NPERIOX(I-2),NPERIO(J-1))  &
                                                  & + ZM(0)*(ZM(14)/ZM_A)/ZPM(NPERIOX(I-2),NPERIO(J-1))
                  PDEP(NPERIOX(I+2),NPERIOY(J-1)) = PDEP(NPERIOX(I+2),NPERIOY(J-1)) &
                                                  & + ZM(0)*(ZM(15)/ZM_A)/ZPM(NPERIOX(I+2),NPERIOY(J-1))
                  PDEP(NPERIOX(I-2),J)            = PDEP(NPERIOX(I-2),J)            &
                                                  & + ZM(0)*(ZM(16)/ZM_A)/ZPM(NPERIOX(I-2),J)
                  PDEP(NPERIOX(I+2),J)            = PDEP(NPERIOX(I+2),J)            &
                                                  & + ZM(0)*(ZM(17)/ZM_A)/ZPM(NPERIOX(I+2),J)
                  PDEP(NPERIOX(I-2),NPERIOY(J+1)) = PDEP(NPERIOX(I-2),NPERIOY(J+1)) &
                                                  & + ZM(0)*(ZM(18)/ZM_A)/ZPM(NPERIOX(I-2),NPERIOY(J+1))
                  PDEP(NPERIOX(I+2),NPERIOY(J+1)) = PDEP(NPERIOX(I+2),NPERIOY(J+1)) &
                                                  & + ZM(0)*(ZM(19)/ZM_A)/ZPM(NPERIOX(I+2),NPERIOY(J+1))
                  PDEP(NPERIOX(I-2),NPERIOY(J+2)) = PDEP(NPERIOX(I-2),NPERIOY(J+2)) &
                                                  & + ZM(0)*(ZM(20)/ZM_A)/ZPM(NPERIOX(I-2),NPERIOY(J+2))
                  PDEP(NPERIOX(I-1),NPERIOY(J+2)) = PDEP(NPERIOX(I-1),NPERIOY(J+2)) &
                                                  & + ZM(0)*(ZM(21)/ZM_A)/ZPM(NPERIOX(I-1),NPERIOY(J+2))
                  PDEP(I,NPERIOY(J+2))            = PDEP(I,NPERIOY(J+2))            &
                                                  & + ZM(0)*(ZM(22)/ZM_A)/ZPM(I,NPERIOY(J+2))
                  PDEP(NPERIOX(I+1),NPERIOY(J+2)) = PDEP(NPERIOX(I+1),NPERIOY(J+2)) &
                                                  & + ZM(0)*(ZM(23)/ZM_A)/ZPM(NPERIOX(I+1),NPERIOY(J+2))
                  PDEP(NPERIOX(I+2),NPERIOY(J+2)) = PDEP(NPERIOX(I+2),NPERIOY(J+2)) &
                                                  & + ZM(0)*(ZM(24)/ZM_A)/ZPM(NPERIOX(I+2),NPERIOY(J+2))
                  PDEP(NPERIOX(I-3),NPERIOY(J-3)) = PDEP(NPERIOX(I-3),NPERIOY(J-3)) &
                                                  & + ZM(0)*(ZM(25)/ZM_A)/ZPM(NPERIOX(I-3),NPERIOY(J-3))
                  PDEP(NPERIOX(I-2),NPERIOY(J-3)) = PDEP(NPERIOX(I-2),NPERIOY(J-3)) &
                                                  & + ZM(0)*(ZM(26)/ZM_A)/ZPM(I,NPERIOY(J-1))
                  PDEP(NPERIOX(I-1),NPERIOY(J-3)) = PDEP(NPERIOX(I-1),NPERIOY(J-3)) &
                                                  & + ZM(0)*(ZM(27)/ZM_A)/ZPM(NPERIOX(I-1),NPERIOY(J-3))
                  PDEP(NPERIOX(I),NPERIOY(J-3))   = PDEP(NPERIOX(I),NPERIOY(J-3))   &
                                                  & + ZM(0)*(ZM(28)/ZM_A)/ZPM(NPERIOX(I-1),NPERIOY(J-3)) 
                  PDEP(NPERIOX(I+1),NPERIOY(J-3)) = PDEP(NPERIOX(I+1),NPERIOY(J-3)) &
                                                  & + ZM(0)*(ZM(29)/ZM_A)/ZPM(NPERIOX(I+1),NPERIOY(J-3))
                  PDEP(NPERIOX(I+2),NPERIOY(J-3)) = PDEP(NPERIOX(I+2),NPERIOY(J-3)) &
                                                  & + ZM(0)*(ZM(30)/ZM_A)/ZPM(NPERIOX(I+2),NPERIOY(J-3))
                  PDEP(NPERIOX(I+3),NPERIOY(J-3)) = PDEP(NPERIOX(I+3),NPERIOY(J-3)) &
                                                  & + ZM(0)*(ZM(31)/ZM_A)/ZPM(NPERIOX(I+3),NPERIOY(J-3))
                  PDEP(NPERIOX(I-3),NPERIOY(J-2)) = PDEP(NPERIOX(I-3),NPERIOY(J-2)) &
                                                  & + ZM(0)*(ZM(32)/ZM_A)/ZPM(NPERIOX(I-3),NPERIOY(J-2))
                  PDEP(NPERIOX(I+3),NPERIOY(J-2)) = PDEP(NPERIOX(I+3),NPERIOY(J-2)) &
                                                  & + ZM(0)*(ZM(33)/ZM_A)/ZPM(NPERIOX(I+3),NPERIOY(J-2))
                  PDEP(NPERIOX(I-3),NPERIOY(J-1)) = PDEP(NPERIOX(I-3),NPERIOY(J-1)) &
                                                  & + ZM(0)*(ZM(34)/ZM_A)/ZPM(NPERIOX(I-3),NPERIOY(J-1))
                  PDEP(NPERIOX(I+3),NPERIOY(J-1)) = PDEP(NPERIOX(I+3),NPERIOY(J-1)) &
                                                  & + ZM(0)*(ZM(35)/ZM_A)/ZPM(NPERIOX(I+3),NPERIOY(J-1))
                  PDEP(NPERIOX(I-3),NPERIOY(J))   = PDEP(NPERIOX(I-3),NPERIOY(J))   &
                                                  & + ZM(0)*(ZM(36)/ZM_A)/ZPM(NPERIOX(I-3),NPERIOY(J))
                  PDEP(NPERIOX(I+3),NPERIOY(J))   = PDEP(NPERIOX(I+3),NPERIOY(J))   &
                                                  & + ZM(0)*(ZM(37)/ZM_A)/ZPM(NPERIOX(I+3),NPERIOY(J))
                  PDEP(NPERIOX(I-3),NPERIOY(J+1)) = PDEP(NPERIOX(I-3),NPERIO(J+1))  &
                                                  & + ZM(0)*(ZM(38)/ZM_A)/ZPM(NPERIOX(I-3),NPERIO(J+1))
                  PDEP(NPERIOX(I+3),NPERIOY(J+1)) = PDEP(NPERIOX(I+3),NPERIOY(J+1)) &
                                                  & + ZM(0)*(ZM(39)/ZM_A)/ZPM(NPERIOX(I+3),NPERIOY(J+1))
                  PDEP(NPERIOX(I-3),NPERIOY(J+2)) = PDEP(NPERIOX(I-3),NPERIOY(J+2)) &
                                                  & + ZM(0)*(ZM(40)/ZM_A)/ZPM(NPERIOX(I-3),NPERIOY(J+2)) 
                  PDEP(NPERIOX(I+3),NPERIOY(J+2)) = PDEP(NPERIOX(I+3),NPERIOY(J+2)) &                                     
                                                  & + ZM(0)*(ZM(41)/ZM_A)/ZPM(NPERIOX(I+3),NPERIOY(J+2))
                  PDEP(NPERIOX(I-3),NPERIOY(J+3)) = PDEP(NPERIOX(I-3),NPERIOY(J+3)) &
                                                  & + ZM(0)*(ZM(42)/ZM_A)/ZPM(NPERIOX(I-3),NPERIOY(J+3))
                  PDEP(NPERIOX(I-2),NPERIOY(J+3)) = PDEP(NPERIOX(I-2),NPERIOY(J+3)) &
                                                  & + ZM(0)*(ZM(43)/ZM_A)/ZPM(NPERIOX(I-2),NPERIOY(J+3))
                  PDEP(NPERIOX(I-1),NPERIOY(J+3)) = PDEP(NPERIOX(I-1),NPERIOY(J+3)) &
                                                  & + ZM(0)*(ZM(44)/ZM_A)/ZPM(NPERIOX(I-1),NPERIOY(J+3))
                  PDEP(NPERIOX(I),NPERIOY(J+3))   = PDEP(NPERIOX(I),NPERIOY(J+3))   &
                                                  & + ZM(0)*(ZM(45)/ZM_A)/ZPM(NPERIOX(I),NPERIOY(J+3))
                  PDEP(NPERIOX(I+1),NPERIOY(J+3)) = PDEP(NPERIOX(I+1),NPERIOY(J+3)) &
                                                  & + ZM(0)*(ZM(46)/ZM_A)/ZPM(NPERIOX(I+1),NPERIOY(J+3))
                  PDEP(NPERIOX(I+2),NPERIOY(J+3)) = PDEP(NPERIOX(I+2),NPERIOY(J+3)) &
                                                  & + ZM(0)*(ZM(47)/ZM_A)/ZPM(NPERIOX(I+2),NPERIOY(J+3))
                  PDEP(NPERIOX(I+3),NPERIOY(J+3)) = PDEP(NPERIOX(I+3),NPERIOY(J+3)) &
                                                  & + ZM(0)*(ZM(48)/ZM_A)/ZPM(NPERIOX(I+3),NPERIOY(J+3))
               END IF
            END IF   
         END IF
      !----------------------------!
      !  DEALING WITH UNDERSHOOTS  !
      !----------------------------!   
      ELSE IF (PDEP(I,J) < ZMIN_LOC(I,J)) THEN
         ZM(0) = (ZMIN_LOC(I,J)-PDEP(I,J))*ZPM(I,J)
         ZM(1) = MAX(PDEP(NPERIOX(I-1),NPERIOY(J-1))-ZMIN_LOC(NPERIOX(I-1),NPERIOY(J-1)),ZERO)*( &
               & ZPM(NPERIOX(I-1),NPERIOY(J-1)))
         ZM(2) = MAX(PDEP(I,NPERIOY(J-1))-ZMIN_LOC(I,NPERIOY(J-1)),ZERO)*(                       &
               & ZPM(I,NPERIOY(J-1)))
         ZM(3) = MAX(PDEP(NPERIOX(I+1),NPERIOY(J-1))-ZMIN_LOC(NPERIOX(I+1),NPERIOY(J-1)),ZERO)*( &
               & ZPM(NPERIOX(I+1),NPERIOY(J-1)))
         ZM(4) = MAX(PDEP(NPERIOX(I-1),J)-ZMIN_LOC(NPERIOX(I-1),J),ZERO)*(                       &
               & ZPM(NPERIOX(I-1),J))
         ZM(5) = MAX(PDEP(NPERIOX(I+1),J)-ZMIN_LOC(NPERIOX(I+1),J),ZERO)*(                       &
               & ZPM(NPERIOX(I+1),J))
         ZM(6) = MAX(PDEP(NPERIOX(I-1),NPERIOY(J+1))-ZMIN_LOC(NPERIOX(I-1),NPERIOY(J+1)),ZERO)*( &
               & ZPM(NPERIOX(I-1),NPERIOY(J+1)))
         ZM(7) = MAX(PDEP(I,NPERIOY(J+1))-ZMIN_LOC(I,NPERIOY(J+1)),ZERO)*(                       &
               & ZPM(I,NPERIOY(J+1)))
         ZM(8) = MAX(PDEP(NPERIOX(I+1),NPERIOY(J+1))-ZMIN_LOC(NPERIOX(I+1),NPERIOY(J+1)),ZERO)*( &
               & ZPM(NPERIOX(I+1),NPERIOY(J+1)))
         ZM_A  = ZERO 
         DO K=1,8
            ZM_A = ZM_A + ZM(K)
         END DO
         IF (ZM_A > ZM(0)) THEN
            PDEP(I,J) = ZMIN_LOC(I,J)
            PDEP(NPERIOX(I-1),NPERIOY(J-1)) = PDEP(NPERIOX(I-1),NPERIOY(J-1)) &
                                            & - ZM(0)*(ZM(1)/ZM_A)/ZPM(NPERIOX(I-1),NPERIOY(J-1))
            PDEP(I,NPERIOY(J-1))            = PDEP(I,NPERIOY(J-1))            &
                                            & - ZM(0)*(ZM(2)/ZM_A)/ZPM(I,NPERIOY(J-1))
            PDEP(NPERIOX(I+1),NPERIOY(J-1)) = PDEP(NPERIOX(I+1),NPERIOY(J-1)) &
                                            & - ZM(0)*(ZM(3)/ZM_A)/ZPM(NPERIOX(I+1),NPERIOY(J-1))
            PDEP(NPERIOX(I-1),J)            = PDEP(NPERIOX(I-1),J)            &
                                            & - ZM(0)*(ZM(4)/ZM_A)/ZPM(NPERIOX(I-1),J) 
            PDEP(NPERIOX(I+1),J)            = PDEP(NPERIOX(I+1),J)            &
                                            & - ZM(0)*(ZM(5)/ZM_A)/ZPM(NPERIOX(I+1),J)
            PDEP(NPERIOX(I-1),NPERIOY(J+1)) = PDEP(NPERIOX(I-1),NPERIOY(J+1)) &
                                            & - ZM(0)*(ZM(6)/ZM_A)/ZPM(NPERIOX(I-1),NPERIOY(J+1))
            PDEP(I,NPERIOY(J+1))            = PDEP(I,NPERIOY(J+1))            &
                                            & - ZM(0)*(ZM(7)/ZM_A)/ZPM(I,NPERIOY(J+1))
            PDEP(NPERIOX(I+1),NPERIOY(J+1)) = PDEP(NPERIOX(I+1),NPERIOY(J+1)) &
                                            & - ZM(0)*(ZM(8)/ZM_A)/ZPM(NPERIOX(I+1),NPERIOY(J+1))
         ELSE
            ZM(9)  = MAX(PDEP(NPERIOX(I-2),NPERIOY(J-2))-ZMIN_LOC(NPERIOX(I-2),NPERIOY(J-2)),ZERO)*( &
                   & ZPM(NPERIOX(I-2),NPERIOY(J-2)))
            ZM(10) = MAX(PDEP(NPERIOX(I-1),NPERIOY(J-2))-ZMIN_LOC(NPERIOX(I-1),NPERIOY(J-2)),ZERO)*( &
                   & ZPM(NPERIOX(I-1),NPERIOY(J-2)))
            ZM(11) = MAX(PDEP(I,NPERIOY(J-2))-ZMIN_LOC(I,NPERIOY(J-2)),ZERO)*(                       &
                   & ZPM(I,NPERIOY(J-2)))
            ZM(12) = MAX(PDEP(NPERIOX(I+1),NPERIOY(J-2))-ZMIN_LOC(NPERIOX(I+1),NPERIOY(J-2)),ZERO)*( &
                   & ZPM(NPERIOX(I+1),NPERIOY(J-2)))
            ZM(13) = MAX(PDEP(NPERIOX(I+2),NPERIOY(J-2))-ZMIN_LOC(NPERIOX(I+2),NPERIOY(J-2)),ZERO)*( &
                   & ZPM(NPERIOX(I+2),NPERIOY(J-2)))
            ZM(14) = MAX(PDEP(NPERIOX(I-2),NPERIOY(J-1))-ZMIN_LOC(NPERIOX(I-2),NPERIOY(J-1)),ZERO)*( &
                   & ZPM(NPERIOX(I-2),NPERIOY(J-1)))
            ZM(15) = MAX(PDEP(NPERIOX(I+2),NPERIOY(J-1))-ZMIN_LOC(NPERIOX(I+2),NPERIOY(J-1)),ZERO)*( &
                   & ZPM(NPERIOX(I+2),NPERIOY(J-1)))
            ZM(16) = MAX(PDEP(NPERIOX(I-2),J)-ZMIN_LOC(NPERIOX(I-2),J),ZERO)*(                       &
                   & ZPM(NPERIOX(I-2),J))
            ZM(17) = MAX(PDEP(NPERIOX(I+2),J)-ZMIN_LOC(NPERIOX(I+2),J),ZERO)*(                       &
                   & ZPM(NPERIOX(I+2),J))
            ZM(18) = MAX(PDEP(NPERIOX(I-2),NPERIOY(J+1))-ZMIN_LOC(NPERIOX(I-2),NPERIOY(J+1)),ZERO)*( &
                   & ZPM(NPERIOX(I-2),NPERIOY(J+1)))
            ZM(19) = MAX(PDEP(NPERIOX(I+2),NPERIOY(J+1))-ZMIN_LOC(NPERIOX(I+2),NPERIOY(J+1)),ZERO)*( &
                   & ZPM(NPERIOX(I+2),NPERIOY(J+1)))
            ZM(20) = MAX(PDEP(NPERIOX(I-2),NPERIOY(J+2))-ZMIN_LOC(NPERIOX(I-2),NPERIOY(J+2)),ZERO)*( &
                   & ZPM(NPERIOX(I-2),NPERIOY(J+2))) 
            ZM(21) = MAX(PDEP(NPERIOX(I-1),NPERIOY(J+2))-ZMIN_LOC(NPERIOX(I-1),NPERIOY(J+2)),ZERO)*( &
                   & ZPM(NPERIOX(I-1),NPERIOY(J+2)))
            ZM(22) = MAX(PDEP(I,NPERIOY(J+2))-ZMIN_LOC(I,NPERIOY(J+2)),ZERO)*(                       &
                   & ZPM(I,NPERIOY(J+2)))
            ZM(23) = MAX(PDEP(NPERIOX(I+1),NPERIOY(J+2))-ZMIN_LOC(NPERIOX(I+1),NPERIOY(J+2)),ZERO)*( &
                   & ZPM(NPERIOX(I+1),NPERIOY(J+2)))
            ZM(24) = MAX(PDEP(NPERIOX(I+2),NPERIOY(J+2))-ZMIN_LOC(NPERIOX(I+2),NPERIOY(J+2)),ZERO)*( &
                   & ZPM(NPERIOX(I+2),NPERIOY(J+2)))
            DO K=9,24
               ZM_A = ZM_A + ZM(K)
            END DO
            IF (ZM_A > ZM(0)) THEN
               PDEP(I,J) = ZMIN_LOC(I,J)
               PDEP(NPERIOX(I-1),NPERIOY(J-1)) = PDEP(NPERIOX(I-1),NPERIOY(J-1)) &
                                               & - ZM(0)*(ZM(1)/ZM_A)/ZPM(NPERIOX(I-1),NPERIOY(J-1))
               PDEP(I,NPERIOY(J-1))            = PDEP(I,NPERIOY(J-1))            &
                                               & - ZM(0)*(ZM(2)/ZM_A)/ZPM(I,NPERIOY(J-1))
               PDEP(NPERIOX(I+1),NPERIOY(J-1)) = PDEP(NPERIOX(I+1),NPERIOY(J-1)) &
                                               & - ZM(0)*(ZM(3)/ZM_A)/ZPM(NPERIOX(I+1),NPERIOY(J-1))
               PDEP(NPERIOX(I-1),J)            = PDEP(NPERIOX(I-1),J)            &
                                               & - ZM(0)*(ZM(4)/ZM_A)/ZPM(NPERIOX(I-1),J) 
               PDEP(NPERIOX(I+1),J)            = PDEP(NPERIOX(I+1),J)            &
                                               & - ZM(0)*(ZM(5)/ZM_A)/ZPM(NPERIOX(I+1),J)
               PDEP(NPERIOX(I-1),NPERIOY(J+1)) = PDEP(NPERIOX(I-1),NPERIOY(J+1)) &
                                               & - ZM(0)*(ZM(6)/ZM_A)/ZPM(NPERIOX(I-1),NPERIOY(J+1))
               PDEP(I,NPERIOY(J+1))            = PDEP(I,NPERIOY(J+1))            &
                                               & - ZM(0)*(ZM(7)/ZM_A)/ZPM(I,NPERIOY(J+1))
               PDEP(NPERIOX(I+1),NPERIOY(J+1)) = PDEP(NPERIOX(I+1),NPERIOY(J+1)) &
                                               & - ZM(0)*(ZM(8)/ZM_A)/ZPM(NPERIOX(I+1),NPERIOY(J+1))
               PDEP(NPERIOX(I-2),NPERIOY(J-2)) = PDEP(NPERIOX(I-2),NPERIOY(J-2)) &
                                               & - ZM(0)*(ZM(9)/ZM_A)/ZPM(NPERIOX(I-2),NPERIOY(J-2))
               PDEP(NPERIOX(I-1),NPERIOY(J-2)) = PDEP(NPERIOX(I-1),NPERIOY(J-2)) &
                                               & - ZM(0)*(ZM(10)/ZM_A)/ZPM(NPERIOX(I-1),NPERIOY(J-2))
               PDEP(I,NPERIOY(J-2))            = PDEP(I,NPERIOY(J-2))            &
                                               & - ZM(0)*(ZM(11)/ZM_A)/ZPM(I,NPERIOY(J-2))
               PDEP(NPERIOX(I+1),NPERIOY(J-2)) = PDEP(NPERIOX(I+1),NPERIOY(J-2)) &
                                               & - ZM(0)*(ZM(12)/ZM_A)/ZPM(NPERIOX(I+1),NPERIOY(J-2))
               PDEP(NPERIOX(I+2),NPERIOY(J-2)) = PDEP(NPERIOX(I+2),NPERIOY(J-2)) &
                                               & - ZM(0)*(ZM(13)/ZM_A)/ZPM(NPERIOX(I+2),NPERIOY(J-2))
               PDEP(NPERIOX(I-2),NPERIOY(J-1)) = PDEP(NPERIOX(I-2),NPERIO(J-1))  &
                                               & - ZM(0)*(ZM(14)/ZM_A)/ZPM(NPERIOX(I-2),NPERIO(J-1))
               PDEP(NPERIOX(I+2),NPERIOY(J-1)) = PDEP(NPERIOX(I+2),NPERIOY(J-1)) &
                                               & - ZM(0)*(ZM(15)/ZM_A)/ZPM(NPERIOX(I+2),NPERIOY(J-1))
               PDEP(NPERIOX(I-2),J)            = PDEP(NPERIOX(I-2),J)            &
                                               & - ZM(0)*(ZM(16)/ZM_A)/ZPM(NPERIOX(I-2),J)
               PDEP(NPERIOX(I+2),J)            = PDEP(NPERIOX(I+2),J)            &
                                               & - ZM(0)*(ZM(17)/ZM_A)/ZPM(NPERIOX(I+2),J)
               PDEP(NPERIOX(I-2),NPERIOY(J+1)) = PDEP(NPERIOX(I-2),NPERIOY(J+1)) &
                                               & - ZM(0)*(ZM(18)/ZM_A)/ZPM(NPERIOX(I-2),NPERIOY(J+1))
               PDEP(NPERIOX(I+2),NPERIOY(J+1)) = PDEP(NPERIOX(I+2),NPERIOY(J+1)) &
                                               & - ZM(0)*(ZM(19)/ZM_A)/ZPM(NPERIOX(I+2),NPERIOY(J+1))
               PDEP(NPERIOX(I-2),NPERIOY(J+2)) = PDEP(NPERIOX(I-2),NPERIOY(J+2)) &
                                               & - ZM(0)*(ZM(20)/ZM_A)/ZPM(NPERIOX(I-2),NPERIOY(J+2))
               PDEP(NPERIOX(I-1),NPERIOY(J+2)) = PDEP(NPERIOX(I-1),NPERIOY(J+2)) &
                                               & - ZM(0)*(ZM(21)/ZM_A)/ZPM(NPERIOX(I-1),NPERIOY(J+2))
               PDEP(I,NPERIOY(J+2))            = PDEP(I,NPERIOY(J+2))            &
                                               & - ZM(0)*(ZM(22)/ZM_A)/ZPM(I,NPERIOY(J+2))
               PDEP(NPERIOX(I+1),NPERIOY(J+2)) = PDEP(NPERIOX(I+1),NPERIOY(J+2)) &
                                               & - ZM(0)*(ZM(23)/ZM_A)/ZPM(NPERIOX(I+1),NPERIOY(J+2))
               PDEP(NPERIOX(I+2),NPERIOY(J+2)) = PDEP(NPERIOX(I+2),NPERIOY(J+2)) &
                                               & - ZM(0)*(ZM(24)/ZM_A)/ZPM(NPERIOX(I+2),NPERIOY(J+2))
            ELSE
               ZM(25) = MAX(PDEP(NPERIOX(I-3),NPERIOY(J-3))-ZMIN_LOC(NPERIOX(I-3),NPERIOY(J-3)),ZERO)*( &
                      & ZPM(NPERIOX(I-3),NPERIOY(J-3)))
               ZM(26) = MAX(PDEP(NPERIOX(I-2),NPERIOY(J-3))-ZMIN_LOC(NPERIOX(I-2),NPERIOY(J-3)),ZERO)*( &
                      & ZPM(NPERIOX(I-2),NPERIOY(J-3)))
               ZM(27) = MAX(PDEP(NPERIOX(I-1),NPERIOY(J-3))-ZMIN_LOC(NPERIOX(I-1),NPERIOY(J-3)),ZERO)*( &
                      & ZPM(NPERIOX(I-1),NPERIOY(J-1)))
               ZM(28) = MAX(PDEP(NPERIOX(I),NPERIOY(J-3))-ZMIN_LOC(NPERIOX(I),NPERIOY(J-3)),ZERO)*(     &
                      & ZPM(NPERIOX(I),NPERIOY(J-1)))
               ZM(29) = MAX(PDEP(NPERIOX(I+1),NPERIOY(J-3))-ZMIN_LOC(NPERIOX(I+1),NPERIOY(J-3)),ZERO)*( &
                      & ZPM(NPERIOX(I+1),NPERIOY(J-2)))
               ZM(30) = MAX(PDEP(NPERIOX(I+2),NPERIOY(J-3))-ZMIN_LOC(NPERIOX(I+2),NPERIOY(J-3)),ZERO)*( &
                      & ZPM(NPERIOX(I+2),NPERIOY(J-3)))
               ZM(31) = MAX(PDEP(NPERIOX(I+3),NPERIOY(J-3))-ZMIN_LOC(NPERIOX(I+3),NPERIOY(J-3)),ZERO)*( &
                      & ZPM(NPERIOX(I+3),NPERIOY(J-3)))
               ZM(32) = MAX(PDEP(NPERIOX(I-3),NPERIOY(J-2))-ZMIN_LOC(NPERIOX(I-3),NPERIOY(J-2)),ZERO)*( &
                      & ZPM(NPERIOX(I-3),NPERIOY(J-2)))
               ZM(33) = MAX(PDEP(NPERIOX(I+3),NPERIOY(J-2))-ZMIN_LOC(NPERIOX(I+3),NPERIOY(J-2)),ZERO)*( &
                      & ZPM(NPERIOX(I+3),NPERIOY(J-2)))
               ZM(34) = MAX(PDEP(NPERIOX(I-3),NPERIOY(J-1))-ZMIN_LOC(NPERIOX(I-3),NPERIOY(J-1)),ZERO)*( &
                      & ZPM(NPERIOX(I-3),NPERIOY(J-1)))
               ZM(35) = MAX(PDEP(NPERIOX(I+3),NPERIOY(J-1))-ZMIN_LOC(NPERIOX(I+3),NPERIOY(J-1)),ZERO)*( &
                      & ZPM(NPERIOX(I+3),NPERIOY(J-1)))
               ZM(36) = MAX(PDEP(NPERIOX(I-3),J)-ZMIN_LOC(NPERIOX(I-3),J),ZERO)*(                       &
                      & ZPM(NPERIOX(I-3),J))
               ZM(37) = MAX(PDEP(NPERIOX(I+3),J)-ZMIN_LOC(NPERIOX(I+3),J),ZERO)*(                       &
                      & ZPM(NPERIOX(I+3),J))
               ZM(38) = MAX(PDEP(NPERIOX(I-3),NPERIOY(J+1))-ZMIN_LOC(NPERIOX(I-3),NPERIOY(J+1)),ZERO)*( &
                      & ZPM(NPERIOX(I-3),NPERIOY(J+1)))
               ZM(39) = MAX(PDEP(NPERIOX(I+3),NPERIOY(J+1))-ZMIN_LOC(NPERIOX(I+3),NPERIOY(J+1)),ZERO)*( &
                      & ZPM(NPERIOX(I+3),NPERIOY(J+1)))
               ZM(40) = MAX(PDEP(NPERIOX(I-3),NPERIOY(J+2))-ZMIN_LOC(NPERIOX(I-3),NPERIOY(J+2)),ZERO)*( &
                      & ZPM(NPERIOX(I-3),NPERIOY(J+2))) 
               ZM(41) = MAX(PDEP(NPERIOX(I+3),NPERIOY(J+2))-ZMIN_LOC(NPERIOX(I+3),NPERIOY(J+2)),ZERO)*( &
                      & ZPM(NPERIOX(I+3),NPERIOY(J+2)))
               ZM(42) = MAX(PDEP(NPERIOX(I-3),NPERIOY(J+3))-ZMIN_LOC(NPERIOX(I-3),NPERIOY(J+3)),ZERO)*( &
                      & ZPM(NPERIOX(I-3),NPERIOY(J+3)))
               ZM(43) = MAX(PDEP(NPERIOX(I-2),NPERIOY(J+3))-ZMIN_LOC(NPERIOX(I-2),NPERIOY(J+3)),ZERO)*( &
                      & ZPM(NPERIOX(I-2),NPERIOY(J+3)))
               ZM(44) = MAX(PDEP(NPERIOX(I-1),NPERIOY(J+3))-ZMIN_LOC(NPERIOX(I-1),NPERIOY(J+3)),ZERO)*( &
                      & ZPM(NPERIOX(I-1),NPERIOY(J+3)))
               ZM(45) = MAX(PDEP(NPERIOX(I),NPERIOY(J+3))-ZMIN_LOC(NPERIOX(I),NPERIOY(J+3)),ZERO)*( &
                      & ZPM(NPERIOX(I),NPERIOY(J+3)))
               ZM(46) = MAX(PDEP(NPERIOX(I+1),NPERIOY(J+3))-ZMIN_LOC(NPERIOX(I+1),NPERIOY(J+3)),ZERO)*(     &
                      & ZPM(NPERIOX(I+1),NPERIOY(J+3)))
               ZM(47) = MAX(PDEP(NPERIOX(I+2),NPERIOY(J+3))-ZMIN_LOC(NPERIOX(I+2),NPERIOY(J+3)),ZERO)*( &
                      & ZPM(NPERIOX(I+1),NPERIOY(J+3)))
               ZM(48) = MAX(PDEP(NPERIOX(I+3),NPERIOY(J+3))-ZMIN_LOC(NPERIOX(I+3),NPERIOY(J+3)),ZERO)*( &
                      & ZPM(NPERIOX(I+3),NPERIOY(J+3)))
               DO K=25,48
                  ZM_A = ZM_A + ZM(K)
               END DO
               IF (ZM_A > ZM(0)) THEN
                  PDEP(I,J) = ZMIN_LOC(I,J)
                  PDEP(NPERIOX(I-1),NPERIOY(J-1)) = PDEP(NPERIOX(I-1),NPERIOY(J-1)) &
                                                  & - ZM(0)*(ZM(1)/ZM_A)/ZPM(NPERIOX(I-1),NPERIOY(J-1))
                  PDEP(I,NPERIOY(J-1))            = PDEP(I,NPERIOY(J-1))            &
                                                  & - ZM(0)*(ZM(2)/ZM_A)/ZPM(I,NPERIOY(J-1))
                  PDEP(NPERIOX(I+1),NPERIOY(J-1)) = PDEP(NPERIOX(I+1),NPERIOY(J-1)) &
                                                  & - ZM(0)*(ZM(3)/ZM_A)/ZPM(NPERIOX(I+1),NPERIOY(J-1))
                  PDEP(NPERIOX(I-1),J)            = PDEP(NPERIOX(I-1),J)            &
                                                  & - ZM(0)*(ZM(4)/ZM_A)/ZPM(NPERIOX(I-1),J) 
                  PDEP(NPERIOX(I+1),J)            = PDEP(NPERIOX(I+1),J)            &
                                                  & - ZM(0)*(ZM(5)/ZM_A)/ZPM(NPERIOX(I+1),J)
                  PDEP(NPERIOX(I-1),NPERIOY(J+1)) = PDEP(NPERIOX(I-1),NPERIOY(J+1)) &
                                                  & - ZM(0)*(ZM(6)/ZM_A)/ZPM(NPERIOX(I-1),NPERIOY(J+1))
                  PDEP(I,NPERIOY(J+1))            = PDEP(I,NPERIOY(J+1))            &
                                                  & - ZM(0)*(ZM(7)/ZM_A)/ZPM(I,NPERIOY(J+1))
                  PDEP(NPERIOX(I+1),NPERIOY(J+1)) = PDEP(NPERIOX(I+1),NPERIOY(J+1)) &
                                                  & - ZM(0)*(ZM(8)/ZM_A)/ZPM(NPERIOX(I+1),NPERIOY(J+1))
                  PDEP(NPERIOX(I-2),NPERIOY(J-2)) = PDEP(NPERIOX(I-2),NPERIOY(J-2)) &
                                                  & - ZM(0)*(ZM(9)/ZM_A)/ZPM(NPERIOX(I-2),NPERIOY(J-2))
                  PDEP(NPERIOX(I-1),NPERIOY(J-2)) = PDEP(NPERIOX(I-1),NPERIOY(J-2)) &
                                                  & - ZM(0)*(ZM(10)/ZM_A)/ZPM(NPERIOX(I-1),NPERIOY(J-2))
                  PDEP(I,NPERIOY(J-2))            = PDEP(I,NPERIOY(J-2))            &
                                                  & - ZM(0)*(ZM(11)/ZM_A)/ZPM(I,NPERIOY(J-2))
                  PDEP(NPERIOX(I+1),NPERIOY(J-2)) = PDEP(NPERIOX(I+1),NPERIOY(J-2)) &
                                                  & - ZM(0)*(ZM(12)/ZM_A)/ZPM(NPERIOX(I+1),NPERIOY(J-2))
                  PDEP(NPERIOX(I+2),NPERIOY(J-2)) = PDEP(NPERIOX(I+2),NPERIOY(J-2)) &
                                                  & - ZM(0)*(ZM(13)/ZM_A)/ZPM(NPERIOX(I+2),NPERIOY(J-2))
                  PDEP(NPERIOX(I-2),NPERIOY(J-1)) = PDEP(NPERIOX(I-2),NPERIO(J-1))  &
                                                  & - ZM(0)*(ZM(14)/ZM_A)/ZPM(NPERIOX(I-2),NPERIO(J-1))
                  PDEP(NPERIOX(I+2),NPERIOY(J-1)) = PDEP(NPERIOX(I+2),NPERIOY(J-1)) &
                                                  & - ZM(0)*(ZM(15)/ZM_A)/ZPM(NPERIOX(I+2),NPERIOY(J-1))
                  PDEP(NPERIOX(I-2),J)            = PDEP(NPERIOX(I-2),J)            &
                                                  & - ZM(0)*(ZM(16)/ZM_A)/ZPM(NPERIOX(I-2),J)
                  PDEP(NPERIOX(I+2),J)            = PDEP(NPERIOX(I+2),J)            &
                                                  & - ZM(0)*(ZM(17)/ZM_A)/ZPM(NPERIOX(I+2),J)
                  PDEP(NPERIOX(I-2),NPERIOY(J+1)) = PDEP(NPERIOX(I-2),NPERIOY(J+1)) &
                                                  & - ZM(0)*(ZM(18)/ZM_A)/ZPM(NPERIOX(I-2),NPERIOY(J+1))
                  PDEP(NPERIOX(I+2),NPERIOY(J+1)) = PDEP(NPERIOX(I+2),NPERIOY(J+1)) &
                                                  & - ZM(0)*(ZM(19)/ZM_A)/ZPM(NPERIOX(I+2),NPERIOY(J+1))
                  PDEP(NPERIOX(I-2),NPERIOY(J+2)) = PDEP(NPERIOX(I-2),NPERIOY(J+2)) &
                                                  & - ZM(0)*(ZM(20)/ZM_A)/ZPM(NPERIOX(I-2),NPERIOY(J+2))
                  PDEP(NPERIOX(I-1),NPERIOY(J+2)) = PDEP(NPERIOX(I-1),NPERIOY(J+2)) &
                                                  & - ZM(0)*(ZM(21)/ZM_A)/ZPM(NPERIOX(I-1),NPERIOY(J+2))
                  PDEP(I,NPERIOY(J+2))            = PDEP(I,NPERIOY(J+2))            &
                                                  & - ZM(0)*(ZM(22)/ZM_A)/ZPM(I,NPERIOY(J+2))
                  PDEP(NPERIOX(I+1),NPERIOY(J+2)) = PDEP(NPERIOX(I+1),NPERIOY(J+2)) &
                                                  & - ZM(0)*(ZM(23)/ZM_A)/ZPM(NPERIOX(I+1),NPERIOY(J+2))
                  PDEP(NPERIOX(I+2),NPERIOY(J+2)) = PDEP(NPERIOX(I+2),NPERIOY(J+2)) &
                                                  & - ZM(0)*(ZM(24)/ZM_A)/ZPM(NPERIOX(I+2),NPERIOY(J+2))
                  PDEP(NPERIOX(I-3),NPERIOY(J-3)) = PDEP(NPERIOX(I-3),NPERIOY(J-3)) &
                                                  & - ZM(0)*(ZM(25)/ZM_A)/ZPM(NPERIOX(I-3),NPERIOY(J-3))
                  PDEP(NPERIOX(I-2),NPERIOY(J-3)) = PDEP(NPERIOX(I-2),NPERIOY(J-3)) &
                                                  & - ZM(0)*(ZM(26)/ZM_A)/ZPM(I,NPERIOY(J-1))
                  PDEP(NPERIOX(I-1),NPERIOY(J-3)) = PDEP(NPERIOX(I-1),NPERIOY(J-3)) &
                                                  & - ZM(0)*(ZM(27)/ZM_A)/ZPM(NPERIOX(I-1),NPERIOY(J-3))
                  PDEP(NPERIOX(I),NPERIOY(J-3))   = PDEP(NPERIOX(I),NPERIOY(J-3))   &
                                                  & - ZM(0)*(ZM(28)/ZM_A)/ZPM(NPERIOX(I-1),NPERIOY(J-3)) 
                  PDEP(NPERIOX(I+1),NPERIOY(J-3)) = PDEP(NPERIOX(I+1),NPERIOY(J-3)) &
                                                  & - ZM(0)*(ZM(29)/ZM_A)/ZPM(NPERIOX(I+1),NPERIOY(J-3))
                  PDEP(NPERIOX(I+2),NPERIOY(J-3)) = PDEP(NPERIOX(I+2),NPERIOY(J-3)) &
                                                  & - ZM(0)*(ZM(30)/ZM_A)/ZPM(NPERIOX(I+2),NPERIOY(J-3))
                  PDEP(NPERIOX(I+3),NPERIOY(J-3)) = PDEP(NPERIOX(I+3),NPERIOY(J-3)) &
                                                  & - ZM(0)*(ZM(31)/ZM_A)/ZPM(NPERIOX(I+3),NPERIOY(J-3))
                  PDEP(NPERIOX(I-3),NPERIOY(J-2)) = PDEP(NPERIOX(I-3),NPERIOY(J-2)) &
                                                  & - ZM(0)*(ZM(32)/ZM_A)/ZPM(NPERIOX(I-3),NPERIOY(J-2))
                  PDEP(NPERIOX(I+3),NPERIOY(J-2)) = PDEP(NPERIOX(I+3),NPERIOY(J-2)) &
                                                  & - ZM(0)*(ZM(33)/ZM_A)/ZPM(NPERIOX(I+3),NPERIOY(J-2))
                  PDEP(NPERIOX(I-3),NPERIOY(J-1)) = PDEP(NPERIOX(I-3),NPERIOY(J-1)) &
                                                  & - ZM(0)*(ZM(34)/ZM_A)/ZPM(NPERIOX(I-3),NPERIOY(J-1))
                  PDEP(NPERIOX(I+3),NPERIOY(J-1)) = PDEP(NPERIOX(I+3),NPERIOY(J-1)) &
                                                  & - ZM(0)*(ZM(35)/ZM_A)/ZPM(NPERIOX(I+3),NPERIOY(J-1))
                  PDEP(NPERIOX(I-3),NPERIOY(J))   = PDEP(NPERIOX(I-3),NPERIOY(J))   &
                                                  & - ZM(0)*(ZM(36)/ZM_A)/ZPM(NPERIOX(I-3),NPERIOY(J))
                  PDEP(NPERIOX(I+3),NPERIOY(J))   = PDEP(NPERIOX(I+3),NPERIOY(J))   &
                                                  & - ZM(0)*(ZM(37)/ZM_A)/ZPM(NPERIOX(I+3),NPERIOY(J))
                  PDEP(NPERIOX(I-3),NPERIOY(J+1)) = PDEP(NPERIOX(I-3),NPERIO(J+1))  &
                                                  & - ZM(0)*(ZM(38)/ZM_A)/ZPM(NPERIOX(I-3),NPERIO(J+1))
                  PDEP(NPERIOX(I+3),NPERIOY(J+1)) = PDEP(NPERIOX(I+3),NPERIOY(J+1)) &
                                                  & - ZM(0)*(ZM(39)/ZM_A)/ZPM(NPERIOX(I+3),NPERIOY(J+1))
                  PDEP(NPERIOX(I-3),NPERIOY(J+2)) = PDEP(NPERIOX(I-3),NPERIOY(J+2)) &
                                                  & - ZM(0)*(ZM(40)/ZM_A)/ZPM(NPERIOX(I-3),NPERIOY(J+2))
                  PDEP(NPERIOX(I+3),NPERIOY(J+2)) = PDEP(NPERIOX(I+3),NPERIOY(J+2)) &
                                                  & - ZM(0)*(ZM(41)/ZM_A)/ZPM(NPERIOX(I+3),NPERIOY(J+2))
                  PDEP(NPERIOX(I-3),NPERIOY(J+3)) = PDEP(NPERIOX(I-3),NPERIOY(J+3)) &
                                                  & - ZM(0)*(ZM(42)/ZM_A)/ZPM(NPERIOX(I-3),NPERIOY(J+3))
                  PDEP(NPERIOX(I-2),NPERIOY(J+3)) = PDEP(NPERIOX(I-2),NPERIOY(J+3)) &
                                                  & - ZM(0)*(ZM(43)/ZM_A)/ZPM(NPERIOX(I-2),NPERIOY(J+3))
                  PDEP(NPERIOX(I-1),NPERIOY(J+3)) = PDEP(NPERIOX(I-1),NPERIOY(J+3)) &
                                                  & - ZM(0)*(ZM(44)/ZM_A)/ZPM(NPERIOX(I-1),NPERIOY(J+3))
                  PDEP(NPERIOX(I),NPERIOY(J+3))   = PDEP(NPERIOX(I),NPERIOY(J+3))   &
                                                  & - ZM(0)*(ZM(45)/ZM_A)/ZPM(NPERIOX(I),NPERIOY(J+3))
                  PDEP(NPERIOX(I+1),NPERIOY(J+3)) = PDEP(NPERIOX(I+1),NPERIOY(J+3)) &
                                                  & - ZM(0)*(ZM(46)/ZM_A)/ZPM(NPERIOX(I+1),NPERIOY(J+3))
                  PDEP(NPERIOX(I+2),NPERIOY(J+3)) = PDEP(NPERIOX(I+2),NPERIOY(J+3)) &
                                                  & - ZM(0)*(ZM(47)/ZM_A)/ZPM(NPERIOX(I+2),NPERIOY(J+3))
                  PDEP(NPERIOX(I+3),NPERIOY(J+3)) = PDEP(NPERIOX(I+3),NPERIOY(J+3)) &
                                                  & - ZM(0)*(ZM(48)/ZM_A)/ZPM(NPERIOX(I+3),NPERIOY(J+3))
               END IF   
            END IF   
         END IF
      END IF   
   ENDDO
ENDDO


END SUBROUTINE ILMC_FILTER 
!*************************************************************************************!
!*************************************************************************************!
!*************************************************************************************!
!*************************************************************************************!
!*************************************************************************************!
!*************************************************************************************!
SUBROUTINE BC_MASS_FIXER(PDEP,PARR,PXWEI,KXLAG,PYWEI,KYLAG,PMASS0)
  
IMPLICIT NONE

REAL(8),          INTENT(INOUT) :: PDEP(GXPT,GYPT)
REAL(8),          INTENT(IN)    :: PARR(GXPT,GYPT)
REAL(8),          INTENT(IN)    :: PMASS0
REAL(8),          INTENT(IN)    :: PXWEI(6,GXPT,GYPT),PYWEI(6,GXPT,GYPT)
INTEGER(8),       INTENT(IN)    :: KXLAG(6,GXPT,GYPT),KYLAG(6,GXPT,GYPT)

INTEGER(8)                      :: I,  J,  K,  L
REAL(8)                         :: ZMASS,ZDM,ZDEP_L,ZEPS
REAL(8)                         :: ZMIN_LOC,ZMAX_LOC,ZWS
REAL(8)                         :: ZW(GXPT,GYPT),ZDEP(GXPT,GYPT)

ZEPS = ABS(100*EPSILON(ONE))

!* Compute Mass error
ZMASS = ZERO
DO I=1,GXPT
   DO J=1,GYPT
      ZMASS = ZMASS + PDEP(I,J)
   ENDDO
ENDDO
ZDM = ZMASS - PMASS0


IF (ABS(ZDM) < ZEPS) THEN
   !* mass conserving coefficient
   DO I=1,GXPT
      DO J=1,GYPT
         ZDEP_L = ZERO
         DO K=5,6
            DO L=5,6
               ZDEP_L = ZDEP_L + PXWEI(K,I,J)*PYWEI(L,I,J)*(   &
                      &   PARR(KXLAG(K,I,J),KYLAG(L,I,J)) )
            ENDDO
         ENDDO
         ZW(I,J) = MAX(ZERO,SIGN(ZDM)*SIGN(PDEP(I,J)-ZDEP_L)*( &
                 & EXP(RQ*LOG(ABS(PDEP(I,J)-ZDEP_L))) ) )
      ENDDO
   ENDDO
   !* compute sum of weights
   ZWS = ZERO
   DO I=1,GXPT
      DO J=1,GYPT
         ZWS = ZWS + ZW(I,J)
      ENDDO
   ENDDO
   !* mass correction
   IF (ZWS > ZEPS) THEN
      DO I=1,GXPT
         DO J=1,GYPT
            PDEP(I,J) = PDEP(I,J) - (ZDM/ZWS)*ZW(I,J)
         ENDDO
      ENDDO
   END IF
   ! shape preservation
   DO I=1,GXPT
      DO J=1,GYPT
         ZMAX_LOC  =  MAX(                            &
                   & PARR(KXLAG(2,I,J),KYLAG(2,I,J)), &
                   & PARR(KXLAG(2,I,J),KYLAG(3,I,J)), &
                   & PARR(KXLAG(3,I,J),KYLAG(2,I,J)), &
                   & PARR(KXLAG(3,I,J),KYLAG(3,I,J))  )
         ZMIN_LOC  =  MIN(                            &
                   & PARR(KXLAG(2,I,J),KYLAG(2,I,J)), &
                   & PARR(KXLAG(2,I,J),KYLAG(3,I,J)), &
                   & PARR(KXLAG(3,I,J),KYLAG(2,I,J)), &
                   & PARR(KXLAG(3,I,J),KYLAG(3,I,J))  )
         ZDEP(I,J) = MAX(ZERO,SIGN(ZDM))*MAX(ZMIN_LOC,PDEP(I,J))  &
                   & + MAX(ZERO,SIGN(-ZDM))*MIN(ZMAX_LOC,PDEP(I,J))
      ENDDO
   ENDDO
   !* Compute Mass error and correction
   ZMASS = ZERO
   DO I=1,GXPT
      DO J=1,GYPT
         ZMASS = ZMASS + ZDEP(I,J)
      ENDDO
   ENDDO
   DO I=1,GXPT
      DO J=1,GYPT
         PDEP(I,J) = (PMASS0/ZMASS)*ZDEP(I,J)
      ENDDO
   ENDDO
END IF

END SUBROUTINE BC_MASS_FIXER 
!*************************************************************************************!
!*************************************************************************************!
!*************************************************************************************!
!*************************************************************************************!

END MODULE MOD_SLAG
!=====================================================================================!
!=====================================================================================!
!=====================================================================================!
