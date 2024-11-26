!#################################  MODELE JOUET 2D  ##################################!
!#                                                                                    #!
!# auteurs : F.Voitus                                                                 #!
!# sujet   : Module servant à faire l'ensemble des calculs dynamiques                 #!
!#                                                                                    #!
!######################################################################################!

!=====================================================================!
!=========================  MODULE PRINCIPAL  ========================!
!=====================================================================!

MODULE MOD_MPDTA
  
  USE MOD_PARAM, ONLY : GXPT,GYPT,ONE,TWO,ZERO,HALF,RDT,RDX,RDY,LSETTLS, &
                      & NPERIOX,NPERIOY,LPERIO,NORD1,NORD2,RFCT,REPS,    &
                      & LPRINTLEV,LFCTM,LMPDATA_FV,LGAUGE
  USE MOD_SETUP, ONLY : MPDATA,RHSVAR
 
CONTAINS

  !****************************************************************************************!
  !****************************************************************************************!
  !****************************************************************************************!
  !****************************************************************************************!
  !****************************************************************************************!
  SUBROUTINE MPDATA_TRANSPORT_SCHEME(RHS_DEP,RHS_ARR,MPDTA, &
             & PM,PXDOT0,PYDOT0,PXDOT9,PYDOT9,PXDOT1,PYDOT1,KSTEP,CPART)
  
  IMPLICIT NONE

  TYPE(RHSVAR),             INTENT(INOUT)  :: RHS_ARR
  TYPE(RHSVAR),             INTENT(INOUT)  :: RHS_DEP
  TYPE(MPDATA),             INTENT(INOUT)  :: MPDTA
  REAL(8),                  INTENT(IN)     :: PM(GXPT,GYPT)  
  REAL(8),                  INTENT(IN)     :: PXDOT0(GXPT,GYPT),PXDOT1(GXPT,GYPT)
  REAL(8),                  INTENT(IN)     :: PXDOT9(GXPT,GYPT),PYDOT9(GXPT,GYPT)
  REAL(8),                  INTENT(IN)     :: PYDOT0(GXPT,GYPT),PYDOT1(GXPT,GYPT)
  INTEGER(8),               INTENT(IN)     :: KSTEP
  CHARACTER(LEN=1),         INTENT(IN)     :: CPART

  REAL(8)                                  :: ZLIN,ZWEI
  REAL(8)                                  :: ZWIND_HOR_MAX
  INTEGER(8)                               :: I,J

  
  ZLIN = ONE
  IF (CPART == 'LI') ZLIN =ZERO
    
  !**************************************************!
  !* Advecting velocity at interfaces at (t+0.5*dt) *!
  !**************************************************!
     
    IF (LSETTLS) THEN
       DO I=1,GXPT
          DO J=1,GYPT
             MPDTA%XDOT(I,J)   = ZLIN*HALF*(                                 &
                               & (1.5d0*PXDOT0(I,J)-0.5d0*PXDOT9(I,J))       &
                               & + (1.5d0*PXDOT0(NPERIOX(I+1),J)             &
                               & -0.5d0*PXDOT9(NPERIOX(I+1),J)) )
             MPDTA%YDOT(I,J)   = ZLIN*HALF*(                                 &
                               & (1.5d0*PYDOT0(I,J)-0.5d0*PYDOT9(I,J))       &
                               & + (1.5d0*PYDOT0(I,NPERIOY(J+1))             &
                               & -0.5d0*PYDOT9(I,NPERIOY(J+1))) )
          END DO
       END DO
    ELSE
       DO I=1,GXPT
          DO J=1,GYPT
             MPDTA%XDOT(I,J)   = ZLIN*HALF*(                                 &
                               &   HALF*(PXDOT0(NPERIOX(I+1),J)+PXDOT0(I,J)) &
                               & + HALF*(PXDOT1(NPERIOX(I+1),J)+PXDOT1(I,J)) )
             MPDTA%YDOT(I,J)   = ZLIN*HALF*(                                 &
                               &   HALF*(PYDOT0(I,J)+PYDOT0(I,NPERIOY(J+1))) &
                               & + HALF*(PYDOT1(I,NPERIOY(J+1))+PYDOT1(I,J)) )
          END DO
       END DO
     ENDIF


     !****************************************************************!
     ! STABILITY CHECK-POINT MAXIMUM OF HORIZONTAL WIND VELOCITY      !
     !****************************************************************!
     ZWIND_HOR_MAX = MAX(MAXVAL(MPDTA%XDOT),MAXVAL(MPDTA%YDOT))       !
                                                                      !
     IF (LPRINTLEV) THEN                                              !
        WRITE(*,*) 'NSTEP : ',KSTEP,         &                        !
                 & '--> MAX.WIND_HOR : ',ZWIND_HOR_MAX                !
     END IF                                                           !
     IF (ZWIND_HOR_MAX .GE. 300.d0) THEN                              !
        WRITE(*,*) 'WIND TOO STRONG EXPLOSION AT TIMESTEP = ',KSTEP   !
        STOP                                                          !
     END IF                                                           !
     !****************************************************************!
     
     !*  mass transport scheme  *!
     CALL MPDATA_MASS(RHS_DEP%M,RHS_ARR%M,MPDTA%XDOT,MPDTA%YDOT) 

     MPDTA%M(0,:,:) = RHS_ARR%M(:,:)
     MPDTA%M(1,:,:) = RHS_DEP%M(:,:)

     !*  mass-weighted tracer transport scheme  *!
     CALL MPDATA_FIELD(RHS_DEP%Q,RHS_ARR%Q,MPDTA%XDOT,MPDTA%YDOT,MPDTA%M)

     
END SUBROUTINE MPDATA_TRANSPORT_SCHEME
!######################################################################################!
!**************************************************************************************!
!**************************************************************************************!
!**************************************************************************************!
!**************************************************************************************!
!**************************************************************************************!
!**************************************************************************************!
!**************************************************************************************!
!**************************************************************************************!
!**************************************************************************************!
!**************************************************************************************!
!**************************************************************************************!
!**************************************************************************************!
!**************************************************************************************!
!**************************************************************************************!
!**************************************************************************************!
!**************************************************************************************!
!######################################################################################!
SUBROUTINE MPDATA_MASS(PM_DEP,PM_ARR,PXDOT,PYDOT)

  IMPLICIT NONE

  REAL(8),                       INTENT(OUT)    :: PM_DEP(GXPT,GYPT)
  REAL(8),                       INTENT(OUT)    :: PM_ARR(GXPT,GYPT)
  REAL(8),                       INTENT(INOUT)  :: PXDOT(GXPT,GYPT)
  REAL(8),                       INTENT(INOUT)  :: PYDOT(GXPT,GYPT)

  REAL(8)                                       :: ZUN(GXPT,GYPT)

  ZUN(:,:) = ONE
  
  CALL MPDATA_ADV_XY(PM_DEP,PM_ARR,PXDOT,PYDOT, &
       & ZUN,ZUN,ZERO,RDT,LDMASS=.TRUE.)     
  
END SUBROUTINE MPDATA_MASS
!######################################################################################!
!**************************************************************************************!
!**************************************************************************************!
!**************************************************************************************!
!######################################################################################!
SUBROUTINE MPDATA_FIELD(PPSI_DEP,PPSI_ARR,PXDOT,PYDOT,PM,PPSI_MIN)

  IMPLICIT NONE

  REAL(8),          INTENT(OUT)    :: PPSI_DEP(GXPT,GYPT)
  REAL(8),          INTENT(IN)     :: PPSI_ARR(GXPT,GYPT)
  REAL(8),          INTENT(INOUT)  :: PXDOT(GXPT,GYPT)
  REAL(8),          INTENT(INOUT)  :: PYDOT(GXPT,GYPT)
  REAL(8),          INTENT(IN)     :: PM(0:1,GXPT,GYPT)
 

  CALL MPDATA_ADV_XY(PPSI_DEP,PPSI_ARR,PXDOT,PYDOT, &
       & PM(1,:,:),PM(0,:,:),PPSI_MIN,RDT)

END SUBROUTINE  MPDATA_FIELD
!#######################################################################################!
!***************************************************************************************!
!***************************************************************************************!
!***************************************************************************************!
!***************************************************************************************!
!***************************************************************************************!
!***************************************************************************************!
!***************************************************************************************!
!***************************************************************************************!
!#######################################################################################!
SUBROUTINE MPDATA_ADV_XY(PPSI_DEP,PPSI_ARR,PXDOT,PYDOT,PM1,PM0,PDT,PPSI_MIN,LDMASS)

  IMPLICIT NONE

  REAL(8),                       INTENT(OUT)    :: PPSI_DEP(GXPT,GYPT)
  REAL(8),                       INTENT(IN)     :: PPSI_ARR(GXPT,GYPT)
  REAL(8),                       INTENT(INOUT)  :: PXDOT(GXPT,GYPT)
  REAL(8),                       INTENT(INOUT)  :: PYDOT(GXPT,GYPT)
  REAL(8),                       INTENT(IN)     :: PM1(GXPT,GYPT)
  REAL(8),                       INTENT(IN)     :: PM0(GXPT,GYPT)
  REAL(8),                       INTENT(IN)     :: PDT,PPSI_MIN
  LOGICAL,OPTIONAL,              INTENT(IN)     :: LDMASS
  
  
  REAL(8)                                       :: YPSI(0:NORD,GXPT,GYPT)
  REAL(8)                                       :: YXDOT(NORD,GXPT,GYPT)
  REAL(8)                                       :: YYDOT(NORD,GXPT,GYPT)
  REAL(8)                                       :: YPSIG(GXPT,GYPT)
  REAL(8)                                       :: ZIN(GXPT,GYPT),ZOUT(GXPT,GYPT)
  REAL(8)                                       :: ZMAX,ZMIN
  REAL(8)                                       :: ZPSI_MIN,ZXDOT_MONO,ZYDOT_MONO
  INTEGER(8)                                    :: I,J,K
  LOGICAL                                       :: LLMASS

  
  IF (PRESENT(LDMASS)) THEN
     LLMASS = LDMASS
  ELSE
     LLMASS = .FALSE.
  END IF   

  IF (LLMASS) THEN
     ZPSI_MIN = ZERO
  ELSE
     ZPSI_MIN = MIN(ZERO,PPSI_MIN,MINVAL(PPSI_ARR))
  END IF   

  DO I=1,GXPT
     DO J=1,GYPT
        YPSI(0,I,J)  = PPSI_ARR(I,J) - ZPSI_MIN
        YXDOT(1,I,J) = PXDOT(I,J)
        YYDOT(1,I,J) = PYDOT(I,J)
     END DO   
  END DO

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !***************************************!
  ! DONOR-CELL UPWIND FIRST-ORDER SCHEME *!
  !***************************************!
  DO I=1,GXPT
     DO J=1,GYPT
        YPSI(1,I,J)       = YPSI(0,I,J)                                                                &
                          & -(ONE/PM0(I,J))*(PDT/RDX)*(                                                &
                          & (MIN(ZERO,YXDOT(1,I,J))*YPSI(0,NPERIOX(I+1),J)                             &
                          & +MAX(ZERO,YXDOT(1,I,J))*YPSI(0,I,J))                                       &
                          & -(MIN(ZERO,YXDOT(1,NPERIOX(I-1),J))*YPSI(0,I,J)                            & 
                          & +MAX(ZERO,YXDOT(1,NPERIOX(I-1),J))*YPSI(0,NPERIOX(I-1),J)) )               &
                          & -(ONE/PM0(I,J))*(PDT/RDY)*(                                                &
                          & (MIN(ZERO,YYDOT(1,I,J))*YPSI(0,I,NPERIOY(J+1))                             &
                          & +MAX(ZERO,YYDOT(1,I,J))*YPSI(0,I,J))                                       &
                          & -(MIN(ZERO,YYDOT(1,I,NPERIOY(J-1)))*YPSI(0,I,J)                            &  
                          & +MAX(ZERO,YYDOT(1,I,NPERIOY(J-1)))*YPSI(0,I,NPERIOY(J-1))) )
        YPSIG(I,J)        = (PM0(I,J)/PM1(I,J))*YPSI(1,I,J)
     END DO   
  END DO
  
  !********************!
  !  CORRECTION-STEP  *!
  !********************!
  IF (.NOT.LGAUGE) THEN
     DO K=1,NORD-1
       IF (LMPDATA_FV) THEN ! MPDATA-FV FORMULATION
         DO I=1,GXPT
            DO J=1,GYPT
               ZDIV_XA    = HALF*(                                                                     & !0.5*div(i,j)
                          & (((YXDOT(K,I,J)*(ABS(YPSIG(NPERIOX(I+1),J))+ABS(YPSIG(I,J)))               &
                          & -YXDOT(K,NPERIOX(I-1),J)*(ABS(YPSIG(I,J))+ABS(YPSIG(NPERIOX(I-1),J))))/RDX &
                          & +(YYDOT(K,I,J)*(ABS(YPSIG(I,NPERIOY(J+1)))+ABS(YPSIG(I,J)))                &
                          & -YYDOT(K,I,NPERIOY(J-1))*(ABS(YPSIG(I,J))+ABS(YPSIG(I,NPERIOY(J-1)))))/RDY)&
                          & /(0.25d0*(ABS(YPSIG(NPERIOX(I+1),J))+ABS(YPSIG(NPERIOX(I-1),J))            &
                          & +ABS(YPSIG(I,NPERIOY(J+1)))+ABS(YPSIG(I,NPERIOY(J-1)))                     &
                          & +(TWO+TWO)*ABS(YPSIG(I,J)))+REPS))                                         &
                          & +(((YXDOT(K,NPERIOX(I+1),J)*(ABS(YPSIG(NPERIOX(I+2),J))                    & !0.5*div(i+1,j)
                          & +ABS(YPSIG(NPERIOX(I+1),J)))-YXDOT(K,I,J)*(ABS(YPSIG(I,J))                 &
                          & +ABS(YPSIG(NPERIOX(I+1),J))))/RDX                                          &
                          & +(YYDOT(K,NPERIOX(I+1),J)*(ABS(YPSIG(NPERIOX(I+1),NPERIOY(J+1)))           &
                          & +ABS(YPSIG(NPERIOX(I+1),J)))-YYDOT(K,NPERIOX(I+1),NPERIOY(J-1))*(          &
                          & ABS(YPSIG(NPERIOX(I+1),J))+ABS(YPSIG(NPERIOX(I+1),NPERIOY(J-1)))))/RDY)    & ! option A : div(i+0.5,j)
                          & /(0.25d0*(ABS(YPSIG(NPERIOX(I+2),J))+ABS(YPSIG(I,J))                       & ! ---> 0.5[div(i,j)+div(i+1,j)]
                          & +ABS(YPSIG(NPERIOX(I+1),NPERIOY(J+1)))                                     &
                          & +ABS(YPSIG(NPERIOX(I+1),NPERIOY(J-1)))                                     &
                          & +(TWO+TWO)*ABS(YPSIG(NPERIOX(I+1),J)))+REPS)) )
               
               ZDIV_XB    = ( (((YXDOT(K,NPERIOX(I+1),J)+YXDOT(K,I,J))*ABS(YPSIG(NPERIOX(I+1),J))      & ! div(i+0.5,j) : option B
                          &  -(YXDOT(K,NPERIOX(I-1),J)+YXDOT(K,I,J))*ABS(YPSIG(NPERIOX(I),J)))/RDX)    & ! interpolate fields at interfaces 
                          &  +0.25d0*(((YYDOT(K,I,J)+YYDOT(K,NPERIOX(I+1),J))*(                        & ! then compute fluxes at (i,j+0.5) and (i,j-0.5)
                          &  ABS(YPSIG(NPERIOX(I+1),NPERIOY(J+1)))+ABS(YPSIG(I,NPERIOY(J+1)))          &
                          & +ABS(YPSIG(NPERIOX(I+1),J))+ABS(YPSIG(I,J))  )                             &
                          & -(YYDOT(K,I,NPERIOY(J-1))+YYDOT(K,NPERIOX(I+1),NPERIOY(J-1)))*(            &
                          &  ABS(YPSIG(NPERIOX(I+1),J))+ABS(YPSIG(I,NPERIOY(J-1)))                     &
                          & +ABS(YPSIG(NPERIOX(I+1),NPERIOY(J-1)))+ABS(YPSIG(I,J))))/RDY) )            &
                          & /(0.2d0*(ABS(YPSIG(NPERIOX(I+1),NPERIOY(J-1)))+ABS(YPSIG(I,NPERIOY(J-1)))  &
                          &  ABS(YPSIG(NPERIOX(I+1),NPERIOY(J+1)))+ABS(YPSIG(I,NPERIOY(J+1)))          &
                          & +(ONE+TWO)*(ABS(YPSIG(NPERIOX(I+1),J))+ABS(YPSIG(I,J))))+REPS)

               ZDIV_XC    = ( (((YXDOT(K,NPERIOX(I+1),J)+YXDOT(K,I,J))*ABS(YPSIG(NPERIOX(I+1),J))      & ! div(i+0.5,j) : option C
                          & -(YXDOT(K,NPERIOX(I-1),J)+YXDOT(K,I,J))*ABS(YPSIG(NPERIOX(I),J)))/RDX)     & ! compute fluxes first then 
                          & +0.5d0*((YYDOT(K,I,J)*(ABS(YPSIG(I,NPERIOY(J+1)))+ABS(YPSIG(I,J)))         & ! interpolate the result at interfaces
                          & +YYDOT(K,NPERIOX(I+1),J)*(ABS(YPSIG(NPERIOX(I+1),NPERIOY(J+1)))            & ! (i,j+0.5) and (i,j-0.5)
                          & +ABS(YPSIG(NPERIOX(I+1),J)))-YYDOT(K,I,NPERIOY(J-1))*(                     &
                          & ABS(YPSIG(I,NPERIOY(J-1)))+ABS(YPSIG(I,J)))                                &
                          & -YYDOT(K,NPERIOX(I+1),NPERIOY(J-1))*(ABS(YPSIG(NPERIOX(I+1),J))            &
                          & +ABS(YPSIG(NPERIOX(I+1),NPERIOY(J-1)))))/RDY) )                            &
                          & /(0.2d0*(ABS(YPSIG(NPERIOX(I+1),NPERIOY(J-1)))+ABS(YPSIG(I,NPERIOY(J-1)))  &
                          &  ABS(YPSIG(NPERIOX(I+1),NPERIOY(J+1)))+ABS(YPSIG(I,NPERIOY(J+1)))          &
                          & +(ONE+TWO)*(ABS(YPSIG(NPERIOX(I+1),J))+ABS(YPSIG(I,J))))+REPS)
               
               YXDOT(K+1,I,J) = ABS(YXDOT(K,I,J))*( (ABS(YPSIG(NPERIOX(I+1),J))-ABS(YPSIG(I,J)))       &  
                          & /((ABS(YPSIG(I,J))+ABS(YPSIG(NPERIOX(I+1),J)))+REPS) )                     &
                          & -PDT*(YXDOT(K,I,J)/(PM0(NPERIOX(I+1),J)+PM0(I,J)))*(                       & 
                          & (ONE-ROPT1)*ZDIV_XA + ROPT1*((ONE-ROPT2)*ZDIV_XB+ROPT2*ZDIV_XC)            &
                          & +HALF*((PM1(NPERIOX(I+1),J)-PM0(NPERIOX(I+1),J)+PM1(I,J)-PM0(I,J))/PDT))


               
               ZDIV_YA    = HALF*(                                                                     & !option A : 0.5*div(i,j)
                          & (((YXDOT(K,I,J)*(ABS(YPSIG(NPERIOX(I+1),J))+ABS(YPSIG(I,J)))               & !---> 0.5[div(i,j)+div(i,j+1)]
                          & -YXDOT(K,NPERIOX(I-1),J)*(ABS(YPSIG(I,J))+ABS(YPSIG(NPERIOX(I-1),J))))/RDX &
                          & +(YYDOT(K,I,J)*(ABS(YPSIG(I,NPERIOY(J+1)))+ABS(YPSIG(I,J)))                &
                          & -YYDOT(K,I,NPERIOY(J-1))*(ABS(YPSIG(I,J))+ABS(YPSIG(I,NPERIOY(J-1)))))/RDY)&
                          & /(0.25d0*(ABS(YPSIG(NPERIOX(I+1),J))+ABS(YPSIG(NPERIOX(I-1),J))            &
                          & +ABS(YPSIG(I,NPERIOY(J+1)))+ABS(YPSIG(I,NPERIOY(J-1)))                     &
                          & +(TWO+TWO)*ABS(YPSIG(I,J)))+REPS))                                         &
                          & +(((YXDOT(K,I,NPERIOX(J+1))*(ABS(YPSIG(NPERIOX(I+1),NPERIOY(J+1)))         & 
                          & +ABS(YPSIG(I,NPERIOY(J+1))))-YXDOT(K,NPERIOX(I-1),NPERIOY(J+1))*(          &
                          & ABS(YPSIG(NPERIOX(I-1),NPERIOY(J+1)))+ABS(YPSIG(I,NPERIOY(J+1)))))/RDX     &
                          & +(YYDOT(K,I,NPERIOY(J+1))*(ABS(YPSIG(I,NPERIOY(J+1)))                      &
                          & +ABS(YPSIG(I,NPERIOY(J+2))))-YYDOT(K,I,J)*(ABS(YPSIG(I,J))                 &
                          & +ABS(YPSIG(I,NPERIOY(J+1)))))/RDY)                                         &
                          & /(0.25d0*(ABS(YPSIG(I,NPERIOY(J+2)))+ABS(YPSIG(I,J))                       &
                          & +ABS(YPSIG(NPERIOX(I+1),NPERIOY(J+1)))                                     &
                          & +ABS(YPSIG(NPERIOX(I-1),NPERIOY(J+1)))                                     &
                          & +(TWO+TWO)*ABS(YPSIG(I,NPERIOY(J+1))))+REPS)) )
               
               ZDIV_YB    = ( (((YYDOT(K,I,NPERIOY(J+1))+YYDOT(K,I,J))*ABS(YPSIG(I,NPERIOY(J+1)))      & ! div(i,j+0.5) : option B
                          &  -(YYDOT(K,I,NPERIOY(J-1))+YYDOT(K,I,J))*ABS(YPSIG(I,J)))/RDY)             & ! interpolated fields then 
                          &  +0.25d0*(((YXDOT(K,I,J)+YXDOT(K,I,NPERIOY(J+1)))*(                        & ! compute fluxes at interfaces
                          &  ABS(YPSIG(NPERIOX(I+1),NPERIOY(J+1)))+ABS(YPSIG(I,NPERIOY(J+1)))          &
                          & +ABS(YPSIG(NPERIOX(I+1),J))+ABS(YPSIG(I,J))  )                             &
                          & -(YXDOT(K,NPERIOX(I-1),J)+YXDOT(K,NPERIOX(I-1),NPERIOY(J+1)))*(            &
                          &  ABS(YPSIG(NPERIOX(I-1),J))+ABS(YPSIG(I,NPERIOY(J+1)))                     &
                          & +ABS(YPSIG(NPERIOX(I-1),NPERIOY(J+1)))+ABS(YPSIG(I,J))))/RDX) )            &
                          & /(0.2d0*(ABS(YPSIG(NPERIOX(I-1),NPERIOY(J+1)))+ABS(YPSIG(NPERIOX(I-1),J))  &
                          &  ABS(YPSIG(NPERIOX(I+1),NPERIOY(J+1)))+ABS(YPSIG(NPERIOX(I+1),J))          &
                          & +(ONE+TWO)*(ABS(YPSIG(I,NPERIOY(J+1)))+ABS(YPSIG(I,J))))+REPS)

                ZDIV_YC   = ( (((YYDOT(K,I,NPERIOY(J+1))+YYDOT(K,I,J))*ABS(YPSIG(I,NPERIOY(J+1)))      & ! div(i,j+0.5) : option C
                          & -(YYDOT(K,I,NPERIOY(J-1))+YYDOT(K,I,J))*ABS(YPSIG(I,J)))/RDY)              & ! compute fluxes then interpolate
                          & +0.5d0*((YXDOT(K,I,J)*(ABS(YPSIG(NPERIOX(I+1),J))+ABS(YPSIG(I,J)))         & ! the resulting fluxes at interfaces
                          & +YXDOT(K,I,NPERIOY(J+1))*(ABS(YPSIG(NPERIOX(I+1),NPERIOY(J+1)))            &
                          & +ABS(YPSIG(I,NPERIOY(J+1))))-YXDOT(K,NPERIOX(I-1),J)*(                     &
                          & ABS(YPSIG(NPERIOX(I-1),J))+ABS(YPSIG(I,J)))                                &
                          & -YXDOT(K,NPERIOX(I-1),NPERIOY(J+1))*(ABS(YPSIG(I,NPERIOY(J+1)))            &
                          & +ABS(YPSIG(NPERIOX(I-1),NPERIOY(J+1)))))/RDX) )                            &
                          & /(0.2d0*(ABS(YPSIG(NPERIOX(I-1),NPERIOY(J+1)))+ABS(YPSIG(NPERIOX(I-1),J))  &
                          &  ABS(YPSIG(NPERIOX(I+1),NPERIOY(J+1)))+ABS(YPSIG(NPERIOX(I+1),J))          &
                          & +(ONE+TWO)*(ABS(YPSIG(I,NPERIOY(J+1)))+ABS(YPSIG(I,J))))+REPS)
               
               YYDOT(K+1,I,J) = ABS(YYDOT(K,I,J))*( (ABS(YPSIG(I,NPERIOY(J+1)))-ABS(YPSIG(I,J)))       &  
                          & /(ABS(YPSIG(I,J))+ABS(YPSIG(I,NPERIOY(J+1)))+REPS) )                       &
                          & -PDT*(YYDOT(K,I,J)/(PM0(I,NPERIOY(J+1))+PM0(I,J)))*(                       &
                          & (ONE-ROPT1)*ZDIV_YA + ROPT1*((ONE-ROPT2)*ZDIV_YB+ROPT2*ZDIV_YC)            &
                          & +HALF*((PM1(I,NPERIOY(J+1))-PM0(I,NPERIOY(J+1))+PM1(I,J)-PM0(I,J))/PDT))
           END DO 
         END DO
       ELSE               ! MPDATA-FD FORMULATION [Cf... Smolarkiewicz and Margolin (1996)]
         DO I=1,GXPT
            DO J=1,GYPT
   
               YXDOT(K+1,I,J) = ( ABS(YXDOT(K,I,J))-(PDT/RDX)*((YXDOT(K,I,J)**2)/                      &
                          & (HALF*(PM0(I,J)+PM0(NPERIOX(I+1),J)))) )*(                                 &
                          & (ABS(YPSIG(NPERIOX(I+1),J))-ABS(YPSIG(I,J)))                               &  
                          & /(ABS(YPSIG(I,J))+ABS(YPSIG(NPERIOX(I+1),J))+REPS) )                       &
                          & -0.25d0*(PDT/RDY)*(YXDOT(K,I,J)/(PM0(I,J)+PM0(NPERIOX(I+1),J)))*(          &
                          & (YYDOT(K,I,J)+YYDOT(K,NPERIOX(I+1),J)                                      &
                          & YYDOT(K,I,NPERIOY(J-1))+YYDOT(K,NPERIOX(I+1),NPERIOY(J-1)))*(              &
                          & (ABS(YPSIG(NPERIOX(I+1),NPERIOY(J+1)))+ABS(YPSIG(I,NPERIOY(J+1)))          &
                          & -ABS(YPSIG(NPERIOX(I+1),NPERIOY(J-1)))-ABS(YPSIG(I,NPERIOY(J-1))))         &
                          & /(ABS(YPSIG(NPERIOX(I+1),NPERIOY(J+1)))+ABS(YPSIG(I,NPERIOY(J+1)))         &
                          & +ABS(YPSIG(NPERIOX(I+1),NPERIOY(J-1)))+ABS(YPSIG(I,NPERIOY(J-1)))+REPS)))  &
                          & -HALF*PDT*(YXDOT(K,I,J)/(PM0(NPERIOX(I+1),J)+PM0(I,J)))*(                  &
                          & ((YXDOT(K,NPERIOX(I+1),J)-YXDOT(K,NPERIOX(I-1),J))/RDX)                    &
                          & +((YYDOT(K,NPERIOX(I+1),J)+YYDOT(K,I,J)                                    &
                          & -YYDOT(K,NPERIOX(I+1),NPERIOY(J-1))-YYDOT(K,I,NPERIOY(J-1)))/RDY)          &
                          & +((PM1(NPERIOX(I+1),J)-PM0(NPERIOX(I+1),J)+PM1(I,J)-PM0(I,J))/PDT) )

               YYDOT(K+1,I,J) = ( ABS(YYDOT(K,I,J))-(PDT/RDY)*((YYDOT(K,I,J)**2)/                      &
                          & (HALF*(PM0(I,J)+PM0(I,NPERIOY(J+1))))) )*(                                 &
                          & (ABS(YPSIG(I,NPERIOY(J+1)))-ABS(YPSIG(I,J)))                               &  
                          & /(ABS(YPSIG(I,J))+ABS(YPSIG(I,NPERIOY(J+1)))+REPS) )                       &
                          & -0.25d0*(PDT/RDX)*(YYDOT(K,I,J)/(PM0(I,J)+PM0(I,NPERIOY(J+1))))*(          &
                          & (YXDOT(K,I,J)+YXDOT(K,I,NPERIOY(J+1))                                      &
                          & YXDOT(K,NPERIOX(I-1),J)+YXDOT(K,NPERIOX(I-1),NPERIOX(J+1)))*(              &
                          & (ABS(YPSIG(NPERIOX(I+1),NPERIOY(J+1)))+ABS(YPSIG(NPERIO(I+1),J))           &
                          & -ABS(YPSIG(NPERIOX(I-1),NPERIOY(J+1)))-ABS(YPSIG(NPERIOX(I-1),J)))         &
                          & /(ABS(YPSIG(NPERIOX(I+1),NPERIOY(J+1)))+ABS(YPSIG(NPERIO(I+1),J))          &
                          & +ABS(YPSIG(NPERIOX(I-1),J))+ABS(YPSIG(NPERIO(I-1),NPERIOY(J+1)))+REPS)))   &
                          & -HALF*PDT*(YYDOT(K,I,J)/(PM0(I,NPERIOY(J+1))+PM0(I,J)))*(                  &
                          & ((YYDOT(K,I,NPERIOY(J+1))-YXDOT(K,I,NPERIOY(J-1)))/RDY)                    &
                          & +((YXDOT(K,I,NPERIOY(J+1))+YXDOT(K,I,J)                                    &
                          & -YXDOT(K,NPERIOX(I-1),NPERIOY(J+1))-YXDOT(K,NPERIOX(I-1),J))/RDX)          &
                          & +((PM1(I,NPERIOY(J+1))-PM0(I,NPERIOY(J+1))+PM1(I,J)-PM0(I,J))/PDT) ) 
           END DO 
         END DO
       END IF   
       !****************************************!
       !* GLOBAL FCT NON-OSCILLATORY TREATMENT *!
       !****************************************!
       IF (LFCTM) THEN
         DO I=1,GXPT
            DO J=1,GYPT

               ZMAX       = MAX(YPSIG(NPERIOX(I+1),J),YPSIG(I,NPERIOY(J+1)),YPSIG(I,J),                &
                          & YPSIG(NPERIOX(I-1),J),YPSIG(I,NPERIOY(J-1)),YPSI(0,NPERIOX(I+1),J),        & 
                          & YPSI(0,I,NPERIOY(J+1)),YPSI(0,I,J),YPSI(0,NPERIOX(I-1),J),                 &
                          & YPSI(0,I,NPERIOY(J-1)))

              ZIN(I,J)    = (ZMAX-YPSIG(I,J))/( (ONE/PM0(I,J))*(                                       & 
                          & (PDT/RDX)*(MAX(ZERO,YXDOT(K+1,NPERIOX(I-1),J))*YPSIG(NPERIOX(I-1),J)       & 
                          & -MIN(ZERO,YXDOT(K+1,I,J))*YPSIG(NPERIOX(I+1),J))                           &
                          & +(PDT/RDY)*(MAX(ZERO,YYDOT(K+1,I,NPERIOY(J-1)))*YPSIG(I,NPERIOY(J-1))      &
                          & -MIN(ZERO,YYDOT(K+1,I,J))*YPSIG(I,NPERIOY(J+1))) ) + REPS  )
              
             ZMIN         = MIN(YPSIG(NPERIOX(I+1),J),YPSIG(I,NPERIOY(J+1)),YPSIG(I,J),                &
                          & YPSIG(NPERIOX(I-1),J),YPSIG(I,NPERIOY(J-1)),YPSI(0,NPERIOX(I+1),J),        & 
                          & YPSI(0,I,NPERIOY(J+1)),YPSI(0,I,J),YPSI(0,NPERIOX(I-1),J),                 &
                          & YPSI(0,I,NPERIOY(J-1)))
              
              ZOUT(I,J)   = (YPSIG(I,J)-ZMIN)/( (ONE/PM0(I,J))*(                                       &
                          & (PDT/RDX)*(MAX(ZERO,YXDOT(K+1,I,J))*YPSIG(I,J)                             &
                          & -MIN(ZERO,YXDOT(K+1,NPERIOX(I-1),J))*YPSIG(I,J))                           &
                          & +(PDT/RDY)*(MAX(ZERO,YYDOT(K+1,I,J))*YPSIG(I,J)                            &
                          & -MIN(ZERO,YYDOT(K+1,I,NPERIOY(J-1)))*YPSIG(I,J)) ) + REPS  )
           END DO   
         END DO

         DO I=1,GXPT
           DO J=1,GYPT
              YXDOT(K+1,I,J) = MIN(ONE,ZOUT(I,J),ZIN(NPERIOX(I+1),J))*MAX(ZERO,YXDOT(K+1,I,J))         &  
                          & + MIN(ONE,ZIN(I,J),ZOUT(NPERIOX(I+1),J))*MIN(ZERO,YXDOT(K+1,I,J))
              YYDOT(K+1,I,J) = MIN(ONE,ZOUT(I,J),ZIN(I,NPERIOY(J+1)))*MAX(ZERO,YYDOT(K+1,I,J))         & 
                          & + MIN(ONE,ZIN(I,J),ZOUT(I,NPERIOY(J+1))*MIN(ZERO,YYDOT(K+1,I,J))  
           END DO   
         END DO
       END IF
       !***********************************************************************************************!
       !* CORRECTED UPWIND SCHEME *!
       !***************************!
       DO I=1,GXPT
          DO J=1,GYPT
             YPSI(K+1,I,J)=  YPSI(K,I,J)-(ONE/PM0(I,J))*(PDT/RDX)*(                                    &
                          &  (MIN(ZERO,YXDOT(K+1,I,J))*YPSIG(NPERIOX(I+1),J)                           &
                          &  +MAX(ZERO,YXDOT(K+1,I,J))*YPSIG(I,J))                                     &
                          & -(MIN(ZERO,YXDOT(K+1,NPERIOX(I-1),J))*YPSIG(I,J)                           &
                          &  +MAX(ZERO,YXDOT(K+1,NPERIOX(I-1),J))*YPSIG(NPERIOX(I-1),J)) )             &
                          & -(ONE/PM0(I,J))*(PDT/RDY)*(                                                &
                          &  (MIN(ZERO,YYDOT(K+1,I,J))*YPSIG(I,NPERIOY(J+1))                           &
                          &  +MAX(ZERO,YYDOT(K+1,I,J))*YPSIG(I,J))                                     &
                          & -(MIN(ZERO,YYDOT(K+1,I,NPERIOY(J-1)))*YPSIG(I,J)                           &
                          &  +MAX(ZERO,YYDOT(K+1,I,NPERIOY(J-1)))*YPSIG(I,NPERIOY(J-1))) )       
          END DO   
       END DO
       DO I=1,GXPT
          DO J=1,GYPT
             YPSIG(I,J)   = (PM0(I,J)/PM1(I,J))*YPSI(K+1,I,J)    
          END DO   
       END DO
     END DO

  ELSE  ! INFINITY GAUGE  OPTION
   
     DO K=1,NORD-1          
        IF (LMPDATA_FV) THEN ! MPDATA-FV FORMULATION
           DO I=1,GXPT
              DO J=1,GYPT

                 ZDIV_XA    = HALF*(((YXDOT(K,I,J)*(YPSIG(NPERIOX(I+1),J)+YPSIG(I,J))                  &
                            & -YXDOT(K,NPERIOX(I-1),J)*(YPSIG(I,J)+YPSIG(NPERIOX(I-1),J)))/RDX         &
                            & +(YYDOT(K,I,J)*(YPSIG(I,NPERIOY(J+1))+YPSIG(I,J))                        &
                            & -YYDOT(K,I,NPERIOY(J-1))*(YPSIG(I,J)+YPSIG(I,NPERIOY(J-1))))/RDY))       & !div(i+0.5,j) : option A
                            & +HALF*(((YXDOT(K,NPERIOX(I+1),J)*(YPSIG(NPERIOX(I+2),J)                  & !0.5*(div(i,j)+div(i+1,j))
                            & +YPSIG(NPERIOX(I+1),J))                                                  &
                            & -YXDOT(K,I,J)*(YPSIG(I,J)+YPSIG(NPERIOX(I+1),J)))/RDX                    &
                            & +(YYDOT(K,NPERIOX(I+1),J)*(YPSIG(NPERIOX(I+1),NPERIOY(J+1))              &
                            & +YPSIG(NPERIOX(I+1),J))                                                  &
                            & -YYDOT(K,NPERIOX(I+1),NPERIOY(J-1))*(YPSIG(NPERIOX(I+1),J)               &
                            & +YPSIG(NPERIOX(I+1),NPERIOY(J-1))))/RDY))
               
                 ZDIV_XB    = ((((YXDOT(K,NPERIOX(I+1),J)+YXDOT(K,I,J))*YPSIG(NPERIOX(I+1),J)          & ! div(i+0.5,j) : option B
                            & -(YXDOT(K,NPERIOX(I-1),J)+YXDOT(K,I,J))*YPSIG(NPERIOX(I),J))/RDX)        & ! interpolate fields --->  
                            & +0.25d0*(((YYDOT(K,I,J)+YYDOT(K,NPERIOX(I+1),J))*(                       & ! compute fluxes at interfaces
                            &  ABS(YPSIG(NPERIOX(I+1),NPERIOY(J+1)))+YPSIG(I,NPERIOY(J+1))             & ! (i,j+0.5) and (i,j-0.5)
                            & +YPSIG(NPERIOX(I+1),J)+YPSIG(I,J)  )                                     &
                            & -(YYDOT(K,I,NPERIOY(J-1))+YYDOT(K,NPERIOX(I+1),NPERIOY(J-1)))*(          &
                            &  YPSIG(NPERIOX(I+1),J)+YPSIG(I,NPERIOY(J-1))                             &
                            & +YPSIG(NPERIOX(I+1),NPERIOY(J-1))+YPSIG(I,J)))/RDY))

                 ZDIV_XC    = ( (((YXDOT(K,NPERIOX(I+1),J)+YXDOT(K,I,J))*YPSIG(NPERIOX(I+1),J)         & ! div(i+0.5,j) : option C
                            & -(YXDOT(K,NPERIOX(I-1),J)+YXDOT(K,I,J))*YPSIG(NPERIOX(I),J))/RDX)        & ! compute fluxes at (i,j+1) and (i,j) 
                            & +0.5d0*((YYDOT(K,I,J)*(YPSIG(I,NPERIOY(J+1))+YPSIG(I,J))                 & ! ---> then interpolate at interfaces
                            & +YYDOT(K,NPERIOX(I+1),J)*(YPSIG(NPERIOX(I+1),NPERIOY(J+1))               & ! (i,j+0.5) and (i,j-0.5)
                            & +YPSIG(NPERIOX(I+1),J))-YYDOT(K,I,NPERIOY(J-1))*(                        &
                            & YPSIG(I,NPERIOY(J-1))+YPSIG(I,J))                                        &
                            & -YYDOT(K,NPERIOX(I+1),NPERIOY(J-1))*(YPSIG(NPERIOX(I+1),J)               &
                            & +YPSIG(NPERIOX(I+1),NPERIOY(J-1))))/RDY) )    
                
                 YXDOT(K+1,I,J) = HALF*ABS(YXDOT(K,I,J))*(YPSIG(NPERIOX(I+1),J)-YPSIG(I,J))            &
                            & -HALF*PDT*(YXDOT(K,I,J)/(PM0(NPERIOX(I+1),J)+PM0(I,J)))*(                &
                            & (ONE-ROPT1)*ZDIV_XA +ROPT1*((ONE-ROPT2)*ZDIV_XB+ROPT2*ZDIV_XC)           &
                            & +HALF*(YPSIG(NPERIOX(I+1),J)+YPSIG(I,J))*(                               &
                            & (PM1(NPERIOX(I+1),J)-PM0(NPERIOX(I+1),J)+PM1(I,J)-PM0(I,J))/PDT) )       

                  
                  ZDIV_YA   = HALF*(((YXDOT(K,I,J)*(YPSIG(NPERIOX(I+1),J)+YPSIG(I,J))                  &
                            & -YXDOT(K,NPERIOX(I-1),J)*(YPSIG(I,J)+YPSIG(NPERIOX(I-1),J)))/RDX         &
                            & +(YYDOT(K,I,J)*(YPSIG(I,NPERIOY(J+1))+YPSIG(I,J))                        &
                            & -YYDOT(K,I,NPERIOY(J-1))*(YPSIG(I,J)+YPSIG(I,NPERIOY(J-1))))/RDY))       & !div(i,j+0.5) : option A
                            & +HALF*(((YXDOT(K,I,NPERIOX(J+1))*(YPSIG(NPERIOX(I+1),NPERIOY(J+1))       & !0.5*(div(i,j+1)+div(i,j))
                            & +YPSIG(I,NPERIOY(J+1)))-YXDOT(K,NPERIOX(I-1),NPERIOY(J+1))*(             &
                            & YPSIG(NPERIOX(I-1),NPERIOY(J+1))+YPSIG(I,NPERIOY(J+1))))/RDX             &
                            & +(YYDOT(K,I,NPERIOY(J+1))*(YPSIG(I,NPERIOY(J+1))+YPSIG(I,NPERIOY(J+2)))  &
                            & -YYDOT(K,I,J)*(YPSIG(I,J)+YPSIG(I,NPERIOY(J+1))))/RDY)) 
               
                  ZDIV_YB   = ( (((YYDOT(K,I,NPERIOY(J+1))+YYDOT(K,I,J))*YPSIG(I,NPERIOY(J+1))         & ! div(i,j+0.5) : option B
                            &  -(YYDOT(K,I,NPERIOY(J-1))+YYDOT(K,I,J))*YPSIG(I,J))/RDY)                &
                            &  +0.25d0*(((YXDOT(K,I,J)+YXDOT(K,I,NPERIOY(J+1)))*(                      &
                            &  YPSIG(NPERIOX(I+1),NPERIOY(J+1))+YPSIG(I,NPERIOY(J+1))                  &
                            & +YPSIG(NPERIOX(I+1),J)+YPSIG(I,J)  )                                     &
                            & -(YXDOT(K,NPERIOX(I-1),J)+YXDOT(K,NPERIOX(I-1),NPERIOY(J+1)))*(          &
                            &  YPSIG(NPERIOX(I-1),J)+YPSIG(I,NPERIOY(J+1))                             &
                            & +YPSIG(NPERIOX(I-1),NPERIOY(J+1))+YPSIG(I,J)))/RDX) )

                  ZDIV_YC   = ( (((YYDOT(K,I,NPERIOY(J+1))+YYDOT(K,I,J))*YPSIG(I,NPERIOY(J+1))         & ! div(i,j+0.5) : option C
                            & -(YYDOT(K,I,NPERIOY(J-1))+YYDOT(K,I,J))*YPSIG(I,J))/RDY)                 &
                            & +0.5d0*((YXDOT(K,I,J)*(YPSIG(NPERIOX(I+1),J)+YPSIG(I,J))                 &
                            & +YXDOT(K,I,NPERIOY(J+1))*(YPSIG(NPERIOX(I+1),NPERIOY(J+1))               &
                            & +YPSIG(I,NPERIOY(J+1)))-YXDOT(K,NPERIOX(I-1),J)*(YPSIG(NPERIOX(I-1),J)   &
                            & +YPSIG(I,J))-YXDOT(K,NPERIOX(I-1),NPERIOY(J+1))*(YPSIG(I,NPERIOY(J+1))   &
                            & +YPSIG(NPERIOX(I-1),NPERIOY(J+1))))/RDX) )                                 
                         
                  YYDOT(K+1,I,J) = HALF*ABS(YYDOT(K,I,J))*(YPSIG(I,NPERIOY(J+1))-YPSIG(I,J))           &
                            & -HALF*PDT*(YYDOT(K,I,J)/(PM0(I,NPERIOY(J+1))+PM0(I,J)))*(                &
                            & (ONE-ROPT1)*ZDIV_YA +ROPT1*((ONE-ROPT2)*ZDIV_YB+ROPT2*ZDIV_YC)           &
                            & +HALF*(YPSIG(I,NPERIOY(J+1))+YPSIG(I,J))*(                               &
                            & (PM1(I,NPERIOY(J+1))-PM0(I,NPERIOY(J+1))+PM1(I,J)-PM0(I,J))/PDT) )      
               END DO        
            END DO
         ELSE                ! MPDATA-FD FORMULATION [Cf Smolarkiewicz and Margolin 1996]
            DO I=1,GXPT
               DO J=1,GYPT
                  YXDOT(K+1,I,J) = HALF*( ABS(YXDOT(K,I,J))-(PDT/RDX)*((YXDOT(K,I,J)**2)/              &
                            & (HALF*(PM0(I,J)+PM0(NPERIOX(I+1),J)))) )*(                               &
                            & (YPSIG(NPERIOX(I+1),J)-YPSIG(I,J)) )                                     &
                            & -0.25d0*(PDT/RDY)*(YXDOT(K,I,J)/(PM0(I,J)+PM0(NPERIOX(I+1),J)))*(        &
                            & (YYDOT(K,I,J)+YYDOT(K,NPERIOX(I+1),J)                                    &
                            & +YYDOT(K,I,NPERIOY(J-1))+YYDOT(K,NPERIOX(I+1),NPERIOY(J-1)))*(           &
                            & 0.25d0*( (YPSIG(NPERIOX(I+1),NPERIOY(J+1))+YPSIG(I,NPERIOY(J+1))         &
                            & -YPSIG(NPERIOX(I+1),NPERIOY(J-1))-YPSIG(I,NPERIOY(J-1))))))              &
                            & -HALF*PDT*YXDOT(K,I,J)*(HALF*(YPSIG(NPERIOX(I+1),J)+YPSIG(I,J)))*(       &
                            & ((YXDOT(K,NPERIOX(I+1),J)-YXDOT(K,NPERIOX(I-1),J))/RDX)                  &
                            & +((YYDOT(K,NPERIOX(I+1),J)+YYDOT(K,I,J)                                  &
                            & -YYDOT(K,NPERIOX(I+1),NPERIOY(J-1))-YYDOT(K,I,NPERIOY(J-1)))/RDY)        &
                            & +((PM1(NPERIOX(I+1),J)-PM0(NPERIOX(I+1),J)+PM1(I,J)-PM0(I,J))/PDT) )     &
                            & /(PM0(NPERIOX(I+1),J)+PM0(I,J)) 

               YYDOT(K+1,I,J) = HALF*( ABS(YYDOT(K,I,J))-(PDT/RDY)*((YYDOT(K,I,J)**2)/                 &
                            & (HALF*(PM0(I,J)+PM0(I,NPERIOY(J+1))))) )*(                               &
                            & (YPSIG(I,NPERIOY(J+1))-YPSIG(I,J)) )                                     &
                            & -0.25d0*(PDT/RDX)*(YYDOT(K,I,J)/(PM0(I,J)+PM0(I,NPERIOY(J+1))))*(        &
                            & (YXDOT(K,I,J)+YXDOT(K,I,NPERIOY(J+1))                                    &
                            & +YXDOT(K,NPERIOX(I-1),J)+YXDOT(K,NPERIOX(I-1),NPERIOX(J+1)))*(           &
                            & 0.25d0*((YPSIG(NPERIOX(I+1),NPERIOY(J+1))+YPSIG(NPERIO(I+1),J)           &
                            & -YPSIG(NPERIOX(I-1),NPERIOY(J+1))-YPSIG(NPERIOX(I-1),J)))))              &
                            & -HALF*PDT*YYDOT(K,I,J)*(HALF*(YPSIG(I,NPERIOY(J+1))+YPSIG(I,J)))*(       &
                            & ((YYDOT(K,I,NPERIOY(J+1))-YXDOT(K,I,NPERIOY(J-1)))/RDY)                  &
                            & +((YXDOT(K,I,NPERIOY(J+1))+YXDOT(K,I,J)                                  &
                            & -YXDOT(K,NPERIOX(I-1),NPERIOY(J+1))-YXDOT(K,NPERIOX(I-1),J))/RDX)        &
                            & +((PM1(I,NPERIOY(J+1))-PM0(I,NPERIOY(J+1))+PM1(I,J)-PM0(I,J))/PDT) )     &
                            & /(PM0(I,NPERIOY(J+1))+PM0(I,J)) 
               END DO 
            END DO
         END IF   
         !* GLOBAL FCT NON-OSCILLATORY TREATMENT *!
         IF (LFCTM) THEN
            DO I=1,GXPT
               DO J=1,GYPT
              

                  ZMAX      = MAX(YPSIG(NPERIOX(I+1),J),YPSIG(I,NPERIOY(J+1)),YPSIG(I,J),              &
                            & YPSIG(NPERIOX(I-1),J),YPSIG(I,NPERIOY(J-1)),YPSI(0,NPERIOX(I+1),J),      & 
                            & YPSI(0,I,NPERIOY(J+1)),YPSI(0,I,J),YPSI(0,NPERIOX(I-1),J),               &
                            & YPSI(0,I,NPERIOY(J-1)))
                 
                  ZIN(I,J)  = (ZMAX-YPSIG(I,J))/( (ONE/PM0(I,J))*(                                     & 
                            & (PDT/RDX)*(MAX(ZERO,YXDOT(K+1,NPERIOX(I-1),J))                           &
                            & -MIN(ZERO,YXDOT(K+1,I,J)))                                               &
                            & +(PDT/RDY)*(MAX(ZERO,YYDOT(K+1,I,NPERIOY(J-1)))                          &
                            & -MIN(ZERO,YYDOT(K+1,I,J))) ) + REPS  )
              
                  ZMIN      = MIN(YPSIG(NPERIOX(I+1),J),YPSIG(I,NPERIOY(J+1)),YPSIG(I,J),              &
                            & YPSIG(NPERIOX(I-1),J),YPSIG(I,NPERIOY(J-1)),YPSI(0,NPERIOX(I+1),J),      & 
                            & YPSI(0,I,NPERIOY(J+1)),YPSI(0,I,J),YPSI(0,NPERIOX(I-1),J),               &
                            & YPSI(0,I,NPERIOY(J-1)))
              
                  ZOUT(I,J) = (YPSIG(I,J)-ZMIN)/( (ONE/PM0(I,J))*(                                     &
                            & (PDT/RDX)*(MAX(ZERO,YXDOT(K+1,I,J))                                      &
                            & -MIN(ZERO,YXDOT(K+1,NPERIOX(I-1),J)))                                    &
                            & +(PDT/RDY)*(MAX(ZERO,YYDOT(K+1,I,J))                                     &
                            & -MIN(ZERO,YYDOT(K+1,I,NPERIOY(J-1)))) ) + REPS  )
                END DO   
             END DO

             DO I=1,GXPT
                DO J=1,GYPT
                   YXDOT(K+1,I,J) = MIN(ONE,ZOUT(I,J),ZIN(NPEIOX(I+1),J))*MAX(ZERO,YXDOT(K+1,I,J))     & 
                            & +MIN(ONE,ZIN(I,J),ZOUT(NPERIOX(I+1),J))*MIN(ZERO,YXDOT(K+1,I,J))
                   YYDOT(K+1,I,J)  = MIN(ONE,ZOUT(I,J),ZIN(I,NPERIOY(J+1)))*MAX(ZERO,YYDOT(K+1,I,J))   & 
                            & +MIN(ONE,ZIN(I,J),ZOUT(I,NPERIOY(J+1)))*MIN(ZERO,YYDOT(K+1,I,J))
                END DO   
             END DO   
         END IF
         !*********************************************************************************************!
         !* CORRECTED UPWIND SCHEME *!
         !***************************!
         DO I=1,GXPT
            DO J=1,GYPT  
               YPSI(K+1,I,J) = YPSI(K,I,J) - (ONE/PM0(I,J))*(                                          &
                            &  (PDT/RDX)*(YXDOT(K+1,I,J)-YXDOT(K+1,NPERIOX(I-1),J))                    &
                            & +(PDT/RDY)*(YYDOT(K+1,I,J)-YYDOT(K+1,I,NPERIOY(J-1)))  )           
               YPSIG(I,J)   = (PM0(I,J)/PM1(I,J))*YPSI(K+1,I,J)  
            END DO   
         END DO
     END DO
  END IF

 
  !****************************************************************************************************!
  DO I=1,GXPT
     DO J=1,GYPT
        PPSI_DEP(I,J)      = YPSIG(I,J) + ZPSI_MIN
        IF (LLMASS) THEN
           PXDOT(I,J)      = ZERO
           PYDOT(I,J)      = ZERO
           DO K=1,NORD2
              PXDOT(I,J)   = PXDOT(I,J) + HALF*YXDOT(K,I,J)*(                                          &
                           & (YPSI(K-1,NPERIOX(I+1),J)+YPSI(K-1,I,J)) )
              PYDOT(I,J)   = PYDOT(I,J) + HALF*YYDOT(K,I,J)*(                                          &
                           & (YPSI(K-1,I,NPERIOY(J+1))+YPSI(K-1,I,J)) )
           END DO   
        END IF 
     END DO   
  END DO 
  
END SUBROUTINE MPDATA_ADV_XY
!#########################################################################################################!
!#########################################################################################################!
!#########################################################################################################!
!#########################################################################################################!
!#########################################################################################################!     
!#########################################################################################################!
!#########################################################################################################!
!#########################################################################################################!
!#########################################################################################################!
END MODULE MOD_MPDTA
