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

  USE MOD_SHARE, ONLY : NXPT,NLEV,RR,RG,RPI,ONE,TWO,ZERO,HALF,RDT,RDX,     &
                      & NPERIO,LPERIO,NORD1,NORD2,REPS,RFCT,NBC,NBL,RVMAX, &
                      & LPRINTLEV,LADV_SPLIT,LADV_PIS,LSETTLS,LFCT_MONO,   &
                      & LMPDATA_FV
  USE MOD_SETUP, ONLY : MPDATA_STRUCT,RHSVAR,GEOMETRY
 
CONTAINS

  !******************************************************************************************************!
  !******************************************************************************************************!
  !******************************************************************************************************!
  !******************************************************************************************************!
  !******************************************************************************************************!
  SUBROUTINE MPDATA_TRACER_TRANSPORT_SCHEME(RHS_DEP,RHS_ARR,MPDTA,PMF0, &
             & PXDOTF1,PXDOTF0,PXDOTF9,PZDOTH1,PZDOTH0,PZDOTH9, &
             & PDETAF,PDETAH,PDELTB,KSTEP,CPART)
  
  IMPLICIT NONE

  TYPE(RHSVAR),         INTENT(OUT)    :: RHS_DEP
  TYPE(RHSVAR),         INTENT(IN)     :: RHS_ARR
  TYPE(MPDATA_STRUCT),  INTENT(INOUT)  :: MPDTA
  REAL(8),              INTENT(IN)     :: PMF0(NLEV,NXPT)  
  REAL(8),              INTENT(IN)     :: PXDOTF0(NLEV,NXPT),PXDOTF9(NLEV,NXPT),PXDOTF1(NLEV,NXPT)
  REAL(8),              INTENT(IN)     :: PZDOTH0(0:NLEV,NXPT),PZDOTH9(0:NLEV,NXPT),PZDOTH1(0:NLEV,NXPT)
  REAL(8),              INTENT(IN)     :: PDELTB(NLEV),PDETAF(NLEV),PDETAH(0:NLEV)
  INTEGER(8),           INTENT(IN)     :: KSTEP
  CHARACTER(LEN=1),     INTENT(IN)     :: CPART

  INTEGER(8)                           :: I,J,K,IS
  REAL(8)                              :: ZLIN,ZVMAX,ZWIND
  

  ZLIN = HALF
  IF (CPART == 'LI') ZLIN =ZERO

  IS = 1
  IF (LADV_SPLIT) IS = 3
  

  !**************************************!
  !* Initialize vertical metric (mf,mh) *!
  !**************************************!
  
  DO I=1,NXPT
     MPDTA%MF(0,:,I) = PMF0(:,I)
  END DO
  
  !**************************************!
  !* Advecting velocity at (t+0.5*dt)   *!
  !**************************************!

  IF (LSETTLS) THEN
     DO I=1,NXPT
        MPDTA%XDOTF(:,I)   =  ZLIN*HALF*( (3.d0*PXDOTF0(:,I)-PXDOTF9(:,I)) &
                           &  + (3.d0*PXDOTF0(:,NPERIO(I+1))-PXDOTF9(:,NPERIO(I+1))) )
        MPDTA%ZDOTH(1,:,I) =  ZLIN*(3.d0*PZDOTH0(:,I)-PZDOTH9(:,I))
     END DO 
  ELSE
     DO I=1,NXPT
        MPDTA%XDOTF(:,I)   =  ZLIN*HALF*((PXDOTF1(:,I)+PXDOTF0(:,I)) &
                           &  + (PXDOTF1(:,NPERIO(I+1))+PXDOTF0(:,NPERIO(I+1))))
        MPDTA%ZDOTH(1,:,I) =  ZLIN*(PZDOTH1(:,I)+PZDOTH0(:,I))
     END DO
  ENDIF

  !********************************!
  !* special treatment at surface *!
  !********************************!
  
  DO I=1,NXPT     
     MPDTA%XDOTS(I)        = ZERO
     DO J=1,NLEV
        MPDTA%XDOTS(I)     = MPDTA%XDOTS(I) + PDELTB(J)*MPDTA%XDOTF(J,I) 
     END DO
  END DO

  !*****************************!
  !*  masses transport scheme  *!
  !*****************************!
  
  IF (LADV_SPLIT) THEN 
     CALL MPDATA_SPLIT_2D_FULL_MASS(MPDTA%MF,    &
          & MPDTA%XDOTF,MPDTA%ZDOTH,PDETAF) 
  ELSE
     CALL MPDATA_UNSPLIT_2D_FULL_MASS(MPDTA%MF,  &
          & MPDTA%XDOTF,MPDTA%ZDOTH(1,:,:),PDETAF)
  END IF

  !************************************************************************!
  IF (LPRINTLEV) THEN                                                      !
     WRITE(*,*) 'NSTEP : ',KSTEP,                        &                 !
              & ' MAX.WIND.HOR : ',MAXVAL(ABS(PXDOTF0)), &                 !
              & ' MAX.WIND.VER : ',MAXVAL(ABS(PZDOTH0))                    !
  END IF                                                                   !
  IF (MAXVAL(PXDOTF0) .GE. RVMAX) THEN                                     !
     DO I=1,NXPT                                                           !
        DO J=1,NLEV                                                        !
             ZWIND=PXDOTF0(J,I)                                            !
             IF (ZWIND .GE. RVMAX) THEN                                    !
               WRITE (*,*) '(I,J) = ','(',I,',',J,')',' .... ',ZWIND       !
            END IF                                                         !
         END DO                                                            !
      END DO                                                               !
     WRITE(*,*) 'WIND TOO STRONG EXPLOSION AT TIMESTEP = ',KSTEP           !
     STOP                                                                  !
  END IF                                                                   !
  !************************************************************************!

  
  CALL MPDATA_SURF_PRESSURE(RHS_DEP%PIS,RHS_ARR%PIS,MPDTA%XDOTS)

  
  IF (LADV_SPLIT) THEN
     CALL MPDATA_SPLIT_2D_FULL_FIELD(RHS_DEP%Q,RHS_ARR%Q,      &
          & MPDTA%XDOTF,MPDTA%ZDOTH,MPDTA%MF,PDETAF)
  ELSE
     CALL MPDATA_UNSPLIT_2D_FULL_FIELD(RHS_DEP%Q,RHS_ARR%Q,    &
          & MPDTA%XDOTF,MPDTA%ZDOTH(1,:,:),MPDTA%MF,PDETAF)
  END IF   

  
END SUBROUTINE MPDATA_TRACER_TRANSPORT_SCHEME
!############################################################################################!
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!############################################################################################!
SUBROUTINE MPDATA_SPLIT_2D_MASS(PMF,PXDOTF,PZDOTH,PDETAF)

  IMPLICIT NONE

  REAL(8),       INTENT(INOUT)  :: PMF(0:3,NLEV,NXPT)
  REAL(8),       INTENT(INOUT)  :: PXDOTF(NLEV,NXPT)
  REAL(8),       INTENT(INOUT)  :: PZDOTH(3,0:NLEV,NXPT)
  REAL(8),       INTENT(IN)     :: PDETAF(NLEV)

  REAL(8)                       :: ZUN(NLEV,NXPT)

  ZUN(:,:) = ONE
  
  CALL ADV_Z(PMF(1,:,:),PMF(0,:,:),PZDOTH,     &
       & ZUN,ZUN,PDETAF,RDT/TWO,1,LDMASS=.TRUE.)
  CALL ADV_X(PMF(2,:,:),PMF(1,:,:),PXDOTF,     &
       & ZUN,ZUN,RDT,1_8,LDMASS=.TRUE.)
  CALL ADV_Z(PMF(3,:,:),PMF(2,:,:),PZDOTH,     &
       & ZUN,ZUN,PDETAF,RDT/TWO,2,LDMASS=.TRUE.)     
  
END SUBROUTINE MPDATA_SPLIT_2D_MASS
!################################################################################################!
!************************************************************************************************!
!************************************************************************************************!
!************************************************************************************************!
!************************************************************************************************!
!################################################################################################!
SUBROUTINE MPDATA_SPLIT_2D_FIELD(PPSI_DEP,PPSI_ARR,PXDOTF,PZDOTH,PMF,PDETAF)

  IMPLICIT NONE

  REAL(8),       INTENT(OUT)    :: PPSI_DEP(NLEV,NXPT)
  REAL(8),       INTENT(IN)     :: PPSI_ARR(NLEV,NXPT)
  REAL(8),       INTENT(INOUT)  :: PXDOTF(NLEV,NXPT)
  REAL(8),       INTENT(INOUT)  :: PZDOTH(3,0:NLEV,NXPT)
  REAL(8),       INTENT(IN)     :: PMF(0:3,NLEV,NXPT)
  REAL(8),       INTENT(IN)     :: PDETAF(NLEV)
  
  REAL(8)                       :: YPSI(3,NLEV,NXPT)
  INTEGER(8)                    :: I,J

  CALL ADV_Z(YPSI(1,:,:),PPSI_ARR,PZDOTH,      &
       & PMF(1,:,:),PMF(0,:,:),PDETAF,RDT/TWO,1)
  CALL ADV_X(YPSI(2,:,:),YPSI(1,:,:),PXDOTF,   &
       & PMF(2,:,:),PMF(1,:,:),RDT,1_8)
  CALL ADV_Z(YPSI(3,:,:),YPSI(2,:,:),PZDOTH,   &
       & PMF(3,:,:),PMF(2,:,:),PDETAF,RDT/TWO,2)

  DO I=1,NXPT
     DO J=1,NLEV
        PPSI_DEP(J,I) = (PMF(0,J,I)/PMF(3,J,I))*YPSI(3,J,I)
     END DO
  END DO

END SUBROUTINE  MPDATA_SPLIT_2D_FIELD
!################################################################################################!
!************************************************************************************************!
!************************************************************************************************!
!************************************************************************************************!
!************************************************************************************************!
!************************************************************************************************!
!************************************************************************************************!
!************************************************************************************************!
!************************************************************************************************!
!************************************************************************************************!
!************************************************************************************************!
!################################################################################################!
SUBROUTINE ADV_X(PPSI_DEP,PPSI_ARR,PXDOT,PM1,PM0,PDT,PPSI_MIN,KTOP,LDMASS)
  
  IMPLICIT NONE

  REAL(8),                       INTENT(OUT)    :: PPSI_DEP(KTOP:NLEV,NXPT)
  REAL(8),                       INTENT(IN)     :: PPSI_ARR(KTOP:NLEV,NXPT)
  REAL(8),                       INTENT(INOUT)  :: PXDOT(KTOP:NLEV,NXPT)
  REAL(8),                       INTENT(IN)     :: PM1(KTOP:NLEV,NXPT)
  REAL(8),                       INTENT(IN)     :: PM0(KTOP:NLEV,NXPT)
  REAL(8),                       INTENT(IN)     :: PDT,PPSI_MIN
  INTEGER(8),                    INTENT(IN)     :: KTOP
  LOGICAL,OPTIONAL,              INTENT(IN)     :: LDMASS
  
  REAL(8)                                       :: YPSI(0:NORD2,KTOP:NLEV,NXPT)
  REAL(8)                                       :: YXDOT(NORD2,KTOP:NLEV,NXPT)
  REAL(8)                                       :: YPSIG(KTOP:NLEV,NXPT)
  REAL(8)                                       :: ZIN(NLEV,NXPT),ZOUT(NLEV,NXPT)
  REAL(8)                                       :: ZPSI_MIN,ZXDOT_MONO
  INTEGER(8)                                    :: I,J,K
  LOGICAL                                       :: LLMASS                               

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IF (PRESENT(LDMASS)) THEN
     LLMASS = LDMASS
  ELSE
     LLMASS = .FALSE.
  END IF

  IF (LLMASS) THEN
     ZPSI_MIN = ZERO
  ELSE   
     ZPSI_MIN = MIN(ZERO,MIN(PPSI_MIN,MINVAL(PPSI_ARR)))
  END IF

  DO I=1,NXPT
     YPSI(0,:,I)  = PPSI_ARR(:,I) - ZPSI_MIN
     YXDOT(1,:,I) = PXDOT(:,I)
  END DO
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  
  !***************************************!
  ! DONOR-CELL UPWIND FIRST-ORDER SCHEME *!
  !***************************************!
  DO I=1,NXPT
     DO J=KTOP,NLEV
        YPSI(1,J,I)          = YPSI(0,J,I) - (ONE/PM0(J,I))*(PDT/RDX)*(                           &
                             & (MIN(ZERO,YXDOT(1,J,I))*YPSI(0,J,NPERIO(I+1))                      &
                             & +MAX(ZERO,YXDOT(1,J,I))*YPSI(0,J,I))                               &
                             & -(MIN(ZERO,YXDOT(1,J,NPERIO(I-1)))*YPSI(0,J,I)                     &
                             & +MAX(ZERO,YXDOT(1,J,NPERIO(I-1)))*YPSI(0,J,NPERIO(I-1))) )
        YPSIG(J,I)           = (PM0(J,I)/PM1(J,I))*YPSI(1,J,I)
     END DO   
  END DO
  
  !************************!
  ! FIRST CORRECTION-STEP *!
  !************************!
  IF (.NOT.LGAUGE) THEN
    DO J=KTOP,NLEV
      DO K=1,NORD-1

        !*****************************************************************************************!
        !* COMPUTE SURFACE PSEUDO-VELOCITIES
        IF (LMPDATA_FV) THEN
           DO I=1,NXPT
              YXDOT(K+1,J,I) = ABS(YXDOT(K,J,I))*( (ABS(YPSIG(J,NPERIO(I+1)))-ABS(YPSIG(J,I)))    &
                             & /(ABS(YPSIG(J,I))+ABS(YPSIG(J,NPERIO(I+1)))+REPS) )                &
                             & -PDT*(YXDOT(K,J,I)/(PM0(J,NPERIO(I+1))+PM0(J,I)))*(                &
                             & (((((YXDOT(K,J,I)+YXDOT(K,J,NPERIO(I+1)))*ABS(YPSIG(J,NPERIO(I+1)))&
                             & -(YXDOT(K,J,NPERIO(I-1))+YXDOT(K,J,I))*ABS(YPSIG(J,I))))/RDX)      &
                             & /((ABS(YPSIG(J,NPERIO(I+1)))+ABS(YPSIG(J,I)))+REPS))               &
                             & +(HALF*(PM1(J,NPERIO(I+1))-PM0(J,NPERIO(I+1))                      &
                             & +PM1(J,I)-PM0(J,I))/PDT) )
           END DO   
        ELSE
           DO I=1,NXPT
              YXDOT(K+1,J,I) = (ABS(YXDOT(K,J,I))-(PDT/RDX)*((YXDOT(K,J,I)**2)                    &
                             & /(HALF*(PM0(J,NPERIO(I+1))+PM0(J,I)))))*(                          &
                             & (ABS(YPSIG(J,NPERIO(I+1)))-ABS(YPSIG(J,I)))                        &
                             & /(ABS(YPSIG(J,I))+ABS(YPSIG(J,NPERIO(I+1)))+REPS) )                &
                             & -HALF*PDT*(YXDOT(K,J,I)/(PM0(J,NPERIO(I+1))+PM0(J,I)))*(           &
                             & ((YXDOT(K,J,NPERIO(I+1))-YXDOT(K,J,NPERIO(I-1)))/RDX)              &
                             & +((PM1(J,NPERIO(I+1))-PM0(J,NPERIO(I+1))+PM1(J,I)-PM0(J,I))/PDT) )
           END DO
        END IF   
        !*****************************************************************************************!
        !* GLOBAL FCT NON-OSCILLATORY TREATMENT
        IF (LFCT_MONO) THEN 
           DO I=1,NXPT
              ZIN(J,I)       = ( MAX( YPSIG(J,NPERIO(I+1)),YPSIG(J,NPERIO(I)),                    &
                             & YPSIG(J,NPERIO(I-1)),YPSI(0,J,NPERIO(I+1)),                        &
                             & YPSI(0,J,I),YPSI(0,J,NPERIO(I-1)))-YPSIG(J,I) )                    &
                             & /( (ONE/PM0(J,I))*(PDT/RDX)*(                                      &
                             & MAX(ZERO,YXDOT(K+1,J,NPERIO(I-1)))*YPSI(K,J,NPERIO(I-1))           &
                             & -MIN(ZERO,YXDOT(K+1,J,I))*YPSI(K,J,NPERIO(I+1)))+REPS  )
              
              ZOUT(J,I)      = -( MIN( YPSIG(J,NPERIO(I+1)),YPSIG(J,I),                           &
                             & YPSIG(J,NPERIO(I-1)),YPSI(0,J,NPERIO(I+1)),                        &
                             & YPSI(0,J,I),YPSI(0,J,NPERIO(I-1)))-YPSIG(J,I) )                    &
                             & /( (ONE/PM0(J,I))*(PDT/RDX)*(                                      &
                             & MAX(ZERO,YXDOT(K+1,J,I))*YPSIG(J,I)                                &
                             & -MIN(ZERO,YXDOT(K+1,J,NPERIO(I-1)))*YPSIG(J,I))+REPS)
           END DO
           DO I=1,NXPT
              YXDOT(K+1,J,I) =  MIN(ONE,ZOUT(J,I),ZIN(J,NPERIO(I+1)))*MAX(ZERO,YXDOT(K+1,J,I))    &
                             & +MIN(ONE,ZIN(J,I),ZOUT(J,NPERIO(I+1)))*MIN(ZERO,YXDOT(K+1,J,I))
           END DO
        END IF   
        !*****************************************************************************************!
        !* CORRECTED UPWIND SCHEME
        DO I=1,NXPT
           YPSI(K+1,J,I)     = YPSI(K,J,I) - (ONE/PM0(J,I))*(PDT/RDX)*(                           &
                             &  (MIN(ZERO,YXDOT(K+1,J,I))*YPSIG(J,NPERIO(I+1))                    &
                             & +MAX(ZERO,YXDOT(K+1,J,NPERIO(I)))*YPSIG(J,I))                      &
                             & -(MIN(ZERO,YXDOT(K+1,J,NPERIO(I-1)))*YPSIG(J,I)                    &
                             & +MAX(ZERO,YXDOT(K+1,J,NPERIO(I-1)))*YPSIG(J,NPERIO(I-1))) )
        END DO
        DO I=1,NXPT
           YPSIG(J,I)        = (PM0(J,I)/PM1(J,I))*YPSI(K+1,J,I)
        END DO
      END DO
    END DO

  ELSE !INFINITY GAUGE
    
    !*************************!
    ! SECOND CORRECTION-STEP *!
    !*************************!
    DO J=KTOP,NLEV 
       DO K=1,NORD-1
          !***************************************************************************************!
          !* COMPUTE SURFACE PSEUDO-VELOCITIES
          IF (LMPDATA_FV) THEN
            DO I=1,NXPT
               YXDOT(K+1,J,I) = HALF*ABS(YXDOT(K,J,I))*( (YPSIG(J,NPERIO(I+1))-YPSIG(J,I)) )      &
                              & -HALF*PDT*(YXDOT(K,J,I)/(PM0(J,NPERIO(I+1))+PM0(J,I)))*(          &
                              & (((((YXDOT(K,J,I)+YXDOT(K,J,NPERIO(I+1)))*YPSIG(J,NPERIO(I+1))    &
                              & -(YXDOT(K,J,NPERIO(I-1))+YXDOT(K,J,I))*YPSIG(J,I)))/RDX))         &
                              & +(HALF*(YPSIG(J,NPERIO(I+1))+YPSIG(J,I))*(                        &
                              & PM1(J,NPERIO(I+1))-PM0(J,NPERIO(I+1))+PM1(J,I)-PM0(J,I))/PDT) )
            END DO   
          ELSE
            DO I=1,NXPT
               YXDOT(K+1,J,I) = HALF*(ABS(YXDOT(K,J,I))-(PDT/RDX)*((YXDOT(K,J,I)**2)              &
                              & /(HALF*(PM0(J,NPERIO(I+1))+PM0(J,I)))))*(                         &
                              & (YPSIG(J,NPERIO(I+1))-YPSIG(J,I)) )                               &
                              & -HALF*PDT*(YXDOT(K,J,I)/(PM0(J,NPERIO(I+1))+PM0(J,I)))*(          &
                              & HALF*(YPSIG(J,NPERIO(I+1))+YPSIG(J,I))*(                          &
                              & ((YXDOT(K,J,NPERIO(I+1))-YXDOT(K,J,NPERIO(I-1)))/RDX)             &
                              & +((PM1(J,NPERIO(I+1))-PM0(J,NPERIO(I+1))+PM1(J,I)-PM0(J,I))/PDT)))
            END DO
          END IF   
          !***************************************************************************************!
          !* GLOBAL FCT NON-OSCILLATORY TREATMENT
          IF (LFCT_MONO) THEN 
             DO I=1,NXPT
        
                ZIN(J,I)      = ( MAX( YPSIG(J,NPERIO(I+1)),YPSIG(J,NPERIO(I)),                   &
                              & YPSIG(J,NPERIO(I-1)),YPSI(0,J,NPERIO(I+1)),                       &
                              & YPSI(0,J,I),YPSI(0,J,NPERIO(I-1)))-YPSIG(J,I) )                   &
                              & /( (ONE/PM0(J,I))*(PDT/RDX)*(                                     &
                              & MAX(ZERO,YXDOT(K+1,J,NPERIO(I-1)))                                &
                              & -MIN(ZERO,YXDOT(K+1,J,I)))+REPS  )
               
                ZOUT(J,I)     = -( MIN( YPSIG(J,NPERIO(I+1)),YPSIG(J,I),                          &
                              & YPSIG(J,NPERIO(I-1)),YPSI(0,J,NPERIO(I+1)),                       &
                              & YPSI(0,J,I),YPSI(0,J,NPERIO(I-1)))-YPSIG(J,I) )                   &
                              & /( (ONE/PM0(J,I))*(PDT/RDX)*(                                     &
                              & MAX(ZERO,YXDOT(K+1,J,I))                                          &
                              & -MIN(ZERO,YXDOT(K+1,J,NPERIO(I-1))))+REPS  )
             END DO
             DO I=1,NXPT
                YXDOT(K+1,J,I) = MIN(ONE,ZOUT(J,I),ZIN(J,NPERIO(I+1)))*MAX(ZERO,YXDOT(K+1,J,I))   &
                               & +MIN(ONE,ZIN(J,I),ZOUT(J,NPERIO(I+1)))*MIN(ZERO,YXDOT(K+1,J,I))
             END DO
          END IF                                                      
          !***************************************************************************************!
          !* CORRECTED UPWIND SCHEME
          DO I=1,NXPT   
             YPSI(K+1,J,I)     = YPSI(K,J,I)-(ONE/PM0(J,I))*(PDT/RDX)*(                           &
                               & YXDOT(K+1,J,I)-YXDOT(K+1,J,NPERIO(I-1))  )
             YPSIG(J,I)        = (PM0(J,I)/PM1(J,I))*YPSI(K+1,J,I)
          END DO
       END DO
    END DO
  END IF  
  !***********************************************************************************************!

  DO I=1,NXPT
     PPSI_DEP(:,I)     = YPSIG(:,I)+ZPSI_MIN   
     IF (LLMASS) THEN
        PXDOT(:,I)     = ZERO
        DO K=1,NORD
           PXDOT(:,I)  = PXDOT(:,I) + HALF*YXDOT(K,:,I)*(                                    &
                       & (YPSI(K-1,:,NPERIO(I+1))+YPSI(K-1,:,I)) )
        END DO
     END IF   
  END DO
  
END SUBROUTINE ADV_X
!################################################################################################!
!************************************************************************************************!
!************************************************************************************************!
!************************************************************************************************!
!************************************************************************************************!
!################################################################################################!
SUBROUTINE ADV_ZF(PPSI_DEP,PPSI_ARR,PZDOT,PM1,PM0,PDETAF,PDETAH,PDT,KS,LDMASS)
  
  IMPLICIT NONE

  REAL(8),                       INTENT(OUT)    :: PPSI_DEP(NLEV,NXPT)
  REAL(8),                       INTENT(IN)     :: PPSI_ARR(NLEV,NXPT)
  REAL(8),                       INTENT(INOUT)  :: PZDOT(3,0:NLEV,NXPT)
  REAL(8),                       INTENT(IN)     :: PM1(NLEV,NXPT),PM0(NLEV,NXPT)
  REAL(8),                       INTENT(IN)     :: PDETAF(NLEV),PDETAH(0:NLEV)
  REAL(8),                       INTENT(IN)     :: PDT
  INTEGER,                       INTENT(IN)     :: KS
  LOGICAL,OPTIONAL,              INTENT(IN)     :: LDMASS

  REAL(8)                                       :: YPSI(0:NORD2,NLEV,NXPT)
  REAL(8)                                       :: YZDOT(NORD2,0:NLEV,NXPT)
  REAL(8)                                       :: YPSIG(NLEV,NXPT)
  REAL(8)                                       :: ZIN(NLEV,NXPT),ZOUT(NLEV,NXPT)
  REAL(8)                                       :: ZPSI_MIN,ZZDOT_MONO
  INTEGER(8)                                    :: I,J,K
  LOGICAL                                       :: LLMASS

  IF (PRESENT(LDMASS)) THEN
     LLMASS = LDMASS
  ELSE
     LLMASS = .FALSE.
  END IF   
  
  ZPSI_MIN = MINVAL(PPSI_ARR)

  DO I=1,NXPT
     YPSI(0,:,I)  = PPSI_ARR(:,I) - ZPSI_MIN
     YZDOT(1,:,I) = PZDOT(KS,:,I) 
  END DO

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !***************************************!
  ! DONOR-CELL UPWIND FIRST-ORDER SCHEME *!
  !***************************************!
  DO I=1,NXPT
     DO J=1,NLEV
        YPSI(1,J,I)          = YPSI(0,J,I) - (ONE/PM0(J,I))*(PDT/PDETAF(J))*(                 &
                             & (MIN(ZERO,YZDOT(1,J,I))*YPSI(0,NBC(J+1),I)                     &
                             & +MAX(ZERO,YZDOT(1,J,I))*YPSI(0,J,I))                           &
                             & -(MIN(ZERO,YZDOT(1,J-1,I))*YPSI(0,J,I)                         &
                             & +MAX(ZERO,YZDOT(1,J-1,I))*YPSI(0,NBC(J-1),I)) )
        YPSIG(J,I)           = (PM0(J,I)/PM1(J,I))*YPSI(1,J,I)
     END DO   
  END DO
  
  !************************!
  ! FIRST CORRECTION-STEP *!
  !************************!
  IF (.NOT.LGAUGE) THEN

    DO I=1,NXPT 
      DO K=1,NORD-1
        
        !*************************************************************************************!
        !* COMPUTE SURFACE PSEUDO-VELOCITIES
        IF (LMPDATA_FV) THEN
           DO J=1,NLEV-1
              YZDOT(K+1,J,I) = ABS(YZDOT(K,J,I))*((ABS(YPSIG(J+1,I))-ABS(YPSIG(J,I)))         &
                             & /(ABS(YPSIG(J,I))+ABS(YPSIG(J+1,I))+REPS) )                    &
                             & -PDT*(YZDOT(K,J,I)/(PM0(J+1,I)+PM0(J,I)))*(                    &
                             & ((((YZDOT(K,J,I)+YZDOT(K,J+1,I))*ABS(YPSIG(J+1,I))             &
                             & -(YZDOT(K,J-1,I)+YZDOT(K,J,I))*ABS(YPSIG(J,I)))/PDETAH(J))     &
                             & /(ABS(YPSIG(J+1,I))+ABS(YPSIG(J,I)) + RESP))                   &
                             & +(HALF*(PM1(J+1,I)-PM0(J+1,I)+PM1(J,I)-PM0(J,I))/PDT) )
           END DO
        ELSE
           DO J=1,NLEV-1
              YZDOT(K+1,J,I) = (ABS(YZDOT(K,J,I))-(PDT/PDETAH(J))*(                           &
                             & ((YZDOT(K,J,I)**2)/(HALF*(PM0(J+1,I)+PM0(J,I))))))*(           &
                             & ((ABS(YPSIG(J+1,I))-ABS(YPSIG(J,I)))                           &
                             & /(ABS(YPSIG(J+1,I))+ABS(YPSIG(J,I))+REPS)) )                   &
                             & -HALF*PDT*(YZDOT(K,J,I)/(PM0(J+1,I)+PM0(J,I)))*(               &
                             & (((PDETAF(J)/PDETAF(J+1))*(YZDOT(K,J+1,I)-YZDOT(K,J,I))        &
                             & +(PDETAF(J+1)/PDETAF(J))*(YZDOT(K,J,I)-YZDOT(K,J-1,I)))        &
                             & /PDETAH(J))+((PM1(J+1,I)-PM0(J+1,I)+PM1(J,I)-PM0(J,I))/PDT))
           END DO
        END IF
        
        YZDOT(K+1,0,I)       = ZERO
        YZDOT(K+1,NLEV,I)    = ZERO
        
        !*************************************************************************************!
        !* GLOBAL FCT NON-OSCILLATORY TREATMENT
        IF (LFCT_MONO) THEN
           DO J=1,NLEV
              
              ZIN(J,I)       = ( MAX( YPSIG(NBC(J+1),I),YPSIG(J,I),                           &
                             & YPSIG(NBC(J-1),I),YPSI(0,NBC(J+1),I),                          &
                             & YPSI(0,J,I),YPSI(0,NBC(J-1),I))-YPSIG(J,I) )                   &
                             & /( (ONE/PM0(J,I))*(PDT/PDETAF(J))*(                            &
                             & MAX(ZERO,YZDOT(K+1,J-1,I))*YPSIG(NBC(J-1),I)                   &
                             & -MIN(ZERO,YZDOT(K+1,J,I))*YPSIG(NBC(J+1),I))+REPS  )
              
              ZOUT(J,I)      = -( MIN( YPSIG(NBC(J+1),I),YPSIG(J,I),                          &
                             & YPSIG(NBC(J-1),I),YPSI(0,NBC(J+1),I),                          &
                             & YPSI(0,J,I),YPSI(0,NBC(J-1),I))-YPSIG(J,I) )                   &
                             & /( (ONE/PM0(J,I))*(PDT/PDETAF(J))*(                            &
                             & MAX(ZERO,YZDOT(K+1,J,I))*YPSIG(J,I)                            &
                             & -MIN(ZERO,YZDOT(K+1,J-1,I))*YPSIG(J,I))+REPS  )
           END DO
           DO J=1,NLEV-1
              YZDOT(K+1,J,I) = MIN(ONE,ZOUT(J,I),ZIN(J+1,I)*MAX(ZERO,YZDOT(K+1,J,I))          &
                             & + MIN(ONE,ZIN(J,I),ZOUT(J+1,I)*MIN(ZERO,YZDOT(K+1,J,I))
           END DO   
        END IF
        !*************************************************************************************!
        !* CORRECTED UPWIND SCHEME
        DO J=1,NLEV   
           YPSI(K+1,J,I)     = YPSI(K,J,I) - (ONE/PM0(J,I))*(PDT/PDETAF(J))*(                 &
                             &  (MIN(ZERO,YZDOT(K+1,J,I))*YPSIG(NBC(J+1),I)                   &
                             &  +MAX(ZERO,YZDOT(K+1,J,I))*YPSIG(J,I))                         &
                             & -(MIN(ZERO,YZDOT(K+1,J-1,I))*YPSIG(J,I)                        &
                             &  +MAX(ZERO,YZDOT(K+1,J-1,I))*YPSIG(NBC(J-1),I)) )
        END DO
        DO J=1,NLEV   
           YPSIG(J,I)        = (PM0(J,I)/PM1(J,I))*YPSI(K+1,J,I)
        END DO
      END DO
    END DO

  ELSE ! INFINITY GAUGE
   
    DO I=1,NXPT
       DO K=1,NORD-1
        !*************************************************************************************!
        !* COMPUTE SURFACE PSEUDO-VELOCITIES
        IF (LMPDATA_FV) THEN
           DO J=1,NLEV-1
              YZDOT(K+1,J,I) = HALF*ABS(YZDOT(K,J,I))*((YPSIG(J+1,I)-YPSIG(J,I)))             &
                             & -HALF*PDT*(YZDOT(K,J,I)/(PM0(J+1,I)+PM0(J,I)))*(               &
                             & ((((YZDOT(K,J,I)+YZDOT(K,J+1,I))*YPSIG(J+1,I)                  &
                             & -(YZDOT(K,J-1,I)+YZDOT(K,J,I))*YPSIG(J,I))/PDETAH(J)))         &
                             & +HALF*(YPSIG(J+1,I)+YPSIG(J,I))*(                              &
                             & ((PM1(J+1,I)-PM0(J+1,I)+PM1(J,I)-PM0(J,I))/PDT)))
           END DO
        ELSE
           DO J=1,NLEV-1
              YZDOT(K+1,J,I) = HALF*(ABS(YZDOT(K,J,I))-(PDT/PDETAH(J))*(                      &
                             & ((YZDOT(K,J,I)**2)/(HALF*(PM0(J+1,I)+PM0(J,I))))))*(           &
                             & ((YPSIG(J+1,I)-YPSIG(J,I))) )                                  &
                             & -HALF*PDT*(YZDOT(K,J,I)/(PM0(J+1,I)+PM0(J,I)))*(               &
                             & HALF*(YPSIG(J+1,I)+YPSIG(J,I))*(                               &
                             & (((PDETAF(J)/PDETAF(J+1))*(YZDOT(K,J+1,I)-YZDOT(K,J,I))        &
                             & +(PDETAF(J+1)/PDETAF(J))*(YZDOT(K,J,I)-YZDOT(K,J-1,I)))        &
                             & /PDETAH(J))+((PM1(J+1,I)-PM0(J+1,I)+PM1(J,I)-PM0(J,I))/PDT)))
           END DO
        END IF

        YZDOT(K+1,0,I)       = ZERO
        YZDOT(K+1,NLEV,I)    = ZERO

        !*************************************************************************************!
        !* GLOBAL FCT NON-OSCILLATORY TREATMENT
        IF (LFCT_MONO) THEN
           DO J=1,NLEV
              
              ZIN(J,I)       = ( MAX( YPSIG(NBC(J+1),I),YPSIG(J,I),                           &
                             & YPSIG(NBC(J-1),I),YPSI(0,NBC(J+1),I),                          &
                             & YPSI(0,J,I),YPSI(0,NBC(J-1),I))-YPSIG(J,I) )                   &
                             & /( (ONE/PM0(J,I))*(PDT/PDETAF(J))*(                            &
                             & MAX(ZERO,YZDOT(K+1,J-1,I))                                     &
                             & -MIN(ZERO,YZDOT(K+1,J,I)))+REPS  )
              
              ZOUT(J,I)      = -( MIN( YPSIG(NBC(J+1),I),YPSIG(J,I),                          &
                             & YPSIG(NBC(J-1),I),YPSI(0,NBC(J+1),I),                          &
                             & YPSI(0,J,I),YPSI(0,NBC(J-1),I))-YPSIG(J,I) )                   &
                             & /( (ONE/PM0(J,I))*(PDT/PDETAF(J))*(                            &
                             & MAX(ZERO,YZDOT(K+1,J,I))                                       &
                             & -MIN(ZERO,YZDOT(K+1,J-1,I)))+REPS  )
           END DO
           DO J=1,NLEV-1
              YZDOT(K+1,J,I) = MIN(ONE,ZOUT(J,I),ZIN(J+1,I)*MAX(ZERO,YZDOT(K+1,J,I))          &
                             & + MIN(ONE,ZIN(J,I),ZOUT(J+1,I)*MIN(ZERO,YZDOT(K+1,J,I))
           END DO   
        END IF
        !*************************************************************************************!
        !* CORRECTED UPWIND SCHEME
           DO J=1,NLEV
              YPSI(K+1,J,I)  = YPSI(K,J,I) - (ONE/PM0(J,I))*(PDT/PDETAF(J))*(                 &
                             & YZDOT(K+1,J,I)-YZDOT(K+1,J-1,I) )
              YPSIG(J,I)     = (PM0(J,I)/PM1(J,I))*YPSI(K+1,J,I)
           END DO      
        END DO
     END DO
  END IF
  
  !*****************************************************
  DO I=1,NXPT
     DO J=1,NLEV
        PPSI_DEP(J,I)        = YPSI(NORD2,J,I)+ZPSI_MIN
     END DO   
 
     IF (LLMASS) THEN
        PZDOT(KS,0,I)        = ZERO
        PZDOT(KS+1,0,I)      = ZERO
        DO J=1,NLEV-1
           PZDOT(KS,J,I)     = HALF*YZDOT(NORD2,J,I)*(                &
                             & HALF*(PPSI_ARR(J+1,I)+PPSI_DEP(J+1,I)) &
                             & +HALF*(PPSI_ARR(J,I)+PPSI_DEP(J,I)) )
           PZDOT(KS+1,J,I)   = YZDOT(NORD2,J,I)
        END DO
        PZDOT(KS,NLEV,I)     = ZERO
        PZDOT(KS+1,NLEV,I)   = ZERO
     END IF
  END DO
  
END SUBROUTINE ADV_ZF
!################################################################################################!
!************************************************************************************************!
!************************************************************************************************!
!************************************************************************************************!
!************************************************************************************************!
!################################################################################################!
SUBROUTINE ADV_ZH(PPSI_DEP,PPSI_ARR,PZDOT,PM1,PM0,PDETAH,PDETAF,PDT,KS,LDMASS)
  
  IMPLICIT NONE

  REAL(8),                       INTENT(OUT)    :: PPSI_DEP(0:NLEV,NXPT)
  REAL(8),                       INTENT(IN)     :: PPSI_ARR(0:NLEV,NXPT)
  REAL(8),                       INTENT(INOUT)  :: PZDOT(3,NLEV,NXPT)
  REAL(8),                       INTENT(IN)     :: PM1(0:NLEV,NXPT)
  REAL(8),                       INTENT(IN)     :: PM0(0:NLEV,NXPT)
  REAL(8),                       INTENT(IN)     :: PDETAH(0:NLEV)
  REAL(8),                       INTENT(IN)     :: PDETAF(NLEV)
  REAL(8),                       INTENT(IN)     :: PDT
  INTEGER,                       INTENT(IN)     :: KS
  LOGICAL,OPTIONAL,              INTENT(IN)     :: LDMASS

  REAL(8)                                       :: YPSI(0:NORD2,0:NLEV,NXPT)
  REAL(8)                                       :: YZDOT(NORD2,0:NLEV+1,NXPT)
  REAL(8)                                       :: YPSIG(0:NLEV,NXPT)
  REAL(8)                                       :: ZIN_UP,ZOUT_UP,ZIN_DOWN,ZOUT_DOWN
  REAL(8)                                       :: ZPSI_MIN
  INTEGER                                       :: I,J,K
  LOGICAL                                       :: LLMASS

  IF (PRESENT(LDMASS)) THEN
     LLMASS = LDMASS
  ELSE
     LLMASS = .FALSE.
  END IF
  
  ZPSI_MIN = MINVAL(PPSI_ARR)
  
  YPSI(0,:,:)       = PPSI_ARR(:,:)-ZPSI_MIN
    
  YZDOT(1,1:NLEV,:) = PZDOT(KS,1:NLEV,:)
  YZDOT(:,0,:)      = ZERO
  YZDOT(:,NLEV+1,:) = ZERO

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !***************************************!
  ! DONOR-CELL UPWIND FIRST-ORDER SCHEME *!
  !***************************************!

  DO I=1,NXPT
     DO J=0,NLEV
        YPSI(1,J,I)          = YPSI(0,J,I) - (ONE/PM0(J,I))*(PDT/PDETAH(J))*(                 &
                             &  (MIN(ZERO,YZDOT(1,J+1,I))*YPSI(0,NBL(J+1),I)                  &
                             & +MAX(ZERO,YZDOT(1,J+1,I))*YPSI(0,J,I))                         &
                             & -(MIN(ZERO,YZDOT(1,J,I))*YPSI(0,J,I)                           &
                             & +MAX(ZERO,YZDOT(1,J,I))*YPSI(0,NBL(J-1),I)) )
        YPSIG(J,I)           = (PM0(J,I)/PM1(J,I))*YPSI(1,J,I)
     END DO
  END DO
  !******************!
  ! CORRECTION-STEP *!
  !******************!
  IF (.NOT.LGAUGE) THEN
    DO I=1,NXPT
       DO K=1,NORD-1
        !*************************************************************************************!
        !* COMPUTE SURFACE PSEUDO-VELOCITIES
          IF (LMPDATA_FV) THEN
             DO J=1,NLEV
                YZDOT(K+1,J,I) = ABS(YZDOT(K,J,I))*((ABS(YPSIG(J,I))-ABS(YPSIG(J-1,I)))       &
                               & /(ABS(YPSIG(J,I))+ABS(YPSIG(J-1,I))+REPS) )                  &
                               & -PDT*(YZDOT(K,J,I)/(PM0(J-1,I)+PM0(J,I)))*(                  &
                               & ((((YZDOT(K,J,I)+YZDOT(K,J+1,I))*ABS(YPSIG(J,I))             &
                               & -(YZDOT(K,J-1,I)+YZDOT(K,J,I))*ABS(YPSIG(J-1,I)))/PDETAF(J)) &
                               & /(ABS(YPSIG(J-1,I))+ABS(YPSIG(J,I)) + RESP))                 &
                               & +(HALF*(PM1(J,I)-PM0(J,I)+PM1(J-1,I)-PM0(J-1,I))/PDT) )
             END DO
          ELSE
             DO J=1,NLEV
                YZDOT(K+1,J,I) = (ABS(YZDOT(K,J,I))-(PDT/PDETAF(J))*(                         &
                               & ((YZDOT(K,J,I)**2)/(HALF*(PM0(J-1,I)+PM0(J,I))))))*(         &
                               & ((ABS(YPSIG(J,I))-ABS(YPSIG(J-1,I)))                         &
                               & /(ABS(YPSIG(J,I))+ABS(YPSIG(J-1,I))+REPS)) )                 &
                               & -HALF*PDT*(YZDOT(K,J,I)/(PM0(J-1,I)+PM0(J,I)))*(             &
                               & (((PDETAF(J)/PDETAH(J))*(YZDOT(K,J+1,I)-YZDOT(K,J,I))        &
                               & +(PDETAF(J)/PDETAH(J-1))*(YZDOT(K,J,I)-YZDOT(K,J-1,I)))      &
                               & /PDETAF(J))+((PM1(J,I)-PM0(J,I)+PM1(J-1,I)-PM0(J-1,I))/PDT))
             END DO
          END IF             

          YZDOT(K+1,0,I)       = ZERO
          YZDOT(K+1,NLEV+1,I)  = ZERO
        
          !***********************************************************************************!
          !* GLOBAL FCT NON-OSCILLATORY OPTION
          IF (LFCT_MONO) THEN
             DO J=0,NLEV
              
                ZIN(J,I)       = ( MAX( YPSIG(NBC(J+1),I),YPSIG(J,I),                         &
                               & YPSIG(NBL(J-1),I),YPSI(0,NBC(J+1),I),                        &
                               & YPSI(0,J,I),YPSI(0,NBL(J-1),I))-YPSIG(J,I) )                 &
                               & /( (ONE/PM0(J,I))*(PDT/PDETAH(J))*(                          &
                               & MAX(ZERO,YZDOT(K+1,J,I))*YPSIG(NBL(J-1),I)                   &
                               & -MIN(ZERO,YZDOT(K+1,J+1,I))*YPSIG(NBC(J+1),I))+REPS  )
              
                ZOUT(J,I)      = -( MIN( YPSIG(NBC(J+1),I),YPSIG(J,I),                        &
                               & YPSIG(NBL(J-1),I),YPSI(0,NBC(J+1),I),                        &
                               & YPSI(0,J,I),YPSI(0,NBL(J-1),I))-YPSIG(J,I) )                 &
                               & /( (ONE/PM0(J,I))*(PDT/PDETAH(J))*(                          &
                               & MAX(ZERO,YZDOT(K+1,J+1,I))*YPSIG(J,I)                        &
                               & -MIN(ZERO,YZDOT(K+1,J,I))*YPSIG(J,I))+REPS  )
             END DO
             DO J=1,NLEV
                YZDOT(K+1,J,I) = MIN(ONE,ZOUT(J-1,I),ZIN(J,I)*MAX(ZERO,YZDOT(K+1,J,I))        &
                               & + MIN(ONE,ZIN(J-1,I),ZOUT(J,I)*MIN(ZERO,YZDOT(K+1,J,I))
             END DO
          END IF          
          !***********************************************************************************!
          !* CORRECTED UPWIND SCHEME
          DO J=0,NLEV   
             YPSI(K+1,J,I)     = YPSI(K,J,I) - (ONE/PM0(J,I))*(PDT/PDETAH(J))*(               &
                               & (MIN(ZERO,YZDOT(K+1,J+1,I))*YPSIG(NBL(J+1),I)                &
                               & +MAX(ZERO,YZDOT(K+1,J+1,I))*YPSIG(J,I))                      &
                               & -(MIN(ZERO,YZDOT(K+1,J,I))*YPSIG(J,I)                        &
                               & +MAX(ZERO,YZDOT(K+1,J,I))*YPSIG(NBL(J-1),I)) )
          END DO
          DO J=0,NLEV   
             YPSIG(J,I)        = (PM0(J,I)/PM1(J,I))*YPSI(K+1,J,I)
          END DO
       END DO
    END DO

  ELSE
  !*************************!
  ! SECOND CORRECTION-STEP *!
  !*************************!

    DO I=1,NXPT
       DO K=1,NORD-1

        !*************************************************************************************!
          !* COMPUTE SURFACE PSEUDO-VELOCITIES
          IF (LMPDATA_FV) THEN
             DO J=1,NLEV
                YZDOT(K+1,J,I) = HALF*ABS(YZDOT(K,J,I))*((YPSIG(J,I)-YPSIG(J-1,I)))           &
                             & -HALF*PDT*(YZDOT(K,J,I)/(PM0(J-1,I)+PM0(J,I)))*(               &
                             & ((((YZDOT(K,J,I)+YZDOT(K,J+1,I))*YPSIG(J,I)                    &
                             & -(YZDOT(K,J-1,I)+YZDOT(K,J,I))*YPSIG(J-1,I))/PDETAF(J)))       &
                             & +HALF*(YPSIG(J-1,I)+YPSIG(J,I))*(                              &
                             & ((PM1(J-1,I)-PM0(J-1,I)+PM1(J,I)-PM0(J,I))/PDT)))
             END DO
          ELSE
             DO J=1,NLEV
                YZDOT(K+1,J,I) = HALF*(ABS(YZDOT(K,J,I))-(PDT/PDETAF(J))*(                    &
                             & ((YZDOT(K,J,I)**2)/(HALF*(PM0(J-1,I)+PM0(J,I))))))*(           &
                             & ((YPSIG(J,I)-YPSIG(J-1,I))) )                                  &
                             & -HALF*PDT*(YZDOT(K,J,I)/(PM0(J-1,I)+PM0(J,I)))*(               &
                             & HALF*(YPSIG(J-1,I)+YPSIG(J,I))*(                               &
                             & (((PDETAF(J)/PDETAH(J))*(YZDOT(K,J+1,I)-YZDOT(K,J,I))          &
                             & +(PDETAF(J)/PDETAH(J-1))*(YZDOT(K,J,I)-YZDOT(K,J-1,I)))        &
                             & /PDETAF(J))+((PM1(J-1,I)-PM0(J-1,I)+PM1(J,I)-PM0(J,I))/PDT)))
             END DO
          END IF
       
          YZDOT(K+1,0,I)      = ZERO
          YZDOT(K+1,NLEV+1,I) = ZERO
        
          !***********************************************************************************!
          !* GLOBAL FCT NON-OSCILLATORY OPTION
          IF (LFCT_MONO) THEN
             DO J=0,NLEV
              
                ZIN(J,I)     = ( MAX( YPSIG(NBC(J+1),I),YPSIG(J,I),                           &
                             & YPSIG(NBL(J-1),I),YPSI(0,NBC(J+1),I),                          &
                             & YPSI(0,J,I),YPSI(0,NBL(J-1),I))-YPSIG(J,I) )                   &
                             & /( (ONE/PM0(J,I))*(PDT/PDETAH(J))*(                            &
                             & MAX(ZERO,YZDOT(K+1,J,I))                                       &
                             & -MIN(ZERO,YZDOT(K+1,J+1,I)))+REPS  )
              
                ZOUT(J,I)    = -( MIN( YPSIG(NBC(J+1),I),YPSIG(J,I),                          &
                             & YPSIG(NBL(J-1),I),YPSI(0,NBC(J+1),I),                          &
                             & YPSI(0,J,I),YPSI(0,NBL(J-1),I))-YPSIG(J,I) )                   &
                             & /( (ONE/PM0(J,I))*(PDT/PDETAH(J))*(                            &
                             & MAX(ZERO,YZDOT(K+1,J+1,I))                                     &
                             & -MIN(ZERO,YZDOT(K+1,J,I)))+REPS  )
             END DO
             DO J=1,NLEV
                YZDOT(K+1,J,I) = MIN(ONE,ZOUT(J-1,I),ZIN(J,I)*MAX(ZERO,YZDOT(K+1,J,I))        &
                             & + MIN(ONE,ZIN(J-1,I),ZOUT(J,I)*MIN(ZERO,YZDOT(K+1,J,I))
             END DO   
          END IF
          !***********************************************************************************!
          !* CORRECTED UPWIND SCHEME
          DO J=0,NLEV
             YPSI(K+1,J,I)     = YPSI(K,J,I) - (ONE/ZM0(J,I))*(PDT/PDETAH(J))*(               &
                               &  YZDOT(K+1,J+1,I)-YZDOT(K+1,J,I) )
             YPSIG(J,I)        = (PM0(J,I)/PM1(J,I))*YPSI(K+1,J,I)
          END DO
       END DO
    END DO
  END IF
 
  !*******************************************************************************************!
  DO I=1,NXPT
     DO J=0,NLEV
        PPSI_DEP(J,I)      = YPSI(NORD2,J+1,I)+ZPSI_MIN
     END DO
     IF (LLMASS) THEN
        DO J=1,NLEV
           PZDOT(KS,J,I)   = HALF*YZDOT(NORD2,J,I)*(                                          &
                           & HALF*(PPSI_ARR(NBL(J+1),I)+PPSI_DEP(J,I))                        &
                           & + HALF*(PPSI_ARR(J,I)+PPSI_DEP(NBL(J+1),I)) )
           PZDOT(KS+1,J,I) = YZDOT(NORD2,J,I)
        END DO
     END IF
  END DO
  
END SUBROUTINE ADV_ZH
!################################################################################################!
!************************************************************************************************!
!************************************************************************************************!
!************************************************************************************************!
!************************************************************************************************!
!************************************************************************************************!
!************************************************************************************************!
!************************************************************************************************!
!************************************************************************************************!
!************************************************************************************************!
!************************************************************************************************!
!************************************************************************************************!
!************************************************************************************************!
!************************************************************************************************!
!************************************************************************************************!
!************************************************************************************************!
!************************************************************************************************!
!************************************************************************************************!
!************************************************************************************************!
!************************************************************************************************!
!************************************************************************************************!
!************************************************************************************************!
!************************************************************************************************!
!************************************************************************************************!
!************************************************************************************************!
!################################################################################################!
SUBROUTINE MPDATA_UNSPLIT_2D_FULL_MASS(PMF,PXDOTF,PZDOTH,PDETAF)

  IMPLICIT NONE

  REAL(8),         INTENT(INOUT)  :: PMF(0:1,NLEV,NXPT)
  REAL(8),         INTENT(INOUT)  :: PXDOTF(NLEV,NXPT)
  REAL(8),         INTENT(INOUT)  :: PZDOTH(0:NLEV,NXPT)
  REAL(8),         INTENT(IN)     :: PDETAF(NLEV)

  REAL(8)                         :: ZUN(NLEV,NXPT)

  ZUN(:,:) = ONE
  
  CALL ADV_XZ(PMF(1,:,:),PMF(0,:,:),PXDOTF,     &
       & PZDOTH,ZUN,ZUN,PDETAF,RDT,LDMASS=.TRUE.)
       
END SUBROUTINE MPDATA_UNSPLIT_2D_FULL_MASS
!################################################################################################!
!************************************************************************************************!
!************************************************************************************************!
!################################################################################################!
SUBROUTINE MPDATA_UNSPLIT_2D_FULL_FIELD(PPSI_DEP,PPSI_ARR,PXDOTF,PZDOTH,PMF,PDETAF)

  IMPLICIT NONE

  REAL(8),    INTENT(OUT)    :: PPSI_DEP(NLEV,NXPT)
  REAL(8),    INTENT(IN)     :: PPSI_ARR(NLEV,NXPT)
  REAL(8),    INTENT(INOUT)  :: PXDOTF(NLEV,NXPT)
  REAL(8),    INTENT(INOUT)  :: PZDOTH(0:NLEV,NXPT)
  REAL(8),    INTENT(IN)     :: PMF(0:1,NLEV,NXPT)
  REAL(8),    INTENT(IN)     :: PDETAF(NLEV)

 
  CALL ADV_XZ(PPSI_DEP,PPSI_ARR,PXDOTF,         &
       & PZDOTH,PMF(1,:,:),PMF(0,:,:),PDETAF,RDT)
  
END SUBROUTINE  MPDATA_UNSPLIT_2D_FULL_FIELD
!######################################################################################################!
!******************************************************************************************************!
!******************************************************************************************************!
!******************************************************************************************************!
!******************************************************************************************************!
!######################################################################################################!
SUBROUTINE ADV_XZ(PPSI_DEP,PPSI_ARR,PXDOT,PZDOT,PM1,PM0,PDETAF,PDETAH,PDT,LDMASS)

  IMPLICIT NONE

  REAL(8),                       INTENT(OUT)    :: PPSI_DEP(NLEV,NXPT)
  REAL(8),                       INTENT(IN)     :: PPSI_ARR(NLEV,NXPT)
  REAL(8),                       INTENT(INOUT)  :: PXDOT(NLEV,NXPT)
  REAL(8),                       INTENT(INOUT)  :: PZDOT(0:NLEV,NXPT)
  REAL(8),                       INTENT(IN)     :: PM1(NLEV,NXPT)
  REAL(8),                       INTENT(IN)     :: PM0(NLEV,NXPT)
  REAL(8),                       INTENT(IN)     :: PDETAF(NLEV)
  REAL(8),                       INTENT(IN)     :: PDETAH(0:NLEV)
  REAL(8),                       INTENT(IN)     :: PDT
  LOGICAL,OPTIONAL,              INTENT(IN)     :: LDMASS
  
  REAL(8)                                       :: YPSI(0:NORD,NLEV,NXPT)
  REAL(8)                                       :: YXDOT(NORD,NLEV,NXPT)
  REAL(8)                                       :: YZDOT(NORD,0:NLEV,NXPT)
  REAL(8)                                       :: YPSIG(NLEV,NXPT)
  REAL(8)                                       :: ZIN(NLEV,NXPT),ZOUT(NLEV,NXPT)
  REAL(8)                                       :: ZMAX,ZMIN
  REAL(8)                                       :: ZPSI_MIN
  REAL(8)                                       :: ZXDOT_MONO,ZZDOT_MONO
  INTEGER(8)                                    :: I,J,K
  LOGICAL                                       :: LLMASS

  IF (PRESENT(LDMASS)) THEN
     LLMASS = LDMASS
  ELSE
     LLMASS = .FALSE.
  END IF   

  ZPSI_MIN = MINVAL(PPSI_ARR)

  DO I=1,NXPT
     YPSI(0,:,I)  = PPSI_ARR(:,I) - ZPSI_MIN
     YXDOT(1,:,I) = PXDOT(:,I)
     YZDOT(1,:,I) = PZDOT(:,I)
  END DO

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !***************************************!
  ! DONOR-CELL UPWIND FIRST-ORDER SCHEME *!
  !***************************************!

  DO I=1,NXPT
     DO J=1,NLEV
        YPSI(1,J,I)       = YPSI(0,J,I)                                                              &
                          & -(ONE/PM0(J,I))*(PDT/RDX)*(                                              &
                          & (MIN(ZERO,YXDOT(1,J,I))*YPSI(0,J,NPERIO(I+1))                            &
                          & +MAX(ZERO,YXDOT(1,J,I))*YPSI(0,J,I))                                     &
                          & -(MIN(ZERO,YXDOT(1,J,NPERIO(I-1)))*YPSI(0,J,I)                           & 
                          & +MAX(ZERO,YXDOT(1,J,NPERIO(I-1)))*YPSI(0,J,NPERIO(I-1))) )               &
                          & -(ONE/PM0(J,I))*(PDT/PDETAF(J))*(                                        &
                          & (MIN(ZERO,YZDOT(1,J,I))*YPSI(0,NBC(J+1),I)                               &
                          & +MAX(ZERO,YZDOT(1,J,I))*YPSI(0,J,I))                                     &
                          & -(MIN(ZERO,YZDOT(1,J-1,I))*YPSI(0,J,I)                                   &
                          & +MAX(ZERO,YZDOT(1,J-1,I))*YPSI(0,NBC(J-1),I)) )
        YPSIG(J,I)        = (PM0(J,I)/PM1(J,I))*YPSI(1,J,I)
     END DO   
  END DO
  
  !************************!
  ! FIRST CORRECTION-STEP *!
  !************************!
  IF (.NOT.LGAUGE) THEN
     DO K=1,NORD-1   
        IF (LMPDATA_FV) THEN
           !... [FV] FORMULATION with DIV_C OPTION (Cf. XY tracer_labo) :
           !    Fluxes are first computed, then interpolated at the (i+0.5,j+0.5), (i-0.5,j+0.5)
           !    and (i+0.5,j+0.5), (i+0.5,j-0.5) 
           DO I=1,NXPT   
              DO J=1,NLEV
                 YXDOT(K+1,J,I) = ABS(YXDOT(K,J,I))*((ABS(YPSIG(J,NPERIO(I+1)))-ABS(YPSIG(J,I)))     &
                          & /(ABS(YPSIG(J,I))+ABS(YPSIG(J,NPERIO(I+1)))+REPS) )                      &
                          & -PDT*(YXDOT(K,J,I)/(PM0(J,NPERIO(I+1))+PM0(J,I)))*(                      &
                          & (((((YXDOT(K,J,I)+YXDOT(K,J,NPERIO(I+1)))*ABS(YPSIG(J,NPERIO(I+1)))      &
                          & -(YXDOT(K,J,NPERIO(I-1))+YXDOT(K,J,I))*ABS(YPSIG(J,I)))/RDX)             &
                          & +(HALF*((PDETAF(J)/PDETAH(J))*((YZDOT(K,J,I)*ABS(YPSIG(NBC(J+1),I)))     &
                          & +(YZDOT(K,J,NPERIO(I+1))*ABS(YPSIG(NBC(J+1),NPERIO(I+1)))))              &  
                          & +(PDETAF(NBC(J+1))/PDETAH(J))*((YZDOT(K,J,I)*ABS(YPSIG(J,I)))            &
                          & +(YZDOT(K,J,NPERIO(I+1))*ABS(YPSIG(J,NPERIO(I+1)))))                     & 
                          & -(PDETAF(NBC(J-1))/PDETAH(J-1))*((YZDOT(K,J-1,I)*ABS(YPSIG(J,I)))        &
                          & +(YZDOT(K,J-1,NPERIO(I+1))*ABS(YPSIG(J,NPERIO(I+1)))))                   &
                          & -(PDETAF(J)/PDETAH(J-1))*((YZDOT(K,J-1,I)*ABS(YPSIG(NBC(J-1),I)))        &
                          & +(YZDOT(K,J-1,NPERIO(I+1))*ABS(YPSIG(NBC(J-1),NPERIO(I+1)))))) ))        &
                          & /((0.2d0*((PDETAF(J)/PDETAH(J))*(ABS(YPSIG(NBC(J+1),I))                  &
                          & +ABS(YPSIG(NBC(J+1),NPERIO(I+1))))+(ONE+(PDETAF(NBC(J+1))/PDETAH(J))     &
                          & +(PDETAF(NBC(J-1))/PDETAH(J-1)))*(ABS(YPSIG(J,I))                        &
                          & +ABS(YPSIG(J,NPERIO(I+1))))+(PDETAF(J)/PDETAH(J-1))*(                    &
                          &  ABS(YPSIG(NBC(J-1),I))+ABS(YPSIG(NBC(J-1),NPERIO(I+1))))))+REPS))       &
                          & +(HALF*(PM1(J,NPERIO(I+1))-PM0(J,NPERIO(I+1))+PM1(J,I)-PM0(J,I))/PDT) ) 

              END DO    
              DO J=1,NLEV-1
                 YZDOT(K+1,J,I) = ABS(YZDOT(K,J,I))*((ABS(YPSIG(J+1,I))-ABS(YPSIG(J,I)))             & 
                          & /(ABS(YPSIG(J,I))+ABS(YPSIG(J+1,I))+REPS) )                              &
                          & -PDT*(YZDOT(K,J,I)/(PM0(J+1,I)+PM0(J,I)))*(                              &
                          & (((((YZDOT(K,J+1,I)+YZDOT(K,J,I))*ABS(YPSIG(J+1,I))                      &  
                          & -(YZDOT(K,J-1,I)+YZDOT(K,J,I))*ABS(YPSIG(J,I))))/PDETAH(J))              &
                          & +(HALF*(((PDETAF(J)/PDETAH(J))*(YXDOT(K,J+1,I)*(                         &
                          &  ABS(YPSIG(J+1,NPERIO(I+1)))+ABS(YPSIG(J+1,I)))                          & 
                          & -YXDOT(K,J+1,NPERIO(I-1))*(ABS(YPSIG(J+1,NPERIO(I-1)))                   &   
                          & +ABS(YPSIG(J+1,I))))+(PDETAF(J+1)/PDETAH(J))*(YXDOT(K,J,I)*(             & 
                          &  ABS(YPSIG(J,NPERIO(I+1)))+ABS(YPSIG(J,I)))-YXDOT(K,J,NPERIO(I-1))*(     &
                          &  ABS(YPSIG(J,NPERIO(I-1)))+ABS(YPSIG(J,I)))))/RDX))                      &
                          & /((0.2d0*((PDETAF(J)/PDETAH(J))*(ABS(YPSIG(J+1,NPERIO(I+1)))             &
                          & +ABS(YPSIG(J+1,NPERIO(I-1))))+(PDETAF(J+1)/PDETAH(J))*(                  &
                          &  ABS(YPSIG(J,NPERIO(I+1)))+ABS(YPSIG(J,NPERIO(I-1))))                    &
                          & +(ONE+TWO*(PDETAF(J)/PDETAH(J)))*ABS(YPSIG(J+1,I))                       &
                          & +(ONE+TWO*(PDETAF(J+1)/PDETAH(J)))*ABS(YPSIG(J,I))))+REPS))              &
                          & +(HALF*(PM1(J+1,I)-PM0(J+1,I)+PM1(J,I)-PM0(J,I))/RDT) )
              END DO
              YZDOT(K+1,0,I)    = ZERO
              YZDOT(K+1,NLEV,I) = ZERO
           END DO
        ELSE                 !... [FD] FORMULATION
           DO I=1,NXPT   
              DO J=1,NLEV
                 YXDOT(K+1,J,I) = (ABS(YXDOT(K,J,I))-(PDT/RDX)*((YXDOT(K,J,I)**2)                    &
                          & /(HALF*(PM0(J,NPERIO(I+1))+PM0(J,I)))))*(                                &
                          & (ABS(YPSIG(J,NPERIO(I+1)))-ABS(YPSIG(J,I)))                              &
                          & /(ABS(YPSIG(J,I))+ABS(YPSIG(J,NPERIO(I+1)))+REPS) )                      &
                          & -0.25d0*(PDT/PDETAF(J))*(YXDOT(K,J,I)/(PM0(J,NPERIO(I+1))+PM0(J,I)))*(   & 
                          & (YZDOT(K,J,NPERIO(I+1))+YZDOT(K,J,I)+YZDOT(K,J-1,NPERIO(I+1))            &
                          & +YZDOT(K,J-1,I))*( (((PDETAF(J)/PDETAH(J))*(                             &
                          & ABS(YPSIG(NBC(J+1),NPERIO(I+1)))+ABS(YPSIG(NBC(J+1),I)))                 & 
                          & +((PDETAF(J+1)/PDETAH(J))-(PDETAF(J-1)/PDETAH(J-1)))*(                   &
                          & ABS(YPSIG(J,NPERIO(I+1)))+ABS(YPSIG(J,I)))                               &
                          & -(PDETAF(J)/PDETAH(J-1))*(ABS(YPSIG(NBC(J-1),NPERIO(I+1)))               &
                          & +ABS(YPSIG(NBC(J-1),I))))/((PDETAF(J)/PDETAH(J))*(                       &
                          & ABS(YPSIG(NBC(J+1),NPERIO(I+1)))+ABS(YPSIG(NBC(J+1),I)))                 & 
                          & +((PDETAF(J+1)/PDETAH(J))+(PDETAF(J-1)/PDETAH(J-1)))*(                   &
                          & ABS(YPSIG(J,NPERIO(I+1)))+ABS(YPSIG(J,I)))                               &
                          & +(PDETAF(J)/PDETAH(J-1))*(ABS(YPSIG(NBC(J-1),NPERIO(I+1)))               &
                          & +ABS(YPSIG(NBC(J-1),I))) + REPS) ))                                      &
                          & -HALF*PDT*(YXDOT(K,J,I)/(PM0(J,NPERIO(I+1))+PM0(J,I)))*(                 &
                          & ((YXDOT(K,J,NPERIO(I+1))-YXDOT(K,J,NPERIO(I-1)))/RDX)                    &
                          & +((YZDOT(K,J,NPERIO(I+1))+YZDOT(K,J,I)                                   &
                          & -YZDOT(K,J-1,NPERIO(I+1))-YZDOT(K,J-1,I))/PDETAF(J))                     &
                          & +((PM1(J,NPERIO(I+1))-PM0(J,NPERIO(I+1))+PM1(J,I)-PM0(J,I) )/PDT) )

              END DO    
              DO J=1,NLEV-1
                 YZDOT(K+1,J,I) = (ABS(YZDOT(K,J,I))-(PDT/PDETAH(J))*((YZDOT(K,J,I)**2)              &
                          & /(HALF*(PM0(J+1,I)+PM0(J,I)))))*(                                        &
                          & (ABS(YPSIG(J+1,I))-ABS(YPSIG(J,I)))                                      &
                          & /(ABS(YPSIG(J,I))+ABS(YPSIG(J+1,I))+REPS) )                              &
                          & -0.25d0*(PDT/RDX)*(YZDOT(K,J,I)/(PM0(J+1,I)+PM0(J,I)))*(                 &
                          & ((PDETAF(J)/PDETAH(J))*(YXDOT(K,J+1,I)+YXDOT(K,J+1,NPERIO(I-1)))         &
                          & +(PDETAF(J+1)/PDETAH(J))*(YXDOT(K,J,I)+YXDOT(K,J,NPERIO(I-1))))*(        &
                          & ((PDETAF(J)/PDETAH(J))*(ABS(YPSIG(J+1,NPERIO(I+1)))                      &
                          & -ABS(YPSIG(J+1,NPERIO(I-1))))+(PDETAF(J+1)/PDETAH(J))*(                  &
                          & ABS(YPSIG(J,NPERIO(I+1)))-ABS(YPSIG(J,NPERIO(I-1)))))                    &
                          & /((PDETAF(J)/PDETAH(J))*(ABS(YPSIG(J+1,NPERIO(I+1)))                     &
                          & +ABS(YPSIG(J+1,NPERIO(I-1))))+(PDETAF(J+1)/PDETAH(J))*(                  &
                          & ABS(YPSIG(J,NPERIO(I+1)))+ABS(YPSIG(J,NPERIO(I-1))))+REPS)))             &
                          & -HALF*PDT*(YZDOT(K,J,I)/(PM0(J+1,I)+PM0(J,I)))*(                         &
                          & ((PDETAF(J)/PDETAH(J))*(YXDOT(K,J+1,I)-YXDOT(K,J+1,NPERIO(I-1)))         &
                          & +(PDETAF(J+1)/PDETAH(J))*(YXDOT(K,J,I)-YXDOT(K,J,NPERIO(I-1)))/RDX)      &
                          & +(((PDETAF(J)/PDETAF(J+1))*(YZDOT(K,J+1,I)-YZDOT(K,J,I))                 &
                          & +(PDETAF(J+1)/PDETAF(J))*(YZDOT(K,J,I)-YZDOT(K,J-1,I)))/PDETAH(J))       &
                          & +((PM1(J+1,I)-PM0(J+1,I)+PM1(J,I)-PM0(J,I))/RDT) )
              END DO
              YZDOT(K+1,0,I)    = ZERO
              YZDOT(K+1,NLEV,I) = ZERO
           END DO
        END IF
        !********************************************************************************************!
        !* GLOBAL FCT NON-OSCILLATORY TREATMENT
        IF (LFCT_MONO) THEN
           DO I=1,NXPT
              DO J=1,NLEV

                 ZMAX        = MAX( YPSIG(J,NPERIO(I+1)),YPSIG(NBC(J+1),I),YPSIG(J,I),               &
                             & YPSIG(J,NPERIO(I-1)),YPSIG(NBC(J-1),I),YPSI(0,NBC(J+1),I),            &
                             & YPSI(0,J,NPERIO(I+1)),YPSI(0,J,I),YPSI(0,NBC(J-1),I),                 &
                             & YPSI(0,J,NPERIO(I-1)))

                 ZIN(J,I)    = ( ZMAX-YPSIG(J,I))/( (ONE/PM0(J,I))*(                                 & 
                             & (PDT/RDX)*(MAX(ZERO,YXDOT(K+1,J,NPERIO(I-1)))*YPSIG(J,NPERIO(I-1))    & 
                             & -MIN(ZERO,YXDOT(K+1,J,I))*YPSIG(J,NPERIO(I+1)))                       &
                             & +(PDT/PDETAF(J))*(MAX(ZERO,YZDOT(K+1,J-1,I))*YPSIG(NBC(J-1),I)        &
                             & -MIN(ZERO,YZDOT(K+1,J,I))*YPSIG(NBC(J+1),I)) ) + REPS  )


                 ZMIN        =  MIN( YPSIG(J,NPERIO(I+1)),YPSIG(NBC(J+1),I),YPSIG(J,I),              &
                             & YPSIG(J,NPERIO(I-1)),YPSIG(NBC(J-1),I),YPSI(0,NBC(J+1),I),            &
                             & YPSI(0,J,NPERIO(I+1)),YPSI(0,J,I),YPSI(0,NBC(J-1),I),                 &
                             & YPSI(0,J,NPERIO(I-1)))
                 
                 ZOUT(J,I)   = (YPSIG(J,I)-ZMIN)/( (ONE/PM0(J,I))*(                                  &
                             & (PDT/RDX)*(MAX(ZERO,YXDOT(K+1,J,I))*YPSIG(J,I)                        &
                             & -MIN(ZERO,YXDOT(K+1,J,I))*YPSIG(J,I))                                 &
                             & +(PDT/PDETAF(J))*(MAX(ZERO,YZDOT(K+1,J,I))*YPSIG(J,I)                 &
                             & -MIN(ZERO,YZDOT(K+1,J-1,I))*YPSIG(J,I)) ) + REPS  )
                   
              END DO
           END DO
           DO I=1,NXPT
              DO J=1,NLEV
                 YXDOT(K+1,J,I) = MIN(ONE,ZOUT(J,I),ZIN(J,NPERIO(I+1)))*MAX(ZERO,YXDOT(K+1,J,I))     & 
                                & +MIN(ONE,ZIN(J,I),ZOUT(J,NPERIO(I+1)))*MIN(ZERO,YXDOT(K+1,J,I))
              END DO
              DO J=1,NLEV-1
                 YZDOT(K+1,J,I) =MIN(ONE,ZOUT(J,I),ZIN(J+1,I))*MAX(ZERO,YZDOT(K+1,J,I))              & 
                                & + MIN(ONE,ZIN(J,I),ZOUT(J+1,I)*MIN(ZERO,YZDOT(K+1,J,I))
              END DO
           END DO      
        END IF     
        !********************************************************************************************!
        !* CORRECTED UPWIND SCHEME
        DO I=1,NXPT
           DO J=1,NLEV
              YPSI(K+1,J,I)  = YPSI(K,J,I)-(ONE/PM0(J,I))*(PDT/RDX)*(                                &
                             &  (MIN(ZERO,YXDOT(K+1,J,I))*YPSIG(J,NPERIO(I+1))                       &
                             & +MAX(ZERO,YXDOT(K+1,J,I))*YPSIG(J,I))                                 &
                             & -(MIN(ZERO,YXDOT(K+1,J,NPERIO(I-1)))*YPSIG(J,I)                       &
                             & +MAX(ZERO,YXDOT(K+1,J,NPERIO(I-1)))*YPSIG(J,NPERIO(I-1))) )           &
                             & -(ONE/PM0(J,I))*(PDT/PDETAF(J))*(                                     &
                             &  (MIN(ZERO,YZDOT(K+1,J,I))*YPSIG(NBC(J+1),I)                          &
                             & +MAX(ZERO,YZDOT(K+1,J,I))*YPSIG(J,I))                                 &
                             & -(MIN(ZERO,YZDOT(K+1,J-1,I))*YPSIG(J,I)                               &
                             & +MAX(ZERO,YZDOT(K+1,J-1,I))*YPSIG(NBC(J-1),I)) )
           END DO   
        END DO
        DO I=1,NXPT
           DO J=1,NLEV
              YPSIG(J,I)     = (PM0(J,I)/PM1(J,I))*YPSI(K+1,J,I)
           END DO   
        END DO
     END DO

  ELSE ! INFINITY GAUGE OPTION
  
     DO K=1,NORD-1
     !***********************************************************************************************!
     !* COMPUTE SURFACE PSEUDO-VELOCITIES
         IF (LMPDATA_FV) THEN
           DO I=1,NXPT
              DO J=1,NLEV
                  YXDOT(K+1,J,I) = HALF*ABS(YXDOT(K,J,I))*((YPSIG(J,NPERIO(I+1))-YPSIG(J,I)))        &
                          & -HALF*PDT*(YXDOT(K,J,I)/(PM0(J,NPERIO(I+1))+PM0(J,I)))*(                 &
                          & (((((YXDOT(K,J,I)+YXDOT(K,J,NPERIO(I+1)))*YPSIG(J,NPERIO(I+1))           &
                          & -(YXDOT(K,J,NPERIO(I-1))+YXDOT(K,J,I))*YPSIG(J,I))/RDX)                  &
                          & +(HALF*((PDETAF(J)/PDETAH(J))*((YZDOT(K,J,I)*YPSIG(NBC(J+1),I))          &
                          & +(YZDOT(K,J,NPERIO(I+1))*YPSIG(NBC(J+1),NPERIO(I+1))))                   &  
                          & +(PDETAF(NBC(J+1))/PDETAH(J))*((YZDOT(K,J,I)*YPSIG(J,I))                 &
                          & +(YZDOT(K,J,NPERIO(I+1))*YPSIG(J,NPERIO(I+1))))                          & 
                          & -(PDETAF(NBC(J-1))/PDETAH(J-1))*((YZDOT(K,J-1,I)*YPSIG(J,I))             &
                          & +(YZDOT(K,J-1,NPERIO(I+1))*YPSIG(J,NPERIO(I+1))))                        &
                          & -(PDETAF(J)/PDETAH(J-1))*((YZDOT(K,J-1,I)*YPSIG(NBC(J-1),I))             &
                          & +(YZDOT(K,J-1,NPERIO(I+1))*YPSIG(NBC(J-1),NPERIO(I+1))))) )))            &
                          & +(HALF*(YPSIG(J,I)+YPSIG(J,NPERIO(I+1)))*(                               &
                          &  PM1(J,NPERIO(I+1))-PM0(J,NPERIO(I+1))+PM1(J,I)-PM0(J,I))/PDT) ) 
              END DO
              DO J=1,NLEV-1
                 YZDOT(K+1,J,I) = HALF*ABS(YZDOT(K,J,I))*(YPSIG(J+1,I)-YPSIG(J,I))                   &
                          & -HALF*PDT*(YZDOT(K,J,I)/(PM0(J+1,I)+PM0(J,I)))*(                         &
                          & (((((YZDOT(K,J+1,I)+YZDOT(K,J,I))*YPSIG(J+1,I)                           &  
                          & -(YZDOT(K,J-1,I)+YZDOT(K,J,I))*YPSIG(J,I)))/PDETAH(J))                   &
                          & +(HALF*(((PDETAF(J)/PDETAH(J))*(YXDOT(K,J+1,I)*(                         &
                          &  YPSIG(J+1,NPERIO(I+1))+YPSIG(J+1,I))                                    & 
                          & -YXDOT(K,J+1,NPERIO(I-1))*(YPSIG(J+1,NPERIO(I-1))                        &   
                          & +YPSIG(J+1,I)))+(PDETAF(J+1)/PDETAH(J))*(YXDOT(K,J,I)*(                  & 
                          &  YPSIG(J,NPERIO(I+1))+YPSIG(J,I))-YXDOT(K,J,NPERIO(I-1))*(               &
                          &  YPSIG(J,NPERIO(I-1))+YPSIG(J,I))))/RDX)))                               &
                          & +(HALF*(YPSIG(J+1,I)+YPSIG(J,I))*(                                       &
                          &  PM1(J+1,I)-PM0(J+1,I)+PM1(J,I)-PM0(J,I))/PDT) )
              END DO
              YZDOT(K+1,0,I)    = ZERO
              YZDOT(K+1,NLEV,I) = ZERO
           END DO
        ELSE
           DO I=1,NXPT   
              DO J=1,NLEV
                 YXDOT(K+1,J,I) = HALF*(ABS(YXDOT(K,J,I))-(PDT/RDX)*((YXDOT(K,J,I)**2)               &
                          & /(HALF*(PM0(J,NPERIO(I+1))+PM0(J,I)))))*(                                &
                          &  YPSIG(J,NPERIO(I+1))-YPSIG(J,I))-(PDT/PDETAF(J))*(                      &
                          &  YXDOT(K,J,I)/(PM0(J,NPERIO(I+1))+PM0(J,I)))*(0.25d0*(                   & 
                          &  (YZDOT(K,J,NPERIO(I+1))+YZDOT(K,J,I)+YZDOT(K,J-1,NPERIO(I+1))           &
                          & +YZDOT(K,J-1,I)))*(0.25d0*((((PDETAF(J)/PDETAH(J))*(                     &
                          &  YPSIG(NBC(J+1),NPERIO(I+1))+YPSIG(NBC(J+1),I))                          & 
                          & +((PDETAF(J+1)/PDETAH(J))-(PDETAF(J-1)/PDETAH(J-1)))*(                   &
                          &  YPSIG(J,NPERIO(I+1))+YPSIG(J,I))-(PDETAF(J)/PDETAH(J-1))*(              &
                          &  YPSIG(NBC(J-1),NPERIO(I+1))+YPSIG(NBC(J-1),I)))))))                     &
                          & -HALF*PDT*(YXDOT(K,J,I)/(PM0(J,NPERIO(I+1))+PM0(J,I)))*(                 &
                          &  HALF*(YPSIG(J,I)+YPSIG(J,NPERIO(I+1)))*(                                &
                          &  ((YXDOT(K,J,NPERIO(I+1))-YXDOT(K,J,NPERIO(I-1)))/RDX)                   &
                          & +((YZDOT(K,J,NPERIO(I+1))+YZDOT(K,J,I)                                   &
                          & -YZDOT(K,J-1,NPERIO(I+1))-YZDOT(K,J-1,I))/PDETAF(J))                     &
                          & +((PM1(J,NPERIO(I+1))-PM0(J,NPERIO(I+1))+PM1(J,I)-PM0(J,I) )/PDT)) )
              END DO    
              DO J=1,NLEV-1
                 YZDOT(K+1,J,I) = HALF*(ABS(YZDOT(K,J,I))-(PDT/PDETAH(J))*((YZDOT(K,J,I)**2)         &
                          & /(HALF*(PM0(J+1,I)+PM0(J,I)))))*(YPSIG(J+1,I)-YPSIG(J,I))                &
                          & -(PDT/RDX)*(YZDOT(K,J,I)/(PM0(J+1,I)+PM0(J,I)))*(                        &
                          & (0.25d0*((PDETAF(J)/PDETAH(J))*(YXDOT(K,J+1,I)                           &
                          & +YXDOT(K,J+1,NPERIO(I-1)))+(PDETAF(J+1)/PDETAH(J))*(                     &
                          & YXDOT(K,J,I)+YXDOT(K,J,NPERIO(I-1)))))*0.25d0*(                          &
                          & ((PDETAF(J)/PDETAH(J))*(YPSIG(J+1,NPERIO(I+1))                           &
                          & -YPSIG(J+1,NPERIO(I-1)))+(PDETAF(J+1)/PDETAH(J))*(                       &
                          &  YPSIG(J,NPERIO(I+1))-YPSIG(J,NPERIO(I-1))))))                           &
                          & -HALF*PDT*(YZDOT(K,J,I)/(PM0(J+1,I)+PM0(J,I)))*(                         &
                          &  HALF*(YPSIG(J+1,I)+YPSIG(J,I))*(((((PDETAF(J)/PDETAH(J))*(              &
                          &  YXDOT(K,J+1,I)-YXDOT(K,J+1,NPERIO(I-1)))+(PDETAF(J+1)/PDETAH(J))*(      &
                          &  YXDOT(K,J,I)-YXDOT(K,J,NPERIO(I-1))))/RDX)+(((PDETAF(J)/PDETAF(J+1))*(  &
                          &  YZDOT(K,J+1,I)-YZDOT(K,J,I))+(PDETAF(J+1)/PDETAF(J))*(                  &
                          &  YZDOT(K,J,I)-YZDOT(K,J-1,I)))/PDETAH(J)))                               &
                          & +((PM1(J+1,I)-PM0(J+1,I)+PM1(J,I)-PM0(J,I))/PDT)) ) 
              END DO
              YZDOT(K+1,0,I)    = ZERO
              YZDOT(K+1,NLEV,I) = ZERO
           END DO
        END IF    
        !***************************************************************************************!
        !* GLOBAL FCT NON-OSCILLATORY TREATMENT
        IF (LFCT_MONO) THEN
           DO I=1,NXPT
              DO J=1,NLEV

                 ZMAX        = MAX( YPSIG(J,NPERIO(I+1)),YPSIG(NBC(J+1),I),YPSIG(J,I),               &
                             & YPSIG(J,NPERIO(I-1)),YPSIG(NBC(J-1),I),YPSI(0,NBC(J+1),I),            &
                             & YPSI(0,J,NPERIO(I+1)),YPSI(0,J,I),YPSI(0,NBC(J-1),I),                 &
                             & YPSI(0,J,NPERIO(I-1)))

                 ZIN(J,I)   = ( ZMAX-YPSIG(J,I))/( (ONE/PM0(J,I))*(                                  & 
                             & (PDT/RDX)*(MAX(ZERO,YXDOT(K+1,J,NPERIO(I-1)))                         & 
                             & -MIN(ZERO,YXDOT(K+1,J,I)))                                            &
                             & +(PDT/PDETAF(J))*(MAX(ZERO,YZDOT(K+1,J-1,I))                          &
                             & -MIN(ZERO,YZDOT(K+1,J,I))) ) + REPS  )


                 ZMIN        =  MIN( YPSIG(J,NPERIO(I+1)),YPSIG(NBC(J+1),I),YPSIG(J,I),              &
                             & YPSIG(J,NPERIO(I-1)),YPSIG(NBC(J-1),I),YPSI(0,NBC(J+1),I),            &
                             & YPSI(0,J,NPERIO(I+1)),YPSI(0,J,I),YPSI(0,NBC(J-1),I),                 &
                             & YPSI(0,J,NPERIO(I-1)))
                 
                 ZOUT(J,I)  = (YPSIG(J,I) - ZMIN)/( (ONE/PM0(J,I))*(                                 &
                             & (PDT/RDX)*(MAX(ZERO,YXDOT(K+1,J,I))                                   &
                             & -MIN(ZERO,YXDOT(K+1,J,I)))                                            &
                             & +(PDT/PDETAF(J))*(MAX(ZERO,YZDOT(K+1,J,I))                            &
                             & -MIN(ZERO,YZDOT(K+1,J-1,I))) ) + REPS  )
                   
              END DO
           END DO
           DO I=1,NXPT
              DO J=1,NLEV
                 YXDOT(K+1,J,I) = MIN(ONE,ZOUT(J,I),ZIN(J,NPERIO(I+1)))*MAX(ZERO,YXDOT(K+1,J,I))     & 
                                & +MIN(ONE,ZIN(J,I),ZOUT(J,NPERIO(I+1)))*MIN(ZERO,YXDOT(K+1,J,I))
              END DO
              DO J=1,NLEV-1
                 YZDOT(K+1,J,I) =MIN(ONE,ZOUT(J,I),ZIN(J+1,I))*MAX(ZERO,YZDOT(K+1,J,I))              & 
                                & + MIN(ONE,ZIN(J,I),ZOUT(J+1,I)*MIN(ZERO,YZDOT(K+1,J,I))
              END DO
           END DO      
        END IF
        !********************************************************************************************!
        !* CORRECTED UPWIND SCHEME
        DO I=1,NXPT   
           DO J=1,NLEV
              YPSI(K+1,J,I) = YPSI(K,J,I) - (ONE/PM0(J,I))*(                                         &
                            &  (PDT/RDX)*(YXDOT(K+1,J,I)-YXDOT(K+1,J,NPERIO(I-1)))                   &
                            & +(PDT/PDETA(J))*(YZDOT(K+1,J,I)-YZDOT(K+1,J-1,I)) )           
              YPSIG(J,I)    = (PM0(J,I)/PM1(J,I))*YPSI(K+1,J,I)
           END DO   
        END DO
     END DO
  END IF
  
  !*******************************!
  DO I=1,NXPT
     DO J=1,NLEV
        PPSI_DEP(J,I)  = (PM0(J,I)/PM1(J,I))*(YPSI(NORD2,J,I)+ZPSI_MIN)
     END DO   
     IF (LLMASS) THEN
        DO J=1,NLEV
           PXDOT(J,I)  = HALF*YXDOT(NORD2,J,I)*(                                                     &
                       & HALF*(PPSI_ARR(J,NPERIO(I+1))+PPSI_DEP(J,I))                                &
                       & + HALF*(PPSI_ARR(J,I)+PPSI_DEP(J,NPERIO(I+1))) )
        END DO
        PZDOT(0,I)     = ZERO
        PZDOT(NLEV,I)  = ZERO
        DO J=1,NLEV-1
           PZDOT(J,I)  = HALF*YZDOT(NORD2,J,I)*(                                                     &
                       & HALF*(PPSI_ARR(J+1,I)+PPSI_DEP(J+1,I))                                      &
                       & + HALF*(PPSI_ARR(J,I)+PPSI_DEP(J,I)) )
        END DO
     END IF     
  END DO 
  
END SUBROUTINE ADV_XZ
!#####################################################################################################!
!*****************************************************************************************************!
!*****************************************************************************************************!
!*****************************************************************************************************!
!*****************************************************************************************************!
!*****************************************************************************************************!
!*****************************************************************************************************!
!#####################################################################################################!
SUBROUTINE MPDATA_SURF_PRESSURE(PPIS_DEP,PPIS_ARR,PXDOTS,PIS_MIN)
  
  IMPLICIT NONE

  REAL(8),                       INTENT(OUT)    :: PPIS_DEP(NLEV:NLEV,NXPT)
  REAL(8),                       INTENT(IN)     :: PPIS_ARR(NLEV:NLEV,NXPT)
  REAL(8),                       INTENT(IN)     :: PXDOTS(NXPT)
  REAL(8),                       INTENT(IN)     :: PIS_MIN

  INTEGER(8)                                    :: I,K
  REAL(8)                                       :: YXDOTS(NORD2,NXPT)
  REAL(8)                                       :: YPIS(0:NORD2,NXPT)
  REAL(8)                                       :: ZIN_UP,ZOUT_UP
  REAL(8)                                       :: ZIN_DOWN,ZOUT_DOWN
  REAL(8)                                       :: ZPIS_MIN,ZXDOTS_FCT(NXPT)
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ZPIS_MIN = MIN(ZERO,PIS_MIN,MINVAL(PPIS_ARR))
  
  DO I=1,NXPT
     YXDOTS(1,I)  = PXDOTS(I)
     YPIS(0,I)    = PPIS_ARR(NLEV,I) - ZPIS_MIN
  END DO
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !***************************************!
  ! DONOR-CELL UPWIND FIRST-ORDER SCHEME *!
  !***************************************!

  DO I=1,NXPT
     YPIS(1,I)        = YPIS(0,I) - (RDT/RDX)*(                                             &
                      &  (MIN(ZERO,YXDOTS(1,NPERIO(I)))*YPIS(0,NPERIO(I+1))                 &
                      &  +MAX(ZERO,YXDOTS(1,NPERIO(I)))*YPIS(0,I))                          &
                      & -(MIN(ZERO,YXDOTS(1,NPERIO(I-1)))*YPIS(0,I)                         &
                      &  +MAX(ZERO,YXDOTS(1,NPERIO(I-1)))*YPIS(0,NPERIO(I-1))) )
  END DO
  
  !************************!
  ! FIRST CORRECTION-STEP *!
  !************************!

  IF (.NOT.LGAUGE) THEN

     DO K=1,NORD1-1
     !**************************************************************************************!
     !* COMPUTE SURFACE PSEUDO-VELOCITIES
        IF (LMPDATA_FV) THEN
           DO I=1,NXPT
              YXDOTS(K+1,I) = HALF*ABS(YXDOTS(K,I))*(                                      & 
                         &  (ABS(YPIS(K,NPERIO(I+1)))-ABS(YPIS(K,I)))                      &
                         &  /(ABS(YPIS(K,I))+ABS(YPIS(K,NPERIO(I+1)))+REPS) )              &
                         & -(RDT/TWO)*YXDOTS(K,I)*HALF*(                                   &
                         &  (HALF*( ( YXDOTS(K,I)*(                                        &
                         & ABS(YPIS(K,NPERIO(I+1)))+ABS(YPIS(K,I)))                        &
                         & -YXDOTS(K,NPERIO(I-1))*(ABS(YPIS(K,NPERIO(I-1)))                &   
                         & +ABS(YPIS(K,I))) )/RDX)/( (0.25d0*(                             &
                         &  ABS(YPIS(K,NPERIO(I+1)))+ABS(YPIS(K,I))                        &
                         & +ABS(YPIS(K,NPERIO(I-1)))+ABS(YPIS(K,I))))+REPS))               &    
                         & + (HALF*( ( YXDOTS(K,NPERIO(I+1))*(                             &
                         & ABS(YPIS(K,NPERIO(I+2)))+ABS(YPIS(K,NPERIO(I+1))))              &
                         & -YXDOTS(K,I)*(ABS(YPIS(K,NPERIO(I+1)))                          &
                         & +ABS(YPIS(K,NPERIO(I)))) )/RDX)/( (0.25d0*(                     &  
                         &  ABS(YPIS(K,NPERIO(I+2)))+ABS(YPIS(K,NPERIO(I+1)))              & 
                         & +ABS(YPIS(K,NPERIO(I+1)))+ABS(YPIS(K,I))))+REPS))   )
           END DO
        ELSE
           DO I=1,NXPT
              YXDOTS(K+1,I) = HALF*ABS(YXDOTS(K,I))*(                                      & 
                         &  (ABS(YPIS(K,NPERIO(I+1)))-ABS(YPIS(K,I)))                      &
                         &  /(ABS(YPIS(K,I))+ABS(YPIS(K,NPERIO(I+1)))+REPS) )              &
                         & -(RDT/TWO)*YXDOTS(K,I)*HALF*(                                   &
                         &  (HALF*( ( YXDOTS(K,I)*(                                        &
                         & ABS(YPIS(K,NPERIO(I+1)))+ABS(YPIS(K,I)))                        &
                         & -YXDOTS(K,NPERIO(I-1))*(ABS(YPIS(K,NPERIO(I-1)))                &   
                         & +ABS(YPIS(K,I))) )/RDX)/( (0.25d0*(                             &
                         &  ABS(YPIS(K,NPERIO(I+1)))+ABS(YPIS(K,I))                        &
                         & +ABS(YPIS(K,NPERIO(I-1)))+ABS(YPIS(K,I))))+REPS))               &    
                         & + (HALF*( ( YXDOTS(K,NPERIO(I+1))*(                             &
                         & ABS(YPIS(K,NPERIO(I+2)))+ABS(YPIS(K,NPERIO(I+1))))              &
                         & -YXDOTS(K,I)*(ABS(YPIS(K,NPERIO(I+1)))                          &
                         & +ABS(YPIS(K,NPERIO(I)))) )/RDX)/( (0.25d0*(                     &  
                         &  ABS(YPIS(K,NPERIO(I+2)))+ABS(YPIS(K,NPERIO(I+1)))              & 
                         & +ABS(YPIS(K,NPERIO(I+1)))+ABS(YPIS(K,I))))+REPS))   )
           END DO
        END IF   
        !**************************************************************************************!
        !* GLOBAL FCT NON-OSCILLATORY TREATMENT
        IF (LFCT_MONO) THEN

           DO I=1,NXPT
              ZIN_DOWN   = ( MAX( YPIS(K,NPERIO(I+1)),YPIS(K,NPERIO(I)),                    &
                         & YPIS(K,NPERIO(I-1)),YPIS(0,NPERIO(I+1)),                         &
                         & YPIS(0,NPERIO(I)),YPIS(0,NPERIO(I-1)))-YPIS(K,I) )               &
                         & /( (RDT/RDX)*(                                                   &
                         & MAX(ZERO,YXDOTS(K+1,NPERIO(I-1)))*YPIS(K,NPERIO(I-1))            &
                         & -MIN(ZERO,YXDOTS(K+1,I))*YPIS(K,NPERIO(I+1)))+REPS)

              ZIN_UP     = ( MAX( YPIS(K,NPERIO(I+2)),YPIS(K,I),                            &
                         & YPIS(K,NPERIO(I+1)),YPIS(0,NPERIO(I+2)),                         &
                         & YPIS(0,I),YPIS(0,NPERIO(I+1)))-YPIS(K,NPERIO(I+1)) )             &
                         & /( (RDT/RDX)*( MAX(ZERO,YXDOTS(K+1,I))*YPIS(K,I)                 &
                         & -MIN(ZERO,YXDOTS(K+1,NPERIO(I+1)))*YPIS(K,NPERIO(I+2)))+REPS)
        
              ZOUT_DOWN  = -( MIN( YPIS(K,NPERIO(I+1)),YPIS(K,NPERIO(I)),                   &
                         & YPIS(K,NPERIO(I-1)),YPIS(0,NPERIO(I+1)),                         &
                         & YPIS(0,NPERIO(I)),YPIS(0,NPERIO(I-1)))-YPIS(K,I) )               &
                         & /( (RDT/RDX)*( MAX(ZERO,YXDOTS(K+1,I))*YPIS(K,I)                 &
                         & -MIN(ZERO,YXDOTS(K+1,NPERIO(I-1)))*YPIS(K,I))+REPS)

              ZOUT_UP    = -( MIN( YPIS(K,NPERIO(I+2)),YPIS(K,I),                           &
                         & YPIS(K,NPERIO(I+1)),YPIS(0,NPERIO(I+2)),                         &
                         & YPIS(0,I),YPIS(0,NPERIO(I+1)))-YPIS(K,NPERIO(I+1)) )             &
                         & /( (RDT/RDX)*(                                                   &
                         & MAX(ZERO,YXDOTS(K+1,NPERIO(I+1)))*YPIS(K,NPERIO(I+1))            &
                         & -MIN(ZERO,YXDOTS(K+1,I))*YPIS(K,NPERIO(I+1)))+REPS)
        
              ZXDOTS_MONO =  MIN(ONE,ZOUT_DOWN,ZIN_UP)*MAX(ZERO,YXDOTS(K+1,I))               &      
                          & +MIN(ONE,ZIN_DOWN,ZOUT_UP)*MIN(ZERO,YXDOTS(K+1,I))              
              YXDOTS(K+1,I) = (ONE-RFCT)*YXDOTS(K+1,I)+RFCT*ZXDOTS_MONO
           END DO                                                                                  
        END IF
     
        !**************************************************************************************!
        !* CORRECTED UPWIND SCHEME
        DO I=1,NXPT   
           YPIS(K+1,I)   = YPIS(K,I) - (RDT/RDX)*(                                          &
                         & (MIN(ZERO,YXDOTS(K+1,NPERIO(I)))*YPIS(K,NPERIO(I+1))             &
                         & +MAX(ZERO,YXDOTS(K+1,NPERIO(I)))*YPIS(K,I))                      &
                         & -(MIN(ZERO,YXDOTS(K+1,NPERIO(I-1)))*YPIS(K,I)                    &
                         & +MAX(ZERO,YXDOTS(K+1,NPERIO(I-1)))*YPIS(K,NPERIO(I-1))) )
        END DO 
     END DO

  ELSE
     
  !*************************!
  ! SECOND CORRECTION-STEP *!
  !*************************!
  
     DO K=1,NORD-1
     !**************************************************************************************!
     !* COMPUTE SURFACE PSEUDO-VELOCITIES
        IF (LMPDATA_FV) THEN 
           DO I=1,NXPT
              YXDOTS(K+1,I) = HALF*ABS(YXDOTS(K,I))*(YPIS(K,NPERIO(I+1))-YPIS(K,I))         &
                         & -(RDT/TWO)*YXDOTS(K,I)*HALF*(                                    &
                         & (HALF*((YXDOTS(K,I)*(YPIS(K,NPERIO(I+1))+YPIS(K,I))              &
                         & -YXDOTS(K,NPERIO(I-1))*(YPIS(K,NPERIO(I-1))+YPIS(K,I)) )         &
                         & /RDX))+(HALF*((YXDOTS(K,NPERIO(I+1))*(YPIS(K,NPERIO(I+2))        &
                         & +YPIS(K,NPERIO(I+1)))-YXDOTS(K,I)*(YPIS(K,NPERIO(I+1))           &
                         & +YPIS(K,NPERIO(I))) )/RDX)) )
           END DO
        ELSE
           DO I=1,NXPT
              YXDOTS(K+1,I) = HALF*ABS(YXDOTS(K,I))*(YPIS(K,NPERIO(I+1))-YPIS(K,I))         &
                         & -(RDT/TWO)*YXDOTS(K,I)*HALF*(                                    &
                         & (HALF*((YXDOTS(K,I)*(YPIS(K,NPERIO(I+1))+YPIS(K,I))              &
                         & -YXDOTS(K,NPERIO(I-1))*(YPIS(K,NPERIO(I-1))+YPIS(K,I)) )         &
                         & /RDX))+(HALF*((YXDOTS(K,NPERIO(I+1))*(YPIS(K,NPERIO(I+2))        &
                         & +YPIS(K,NPERIO(I+1)))-YXDOTS(K,I)*(YPIS(K,NPERIO(I+1))           &
                         & +YPIS(K,NPERIO(I))) )/RDX)) )
           END DO
        END IF   
        !**************************************************************************************!
        !* GLOBAL FCT NON-OSCILLATORY TREATMENT
        IF (LFCT_MONO) THEN
           DO I=1,NXPT
              ZIN_DOWN   = ( MAX( YPIS(K,NPERIO(I+1)),YPIS(K,NPERIO(I)),                    &
                         & YPIS(K,NPERIO(I-1)),YPIS(0,NPERIO(I+1)),                         &
                         & YPIS(0,NPERIO(I)),YPIS(0,NPERIO(I-1)))-YPIS(K,I) )               &
                         & /( (RDT/RDX)*(                                                   &
                         & MAX(ZERO,YXDOTS(K+1,NPERIO(I-1)))*YPIS(K,NPERIO(I-1))            &
                         & -MIN(ZERO,YXDOTS(K+1,I))*YPIS(K,NPERIO(I+1)))+REPS)

              ZIN_UP     = ( MAX( YPIS(K,NPERIO(I+2)),YPIS(K,I),                            &
                         & YPIS(K,NPERIO(I+1)),YPIS(0,NPERIO(I+2)),                         &
                         & YPIS(0,I),YPIS(0,NPERIO(I+1)))-YPIS(K,NPERIO(I+1)) )             &
                         & /( (RDT/RDX)*( MAX(ZERO,YXDOTS(K+1,I))*YPIS(K,I)                 &
                         & -MIN(ZERO,YXDOTS(K+1,NPERIO(I+1)))*YPIS(K,NPERIO(I+2)))+REPS)
        
              ZOUT_DOWN  = -( MIN( YPIS(K,NPERIO(I+1)),YPIS(K,NPERIO(I)),                   &
                         & YPIS(K,NPERIO(I-1)),YPIS(0,NPERIO(I+1)),                         &
                         & YPIS(0,NPERIO(I)),YPIS(0,NPERIO(I-1)))-YPIS(K,I) )               &
                         & /( (RDT/RDX)*( MAX(ZERO,YXDOTS(K+1,I))*YPIS(K,I)                 &
                         & -MIN(ZERO,YXDOTS(K+1,NPERIO(I-1)))*YPIS(K,I))+REPS)

              ZOUT_UP    = -( MIN( YPIS(K,NPERIO(I+2)),YPIS(K,I),                           &
                         & YPIS(K,NPERIO(I+1)),YPIS(0,NPERIO(I+2)),                         &
                         & YPIS(0,I),YPIS(0,NPERIO(I+1)))-YPIS(K,NPERIO(I+1)) )             &
                         & /( (RDT/RDX)*(                                                   &
                         & MAX(ZERO,YXDOTS(K+1,NPERIO(I+1)))*YPIS(K,NPERIO(I+1))            &
                         & -MIN(ZERO,YXDOTS(K+1,I))*YPIS(K,NPERIO(I+1)))+REPS)
        
              ZXDOTS_FCT(I) =  MIN(ONE,ZOUT_DOWN,ZIN_UP)*MAX(ZERO,YXDOTS(K+1,I))            &      
                            & +MIN(ONE,ZIN_DOWN,ZOUT_UP)*MIN(ZERO,YXDOTS(K+1,I))              
           END DO
           YXDOTS(K+1,:) = (ONE-RFCT)*YXDOTS(K+1,:)+RFCT*ZXDOTS_FCT(:)
        END IF
     
        !************************************************************************************!
        !* CORRECTED UPWIND SCHEME
        DO I=1,NXPT   
           YPIS(K+1,I)   = YPIS(K,I) - (RDT/RDX)*( YXDOTS(K+1,I)-YXDOTS(K+1,NPERIO(I-1)) )
        END DO
     END DO
     
  END IF
  
  DO I=1,NXPT
     PPIS_DEP(NLEV,I) = YPIS(NORD2,I) + ZPIS_MIN 
  END DO

  
END SUBROUTINE MPDATA_SURF_PRESSURE
!################################################################################################!
!################################################################################################!
!################################################################################################!     
!################################################################################################!
!################################################################################################!
!################################################################################################!
!################################################################################################!
!################################################################################################!
!################################################################################################!     
!################################################################################################!
!################################################################################################!
!################################################################################################!
!################################################################################################!
!************************************************************************************************!
END MODULE MOD_MPDTA
