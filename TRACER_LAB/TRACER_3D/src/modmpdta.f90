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

  USE MPI
  USE OMP_LIB
  USE MOD_PARAMETER, ONLY : NXPT,NLEV,RR,RG,RCP,RCV,RPI,ONE,TWO,ZERO,HALF,RDT,RDX, &
                          & LSETTLS,NPERIO,LPERIO,NORD1,NORD2,RFCT,REPS,NBC,NBL,   &
                          & LPRINTLEV,RVMAX,LADV_PIS2D,LFCT_MONO
  USE MOD_STRUCTURE, ONLY : MPDATA,RHS_LOC,RHS_GLB,GEOMETRY
  USE MOD_COMM,      ONLY : TRANSGLB_RHS,WIND_GLB_TRANSFER,TRANSLOC_RHS,TRANSGLB
 
CONTAINS

  !****************************************************************************************!
  !****************************************************************************************!
  !****************************************************************************************!
  !****************************************************************************************!
  !****************************************************************************************!
  SUBROUTINE MPDATA_TRANSPORT_SCHEME(RHSL_DEP,RHSL_ARR,RHSG_DEP,RHSG_ARR,   &
             & MPDTA,PMF,PMH,PXDOTF0,PYDOTF0,PZDOTF0,PZDOTH0,PETAF,PDETAF,  &
             & PETAH,PDETAH,PDELTB,KITER,KDIM,KSTEP,ABS_GLB,ORD_GLB,MYPROC, &
             & CODE,CPART)
  
  IMPLICIT NONE

  TYPE(RHS_LOC),            INTENT(OUT)    :: RHSL_DEP
  TYPE(RHS_LOC),            INTENT(IN)     :: RHSL_ARR
  TYPE(RHS_GLB),            INTENT(INOUT)  :: RHSG_ARR
  TYPE(RHS_GLB),            INTENT(INOUT)  :: RHSG_DEP
  TYPE(MPDATA),             INTENT(INOUT)  :: MPDTA
  REAL(8),                  INTENT(IN)     :: PMF(NLEV,KDIM),PMH(0:NLEV,KDIM)  
  REAL(8),                  INTENT(IN)     :: PXDOTF0(NLEV,KDIM)
  REAL(8),                  INTENT(IN)     :: PYDOTF0(NLEV,KDIM)
  REAL(8),                  INTENT(IN)     :: PZDOTF0(NLEV,KDIM)
  REAL(8),                  INTENT(IN)     :: PZDOTH0(0:NLEV,KDIM)
  REAL(8),                  INTENT(IN)     :: PDETAF(NLEV),PETAF(NLEV)
  REAL(8),                  INTENT(IN)     :: PDETAH(0:NLEV),PETAH(0:NLEV)
  REAL(8),                  INTENT(IN)     :: PDELTB(NLEV)
  INTEGER(8),               INTENT(IN)     :: KITER,KSTEP,KDIM
  INTEGER,                  INTENT(INOUT)  :: CODE
  INTEGER,                  INTENT(IN)     :: MYPROC
  INTEGER(8),               INTENT(IN)     :: ABS_GLB(0:NB_PROCS-1)
  INTEGER(8),               INTENT(IN)     :: ORD_GLB(0:NB_PROCS-1)
  CHARACTER(LEN=1),         INTENT(IN)     :: CPART

  REAL(8), DIMENSION(NLEV,GXPT,GYPT)       :: ZXDOTF0,ZYDOTF0
  REAL(8)                                  :: ZLIN,ZWEI
  REAL(8)                                  :: ZWIND_HOR_MAX,ZWIND_VER_MAX
  INTEGER(8)                               :: I,J,K,L

  
  ZLIN = ONE
  IF (CPART == 'LI') ZLIN =ZERO

  IF (KITER == 0) THEN
     
     CALL TRANSGLB_RHS(RHSG_ARR,RHSL_ARR,      &
          & ABS_GLB,ORD_GLB,MYPROC,CODE,1)
     
     CALL TRANSGLB(MPDTA%MF(0,:,:,:),PMF,1,    &
          & ABS_GLB,ORD_GLB,MYPROC,CODE,10)
     CALL TRANSGLB(MPDTA%MH(0,:,:,:),PMH,0,    &
          & ABS_GLB,ORD_GLB,MYPROC,CODE,11)
     CALL TRANSGLB(MPDTA%ZDOTH0,PZDOTH0,0,     &
          & ABS_GLB,ORD_GLB,MYPROC,CODE,12)
     
     CALL WIND_GLB_TRANSFER(ZXDOTF,ZYDOTF,     &
          & MPDTA%ZDOTF0,PXDOTF0,PYDOTF0,      &
          &PZDOTF0,1,ABS_GLB,ORD_GLB,MYPROC,   &
          & CODE,20)
  ELSE

     CALL TRANSGLB(MPDTA%ZDOTH(1,:,:,:),       &
          & PZDOTH0,0,ABS_GLB,ORD_GLB,MYPROC,  &
          & CODE,12)
     
     CALL WIND_GLB_TRANSFER(ZXDOTF,ZYDOTF,     &
          & MPDTA%ZDOTF(1,:,:,:),PXDOTF0,      &
          & PYDOTF0,PZDOTF0,1,ABS_GLB,ORD_GLB, &
          & MYPROC,CODE,20)
     
  ENDIF

  !*********************************************************************
  IF (MYPROC == 0) THEN
    
     !**************************************************!
     !* Advecting velocity at interfaces at (t+0.5*dt) *!
     !**************************************************!
     
     IF (KITER == 0) THEN
        DO I=1,GXPT
           DO J=1,GYPT
              MPDTA%ZDOTF(1,:,I,J) = ZLIN*MPDTA%ZDOTF0(:,I,J) 
              MPDTA%ZDOTH(1,:,I,J) = ZLIN*MPDTA%ZDOTH0(:,I,J)
              MPDTA%XDOTF0(:,I,J)  = ZLIN*HALF*(ZXDOTF(:,I,J)         &
                                   & + ZXDOTF(:,NPERIOX(I+1),J) )
              MPDTA%YDOTF0(:,I,J)  = ZLIN*HALF*(ZYDOTF(:,I,J)         &
                                   & + ZYDOTF(:,I,NPERIOY(J+1)) )
              MPDTA%XDOTF(:,I,J)   = ZLIN*MPDTA%XDOTF0(:,I,J)
              MPDTA%YDOTF(:,I,J)   = ZLIN*MPDTA%YDOTF0(:,I,J)
              
           END DO
        END DO
           
     ELSE
        DO I=1,GXPT
           DO J=1,GYPT
              MPDTA%XDOTF(:,I,J)   = ZLIN*HALF*( MPDTA%XDOTF0(:,I,J)  &
                                   & + HALF*(ZXDOTF(:,NPERIOX(I+1),J) &
                                   & + ZXDOTF(:,I,J)) )
              MPDTA%YDOTF(:,I,J)   = ZLIN*HALF*( MPDTA%YDOTF0(:,I,J)  &
                                   & + HALF*(ZYDOTF(:,I,NPERIOY(J+1)) &
                                   & + ZYDOTF(:,I,J)) )
              MPDTA%ZDOTF(1,:,I,J) = ZLIN*( MPDTA%ZDOTF(1,:,I,J)      &
                                   & + MPDTA%ZDOTF0(:,I,J) )/TWO
              MPDTA%ZDOTH(1,:,I,J) = ZLIN*( MPDTA%ZDOTH(1,:,I,J)      &
                                   & + MPDTA%ZDOTH0(:,I,J) )/TWO
           END DO
        END DO
     ENDIF

     !************************************************!
     !* special treatment at half-levels and surface *!
     !************************************************!
  
     DO I=1,GXPT
        DO J=1,GYPT

           MPDTA%XDOTH(0,I,J)      = MPDTA%XDOTF(1,I,J)
           MPDTA%YDOTH(0,I,J)      = MPDTA%YDOTF(1,I,J)
           DO K=1,NLEV-1
              ZWEI                 = ( (PETAH(K)-PETAF(K))         &
                                   & /(PETAF(K+1)-PETAF(K)) ) 
              MPDTA%XDOTH(K,I,J)   = MPDTA%XDOTF(K,I,J)            &
                                   & + ZWEI*( MPDTA%XDOTF(K+1,I,J) &
                                   & - MPDTA%XDOTF(K,I,J) )
              MPDTA%YDOTH(K,I,J)   = MPDTA%YDOTF(K,I,J)            &
                                   & + ZWEI*( MPDTA%YDOTF(K+1,I,J) &
                                   & - MPDTA%YDOTF(K,I,J) )
           ENDDO
           MPDTA%XDOTH(NLEV,I,J)   = MPDTA%XDOTF(NLEV,I,J)
           MPDTA%YDOTH(NLEV,I,J)   = MPDTA%YDOTF(NLEV,I,J)
          
          
           MPDTA%XDOTS(I,J)        = ZERO
           MPDTA%YDOTS(I,J)        = ZERO
           DO K=1,NLEV
              MPDTA%XDOTS(I,J)     = MPDTA%XDOTS(I,J)              &
                                   & + PDELTB(K)*MPDTA%XDOTF(K,J,I)
              MPDTA%YDOTS(I,J)     = MPDTA%YDOTS(I,J)              &
                                   & + PDELTB(K)*MPDTA%YDOTF(K,J,I)
           END DO
           
        END DO   
     END DO
  
     !*****************************!
     !*  masses transport scheme  *!
     !*****************************!
  
     CALL MPDATA_FULL_MASS(MPDTA%MF,MPDTA%XDOTF,MPDTA%YDOTF, &
          & MPDTA%ZDOTH,PDETAF) 
     CALL MPDATA_HALF_MASS(MPDTA%MH,MPDTA%XDOTH,MPDTA%YDOTH, &
          & MPDTA%ZDOTF,PDETAH)
     
     !****************************************************************!
     ! STABILITY CHECK-POINT MAXIMUM OF HORIZONTAL WIND VELOCITY      !
     !****************************************************************!
     ZWIND_HOR_MAX = MAX(MAXVAL(MPDTA%XDOTF0),MAXVAL(MPDTA%YDOTF0))   !
     ZWIND_VER_MAX = MAXVAL(MPDTA%ZDOTH0)                             !
                                                                      !
     IF (LPRINTLEV) THEN                                              !
        IF (KITER == 0) WRITE(*,*) 'NSTEP : ',KSTEP,         &        !
                      & '--> MAX.WIND_HOR : ',ZWIND_HOR_MAX, &        !
                      & '--> MAX.WIND_VER : ',ZWIND_VER_MAX           !
     END IF                                                           !
     IF (ZWIND_HOR_MAX .GE. 300.d0) THEN                              !
        WRITE(*,*) 'WIND TOO STRONG EXPLOSION AT TIMESTEP = ',KSTEP   !
        STOP                                                          !
     END IF                                                           !
     !****************************************************************!

     
     IF (LADV_PIS2D) THEN
        CALL MPDATA_SURF_PRESSURE(RHS_DEP%PIS,RHS_ARR%PIS,     &
             & MPDTA%XDOTS,MPDTA%YDOTS)
     ELSE 
        DO I=1,GXPT
           DO J=1,GYPT
              RHS_DEP%PIS(I,J)    = ZERO
              DO K=1,NLEV
                 RHS_DEP%PIS(I,J) = RHS_DEP%PIS(I,J)           &
                                  & + MPDTA%MF(3,K,I,J)*PDETAF(K)
              END DO
           END DO   
        END DO   
     END IF
  
     CALL MPDATA_2D_FULL_FIELD(RHS_DEP%T,RHS_ARR%T,            &
          & MPDTA%XDOTF,MPDTA%YDOTF,MPDTA%ZDOTH,MPDTA%MF,PDETAF)
     CALL MPDATA_2D_HALF_FIELD(RHS_DEP%Q,RHS_ARR%Q,            &
          & MPDTA%XDOTH,MPDTA%YDOTH,MPDTA%ZDOTF,MPDTA%MH,PDETAH)

  END IF 

  !*************************************!
  !waiting for the end of myproc=0 task !
  CALL MPI_BARRIER(MPI_COMM_WORLD,CODE) !
  !*************************************!
   
  CALL TRANSLOC_RHS(RHSL_DEP,RHSG_DEP,ABS_GLB,ORD_GLB,MYPROC,CODE,1)
  
END SUBROUTINE MPDATA_TRANSPORT_SCHEME
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
!################################################################################################!
SUBROUTINE MPDATA_FULL_MASS(PMF,PXDOTF,PYDOTF,PZDOTH,PDETAF)

  IMPLICIT NONE

  REAL(8),                       INTENT(INOUT)  :: PMF(0:3,NLEV,GXPT,GYPT)
  REAL(8),                       INTENT(INOUT)  :: PXDOTF(NLEV,GXPT,GYPT)
  REAL(8),                       INTENT(INOUT)  :: PYDOTF(NLEV,GXPT,GYPT)
  REAL(8),                       INTENT(INOUT)  :: PZDOTH(3,0:NLEV,GXPT,GYPT)
  REAL(8),                       INTENT(IN)     :: PDETAF(NLEV)

  REAL(8)                                       :: ZUN(NLEV,GXPT,GYPT)

  ZUN(:,:,:) = ONE
  
  CALL ADV_ZF(PMF(1,:,:,:),PMF(0,:,:,:),PZDOTH, &
       & ZUN,ZUN,PDETAF,RDT/TWO,1,LDMASS=.TRUE.)
  CALL ADV_XY(PMF(2,:,:,:),PMF(1,:,:,:),PXDOTF, &
       & PYDOTF,ZUN,ZUN,RDT,1_8,LDMASS=.TRUE.)
  CALL ADV_ZF(PMF(3,:,:,:),PMF(2,:,:,:),PZDOTH, &
       & ZUN,ZUN,PDETAF,RDT/TWO,2,LDMASS=.TRUE.)     
  
END SUBROUTINE MPDATA_FULL_MASS
!################################################################################################!
!************************************************************************************************!
!************************************************************************************************!
!************************************************************************************************!
!################################################################################################!
SUBROUTINE MPDATA_FULL_FIELD(PPSI_DEP,PPSI_ARR,PXDOTF,PYDOTF,PZDOTH,PMF,PDETAF)

  IMPLICIT NONE

  REAL(8),          INTENT(OUT)    :: PPSI_DEP(NLEV,GXPT,GYPT)
  REAL(8),          INTENT(IN)     :: PPSI_ARR(NLEV,GXPT,GYPT)
  REAL(8),          INTENT(INOUT)  :: PXDOTF(NLEV,GXPT,GYPT)
  REAL(8),          INTENT(INOUT)  :: PYDOTF(NLEV,GXPT,GYPT)
  REAL(8),          INTENT(INOUT)  :: PZDOTH(3,0:NLEV,GXPT,GYPT)
  REAL(8),          INTENT(IN)     :: PMF(0:3,NLEV,GXPT,GYPT)
  REAL(8),          INTENT(IN)     :: PDETAF(NLEV)
  
  REAL(8)                          :: YPSI(3,NLEV,GXPT,GYPT)
  INTEGER(8)                       :: I,J,K

  CALL ADV_ZF(YPSI(1,:,:,:),PPSI_ARR,PZDOTH,       &
       & PMF(1,:,:,:),PMF(0,:,:,:),PDETAF,RDT/TWO,1)
  CALL ADV_XY(YPSI(2,:,:,:),YPSI(1,:,:,:),PXDOTF,  &
       & PYDOTF,PMF(2,:,:,:),PMF(1,:,:,:),RDT,1_8)
  CALL ADV_ZF(YPSI(3,:,:,:),YPSI(2,:,:,:),PZDOTH,  &
       & PMF(3,:,:,:),PMF(2,:,:,:),PDETAF,RDT/TWO,2)

  DO I=1,GXPT
     DO J=1,GYPT
        DO K=1,NLEV
           PPSI_DEP(K,I,J) = &
                & (PMF(0,K,I,J)/PMF(3,K,I,J))*YPSI(3,K,I,J)
        END DO
     END DO
  END DO

END SUBROUTINE  MPDATA_FULL_FIELD
!################################################################################################!
!************************************************************************************************!
!************************************************************************************************!
!************************************************************************************************!
!################################################################################################!
SUBROUTINE MPDATA_HALF_MASS(PMH,PXDOTH,PYDOTH,PZDOTF,PDETAH)

  IMPLICIT NONE

  REAL(8),          INTENT(INOUT)  :: PMH(0:3,0:NLEV,GXPT,GYPT)
  REAL(8),          INTENT(INOUT)  :: PXDOTH(0:NLEV,GXPT,GYPT)
  REAL(8),          INTENT(INOUT)  :: PYDOTH(0:NLEV,GXPT,GYPT)
  REAL(8),          INTENT(INOUT)  :: PZDOTF(3,NLEV,GXPT,GYPT)
  REAL(8),          INTENT(IN)     :: PDETAH(0:NLEV)

  REAL(8)                                       :: ZUN(0:NLEV,GXPT,GYPT)

  ZUN(:,:,:) = ONE
  
  CALL ADV_ZH(PMH(1,:,:,:),PMH(0,:,:,:),PZDOTF,  &
       & ZUN,ZUN,PDETAH,RDT/TWO,1,LDMASS=.TRUE.)
  CALL ADV_XY(PMH(2,:,:,:),PMH(1,:,:,:),PXDOTH,  &
       & PYDOTH,ZUN,ZUN,RDT,0_8,LDMASS=.TRUE.)
  CALL ADV_ZH(PMH(3,:,:,:),PMH(2,:,:,:),PZDOTF,  &
       & ZUN,ZUN,PDETAH,RDT/TWO,2,LDMASS=.TRUE.)
  
END SUBROUTINE MPDATA_HALF_MASS
!################################################################################################!
!************************************************************************************************!
!************************************************************************************************!
!************************************************************************************************!
!################################################################################################!
SUBROUTINE MPDATA_HALF_FIELD(PPSI_DEP,PPSI_ARR,PXDOTH,PYDOTH,PZDOTF,PMH,PDETAH)

  IMPLICIT NONE

  REAL(8),        INTENT(OUT)    :: PPSI_DEP(0:NLEV,GXPT,GYPT)
  REAL(8),        INTENT(IN)     :: PPSI_ARR(0:NLEV,GXPT,GYPT)
  REAL(8),        INTENT(INOUT)  :: PXDOTH(0:NLEV,GXPT,GYPT)
  REAL(8),        INTENT(INOUT)  :: PYDOTH(0:NLEV,GXPT,GYPT)
  REAL(8),        INTENT(INOUT)  :: PZDOTF(3,NLEV,GXPT,GYPT)
  REAL(8),        INTENT(IN)     :: PMH(0:3,0:NLEV,GXPT,GYPT)
  REAL(8),        INTENT(IN)     :: PDETAH(0:NLEV)
  
  REAL(8)                        :: YPSI(3,0:NLEV,GXPT,GYPT)
  INTEGER(8)                     :: I,J

  CALL ADV_ZH(YPSI(1,:,:,:),PPSI_ARR,PZDOTF,      &
       & PMH(1,:,:,:),PMH(0,:,:,:),PDETAH,RDT/TWO,1)
  CALL ADV_XY(YPSI(2,:,:,:),YPSI(1,:,:,:),PXDOTH, &
       & PYDOTH,PMH(2,:,:,:),PMH(1,:,:,:),RDT,0_8)
  CALL ADV_ZH(YPSI(3,:,:,:),YPSI(2,:,:,:),PZDOTF, &
       & PMH(3,:,:,:),PMH(2,:,:,:),PDETAH,RDT/TWO,2)     

  DO I=1,GXPT
     DO J=1,GYPT
        DO J=0,NLEV
           PPSI_DEP(K,I,J) =  &
               & (PMH(0,K,I,J)/PMH(3,K,I,J))*YPSI(3,K,I,J)
        END DO
     END DO 
  END DO
  
END SUBROUTINE  MPDATA_HALF_FIELD
!################################################################################################!
!************************************************************************************************!
!************************************************************************************************!
!************************************************************************************************!
!************************************************************************************************!
!************************************************************************************************!
!************************************************************************************************!
!************************************************************************************************!
!************************************************************************************************!
!################################################################################################!
SUBROUTINE ADV_ZF(PPSI_DEP,PPSI_ARR,PZDOT,PM1,PM0,PDETA,PDT,KS,LDMASS)
  
  IMPLICIT NONE

  REAL(8),                       INTENT(OUT)    :: PPSI_DEP(NLEV,GXPT,GYPT)
  REAL(8),                       INTENT(IN)     :: PPSI_ARR(NLEV,GXPT,GYPT)
  REAL(8),                       INTENT(INOUT)  :: PZDOT(3,0:NLEV,GXPT,GYPT)
  REAL(8),                       INTENT(IN)     :: PM1(NLEV,GXPT,GYPT)
  REAL(8),                       INTENT(IN)     :: PM0(NLEV,GXPT,GYPT)
  REAL(8),                       INTENT(IN)     :: PDETA(NLEV)
  REAL(8),                       INTENT(IN)     :: PDT
  INTEGER,                       INTENT(IN)     :: KS
  LOGICAL,OPTIONAL,              INTENT(IN)     :: LDMASS

  REAL(8)                                       :: YPSI(0:NORD2,NLEV,GXPT,GYPT)
  REAL(8)                                       :: YZDOT(NORD2,0:NLEV,GXPT,GYPT)
  REAL(8)                                       :: YPSIG(NLEV,GXPT,GYPT),
  REAL(8)                                       :: ZZDOT_MONO,ZPSI_MIN
  REAL(8)                                       :: ZIN_DOWN,ZOUT_DOWN,ZIN_UP,ZOUT_UP
  INTEGER(8)                                    :: IX,IY,J,K
  LOGICAL                                       :: LLMASS

  IF (PRESENT(LDMASS)) THEN
     LLMASS = LDMASS
  ELSE
     LLMASS = .FALSE.
  END IF   
  
  ZPSI_MIN = MINVAL(PPSI_ARR)

  DO IX=1,GXPT
     DO IY=1,GYPT
        YPSI(0,:,IX,IY)  = PPSI_ARR(:,IX,IY) - ZPSI_MIN
        YZDOT(1,:,IX,IY) = PZDOT(KS,:,IX,IY)
     END DO   
  END DO
 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !***************************************!
  ! DONOR-CELL UPWIND FIRST-ORDER SCHEME *!
  !***************************************!

  DO IX=1,GXPT
     DO IY=1,GYPT
        DO J=1,NLEV
           YPSI(1,J,IX,IY) = YPSI(0,J,IX,IY) - (ONE/PM0(J,IX,IY))*(PDT/PDETA(J))*(   &
                           & (MIN(ZERO,YZDOT(1,J,IX,IY))*YPSI(0,NBC(J+1),IX,IY)      &
                           & +MAX(ZERO,YZDOT(1,J,IX,IY))*YPSI(0,J,IX,IY))            &
                           & -(MIN(ZERO,YZDOT(1,J-1,IX,IY))*YPSI(0,J,IX,IY)          &
                           & +MAX(ZERO,YZDOT(1,J-1,IX,IY))*YPSI(0,NBC(J-1),IX,IY)) )
           YPSIG(J,IX,IY)  = (PM0(J,IX,IY)/PM1(J,IX,IY))*YPSI(1,J,IX,IY)
        END DO  
     END DO   
  END DO
  
  !************************!
  ! FIRST CORRECTION-STEP *!
  !************************!

  DO K=1,NORD1-1
     
     DO IX=1,GXPT
        DO IY=1,GYPT

           !*****************************************************************************************!
           !* COMPUTE CORRECTED PSEUDO-VELOCITIES   
           DO J=1,NLEV-1
           YZDOT(K+1,J,IX,IY)    = HALF*ABS(YZDOT(K,J,IX,IY))*(                                      &
                                 & (ABS(YPSIG(J+1,IX,IY))-ABS(YPSIG(J,IX,IY)))                       &
                                 & /(ABS(YPSIG(J,IX,IY))+ABS(YPSIG(J+1,IX,IY))+REPS) )               &
                                 & -(PDT/TWO)*YZDOT(K,J,IX,IY)*( (HALF*( ( YZDOT(K,J,IX,IY)*(        &
                                 & ABS(YPSIG(J+1,IX,IY))+ABS(YPSIG(J,IX,IY)))                        &
                                 & -YZDOT(K,J-1,IX,IY)*(ABS(YPSIG(NBC(J-1),IX,IY))                   &
                                 & +ABS(YPSIG(J,IX,IY))) )/PDETA(J))/( (0.25d0*(                     &
                                 &  ABS(YPSIG(J+1,IX,IY))+ABS(YPSIG(J,IX,IY))                        &
                                 & +ABS(YPSIG(NBC(J-1),IX,IY))+ABS(YPSIG(J,IX,IY))))+REPS))          &
                                 & + (HALF*( ( YZDOT(K,J+1,IX,IY)*(                                  &
                                 & ABS(YPSIG(NBC(J+2),IX,IY))+ABS(YPSIG(J+1,IX,IY)))                 &
                                 & -YZDOT(K,J,IX,IY)*(ABS(YPSIG(J+1,IX,IY))                          &
                                 & +ABS(YPSIG(J,IX,IY))) )/PDETA(J+1))/( (0.25d0*(                   &
                                 &  ABS(YPSIG(NBC(J+2),IX,IY))+ABS(YPSIG(J+1,IX,IY))                 &
                                 & +ABS(YPSIG(J+1,IX,IY))+ABS(YPSIG(J,IX,IY))))+REPS))               &
                                 & +((PM1(J+1,IX,IY)-PM0(J+1,IX,IY)+PM1(J,IX,IY)-PM0(J,IX,IY))/PDT)) &
                                 & /(PM0(J+1,IX,IY)+PM0(J,IX,IY))
           END DO

           YZDOT(K+1,0,IX,IY)    = ZERO
           YZDOT(K+1,NLEV,IX,IY) = ZERO
        
           !********************************************************************************************!
           !* GLOBAL FCT NON-OSCILLATORY TREATMENT
           IF (LFCT_MONO) THEN   
              DO J=1,NLEV-1

                 ZIN_DOWN           = ( MAX( YPSI(K,J+1,IX,IY),YPSI(K,J,IX,IY),                  &
                                    & YPSI(K,NBC(J-1),IX,IY),YPSI(0,J+1,IX,IY),                  &
                                    & YPSI(0,J,IX,IY),YPSI(0,NBC(J-1),IX,IY))-YPSI(K,J,IX,IY) )  &
                                    & /( (ONE/PM0(J,IX,IY))*(PDT/PDETA(J))*(                     &
                                    & MAX(ZERO,YZDOT(K+1,J-1,IX,IY))*YPSI(K,NBC(J-1),IX,IY)      &
                                    & -MIN(ZERO,YZDOT(K+1,J,IX,IY))*YPSI(K,J+1,IX,IY))+REPS  )

                 ZIN_UP             = ( MAX( YPSI(K,NBC(J+2),IX,IY),YPSI(K,J+1,IX,IY),           &
                                    & YPSI(K,J,IX,IY),YPSI(0,NBC(J+2),IX,IY),                    &
                                    & YPSI(0,J,IX,IY),YPSI(0,J+1,IX,IY))-YPSI(K,J+1,IX,IY) )     &
                                    & /( (ONE/PM0(J+1,IX,IY))*(PDT/PDETA(J+1))*(                 &
                                    & MAX(ZERO,YZDOT(K+1,J,IX,IY))*YPSI(K,J,IX,IY)               &
                                    & -MIN(ZERO,YZDOT(K+1,J+1,IX,IY))*YPSI(K,NBC(J+2),IX,IY))    &
                                    & +REPS )
           
                 ZOUT_DOWN          = -( MIN( YPSI(K,J+1,IX,IY),YPSI(K,J,IX,IY),                 &
                                    & YPSI(K,NBC(J-1),IX,IY),YPSI(0,J+1,IX,IY),                  &
                                    & YPSI(0,J,I),YPSI(0,NBC(J-1),IX,IY))-YPSI(K,J,IX,IY) )      &
                                    & /( (ONE/PM0(J,IX,IY))*(PDT/PDETA(J))*(                     &
                                    & MAX(ZERO,YZDOT(K+1,J,IX,IY))*YPSI(K,J,IX,IY)               &
                                    & -MIN(ZERO,YZDOT(K+1,J-1,IX,IY))*YPSI(K,J,IX,IY))+REPS  )

                 ZOUT_UP            = -( MIN( YPSI(K,J+1,IX,IY),YPSI(K,J,IX,IY),                 &
                                    & YPSI(K,NBC(J+2),IX,IY),YPSI(0,J+1,IX,IY),                  &
                                    & YPSI(0,J,I),YPSI(0,NBC(J+2),IX,IY))-YPSI(K,J+1,IX,IY) )    &
                                    & /( (ONE/PM0(J+1,IX,IY))*(PDT/PDETA(J+1))*(                 &
                                    & MAX(ZERO,YZDOT(K+1,J+1,IX,IY))*YPSI(K,J+1,IX,IY)           &
                                    & -MIN(ZERO,YZDOT(K+1,J,IX,IY))*YPSI(K,J+1,IX,IY))+REPS  )
           
                 ZZDOT_MONO         = MIN(ONE,ZOUT_DOWN,ZIN_UP)*MAX(ZERO,YZDOT(K+1,J,IX,IY))     &
                                    & + MIN(ONE,ZIN_DOWN,ZOUT_UP)*MIN(ZERO,YZDOT(K+1,J,IX,IY))
           
                 YZDOT(K+1,J,IX,IY) = (ONE-RFCT)*YZDOT(K+1,J,IX,IY)+RFCT*ZZDOT_MONO
                 
              END DO
           END IF
        
           !*************************************************************************************!
           !* CORRECTED UPWIND SCHEME
           DO J=1,NLEV   
              YPSI(K+1,J,IX,IY)  = YPSI(K,J,IX,IY) - (ONE/PM0(J,IX,IY))*(PDT/PDETA(J))*(         &
                                 &  (MIN(ZERO,YZDOT(K+1,J,IX,IY))*YPSI(K,NBC(J+1),IX,IY)         &
                                 &  +MAX(ZERO,YZDOT(K+1,J,IX,IY))*YPSI(K,J,IX,IY))               &
                                 & -(MIN(ZERO,YZDOT(K+1,J-1,IX,IY))*YPSI(K,J,IX,IY)              &
                                 &  +MAX(ZERO,YZDOT(K+1,J-1,IX,IY))*YPSI(K,NBC(J-1),IX,IY)) )
              YPSIG(J,IX,IY)     = (PM0(J,IX,IY)/PM1(J,IX,IY))*YPSI(K+1,J,IX,IY)
           END DO
           
        END DO
     END DO
     
  END DO

  
  !*************************!
  ! SECOND CORRECTION-STEP *!
  !*************************!

  DO K=NORD1,NORD2-1
     
     DO IX=1,NXPT
        DO IY=1,GYPT
           
           !**********************************************************************************************!
           !* COMPUTE CORRECTED PSEUDO-VELOCITIES
           DO J=1,NLEV-1
              YZDOT(K+1,J,IX,IY)  = HALF*ABS(YZDOT(K,J,IX,IY))*(YPSIG(J+1,IX,IY)-YPSIG(J,IX,IY) )          &
                                  & -(PDT/TWO)*YZDOT(K,J,IX,IY)*(                                          & 
                                  & (HALF*((YZDOT(K,J,IX,IY)*(YPSIG(J+1,IX,IY)+YPSIG(J,IX,IY))             &
                                  & -YZDOT(K,J-1,IX,IY)*(YPSIG(NBC(J-1),IX,IY)+YPSIG(J,IX,IY)))/PDETA(J))) &
                                  & +(HALF*((YZDOT(K,J+1,IX,IY)*(YPSIG(NBC(J+2),I)+YPSIG(J+1,IX,IY))       &
                                  & -YZDOT(K,J,IX,IY)*(YPSIG(J+1,IX,IY)+YPSIG(J,IX,IY)))/PDETA(J)))        &
                                  & +HALF*(YPSIG(J+1,IX,IY)+YPSIG(J,IX,IY))                                &
                                  & *((PM1(J+1,IX,IY)-PM0(J+1,IX,IY)+PM1(J,IX,IY)-PM0(J,IX,IY))/PDT) )     &
                                  & /(PM0(J+1,IX,IY)+PM0(J,IX,IY))
           END DO
        
           YZDOT(K+1,0,IX,IY)     = ZERO
           YZDOT(K+1,NLEV,IX,IY)  = ZERO

           
           !********************************************************************************************!
           !* GLOBAL FCT NON-OSCILLATORY TREATMENT
           IF (LFCT_MONO) THEN
              DO J=1,NLEV-1

                 ZIN_DOWN           = ( MAX( YPSI(K,J+1,IX,IY),YPSI(K,J,IX,IY),                  &
                                    & YPSI(K,NBC(J-1),IX,IY),YPSI(0,J+1,IX,IY),                  &
                                    & YPSI(0,J,IX,IY),YPSI(0,NBC(J-1),IX,IY))-YPSI(K,J,IX,IY) )  &
                                    & /( (ONE/PM0(J,IX,IY))*(PDT/PDETA(J))*(                     &
                                    & MAX(ZERO,YZDOT(K+1,J-1,IX,IY))*YPSI(K,NBC(J-1),IX,IY)      &
                                    & -MIN(ZERO,YZDOT(K+1,J,IX,IY))*YPSI(K,J+1,IX,IY))+REPS  )
 
                 ZIN_UP             = ( MAX( YPSI(K,NBC(J+2),IX,IY),YPSI(K,J+1,IX,IY),           &
                                    & YPSI(K,J,IX,IY),YPSI(0,NBC(J+2),IX,IY),                    &
                                    & YPSI(0,J,IX,IY),YPSI(0,J+1,IX,IY))-YPSI(K,J+1,IX,IY) )     &
                                    & /( (ONE/PM0(J+1,IX,IY))*(PDT/PDETA(J+1))*(                 &
                                    & MAX(ZERO,YZDOT(K+1,J,IX,IY))*YPSI(K,J,IX,IY)               &
                                    & -MIN(ZERO,YZDOT(K+1,J+1,IX,IY))*YPSI(K,NBC(J+2),IX,IY))    &
                                    & +REPS )
           
                 ZOUT_DOWN          = -( MIN( YPSI(K,J+1,IX,IY),YPSI(K,J,IX,IY),                 &
                                    & YPSI(K,NBC(J-1),IX,IY),YPSI(0,J+1,IX,IY),                  &
                                    & YPSI(0,J,I),YPSI(0,NBC(J-1),IX,IY))-YPSI(K,J,IX,IY) )      &
                                    & /( (ONE/PM0(J,IX,IY))*(PDT/PDETA(J))*(                     &
                                    & MAX(ZERO,YZDOT(K+1,J,IX,IY))*YPSI(K,J,IX,IY)               &
                                    & -MIN(ZERO,YZDOT(K+1,J-1,IX,IY))*YPSI(K,J,IX,IY))+REPS  )

                 ZOUT_UP            = -( MIN( YPSI(K,J+1,IX,IY),YPSI(K,J,IX,IY),                 &
                                    & YPSI(K,NBC(J+2),IX,IY),YPSI(0,J+1,IX,IY),                  &
                                    & YPSI(0,J,I),YPSI(0,NBC(J+2),IX,IY))-YPSI(K,J+1,IX,IY) )    &
                                    & /( (ONE/PM0(J+1,IX,IY))*(PDT/PDETA(J+1))*(                 &
                                    & MAX(ZERO,YZDOT(K+1,J+1,IX,IY))*YPSI(K,J+1,IX,IY)           &
                                    & -MIN(ZERO,YZDOT(K+1,J,IX,IY))*YPSI(K,J+1,IX,IY))+REPS  )
           
                 ZZDOT_MONO         = MIN(ONE,ZOUT_DOWN,ZIN_UP)*MAX(ZERO,YZDOT(K+1,J,IX,IY))     &
                                    & + MIN(ONE,ZIN_DOWN,ZOUT_UP)*MIN(ZERO,YZDOT(K+1,J,IX,IY))
           
                 YZDOT(K+1,J,IX,IY) = (ONE-RFCT)*YZDOT(K+1,J,IX,IY)+RFCT*ZZDOT_MONO
           
              END DO
           END IF

           !***********************************************************************************!
           !* CORRECTED UPWIND SCHEME
           DO J=1,NLEV
              YPSI(K+1,J,IX,IY)     = YPSI(K,J,IX,IY) - (ONE/PM0(J,IX,IY))*(PDT/PDETA(J))*(     &
                                    & YZDOT(K+1,J,IX,IY)-YZDOT(K+1,J-1,IX,IY) )
              YPSIG(J,IX,IY)        = (PM0(J,IX,IY)/PM1(J,IX,IY))*YPSI(K+1,J,IX,IY)
           END DO
           
        END DO
     END DO
     
  END DO
  !*****************************************************************************************!
  !*****************************************************************************************!

  DO IX=1,GXPT
     DO IY=1,GYPT

        PPSI_DEP(:,IX,IY)         = YPSI(NORD2,:,IX,IY)+ZPSI_MIN

        IF (LLMASS) THEN
           PZDOT(KS,0,IX,IY)      = ZERO
           PZDOT(KS+1,0,IX,IY)    = ZERO
           DO J=1,NLEV-1
              PZDOT(KS,J,IX,IY)   = HALF*YZDOT(NORD2,J,IX,IY)*(                    &
                                  & HALF*(PPSI_ARR(J+1,IX,IY)+PPSI_DEP(J+1,IX,IY)) &
                                  & + HALF*(PPSI_ARR(J,IX,IY)+PPSI_DEP(J,IX,IY))   )
              PZDOT(KS+1,J,IX,IY) = YZDOT(NORD2,J,IX,IY)
           END DO
           PZDOT(KS,NLEV,IX,IY)   = ZERO
           PZDOT(KS+1,NLEV,IX,IY) = ZERO
        END IF
        
     END DO   
  END DO
    
  
END SUBROUTINE ADV_ZF
!################################################################################################!
!************************************************************************************************!
!************************************************************************************************!
!************************************************************************************************!
!************************************************************************************************!
!################################################################################################!
SUBROUTINE ADV_ZH(PPSI_DEP,PPSI_ARR,PZDOT,PM1,PM0,PDETA,PDT,KS,LDMASS)
  
  IMPLICIT NONE

  REAL(8),                       INTENT(OUT)    :: PPSI_DEP(0:NLEV,GXPT,GYPT)
  REAL(8),                       INTENT(IN)     :: PPSI_ARR(0:NLEV,GXPT,GYPT)
  REAL(8),                       INTENT(INOUT)  :: PZDOT(3,NLEV,GXPT,GYPT)
  REAL(8),                       INTENT(IN)     :: PM1(0:NLEV,GXPT,GYPT)
  REAL(8),                       INTENT(IN)     :: PM0(0:NLEV,GXPT,GYPT)
  REAL(8),                       INTENT(IN)     :: PDETA(0:NLEV)
  REAL(8),                       INTENT(IN)     :: PDT
  INTEGER,                       INTENT(IN)     :: KS
  LOGICAL,OPTIONAL,              INTENT(IN)     :: LDMASS

  REAL(8)                                       :: YPSI(0:NORD2,NLEV+1,GXPT,GYPT)
  REAL(8)                                       :: YZDOT(NORD2,0:NLEV+1,GXPT,GYPT)
  REAL(8)                                       :: YPSIG(NLEV+1,GXPT,GYPT)
  REAL(8)                                       :: ZM0(NLEV+1,GXPT,GYPT)
  REAL(8)                                       :: ZM1(NLEV+1,GXPT,GYPT)
  REAL(8)                                       :: ZDETA(NLEV+1)
  REAL(8)                                       :: ZPSI_MIN,ZZDOT_MONO
  REAL(8)                                       :: ZIN_DOWN,ZIN_UP,ZOUT_DOWN,ZOUT_UP
  INTEGER(8)                                    :: IX,IY,J,K
  LOGICAL                                       :: LLMASS

  IF (PRESENT(LDMASS)) THEN
     LLMASS = LDMASS
  ELSE
     LLMASS = .FALSE.
  END IF
  
  ZPSI_MIN = MINVAL(PPSI_ARR)

  !IF (KTOP==0)
  DO J=1,NLEV+1
     ZDETA(J)       = PDETA(J-1)
     ZM0(J,:,:)     = PM0(J-1,:,:)
     ZM1(J,:,:)     = PM1(J-1,:,:)
     YPSI(0,J,:,:)  = PPSI_ARR(J-1,:,:)-ZPSI_MIN
  END DO
     
  DO J=1,NLEV
     YZDOT(1,J,:,:) = PZDOT(KS,J,:,:)
  END DO   

  YZDOT(:,0,:,:)      = ZERO
  YZDOT(:,NLEV+1,:,:) = ZERO

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !***************************************!
  ! DONOR-CELL UPWIND FIRST-ORDER SCHEME *!
  !***************************************!

  DO I=1,GXPT
     DO IY=1,GYPT
        DO J=1,NLEV+1
           YPSI(1,J,IX,IY) = YPSI(0,J,IX,IY) - (ONE/ZM0(J,IX,IY))*(PDT/ZDETA(J))*(          &
                           &  (MIN(ZERO,YZDOT(1,J,IX,IY))*YPSI(0,NBL(J+1),IX,IY)            &
                           & +MAX(ZERO,YZDOT(1,J,IX,IY))*YPSI(0,J,IX,IY))                   &
                           & -(MIN(ZERO,YZDOT(1,J-1,IX,IY))*YPSI(0,J,IX,IY)                 &
                           & +MAX(ZERO,YZDOT(1,J-1,IX,IY))*YPSI(0,NBL(J-1),IX,IY)) )
           YPSIG(J,IX,IY)  = (ZM0(J,IX,IY)/ZM1(J,IX,IY))*YPSI(1,J,IX,IY)
     END DO
  END DO
  
  !************************!
  ! FIRST CORRECTION-STEP *!
  !************************!

  DO K=1,NORD1-1
     
     DO IX=1,GXPT
        DO IY=1,GYPT
        
           !********************************************************************************************!
           !* COMPUTE SURFACE PSEUDO-VELOCITIES   
           DO J=1,NLEV
              YZDOT(K+1,J,IX,IY)    = HALF*ABS(YZDOT(K,J,IX,IY))*(                                     &
                                    & (ABS(YPSIG(NBL(J+1),IX,IY))-ABS(YPSIG(J,IX,IY)))                 &
                                    & /(ABS(YPSIG(J,IX,IY))+ABS(YPSIG(NBL(J+1),IX,IY))+REPS) )         &
                                    & -(PDT/TWO)*YZDOT(K,J,IX,IY)*( (HALF*( ( YZDOT(K,J,IX,IY)*(       &
                                    & ABS(YPSIG(NBL(J+1),IX,IY))+ABS(YPSIG(J,IX,IY)))                  &
                                    & -YZDOT(K,J-1,IX,IY)*(ABS(YPSIG(NBL(J-1),IX,IY))                  &
                                    & +ABS(YPSIG(J,IX,IY))) )/ZDETA(J))/( (0.25d0*(                    &
                                    &  ABS(YPSIG(NBL(J+1),IX,IY))+ABS(YPSIG(J,IX,IY))                  &
                                    & +ABS(YPSIG(NBL(J-1),IX,IY))+ABS(YPSIG(J,IX,IY))))+REPS))         &
                                    & + (HALF*( ( YZDOT(K,J+1,IX,IY)*(                                 &
                                    & ABS(YPSIG(NBL(J+2),IX,IY))+ABS(YPSIG(NBL(J+1),IX,IY)))           &
                                    & -YZDOT(K,J,I)*(ABS(YPSIG(NBL(J+1),IX,IY))                        &
                                    & +ABS(YPSIG(J,IX,IY))))/ZDETA(J+1))/( (0.25d0*(                   &
                                    &  ABS(YPSIG(NBL(J+2),IX,IY))+ABS(YPSIG(NBL(J+1),IX,IY))           &
                                    & +ABS(YPSIG(NBL(J+1),IX,IY))+ABS(YPSIG(J,IX,IY))))+REPS))         &
                                    & +((ZM1(NBL(J+1),IX,IY)-ZM0(NBL(J+1),IX,IY)                       &
                                    & +ZM1(J,IX,IY)-ZM0(J,IX,IY))/PDT))                                &
                                    & /(ZM0(NBL(J+1),IX,IY)+ZM0(J,IX,IY))
           END DO

           YZDOT(K+1,0,IX,IY)       = ZERO
           YZDOT(K+1,NLEV+1,IX,IY)  = ZERO
        
           !*******************************************************************************************!
           !* GLOBAL FCT NON-OSCILLATORY OPTION
           IF (LFCT_MONO) THEN
              DO J=1,NLEV
           
                 ZIN_DOWN           = ( MAX( YPSI(K,NBL(J+1),IX,IY),YPSI(K,J,IX,IY),                   &  
                                    & YPSI(K,NBL(J-1),IX,IY),YPSI(0,NBL(J+1),IX,IY),                   & 
                                    & YPSI(0,J,IX,IY),YPSI(0,NBL(J-1),IX,IY))-YPSI(K,J,IX,IY) )        &
                                    & /( (ONE/ZM0(J,IX,IY))*(PDT/ZDETA(J))*(                           &
                                    & MAX(ZERO,YZDOT(K+1,J-1,IX,IY))*YPSI(K,NBL(J-1),IX,IY)            &
                                    & -MIN(ZERO,YZDOT(K+1,J,IX,IY))*YPSI(K,NBL(J+1),IX,IY))+REPS )

                 ZIN_UP             = ( MAX( YPSI(K,NBL(J+2),IX,IY),YPSI(K,J,IX,IY),                   &  
                                    & YPSI(K,NBL(J+1),IX,IY),YPSI(0,NBL(J+2),IX,IY),                   & 
                                    & YPSI(0,J,IX,IY),YPSI(0,NBL(J+1),IX,IY))-YPSI(K,NBL(J+1),IX,IY) ) &
                                    & /( (ONE/ZM0(J+1,IX,IY))*(PDT/ZDETA(J+1))*(                       &
                                    & MAX(ZERO,YZDOT(K+1,J,IX,IY))*YPSI(K,J,IX,IY)                     &
                                    & -MIN(ZERO,YZDOT(K+1,NBL(J+1),IX,IY))*YPSI(K,NBL(J+2),IX,IY))     &
                                    & + REPS )
              
                 ZOUT_DOWN          = -( MIN( YPSI(K,NBL(J+1),IX,IY),YPSI(K,J,IX,IY),                  &
                                    & YPSI(K,NBL(J-1),IX,IY),YPSI(0,NBL(J+1),IX,IY),                   & 
                                    & YPSI(0,J,IX,IY),YPSI(0,NBL(J-1),IX,IY))-YPSI(K,J,IX,IY) )        &
                                    & /( (ONE/ZM0(J,IX,IY))*(PDT/ZDETA(J))*(                           &
                                    & MAX(ZERO,YZDOT(K+1,J,IX,IY))*YPSI(K,J,IX,IY)                     &
                                    & -MIN(ZERO,YZDOT(K+1,J-1,IX,IY))*YPSI(K,J,IX,IY))+REPS )

                 ZOUT_UP            = -( MIN( YPSI(K,NBL(J+2),IX,IY),YPSI(K,J,IX,IY),                  &
                                    & YPSI(K,NBL(J+1),IX,IY),YPSI(0,NBL(J+2),IX,IY),                   & 
                                    & YPSI(0,J,IX,IY),YPSI(0,NBL(J+1),IX,IY))-YPSI(K,NBL(J+1),IX,IY) ) &
                                    & /( (ONE/ZM0(J+1,IX,IY))*(PDT/ZDETA(J+1))*(                       &
                                    & MAX(ZERO,YZDOT(K+1,NBL(J+1),IX,IY))*YPSI(K,NBL(J+1),IX,IY)       &
                                    & -MIN(ZERO,YZDOT(K+1,J,IX,IY))*YPSI(K,NBL(J+1),IX,IY))+REPS )
              
                 ZZDOT_MONO         = MIN(ONE,ZOUT_DOWN,ZIN_UP)*MAX(ZERO,YZDOT(K+1,J,IX,IY))           &
                                    & +MIN(ONE,ZIN_DOWN,ZOUT_UP)*MIN(ZERO,YZDOT(K+1,J,IX,IY))
              
                 YZDOT(K+1,J,IX,IY) = (ONE-RFCT)*YZDOT(K+1,J,IX,IY)+RFCT*ZZDOT_MONO
              
              END DO
           END IF

           !*******************************************************************************************!
           !* CORRECTED UPWIND SCHEME
           DO J=1,NLEV+1   
              YPSI(K+1,J,IX,IY)     = YPSI(K,J,IX,IY) - (ONE/ZM0(J,IX,IY))*(PDT/ZDETA(J))*(            &
                                    & (MIN(ZERO,YZDOT(K+1,J,IX,IY))*YPSI(K,NBL(J+1),IX,IY)             &
                                    & +MAX(ZERO,YZDOT(K+1,J,IX,IY))*YPSI(K,J,IX,IY))                   &
                                    & -(MIN(ZERO,YZDOT(K+1,J-1,IX,IY))*YPSI(K,J,IX,IY)                 &
                                    & +MAX(ZERO,YZDOT(K+1,J-1,IX,IY))*YPSI(K,NBL(J-1),IX,IY)) )
              YPSIG(J,IX,IY)        = (ZM0(J,IX,IY)/ZM1(J,IX,IY))*YPSI(K+1,J,IX,IY)
        END DO
        
     END DO
  END DO

  !*************************!
  ! SECOND CORRECTION-STEP *!
  !*************************!

  DO K=NORD1,NORD2-1
     
     DO IX=1,GXPT
        DO IY=1,GYPT

           !**********************************************************************************************!
           !* COMPUTE SURFACE PSEUDO-VELOCITIES
           DO J=1,NLEV
              YZDOT(K+1,J,IX,IY)    = HALF*ABS(YZDOT(K,J,IX,IY))*(YPSIG(NBL(J+1),IX,IY)-YPSIG(J,IX,IY))   & 
                                    & -(PDT/TWO)*YZDOT(K,J,IX,IY)*(                                       &
                                    & (HALF*((YZDOT(K,J,IX,IY)*(YPSIG(NBL(J+1),IX,IY)+YPSIG(J,IX,IY))     &
                                    & -YZDOT(K,J-1,IX,IY)*(                                               &
                                    & YPSIG(NBL(J-1),IX,IY)+YPSIG(J,IX,IY)) )/ZDETA(J)))                  &
                                    & +(HALF*((YZDOT(K,J+1,IX,IY)*(YPSIG(NBL(J+2),IX,IY)                  &
                                    & +YPSIG(NBL(J+1),IX,IY))                                             &
                                    & -YZDOT(K,J,IX,IY)*(YPSIG(NBL(J+1),IX,IY)+YPSIG(J,IX,IY)))           &
                                    & /ZDETA(J+1)))                                                       &
                                    & +HALF*(YPSIG(NBL(J+1),IX,IY)+YPSIG(J,IX,IY))                        &
                                    & *((ZM1(NBL(J+1),IX,IY)-ZM0(NBL(J+1),IX,IY)                          &
                                    & +ZM1(J,IX,IY)-ZM0(J,IX,IY))/PDT))/(ZM0(NBL(J+1),IX,IY)+ZM0(J,IX,IY))
           END DO
        
           YZDOT(K+1,0,I)           = ZERO
           YZDOT(K+1,NLEV+1,I)      = ZERO

           !**********************************************************************************************!
           !* GLOBAL FCT NON-OSCILLATORY OPTION
           IF (LFCT_MONO) THEN
              DO J=1,NLEV
           
                 ZIN_DOWN           = ( MAX( YPSI(K,NBL(J+1),IX,IY),YPSI(K,J,IX,IY),                      &  
                                    & YPSI(K,NBL(J-1),IX,IY),YPSI(0,NBL(J+1),IX,IY),                      & 
                                    & YPSI(0,J,IX,IY),YPSI(0,NBL(J-1),IX,IY))-YPSI(K,J,IX,IY) )           &
                                    & /( (ONE/ZM0(J,IX,IY))*(PDT/ZDETA(J))*(                              &
                                    & MAX(ZERO,YZDOT(K+1,J-1,IX,IY))*YPSI(K,NBL(J-1),IX,IY)               &
                                    & -MIN(ZERO,YZDOT(K+1,J,IX,IY))*YPSI(K,NBL(J+1),IX,IY))+REPS )

                 ZIN_UP             = ( MAX( YPSI(K,NBL(J+2),IX,IY),YPSI(K,J,IX,IY),                      &  
                                    & YPSI(K,NBL(J+1),IX,IY),YPSI(0,NBL(J+2),IX,IY),                      & 
                                    & YPSI(0,J,IX,IY),YPSI(0,NBL(J+1),IX,IY))-YPSI(K,NBL(J+1),IX,IY) )    &
                                    & /( (ONE/ZM0(NBL(J+1),IX,IY))*(PDT/ZDETA(NBL(J+1)))*(                &
                                    & MAX(ZERO,YZDOT(K+1,J,IX,IY))*YPSI(K,J,IX,IY)                        &
                                    & -MIN(ZERO,YZDOT(K+1,NBL(J+1),IX,IY))*YPSI(K,NBL(J+2),IX,IY))        &
                                    & + REPS )
              
                 ZOUT_DOWN          = -( MIN( YPSI(K,NBL(J+1),IX,IY),YPSI(K,J,IX,IY),                     &
                                    & YPSI(K,NBL(J-1),IX,IY),YPSI(0,NBL(J+1),IX,IY),                      &  
                                    & YPSI(0,J,IX,IY),YPSI(0,NBL(J-1),IX,IY))-YPSI(K,J,IX,IY) )           &
                                    & /( (ONE/ZM0(J,IX,IY))*(PDT/ZDETA(J))*(                              &
                                    & MAX(ZERO,YZDOT(K+1,J,IX,IY))*YPSI(K,J,IX,IY)                        &
                                    & -MIN(ZERO,YZDOT(K+1,J-1,IX,IY))*YPSI(K,J,IX,IY))+REPS )

                 ZOUT_UP            = -( MIN( YPSI(K,NBL(J+2),IX,IY),YPSI(K,J,IX,IY),                     &
                                    & YPSI(K,NBL(J+1),IX,IY),YPSI(0,NBL(J+2),IX,IY),                      & 
                                    & YPSI(0,J,IX,IY),YPSI(0,NBL(J+1),IX,IY))-YPSI(K,NBL(J+1),IX,IY) )    &
                                    & /( (ONE/ZM0(NBL(J+1),IX,IY))*(PDT/ZDETA(NBL(J+1)))*(                &
                                    & MAX(ZERO,YZDOT(K+1,NBL(J+1),IX,IY))*YPSI(K,NBL(J+1),IX,IY)          &
                                    & -MIN(ZERO,YZDOT(K+1,J,IX,IY))*YPSI(K,NBL(J+1),IX,IY))+REPS )
              
                 ZZDOT_MONO         = MIN(ONE,ZOUT_DOWN,ZIN_UP)*MAX(ZERO,YZDOT(K+1,J,IX,IY))              &
                                    & +MIN(ONE,ZIN_DOWN,ZOUT_UP)*MIN(ZERO,YZDOT(K+1,J,IX,IY))
              
                 YZDOT(K+1,J,IX,IY) = (ONE-RFCT)*YZDOT(K+1,J,IX,IY)+RFCT*ZZDOT_MONO
              
              END DO
           END IF
           
           !**********************************************************************************************!
           !* CORRECTED UPWIND SCHEME
           DO J=1,NLEV+1
              YPSI(K+1,J,IX,IY)   = YPSI(K,J,IX,IY) - (ONE/ZM0(J,IX,IY))*(PDT/ZDETA(J))*(                 &
                                  &    YZDOT(K+1,J,IX,IY)-YZDOT(K+1,J-1,IX,IY) )
              YPSIG(J,IX,IY)      = (ZM0(J,IX,IY)/ZM1(J,IX,IY))*YPSI(K+1,J,IX,IY)
           END DO
           
        END DO
     END DO
  END DO

  !*************************************************************************************!
  DO IX=1,GXPT
     DO IY=1,GYPT
        
        DO J=0,NLEV
           PPSI_DEP(J,IX,IY)      = YPSI(NORD2,J+1,IX,IY)+ZPSI_MIN
        END DO

        IF (LLMASS) THEN
           DO J=1,NLEV
              PZDOT(KS,J,IX,IY)   = HALF*YZDOT(NORD2,J,IX,IY)*(                         &
                                  & HALF*(PPSI_ARR(NBL(J+1),IX,IY)+PPSI_DEP(J,IX,IY))   &
                                  & + HALF*(PPSI_ARR(J,IX,IY)+PPSI_DEP(NBL(J+1),IX,IY)) )
              PZDOT(KS+1,J,IX,IY) = YZDOT(NORD2,J,IX,IY)
           END DO
        END IF

     END DO
   END DO
     
END SUBROUTINE ADV_ZH
!################################################################################################!
!************************************************************************************************!
!************************************************************************************************!
!************************************************************************************************!
!************************************************************************************************!
!################################################################################################!
SUBROUTINE ADV_XY(PPSI_DEP,PPSI_ARR,PXDOT,PYDOT,PM1,PM0,PDT,KTOP,LDMASS)

  IMPLICIT NONE

  REAL(8),                       INTENT(OUT)    :: PPSI_DEP(KTOP:NLEV,GXPT,GYPT)
  REAL(8),                       INTENT(IN)     :: PPSI_ARR(KTOP:NLEV,GXPT,GYPT)
  REAL(8),                       INTENT(INOUT)  :: PXDOT(KTOP:NLEV,GXPT,GYPT)
  REAL(8),                       INTENT(INOUT)  :: PYDOT(KTOP:NLEV,GXPT,GYPT)
  REAL(8),                       INTENT(IN)     :: PM1(KTOP:NLEV,GXPT,GYPT)
  REAL(8),                       INTENT(IN)     :: PM0(KTOP:NLEV,GXPT,GYPT)
  REAL(8),                       INTENT(IN)     :: PDT
  INTEGER(8),                    INTENT(IN)     :: KTOP
  LOGICAL,OPTIONAL,              INTENT(IN)     :: LDMASS
  
  
  REAL(8)                                       :: YPSI(0:NORD2,KTOP:NLEV,GXPT,GYPT)
  REAL(8)                                       :: YXDOT(NORD2,KTOP:NLEV,GXPT,GYPT)
  REAL(8)                                       :: YYDOT(NORD2,KTOP:NLEV,GXPT,GYPT)
  REAL(8)                                       :: YPSIG(KTOP:NLEV,GXPT,GYPT)
  REAL(8)                                       :: ZINX_UP,ZOUTX_UP,ZINX_DOWN,ZOUTX_DOWN
  REAL(8)                                       :: ZINY_UP,ZOUTY_UP,ZINY_DOWN,ZOUTY_DOWN
  REAL(8)                                       :: ZPSI_MIN,ZXDOT_MONO,ZYDOT_MONO
  INTEGER(8)                                    :: I,J,K,L
  LOGICAL                                       :: LLMASS

  
  IF (PRESENT(LDMASS)) THEN
     LLMASS = LDMASS
  ELSE
     LLMASS = .FALSE.
  END IF   

  ZPSI_MIN = MINVAL(PPSI_ARR)

  DO I=1,GXPT
     DO J=1,GYPT
        YPSI(0,:,I,J)  = PPSI_ARR(:,I,J) - ZPSI_MIN
        YXDOT(1,:,I,J) = PXDOT(:,I,J)
        YYDOT(1,:,I,J) = PYDOT(:,I,J)
     END DO   
  END DO


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !***************************************!
  ! DONOR-CELL UPWIND FIRST-ORDER SCHEME *!
  !***************************************!

  DO I=1,GXPT
     DO J=1,GYPT

        DO J=KTOP,NLEV
           YPSI(1,J,I)    = YPSI(0,L,I,J)                                                            &
                          & -(ONE/PM0(L,I,J))*(PDT/RDX)*(                                            &
                          & (MIN(ZERO,YXDOT(1,L,I,J))*YPSI(0,L,NPERIOX(I+1),J)                       &
                          & +MAX(ZERO,YXDOT(1,L,I,J))*YPSI(0,L,I,J))                                 &
                          & -(MIN(ZERO,YXDOT(1,L,NPERIOX(I-1),J))*YPSI(0,L,I,J)                      & 
                          & +MAX(ZERO,YXDOT(1,L,NPERIOX(I-1),J))*YPSI(0,L,NPERIOX(I-1),J)) )         &
                          & -(ONE/PM0(L,I,J))*(PDT/RDY)*(                                            &
                          & (MIN(ZERO,YYDOT(1,L,I,J))*YPSI(0,L,I,NPERIOY(J+1))                       &
                          & +MAX(ZERO,YYDOT(1,L,I,J))*YPSI(0,L,I,J))                                 &
                          & -(MIN(ZERO,YYDOT(1,L,I,NPERIOY(J-1)))*YPSI(0,L,I,J)                      & 
                          & +MAX(ZERO,YYDOT(1,L,I,NPERIOY(J-1)))*YPSI(0,L,I,NPERIOY(J-1))) )
           
           YPSIG(L,I,J)   = (PM0(L,I,J)/PM1(L,I,J))*YPSI(1,L,I,J)
        END DO
        
     END DO   
  END DO
  
  !************************!
  ! FIRST CORRECTION-STEP *!
  !************************!

  DO K=1,NORD1-1   

     !***********************************************************************************************!
     !* COMPUTE SURFACE PSEUDO-VELOCITIES *!
     !*************************************!
     DO I=1,GXPT
        DO J=1,GYPT
           DO L=KTOP,NLEV   

              YXDOT(K+1,L,I,J) = HALF*ABS(YXDOT(K,L,I,J))*(                                          &
                          & (ABS(YPSIG(L,NPERIO(I+1),J))-ABS(YPSIG(L,I,J)))                          &  
                          & /(ABS(YPSIG(L,I,J))+ABS(YPSIG(L,NPERIO(I+1),J))+REPS) )                  &
                          & -(PDT/TWO)*YXDOT(K,L,I,J)*( (                                            &
                          & ( HALF*((YXDOT(K,L,I,J)*(ABS(YPSIG(L,NPERIO(I+1),J))+ABS(YPSIG(L,I,J)))  & ! div(i,j)
                          & -YXDOT(K,L,NPERIO(I-1),J)*(ABS(YPSIG(L,NPERIO(I-1),J))                   &
                          & +ABS(YPSIG(L,I,J))))/RDX)                                                &
                          & + HALF*((YYDOT(K,L,I,J)*(ABS(YPSIG(L,I,NPERIO(J+1)))+ABS(YPSIG(L,I,J)))  & 
                          & -YYDOT(K,L,I,NPERIO(J-1))*(ABS(YPSIG(L,I,NPERIO(J-1)))                   &
                          & +ABS(YPSIG(L,I,J))))/RDY) )                                              &
                          & /((0.125d0*(ABS(YPSIG(L,NPERIO(I+1),J))+4.d0*ABS(YPSIG(L,I,J))           &
                          & +ABS(YPSIG(L,NPERIO(I-1),J))+ABS(YPSIG(L,I,NPERIO(J+1)))                 &
                          & +ABS(YPSIG(L,I,NPERIO(J-1)))))+REPS) )                                   &
                          & +( (HALF*((YXDOT(K,L,NPERIO(I+1),J)*(                                    & ! div(i+1,j)
                          & ABS(YPSIG(L,NPERIO(I+2),J))+ABS(YPSIG(L,NPERIO(I+1),J)))                 &
                          & -YXDOT(K,L,I,J)*(ABS(YPSIG(L,NPERIO(I+1),J))+ABS(YPSIG(L,I,J))))/RDX)    &
                          & +HALF*((YYDOT(K,L,NPERIO(I+1),J)*(                                       &  
                          & ABS(YPSIG(L,NPERIO(I+1),NPERIO(J+1)))+ABS(YPSIG(L,NPERIO(I+1),J)))       & 
                          & -YYDOT(K,L,NPERIO(I+1),NPERIO(J-1))*(ABS(YPSIG(L,NPERIO(I+1),J))         &
                          & +ABS(YPSIG(L,NPERIO(I+1),NPERIO(J-1)))))/RDY) )                          & 
                          & /((0.125d0*(ABS(YPSIG(L,NPERIO(I+2),J))+ABS(YPSIG(L,I,J))                & 
                          & +4.d0*ABS(YPSIG(L,NPERIO(I+1),J))+ABS(YPSIG(L,NPERIO(I+1),NPERIO(J+1)))  & 
                          & +ABS(YPSIG(L,NPERIO(I+1),NPERIO(J-1)))))+REPS) )                         & 
                          & +((PM1(L,NPERIO(I+1),J)-PM0(L,NPERIO(I+1),J)                             & ! <DG/Dt>(i+0.5,j)
                          & + PM1(L,I,J)-PM0(L,I,J) )/PDT) )/(PM0(L,NPERIO(I+1),J)+PM0(L,I,J))

              YYDOT(K+1,L,I,J) = HALF*ABS(YYDOT(K,L,I,J))*(                                          &
                          & (ABS(YPSIG(L,I,NPERIO(J+1)))-ABS(YPSIG(L,I,J)))                          & ! 
                          & /(ABS(YPSIG(L,I,J))+ABS(YPSIG(L,I,NPERIO(J+1)))+REPS) )                  &
                          & -(PDT/TWO)*YYDOT(K,L,I,J)*( (                                            &
                          & ( HALF*((YXDOT(K,L,I,J)*(ABS(YPSIG(L,NPERIO(I+1),J))+ABS(YPSIG(L,I,J)))  & ! div(i,j)
                          & -YXDOT(K,L,NPERIO(I-1),J)*(ABS(YPSIG(L,NPERIO(I-1),J))                   &
                          & +ABS(YPSIG(L,I,J))))/RDX)                                                &
                          & + HALF*((YYDOT(K,L,I,J)*(ABS(YPSIG(L,I,NPERIO(J+1)))+ABS(YPSIG(L,I,J)))  & 
                          & -YYDOT(K,L,I,NPERIO(J-1))*(ABS(YPSIG(L,I,NPERIO(J-1)))                   &
                          & +ABS(YPSIG(L,I,J))))/RDY) )                                              &
                          & /((0.125d0*(ABS(YPSIG(L,NPERIO(I+1),J))+4.d0*ABS(YPSIG(L,I,J))           &
                          & +ABS(YPSIG(L,NPERIO(I-1),J))+ABS(YPSIG(L,I,NPERIO(J+1)))                 &
                          & +ABS(YPSIG(L,I,NPERIO(J-1)))))+REPS) )                                   &
                          & +( ( HALF*((YXDOT(K,L,I,NPERIO(J+1))*(                                   & ! div(i,j+1)
                          & ABS(YPSIG(L,NPERIO(I+1),NPERIO(J+1)))+ABS(YPSIG(L,I,NPERIO(J+1))))       &
                          & -YXDOT(K,L,,NPERIO(I-1),NPERIO(J+1))*(ABS(YPSIG(L,I,NPERIO(J+1)))        &
                          & +ABS(YPSIG(L,NPERIO(I-1),NPERIO(J+1)))))/RDX)                            &
                          & +HALF*((YYDOT(K,L,I,NPERIO(J+1))*(                                       &  
                          & ABS(YPSIG(L,I,NPERIO(J+2)))+ABS(YPSIG(L,I,NPERIO(J+1))))                 & 
                          & -YYDOT(K,L,I,J)*(ABS(YPSIG(L,I,NPERIO(J+1)))+ABS(YPSIG(L,I,J))))/RDY) )  & 
                          & /((0.125d0*(ABS(YPSIG(L,NPERIO(I+1),NPERIO(J+1)))+ABS(YPSIG(L,I,J))      &
                          & +ABS(YPSIG(L,NPERIO(I-1),NPERIO(J+1)))+4.d0*ABS(YPSIG(L,I,NPERIO(J+1)))  &
                          & +ABS(YPSIG(L,I,NPERIO(J+2)))))+REPS) )                                   & 
                          & +((PM1(L,I,NPERIO(J+1))-PM0(L,I,NPERIO(J+1))                             & ! <DG/Dt>(i,j+0.5)
                          & + PM1(L,I,J)-PM0(L,I,J) )/PDT) )/(PM0(L,I,NPERIO(J+1))+PM0(L,I,J))
           END DO
        END DO 
     END DO

     !********************************************************************************************************!
     !* GLOBAL FCT NON-OSCILLATORY TREATMENT *!
     !****************************************!
     IF (LFCT_MONO) THEN
        DO I=1,GXPT
           DO J=1,GYPT
              DO L=KTOP,NLEV
                 
                 !***************!
                 !* X-DIRECTION *!
                 !***************!

                 ZINX_DOWN        = ( MAX(YPSI(K,L,NPERIO(I+1),J),YPSI(K,L,I,J),YPSI(K,L,NPERIO(I-1),J),      &
                                  & YPSI(0,L,NPERIO(I+1),J),YPSI(0,L,I,J),YPSI(0,L,NPERIO(I-1),J))            &
                                  & - YPSI(K,L,I,J) )/( (ONE/PM0(L,I,J))*(                                    & 
                                  & (PDT/RDX)*(MAX(ZERO,YXDOT(K+1,L,NPERIO(I-1),J))*YPSI(K,L,NPERIO(I-1),J)   & 
                                  & -MIN(ZERO,YXDOT(K+1,L,I,J))*YPSI(K,L,NPERIO(I+1),J))                      &
                                  & +(PDT/RDY)*(MAX(ZERO,YYDOT(K+1,L,I,NPERIO(J-1)))*YPSI(K,L,I,NPERIO(J-1))  &
                                  & -MIN(ZERO,YYDOT(K+1,L,I,J))*YPSI(K,L,I,NPERIO(J+1))) ) + REPS  )
              
                 ZINX_UP          = ( MAX(YPSI(K,L,NPERIO(I+2),J),YPSI(K,L,NPERIO(I+1),J),YPSI(K,L,I,J),      &
                                  & YPSI(0,L,NPERIO(I+2),J),YPSI(0,L,I,J),YPSI(0,L,NPERIO(I+1),J))            &
                                  & - YPSI(K,L,NPERIO(I+1),J) )/( (ONE/PM0(L,NPERIO(I+1),J))*(                & 
                                  & (PDT/RDX)*(MAX(ZERO,YXDOT(K+1,L,I,J))*YPSI(K,L,I,J)                       & 
                                  & -MIN(ZERO,YXDOT(K+1,L,NPERIO(I+1),J))*YPSI(K,L,NPERIO(I+2),J))            &
                                  & +(PDT/RDY)*(MAX(ZERO,YYDOT(K+1,L,NPERIO(I+1),NPERIO(J-1)))                &
                                  & *YPSI(K,L,NPERIO(I+1),NPERIO(J-1))-MIN(ZERO,YYDOT(K+1,L,NPERIO(I+1),J))   &
                                  & *YPSI(K,L,NPERIO(I+1),NPERIO(J+1))) ) + REPS  )
              
                 ZOUTX_DOWN       = -( MIN( YPSI(K,L,NPERIO(I+1),J),YPSI(K,L,I,J),YPSI(K,L,NPERIO(I-1),J),    &
                                  & YPSI(0,L,NPERIO(I+1),J),YPSI(0,L,I,J),YPSI(0,L,NPERIO(I-1),J))            &
                                  & - YPSI(K,L,I,J) )/( (ONE/PM0(L,I,J))*(                                    &
                                  &  (PDT/RDX)*(MAX(ZERO,YXDOT(K+1,L,I,J))*YPSI(K,L,I,J)                      &
                                  & -MIN(ZERO,YXDOT(K+1,L,NPERIO(I-1),J))*YPSI(K,L,I,J))                      &
                                  & +(PDT/RDY)*(MAX(ZERO,YYDOT(K+1,L,I,J))*YPSI(K,L,I,J)                      &
                                  & -MIN(ZERO,YYDOT(K+1,L,I,NPERIO(J-1)))*YPSI(K,L,I,J)) ) + REPS  )

                 ZOUTX_UP         = -( MIN( YPSI(K,L,NPERIO(I+2),J),YPSI(K,L,I,J),YPSI(K,L,NPERIO(I+1),J),    &
                                  & YPSI(0,L,NPERIO(I+2),J),YPSI(0,L,I,J),YPSI(0,L,NPERIO(I+1),J))            &
                                  & - YPSI(K,L,NPERIO(I+1),J) )/( (ONE/PM0(L,NPERIO(I+1),J))*(                &  
                                  &  (PDT/RDX)*(MAX(ZERO,YXDOT(K+1,L,NPERIO(I+1),J))*YPSI(K,L,NPERIO(I+1),J)  &
                                  & -MIN(ZERO,YXDOT(K+1,L,I,J))*YPSI(K,L,NPERIO(I+1),J))                      &
                                  & +(PDT/RDY)*(MAX(ZERO,YYDOT(K+1,L,NPERIO(I+1),J))*YPSI(K,L,NPERIO(I+1),J)  &
                                  & -MIN(ZERO,YYDOT(K+1,L,NPERIO(I+1),NPERIO(J-1)))*YPSI(K,L,NPERIO(I+1),J))) &
                                  & + REPS  )

                 ZXDOT_MONO       = YXDOT(K+1,L,I,J)
                 YXDOT(K+1,L,I,J) = MIN(ONE,ZOUTX_DOWN,ZINX_UP)*MAX(ZERO,ZXDOT_MONO)                          & 
                                  & + MIN(ONE,ZINX_DOWN,ZOUTX_UP)*MIN(ZERO,ZXDOT_MONO)
              

                 !***************!
                 !* Y-DIRECTION *!
                 !***************!
                 
                 ZINY_DOWN        = ( MAX( YPSI(K,L,I,NPERIO(J+1)),YPSI(K,L,I,J),YPSI(K,L,I,NPERIO(J-1)),     &
                                  & YPSI(0,L,I,NPERIO(J+1)),YPSI(0,L,I,J),YPSI(0,L,I,NPERIO(J-1)))            &
                                  & - YPSI(K,L,I,J) )/( (ONE/PM0(L,I,J))*(                                    & 
                                  & (PDT/RDX)*(MAX(ZERO,YXDOT(K+1,L,NPERIO(I-1),J))*YPSI(K,L,NPERIO(I-1),J)   & 
                                  & -MIN(ZERO,YXDOT(K+1,L,I,J))*YPSI(K,L,NPERIO(I+1),J))                      &
                                  & +(PDT/RDY)*(MAX(ZERO,YYDOT(K+1,L,I,NPERIO(J-1)))*YPSI(K,L,I,NPERIO(J-1))  &
                                  & -MIN(ZERO,YYDOT(K+1,L,I,J))*YPSI(K,L,I,NPERIO(J+1))) ) + REPS  )
              
                 ZINY_UP          = ( MAX( YPSI(K,L,I,NPERIO(J+2)),YPSI(K,L,I,J),YPSI(K,L,I,NPERIO(J+1)),     &
                                  & YPSI(0,L,I,NPERIO(J+2)),YPSI(0,L,I,J),YPSI(0,L,I,NPERIO(J+1)))            &
                                  & - YPSI(K,L,I,NPERIO(J+1)) )/( (ONE/PM0(L,I,NPERIO(J+1)))*(                & 
                                  & (PDT/RDX)*(MAX(ZERO,YXDOT(K+1,L,NPERIO(I-1),NPERIO(J+1)))                 &
                                  & *YPSI(K,L,NPERIO(I-1),NPERIO(J+1))                                        & 
                                  & -MIN(ZERO,YXDOT(K+1,L,I,NPERIO(J+1)))*YPSI(K,L,NPERIO(I+1),NPERIO(J+1)))  &
                                  & +(PDT/RDY)*(MAX(ZERO,YYDOT(K+1,L,I,J))*YPSI(K,L,I,J)                      &
                                  & -MIN(ZERO,YYDOT(K+1,L,I,NPERIO(J+1)))*YPSI(K,L,I,NPERIO(J+2))) ) + REPS  )
        
                 ZOUTY_DOWN       = -( MIN( YPSI(K,L,I,NPERIO(J+1)),YPSI(K,L,I,J),YPSI(K,L,I,NPERIO(J-1)),    &
                                  & YPSI(0,L,I,NPERIO(J+1)),YPSI(0,L,I,J),YPSI(0,L,I,NPERIO(J-1)))            &
                                  & - YPSI(K,L,I,J) )/( (ONE/PM0(L,I,J))*(                                    &
                                  &  (PDT/RDX)*(MAX(ZERO,YXDOT(K+1,L,I,J))*YPSI(K,L,I,J)                      &
                                  & -MIN(ZERO,YXDOT(K+1,L,NPERIO(I-1),J))*YPSI(K,L,I,J))                      &
                                  & +(PDT/RDY)*(MAX(ZERO,YYDOT(K+1,L,I,J))*YPSI(K,L,I,J)                      &
                                  & -MIN(ZERO,YYDOT(K+1,L,I,NPERIO(J-1)))*YPSI(K,L,I,J)) ) + REPS  )

                 ZOUTY_UP         = -( MIN( YPSI(K,L,I,NPERIO(J+2)),YPSI(K,L,I,J),YPSI(K,L,I,NPERIO(J+1)),    &
                                  & YPSI(0,L,I,NPERIO(J+2)),YPSI(0,L,I,J),YPSI(0,L,I,NPERIO(J+1)))            &
                                  & - YPSI(K,L,I,NPERIO(J+1)) )/( (ONE/PM0(L,I,NPERIO(J+1)))*(                &
                                  &  (PDT/RDX)*(MAX(ZERO,YXDOT(K+1,L,I,NPERIO(J+1)))*YPSI(K,L,I,NPERIO(J+1))  &
                                  & -MIN(ZERO,YXDOT(K+1,L,NPERIO(I-1),NPERIO(J+1)))*YPSI(K,L,I,NPERIO(J+1)))  &
                                  & +(PDT/RDY)*(MAX(ZERO,YYDOT(K+1,L,I,NPERIO(J+1)))*YPSI(K,L,I,NPERIO(J+1))  &
                                  & -MIN(ZERO,YYDOT(K+1,L,I,J))*YPSI(K,L,I,NPERIO(J+1))) ) + REPS  )

                 ZYDOT_MONO       = YYDOT(K+1,L,I,J)
                 YYDOT(K+1,L,I,J) = MIN(ONE,ZOUTY_DOWN,ZINY_UP)*MAX(ZERO,ZYDOT_MONO)                          & 
                                  & +MIN(ONE,ZINY_DOWN,ZOUTY_UP)*MIN(ZERO,ZYDOT_MONO)
            

              END DO
           END DO   
        END DO
     END IF
                        
     !****************************************************************************************************!
     !* CORRECTED UPWIND SCHEME *!
     !***************************!
     DO I=1,GXPT
        DO J=1,GYPT
           DO L=KTOP,NLEV
              
              YPSI(K+1,L,I,J) = YPSI(K,L,I,J)-(ONE/PM0(L,I,J))*(PDT/RDX)*(                               &
                              &  (MIN(ZERO,YXDOT(K+1,L,I,J))*YPSI(K,L,NPERIO(I+1),J)                     &
                              &  +MAX(ZERO,YXDOT(K+1,L,I,J))*YPSI(K,L,I,J))                              &
                              & -(MIN(ZERO,YXDOT(K+1,L,NPERIO(I-1),J))*YPSI(K,K,I,J)                     &
                              &  +MAX(ZERO,YXDOT(K+1,L,NPERIO(I-1),J))*YPSI(K,L,NPERIO(I-1),J)) )        &
                              & -(ONE/PM0(L,I,J))*(PDT/RDY)*(                                            &
                              &  (MIN(ZERO,YYDOT(K+1,L,I,J))*YPSI(K,L,I,NPERIO(J+1))                     &
                              &  +MAX(ZERO,YYDOT(K+1,L,I,J))*YPSI(K,L,I,J))                              &
                              & -(MIN(ZERO,YYDOT(K+1,L,I,NPERIO(J-1)))*YPSI(K,L,I,J)                     &
                              &  +MAX(ZERO,YYDOT(K+1,L,I,NPERIO(J-1)))*YPSI(K,L,I,NPERIO(J-1))) )
                          
              YPSIG(L,I,J)    = (PM0(L,I,J)/PM1(L,I,J))*YPSI(K+1,L,I,J)
              
           END DO   
        END DO   
     END DO

  END DO
  
  !*************************!
  ! SECOND CORRECTION-STEP *!
  !*************************!

  DO K=NORD1,NORD2-1
     
     !***********************************************************************************************!
     !* COMPUTE PSEUDO-VELOCITIES *!
     !*****************************!
     DO I=1,GXPT
        DO J=1,GYPT
           DO L=KTOP,NLEV
              
              YXDOT(K+1,L,I,J) = HALF*ABS(YXDOT(K,L,I,J))*(YPSIG(L,NPERIO(I+1),J)-YPSIG(L,I,J))      &
                          & -(PDT/TWO)*YXDOT(K,L,I,J)*(                                              &
                          & ( (HALF*((YXDOT(K,L,I,J)*(YPSIG(L,NPERIO(I+1),J)                         & ! div(i,j)
                          & +YPSIG(L,I,J))-YXDOT(K,L,NPERIO(I-1),J)*(                                &
                          & YPSIG(L,NPERIO(I-1),J)+YPSIG(L,I,J)) )/RDX))                             &
                          & +(HALF*((YYDOT(K,L,I,J)*(YPSIG(L,I,NPERIO(J+1))                          & 
                          & +YPSIG(L,I,J))-YYDOT(K,L,I,NPERIO(J-1))*(                                &
                          & YPSIG(L,I,NPERIO(J-1))+YPSIG(L,I,J)) )/RDY)) )                           &
                          & +( (HALF*((YXDOT(K,L,NPERIO(I+1),J)*(                                    & ! div(i+1,j)
                          & YPSIG(L,NPERIO(I+2),J)+YPSIG(L,NPERIO(I+1),J) )                          &
                          & -YXDOT(K,L,I,J)*(YPSIG(L,NPERIO(I+1),J)+YPSIG(L,I,J)) )/RDX))            &
                          & +(HALF*((YYDOT(K,L,NPERIO(I+1),J)*(                                      &  
                          & YPSIG(L,NPERIO(I+1),NPERIO(J+1))+YPSIG(L,NPERIO(I+1),J))                 & 
                          & -YYDOT(K,L,NPERIO(I+1),NPERIO(J-1))*(YPSIG(L,NPERIO(I+1),J)              &
                          & +YPSIG(L,NPERIO(I+1),NPERIO(J-1))))/RDY)) )                              &
                          & +HALF*(YPSIG(L,NPERIO(I+1),J)+YPSIG(L,I,J))                              &
                          & *((PM1(L,NPERIO(I+1),J)-PM0(L,NPERIO(I+1),J)                             & ! <DG/DT>(i+0.5,j)
                          & + PM1(L,I,J)-PM0(L,I,J) )/PDT) )/(PM0(L,NPERIO(I+1),J)+PM0(L,I,J))

              YYDOT(K+1,L,I,J) = HALF*ABS(YYDOT(K,L,I,J))*(YPSIG(L,I,NPERIO(J+1))-YPSIG(L,I,J))      &
                          & -(PDT/TWO)*YYDOT(K,L,I,J)*(                                              &
                          & ( (HALF*((YXDOT(K,L,I,J)*(YPSIG(L,NPERIO(I+1),J)                         & ! div(i,j)
                          & +YPSIG(L,I,J))-YXDOT(K,L,NPERIO(I-1),J)*(                                &
                          & YPSIG(L,NPERIO(I-1),J)+YPSIG(L,I,J)) )/RDX))                             &
                          & +(HALF*((YYDOT(K,L,I,J)*(YPSIG(L,I,NPERIO(J+1))                          & 
                          & +YPSIG(L,I,J))-YYDOT(K,L,I,NPERIO(J-1))*(                                &
                          & YPSIG(L,I,NPERIO(J-1))+YPSIG(L,I,J)) )/RDY)) )                           &
                          & +( (HALF*((YXDOT(K,L,I,NPERIO(J+1))*(                                    & ! div(i,j+1)
                          & YPSIG(L,NPERIO(I+1),NPERIO(J+1))+YPSIG(L,I,NPERIO(J+1)) )                &
                          & -YXDOT(K,L,NPERIO(I-1),NPERIO(J+1))*(YPSIG(L,I,NPERIO(J+1))              &
                          & +YPSIG(L,NPERIO(I-1),NPERIO(J+1))) )/RDX))                               &
                          & +(HALF*((YYDOT(K,L,I,NPERIO(J+1))*(YPSIG(L,I,NPERIO(J+2))                & 
                          & +YPSIG(L,I,NPERIO(J+1)))-YYDOT(K,L,I,J)*(YPSIG(L,I,NPERIO(J+1))          &
                          & +YPSIG(L,I,J)))/RDY)) )                                                  &
                          & +HALF*(YPSIG(L,I,NPERIO(J+1))+YPSIG(L,I,J))                              &
                          & *((PM1(L,I,NPERIO(J+1))-PM0(L,I,NPERIO(J+1))                             & ! <DG/DT>(i,j+0.5)
                          & + PM1(L,I,J)-PM0(L,I,J) )/PDT) )/(PM0(L,I,NPERIO(J+1))+PM0(L,I,J))

           END DO
        END DO        
     END DO
     
     !**********************************************************************************************************!
     !* GLOBAL FCT NON-OSCILLATORY TREATMENT *!
     !****************************************!
     IF (LFCT_MONO) THEN
        DO I=1,GXPT
           DO J=1,GYPT
              DO L=KTOP,NLEV

                 !***************!
                 !* X-DIRECTION *!
                 !***************!
                 
                 ZINX_DOWN        = ( MAX(YPSI(K,L,NPERIO(I+1),J),YPSI(K,L,I,J),YPSI(K,L,NPERIO(I-1),J),      &
                                  & YPSI(0,L,NPERIO(I+1),J),YPSI(0,L,I,J),YPSI(0,L,NPERIO(I-1),J))            &
                                  & - YPSI(K,L,I,J) )/( (ONE/PM0(L,I,J))*(                                    & 
                                  & (PDT/RDX)*(MAX(ZERO,YXDOT(K+1,L,NPERIO(I-1),J))*YPSI(K,L,NPERIO(I-1),J)   & 
                                  & -MIN(ZERO,YXDOT(K+1,L,I,J))*YPSI(K,L,NPERIO(I+1),J))                      &
                                  & +(PDT/RDY)*(MAX(ZERO,YYDOT(K+1,L,I,NPERIO(J-1)))*YPSI(K,L,I,NPERIO(J-1))  &
                                  & -MIN(ZERO,YYDOT(K+1,L,I,J))*YPSI(K,L,I,NPERIO(J+1))) ) + REPS  )
              
                 ZINX_UP          = ( MAX(YPSI(K,L,NPERIO(I+2),J),YPSI(K,L,NPERIO(I+1),J),YPSI(K,L,I,J),      &
                                  & YPSI(0,L,NPERIO(I+2),J),YPSI(0,L,I,J),YPSI(0,L,NPERIO(I+1),J))            &
                                  & - YPSI(K,L,NPERIO(I+1),J) )/( (ONE/PM0(L,NPERIO(I+1),J))*(                & 
                                  & (PDT/RDX)*(MAX(ZERO,YXDOT(K+1,L,I,J))*YPSI(K,L,I,J)                       & 
                                  & -MIN(ZERO,YXDOT(K+1,L,NPERIO(I+1),J))*YPSI(K,L,NPERIO(I+2),J))            &
                                  & +(PDT/RDY)*(MAX(ZERO,YYDOT(K+1,L,NPERIO(I+1),NPERIO(J-1)))                &
                                  & *YPSI(K,L,NPERIO(I+1),NPERIO(J-1))-MIN(ZERO,YYDOT(K+1,L,NPERIO(I+1),J))   &
                                  & *YPSI(K,L,NPERIO(I+1),NPERIO(J+1))) ) + REPS  )
              
                 ZOUTX_DOWN       = -( MIN( YPSI(K,L,NPERIO(I+1),J),YPSI(K,L,I,J),YPSI(K,L,NPERIO(I-1),J),    &
                                  & YPSI(0,L,NPERIO(I+1),J),YPSI(0,L,I,J),YPSI(0,L,NPERIO(I-1),J))            &
                                  & - YPSI(K,L,I,J) )/( (ONE/PM0(L,I,J))*(                                    &
                                  &  (PDT/RDX)*(MAX(ZERO,YXDOT(K+1,L,I,J))*YPSI(K,L,I,J)                      &
                                  & -MIN(ZERO,YXDOT(K+1,L,NPERIO(I-1),J))*YPSI(K,L,I,J))                      &
                                  & +(PDT/RDY)*(MAX(ZERO,YYDOT(K+1,L,I,J))*YPSI(K,L,I,J)                      &
                                  & -MIN(ZERO,YYDOT(K+1,L,I,NPERIO(J-1)))*YPSI(K,L,I,J)) ) + REPS  )

                 ZOUTX_UP         = -( MIN( YPSI(K,L,NPERIO(I+2),J),YPSI(K,L,I,J),YPSI(K,L,NPERIO(I+1),J),    &
                                  & YPSI(0,L,NPERIO(I+2),J),YPSI(0,L,I,J),YPSI(0,L,NPERIO(I+1),J))            &
                                  & - YPSI(K,L,NPERIO(I+1),J) )/( (ONE/PM0(L,NPERIO(I+1),J))*(                &  
                                  &  (PDT/RDX)*(MAX(ZERO,YXDOT(K+1,L,NPERIO(I+1),J))*YPSI(K,L,NPERIO(I+1),J)  &
                                  & -MIN(ZERO,YXDOT(K+1,L,I,J))*YPSI(K,L,NPERIO(I+1),J))                      &
                                  & +(PDT/RDY)*(MAX(ZERO,YYDOT(K+1,L,NPERIO(I+1),J))*YPSI(K,L,NPERIO(I+1),J)  &
                                  & -MIN(ZERO,YYDOT(K+1,L,NPERIO(I+1),NPERIO(J-1)))*YPSI(K,L,NPERIO(I+1),J))) &
                                  & + REPS  )

                 ZXDOT_MONO       = YXDOT(K+1,L,I,J)
                 YXDOT(K+1,L,I,J) = MIN(ONE,ZOUTX_DOWN,ZINX_UP)*MAX(ZERO,ZXDOT_MONO)                          & 
                                  & + MIN(ONE,ZINX_DOWN,ZOUTX_UP)*MIN(ZERO,ZXDOT_MONO)
              

                 !***************!
                 !* Y-DIRECTION *!
                 !***************!
              
                 ZINY_DOWN        = ( MAX( YPSI(K,L,I,NPERIO(J+1)),YPSI(K,L,I,J),YPSI(K,L,I,NPERIO(J-1)),    &
                                  & YPSI(0,L,I,NPERIO(J+1)),YPSI(0,L,I,J),YPSI(0,L,I,NPERIO(J-1)))           &
                                  & - YPSI(K,L,I,J) )/( (ONE/PM0(L,I,J))*(                                   & 
                                  & (PDT/RDX)*(MAX(ZERO,YXDOT(K+1,L,NPERIO(I-1),J))*YPSI(K,L,NPERIO(I-1),J)  & 
                                  & -MIN(ZERO,YXDOT(K+1,L,I,J))*YPSI(K,L,NPERIO(I+1),J))                     &
                                  & +(PDT/RDY)*(MAX(ZERO,YYDOT(K+1,L,I,NPERIO(J-1)))*YPSI(K,L,I,NPERIO(J-1)) &
                                  & -MIN(ZERO,YYDOT(K+1,L,I,J))*YPSI(K,L,I,NPERIO(J+1))) ) + REPS  )
              
                 ZINY_UP          = ( MAX( YPSI(K,L,I,NPERIO(J+2)),YPSI(K,L,I,J),YPSI(K,L,I,NPERIO(J+1)),    &
                                  & YPSI(0,L,I,NPERIO(J+2)),YPSI(0,L,I,J),YPSI(0,L,I,NPERIO(J+1)))           &
                                  & - YPSI(K,L,I,NPERIO(J+1)) )/( (ONE/PM0(L,I,NPERIO(J+1)))*(               & 
                                  & (PDT/RDX)*(MAX(ZERO,YXDOT(K+1,L,NPERIO(I-1),NPERIO(J+1)))                &
                                  & *YPSI(K,L,NPERIO(I-1),NPERIO(J+1))  & 
                                  & -MIN(ZERO,YXDOT(K+1,L,I,NPERIO(J+1)))*YPSI(K,L,NPERIO(I+1),NPERIO(J+1))) &
                                  & +(PDT/RDY)*(MAX(ZERO,YYDOT(K+1,L,I,J))*YPSI(K,L,I,J)                     &
                                  & -MIN(ZERO,YYDOT(K+1,L,I,NPERIO(J+1)))*YPSI(K,L,I,NPERIO(J+2))) ) + REPS  )
        
                 ZOUTY_DOWN       = -( MIN( YPSI(K,L,I,NPERIO(J+1)),YPSI(K,L,I,J),YPSI(K,L,I,NPERIO(J-1)),   &
                                  & YPSI(0,L,I,NPERIO(J+1)),YPSI(0,L,I,J),YPSI(0,L,I,NPERIO(J-1)))           &
                                  & - YPSI(K,L,I,J) )/( (ONE/PM0(L,I,J))*(                                   &
                                  &  (PDT/RDX)*(MAX(ZERO,YXDOT(K+1,L,I,J))*YPSI(K,L,I,J)                     &
                                  & -MIN(ZERO,YXDOT(K+1,L,NPERIO(I-1),J))*YPSI(K,L,I,J))                     &
                                  & +(PDT/RDY)*(MAX(ZERO,YYDOT(K+1,L,I,J))*YPSI(K,L,I,J)                     &
                                  & -MIN(ZERO,YYDOT(K+1,L,I,NPERIO(J-1)))*YPSI(K,L,I,J)) ) + REPS  )

                 ZOUTY_UP         = -( MIN( YPSI(K,L,I,NPERIO(J+2)),YPSI(K,L,I,J),YPSI(K,L,I,NPERIO(J+1)),   &
                                  & YPSI(0,L,I,NPERIO(J+2)),YPSI(0,L,I,J),YPSI(0,L,I,NPERIO(J+1)))           &
                                  & - YPSI(K,L,I,NPERIO(J+1)) )/( (ONE/PM0(L,I,NPERIO(J+1)))*(               &
                                  &  (PDT/RDX)*(MAX(ZERO,YXDOT(K+1,L,I,NPERIO(J+1)))*YPSI(K,L,I,NPERIO(J+1)) &
                                  & -MIN(ZERO,YXDOT(K+1,L,NPERIO(I-1),NPERIO(J+1)))*YPSI(K,L,I,NPERIO(J+1))) &
                                  & +(PDT/RDY)*(MAX(ZERO,YYDOT(K+1,L,I,NPERIO(J+1)))*YPSI(K,L,I,NPERIO(J+1)) &
                                  & -MIN(ZERO,YYDOT(K+1,L,I,J))*YPSI(K,L,I,NPERIO(J+1))) ) + REPS  )

                 ZYDOT_MONO       = YYDOT(K+1,L,I,J)
                 YYDOT(K+1,L,I,J) = MIN(ONE,ZOUTY_DOWN,ZINY_UP)*MAX(ZERO,ZYDOT_MONO)                         & 
                                  & +MIN(ONE,ZINY_DOWN,ZOUTY_UP)*MIN(ZERO,ZYDOT_MONO)
            

              END DO
           END DO   
        END DO
     END IF

     !************************************************************************************************!
     !* CORRECTED UPWIND SCHEME *!
     !***************************!
     DO I=1,GXPT
        DO J=1,GYPT
           DO J=KTOP,NLEV
              YPSI(K+1,L,I,J)  = YPSI(K,J,I) - (ONE/PM0(J,I))*(                                     &
                               &  (PDT/RDX)*(YXDOT(K+1,L,I,J)-YXDOT(K+1,L,NPERIO(I-1),J))           &
                               & +(PDT/RDY)*(YYDOT(K+1,L,I,J)-YYDOT(K+1,L,I,NPERIO(J-1)))  )           
              YPSIG(L,I,J)     = (PM0(L,I,J)/PM1(L,I,J))*YPSI(K+1,L,I,J)
           END DO   
        END DO   
     END DO
     
  END DO
 
  !***************************************************************************************************!
  DO I=1,GXPT
     DO J=1,GYPT
        
        DO L=KTOP,NLEV
           PPSI_DEP(L,I,J)  = (PM0(L,I,J)/PM1(L,I,J))*(YPSI(NORD2,L,I,J)+ZPSI_MIN)
        END DO
  
        IF (LLMASS) THEN
           PXDOT(:,I,J) = HALF*YXDOT(NORD2,:,I,J)*(                          &
                        & HALF*(PPSI_ARR(:,NPERIO(I+1),J)+PPSI_DEP(:,I,J))   &
                        & + HALF*(PPSI_ARR(:,I,J)+PPSI_DEP(:,NPERIO(I+1),J)) )
           PYDOT(:,I,J) = HALF*YYDOT(NORD2,:,I,J)*(                          &
                        & HALF*(PPSI_ARR(:,I,NPERIO(J+1))+PPSI_DEP(:,I,J))   &
                        & + HALF*(PPSI_ARR(:,I,J)+PPSI_DEP(:,I,NPERIO(J+1))) )
        END IF
        
     END DO   
  END DO 
  
END SUBROUTINE ADV_XY
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
SUBROUTINE MPDATA_SURF_PRESSURE(PPIS_DEP,PPIS_ARR,PXDOTS,PYDOTS)
  
  IMPLICIT NONE

  REAL(8),                       INTENT(OUT)    :: PPIS_DEP(GXPT,GYPT)
  REAL(8),                       INTENT(IN)     :: PPIS_ARR(GXPT,GYPT)
  REAL(8),                       INTENT(IN)     :: PXDOTS(GXPT,GYPT)
  REAL(8),                       INTENT(IN)     :: PYDOTS(GXPT,GYPT)

  INTEGER(8)                                    :: I,J,K
  REAL(8)                                       :: YXDOTS(NORD2,GXPT,GYPT)
  REAL(8)                                       :: YYDOTS(NORD2,GXPT,GYPT)
  REAL(8)                                       :: YPIS(0:NORD2,GXPT,GYPT)
  REAL(8)                                       :: ZIN_DOWN,ZIN_UP,ZOUT_DOWN,ZOUT_UP
  REAL(8)                                       :: ZPIS_MIN,ZXDOTS_MONO,ZYDOTS_MONO

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ZPIS_MIN = MINVAL(PPIS_ARR)
  
  DO I=1,GXPT
     DO J=1,GYPT
        YXDOTS(1,I,J) = PXDOTS(I,J)
        YYDOTS(1,I,J) = PYDOTS(I,J)
        YPIS(0,I,J)   = PPIS_ARR(I,J) - ZPIS_MIN
     END DO   
  END DO
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !***************************************!
  ! DONOR-CELL UPWIND FIRST-ORDER SCHEME *!
  !***************************************!

  DO I=1,GXPT
     DO J=1,GYPT
        YPIS(1,I,J)   = YPIS(0,I,J)                                                         &
                      & - (RDT/RDX)*( (MIN(ZERO,YXDOTS(1,I,J))*YPIS(0,NPERIO(I+1),J)        &
                      &  +MAX(ZERO,YXDOTS(1,I,J))*YPIS(0,I,J))                              &
                      & -(MIN(ZERO,YXDOTS(1,NPERIO(I-1),J))*YPIS(0,I,J)                     &
                      &  +MAX(ZERO,YXDOTS(1,NPERIO(I-1),J))*YPIS(0,NPERIO(I-1),J)) )        &
                      & - (RDT/RDY)*( (MIN(ZERO,YYDOTS(1,I,J))*YPIS(0,I,NPERIO(J+1))        &
                      &  +MAX(ZERO,YYDOTS(1,I,J))*YPIS(0,I,J))                              &
                      & -(MIN(ZERO,YYDOTS(1,I,NPERIO(J-1)))*YPIS(0,I,J)                     &
                      &  +MAX(ZERO,YYDOTS(1,I,NPERIO(J-1)))*YPIS(0,I,NPERIO(J-1))) )
     END DO   
  END DO
  
  !************************!
  ! FIRST CORRECTION-STEP *!
  !************************!

  DO K=1,NORD1-1
     
     !**************************************************************************************!
     !* COMPUTE SURFACE PSEUDO-VELOCITIES
     DO I=1,GXPT
        DO J=1,GYPT
           
           YXDOTS(K+1,I,J) = HALF*ABS(YXDOTS(K,I,J))*(                                      & 
                      &  (ABS(YPIS(K,NPERIO(I+1),J))-ABS(YPIS(K,I,J)))                      &
                      &  /(ABS(YPIS(K,I,J))+ABS(YPIS(K,NPERIO(I+1),J))+REPS) )              &
                      & -(RDT/TWO)*YXDOTS(K,I,J)*( ( (HALF*( ( YXDOTS(K,I,J)*(              &
                      & ABS(YPIS(K,NPERIO(I+1),J))+ABS(YPIS(K,I),J))                        &
                      & -YXDOTS(K,NPERIO(I-1),J)*(ABS(YPIS(K,NPERIO(I-1),J))                &   
                      & +ABS(YPIS(K,I,J))) )/RDX)) + (HALF*( ( YYDOTS(K,I,J)*(              &
                      & ABS(YPIS(K,I,NPERIO(J+1)))+ABS(YPIS(K,I,J))                         &
                      & -YYDOTS(K,I,NPERIO(J-1))*(ABS(YPIS(K,I,NPERIO(J-1)))                &   
                      & +ABS(YPIS(K,I,J))) )/RDY)) )/( (0.125d0*(                           &
                      &  ABS(YPIS(K,NPERIO(I+1),J))+ABS(YPIS(K,NPERIO(I-1),J))              &
                      & +ABS(YPIS(K,I,NPERIO(J+1)))+ABS(YPIS(K,I,NPERIO(J-1)))              &
                      & +4.d0*ABS(YPIS(K,I,J))))+REPS))                                     &    
                      & +( ((HALF*( ( YXDOTS(K,NPERIO(I+1),J)*(                             &
                      & ABS(YPIS(K,NPERIO(I+2),J))+ABS(YPIS(K,NPERIO(I+1),J)))              &
                      & -YXDOTS(K,I,J)*(ABS(YPIS(K,NPERIO(I+1),J))+ABS(YPIS(K,I,J))))/RDX)) &
                      & + (HALF*( (YYDOTS(K,NPERIO(I+1),J)*(                                &
                      & ABS(YPIS(K,NPERIO(I+1),NPERIO(J+1)))+ABS(YPIS(K,NPERIO(I+1),J)))    &
                      & -YYDOTS(K,NPERIO(I+1),J)*(ABS(YPIS(K,NPERIO(I+1),J))                &
                      & +ABS(YPIS(K,NPERIO(I+1),NPERIO(J-1)))) )/RDY)) )/( (0.25d0*(        &  
                      & ABS(YPIS(K,NPERIO(I+2),J))+ABS(YPIS(K,NPERIO(I+1),NPERIO(J+1)))     & 
                      & +4.d0*ABS(YPIS(K,NPERIO(I+1),J))+ABS(YPIS(K,NPERIO(I+1),J))         &
                      & ABS(YPIS(K,NPERIO(I+1),NPERIO(J-1)))+ABS(YPIS(K,I,J))))+REPS))   )

           YYDOTS(K+1,I,J) = HALF*ABS(YYDOTS(K,I,J))*(                                      & 
                      &  (ABS(YPIS(K,I,NPERIO(J+1)))-ABS(YPIS(K,I,J)))                      &
                      &  /(ABS(YPIS(K,I,J))+ABS(YPIS(K,I,NPERIO(J+1)))+REPS) )              &
                      & -(RDT/TWO)*YYDOTS(K,I,J)*( ( (HALF*( ( YXDOTS(K,I,J)*(              &
                      & ABS(YPIS(K,NPERIO(I+1),J))+ABS(YPIS(K,I),J))                        &
                      & -YXDOTS(K,NPERIO(I-1),J)*(ABS(YPIS(K,NPERIO(I-1),J))                &   
                      & +ABS(YPIS(K,I,J))) )/RDX)) + (HALF*( ( YYDOTS(K,I,J)*(              &
                      & ABS(YPIS(K,I,NPERIO(J+1)))+ABS(YPIS(K,I,J))                         &
                      & -YYDOTS(K,I,NPERIO(J-1))*(ABS(YPIS(K,I,NPERIO(J-1)))                &   
                      & +ABS(YPIS(K,I,J))) )/RDY)) )/( (0.125d0*(                           &
                      &  ABS(YPIS(K,NPERIO(I+1),J))+ABS(YPIS(K,NPERIO(I-1),J))              &
                      & +ABS(YPIS(K,I,NPERIO(J+1)))+ABS(YPIS(K,I,NPERIO(J-1)))              &
                      & +4.d0*ABS(YPIS(K,I,J))))+REPS))                                     &    
                      & +( ((HALF*( ( YXDOTS(K,I,NPERIO(J+1))*(                             &
                      & ABS(YPIS(K,NPERIO(I+1),NPERIO(J+1)))+ABS(YPIS(K,I,NPERIO(J+1))))    &
                      & -YXDOTS(K,NPERIO(I-1),NPERIO(J+1))*(ABS(YPIS(K,I,NPERIO(J+1)))      &
                      & +ABS(YPIS(K,NPERIO(I-1),NPERIO(J+1)))))/RDX))                       &
                      & +(HALF*( (YYDOTS(K,I,NPERIO(J+1))*(                                 &
                      & ABS(YPIS(K,I,NPERIO(J+2)))+ABS(YPIS(K,I,NPERIO(J+1))))              &
                      & -YYDOTS(K,I,J)*(ABS(YPIS(K,I,NPERIO(J+1)))                          &
                      & +ABS(YPIS(K,I,J))))/RDY)))/( (0.125d0*(ABS(YPIS(K,I,NPERIO(J+2)))   &
                      & +ABS(YPIS(K,NPERIO(I-1),NPERIO(J+1)))                               &
                      & +4.d0*ABS(YPIS(K,I,NPERIO(J+1)))                                    &
                      & +ABS(YPIS(K,NPERIO(I+1),NPERIO(J+1)))+ABS(YPIS(K,I,J))))+REPS))  )
           
        END DO   
     END DO

     !**************************************************************************************!
     !* GLOBAL FCT NON-OSCILLATORY TREATMENT
     IF (LFCT_MONO) THEN
        DO I=1,GXPT
           DO J=1,GYPT

              !***************!
              !* X-DIRECTION *!
              !***************!
           
              ZIN_DOWN        = ( MAX( YPIS(K,NPERIO(I+1),J),YPIS(K,I,J),YPIS(K,NPERIO(I-1),J),         &
                              & YPIS(0,NPERIO(I+1),J),YPIS(0,I,J),YPIS(0,NPERIO(I-1),J))-YPIS(K,I,J) )  &
                              & /( (RDT/RDX)*(MAX(ZERO,YXDOTS(K+1,NPERIO(I-1),J))*YPIS(K,NPERIO(I-1),J) &
                              & -MIN(ZERO,YXDOTS(K+1,I,J))*YPIS(K,NPERIO(I+1),J))                       &
                              & +(RDT/RDY)*(MAX(ZERO,YYDOTS(K+1,I,NPERIO(J-1)))*YPIS(K,I,NPERIO(J-1))   &
                              & -MIN(ZERO,YYDOTS(K+1,I,J))*YPIS(K,I,NPERIO(J+1))) + REPS)
        
              ZOUT_DOWN       = -( MIN( YPIS(K,NPERIO(I+1),J),YPIS(K,I,J),YPIS(K,NPERIO(I-1),J),        &
                              & YPIS(0,NPERIO(I+1),J),YPIS(0,I,J),YPIS(0,NPERIO(I-1),J))-YPIS(K,I,J) )  &
                              & /( (RDT/RDX)*(MAX(ZERO,YXDOTS(K+1,I,J))*YPIS(K,I,J)                     &
                              & -MIN(ZERO,YXDOTS(K+1,NPERIO(I-1),J))*YPIS(K,I,J))                       &
                              & (RDT/RDY)*(MAX(ZERO,YYDOTS(K+1,I,J))*YPIS(K,I,J)                        &
                              & -MIN(ZERO,YYDOTS(K+1,I,NPERIO(J-1)))*YPIS(K,I,J)) + REPS)

              ZIN_UP          = ( MAX( YPIS(K,NPERIO(I+2),J),YPIS(K,I,J),YPIS(K,NPERIO(I+1),J),         &
                              & YPIS(0,NPERIO(I+2),J),YPIS(0,I,J),YPIS(0,NPERIO(I+1),J))                &
                              & - YPIS(K,NPERIO(I+1),J) )                                               &
                              & /( (RDT/RDX)*(MAX(ZERO,YXDOTS(K+1,I,J))*YPIS(K,I,J)                     &
                              & -MIN(ZERO,YXDOTS(K+1,NPERIO(I+1),J))*YPIS(K,NPERIO(I+2),J))             &
                              & +(RDT/RDY)*(MAX(ZERO,YYDOTS(K+1,NPERIO(I+1),NPERIO(J-1)))               &
                              & *YPIS(K,NPERIO(I+1),NPERIO(J-1))                                        &
                              & -MIN(ZERO,YYDOTS(K+1,NPERIO(I+1),J))*YPIS(K,NPERIO(I+1),NPERIO(J+1)))   &
                              & + REPS)
        
              ZOUT_UP         = -( MIN( YPIS(K,NPERIO(I+1),J),YPIS(K,I,J),YPIS(K,NPERIO(I+1),J),        &
                              & YPIS(0,NPERIO(I+1),J),YPIS(0,I,J),YPIS(0,NPERIO(I+1),J))                &
                              & - YPIS(K,NPERIO(I+1),J) )                                               &
                              & /( (RDT/RDX)*(MAX(ZERO,YXDOTS(K+1,NPERIO(I+1),J))*YPIS(K,NPERIO(I+1),J) &
                              & -MIN(ZERO,YXDOTS(K+1,I,J))*YPIS(K,NPERIO(I+1),J))                       &
                              & (RDT/RDY)*(MAX(ZERO,YYDOTS(K+1,NPERIO(I+1),J))*YPIS(K,NPERIO(I+1),J)    &
                              & -MIN(ZERO,YYDOTS(K+1,NPERIO(I+1),NPERIO(J-1)))*YPIS(K,NPERIO(I+1),J))   &
                              & + REPS)

              ZXDOTS_MONO     = MIN(ONE,ZOUT_DOWN,ZIN_UP)*MAX(ZERO,YXDOTS(K+1,I,J))                     &      
                              & +MIN(ONE,ZIN_DOWN,ZOUT_UP)*MIN(ZERO,YXDOTS(K+1,I,J))
              YXDOTS(K+1,I,J) = (ONE-RFCT)*YXDOTS(K+1,I,J)+RFCT*ZXDOTS_MONO


              !***************!
              !* Y-DIRECTION *!
              !***************!
           
              ZIN_DOWN        = ( MAX( YPIS(K,I,NPERIO(J+1)),YPIS(K,I,J),YPIS(K,I,NPERIO(J-1)),         &
                              & YPIS(0,I,NPERIO(J+1)),YPIS(0,I,J),YPIS(0,I,NPERIO(J-1)))-YPIS(K,I,J) )  &
                              & /( (RDT/RDX)*(MAX(ZERO,YXDOTS(K+1,NPERIO(I-1),J))*YPIS(K,NPERIO(I-1),J) &
                              & -MIN(ZERO,YXDOTS(K+1,I,J))*YPIS(K,NPERIO(I+1),J))                       &
                              & +(RDT/RDY)*(MAX(ZERO,YYDOTS(K+1,I,NPERIO(J-1)))*YPIS(K,I,NPERIO(J-1))   &
                              & -MIN(ZERO,YYDOTS(K+1,I,J))*YPIS(K,I,NPERIO(J+1))) + REPS)
        
              ZOUT_DOWN       = -( MIN( YPIS(K,I,NPERIO(J+1)),YPIS(K,I,J),YPIS(K,I,NPERIO(J-1)),        &
                              & YPIS(0,I,NPERIO(J+1)),YPIS(0,I,J),YPIS(0,I,NPERIO(J-1)))-YPIS(K,I,J) )  &
                              & /( (RDT/RDX)*(MAX(ZERO,YXDOTS(K+1,I,J))*YPIS(K,I,J)                     &
                              & -MIN(ZERO,YXDOTS(K+1,NPERIO(I-1),J))*YPIS(K,I,J))                       &
                              & (RDT/RDY)*(MAX(ZERO,YYDOTS(K+1,I,J))*YPIS(K,I,J)                        &
                              & -MIN(ZERO,YYDOTS(K+1,I,NPERIO(J-1)))*YPIS(K,I,J)) + REPS)

              ZIN_UP          = ( MAX( YPIS(K,I,NPERIO(J+2)),YPIS(K,I,J),YPIS(K,I,NPERIO(J+1)),         &
                              & YPIS(0,I,NPERIO(J+2)),YPIS(0,I,J),YPIS(0,I,NPERIO(J+1)))                &
                              & - YPIS(K,I,NPERIO(J+1)) )                                               &
                              & /( (RDT/RDX)*(MAX(ZERO,YXDOTS(K+1,NPERIO(I-1),NPERIO(J+1)))             &
                              & *YPIS(K,NPERIO(I-1),NPERIO(J+1))                                        &
                              & -MIN(ZERO,YXDOTS(K+1,I,NPERIO(J+1)))*YPIS(K,NPERIO(I+1),NPERIO(J+1)))   &
                              & +(RDT/RDY)*(MAX(ZERO,YYDOTS(K+1,I,J))*YPIS(K,I,J)                       &
                              & -MIN(ZERO,YYDOTS(K+1,I,NPERIO(J+1)))*YPIS(K,I,NPERIO(J+2))) + REPS)
        
              ZOUT_UP         = -( MIN( YPIS(K,NPERIO(I+2),J),YPIS(K,I,J),YPIS(K,NPERIO(I+1),J),        &
                              & YPIS(0,NPERIO(I+2),J),YPIS(0,I,J),YPIS(0,NPERIO(I+1),J))                &
                              & - YPIS(K,I,NPERIO(J+1)) )                                               &
                              & /( (RDT/RDX)*(MAX(ZERO,YXDOTS(K+1,I,NPERIO(J+1)))*YPIS(K,I,NPERIO(J+1)) &
                              & -MIN(ZERO,YXDOTS(K+1,NPERIO(I-1),NPERIO(J+1)))*YPIS(K,I,NPERIO(J+1)))   &
                              & (RDT/RDY)*(MAX(ZERO,YYDOTS(K+1,I,NPERIO(J+1)))*YPIS(K,I,NPERIO(J+1))    &
                              & -MIN(ZERO,YYDOTS(K+1,I,J))*YPIS(K,I,NPERIO(J+1))) + REPS)

              ZYDOTS_MONO     = MIN(ONE,ZOUT_DOWN,ZIN_UP)*MAX(ZERO,YYDOTS(K+1,I,J))                     &      
                              & +MIN(ONE,ZIN_DOWN,ZOUT_UP)*MIN(ZERO,YYDOTS(K+1,I,J))
              YYDOTS(K+1,I,J) = (ONE-RFCT)*YYDOTS(K+1,I,J)+RFCT*ZYDOTS_MONO

           END DO
        END DO                                                                                     
     END IF
     
     !******************************************************************************************!
     !* CORRECTED UPWIND SCHEME
     DO I=1,NXPT
        DO J=1,GYPT
           YPIS(K+1,I,J) = YPIS(K,I,J)
                         & -(RDT/RDX)*( (MIN(ZERO,YXDOTS(K+1,I,J))*YPIS(K,NPERIO(I+1),J)       &
                         &  +MAX(ZERO,YXDOTS(K+1,I,J))*YPIS(K,I,J))                            &
                         & -(MIN(ZERO,YXDOTS(K+1,NPERIO(I-1),J))*YPIS(K,I,J)                   &
                         &  +MAX(ZERO,YXDOTS(K+1,NPERIO(I-1),J))*YPIS(K,NPERIO(I-1),J)) )      &
                         & -(RDT/RDY)*( (MIN(ZERO,YYDOTS(K+1,I,J))*YPIS(K,I,NPERIO(J+1))       &
                         &  +MAX(ZERO,YYDOTS(K+1,I,J))*YPIS(K,I,J))                            &
                         & -(MIN(ZERO,YYDOTS(K+1,I,NPERIO(J-1)))*YPIS(K,I,J)                   &
                         &  +MAX(ZERO,YYDOTS(K+1,I,NPERIO(J-1)))*YPIS(K,I,NPERIO(J-1))) )
        END DO                 
     END DO
     
  END DO

  
  !*************************!
  ! SECOND CORRECTION-STEP *!
  !*************************!

  DO K=NORD1,NORD2-1
     
     !**************************************************************************************!
     !* COMPUTE SURFACE PSEUDO-VELOCITIES
     DO I=1,GXPT
        DO J=1,GYPT
           YXDOTS(K+1,I,J) = HALF*ABS(YXDOTS(K,I,J))*(YPIS(K,NPERIO(I+1),J)-YPIS(K,I,J))         &
                           & -(RDT/TWO)*YXDOTS(K,I,J)*HALF*(                                     &
                           & (HALF*((YXDOTS(K,I,J)*(YPIS(K,NPERIO(I+1),J)+YPIS(K,I,J))           &
                           & -YXDOTS(K,NPERIO(I-1),J)*(YPIS(K,NPERIO(I-1),J)+YPIS(K,I,J)) )      &
                           & /RDX))                                                              &
                           & +(HALF*((YYDOTS(K,I,J)*(YPIS(K,I,NPERIO(J+1))+YPIS(K,I,J))          &
                           & -YYDOTS(K,I,NPERIO(J-1))*(YPIS(K,I,NPERIO(J-1))+YPIS(K,I,J)) )      &
                           & /RDY))                                                              &
                           & +(HALF*((YXDOTS(K,NPERIO(I+1),J)*(YPIS(K,NPERIO(I+2),J)             &
                           & +YPIS(K,NPERIO(I+1),J))-YXDOTS(K,I,J)*(YPIS(K,NPERIO(I+1),J)        &
                           & +YPIS(K,I,J)) )/RDX))                                               &
                           & +(HALF*((YYDOTS(K,NPERIO(I+1),J)*(YPIS(K,NPERIO(I+1),NPERIO(J+1))   &
                           & +YPIS(K,NPERIO(I+1),J))                                             &
                           & -YYDOTS(K,NPERIO(I+1),NPERIO(J-1))*(YPIS(K,NPERIO(I+1),J)           &
                           & +YPIS(K,NPERIO(I+1),NPERIO(J-1))) )/RDY)) )
           
           YYDOTS(K+1,I,J) = HALF*ABS(YYDOTS(K,I,J))*(YPIS(K,I,NPERIO(J+1))-YPIS(K,I,J))         &
                           & -(RDT/TWO)*YYDOTS(K,I,J)*HALF*(                                     &
                           & (HALF*((YXDOTS(K,I,J)*(YPIS(K,NPERIO(I+1),J)+YPIS(K,I,J))           &
                           & -YXDOTS(K,NPERIO(I-1),J)*(YPIS(K,NPERIO(I-1),J)+YPIS(K,I,J)) )      &
                           & /RDX))                                                              &
                           & +(HALF*((YYDOTS(K,I,J)*(YPIS(K,I,NPERIO(J+1))+YPIS(K,I,J))          &
                           & -YYDOTS(K,I,NPERIO(J-1))*(YPIS(K,I,NPERIO(J-1))+YPIS(K,I,J)) )      &
                           & /RDY))                                                              &
                           & +(HALF*((YXDOTS(K,I,NPERIO(J+1))*(YPIS(K,NPERIO(I+1),NPERIO(J+1))   &
                           & +YPIS(K,I,NPERIO(J+1)))                                             &
                           & -YXDOTS(K,NPERIO(I-1),NPERIO(J+1))*(YPIS(K,NPERIO(I-1),NPERIO(J+1)) &
                           & +YPIS(K,I,NPERIO(J+1))) )/RDX))                                     &
                           & +(HALF*((YYDOTS(K,I,NPERIO(J+1))*(YPIS(K,I,NPERIO(J+2))             &
                           & +YPIS(K,I,NPERIO(J+1)))-YYDOTS(K,I,J)*(YPIS(K,I,NPERIO(J+1))        &
                           & +YPIS(K,I,J)) )/RDY)) )
        END DO   
     END DO

     !**************************************************************************************!
     !* GLOBAL FCT NON-OSCILLATORY TREATMENT
     IF (LFCT_MONO) THEN
        DO I=1,GXPT
           DO J=1,GYPT

              !***************!
              !* X-DIRECTION *!
              !***************!
           
              ZIN_DOWN        = ( MAX( YPIS(K,NPERIO(I+1),J),YPIS(K,I,J),YPIS(K,NPERIO(I-1),J),         &
                              & YPIS(0,NPERIO(I+1),J),YPIS(0,I,J),YPIS(0,NPERIO(I-1),J))-YPIS(K,I,J) )  &
                              & /( (RDT/RDX)*(MAX(ZERO,YXDOTS(K+1,NPERIO(I-1),J))*YPIS(K,NPERIO(I-1),J) &
                              & -MIN(ZERO,YXDOTS(K+1,I,J))*YPIS(K,NPERIO(I+1),J))                       &
                              & +(RDT/RDY)*(MAX(ZERO,YYDOTS(K+1,I,NPERIO(J-1)))*YPIS(K,I,NPERIO(J-1))   &
                              & -MIN(ZERO,YYDOTS(K+1,I,J))*YPIS(K,I,NPERIO(J+1))) + REPS)
        
              ZOUT_DOWN       = -( MIN( YPIS(K,NPERIO(I+1),J),YPIS(K,I,J),YPIS(K,NPERIO(I-1),J),        &
                              & YPIS(0,NPERIO(I+1),J),YPIS(0,I,J),YPIS(0,NPERIO(I-1),J))-YPIS(K,I,J) )  &
                              & /( (RDT/RDX)*(MAX(ZERO,YXDOTS(K+1,I,J))*YPIS(K,I,J)                     &
                              & -MIN(ZERO,YXDOTS(K+1,NPERIO(I-1),J))*YPIS(K,I,J))                       &
                              & (RDT/RDY)*(MAX(ZERO,YYDOTS(K+1,I,J))*YPIS(K,I,J)                        &
                              & -MIN(ZERO,YYDOTS(K+1,I,NPERIO(J-1)))*YPIS(K,I,J)) + REPS)

              ZIN_UP          = ( MAX( YPIS(K,NPERIO(I+2),J),YPIS(K,I,J),YPIS(K,NPERIO(I+1),J),         &
                              & YPIS(0,NPERIO(I+2),J),YPIS(0,I,J),YPIS(0,NPERIO(I+1),J))                &
                              & - YPIS(K,NPERIO(I+1),J) )                                               &
                              & /( (RDT/RDX)*(MAX(ZERO,YXDOTS(K+1,I,J))*YPIS(K,I,J)                     &
                              & -MIN(ZERO,YXDOTS(K+1,NPERIO(I+1),J))*YPIS(K,NPERIO(I+2),J))             &
                              & +(RDT/RDY)*(MAX(ZERO,YYDOTS(K+1,NPERIO(I+1),NPERIO(J-1)))               &
                              & *YPIS(K,NPERIO(I+1),NPERIO(J-1))                                        &
                              & -MIN(ZERO,YYDOTS(K+1,NPERIO(I+1),J))*YPIS(K,NPERIO(I+1),NPERIO(J+1)))   &
                              & + REPS)
        
              ZOUT_UP         = -( MIN( YPIS(K,NPERIO(I+1),J),YPIS(K,I,J),YPIS(K,NPERIO(I+1),J),        &
                              & YPIS(0,NPERIO(I+1),J),YPIS(0,I,J),YPIS(0,NPERIO(I+1),J))                &
                              & - YPIS(K,NPERIO(I+1),J) )                                               &
                              & /( (RDT/RDX)*(MAX(ZERO,YXDOTS(K+1,NPERIO(I+1),J))*YPIS(K,NPERIO(I+1),J) &
                              & -MIN(ZERO,YXDOTS(K+1,I,J))*YPIS(K,NPERIO(I+1),J))                       &
                              & (RDT/RDY)*(MAX(ZERO,YYDOTS(K+1,NPERIO(I+1),J))*YPIS(K,NPERIO(I+1),J)    &
                              & -MIN(ZERO,YYDOTS(K+1,NPERIO(I+1),NPERIO(J-1)))*YPIS(K,NPERIO(I+1),J))   &
                              & + REPS)

              ZXDOTS_MONO     = MIN(ONE,ZOUT_DOWN,ZIN_UP)*MAX(ZERO,YXDOTS(K+1,I,J))                     &      
                              & +MIN(ONE,ZIN_DOWN,ZOUT_UP)*MIN(ZERO,YXDOTS(K+1,I,J))
              YXDOTS(K+1,I,J) = (ONE-RFCT)*YXDOTS(K+1,I,J)+RFCT*ZXDOTS_MONO


              !***************!
              !* Y-DIRECTION *!
              !***************!
           
              ZIN_DOWN        = ( MAX( YPIS(K,I,NPERIO(J+1)),YPIS(K,I,J),YPIS(K,I,NPERIO(J-1)),         &
                              & YPIS(0,I,NPERIO(J+1)),YPIS(0,I,J),YPIS(0,I,NPERIO(J-1)))-YPIS(K,I,J) )  &
                              & /( (RDT/RDX)*(MAX(ZERO,YXDOTS(K+1,NPERIO(I-1),J))*YPIS(K,NPERIO(I-1),J) &
                              & -MIN(ZERO,YXDOTS(K+1,I,J))*YPIS(K,NPERIO(I+1),J))                       &
                              & +(RDT/RDY)*(MAX(ZERO,YYDOTS(K+1,I,NPERIO(J-1)))*YPIS(K,I,NPERIO(J-1))   &
                              & -MIN(ZERO,YYDOTS(K+1,I,J))*YPIS(K,I,NPERIO(J+1))) + REPS)
        
              ZOUT_DOWN       = -( MIN( YPIS(K,I,NPERIO(J+1)),YPIS(K,I,J),YPIS(K,I,NPERIO(J-1)),        &
                              & YPIS(0,I,NPERIO(J+1)),YPIS(0,I,J),YPIS(0,I,NPERIO(J-1)))-YPIS(K,I,J) )  &
                              & /( (RDT/RDX)*(MAX(ZERO,YXDOTS(K+1,I,J))*YPIS(K,I,J)                     &
                              & -MIN(ZERO,YXDOTS(K+1,NPERIO(I-1),J))*YPIS(K,I,J))                       &
                              & (RDT/RDY)*(MAX(ZERO,YYDOTS(K+1,I,J))*YPIS(K,I,J)                        &
                              & -MIN(ZERO,YYDOTS(K+1,I,NPERIO(J-1)))*YPIS(K,I,J)) + REPS)

              ZIN_UP          = ( MAX( YPIS(K,I,NPERIO(J+2)),YPIS(K,I,J),YPIS(K,I,NPERIO(J+1)),         &
                              & YPIS(0,I,NPERIO(J+2)),YPIS(0,I,J),YPIS(0,I,NPERIO(J+1)))                &
                              & - YPIS(K,I,NPERIO(J+1)) )                                               &
                              & /( (RDT/RDX)*(MAX(ZERO,YXDOTS(K+1,NPERIO(I-1),NPERIO(J+1)))             &
                              & *YPIS(K,NPERIO(I-1),NPERIO(J+1))                                        &
                              & -MIN(ZERO,YXDOTS(K+1,I,NPERIO(J+1)))*YPIS(K,NPERIO(I+1),NPERIO(J+1)))   &
                              & +(RDT/RDY)*(MAX(ZERO,YYDOTS(K+1,I,J))*YPIS(K,I,J)                       &
                              & -MIN(ZERO,YYDOTS(K+1,I,NPERIO(J+1)))*YPIS(K,I,NPERIO(J+2))) + REPS)
        
              ZOUT_UP         = -( MIN( YPIS(K,NPERIO(I+2),J),YPIS(K,I,J),YPIS(K,NPERIO(I+1),J),        &
                              & YPIS(0,NPERIO(I+2),J),YPIS(0,I,J),YPIS(0,NPERIO(I+1),J))                &
                              & - YPIS(K,I,NPERIO(J+1)) )                                               &
                              & /( (RDT/RDX)*(MAX(ZERO,YXDOTS(K+1,I,NPERIO(J+1)))*YPIS(K,I,NPERIO(J+1)) &
                              & -MIN(ZERO,YXDOTS(K+1,NPERIO(I-1),NPERIO(J+1)))*YPIS(K,I,NPERIO(J+1)))   &
                              & (RDT/RDY)*(MAX(ZERO,YYDOTS(K+1,I,NPERIO(J+1)))*YPIS(K,I,NPERIO(J+1))    &
                              & -MIN(ZERO,YYDOTS(K+1,I,J))*YPIS(K,I,NPERIO(J+1))) + REPS)

              ZYDOTS_MONO     = MIN(ONE,ZOUT_DOWN,ZIN_UP)*MAX(ZERO,YYDOTS(K+1,I,J))                     &      
                              & +MIN(ONE,ZIN_DOWN,ZOUT_UP)*MIN(ZERO,YYDOTS(K+1,I,J))
              YYDOTS(K+1,I,J) = (ONE-RFCT)*YYDOTS(K+1,I,J)+RFCT*ZYDOTS_MONO

           END DO
        END DO                                                                                     
     END IF
     
     !************************************************************************************!
     !* CORRECTED UPWIND SCHEME *!
     !***************************!
     DO I=1,GXPT
        DO J=1,GYPT
           YPIS(K+1,I,J) = YPIS(K,I,J)                                              &
                         & - (RDT/RDX)*(YXDOTS(K+1,I,J)-YXDOTS(K+1,NPERIO(I-1),J) ) &
                         & - (RDT/RDY)*(YYDOTS(K+1,I,J)-YYDOTS(K+1,I,NPERIO(J-1)) )
        END DO
     END DO
     
  END DO
  
  PPIS_DEP(:,:) = YPIS(NORD2,:,:) + ZPIS_MIN 

  
END SUBROUTINE MPDATA_SURF_PRESSURE
!################################################################################################!
!************************************************************************************************!
!************************************************************************************************!
!************************************************************************************************! 
!************************************************************************************************!
!************************************************************************************************!
!************************************************************************************************!
!************************************************************************************************!
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
