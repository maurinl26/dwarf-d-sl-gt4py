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
  
  USE MPI
  USE OMP_LIB
  USE MOD_STRUCTURE, ONLY : SMILAG,GEOMETRY,RHS_GLB,RHS_LOC
  USE MOD_PARAMETER, ONLY : GXPT,GYPT,NLEV,RR,RG,RCP,RCV,RPI,ONE,TWO,ZERO,HALF,RDT, &
                          & RDX,RDY,NSLTER,NPERIOX,NPERIOY,LPERIO,LPRINTLEV,NB_PROCS,NDIM
  USE MOD_COMM,  ONLY     : TRANSGLB_RHS,WIND_GLB_TRANSFER,TRANSLOC_RHS
 
CONTAINS

  !********************************************************************************************!
  !********************************************************************************************!
  !********************************************************************************************!
  !*******************************  LISTE DES ROUTINES ****************************************!
  !********************************************************************************************!  
  !********************************************************************************************!
  !********************************************************************************************!
  SUBROUTINE SMILAG_TRANSPORT_SCHEME(RHSL_DEP,RHSL_ARR,RHSG_DEP,RHSG_ARR,   &
             & SLSTRUC,PXDOT0,PYDOT0,PZDOT0,PETAF,PETAH,PDELTB,KITER,KSTEP, &
             & KDIM,ABS_GLB,ORD_GLB,MYPROC,CODE,CD_ADV)
    
  IMPLICIT NONE

  TYPE(SMILAG),                        INTENT(INOUT)  :: SLSTRUC
  TYPE(RHS_LOC),                       INTENT(INOUT)  :: RHSL_DEP
  TYPE(RHS_LOC),                       INTENT(IN)     :: RHSL_ARR
  TYPE(RHS_GLB),                       INTENT(INOUT)  :: RHSG_ARR
  TYPE(RHS_GLB),                       INTENT(INOUT)  :: RHSG_DEP
  REAL(8),                             INTENT(IN)     :: PXDOT0(NLEV,KDIM)
  REAL(8),                             INTENT(IN)     :: PYDOT0(NLEV,KDIM)
  REAL(8),                             INTENT(IN)     :: PZDOT0(NLEV,KDIM)
  REAL(8),                             INTENT(IN)     :: PETAF(NLEV)
  REAL(8),                             INTENT(IN)     :: PETAH(0:NLEV)
  REAL(8),                             INTENT(IN)     :: PDELTB(NLEV)
  INTEGER(8),                          INTENT(IN)     :: KITER,KSTEP,KDIM
  INTEGER,                             INTENT(INOUT)  :: CODE
  INTEGER,                             INTENT(IN)     :: MYPROC
  INTEGER(8), DIMENSION(0:NB_PROCS-1), INTENT(IN)     :: ABS_GLB,ORD_GLB
  CHARACTER(LEN=1),                    INTENT(IN)     :: CD_ADV

  INTEGER(8)                                          :: I,J,K
  REAL(8)                                             :: ZLIN,ZWIND_HOR_MAX,ZWIND_VER_MAX

  ZLIN = ONE
  IF (CD_ADV == 'LI') ZLIN =ZERO
  
  IF (KITER == 0) THEN
     
     CALL TRANSGLB_RHS(RHSG_ARR,RHSL_ARR,ABS_GLB,ORD_GLB,MYPROC,CODE,1)
     CALL WIND_GLB_TRANSFER(SLSTRUC%XDOT,SLSTRUC%YDOT,SLSTRUC%ZDOT,    &
          & PXDOT0,PYDOT0,PZDOT0,1,ABS_GLB,ORD_GLB,MYPROC,CODE,10)

  ELSE

     CALL WIND_GLB_TRANSFER(SLSTRUC%XDOTF,SLSTRUC%YDOTF,SLSTRUC%ZDOTF, &
          & PXDOT0,PYDOT0,PZDOT0,1,ABS_GLB,ORD_GLB,MYPROC,CODE,10)     
  ENDIF

  
  IF (MYPROC == 0) THEN

     IF (KITER == 0) THEN
        SLSTRUC%XDOTF(:,:,:) = SLSTRUC%XDOT(:,:,:)
        SLSTRUC%YDOTF(:,:,:) = SLSTRUC%YDOT(:,:,:)
        SLSTRUC%ZDOTF(:,:,:) = SLSTRUC%ZDOT(:,:,:)
     END IF   
     
     !**************************************************************!
     ! STABILITY CHECK-POINT MAXIMUM OF HORIZONTAL WIND VELOCITY    !
     !**************************************************************!
     ZWIND_HOR_MAX = MAX(MAXVAL(SLSTRUC%XDOT),MAXVAL(SLSTRUC%YDOT)) !
     ZWIND_VER_MAX = MAXVAL(SLSTRUC%ZDOT)                           !
                                                                    !
     IF (LPRINTLEV) THEN                                            !
        IF (KITER == 0) WRITE(*,*) 'NSTEP : ',KSTEP,         &      !
                      & '--> MAX.WIND_HOR : ',ZWIND_HOR_MAX, &      !
                      & '--> MAX.WIND_VER : ',ZWIND_VER_MAX         !
     END IF                                                         !
     IF (ZWIND_HOR_MAX .GE. 300.d0) THEN                            !
        WRITE(*,*) 'WIND TOO STRONG EXPLOSION AT TIMESTEP = ',KSTEP !
        STOP                                                        !
     END IF                                                         !
     !**************************************************************!

     
     !********************************************************!
     ! RESEARCH OF DEPARTURE POINTS AND INTERPOLATION WEIGHTS !   
     !********************************************************!
     
     CALL LARCINAF(SLSTRUC,PETAF)
     CALL LARCINAH(SLSTRUC,PETAH,PETAF,PDELTB)

     !********************************************************!
     ! SL INTERPOLATIONS AT DEPARTURE POINTS                  !
     !********************************************************!
     IF (CD_ADV == 'LI') THEN 
        RHSG_DEP%T(:,:,:) = RHSG_ARR%T(:,:,:)
        RHSG_DEP%Q(:,:,:) = RHSG_ARR%Q(:,:,:)
        RHSG_DEP%PIS(:,:) = RHSG_ARR%PIS(:,:)
     ELSE 
        CALL LARCINB(RHSG_DEP%T,RHSG_ARR%T,RHSG_DEP%Q,RHSG_ARR%Q,         &
             & RHSG_DEP%PIS,RHSG_ARR%PIS,SLSTRUC%XWEI,SLSTRUC%NXLAG,      &
             & SLSTRUC%YWEI,SLSTRUC%NYLAG,SLSTRUC%ZWEI,SLSTRUC%NZLAG,     &
             & SLSTRUC%XWEIH,SLSTRUC%NXLAGH,SLSTRUC%YWEIH,SLSTRUC%NYLAGH, &
             & SLSTRUC%ZWEIH,SLSTRUC%NZLAGH,SLSTRUC%XWEIS,SLSTRUC%NXLAGS, &
             & SLSTRUC%YWEIS,SLSTRUC%NYLAGS)
     END IF
     
  END IF

  !*************************************!
  !waiting for the end of myproc=0 task !
  CALL MPI_BARRIER(MPI_COMM_WORLD,CODE) !
  !*************************************!
   
  CALL TRANSLOC_RHS(RHSL_DEP,RHSG_DEP,ABS_GLB,ORD_GLB,MYPROC,CODE,1)

  END SUBROUTINE SMILAG_TRANSPORT_SCHEME
  !################################################################################################
  !################################################################################################
  !################################################################################################
  !################################################################################################  
  SUBROUTINE LARCINAF(SLSTRUC,PETA)  

    IMPLICIT NONE
               
    TYPE(SL_STRUCT),                    INTENT(INOUT) :: SLSTRUC
    REAL(8), DIMENSION(NLEV),           INTENT(IN)    :: PETA

    REAL(8), DIMENSION(NLEV,GXPT,GYPT)                :: XDOT_DEP,YDOT_DEP,ZDOT_DEP
    INTEGER(8)                                        :: JITER

    XDOT_DEP(:,:,:) = SLSTRUC%XDOT(:,:,:)
    YDOT_DEP(:,:,:) = SLSTRUC%YDOT(:,:,:)
    ZDOT_DEP(:,:,:) = SLSTRUC%ZDOT(:,:,:) 

    !Compute first guess SL trajectories
    CALL ELARCHE(SLSTRUC%XTRAJ,SLSTRUC%YTRAJ,SLSTRUC%ZTRAJ,    &
         & SLSTRUC%XDOTF,SLSTRUC%YDOTF,SLSTRUC%ZDOTF,          &
         & XDOT_DEP,YDOT_DEP,ZDOT_DEP,PETA,1_8)

    DO JITER=1,NSLTER
       !Compute SL weights
       CALL ELASCAW(SLSTRUC%XTRAJ,SLSTRUC%YTRAJ,SLSTRUC%ZTRAJ, &
            & SLSTRUC%NXLAG,SLSTRUC%NYLAG,SLSTRUC%NZLAG,       &
            & SLSTRUC%XWEI,SLSTRUC%YWEI,SLSTRUC%ZWEI,PETA,1_8)
       !Interpolations velocity components to Departure points
       CALL ELARMES(XDOT_DEP,YDOT_DEP,ZDOT_DEP,SLSTRUC%XDOT,   &
            & SLSTRUC%YDOT,SLSTRUC%ZDOT,SLSTRUC%NXLAG,         &
            & SLSTRUC%NYLAG,SLSTRUC%NZLAG,SLSTRUC%XWEI,        &
            & SLSTRUC%YWEI,SLSTRUC%ZWEI,1_8)
       !Compute SL trajectories
       CALL ELARCHE(SLSTRUC%XTRAJ,SLSTRUC%YTRAJ,SLSTRUC%ZTRAJ, &
            & SLSTRUC%XDOTF,SLSTRUC%YDOTF,SLSTRUC%ZDOTF,       &
            & XDOT_DEP,YDOT_DEP,ZDOT_DEP,PETA,1_8)
    END DO

    !Compute new SL weights for Sources
    CALL ELASCAW(SLSTRUC%XTRAJ,SLSTRUC%YTRAJ,SLSTRUC%ZTRAJ,    &
         & SLSTRUC%NXLAG,SLSTRUC%NYLAG,SLSTRUC%NZLAG,          &
         & SLSTRUC%XWEI,SLSTRUC%YWEI,SLSTRUC%ZWEI,PETA,1_8)

  END SUBROUTINE LARCINAF
  !*******************************************************************************!
  !*******************************************************************************!
  SUBROUTINE LARCINA_OTH(SLSTRUC,PETAH,PETAF,PDELTB)

    IMPLICIT NONE

    TYPE(SL_STRUCT),   INTENT(INOUT)  :: SLSTRUC
    REAL(8),           INTENT(IN)     :: PETAH(0:NLEV)
    REAL(8),           INTENT(IN)     :: PETAF(NLEV),PDELTB(NLEV)

    INTEGER(8)                        :: II ,JJ, KK
    REAL(8)                           :: ZETAH(0:NLEV)
    REAL(8)                           :: ZXTRAJH(0:NLEV,GXPT,GYPT)
    REAL(8)                           :: ZYTRAJH(0:NLEV,GXPT,GYPT)
    REAL(8)                           :: ZZTRAJH(0:NLEV,GXPT,GYPT)
    REAL(8)                           :: ZXTRAJS(GXPT,GYPT)
    REAL(8)                           :: ZYTRAJS(GXPT,GYPT)
    REAL(8)                           :: ZTRAJZ,ZWEI

    DO II=1,GXPT
       DO JJ=1,GYPT
          
          ZXTRAJH(0,II,JJ)     = SLSTRUC%XTRAJ(1,II,JJ)
          ZYTRAJH(0,II,JJ)     = SLSTRUC%YTRAJ(1,II,JJ)
          ZZTRAJH(0,II,JJ)     = PETAH(0)
          ZETAH(0)             = PETAH(0)

          DO KK=1,NLEV-1
             ZWEI              = ( (PETAH(KK)-PETAF(KK)) &
                               & /(PETAF(KK+1)-PETAF(KK)) )
             ZETAH(KK)         = PETAH(KK)  !PETAF(KK)   &
                               !& + ZWEI*(PETAF(KK+1)-PETAF(KK)) 
             ZXTRAJH(KK,II,JJ) = SLSTRUC%XTRAJ(KK,II,JJ) &
                               & + ZWEI*( SLSTRUC%XTRAJ(KK+1,II,JJ) &
                               & - SLSTRUC%XTRAJ(KK,II,JJ) )
             ZYTRAJH(KK,II,JJ) = SLSTRUC%YTRAJ(KK,II,JJ) &
                               & + ZWEI*( SLSTRUC%YTRAJ(KK+1,II,JJ) &
                               & - SLSTRUC%YTRAJ(KK,II,JJ) )
             ZTRAJZ            = SLSTRUC%ZTRAJ(KK,II,JJ) &
                               & + ZWEI*( SLSTRUC%ZTRAJ(KK+1,II,JJ) &
                               & - SLSTRUC%ZTRAJ(KK,II,JJ) )
             ZZTRAJH(KK,II,JJ) = MAX(PETAH(0),MIN(PETAH(NLEV),ZTRAJZ))
          ENDDO

          ZXTRAJH(NLEV,II,JJ)  = SLSTRUC%XTRAJ(NLEV,II,JJ)
          ZYTRAJH(NLEV,II,JJ)  = SLSTRUC%YTRAJ(NLEV,II,JJ)
          ZZTRAJH(NLEV,II,JJ)  = PETAH(NLEV)
          ZETAH(NLEV)          = PETAH(NLEV)


          ZXTRAJS(II,JJ)       = ZERO
          ZYTRAJS(II,JJ)       = ZERO
          DO KK=1,NLEV
             ZXTRAJS(II,JJ)    = ZXTRAJS(II,JJ) &
                               & + PDELTB(KK)*SLSTRUC%XTRAJ(KK,II,JJ)
             ZYTRAJS(II,JJ)    = ZYTRAJS(II,JJ) &
                               & + PDELTB(KK)*SLSTRUC%YTRAJ(KK,II,JJ)
          END DO
          
       ENDDO
    ENDDO
        
    CALL ELASCAW(ZXTRAJH,ZYTRAJH,ZZTRAJH,SLSTRUC%NXLAGH, &
         & SLSTRUC%NYLAGH,SLSTRUC%NZLAGH,SLSTRUC%XWEIH,  &
         & SLSTRUC%YWEIH,SLSTRUC%ZWEIH,ZETAH,0_8)

    CALL ELASCAW_SURF(ZXTRAJS,ZYTRAJS,SLSTRUC%NXLAGS,    &
         & SLSTRUC%NYLAGS,SLSTRUC%XWEIS,SLSTRUC%YWEIS)
    
  END SUBROUTINE LARCINA_OTH
  !*******************************************************************************!
  !*******************************************************************************!
  !*******************************************************************************!
  !*******************************************************************************!
  !*******************************************************************************!
  !*******************************************************************************!
  !*******************************************************************************!
  SUBROUTINE ELARCHE(PXTRAJ,PYTRAJ,PZTRAJ,PXDOT_ARR,PYDOT_ARR,PZDOT_ARR, &
                    & PXDOT_DEP,PYDOT_DEP,PZDOT_DEP,PETA,KTOP)
    
    IMPLICIT NONE

    REAL(8),DIMENSION(KTOP:NLEV,GXPT,GYPT), INTENT(OUT) :: PXTRAJ,PYTRAJ,PZTRAJ      
    REAL(8),DIMENSION(KTOP:NLEV,GXPT,GYPT), INTENT(IN)  :: PXDOT_ARR,PYDOT_ARR,PZDOT_ARR
    REAL(8),DIMENSION(KTOP:NLEV,GXPT,GYPT), INTENT(IN)  :: PXDOT_DEP,PYDOT_DEP,PZDOT_DEP
    REAL(8),DIMENSION(KTOP:NLEV),           INTENT(IN)  :: PETA
    INTEGER(8),                             INTENT(IN)  :: KTOP

    INTEGER(8)                                          :: II,  JJ,  KK
    REAL(8)                                             :: ZVETAOX, ZVETAON, ZTRAJ_INIT

    ZVETAON=PETA(KTOP)
    ZVETAOX=PETA(NLEV)

    ! COMPUTE SL TRAJECTORIES AT DEPARTURE POINT
    DO II=1,GXPT
       DO JJ=1,GYPT   
          DO KK=KTOP,NLEV
             PXTRAJ(KK,II,JJ) = REAL(II-1,8)*RDX    &
                              & -(RDT/TWO)*(PXDOT_ARR(KK,II,JJ) &
                              & + PXDOT_DEP(KK,II,JJ))
             PYTRAJ(KK,II,JJ) = REAL(JJ-1,8)*RDY    &
                              & -(RDT/TWO)*(PYDOT_ARR(KK,II,JJ) &
                              & + PYDOT_DEP(KK,II,JJ))
             ZTRAJ_INIT       = PETA(JJ)            &
                              & -(RDT/TWO)*(PZDOT_ARR(KK,II,JJ) &
                              & + PZDOT_DEP(KK,II,JJ))
             PZTRAJ(KK,II,JJ) = MAX(ZVETAON,MIN(ZVETAOX,ZTRAJ_INIT))
          ENDDO
       ENDDO
    ENDDO

END SUBROUTINE ELARCHE
!************************************************************************************!
!************************************************************************************!
!************************************************************************************!
!************************************************************************************!
!************************************************************************************!
SUBROUTINE ELASCAW(PXTRAJ,PYTRAJ,PZTRAJ,KXLAG,KYLAG,KZLAG,PXWEI,PYWEI,PZWEI,PETA,KTOP)

IMPLICIT NONE

REAL(8),   DIMENSION(4,KTOP:NLEV,GXPT,GYPT), INTENT(OUT) :: PXWEI,PYWEI,PZWEI
INTEGER(8),DIMENSION(4,KTOP:NLEV,GXPT,GYPT), INTENT(OUT) :: KXLAG,KYLAG,KZLAG
REAL(8),   DIMENSION(KTOP:NLEV,GXPT,GYPT),   INTENT(IN)  :: PXTRAJ,PYTRAJ,PZTRAJ
REAL(8),   DIMENSION(KTOP:NLEV),             INTENT(IN)  :: PETA
INTEGER(8),                                  INTENT(IN)  :: KTOP

INTEGER(8)    :: II   ,JJ   ,KK   ,LL   ,IP   ,JP ,KP ,K0            
REAL(8)       :: ZXP1 ,ZX0  ,ZXM1 ,ZXM2 ,ZXD 
REAL(8)       :: ZYP1 ,ZY0  ,ZYM1 ,ZYM2 ,ZYD    
REAL(8)       :: ZZP1 ,ZZ0  ,ZZM1 ,ZZM2 ,ZZD      
REAL(8)       :: ZEPS ,ZALX ,ZALY ,ZALZ ,ZSIGN ,ZDIFF

ZEPS = ABS(EPSILON(ONE))

!****************************************!
! COMPUTE WEIGHTS FOR SL INTERPOLATIONS *!
!****************************************!
       
  DO II=1,GXPT   
     DO JJ=1,GYPT
        DO KK=KTOP,NLEV
         
           !**************
           ! X-DIRECTION !
           !**************
           
           ZXD   = PXTRAJ(KK,II,JJ) 
           ZALX  = REAL(II-1,8)-(ZXD/RDX)
           IF (ABS(ZALX).LE.ZEPS) ZALX = ZERO
           !IP = INT(ZAL+MIN(ZERO,SIGN(ONE,ZAL)))
           IP = FLOOR(ZALX)
           DO LL=1,4
              KXLAG(LL,KK,II,JJ) = NPERIOX(II-IP+2-LL) 
           END DO
           ! cubic interpolation between II-IP-2, II-IP-1, II-IP and II-IP+1
           ZXP1=REAL(II-IP+1-1,8)*RDX
           ZX0 =REAL(II-IP+0-1,8)*RDX
           ZXM1=REAL(II-IP-1-1,8)*RDX
           ZXM2=REAL(II-IP-2-1,8)*RDX
           PXWEI(1,KK,II,JJ)=((ZXD-ZX0)/(ZXP1-ZX0))*((ZXD-ZXM1)/(ZXP1-ZXM1))*((ZXD-ZXM2)/(ZXP1-ZXM2))
           PXWEI(2,KK,II,JJ)=((ZXD-ZXP1)/(ZX0-ZXP1))*((ZXD-ZXM1)/(ZX0-ZXM1))*((ZXD-ZXM2)/(ZX0-ZXM2))
           PXWEI(3,KK,II,JJ)=((ZXD-ZXP1)/(ZXM1-ZXP1))*((ZXD-ZX0)/(ZXM1-ZX0))*((ZXD-ZXM2)/(ZXM1-ZXM2))
           PXWEI(4,KK,II,JJ)=((ZXD-ZXP1)/(ZXM2-ZXP1))*((ZXD-ZX0)/(ZXM2-ZX0))*((ZXD-ZXM1)/(ZXM2-ZXM1))

           !**************
           ! Y-DIRECTION !
           !**************

           ZYD   = PYTRAJ(KK,II,JJ) 
           ZALY  = REAL(JJ-1,8)-(ZYD/RDY)
           IF (ABS(ZALY).LE.ZEPS) ZALY = ZERO
           !JP = INT(ZALY+MIN(ZERO,SIGN(ONE,ZALY)))
           JP = FLOOR(ZALY)
           DO LL=1,4
              KYLAG(LL,KK,II,JJ) = NPERIOY(JJ-JP+2-LL) 
           END DO
           ! cubic interpolation between JJ-JP-2, JJ-JP-1, JJ-JP and JJ-JP+1
           ZYP1=REAL(JJ-JP+1-1,8)*RDY
           ZY0 =REAL(JJ-JP+0-1,8)*RDY
           ZYM1=REAL(JJ-JP-1-1,8)*RDY
           ZYM2=REAL(JJ-JP-2-1,8)*RDY
           PYWEI(1,KK,II,JJ)=((ZYD-ZY0)/(ZYP1-ZY0))*((ZYD-ZYM1)/(ZYP1-ZYM1))*((ZYD-ZYM2)/(ZYP1-ZYM2))
           PYWEI(2,KK,II,JJ)=((ZYD-ZYP1)/(ZY0-ZYP1))*((ZYD-ZYM1)/(ZY0-ZYM1))*((ZYD-ZYM2)/(ZY0-ZYM2))
           PYWEI(3,KK,II,JJ)=((ZYD-ZYP1)/(ZYM1-ZYP1))*((ZYD-ZY0)/(ZYM1-ZY0))*((ZYD-ZYM2)/(ZYM1-ZYM2))
           PYWEI(4,KK,II,JJ)=((ZYD-ZYP1)/(ZYM2-ZYP1))*((ZYD-ZY0)/(ZYM2-ZY0))*((ZYD-ZYM1)/(ZYM2-ZYM1))

           !*************!
           ! Z-DIRECTION !
           !*************!

           ZZD   = PZTRAJ(KK,II,JJ) 
           ZALZ  = PETA(KK)-PZTRAJ(KK,II,JJ)
           IF (ABS(ZALZ).LT.ZEPS) THEN
              ZALZ = ZERO
              KP   = KK
           ELSE
              KP = KK
              ZSIGN = MAX(ZERO,SIGN(ONE,ZALZ))
              IF (ZSIGN == ABS(ONE)) THEN
                 DO K0=KK,KTOP+1,-1
                    ZDIFF = MAX(ZERO,SIGN(ONE,(ZZD-PETA(K0))*(ZZD-PETA(K0-1))))
                    IF (ZDIFF==ZERO) THEN
                      KP=K0
                      EXIT
                    END IF
                    KP=KTOP
                 END DO
              ELSE
                KP = KK
                DO K0=KK+1,NLEV
                   ZDIFF = MAX(ZERO,SIGN(ONE,(ZZD-PETA(K0))*(ZZD-PETA(K0-1))))
                   IF (ZDIFF==ZERO) THEN
                      KP=K0
                      EXIT   
                   END IF
                   KP = NLEV
                END DO  
              END IF
           END IF
           DO LL=1,4
              KZLAG(LL,KK,II,JJ) = MAX(KTOP,MIN(NLEV,KP+2-LL)) 
           END DO
           ! cubic interpolation between KK-KP-2, KK-KP-1, KK-KP and KK-KP+1
           IF ((KP+1.LE.NLEV).AND.(KP-2.GE.KTOP)) THEN
             ZZP1=PETA(KP+1)
             ZZ0 =PETA(KP+0)
             ZZM1=PETA(KP-1)
             ZZM2=PETA(KP-2)
             PZWEI(1,KK,II,JJ)=((ZZD-ZZ0)/(ZZP1-ZZ0))*((ZZD-ZZM1)/(ZZP1-ZZM1))*((ZZD-ZZM2)/(ZZP1-ZZM2))
             PZWEI(2,KK,II,JJ)=((ZZD-ZZP1)/(ZZ0-ZZP1))*((ZZD-ZZM1)/(ZZ0-ZZM1))*((ZZD-ZZM2)/(ZZ0-ZZM2))
             PZWEI(3,KK,II,JJ)=((ZZD-ZZP1)/(ZZM1-ZZP1))*((ZZD-ZZ0)/(ZZM1-ZZ0))*((ZZD-ZZM2)/(ZZM1-ZZM2))
             PZWEI(4,KK,II,JJ)=((ZZD-ZZP1)/(ZZM2-ZZP1))*((ZZD-ZZ0)/(ZZM2-ZZ0))*((ZZD-ZZM1)/(ZZM2-ZZM1)) 
           ELSEIF (KP-1.EQ.KTOP) THEN
           !  ! quadratic interpolation between KK-KP-1, KK-KP and KK-KP+1 at x=0
             ZZP1=PETA(KP+1)
             ZZ0 =PETA(KP+0)
             ZZM1=PETA(KP-1)
             PZWEI(1,KK,II,JJ)=((ZZD-ZZ0)/(ZZP1-ZZ0))*((ZZD-ZZM1)/(ZZP1-ZZM1))
             PZWEI(2,KK,II,JJ)=((ZZD-ZZP1)/(ZZ0-ZZP1))*((ZZD-ZZM1)/(ZZ0-ZZM1))
             PZWEI(3,KK,II,JJ)=((ZZD-ZZP1)/(ZZM1-ZZP1))*((ZZD-ZZ0)/(ZZM1-ZZ0))
             PZWEI(4,KK,II,JJ)=ZERO
           ELSEIF (KP.LE.KTOP) THEN
             ! Trajectory truncation at x=L
             PZWEI(1,KK,II,JJ)=ZERO
             PZWEI(2,KK,II,JJ)=ONE
             PZWEI(3,KK,II,JJ)=ZERO
             PZWEI(4,KK,II,JJ)=ZERO    
           ELSEIF (KP.EQ.NLEV) THEN
             ! quadratic interpolation between JJ-JP-2, JJ-JP-1 and JJ-JP at x=L 
             ZZ0 =PETA(KP+0)
             ZZM1=PETA(KP-1)
             ZZM2=PETA(KP-2)
             PZWEI(1,KK,II,JJ)=ZERO
             PZWEI(2,KK,II,JJ)=((ZZD-ZZM1)/(ZZ0-ZZM1))*((ZZD-ZZM2)/(ZZ0-ZZM2))
             PZWEI(3,KK,II,JJ)=((ZZD-ZZ0)/(ZZM1-ZZ0))*((ZZD-ZZM2)/(ZZM1-ZZM2))
             PZWEI(4,KK,II,JJ)=((ZZD-ZZ0)/(ZZM2-ZZ0))*((ZZD-ZZM1)/(ZZM2-ZZM1))  
           ELSEIF (KP-1.GE.NLEV) THEN
             ! Trajectory truncation at x=L
             PZWEI(1,KK,II,JJ)=ZERO
             PZWEI(2,KK,II,JJ)=ZERO
             PZWEI(3,KK,II,JJ)=ONE
             PZWEI(4,KK,II,JJ)=ZERO
           ENDIF
          
        END DO  
     END DO
  END DO

END SUBROUTINE ELASCAW
!****************************************************************************************!
!****************************************************************************************!
!****************************************************************************************!
!****************************************************************************************!
!****************************************************************************************!
!****************************************************************************************!
!****************************************************************************************!
SUBROUTINE ELASCAW_SURF(PXTRAJ,PYTRAJ,KXLAG,KYLAG,PXWEI,PYWEI)

IMPLICIT NONE

REAL(8),   DIMENSION(4,GXPT,GYPT), INTENT(OUT) :: PXWEI,PYWEI,PZWEI
INTEGER(8),DIMENSION(4,GXPT,GYPT), INTENT(OUT) :: KXLAG,KYLAG,KZLAG
REAL(8),   DIMENSION(GXPT,GYPT),   INTENT(IN)  :: PXTRAJ,PYTRAJ,PZTRAJ

INTEGER(8)                     :: II       ,JJ       ,LL
INTEGER(8)                     :: IP       ,JP             
REAL(8)                        :: ZXP1     ,ZX0      ,ZXM1     ,ZXM2     ,ZXD 
REAL(8)                        :: ZYP1     ,ZY0      ,ZYM1     ,ZYM2     ,ZYD      
REAL(8)                        :: ZEPS     ,ZALX     ,ZALY

ZEPS = ABS(EPSILON(ONE))

!****************************************!
! COMPUTE WEIGHTS FOR SL INTERPOLATIONS *!
!****************************************!
       
  DO II=1,GXPT   
     DO JJ=1,GYPT
     
           !**************
           ! X-DIRECTION !
           !**************
           ZXD   = PXTRAJ(II,JJ) 
           ZALX  = REAL(II-1,8)-(ZXD/RDX)
           IF (ABS(ZALX).LT.ZEPS) ZALX = ZERO
           !IP = INT(ZAL+MIN(ZERO,SIGN(ONE,ZAL)))
           IP = FLOOR(ZALX)
           DO LL=1,4
              KXLAG(LL,II,JJ) = NPERIOX(II-IP+2-LL) 
           END DO
           ! cubic interpolation between II-IP-2, II-IP-1, II-IP and II-IP+1
           ZXP1=REAL(II-IP+1-1,8)*RDX
           ZX0 =REAL(II-IP+0-1,8)*RDX
           ZXM1=REAL(II-IP-1-1,8)*RDX
           ZXM2=REAL(II-IP-2-1,8)*RDX
           PXWEI(1,II,JJ)=((ZXD-ZX0)/(ZXP1-ZX0))*((ZXD-ZXM1)/(ZXP1-ZXM1))*((ZXD-ZXM2)/(ZXP1-ZXM2))
           PXWEI(2,II,JJ)=((ZXD-ZXP1)/(ZX0-ZXP1))*((ZXD-ZXM1)/(ZX0-ZXM1))*((ZXD-ZXM2)/(ZX0-ZXM2))
           PXWEI(3,II,JJ)=((ZXD-ZXP1)/(ZXM1-ZXP1))*((ZXD-ZX0)/(ZXM1-ZX0))*((ZXD-ZXM2)/(ZXM1-ZXM2))
           PXWEI(4,II,JJ)=((ZXD-ZXP1)/(ZXM2-ZXP1))*((ZXD-ZX0)/(ZXM2-ZX0))*((ZXD-ZXM1)/(ZXM2-ZXM1))

           !**************
           ! Y-DIRECTION !
           !**************
           ZYD   = PYTRAJ(II,JJ) 
           ZALY  = REAL(JJ-1,8)-(ZYD/RDY)
           IF (ABS(ZALY).LT.ZEPS) ZALY = ZERO
           !JP = INT(ZALY+MIN(ZERO,SIGN(ONE,ZALY)))
           JP = FLOOR(ZALY)
           DO LL=1,4
              KYLAG(LL,II,JJ) = NPERIOY(JJ-JP+2-LL) 
           END DO
           ! cubic interpolation between JJ-JP-2, JJ-JP-1, JJ-JP and JJ-JP+1
           ZYP1=REAL(JJ-JP+1-1,8)*RDY
           ZY0 =REAL(JJ-JP+0-1,8)*RDY
           ZYM1=REAL(JJ-JP-1-1,8)*RDY
           ZYM2=REAL(JJ-JP-2-1,8)*RDY
           PYWEI(1,II,JJ)=((ZYD-ZY0)/(ZYP1-ZY0))*((ZYD-ZYM1)/(ZYP1-ZYM1))*((ZYD-ZYM2)/(ZYP1-ZYM2))
           PYWEI(2,II,JJ)=((ZYD-ZYP1)/(ZY0-ZYP1))*((ZYD-ZYM1)/(ZY0-ZYM1))*((ZYD-ZYM2)/(ZY0-ZYM2))
           PYWEI(3,II,JJ)=((ZYD-ZYP1)/(ZYM1-ZYP1))*((ZYD-ZY0)/(ZYM1-ZY0))*((ZYD-ZYM2)/(ZYM1-ZYM2))
           PYWEI(4,II,JJ)=((ZYD-ZYP1)/(ZYM2-ZYP1))*((ZYD-ZY0)/(ZYM2-ZY0))*((ZYD-ZYM1)/(ZYM2-ZYM1))
           
     END DO
  END DO

END SUBROUTINE ELASCAW_SURF
!****************************************************************************************!
!****************************************************************************************!
!****************************************************************************************!
!****************************************************************************************!
SUBROUTINE ELARMES(PXDOT_DEP,PYDOT_DEP,PZDOT_DEP,PXDOT_ARR,PYDOT_ARR,PZDOT_ARR, &
                  & KXLAG,KYLAG,KZLAG,PXWEI,PYWEI,PZWEI,KTOP)

IMPLICIT NONE

REAL(8),DIMENSION(KTOP:NLEV,GXPT,GYPT),      INTENT(OUT) :: PXDOT_DEP ,PYDOT_DEP  ,PZDOT_DEP
REAL(8),DIMENSION(KTOP:NLEV,GXPT,GYPT),      INTENT(IN)  :: PXDOT_ARR ,PYDOT_ARR  ,PZDOT_ARR
REAL(8),DIMENSION(4,KTOP:NLEV,GXPT,GYPT),    INTENT(IN)  :: PXWEI     ,PYWEI      ,PZWEI
INTEGER(8),DIMENSION(4,KTOP:NLEV,GXPT,GYPT), INTENT(IN)  :: KXLAG     ,KYLAG      ,KZLAG
INTEGER(8),                                  INTENT(IN)  :: KTOP

INTEGER(8)                    :: KK  ,II  ,JJ
INTEGER                       :: K   ,I   ,J


DO II=1,GXPT
   DO JJ=1,GYPT
      DO KK=KTOP,NLEV
         PXDOT_DEP(KK,II,JJ) = ZERO
         PYDOT_DEP(KK,II,JJ) = ZERO
         PZDOT_DEP(KK,II,JJ) = ZERO
         DO K=1,4
            DO I=1,4
               DO J=1,4
                  PXDOT_DEP(KK,II,JJ) = PXDOT_DEP(KK,II,JJ) &
                       & + PZWEI(K,KK,II,JJ)*PXWEI(I,KK,II,JJ)*PYWEI(J,KK,II,JJ) &
                       & *PXDOT_ARR(KZLAG(K,KK,II,JJ),KXLAG(I,KK,II,JJ),KYLAG(J,KK,II,JJ))
                  PYDOT_DEP(KK,II,JJ) = PYDOT_DEP(KK,II,JJ) &
                       & + PZWEI(K,KK,II,JJ)*PXWEI(I,KK,II,JJ)*PYWEI(J,KK,II,JJ) &
                       & *PYDOT_ARR(KZLAG(K,KK,II,JJ),KXLAG(I,KK,II,JJ),KYLAG(J,KK,II,JJ))
                  PZDOT_DEP(KK,II,JJ) = PZDOT_DEP(KK,II,JJ) &
                       & + PZWEI(K,KK,II,JJ)*PXWEI(I,KK,II,JJ)*PYWEI(J,KK,II,JJ) &
                       & *PZDOT_ARR(KZLAG(K,KK,II,JJ),KXLAG(I,KK,II,JJ),KYLAG(J,KK,II,JJ))
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE ELARMES
!**********************************************************************************************!
!**********************************************************************************************!
!**********************************************************************************************!
!**********************************************************************************************!
!##############################################################################################!
!##############################################################################################!
!**********************************************************************************************!
!**********************************************************************************************!
!**********************************************************************************************!
!**********************************************************************************************!
SUBROUTINE LARCINB(RT_DEP,RT_ARR,RQ_DEP,RQ_ARR,RP_DEP,RP_ARR,PXWEI,KXLAG,PYWEI,KYLAG,PZWEI,KZLAG, &
           & PXWEIH,KXLAGH,PYWEIH,KYLAGH,PZWEIH,KZLAGH,PXWEIS,KXLAGS,PYWEIS,KYLAGS)
  
IMPLICIT NONE

REAL(8),          INTENT(OUT) :: RT_DEP(NLEV,GXPT,GYPT),RQ_DEP(0:NLEV,GXPT,GYPT),RP_DEP(GXPT,GYPT)
REAL(8),          INTENT(IN)  :: RT_ARR(NLEV,GXPT,GYPT),RQ_ARR(0:NLEV,GXPT,GYPT),RP_ARR(GXPT,GYPT)
REAL(8),          INTENT(IN)  :: PXWEI(4,NLEV,GXPT,GYPT),PZWEI(4,NLEV,GXPT,GYPT),PYWEI(4,NLEV,GXPT,GYPT)
REAL(8),          INTENT(IN)  :: PXWEIH(4,0:NLEV,GXPT,GYPT),PZWEIH(4,0:NLEV,GXPT,GYPT)
REAL(8),          INTENT(IN)  :: PXWEIS(4,GXPT,GYPT),PYWEIS(4,GXPT,GYPT),PYWEIH(4,0:NLEV,GXPT,GYPT)
INTEGER(8),       INTENT(IN)  :: KXLAG(4,NLEV,GXPT,GYPT),KZLAG(4,NLEV,GXPT,GYPT),KYLAG(4,NLEV,GXPT,GYPT)
INTEGER(8),       INTENT(IN)  :: KXLAGH(4,0:NLEV,GXPT,GYPT),KZLAGH(4,0:NLEV,GXPT,GYPT)
INTEGER(8),       INTENT(IN)  :: KXLAGS(4,GXPT,GYPT),KYLAGS(4,GXPT,GYPT),KYLAGH(4,0:NLEV,GXPT,GYPT)

INTEGER(8)                    :: KK ,II ,JJ
INTEGER                       :: K  ,I  ,J

DO II=1,GXPT
   DO JJ=1,GYPT
      
      DO KK=1,NLEV
         RT_DEP(KK,II,JJ) = ZERO
         RQ_DEP(KK,II,JJ) = ZERO
         DO I=1,4
            DO J=1,4
               DO K=1,4
                  RT_DEP(KK,II,JJ) = RT_DEP(KK,II,JJ) &
                    & + PZWEI(K,KK,II,JJ)*PXWEI(I,KK,II,JJ)*PYWEI(J,KK,II,JJ)* &
                    & RT_ARR(KZLAG(K,KK,II,JJ),KXLAG(I,KK,II,JJ),KYLAG(J,KK,II,JJ))
                  RQ_DEP(KK,II,JJ) = RQ_DEP(KK,II,JJ) &
                    & + PZWEIH(K,KK,II,JJ)*PXWEIH(I,KK,II,JJ)*PYWEIH(J,KK,II,JJ)* &
                    & RQ_ARR(KZLAGH(K,KK,II,JJ),KXLAGH(I,KK,II,JJ),KYLAGH(J,KK,II,JJ))
               ENDDO
            ENDDO   
         ENDDO
      ENDDO

      RP_DEP(II,JJ)   = ZERO
      RQ_DEP(0,II,JJ) = ZERO
      DO I=1,4
         DO J=1,4
            RP_DEP(II,JJ) = RP_DEP(II,JJ) &
                 & + PXWEIS(I,II,JJ)*PYWEIS(J,II,JJ)* &
                 & RP_ARR(KXLAGS(I,II,JJ),KYLAGS(J,II,JJ))
            DO K=1,4
               RQ_DEP(0,II,JJ) = RQ_DEP(0,II,JJ) &
                 & + PZWEIH(K,0,II,JJ)*PXWEIH(I,0,II,JJ)*PYWEIH(J,0,II,JJ)* &
                 & RQ_ARR(KZLAGH(K,0,II,JJ),KXLAGH(I,0,II,JJ),KYLAGH(J,0,II,JJ))
            ENDDO
         ENDDO   
      ENDDO
      
   ENDDO
ENDDO

END SUBROUTINE LARCINB 
!*******************************************************************************!
!*******************************************************************************!
!*******************************************************************************!
!*******************************************************************************!
!*******************************************************************************!
!*******************************************************************************!
!*******************************************************************************!
!*******************************************************************************!

END MODULE MOD_SLAG
!===============================================================================!
!===============================================================================!
!===============================================================================!
