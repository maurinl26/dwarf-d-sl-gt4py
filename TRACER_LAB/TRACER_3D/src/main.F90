!##################   PRINCIPALE   ###################!
!#                                                   #!
!# auteur : F.Voitus                                 #!
!#                                                   #!
!#####################################################!

!=====================================================!
!================  MODULE PRINCIPAL  =================!
!=====================================================!
PROGRAM MAIN

  USE MPI
  USE MOD_PARAMETER
  USE MOD_STRUCTURE
  USE MOD_LINALG
  USE MOD_FFT
  USE MOD_DIAGNOSTIC
  USE MOD_VIEW
  USE MOD_SLAG
  USE MOD_STEPO

  IMPLICIT NONE

  INTEGER                              :: CODE,MYPROC
  INTEGER                              :: TEMPS,JOUR
  INTEGER(8)                           :: DIM

  TYPE(GEOMETRY)                       :: GEO
  TYPE(DOMAIN)                         :: DOM
  TYPE(BACKGROUND)                     :: STB

  TYPE(VARPROG)                        :: ST0
  TYPE(VARGLOB)                        :: STG
  TYPE(SMILAG)                         :: YSLAG
  TYPE(MPDATA)                         :: YMDTA
  TYPE(RHS_LOC)                        :: RHSL_A,RHSL_D,RHSL_M
  TYPE(RHS_GLB)                        :: RHSG_A,RHSG_D

  INTEGER(8),DIMENSION(:),ALLOCATABLE  :: ABS_GLOBAL,ORD_GLOBAL
  INTEGER(8)                           :: NSTEP,NITER
  REAL(8)                              :: PDT

  
  CALL MPI_INIT(CODE)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NB_PROCS,CODE)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYPROC,CODE)  
  CALL DEF_MPI_DIMENSIONING(ABS_GLOBAL,ORD_GLOBAL,NB_PROCS,MYPROC)
  
  CALL DEF_CONFIG_PARAMETER
  CALL SUALLOC_MPI(ST0,RHSL_M,RHSL_A,RHSL_D,GEO,NDIM)

  ! Initialisation 
  CALL MODEL_INITIALIZATION(ST0,STB,STG,DOM,GEO, &
       & ABS_GLOBAL,ORD_GLOBAL,NDIM,MYPROC,CODE)
  
  IF (MYPROC .EQ. 0) THEN
     WRITE(*,*) "Integration... \n"
     WRITE(*,*)
     CALL SUALLOC_GLB(YSLAG,YMDTA,RHSG_A,RHSG_D)
     CALL SCREEN_START(TEMPS,JOUR,DOM,STB)
  END IF
  
  CALL MPI_BARRIER(MPI_COMM_WORLD,CODE) 

  ! Dynamic Loop
  DO NSTEP = 1, NTIMESTEP
     
    CALL RHS_MDL(RHSL_M,ST0,GEO)
    CALL RHS_EXP(RHSL_A,ST0,RHSL_M,GEO,PDT/TWO) 

    DO NITER = 0,NSITER

       IF (LSLAG) THEN
          CALL SMILAG_TRANSPORT_SCHEME(RHSL_D,RHSL_A,STG,RHSG_D,RHSG_A,YSLAG, &
               & ST0%U,ST0%V,ST0%ETADOTF,GEO%ETAF,GEO%ETAH,GEO%DELTB,NITER,   &
               & NSTEP,NDIM,ABS_GLOBAL,ORD_GLOBAL,MYPROC,CODE,C_ADVECT)
       ELSE 
          CALL MPDATA_TRANSPORT_SCHEME(RHSL_D,RHSL_A,STG,RHSG_D,RHSG_A,YMDTA, &
               & ST%MF,ST%MH,ST0%U,ST0%V,ST0%ETADOTH,GEO%ETAF,GEO%DETAF,      &
               & GEO%ETAH,GEO%DETAH,GEO%DELTB,GEO%EPS,NITER,NSTEP,NDIM,       &
               & ABS_GLOBAL,ORD_GLOBAL,MYPROC,CODE,C_ADVECT)
       END IF   
       
       CALL RHS_MDL(RHSL_M,ST0,GEO)
       CALL RHS_IMP(ST0,RHSL_M,RHSL_D,GEO,STB,PDT/TWO)
       CALL SPC_DER(ST0,STG,STB,DOM,ABS_GLOBAL,ORD_GLOBAL,MYPROC,CODE,NSTEP)
       
       CALL UPDATE_METRIC_TERMS(ST0,GEO)

    END DO
       
    IF (MYPROC .EQ. 0) THEN
       CALL SCREEN(TEMPS,RHSG,DOM,NSTEP)
    END IF

    
  END DO
 
  !CALL MPI_BARRIER(MPI_COMM_WORLD,CODE)
  !THIS MPI_BARRIER MIGHT BE USELESS
  
  CALL FREEMEM_MPI(ST0,RHS_M,RHS_A,RHS_D,GEO)
  
  IF (MYPROC .EQ. 0) THEN
     CALL SCREEN_END(TEMPS,JOUR)
     CALL FREEMEM_GBL(STG,YSLAG,YMDTA,DOM)
     DEALLOCATE(ABS_GLOBAL,ORD_GLOBAL)
  END IF
  
  CALL MPI_FINALIZE(CODE)

END PROGRAM MAIN

!=====================================================!
!=====================================================!
!=====================================================!
