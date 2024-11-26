!##################   PRINCIPALE   ###################!
!#                                                   #!
!# auteur : F.Voitus                                 #!
!#                                                   #!
!#####################################################!

!=====================================================!
!================  MODULE PRINCIPAL  =================!
!=====================================================!
PROGRAM MAIN

  USE MOD_PARAM
  USE MOD_SETUP
  USE MOD_FFT
  USE MOD_DIAG
  USE MOD_VIEW
  USE MOD_SLAG
  USE MOD_MPDTA
  USE MOD_STEPO

  IMPLICIT NONE

  INTEGER                              :: TEMPS,JOUR

  TYPE(DOMAIN)                         :: DOM
  TYPE(PROGVAR)                        :: ST0,ST9,ST1
  TYPE(SMILAG)                         :: YSLAG
  TYPE(MPDATA)                         :: YMDTA
  TYPE(RHSVAR)                         :: RHSG_A,RHSG_D

  INTEGER(8)                           :: NSTEP,NITER
  REAL(8)                              :: PDT


  !* Allocation
  CALL SUALLOC_PART1(ST0,ST1,ST9)
  CALL SUALLOC_PART2(YSLAG,YMDTA,RHSG_A,RHSG_D)
  
  !* Initialisation
  CALL MODEL_INITIALIZATION(ST0,DOM)
  CALL OUTPUT_INITIAL_FIELDS(ST0,DOM)
  
  WRITE(*,*) "Integration... \n"
  WRITE(*,*)
  CALL SCREEN_START(TEMPS,JOUR,DOM,ST0)

  CALL SWAPP(ST9,ST0)
  
  ! Dynamic Loop
  DO NSTEP = 1, NTIMESTEP
     
    CALL RHS_EXP(RHSG_A,ST0) 
    CALL SET_VELOCITY(ST1,DOM,NSTEP) 

    DO NITER = 0,NSITER
      IF (LSLAG) THEN
         CALL SMILAG_TRANSPORT_SCHEME(RHSG_D,RHSG_A,YSLAG,  &
              & ST9%U,ST9%V,ST0%U,ST0%V,ST0%DUDX,ST0%DVDY,  &
              & ST1%U,ST1%V,ST1%DUDX,ST1%DVDY,              &
              & NSTEP,NITER,C_ADVECT)
      ELSE 
         CALL MPDATA_TRANSPORT_SCHEME(RHSG_D,RHSG_A,YMDTA,  &
              & ST0%M,ST0%U,ST0%V,ST9%U,ST9%V,ST1%U,ST1%V,  &
              & NSTEP,NITER,C_ADVECT)
      END IF   
    END DO
      
    CALL SWAPP(ST9,ST0)
    CALL RHS_IMP(ST0,ST1,RHSG_D)

    IF (MOD(NSTEP,NFREQOUTPUT) == 0) THEN
       CALL OUTPUT_CURRENT_FIELDS(ST0,DOM,NSTEP)
    END IF
    CALL SCREEN(TEMPS,ST0,DOM,NSTEP)
    
  END DO
 
  CALL SCREEN_END(TEMPS,JOUR)
  CALL FREEMEM(ST0,ST9,ST1,YSLAG,YMDTA,RHSG_A,RHSG_D)

END PROGRAM MAIN

!=====================================================!
!=====================================================!
!=====================================================!
