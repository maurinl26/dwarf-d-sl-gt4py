!#####################################################!
!#                                                   #!
!# auteur : F.Voitus                                 #!
!# sujet  : Organisation de la boucle temporelle     #!
!#                                                   #!
!#####################################################!

!======================================================!
!===================== MODULE STEPO ===================!
!======================================================!

MODULE MOD_STEPO

  USE OMP_LIB
  USE MPI
  USE MOD_PARAMETER, ONLY :  NLEV,LXPT,LYPT,HALO,GXPT,GYPT,NPARTX,NPARTY,TRANCHEX,TRANCHEY, &
                          &  RDZ,RR,RCV,RCP,RG,RPI,RKAPPA,ZERO,ONE,TWO,RLX,RLY,RDX,RDY,RXC, &
                          &  RYC,NB_PROCS,NFREQOUTPUT,NTIMESTEP,NTAG,C_OPTION
  USE MOD_STRUCTURE, ONLY :  GEOMETRY,BACKGROUND,VARPROG,RHS_LOC,VARGLOB,DOMAIN,MPDATA,SMILAG

CONTAINS

  !********************************************************************************************
  !********************************************************************************************
  !********************************************************************************************
  !********************************************************************************************
  SUBROUTINE MODEL_INITIALIZATION(ST,STB,STG,YSLAG,YMDTA,DOM,GEO,ABS_GLB,ORD_GLB,KDIM,MYPROC,CODE)
    
    USE MOD_STRUCTURE,  ONLY : DEF_DOMAIN_GEOMETRY,DEF_BACKGROUND_STATE,  &
                             & DEF_A_B_FUNCTIONS,DEF_INITIAL_PERTURBATION,SUALLOC_GBL
    USE MOD_VIEW,       ONLY : OUTPUT_INITIAL_FIELDS
    USE MOD_COMM,       ONLY : TRANSVAR_LOC,TRANSGEO_LOC

    IMPLICIT NONE

    TYPE(VARPROG),                    INTENT(INOUT) :: ST
    TYPE(BACKGROUND),                 INTENT(INOUT) :: STB
    TYPE(VARGLOB),                    INTENT(INOUT) :: STG
    TYPE(SMILAG),                     INTENT(INOUT) :: YSLAG
    TYPE(MPDATA),                     INTENT(INOUT) :: YMDTA
    TYPE(DOMAIN),                     INTENT(INOUT) :: DOM
    TYPE(GEOMETRY),                   INTENT(INOUT) :: GEO
    INTEGER(8), DIMENSION(0:NB_PROCS-1), INTENT(IN) :: ABS_GLB,ORD_GLB
    INTEGER(8),                       INTENT(IN)    :: KDIM
    INTEGER,                          INTENT(IN)    :: MYPROC
    INTEGER,                          INTENT(INOUT) :: CODE

    INTEGER(8)                                      :: NY,NX,JL,NL
    INTEGER                                         :: NPROC
    INTEGER, DIMENSION(MPI_STATUS_SIZE)             :: STATUS
    

    CALL DEF_BACKGROUND_STATE(STB)
    CALL DEF_A_B_FUNCTIONS(GEO%AH,GEO%BH,GEO%ETAF,        &
         & GEO%ETAH,GEO%DETAF,GEO%DETAH,GEO%DELTB,C_OPTION)

    IF (MYPROC == 0) THEN
       CALL SUALLOC_GBL(RHSG,YSLAG,YMDTA,DOM)
       WRITE(*,*) " "
       WRITE(*,*) "Initialisation... \n"
       CALL DEF_DOMAIN_GEOMETRY(DOM,GEO%AH,GEO%BH)
       CALL DEF_INITIAL_PERTURBATION(STG,STB,DOM)
       CALL OUTPUT_INITIAL_FIELDS(STG,STB,DOM)
    END IF
 
    !*************************************!
    !waiting for the end of myproc=0 task !
    CALL MPI_BARRIER(MPI_COMM_WORLD,CODE) !
    !*************************************!
    
    CALL TRANSGEO_LOC(GEO,DOM,ABS_GLB,ORD_GLB,MYPROC,CODE)
    CALL TRANSVAR_LOC(ST,STG,ABS_GLB,ORD_GLB,MYPROC,CODE)
    
    CALL UPDATE_METRIC_TERMS(ST,GEO)

  END SUBROUTINE MODEL_INITIALIZATION
  !*****************************************************************************************
  !*****************************************************************************************
  !*****************************************************************************************
  SUBROUTINE UPDATE_METRIC_TERMS(ST,GEO)
   
    USE MOD_DIAGNOSTIC, ONLY : DIAG_MASS_METRIC,DIAG_ETADOT, &
                             & DIAG_UF2UH,DIAG_VERTICAL_METRICS
    
    IMPLICIT NONE

    TYPE(VARPROG),    INTENT(INOUT) :: ST
    TYPE(GEOMETRY),   INTENT(INOUT) :: GEO

    INTEGER(8)                      :: I,NX,NY 
        
    DO NY = 1, LYPT
       DO NX = 1, LXPT
         I = (HALO+NY-1)*TRANCHEX + HALO + NX

         CALL DIAG_VERTICAL_METRICS(GEO%PIH(:,I),GEO%PIF(:,I),            &
              & GEO%DELTAPI(:,I),GEO%DELTA(:,I),GEO%ALPHA(:,I),           &
              & GEO%EPS(:,I),GEO%AH,GEO%BH,ST%PIS(I))
         CALL DIAG_MASS_METRIC(ST%MF(:,I),ST%MH(:,I),GEO%DELTAPI(:,I),    &
              & GEO%PIF(:,I),GEO%DETAF,GEO%PIH(:,I),GEO%DETAH)
         CALL DIAG_ETADOT(ST%ETADOTF(:,I),ST%ETADOTH(:,I),ST%MF(:,I),     &
              & ST%MH(:,I),ST%U(:,I),ST%V(:,I),ST%DIV(:,I),ST%DPISDX(I),  &
              & ST%DPISDY(I),GEO%DELTAPI(:,I),GEO%DELTB,GEO%BH)
         CALL DIAG_UF2UH(ST%UH(:,I),ST%VH(:,I),ST%U(:,I),ST%V(:,I),       &
              & GEO%EPS(:,I))
         
      END DO  
    END DO
    
  END SUBROUTINE UPDATE_METRIC_TERMS
  !*************************************************************************************!
  !*************************************************************************************!
  !*************************************************************************************!
  !*************************************************************************************!
  SUBROUTINE RHS_MDL(RHS,ST,GEO)

    USE MOD_DIAGNOSTIC,  ONLY : DIAG_PISDOT

    IMPLICIT NONE

    TYPE(RHS_LOC),       INTENT(INOUT) :: RHS
    TYPE(VARPROG),       INTENT(IN)    :: ST
    TYPE(GEOMETRY),      INTENT(IN)    :: GEO

    INTEGER(8)                         :: NX,NY,NL,J
    REAL(8)                            :: ZNDIV
    
    DO NY = 1, LYPT
       DO NX = 1, LXPT
          NL = (HALO+NY-1)*TRANCHEX + HALO + NX

          CALL DIAG_PISDOT(ZNDIV,ST%U(:,NL),ST%V(:,NL),ST%DIV(:,NL), &
               & GEO%PIS(NL),GEO%PIS_REF(NL),GEO%DPISDX_REF(NL),     &
               & GEO%DPISDY_REF(NL),GEO%DELTB,GEO%AH)
        
          RHS%PIS(NL)     = - ZNDIV
          RHS%T(:,NL)     = ZERO
          RHS%Q(:,NL)     = ZERO
             
       END DO
    END DO

  END SUBROUTINE RHS_MDL
  !********************************************************************!
  !********************************************************************!
  !********************************************************************!
  !********************************************************************!
  SUBROUTINE RHS_EXP(RHS,ST0,RHSM0,GEO,PDT)

    USE MOD_PARAMETER, ONLY : RPIS2D
    
    IMPLICIT NONE

    TYPE(RHS_LOC),    INTENT(INOUT) :: RHS
    TYPE(VARPROG),    INTENT(IN)    :: ST0
    TYPE(RHS_LOC),    INTENT(IN)    :: RHSM0
    TYPE(GEO_NHNL),   INTENT(IN)    :: GEO
    REAL(8),          INTENT(IN)    :: PDT    

    INTEGER(8)                      :: NX,NY,NL,JL
    
    DO NY = 1,LYPT
       DO NX = 1,LXPT
          NL = (HALO+NY-1)*TRANCHEX + HALO + NX 

          RHS%Q(:,NL)  = ST0%Q(:,NL)    
          RHS%T(:,NL)  = ST0%T(:,NL) 
          RHS%PIS(NL)  = RPIS2D*(ST%PIS(NL)-GEO%PIS_REF(NL)) &
                       & + RPIS2D*PDT*RHS0M%PIS(NL)
       END DO
    END DO
    
  END SUBROUTINE RHS_EXP
  !********************************************************************!
  !********************************************************************!
  !********************************************************************!
  !********************************************************************!
  SUBROUTINE RHS_IMP(ST,RHSM1,RHSDEP,GEO,STB,PDT)

    USE MOD_PARAMETER, ONLY : RS,RPIS2D
    
    IMPLICIT NONE

    TYPE(VARPROG),    INTENT(OUT)   :: ST
    TYPE(GEO_NHNL),   INTENT(INOUT) :: GEO
    TYPE(RHS_LOC),    INTENT(IN)    :: RHSM1,RHSDEP
    TYPE(BACKGROUND), INTENT(IN)    :: STB
    REAL(8),          INTENT(IN)    :: PDT    

    INTEGER(8)                      :: NX,NY,NL,JL
    
    DO NY = 1, LYPT
       DO NX = 1, LXPT
          NL =  (HALO+NY-1)*TRANCHEX + HALO + NX
          
          ST%T(:,NL)     = RHSDEP%T(:,NL)  
          ST%Q(:,NL)     = RHSDEP%Q(:,NL)
          ST%PIS(NL)     = RHSDEP%PIS(NL)            &
                         & + RPIS2D*(GEO%PIS_REF(NL)+PDT*RHSM1%PIS(NL))
          DO JL=1,NLEV
             ST%V(JL,NL) = STB%V(JL)
             ST%U(JL,NL) = STB%U(JL)*(GEO%ETAH(JL)-GEO%ETAH(JL-1)) &
                         & /( ((GEO%AH(JL)-GEO%AH(JL-1))/RP00)     &
                         & +(RS*(ST%PIS(NL)/RP00)+(ONE-RS))*(GEO%BH(JL)-GEO%BH(JL-1)) )
          END DO   
       END DO
    END DO
    
  END SUBROUTINE RHS_IMP
  !************************************************************************************!
  !************************************************************************************!
  !************************************************************************************!
  !************************************************************************************!
  SUBROUTINE SWAPP(ST9,ST0)

    IMPLICIT NONE

    TYPE(VARPROG),  INTENT(INOUT)  :: ST9
    TYPE(VARPROG),  INTENT(IN)     :: ST0

    INTEGER(8)                     :: NX,NY,NL

    DO NY = 1, LYPT
       DO NX = 1, LXPT
          NL =  (HALO+NY-1)*TRANCHEX + HALO + NX

          ST9%U(:,NL)       = ST0%U(:,NL)    
          ST9%V(:,NL)       = ST0%V(:,NL)        
          ST9%T(:,NL)       = ST0%T(:,NL)        
          ST9%Q(:,NL)       = ST0%Q(:,NL)      
          ST9%DIV(:,NL)     = ST0%DIV(:,NL)          
          ST9%PIS(NL)       = ST0%PIS(NL)

          ST9%MF(:,NL)      = ST0%MF(:,NL)     
          ST9%MH(:,NL)      = ST0%MH(:,NL)
          
          ST9%UH(:,Nl)      = ST0%UH(:,NL)     
          ST9%VH(:,NL)      = ST0%VH(:,NL)
          ST9%ETADOTF(:,NL) = ST0%ETADOTF(:,NL)        
          ST9%ETADOTH(:,NL) = ST0%ETADOTH(:,NL)
          
       END DO
    END DO

  END SUBROUTINE SWAPP
  !*********************************************************************************!
  !*********************************************************************************!
  !*********************************************************************************!
  !*********************************************************************************!
  SUBROUTINE SPC_DERV(STL,STG,STB,DOM,ABS_GLB,ORD_GLB,MYPROC,CODE,KSTEP)

    USE MOD_FFT,       ONLY : FFTDIR_UV,FFTDIV
    USE MOD_COMM,      ONLY : TRANSVAR_GLB,TRANSVAR_LOC
    USE MOD_VIEW,      ONLY : OUTPUT_CURRENT_FIELDS
    
    IMPLICIT NONE

    TYPE(VARGLOB),                       INTENT(INOUT) :: STG
    TYPE(VARPROG),                       INTENT(INOUT) :: STL
    TYPE(BACKGROUND),                    INTENT(IN)    :: STB
    TYPE(DOMAIN),                        INTENT(IN)    :: DOM
    INTEGER,                             INTENT(INOUT) :: CODE
    INTEGER,                             INTENT(IN)    :: MYPROC
    INTEGER(8), DIMENSION(0:NB_PROCS-1), INTENT(IN)    :: ABS_GLB,ORD_GLB
    INTEGER(8),                          INTENT(IN)    :: KSTEP
    
    INTEGER(8)                                         :: JL,NY,NL
    INTEGER                                            :: NPROC
    INTEGER, DIMENSION(MPI_STATUS_SIZE)                :: STATUS
    LOGICAL                                            :: LPRINTLEV

    
    !* TRANSFERT TO GLOBAL ARRAYS ON (MYPROC = 0)
    CALL TRANSVAR_GLB(STG,STL,ABS_GLB,ORD_GLB,MYPROC,CODE)

    !* FFT TRANSFORM AND DERIVATIVES
    IF (MYPROC .EQ. 0) THEN

       DO JL = 1,NLEV
          CALL FFTDIR_UV(STG%U(JL,:,:),STG%V(JL,:,:))
          CALL FFTDIV(STG%U(JL,:,:),STG%V(JL,:,:),STG%DIV(JL,:,:),LALIAZ=LRALIAZ)
       END DO
       
       LPRINTLEV = (MOD(KSTEP,NFREQOUTPUT).EQ.0).OR.((KSTEP==NTIMESTEP) &
                 & .AND.(MOD(NTIMESTEP,NFREQOUTPUT).NE.0))

       IF (LPRINTLEV) THEN
          CALL OUTPUT_CURRENT_FIELDS(STG,STB,DOM,KSTEP)
       END IF
       
    END IF

    !* TRANSFERT TO LOCAL MPI TASKS
    CALL TRANSVAR_LOC(STL,STG,ABS_GLB,ORD_GLB,MYPROC,CODE)
    
  END SUBROUTINE SPC_DERV
  !***************************************************************************************!
  !***************************************************************************************!
  !***************************************************************************************!
  !***************************************************************************************!

END MODULE MOD_STEPO
!=====================================================!
!=====================================================!
!=====================================================!
