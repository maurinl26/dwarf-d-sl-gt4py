!#####################  DERIVATIVE ####################!
!#                                                    #!
!# auteur : F.Voitus                                  #!
!# sujet  : Finite difference derivatives             #!
!#                                                    #!
!######################################################!

!======================================================!
!=================== MODULE DERIVATIVE ================!
!======================================================!

MODULE MOD_GRAD

  USE MOD_SHARE, ONLY : NXPT,NLEV,RLX,TWO,ONE,ZERO,HALF,RPI, &
                      & RFD_WEIGHT,NPERIO,HALO
  USE MOD_FFT,   ONLY : SPEC_GRAD,SPEC_DIV
  
CONTAINS

  !###########################################################################!
  !***************************************************************************!
  !***************************** LISTE DES ROUTINES **************************!
  !***************************************************************************!
  !###########################################################################!
  !***************************************************************************!
  !***************************************************************************!
  !***************************************************************************!
  !***************************************************************************!
  !***************************************************************************!
  !***************************************************************************!
  !***************************************************************************!
  !###########################################################################!
  SUBROUTINE GRAD_1D(PDVDX,PVAR)

    IMPLICIT NONE

    REAL(8), DIMENSION(NXPT), INTENT(IN)    :: PVAR
    REAL(8), DIMENSION(NXPT), INTENT(INOUT) :: PDVDX

    INTEGER(8)                              :: I,N

    DO I=1,NXPT
       PDVDX(I) = ZERO
       DO N=1,HALO
             PDVDX(I) = PDVDX(I)  &
                      & + RFD_WEIGHT(N)*PVAR(NPERIO(I+N)) &
                      & - RFD_WEIGHT(N)*PVAR(NPERIO(I-N))
       END DO   
    END DO
    
  END SUBROUTINE GRAD_1D
  !***********************************************************!
  !***********************************************************!
  !***********************************************************!
  !***********************************************************!
  SUBROUTINE GRAD(PDVDX,PVAR,KTOP)

    IMPLICIT NONE

    REAL(8), DIMENSION(KTOP:NLEV,NXPT), INTENT(IN)    :: PVAR
    REAL(8), DIMENSION(KTOP:NLEV,NXPT), INTENT(INOUT) :: PDVDX
    INTEGER(8),                         INTENT(IN)    :: KTOP

    INTEGER(8)                                        :: I,N,J

    DO J=KTOP,NLEV
       DO I=1,NXPT
          PDVDX(J,I) = ZERO
          DO N=1,HALO
             PDVDX(J,I) = PDVDX(J,I)                          &
                        & + RFD_WEIGHT(N)*PVAR(J,NPERIO(I+N)) &
                        & - RFD_WEIGHT(N)*PVAR(J,NPERIO(I-N))
          END DO
       END DO   
    END DO

  END SUBROUTINE GRAD
  !***********************************************************!
  !***********************************************************!
  !***********************************************************!
  !***********************************************************!
  !***********************************************************!
  !***********************************************************!
  !***********************************************************!
  !***********************************************************!
  SUBROUTINE DIV(PDIV,PU)

    IMPLICIT NONE

    REAL(8), DIMENSION(NLEV,NXPT),   INTENT(INOUT) :: PDIV
    REAL(8), DIMENSION(NLEV,NXPT),   INTENT(IN)    :: PU
    
    INTEGER(8)                                     :: I,N,J

    DO J=1,NLEV
       DO I=1,NXPT
          PDIV(J,I)    = ZERO
          DO N=1,HALO
             PDIV(J,I) = PDIV(J,I)                           &
                       & + RFD_WEIGHT(N)*( PU(J,NPERIO(I+N)) &
                       & - PU(J,NPERIO(I-N)) )
          END DO
       END DO
    END DO
    
  END SUBROUTINE DIV
  !***********************************************************!
  !***********************************************************!
  !***********************************************************!
  !***********************************************************!
  SUBROUTINE DIVU(PDIV,PU,LSPEC,LALIAZ)

    IMPLICIT NONE

    REAL(8), DIMENSION(NLEV,NXPT),   INTENT(INOUT) :: PDIV
    REAL(8), DIMENSION(NLEV,NXPT),   INTENT(INOUT) :: PU
    LOGICAL, OPTIONAL,               INTENT(IN)    :: LSPEC
    LOGICAL, OPTIONAL,               INTENT(IN)    :: LALIAZ
    
    LOGICAL                                        :: LL_ALIAZ
    LOGICAL                                        :: LL_SPEC

     IF (PRESENT(LALIAZ)) THEN
       LL_ALIAZ = LALIAZ
    ELSE
       LL_ALIAZ  = .FALSE.
    END IF
    
    IF (PRESENT(LSPEC)) THEN
       LL_SPEC = LSPEC
    ELSE
       LL_SPEC = .FALSE.
    END IF

    IF (LL_SPEC) THEN
       CALL SPEC_DIV(PDIV,PU,LL_ALIAZ)   
    ELSE
       CALL DIV(PDIV,PU)
    END IF
       
  END SUBROUTINE DIVU
  !****************************************************************!
  !****************************************************************!
  !****************************************************************!
  !****************************************************************!
  SUBROUTINE DIV_MV(PDIV,PU,PM)

    IMPLICIT NONE

    REAL(8), DIMENSION(NLEV,NXPT),   INTENT(INOUT) :: PDIV
    REAL(8), DIMENSION(NLEV,NXPT),   INTENT(IN)    :: PM,PU
    
    INTEGER(8)                                     :: I,N,J

    DO J=1,NLEV
       DO I=1,NXPT
          PDIV(J,I)    = ZERO
          DO N=1,HALO
             PDIV(J,I) = PDIV(J,I) + RFD_WEIGHT(N)*((            &
                       & PM(J,NPERIO(I+N))*PU(J,NPERIO(I+N))     &
                       & - PM(J,NPERIO(I-N))*PU(J,NPERIO(I-N)) ) &
                       & /PM(J,I))
          END DO
       END DO
    END DO
    
  END SUBROUTINE DIV_MV
  !***************************************************************!
  !***************************************************************!
  !***************************************************************!
  !***************************************************************!
  SUBROUTINE GRADP(PDVDX,PVAR,KTOP,LSPEC,LALIAZ)

    IMPLICIT NONE

    REAL(8), DIMENSION(KTOP:NLEV,NXPT), INTENT(INOUT) :: PVAR
    REAL(8), DIMENSION(KTOP:NLEV,NXPT), INTENT(INOUT) :: PDVDX
    INTEGER(8),                         INTENT(IN)    :: KTOP
    LOGICAL, OPTIONAL,                  INTENT(IN)    :: LSPEC
    LOGICAL, OPTIONAL,                  INTENT(IN)    :: LALIAZ                                
    LOGICAL                                           :: LL_ALIAZ,LL_SPEC

    IF (PRESENT(LALIAZ)) THEN
       LL_ALIAZ = LALIAZ
    ELSE
       LL_ALIAZ = .FALSE.
    END IF
    
    IF (PRESENT(LSPEC)) THEN
       LL_SPEC = LSPEC
    ELSE
       LL_SPEC = .FALSE.
    END IF
    
    IF (LL_SPEC) THEN
       CALL SPEC_GRAD(PVAR,PDVDX,KTOP,LL_ALIAZ)
    ELSE   
       CALL GRAD(PDVDX,PVAR,KTOP)
    END IF
    
  END SUBROUTINE GRADP
  !*******************************************************************************!
  !*******************************************************************************!
  !*******************************************************************************!
  !*******************************************************************************!
  !*******************************************************************************!
  !*******************************************************************************!
  
END MODULE MOD_GRAD

!=====================================================!
!=====================================================!
!=====================================================!
