!################   CALCUL DES RHS   #################!
!#                                                   #!
!# auteur : C.Colavolpe & Voitus F.                  #!
!# sujet  : Discrétisation temporelle                #!
!#                                                   #!
!#####################################################!

!=====================================================!
!====================  MODULE DIAG ===================!
!=====================================================!

MODULE MOD_DIAG

  USE MOD_PARAM, ONLY : GXPT,GYPT,RG,RPSUR,ZERO,ONE,TWO,HALF
  
CONTAINS

 !*****************************************************!
 !***************  LISTE DES ROUTINES  ****************!
 !*****************************************************!
  
 !*******************************************************************************!
 !*******************************************************************************!
 SUBROUTINE SP_NORMS(P,PNORM,CNORM_TYPE)
   
   IMPLICIT NONE

   REAL(8), DIMENSION(GXPT,GYPT),       INTENT(IN)    :: P   
   REAL(8),                             INTENT(OUT)   :: PNORM
   CHARACTER(LEN=2),                    INTENT(IN)    :: CNORM_TYPE
   INTEGER(8)                                         :: NX,NY
   REAL(8)                                            :: ZNORM, ZGPTOT

   ZGPTOT = REAL(GYPT,8)*REAL(GXPT,8)

   IF (CNORM_TYPE == 'L2') THEN
     ZNORM = 0.d0
       DO NX=1,GXPT
         DO NY =1,GYPT
            ZNORM = ZNORM + P(NX,NY)**2
         END DO
       END DO
     PNORM = DSQRT(ZNORM/ZGPTOT)
   ELSE IF (CNORM_TYPE == 'L1') THEN
     ZNORM = 0.d0
       DO NX=1,GXPT
         DO NY =1,GYPT
            ZNORM = ZNORM + ABS(P(NX,NY))
         END DO
       END DO  
     PNORM = (ZNORM/ZGPTOT)
   ELSE
     PNORM = MAXVAL(ABS(P))   
   END IF

  END SUBROUTINE SP_NORMS
  !*******************************************************************************!
  !*******************************************************************************!
  !*******************************************************************************!

  
END MODULE MOD_DIAG

!=====================================================!
!=====================================================!
!=====================================================!
