!######################   FFT  #######################!
!#                                                   #!
!# auteur : F.Voitus                                 #!
!# sujet  : Les transformées de Fourier              #!
!#                                                   #!
!#####################################################!

!=====================================================!
!===================  MODULE FFT  ====================!
!=====================================================!

MODULE MOD_FFT

  USE MOD_SHARE, ONLY : NXPT,NLEV,RLX,RPI,TWO,ONE,ZERO,NSMAX
  USE GFT
  
CONTAINS

  !*****************************************************!
  !***************  LISTE DES ROUTINES  ****************!
  !*****************************************************!

   !**********************  Transformee de Fourier d'un RHS  **********************!
  SUBROUTINE FFTDIR(PX,KTOP)

    IMPLICIT NONE

    REAL(8), DIMENSION(KTOP:NLEV,NXPT), INTENT(INOUT) :: PX
    INTEGER(8),                         INTENT(IN)    :: KTOP
    INTEGER, PARAMETER                                :: N = NXPT, SI = +1
    INTEGER                                           :: I,J
    REAL(kind=GFT_Prec),PARAMETER                     :: SC = 1.0_GFT_Prec
    COMPLEX(kind=GFT_Prec), DIMENSION(0:NXPT,1)       :: XP,YP
    TYPE(GFT_CC)                                      :: CC
    
    !... Set size of FFT and cc object
    CALL GFT_set_fft(NX=N, NY=1, FFT=CC)
    !... Initialize Array to be Transformed
    DO J = KTOP, NLEV
       DO I = 0, NXPT-1
          XP(I,1) = PX(J,I+1)
       END DO
       !... Forward FFT
       CALL GFT_do_fft(FFT=CC, isign=SI, scale=SC, c_in=XP, c_out=YP)
       PX(J,1)    = DBLE(REAL(YP(0,1)))
       PX(J,2)    = DBLE(REAL(YP(NXPT/2,1)))
       DO I = 1, NXPT/2-1
          PX(J,2*I+1) = DBLE(REAL(YP(I,1)))
          PX(J,2*I+2) = DBLE(AIMAG(YP(I,1)))
       END DO
    END DO
    CALL GFT_end_fft(FFT=CC)

  END SUBROUTINE FFTDIR
 !**************************************************************
  SUBROUTINE FFTINV(PX,KTOP,LALIAZ)

    IMPLICIT NONE

    REAL(8), DIMENSION(KTOP:NLEV,NXPT), INTENT(INOUT) :: PX
    INTEGER(8),                         INTENT(IN)    :: KTOP
    LOGICAL,                            INTENT(IN)    :: LALIAZ
    INTEGER, PARAMETER                                :: N = NXPT, SI = -1
    INTEGER                                           :: I,J 
    COMPLEX(8), PARAMETER                             :: ZI = (0.d0,1.d0)
    REAL(kind=GFT_Prec), PARAMETER                    :: SC = 1.0_GFT_Prec/REAL(NXPT, kind=GFT_Prec)
    COMPLEX(kind=GFT_Prec), DIMENSION(0:NXPT,1)       :: XP,YP
    TYPE(GFT_CC)                                      :: CC
    
    IF (LALIAZ) THEN
      DO I = NSMAX, NXPT/2-1
         PX(:,2*I+1)  = ZERO
         PX(:,2*I+2)  = ZERO
      END DO
    ENDIF

    !... Set size of FFT and cc object
    CALL GFT_set_fft(NX=N, NY=1, FFT=CC)

    !... Initialize Array to be Transformed
    DO J = KTOP, NLEV

       !... Backward FFT
       YP(0,1)         = PX(J,1)
       YP(NXPT/2,1)    = PX(J,2)
       DO I = 1, NXPT/2-1
          YP(I,1)         = PX(J,2*I+1) + ZI*PX(J,2*I+2)
          YP(NXPT-I,1)    = PX(J,2*I+1) - ZI*PX(J,2*I+2)
       END DO

       CALL GFT_do_fft(FFT=CC, isign=SI, scale=SC, c_in=YP, c_out=XP)

       !... Free FFT object
       DO I = 0, NXPT-1
          PX(J,I+1) = REAL(XP(I,1),8)
       END DO

    END DO

    CALL GFT_end_fft(FFT=CC)

  END SUBROUTINE FFTINV
  !*********************************************************
  SUBROUTINE FFTDER(PVAR,PDVDX,KTOP,LALIAZ)

    IMPLICIT NONE

    REAL(8), DIMENSION(KTOP:NLEV,NXPT), INTENT(INOUT) :: PVAR
    REAL(8), DIMENSION(KTOP:NLEV,NXPT), INTENT(OUT)   :: PDVDX
    INTEGER(8),                         INTENT(IN)    :: KTOP
    LOGICAL,                            INTENT(IN)    :: LALIAZ                            
    INTEGER, PARAMETER                                :: N = NXPT, SI = -1
    INTEGER                                           :: I,J
    REAL(8)                                           :: ZK
    COMPLEX(8), PARAMETER                             :: ZI = (0.d0,1.d0)
    REAL(kind=GFT_Prec), PARAMETER                    :: SC = 1.0_GFT_Prec/REAL(NXPT, kind=GFT_Prec)
    COMPLEX(kind=GFT_Prec), DIMENSION(0:NXPT,1)       :: XP,YP,XPDV,YPDV
    TYPE(GFT_CC)                                      :: CC
    
    
    IF (LALIAZ) THEN
      DO I = NSMAX, NXPT/2-1
         PVAR(:,2*I+1) = ZERO
         PVAR(:,2*I+2) = ZERO
      END DO
    ENDIF

    DO I = 0, NXPT/2-1
      ZK = REAL(I,8)*RPI*TWO/RLX
      DO J = KTOP, NLEV
        PDVDX(J,2*I+1) =  ZK*PVAR(J,2*I+2)
        PDVDX(J,2*I+2) = -ZK*PVAR(J,2*I+1)
      END DO
    END DO

    !... Set size of FFT and cc object
    CALL GFT_set_fft(NX=N, NY=1, FFT=CC)

    !... Initialize Array to be Transformed
    DO J = KTOP, NLEV
       !... Backward FFT
       YP(0,1)         = PVAR(J,1)
       YP(NXPT/2,1)    = PVAR(J,2)
       YPDV(0,1)       = PDVDX(J,1)
       YPDV(NXPT/2,1)  = PDVDX(J,2)
       DO I = 1, NXPT/2-1
          YP(I,1)         = PVAR(J,2*I+1)  + ZI*PVAR(J,2*I+2)
          YP(NXPT-I,1)    = PVAR(J,2*I+1)  - ZI*PVAR(J,2*I+2)
          YPDV(I,1)       = PDVDX(J,2*I+1) + ZI*PDVDX(J,2*I+2)
          YPDV(NXPT-I,1)  = PDVDX(J,2*I+1) - ZI*PDVDX(J,2*I+2)
       END DO

       CALL GFT_do_fft(FFT=CC, isign=SI, scale=SC, c_in=YP, c_out=XP)
       CALL GFT_do_fft(FFT=CC, isign=SI, scale=SC, c_in=YPDV, c_out=XPDV)

       !... Free FFT object
       DO I = 0, NXPT-1
          PVAR(J,I+1)  = REAL(XP(I,1),8)
          PDVDX(J,I+1) = REAL(XPDV(I,1),8)
       END DO
    END DO
    
    CALL GFT_end_fft(FFT=CC)

  END SUBROUTINE FFTDER
  !***********************************************************************************!
  !***********************************************************************************!
  !***********************************************************************************!
  !***********************************************************************************!
  !***********************************************************************************!
  !***********************************************************************************!
  !***********************************************************************************!
  !***********************************************************************************!
  !***********************************************************************************!
  !***********************************************************************************!
  !***********************************************************************************!
  !***********************************************************************************!
  !***********************************************************************************!
  !***********************************************************************************!
  !***********************************************************************************!
  SUBROUTINE SPEC_GRAD(PVAR,PDVDX,KTOP,LALIAZ)
    
    IMPLICIT NONE

    REAL(8), DIMENSION(KTOP:NLEV,NXPT), INTENT(INOUT) :: PVAR
    REAL(8), DIMENSION(KTOP:NLEV,NXPT), INTENT(OUT)   :: PDVDX
    INTEGER(8),                         INTENT(IN)    :: KTOP
    LOGICAL, OPTIONAL,                  INTENT(IN)    :: LALIAZ
    
    INTEGER, PARAMETER                                :: N = NXPT, SI = -1
    INTEGER                                           :: I,J
    REAL(8)                                           :: ZK
    COMPLEX(8), PARAMETER                             :: ZI = (0.d0,1.d0)
    REAL(kind=GFT_Prec), PARAMETER                    :: SC = 1.0_GFT_Prec/REAL(NXPT, kind=GFT_Prec)
    COMPLEX(kind=GFT_Prec), DIMENSION(0:NXPT,1)       :: XP,YP,XPDV,YPDV
    TYPE(GFT_CC)                                      :: CC
    LOGICAL                                           :: LLFDIR,LLALIAZ

    IF (PRESENT(LALIAZ)) THEN
       LLALIAZ = LALIAZ
    ELSE
       LLALIAZ = .FALSE.
    END IF   

    CALL FFTDIR(PVAR,KTOP)

    IF (LLALIAZ) THEN
      DO I = NSMAX, NXPT/2-1
         PVAR(:,2*I+1) = ZERO
         PVAR(:,2*I+2) = ZERO
      END DO
    ENDIF
    
    
    !... Set size of FFT and cc object
    CALL GFT_set_fft(NX=N, NY=1, FFT=CC)

    !... Initialize Array to be Transformed
    DO J = KTOP, NLEV

       DO I = 0, NXPT/2-1
          ZK = REAL(I,8)*RPI*TWO/RLX
          PDVDX(J,2*I+1) =  ZK*PVAR(J,2*I+2)
          PDVDX(J,2*I+2) = -ZK*PVAR(J,2*I+1)
       END DO
       
       !... Backward FFT
       YP(0,1)            = PVAR(J,1)
       YP(NXPT/2,1)       = PVAR(J,2)
       YPDV(0,1)          = PDVDX(J,1)
       YPDV(NXPT/2,1)     = PDVDX(J,2)
       DO I = 1, NXPT/2-1
          YP(I,1)         = PVAR(J,2*I+1)  + ZI*PVAR(J,2*I+2)
          YP(NXPT-I,1)    = PVAR(J,2*I+1)  - ZI*PVAR(J,2*I+2)
          YPDV(I,1)       = PDVDX(J,2*I+1) + ZI*PDVDX(J,2*I+2)
          YPDV(NXPT-I,1)  = PDVDX(J,2*I+1) - ZI*PDVDX(J,2*I+2)
       END DO

       CALL GFT_do_fft(FFT=CC, isign=SI, scale=SC, c_in=YP, c_out=XP)
       CALL GFT_do_fft(FFT=CC, isign=SI, scale=SC, c_in=YPDV, c_out=XPDV)

       !... Free FFT object
       DO I = 0, NXPT-1
          PVAR(J,I+1)  = REAL(XP(I,1),8)
          PDVDX(J,I+1) = REAL(XPDV(I,1),8)
       END DO
       
    END DO
    
    CALL GFT_end_fft(FFT=CC)

  END SUBROUTINE SPEC_GRAD
  !*******************************************************************************!
  !*******************************************************************************!
  !*******************************************************************************!
  !*******************************************************************************!
  SUBROUTINE SPEC_DIV(PDIV,PU,LALIAZ)
    
    IMPLICIT NONE

    REAL(8), DIMENSION(NLEV,NXPT), INTENT(INOUT) :: PDIV
    REAL(8), DIMENSION(NLEV,NXPT), INTENT(INOUT) :: PU
    LOGICAL, OPTIONAL,             INTENT(IN)    :: LALIAZ
    
    INTEGER, PARAMETER                           :: N = NXPT, SI = -1
    INTEGER                                      :: I,J
    REAL(8)                                      :: ZK
    COMPLEX(8), PARAMETER                        :: ZI = (0.d0,1.d0)
    REAL(kind=GFT_Prec), PARAMETER               :: SC = 1.0_GFT_Prec/REAL(NXPT, kind=GFT_Prec)
    COMPLEX(kind=GFT_Prec), DIMENSION(0:NXPT,1)  :: XU,YU,XDIV,YDIV
    TYPE(GFT_CC)                                 :: CC
    LOGICAL                                      :: LLALIAZ

    
    IF (PRESENT(LALIAZ)) THEN
       LLALIAZ = LALIAZ
    ELSE
       LLALIAZ = .FALSE.
    END IF

    CALL FFTDIR(PU,1_8)
    
    IF (LLALIAZ) THEN
      DO I = NSMAX, NXPT/2-1
         PU(:,2*I+1) = ZERO
         PU(:,2*I+2) = ZERO
      END DO
    END IF
    
    !... Set size of FFT and cc object
    CALL GFT_set_fft(NX=N, NY=1, FFT=CC)

    !... Initialize Array to be Transformed
    DO J = 1, NLEV
        
       DO I = 0, NXPT/2-1
          ZK = REAL(I,8)*RPI*TWO/RLX
          PDIV(J,2*I+1)   =  ZK*PU(J,2*I+2)
          PDIV(J,2*I+2)   = -ZK*PU(J,2*I+1)
       END DO
       
       !... Backward FFT
       YU(0,1)            = PU(J,1)
       YU(NXPT/2,1)       = PU(J,2)
       YDIV(0,1)          = PDIV(J,1)
       YDIV(NXPT/2,1)     = PDIV(J,2)
       DO I = 1, NXPT/2-1
          YU(I,1)         = PU(J,2*I+1)   + ZI*PU(J,2*I+2)
          YU(NXPT-I,1)    = PU(J,2*I+1)   - ZI*PU(J,2*I+2)
          YDIV(I,1)       = PDIV(J,2*I+1) + ZI*PDIV(J,2*I+2)
          YDIV(NXPT-I,1)  = PDIV(J,2*I+1) - ZI*PDIV(J,2*I+2)
       END DO

       CALL GFT_do_fft(FFT=CC, isign=SI, scale=SC, c_in=YU, c_out=XU)
       CALL GFT_do_fft(FFT=CC, isign=SI, scale=SC, c_in=YDIV, c_out=XDIV)

       !...  FFT inverse  
       DO I = 0, NXPT-1
          PU(J,I+1)   = REAL(XU(I,1),8)
          PDIV(J,I+1) = REAL(XDIV(I,1),8) 
       END DO
       
    END DO

    !... Free FFT object
    CALL GFT_end_fft(FFT=CC)         
    
    
  END SUBROUTINE SPEC_DIV
  !*******************************************************************************!
  !*******************************************************************************!
  !*******************************************************************************!
  !*******************************************************************************!
  
END MODULE MOD_FFT

!=====================================================!
!=====================================================!
!=====================================================!
