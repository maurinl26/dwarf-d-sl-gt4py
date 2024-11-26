!######################   FFT  #######################!
!#                                                   #!
!# auteur : C.Colavolpe & F.Voitus                   #!
!# sujet  : Les transformées de Fourier              #!
!#                                                   #!
!#####################################################!

!=====================================================!
!===================  MODULE FFT  ====================!
!=====================================================!

MODULE MOD_FFT

  USE GFT
  USE MOD_PARAM, ONLY : GXPT,GYPT,NSMAX,MSMAX,RPI,RDX,RDY,RLX,RLY,&
                      & ZERO,ONE,TWO

CONTAINS

  !*****************************************************!
  !***************  LISTE DES ROUTINES  ****************!
  !*****************************************************!

  !*********  Transformée de Fourier directe  **********!
  SUBROUTINE FFTDIR_2D(XINOUT)

    IMPLICIT NONE

    REAL(8), DIMENSION(GXPT,GYPT), INTENT(INOUT) :: XINOUT

    CALL FFTDIR_X(XINOUT)
    CALL FFTDIR_Y(XINOUT)

  END SUBROUTINE FFTDIR_2D
  !*****************************************************!
  !*****************************************************!
  SUBROUTINE FFTDIR_X(XINOUT)
    
    IMPLICIT NONE

    REAL(8), DIMENSION(GXPT,GYPT), INTENT(INOUT) :: XINOUT
    INTEGER, PARAMETER                           :: SI = +1
    INTEGER                                      :: N
    INTEGER                                      :: IX,IY
    REAL(kind=GFT_Prec),PARAMETER                :: SC = 1.0_GFT_Prec
    COMPLEX(kind=GFT_Prec), DIMENSION(0:GXPT,1)  :: X,HATX
    TYPE(GFT_CC)                                 :: CC

    N = GXPT
    CALL GFT_set_fft(NX=N, NY=1, FFT=CC)

    DO IY = 1, GYPT

       DO IX = 0, GXPT-1
          X(IX,1) = XINOUT(IX+1,IY)
       END DO

       CALL GFT_do_fft(FFT=CC, isign=SI, scale=SC, c_in=X, c_out=HATX)

       XINOUT(1,IY) = DBLE(REAL(HATX(0,1)))
       XINOUT(2,IY) = DBLE(REAL(HATX(GXPT/2,1)))
       DO IX = 1, GXPT/2-1
          XINOUT(2*IX+1,IY) = DBLE(REAL(HATX(IX,1)))
          XINOUT(2*IX+2,IY) = DBLE(AIMAG(HATX(IX,1)))
       END DO

    END DO

    CALL GFT_end_fft(FFT=CC)

  END SUBROUTINE FFTDIR_X
  !*****************************************************!
   !*****************************************************!
  SUBROUTINE FFTDIR_Y(XINOUT)

    IMPLICIT NONE

    REAL(8), DIMENSION(GXPT,GYPT), INTENT(INOUT) :: XINOUT
    INTEGER, PARAMETER                           :: SI = +1
    INTEGER                                      :: N
    INTEGER                                      :: IX,IY
    REAL(kind=GFT_Prec),PARAMETER                :: SC = 1.0_GFT_Prec
    COMPLEX(kind=GFT_Prec), DIMENSION(0:GYPT,1)  :: X,HATX
    TYPE(GFT_CC)                                 :: CC

    N = GYPT
    CALL GFT_set_fft(NX=N, NY=1, FFT=CC)

    DO IX = 1, GXPT

       DO IY = 0, GYPT-1
          X(IY,1) = XINOUT(IX,IY+1)
       END DO

       CALL GFT_do_fft(FFT=CC, isign=SI, scale=SC, c_in=X, c_out=HATX)

       XINOUT(IX,1) = DBLE(REAL(HATX(0,1)))
       XINOUT(IX,2) = DBLE(REAL(HATX(GYPT/2,1)))
       DO IY = 1, GYPT/2-1
          XINOUT(IX,2*IY+1) = DBLE(REAL(HATX(IY,1)))
          XINOUT(IX,2*IY+2) = DBLE(AIMAG(HATX(IY,1)))
       END DO

    END DO

    CALL GFT_end_fft(FFT=CC)

  END SUBROUTINE FFTDIR_Y
  !*****************************************************!
  !*****************************************************!
  !*****************************************************!
  !*****************************************************!
  SUBROUTINE FFTDIR_UV(U,V)

    IMPLICIT NONE

    REAL(8), DIMENSION(GXPT,GYPT), INTENT(INOUT) :: U,V
    INTEGER, PARAMETER                           :: SI = +1
    INTEGER                                      :: N
    INTEGER                                      :: IX,IY
    REAL(kind=GFT_Prec),PARAMETER                :: SC = 1.0_GFT_Prec
    COMPLEX(kind=GFT_Prec), DIMENSION(0:GXPT,1)  :: XU,HATU,XV,HATV
    COMPLEX(kind=GFT_Prec), DIMENSION(0:GYPT,1)  :: YU,HOTU,YV,HOTV
    TYPE(GFT_CC)                                 :: CC

    N = GXPT
    CALL GFT_set_fft(NX=N, NY=1, FFT=CC)

    DO IY = 1, GYPT

       DO IX = 0, GXPT-1
          XU(IX,1) = U(IX+1,IY)
          XV(IX,1) = V(IX+1,IY)
       END DO

       CALL GFT_do_fft(FFT=CC, isign=SI, scale=SC, c_in=XU, c_out=HATU)
       CALL GFT_do_fft(FFT=CC, isign=SI, scale=SC, c_in=XV, c_out=HATV)

       U(1,IY) = DBLE(REAL(HATU(0,1)))
       U(2,IY) = DBLE(REAL(HATU(GXPT/2,1)))
       V(1,IY) = DBLE(REAL(HATV(0,1)))
       V(2,IY) = DBLE(REAL(HATV(GXPT/2,1)))
       DO IX = 1, GXPT/2-1
          U(2*IX+1,IY) = DBLE(REAL(HATU(IX,1)))
          U(2*IX+2,IY) = DBLE(AIMAG(HATU(IX,1)))
          V(2*IX+1,IY) = DBLE(REAL(HATV(IX,1)))
          V(2*IX+2,IY) = DBLE(AIMAG(HATV(IX,1)))
       END DO

    END DO

    CALL GFT_end_fft(FFT=CC)
    
    !*******************************************************************!
    !*******************************************************************!

    N = GYPT
    CALL GFT_set_fft(NX=N, NY=1, FFT=CC)

    DO IX = 1, GXPT

       DO IY = 0, GYPT-1
          YU(IY,1) = U(IX,IY+1)
          YV(IY,1) = V(IX,IY+1)
       END DO

       CALL GFT_do_fft(FFT=CC, isign=SI, scale=SC, c_in=YU, c_out=HOTU)
       CALL GFT_do_fft(FFT=CC, isign=SI, scale=SC, c_in=YV, c_out=HOTV)
       
       U(IX,1) = DBLE(REAL(HOTU(0,1)))
       U(IX,2) = DBLE(REAL(HOTU(GYPT/2,1)))
       V(IX,1) = DBLE(REAL(HOTV(0,1)))
       V(IX,2) = DBLE(REAL(HOTV(GYPT/2,1)))
       DO IY = 1, GYPT/2-1
          U(IX,2*IY+1) = DBLE(REAL(HOTU(IY,1)))
          U(IX,2*IY+2) = DBLE(AIMAG(HOTU(IY,1)))
          V(IX,2*IY+1) = DBLE(REAL(HOTV(IY,1)))
          V(IX,2*IY+2) = DBLE(AIMAG(HOTV(IY,1)))
       END DO

    END DO

    CALL GFT_end_fft(FFT=CC)
    

  END SUBROUTINE FFTDIR_UV
  !***************************************************************************!
  !***************************************************************************!
  !***************************************************************************!

  !*****************************************************!
  !*****************************************************!
  !*********  Transformée de Fourier inverse  **********!
  SUBROUTINE FFTINV_2D(XINOUT,LALIAZ)

    IMPLICIT NONE

    REAL(8), DIMENSION(GXPT,GYPT), INTENT(INOUT) :: XINOUT
    LOGICAL, OPTIONAL,             INTENT(IN)    :: LALIAZ

    LOGICAL                                      :: LL_ALIASING

    IF (PRESENT(LALIAZ)) THEN
       LL_ALIASING = LALIAZ
    ELSE
       LL_ALIASING = .FALSE.
    END IF   

    IF (LL_ALIASING) THEN
      CALL ELLIPTIC_TRUNC(XINOUT) 
    END IF

    CALL FFTINV_Y(XINOUT)
    CALL FFTINV_X(XINOUT)


  END SUBROUTINE FFTINV_2D
  !*****************************************************!
  !*****************************************************!
  SUBROUTINE FFTINV_X(XINOUT)

    IMPLICIT NONE

    REAL(8), DIMENSION(GXPT,GYPT), INTENT(INOUT) :: XINOUT
    INTEGER, PARAMETER                           :: SI=-1
    INTEGER                                      :: N
    INTEGER                                      :: IX,IY
    REAL(kind=GFT_Prec)                          :: SC
    COMPLEX(kind=GFT_Prec), DIMENSION(0:GXPT,1)  :: HATX,X
    COMPLEX(8), PARAMETER                        :: ZI = (0.d0,1.d0)
    TYPE(GFT_CC)                                 :: CC

    N = GXPT
    SC = 1.0_GFT_Prec/REAL(GXPT, kind=GFT_Prec)
    CALL GFT_set_fft(NX=N, NY=1, FFT=CC)

    DO IY = 1, GYPT

       HATX(0,1)      = XINOUT(1,IY)
       HATX(GXPT/2,1) = XINOUT(2,IY)
       DO IX = 1, GXPT/2-1
          HATX(IX,1)      = XINOUT(2*IX+1,IY) + ZI*XINOUT(2*IX+2,IY)
          HATX(GXPT-IX,1) = XINOUT(2*IX+1,IY) - ZI*XINOUT(2*IX+2,IY)
       END DO

       CALL GFT_do_fft(FFT=CC, isign=SI, scale=SC, c_in=HATX, c_out=X)

       DO IX = 0, GXPT-1
          XINOUT(IX+1,IY) = REAL(X(IX,1),8)
       END DO

    END DO

    CALL GFT_end_fft(FFT=CC)

  END SUBROUTINE FFTINV_X
  !******************************************************!
  !******************************************************!
  SUBROUTINE FFTINV_Y(XINOUT)

    IMPLICIT NONE

    REAL(8), DIMENSION(GXPT,GYPT), INTENT(INOUT) :: XINOUT
    INTEGER, PARAMETER                           :: SI=-1
    INTEGER                                      :: N
    INTEGER                                      :: IX,IY
    REAL(kind=GFT_Prec)                          :: SC
    COMPLEX(kind=GFT_Prec), DIMENSION(0:GYPT,1)  :: HATX,X
    COMPLEX(8), PARAMETER                        :: ZI = (0.d0,1.d0)
    TYPE(GFT_CC)                                 :: CC

    N = GYPT
    SC = 1.0_GFT_Prec/REAL(GYPT, kind=GFT_Prec)
    CALL GFT_set_fft(NX=N, NY=1, FFT=CC)

    DO IX = 1,GXPT

       HATX(0,1)      = XINOUT(IX,1)
       HATX(GYPT/2,1) = XINOUT(IX,2)
       DO IY = 1, GYPT/2-1
          HATX(IY,1)      = XINOUT(IX,2*IY+1) + ZI*XINOUT(IX,2*IY+2)
          HATX(GYPT-IY,1) = XINOUT(IX,2*IY+1) - ZI*XINOUT(IX,2*IY+2)
       END DO

       CALL GFT_do_fft(FFT=CC, isign=SI, scale=SC, c_in=HATX, c_out=X)

       DO IY = 0, GYPT-1
          XINOUT(IX,IY+1) = REAL(X(IY,1),8)
       END DO

    END DO

    CALL GFT_end_fft(FFT=CC)

  END SUBROUTINE FFTINV_Y
  !*************************************************************************!
  !*************************************************************************!
  !*************************************************************************!
  !*************************************************************************!
  SUBROUTINE FFTINV_UV(U,V)

    IMPLICIT NONE

    REAL(8), DIMENSION(GXPT,GYPT), INTENT(INOUT) :: U,V
    INTEGER, PARAMETER                           :: SI=-1
    INTEGER                                      :: N
    INTEGER                                      :: IX,IY
    REAL(kind=GFT_Prec)                          :: SC
    COMPLEX(kind=GFT_Prec), DIMENSION(0:GXPT,1)  :: HATU,XU,HATV,XV
    COMPLEX(kind=GFT_Prec), DIMENSION(0:GYPT,1)  :: HOTU,YU,HOTV,YV
    COMPLEX(8), PARAMETER                        :: ZI = (0.d0,1.d0)
    TYPE(GFT_CC)                                 :: CC

    
    N  = GXPT
    SC = 1.0_GFT_Prec/REAL(GXPT, kind=GFT_Prec)
    CALL GFT_set_fft(NX=N, NY=1, FFT=CC)

    DO IY = 1, GYPT

       HATU(0,1)      = U(1,IY)
       HATU(GXPT/2,1) = U(2,IY)
       HATV(0,1)      = V(1,IY)
       HATV(GXPT/2,1) = V(2,IY)
       
       DO IX = 1, GXPT/2-1
          HATU(IX,1)      = U(2*IX+1,IY) + ZI*U(2*IX+2,IY)
          HATU(GXPT-IX,1) = U(2*IX+1,IY) - ZI*U(2*IX+2,IY)
          HATV(IX,1)      = V(2*IX+1,IY) + ZI*V(2*IX+2,IY)
          HATV(GXPT-IX,1) = V(2*IX+1,IY) - ZI*V(2*IX+2,IY)
       END DO

       CALL GFT_do_fft(FFT=CC, isign=SI, scale=SC, c_in=HATU, c_out=XU)
       CALL GFT_do_fft(FFT=CC, isign=SI, scale=SC, c_in=HATV, c_out=XV)
       
       DO IX = 0, GXPT-1
          U(IX+1,IY) = REAL(XU(IX,1),8)
          V(IX+1,IY) = REAL(XV(IX,1),8)
       END DO

    END DO

    CALL GFT_end_fft(FFT=CC)

    !*******************************************************************

    N  = GYPT
    SC = 1.0_GFT_Prec/REAL(GYPT, kind=GFT_Prec)
    CALL GFT_set_fft(NX=N, NY=1, FFT=CC)

    DO IX = 1,GXPT

       HOTU(0,1)      = U(IX,1)
       HOTU(GYPT/2,1) = U(IX,2)
       HOTV(0,1)      = V(IX,1)
       HOTV(GYPT/2,1) = V(IX,2)
       DO IY = 1, GYPT/2-1
          HOTU(IY,1)      = U(IX,2*IY+1) + ZI*U(IX,2*IY+2)
          HOTU(GYPT-IY,1) = U(IX,2*IY+1) - ZI*U(IX,2*IY+2)
          HOTV(IY,1)      = V(IX,2*IY+1) + ZI*V(IX,2*IY+2)
          HOTV(GYPT-IY,1) = V(IX,2*IY+1) - ZI*V(IX,2*IY+2)
       END DO

       CALL GFT_do_fft(FFT=CC, isign=SI, scale=SC, c_in=HOTU, c_out=YU)
       CALL GFT_do_fft(FFT=CC, isign=SI, scale=SC, c_in=HOTV, c_out=YV)

       DO IY = 0, GYPT-1
          U(IX,IY+1) = REAL(YU(IY,1),8)
          V(IX,IY+1) = REAL(YV(IY,1),8)
       END DO

    END DO

    CALL GFT_end_fft(FFT=CC)

  END SUBROUTINE FFTINV_UV
  !************************************************************************!
  !************************************************************************!
  !************************************************************************!

  
  !**************************************************************!
  !**************************************************************!
  !**************************************************************!
  SUBROUTINE FFTDER_DXDY(P,DPDX,DPDY,LALIAZ)

    IMPLICIT NONE

    REAL(8), DIMENSION(GXPT,GYPT), INTENT(INOUT) :: P
    REAL(8), DIMENSION(GXPT,GYPT), INTENT(OUT)   :: DPDX, DPDY
    LOGICAL, OPTIONAL,             INTENT(IN)    :: LALIAZ

    LOGICAL                                      :: LL_ALIASING

    IF (PRESENT(LALIAZ)) THEN
       LL_ALIASING = LALIAZ
    ELSE
       LL_ALIASING = .FALSE.
    END IF   

    IF (LL_ALIASING) THEN
      CALL ELLIPTIC_TRUNC(P) 
    END IF
   
    CALL FFTDER_DX(DPDX,P)
    CALL FFTDER_DY(DPDY,P)

    CALL FFTINV_2D(P)
    CALL FFTINV_2D(DPDX)
    CALL FFTINV_2D(DPDY)
    
  END SUBROUTINE FFTDER_DXDY
  !*****************************************************************!
  !*****************************************************************!
  !*****************************************************************!
  SUBROUTINE FFTDER_DX(PDXDX,PX)

  IMPLICIT NONE
  
  REAL(8), DIMENSION(GXPT,GYPT), INTENT(OUT) :: PDXDX
  REAL(8), DIMENSION(GXPT,GYPT), INTENT(IN)  :: PX
  INTEGER                                    :: I
  REAL(8)                                    :: ZK

  DO I = 0, GXPT/2-1
     ZK = REAL(I,8)*RPI*TWO/RLX
     PDXDX(2*I+1,:) =  ZK*PX(2*I+2,:)
     PDXDX(2*I+2,:) = -ZK*PX(2*I+1,:)
  END DO

  END SUBROUTINE FFTDER_DX
  !******************************************************!
  !******************************************************!
  SUBROUTINE FFTDER_DY(PDXDY,PY)

  IMPLICIT NONE
  
  REAL(8), DIMENSION(GXPT,GYPT), INTENT(OUT) :: PDXDY
  REAL(8), DIMENSION(GXPT,GYPT), INTENT(IN)  :: PY
  INTEGER                                    :: J
  REAL(8)                                    :: ZK

  DO J = 0, GYPT/2-1
     ZK = REAL(J,8)*RPI*TWO/RLY
     PDXDY(:,2*J+1) =  ZK*PY(:,2*J+2)
     PDXDY(:,2*J+2) = -ZK*PY(:,2*J+1)
  END DO

  END SUBROUTINE FFTDER_DY
  !******************************************************!
  !******************************************************!
  !******************************************************!
  
  !*****************************************************!
  !*****************************************************!
  SUBROUTINE ELLIPTIC_TRUNC(P)

    IMPLICIT NONE

    REAL(8), DIMENSION(GXPT,GYPT), INTENT(INOUT) :: P

    INTEGER   :: MY_TRUNC, NX
   
    P(NSMAX+1:GXPT,:) = ZERO
    P(:,MSMAX+1:GYPT) = ZERO

    DO NX=1,NSMAX
       MY_TRUNC = INT(MSMAX*DSQRT(ONE-((REAL(NX,8)/REAL(NSMAX,8))**2)))
       P(NX,MY_TRUNC+1:MSMAX) = ZERO
    END DO    
    
  END SUBROUTINE ELLIPTIC_TRUNC
  !****************************************************************************!
  !****************************************************************************!
  !****************************************************************************!
  SUBROUTINE FFTINV_TRUNC(P,LALIAZ)

    IMPLICIT NONE

    REAL(8), DIMENSION(GXPT,GYPT), INTENT(INOUT) :: P
    LOGICAL, OPTIONAL,             INTENT(IN)    :: LALIAZ

    LOGICAL                                      :: LL_ALIASING

    IF (PRESENT(LALIAZ)) THEN
       LL_ALIASING = LALIAZ
    ELSE
       LL_ALIASING = .FALSE.
    END IF   

    IF (LL_ALIASING) THEN
      CALL ELLIPTIC_TRUNC(P) 
    END IF

    CALL FFTINV_2D(P)
    
  END SUBROUTINE FFTINV_TRUNC
  !****************************************************************************!
  !****************************************************************************!
  !****************************************************************************!
  !****************************************************************************!
  SUBROUTINE FFTDIV(U,V,DIV,LALIAZ)

    IMPLICIT NONE

    REAL(8), DIMENSION(GXPT,GYPT), INTENT(INOUT) :: U,V
    REAL(8), DIMENSION(GXPT,GYPT), INTENT(INOUT) :: DIV
    LOGICAL, OPTIONAL,             INTENT(IN)    :: LALIAZ

    INTEGER(8)                                   :: NX, NY
    REAL(8)                                      :: ZX, ZY
    LOGICAL                                      :: LL_ALIASING

    IF (PRESENT(LALIAZ)) THEN
       LL_ALIASING = LALIAZ
    ELSE
       LL_ALIASING = .FALSE.
    END IF   
    
    IF (LL_ALIASING) THEN
       CALL ELLIPTIC_TRUNC(U)
       CALL ELLIPTIC_TRUNC(V)
    END IF
    
    DO NX = 0,GXPT/2-1
       ZX = REAL(NX,8)*RPI*TWO/RLX
       DO NY = 0,GYPT/2-1
          ZY = REAL(NY,8)*RPI*TWO/RLY
          
          DIV(2*NX+1,2*NY+1) = (+ZX)*U(2*NX+2,2*NY+1)+(+ZY)*V(2*NX+1,2*NY+2) 
          DIV(2*NX+2,2*NY+1) = (-ZX)*U(2*NX+1,2*NY+1)+(+ZY)*V(2*NX+2,2*NY+2) 
          DIV(2*NX+1,2*NY+2) = (+ZX)*U(2*NX+2,2*NY+2)+(-ZY)*V(2*NX+1,2*NY+1)     
          DIV(2*NX+2,2*NY+2) = (-ZX)*U(2*NX+1,2*NY+2)+(-ZY)*V(2*NX+2,2*NY+1)
          
       END DO
    END DO   
   
    CALL FFTINV_UV(U,V) 
    CALL FFTINV_2D(DIV)

  END SUBROUTINE FFTDIV
  !***************************************************************************!
  !***************************************************************************!
  !***************************************************************************!
  !***************************************************************************!
  
END MODULE MOD_FFT

!=====================================================!
!=====================================================!
!=====================================================!
