!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! -*- Mode: F90 -*- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! GFT_do.f90 --- GFT forward/backward FFT routines
!! 
!! Auteur          : Jalel Chergui (CNRS/IDRIS) <Jalel.Chergui@idris.fr>
!! Créé le         : Tue Feb 19 10:25:04 2002
!! Dern. mod. par  : Jalel Chergui (CNRS/IDRIS) <Jalel.Chergui@idris.fr>
!! Dern. mod. le   : Thu Jul 11 17:07:50 2002
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Permission is granted to copy and distribute this file or modified
! versions of this file for no fee, provided the copyright notice and
! this permission notice are preserved on all copies.
! Copyright (C) Februry 2002, CNRS/IDRIS, Jalel.Chergui@idris.fr.
!
MODULE GFT_do
  USE GFT_common
  USE GFT_JMFFT

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: GFT_do_fft
  INTERFACE GFT_do_fft
     MODULE PROCEDURE GFT_do_cc_1d,  GFT_do_cc_2d,  GFT_do_cc_3d, &
                      GFT_do_cr_1d,  GFT_do_cr_2d,  GFT_do_cr_3d, &
                      GFT_do_rc_1d,  GFT_do_rc_2d,  GFT_do_rc_3d, &
                      GFT_do_mcc_1d, GFT_do_mcr_1d, GFT_do_mrc_1d
  END INTERFACE GFT_do_fft

  CONTAINS

  SUBROUTINE GFT_do_cc_1d(FFT, isign, scale, c_in, c_out, code)

  IMPLICIT NONE

  !... Input dummy Parameters
  TYPE(GFT_CC), INTENT(inout)                       :: FFT
  INTEGER, INTENT(in)                               :: isign
  REAL(kind=GFT_prec), INTENT(in)                   :: scale
  COMPLEX(kind=GFT_prec), INTENT(in), DIMENSION(0:) :: c_in

  !... Output dummy Parameters
  COMPLEX(kind=GFT_prec), INTENT(out), DIMENSION(0:) :: c_out
  INTEGER, OPTIONAL, INTENT(out)                     :: code

  !... Local variables
  REAL(kind=GFT_prec), DIMENSION(0:2*SIZE(c_in)-1)  :: x
  REAL(kind=GFT_prec), DIMENSION(0:2*SIZE(c_out)-1) :: y
  INTEGER                                        :: i, ioff, nfact
  INTEGER, DIMENSION(0:99)                       :: fact
  CHARACTER(len=*), PARAMETER                    :: spname="GFT_do_cc_1d"

  !... Check conditions
  IF (.NOT.FFT%Init) CALL GFT_error(spname//": Error: GFT_init routine has not been called prior to this routine")
  IF ( SIZE(c_in) < FFT%Nx) CALL GFT_error(spname//":  Error: Size(C_IN) < Nx")
  IF ( SIZE(c_out) < FFT%Nx) CALL GFT_error(spname//":  Error: Size(C_OUT) < Nx")
  IF (isign /=-1 .AND. isign /= 1) CALL GFT_error(spname//": Error: sign has invalid value")

  FFT%Isign = isign
  FFT%Scale = scale

!dir$ ivdep
!ocl novrec
!cdir nodep
  DO i = 0,FFT%Nx-1
    x(2*i)   = REAL(c_in(i),KIND=GFT_prec)
    x(2*i+1) = AIMAG(c_in(i))
  END DO
  nfact = NINT(FFT%Table(0))
  fact(0:nfact-1) = NINT(FFT%Table(0:nfact-1))

  ! On copie le tableau d'entree dans le tableau de travail
  ! On en profite pour premultiplier et pour tenir compte du signe
  DO i = 0,FFT%Nx-1
    FFT%Work(i)       =       scale* x(2*i)
    FFT%Work(FFT%Nx+i) = isign*scale* x(2*i+1)
  END DO
  ioff = 0

  ! On appelle le sous-programme principal
  CALL jmccm1d(1,FFT%Nx,fact,100,0,FFT%Table,FFT%TableSize,100+0,FFT%Work,FFT%WorkSize,ioff)

  ! On recopie dans le tableau d'arrivee
  DO i = 0,FFT%Nx-1
    y(2*i)   =         FFT%Work(ioff  +i)
    y(2*i+1) = isign * FFT%Work(ioff+FFT%Nx+i)
  END DO

!dir$ ivdep
!ocl novrec
!cdir nodep
  DO i = 0,FFT%Nx-1
     c_out(i) = CMPLX(y(2*i),y(2*i+1),KIND=GFT_prec)
  END DO

  IF( PRESENT(code) ) code=0
END SUBROUTINE GFT_do_cc_1d

SUBROUTINE GFT_do_cc_2d(FFT, isign, scale, c_in, c_out, code)
  IMPLICIT NONE

  !... Input dummy Parameters
  TYPE(GFT_CC), INTENT(inout)                          :: FFT
  INTEGER, INTENT(in)                                  :: isign
  REAL(kind=GFT_prec), INTENT(in)                      :: scale
  COMPLEX(kind=GFT_prec), INTENT(in), DIMENSION(0:,0:) :: c_in

  !... Output dummy Parameters
  COMPLEX(kind=GFT_prec), INTENT(out), DIMENSION(0:,0:) :: c_out
  INTEGER, OPTIONAL, INTENT(out)                        :: code

  !... Local variables
  REAL(kind=GFT_prec), DIMENSION(0:2*SIZE(c_in)-1)  :: x
  REAL(kind=GFT_prec), DIMENSION(0:2*SIZE(c_out)-1) :: y
  INTEGER                     :: i, j, ioff
  INTEGER                     :: nfact, mfact
  INTEGER, DIMENSION(0:99)    :: fact
  INTEGER                     :: ideb, ifin, jdeb, jfin, n_temp, m_temp, nwork_temp
  LOGICAL                     :: debut, fin
  CHARACTER(len=*), PARAMETER :: spname="GFT_do_cc_2d"

  !... Check conditions
  IF (.NOT.FFT%Init) CALL GFT_error(spname//": Error: GFT_init must be called prior to this routine")
  FFT%Ldx1  = SIZE(c_in ,dim=1)
  FFT%Ldy1  = SIZE(c_out,dim=1)
  IF ( SIZE(c_in) < FFT%Nx*FFT%Ny ) CALL GFT_error(spname//":  Error: Size(C_IN) < Nx*Ny")
  IF ( SIZE(c_out) < FFT%Nx*FFT%Ny ) CALL GFT_error(spname//":  Error: Size(C_OUT) < Nx*Ny")
  IF ( FFT%Ldx1 < FFT%Nx ) CALL GFT_error(spname//":  Error: Size(C_IN,dim=1) < Nx")
  IF ( FFT%Ldy1 < FFT%Nx ) CALL GFT_error(spname//":  Error: Size(C_OUT,dim=1) < Nx")
  IF (isign /=-1 .AND. isign /= 1) CALL GFT_error(spname//": Error: isign value is invalid")

  FFT%Isign = isign
  FFT%Scale = scale

  DO j = 0, FFT%Ny-1
!dir$ ivdep
!ocl novrec
!cdir nodep
    DO i = 0,FFT%Nx-1
      x(2*i  +j*2*FFT%Ldx1) = REAL(c_in(i,j),KIND=GFT_prec)
      x(2*i+1+j*2*FFT%Ldx1) = AIMAG(c_in(i,j))
    END DO
  END DO

  nfact = NINT(FFT%Table(0))
  mfact = NINT(FFT%Table(nfact)) + nfact
  fact(0:mfact-1) = NINT(FFT%Table(0:mfact-1))

  ! On fait les T.F. sur la premiere dimension en tronconnant sur la deuxieme
  debut = .TRUE.
  DO

    ! Tronconnage
    ! Note : on met npair a .true. car il n'y a pas de restriction dans ce cas
    CALL jmdecoup(FFT%Ny,4*FFT%Nx,FFT%WorkSize,debut,.TRUE.,m_temp,jdeb,jfin,nwork_temp,fin)

    ! On copie le tableau d'entree dans le tableau de travail
    ! On en profite pour premultiplier et pour tenir compte du signe
    ! Note : On copie en transposant
    DO i = 0,FFT%Nx-1
!dir$ ivdep
!ocl novrec
!cdir nodep
      DO j = jdeb,jfin
        FFT%Work(j-jdeb+m_temp*i)         =       scale*x(2*i  +j*2*FFT%Ldx1)
        FFT%Work(j-jdeb+m_temp*(FFT%Nx+i)) = isign*scale*x(2*i+1+j*2*FFT%Ldx1)
      END DO
    END DO
    ioff = 0

    ! Attention : ioff1 est peut-etre modifie en sortie
    CALL jmccm1d(m_temp,FFT%Nx,fact,100,0    ,FFT%Table,FFT%TableSize,100+0  ,FFT%Work,nwork_temp,ioff)

    ! On recopie dans le tableau d'arrivee
    DO i = 0,FFT%Nx-1
!dir$ ivdep
!ocl novrec
!cdir nodep
      DO j = jdeb,jfin
        y(2*i  +j*2*FFT%Ldy1) = FFT%Work(ioff+j-jdeb+m_temp*i)
        y(2*i+1+j*2*FFT%Ldy1) = FFT%Work(ioff+j-jdeb+m_temp*(FFT%Nx+i))
      END DO
    END DO

    ! A-t-on fini ?
    IF (fin) THEN
      EXIT
    ELSE
      debut = .FALSE.
      CYCLE
    END IF

  END DO

  ! On fait les T.F. sur l'autre dimension
  debut = .TRUE.
  DO

    ! Tronconnage
    CALL jmdecoup(FFT%Nx,4*FFT%Ny,FFT%WorkSize,debut,.TRUE.,n_temp,ideb,ifin,nwork_temp,fin)

    ! On copie
    DO j = 0,FFT%Ny-1
!dir$ ivdep
!ocl novrec
!cdir nodep
      DO i = ideb,ifin
        FFT%Work(i-ideb+n_temp*j)         = y(2*i  +j*2*FFT%Ldy1)
        FFT%Work(i-ideb+n_temp*(FFT%Ny+j)) = y(2*i+1+j*2*FFT%Ldy1)
      END DO
    END DO
    ioff = 0

    CALL jmccm1d(n_temp,FFT%Ny,fact,100,nfact,FFT%Table,FFT%TableSize,100+2*FFT%Nx,FFT%Work,nwork_temp,ioff)

    ! On recopie dans le tableau d'arrivee
    DO j = 0,FFT%Ny-1
!dir$ ivdep
!ocl novrec
!cdir nodep
      DO i = ideb,ifin
        y(2*i  +j*2*FFT%Ldy1) =       FFT%Work(ioff+i-ideb+n_temp*j)
        y(2*i+1+j*2*FFT%Ldy1) = isign*FFT%Work(ioff+i-ideb+n_temp*(FFT%Ny+j))
      END DO
    END DO

    ! A-t-on fini ?
    IF (fin) THEN
      EXIT
    ELSE
      debut = .FALSE.
      CYCLE
    END IF

  END DO

  DO j = 0,FFT%Ny-1
!dir$ ivdep
!ocl novrec
!cdir nodep
     DO i = 0,FFT%Nx-1
        c_out(i,j) = CMPLX(y(2*i+j*2*FFT%Ldy1), y(2*i+1+j*2*FFT%Ldy1),KIND=GFT_prec)
     END DO
  END DO

  IF( PRESENT(code) ) code=0
END SUBROUTINE GFT_do_cc_2d

SUBROUTINE GFT_do_cc_3d(FFT, isign, scale, c_in, c_out, code)
  IMPLICIT NONE

  !... Input dummy Parameters
  TYPE(GFT_CC), INTENT(inout)                             :: FFT
  INTEGER, INTENT(in)                                     :: isign
  REAL(kind=GFT_prec), INTENT(in)                         :: scale
  COMPLEX(kind=GFT_prec), INTENT(in), DIMENSION(0:,0:,0:) :: c_in

  !... Output dummy Parameters
  COMPLEX(kind=GFT_prec), INTENT(out), DIMENSION(0:,0:,0:) :: c_out
  INTEGER, OPTIONAL, INTENT(out)                           :: code

  !... Local variables
  INTEGER                     :: i, j, k, ioff
  INTEGER                     :: nfact, mfact, lfact
  INTEGER, DIMENSION(0:99)    :: fact
  INTEGER                     :: ideb, ifin, i1, i2, jdeb, jfin, j1, j2, kdeb, kfin
  INTEGER                     :: nwork_temp, nmtemp, mltemp, nltemp, iwork
  LOGICAL                     :: debut, fini
  CHARACTER(len=*), PARAMETER :: spname="GFT_do_cc_3d"
  REAL(kind=8), DIMENSION(0:2*SIZE(c_in)-1)  :: x
  REAL(kind=8), DIMENSION(0:2*SIZE(c_out)-1) :: y

  !... Check conditions
  IF (.NOT.FFT%Init) CALL GFT_error(spname//": Error: GFT_init must be called prior to this routine")
  FFT%Ldx1  = SIZE(c_in,dim=1)
  FFT%Ldx2  = SIZE(c_in,dim=2)
  FFT%Ldy1  = SIZE(c_out,dim=1)
  FFT%Ldy2  = SIZE(c_out,dim=2)
  IF ( SIZE(c_in) < FFT%Nx*FFT%Ny*FFT%Nz) CALL GFT_error(spname//":  Error: Size(C_IN) < Nx*Ny*Nz")
  IF ( SIZE(c_out) < FFT%Nx*FFT%Ny*FFT%Nz) CALL GFT_error(spname//":  Error: Size(C_OUT) < Nx*Ny*Nz")
  IF ( FFT%Ldx1 < FFT%Nx ) CALL GFT_error(spname//":  Error: size(C_IN,dim=1) < Nx")
  IF ( FFT%Ldx2 < FFT%Ny ) CALL GFT_error(spname//":  Error: size(C_IN,dim=2) < Ny")
  IF ( FFT%Ldy1 < FFT%Nx ) CALL GFT_error(spname//":  Error: size(C_OUT,dim=1) < Nx")
  IF ( FFT%Ldy2 < FFT%Ny ) CALL GFT_error(spname//":  Error: size(C_OUT,dim=2) < Ny")
  IF (isign /=-1 .AND. isign /= 1) CALL GFT_error(spname//": Error: sign has an invalid value")


  FFT%Isign = isign
  FFT%Scale = scale

  DO k = 0, FFT%Nz-1
    DO j = 0, FFT%Ny-1
!dir$ ivdep
!ocl novrec
!cdir nodep
      DO i = 0,FFT%Nx-1
        x(2*i  +2*FFT%Ldx1*j+2*FFT%Ldx1*FFT%Ldx2*k) = REAL(c_in(i,j,k),KIND=GFT_prec)
        x(2*i+1+2*FFT%Ldx1*j+2*FFT%Ldx1*FFT%Ldx2*k) = AIMAG(c_in(i,j,k))
      END DO
    END DO
  END DO

  nfact           = NINT(FFT%Table(0))
  mfact           = NINT(FFT%Table(nfact)) + nfact
  lfact           = NINT(FFT%Table(mfact)) + mfact
  fact(0:lfact-1) = NINT(FFT%Table(0:lfact-1))

  ! On fait les T.F. sur la troisieme dimension en tronconnant sur la premiere
  ! et la deuxieme
  debut = .TRUE.
  fini  = .FALSE.
  DO WHILE (.NOT.fini)

    ! Tronconnage
    ! Note : on met npair a .true. car il n'y a pas de restriction dans ce cas
    CALL jmdecoup3(FFT%Nx,FFT%Ny,4*FFT%Nz,FFT%WorkSize,debut,.TRUE.,ideb,ifin,jdeb,jfin,nmtemp,nwork_temp,fini)
    debut = .FALSE.

    ! On copie le tableau d'entree dans le tableau de travail
    ! On en profite pour premultiplier et pour tenir compte du signe
    ! On prend garde a la gestion des extremites
    DO k = 0,FFT%Nz-1
      iwork = 0
      DO j = jdeb,jfin
        i1 = 0
        i2 = FFT%Nx-1
        IF (j == jdeb) i1 = ideb
        IF (j == jfin) i2 = ifin
!dir$ ivdep
!ocl novrec
!cdir nodep
        DO i = i1,i2
          FFT%Work(             iwork+k*nmtemp) =       scale*x(2*i  +2*FFT%Ldx1*j+2*FFT%Ldx1*FFT%Ldx2*k)
          FFT%Work(nwork_temp/4+iwork+k*nmtemp) = isign*scale*x(2*i+1+2*FFT%Ldx1*j+2*FFT%Ldx1*FFT%Ldx2*k)
          iwork = iwork+1
        END DO
      END DO
    END DO

    ! On fait les T.F. sur la troisieme dimension
    ioff = 0
    CALL jmccm1d(nmtemp,FFT%Nz,fact,100,mfact,FFT%Table,FFT%TableSize,100+2*(FFT%Nx+FFT%Ny),FFT%Work,nwork_temp,ioff)

    ! On recopie dans le tableau d'arrivee
    DO k = 0,FFT%Nz-1
      iwork = 0
      DO j = jdeb,jfin
        i1 = 0
        i2 = FFT%Nx-1
        IF (j == jdeb) i1 = ideb
        IF (j == jfin) i2 = ifin
!dir$ ivdep
!ocl novrec
!cdir nodep
        DO i = i1,i2
          y(2*i  +2*FFT%Ldy1*j+2*FFT%Ldy1*FFT%Ldy2*k) = FFT%Work(ioff             +iwork+k*nmtemp)
          y(2*i+1+2*FFT%Ldy1*j+2*FFT%Ldy1*FFT%Ldy2*k) = FFT%Work(ioff+nwork_temp/4+iwork+k*nmtemp)
          iwork = iwork+1
        END DO
      END DO
    END DO

  END DO

  ! On fait les T.F. sur la deuxieme dimension en tronconnant sur la premiere
  ! et la troisieme
  debut = .TRUE.
  fini  = .FALSE.
  DO WHILE (.NOT.fini)

    ! Tronconnage
    CALL jmdecoup3(FFT%Nx,FFT%Nz,4*FFT%Ny,FFT%WorkSize,debut,.TRUE.,ideb,ifin,kdeb,kfin,nltemp,nwork_temp,fini)
    debut = .FALSE.

    ! On copie le tableau d'entree dans le tableau de travail
    ! On prend garde a la gestion des extremites
    DO j = 0,FFT%Ny-1
      iwork = 0
      DO k = kdeb,kfin
        i1 = 0
        i2 = FFT%Nx-1
        IF (k == kdeb) i1 = ideb
        IF (k == kfin) i2 = ifin
!dir$ ivdep
!ocl novrec
!cdir nodep
        DO i = i1,i2
          FFT%Work(             iwork+j*nltemp) = y(2*i  +2*FFT%Ldy1*j+2*FFT%Ldy1*FFT%Ldy2*k)
          FFT%Work(nwork_temp/4+iwork+j*nltemp) = y(2*i+1+2*FFT%Ldy1*j+2*FFT%Ldy1*FFT%Ldy2*k)
          iwork = iwork+1
        END DO
      END DO
    END DO

    ! On fait les T.F. sur la deuxieme dimension
    ioff = 0
    CALL jmccm1d(nltemp,FFT%Ny,fact,100,nfact,FFT%table,FFT%TableSize,100+2*FFT%Nx    ,FFT%Work,nwork_temp,ioff)

    ! On recopie dans le tableau d'arrivee
    DO j = 0,FFT%Ny-1
      iwork = 0
      DO k = kdeb,kfin
        i1 = 0
        i2 = FFT%Nx-1
        IF (k == kdeb) i1 = ideb
        IF (k == kfin) i2 = ifin
!dir$ ivdep
!ocl novrec
!cdir nodep
        DO i = i1,i2
          y(2*i  +2*FFT%Ldy1*j+2*FFT%Ldy1*FFT%Ldy2*k) = FFT%Work(ioff             +iwork+j*nltemp)
          y(2*i+1+2*FFT%Ldy1*j+2*FFT%Ldy1*FFT%Ldy2*k) = FFT%Work(ioff+nwork_temp/4+iwork+j*nltemp)
          iwork = iwork+1
        END DO
      END DO
    END DO

  END DO

  ! On fait les T.F. sur la premiere dimension en tronconnant sur la deuxieme
  ! et la troisieme
  debut = .TRUE.
  fini  = .FALSE.
  DO WHILE (.NOT.fini)

    ! Tronconnage
    CALL jmdecoup3(FFT%Ny,FFT%Nz,4*FFT%Nx,FFT%WorkSize,debut,.TRUE.,jdeb,jfin,kdeb,kfin,mltemp,nwork_temp,fini)
    debut = .FALSE.

    ! On copie le tableau d'entree dans le tableau de travail
    ! On prend garde a la gestion des extremites
    DO i = 0,FFT%Nx-1
      iwork = 0
      DO k = kdeb,kfin
        j1 = 0
        j2 = FFT%Ny-1
        IF (k == kdeb) j1 = jdeb
        IF (k == kfin) j2 = jfin
!dir$ ivdep
!ocl novrec
!cdir nodep
        DO j = j1,j2
          FFT%Work(             iwork+i*mltemp) = y(2*i  +2*FFT%Ldy1*j+2*FFT%Ldy1*FFT%Ldy2*k)
          FFT%Work(nwork_temp/4+iwork+i*mltemp) = y(2*i+1+2*FFT%Ldy1*j+2*FFT%Ldy1*FFT%Ldy2*k)
          iwork = iwork+1
        END DO
      END DO
    END DO

    ! On fait les T.F. sur la premiere dimension
    ioff = 0
    CALL jmccm1d(mltemp,FFT%Nx,fact,100,0    ,FFT%Table,FFT%TableSize,100+0      ,FFT%Work,nwork_temp,ioff)

    ! On recopie dans le tableau d'arrivee en tenant compte du signe
    DO i = 0,FFT%Nx-1
      iwork = 0
      DO k = kdeb,kfin
        j1 = 0
        j2 = FFT%Ny-1
        IF (k == kdeb) j1 = jdeb
        IF (k == kfin) j2 = jfin
!dir$ ivdep
!ocl novrec
!cdir nodep
        DO j = j1,j2
          y(2*i  +2*FFT%Ldy1*j+2*FFT%Ldy1*FFT%Ldy2*k) =       FFT%Work(ioff             +iwork+i*mltemp)
          y(2*i+1+2*FFT%Ldy1*j+2*FFT%Ldy1*FFT%Ldy2*k) = isign*FFT%Work(ioff+nwork_temp/4+iwork+i*mltemp)
          iwork = iwork+1
        END DO
      END DO
    END DO

  END DO


  DO k = 0, FFT%Nz-1
    DO j = 0, FFT%Ny-1
!dir$ ivdep
!ocl novrec
!cdir nodep
      DO i = 0, FFT%Nx-1 
         c_out(i,j,k) = CMPLX(y(2*i+2*FFT%Ldy1*j+2*FFT%Ldy1*FFT%Ldy2*k),y(2*i+1+2*FFT%Ldy1*j+2*FFT%Ldy1*FFT%Ldy2*k),KIND=GFT_prec)
      END DO
    END DO
  END DO

  IF( PRESENT(code) ) code = 0
END SUBROUTINE GFT_do_cc_3d

SUBROUTINE GFT_do_cr_1d(FFT, isign, scale, c_in, r_out, code)
  IMPLICIT NONE

  !... Input dummy Parameters
  TYPE(GFT_RCR), INTENT(inout)                      :: FFT
  INTEGER, INTENT(in)                               :: isign
  REAL(kind=GFT_prec), INTENT(in)                   :: scale
  COMPLEX(kind=GFT_prec), INTENT(in), DIMENSION(0:) :: c_in

  !... Output dummy Parameters
  REAL(kind=GFT_prec), INTENT(out), DIMENSION(0:) :: r_out
  INTEGER, OPTIONAL, INTENT(out)                  :: code

  ! Local variables
  REAL(kind=GFT_prec), DIMENSION(0:2*(SIZE(r_out)/2)+1) :: x
  INTEGER                             :: i, nfact
  INTEGER, DIMENSION(0:99)            :: fact
  INTEGER                             :: dimx, debx, incx, jumpx
  INTEGER                             :: dimy, deby, incy, jumpy
  CHARACTER(len=*), PARAMETER         :: spname = "GFT_do_cr_1d"

  !... Check conditions
  IF (.NOT.FFT%Init) CALL GFT_error(spname//": Error: GFT_init must be called prior to this routine")
  FFT%Ldx1 = SIZE(c_in)
  FFT%Ldy1 = SIZE(r_out)
  IF ( FFT%Ldx1 < FFT%Nx/2+1) CALL GFT_error(spname//":  Error: Size(C_IN) < Nx/2+1")
  IF ( FFT%Ldy1 < FFT%Nx) CALL GFT_error(spname//":  Error: Size(R_OUT) < Nx")
  IF (isign /=-1 .AND. isign /= 1) CALL GFT_error(spname//": Error: sign has an invalid value")

  FFT%Isign = isign
  FFT%Scale = scale

  DO i = 0,FFT%Nx/2
    x(2*i)   = REAL(c_in(i),KIND=GFT_prec)
    x(2*i+1) = AIMAG(c_in(i))
  END DO
  nfact = NINT(FFT%Table(0))
  fact(0:nfact-1) = NINT(FFT%Table(0:nfact-1))

  !... Do the job (don't forget that R_out is Y)
  dimx = 2*(FFT%Nx/2)+2 ; debx = 0 ; incx = 1 ; jumpx = 0
  dimy = FFT%Nx         ; deby = 0 ; incy = 1 ; jumpy = 0
  CALL jmcsm1dxy(1,FFT%Nx,fact,100,0,FFT%Table,FFT%TableSize,100+0,    &
                 FFT%Work,FFT%WorkSize,x,dimx,debx,incx,jumpy,r_out,dimy, &
                 deby,incy,jumpy,isign,scale)

  IF( PRESENT(code) ) code = 0
END SUBROUTINE GFT_do_cr_1d

SUBROUTINE GFT_do_cr_2d(FFT, isign, scale, c_in, r_out, code)
  IMPLICIT NONE

  !... Input dummy Parameters
  TYPE(GFT_RCR), INTENT(inout)                         :: FFT
  INTEGER, INTENT(in)                                  :: isign
  REAL(kind=GFT_prec), INTENT(in)                      :: scale
  COMPLEX(kind=GFT_prec), INTENT(in), DIMENSION(0:,0:) :: c_in

  !... Output dummy Parameters
  REAL(kind=GFT_prec), INTENT(out), DIMENSION(0:,0:) :: r_out
  INTEGER, OPTIONAL, INTENT(out)                     :: code

  !... Local variables
  REAL(kind=GFT_prec), DIMENSION(0:2*SIZE(c_in)-1) :: x
  REAL(kind=GFT_prec), DIMENSION(0:SIZE(r_out)-1)  :: y
  INTEGER                     :: i, j, ioff, nfact, mfact
  INTEGER, DIMENSION(0:99)    :: fact
  INTEGER                     :: ideb, ifin, jdeb, jfin
  INTEGER                     :: n_temp, m_temp, nwork_temp
  LOGICAL                     :: debut, fin
  INTEGER                     :: dimy, deby, incy, jumpy
  INTEGER                     :: signe
  REAL(kind=GFT_prec)         :: scale_temp
  CHARACTER(len=*), PARAMETER :: spname = "GFT_do_cr_2d"

  !... Check conditions
  IF (.NOT.FFT%Init) CALL GFT_error(spname//": Error: GFT_init must be called prior to this routine")
  FFT%Ldx1 = SIZE(c_in ,dim=1)
  FFT%Ldy1 = SIZE(r_out,dim=1)
  IF ( FFT%Ldx1 < FFT%Nx/2+1 ) CALL GFT_error(spname//":  Error: Size(C_IN,dim=1) < Nx/2+1")
  IF ( FFT%Ldy1 < FFT%Nx+2 ) CALL GFT_error(spname//":  Error: Size(R_OUT,dim=1) < Nx+2")
  IF ( SIZE(r_out) < (FFT%Nx+2)*FFT%Ny ) CALL GFT_error(spname//":  Error: Size(R_OUT) < (Nx+2)*Ny")
  IF ( SIZE(c_in) < (FFT%Nx/2+1)*FFT%Ny ) CALL GFT_error(spname//":  Error: Size(C_IN) < (Nx/2+1)*Ny")
  IF (isign /=-1 .AND. isign /= 1) CALL GFT_error(spname//": Error: isign has an invalid value")

  FFT%Isign = isign
  FFT%Scale = scale

  DO j = 0, FFT%Ny-1
!dir$ ivdep
!ocl novrec
!cdir nodep
     DO i = 0,FFT%Nx/2
        x(2*i  +2*FFT%Ldx1*j) = REAL(c_in(i,j),KIND=GFT_prec)
        x(2*i+1+2*FFT%Ldx1*j) = AIMAG(c_in(i,j))
     END DO
  END DO

  nfact = NINT(FFT%Table(0))
  mfact = NINT(FFT%Table(nfact)) + nfact
  fact(0:mfact-1) = NINT(FFT%Table(0:mfact-1))

  ! On fait les T.F. sur la premiere dimension en tronconnant sur la deuxieme
  debut = .TRUE.
  DO

    ! Tronconnage
    CALL jmdecoup(FFT%Nx/2+1,4*FFT%Ny,FFT%WorkSize,debut,FFT%even_y,n_temp,ideb,ifin,nwork_temp,fin)

    ! On copie le tableau d'entree dans le tableau de travail sans permuter
    ! les dimensions (n en premier) pour faire d'abord la tf sur m
    ! On en profite pour premultiplier et pour tenir compte du signe
    DO j = 0,FFT%Ny-1
!dir$ ivdep
!ocl novrec
!cdir nodep
      DO i = ideb,ifin
        FFT%Work(             n_temp*j+i-ideb) =       scale*x(2*i  +2*FFT%Ldx1*j)
        FFT%Work(nwork_temp/4+n_temp*j+i-ideb) = isign*scale*x(2*i+1+2*FFT%Ldx1*j)
      END DO
    END DO
    ioff = 0

    ! On fait la FFT complexe -> complexe sur la deuxieme dimension (m)
    CALL jmccm1d(n_temp,FFT%Ny,fact,100,nfact,FFT%Table,FFT%TableSize,100+2*FFT%Nx,FFT%Work,nwork_temp,ioff)

    ! On recopie dans le tableau d'arrivee
    DO j = 0,FFT%Ny-1
!dir$ ivdep
!ocl novrec
!cdir nodep
      DO i = ideb,ifin
         y(2*i  +FFT%Ldy1*j) = FFT%Work(ioff+             n_temp*j+i-ideb)
         y(2*i+1+FFT%Ldy1*j) = FFT%Work(ioff+nwork_temp/4+n_temp*j+i-ideb)
      END DO
    END DO

    ! A-t-on fini ?
    IF (fin) THEN
      EXIT
    ELSE
      debut = .FALSE.
      CYCLE
    END IF

  END DO

  ! On fait les T.F. sur l'autre dimension
  debut = .TRUE.
  DO

    ! Tronconnage
    CALL jmdecoup(FFT%Ny,2*FFT%Nx,FFT%WorkSize,debut,FFT%even_x,m_temp,jdeb,jfin,nwork_temp,fin)

    ! On fait la FFT complexe -> reel sur le premiere dimension (n)
    dimy = FFT%Ldy1*FFT%Ny   ; deby = jdeb*FFT%Ldy1   ; incy = 1 ; jumpy = FFT%Ldy1
    signe = 1
    scale_temp = REAL(1,kind=GFT_prec)
    CALL jmcsm1dxy(m_temp,FFT%Nx,fact,100,0,FFT%Table,FFT%TableSize,100+0, &
                   FFT%Work,nwork_temp, y,dimy,deby,incy,jumpy,y,dimy,deby, &
                   incy,jumpy,signe,scale_temp)

    ! A-t-on fini ?
    IF (fin) THEN
      EXIT
    ELSE
      debut = .FALSE.
      CYCLE
    END IF

  END DO

    DO j = 0,FFT%Ny-1
!dir$ ivdep
!ocl novrec
!cdir nodep
      DO i = 0, FFT%Nx-1
         r_out(i,j) = y(i+FFT%Ldy1*j)
      END DO
    END DO

  IF( PRESENT(code) ) code = 0
END SUBROUTINE GFT_do_cr_2d

SUBROUTINE GFT_do_cr_3d(FFT, isign, scale, c_in, r_out, code)
  IMPLICIT NONE

  !... Input dummy Parameters
  TYPE(GFT_RCR), INTENT(inout)                            :: FFT
  INTEGER, INTENT(in)                                     :: isign
  REAL(kind=GFT_prec), INTENT(in)                         :: scale
  COMPLEX(kind=GFT_prec), INTENT(in), DIMENSION(0:,0:,0:) :: c_in

  !... Output dummy Parameters
  REAL(kind=GFT_prec), INTENT(out), DIMENSION(0:,0:,0:) :: r_out
  INTEGER, OPTIONAL, INTENT(out)                        :: code

  ! Local Variables
  REAL(kind=GFT_prec), DIMENSION(0:2*SIZE(c_in)-1) :: x
  REAL(kind=GFT_prec), DIMENSION(0:SIZE(r_out)-1)  :: y
  INTEGER                  :: i, j, k, ioff
  INTEGER                  :: nfact, mfact, lfact
  INTEGER, DIMENSION(0:99) :: fact
  INTEGER                  :: ideb, ifin, jdeb, jfin, kdeb, kfin, i1, i2, j1, j2
  INTEGER                  :: nltemp, nmtemp, mltemp, nwork_temp, iwork
  LOGICAL                  :: debut, fini
  CHARACTER(len=*), PARAMETER :: spname = "GFT_do_cr_3d"

  ! Check conditions
  IF (.NOT.FFT%Init) CALL GFT_error(spname//": Error: GFT_init must be called prior to this routine")
  FFT%Ldx1 = SIZE(c_in ,dim=1)
  FFT%Ldx2 = SIZE(c_in ,dim=2)
  FFT%Ldy1 = SIZE(r_out,dim=1)
  FFT%Ldy2 = SIZE(r_out,dim=2)
  IF (FFT%Ldx1 < FFT%Nx/2+1) CALL GFT_error(spname//":  Error: Size(C_IN,dim=1) < Nx/2+1")
  IF (FFT%Ldy1 < FFT%Nx+2  ) CALL GFT_error(spname//":  Error: Size(R_OUT,dim=1) < Nx+2")
  IF (FFT%Ldx2 < FFT%Ny    ) CALL GFT_error(spname//":  Error: Size(C_IN,dim=2) < Ny")
  IF (FFT%Ldy2 < FFT%Ny    ) CALL GFT_error(spname//":  Error: Size(R_OUT,dim=2) < Ny")
  IF ( SIZE(r_out) < (FFT%Nx+2)*FFT%Ny*FFT%Nz ) &
                     CALL GFT_error(spname//":  Error: Size(R_OUT) < (Nx+2)*Ny*Nz")
  IF ( SIZE(c_in) < (FFT%Nx/2+1)*FFT%Ny*FFT%Nz ) &
                     CALL GFT_error(spname//":  Error: Size(C_IN) < (Nx/2+1)*Ny*Nz")
  IF (isign /=-1 .AND. isign /= 1) &
                     CALL GFT_error(spname//": Error: isign has an invalid value")

  FFT%Isign = isign
  FFT%Scale = scale

  DO k = 0, FFT%Nz-1
    DO j = 0, FFT%Ny-1
!dir$ ivdep
!ocl novrec
!cdir nodep
       DO i = 0,FFT%Nx/2
         x(2*i  +2*FFT%Ldx1*j+2*FFT%Ldx1*FFT%Ldx2*k) = REAL(c_in(i,j,k),KIND=GFT_prec)
         x(2*i+1+2*FFT%Ldx1*j+2*FFT%Ldx1*FFT%Ldx2*k) = AIMAG(c_in(i,j,k))
       END DO
    END DO
  END DO

  nfact = NINT(FFT%Table(0))
  mfact = NINT(FFT%Table(nfact)) + nfact
  lfact = NINT(FFT%Table(mfact)) + mfact
  fact(0:lfact-1) = NINT(FFT%Table(0:lfact-1))

  ! On fait les T.F. sur la troisieme dimension en tronconnant sur la premiere
  ! et la deuxieme
  debut = .TRUE.
  fini  = .FALSE.
  DO WHILE (.NOT.fini)

    ! Tronconnage
    ! Note : on met npair a .true. car il n'y a pas de restriction dans ce cas
    CALL jmdecoup3(FFT%Nx/2+1,FFT%Ny,4*FFT%Nz,FFT%WorkSize,debut,.TRUE.,ideb,ifin,jdeb,jfin,nmtemp,nwork_temp,fini)
    debut = .FALSE.

    ! On copie le tableau d'entree dans le tableau de travail
    ! On en profite pour premultiplier et pour tenir compte du signe
    ! On prend garde a la gestion des extremites
    DO k = 0,FFT%Nz-1
      iwork = 0
      DO j = jdeb,jfin
        i1 = 0
        i2 = FFT%Nx/2
        IF (j == jdeb) i1 = ideb
        IF (j == jfin) i2 = ifin
!dir$ ivdep
!ocl novrec
!cdir nodep
        DO i = i1,i2
          FFT%Work(             iwork+k*nmtemp) = &
          &       scale*x(2*i  +2*FFT%Ldx1*j+2*FFT%Ldx1*FFT%Ldx2*k)
          FFT%Work(nwork_temp/4+iwork+k*nmtemp) = &
          & isign*scale*x(2*i+1+2*FFT%Ldx1*j+2*FFT%Ldx1*FFT%Ldx2*k)
          iwork = iwork+1
        END DO
      END DO
    END DO

    ! On fait les T.F. sur la troisieme dimension
    ioff = 0
    CALL jmccm1d(nmtemp,FFT%Nz,fact,100,mfact,FFT%Table,FFT%TableSize,100+2*(FFT%Nx+FFT%Ny),FFT%Work,nwork_temp,ioff)

    ! On recopie dans le tableau d'arrivee
    DO k = 0,FFT%Nz-1
      iwork = 0
      DO j = jdeb,jfin
        i1 = 0
        i2 = FFT%Nx/2
        IF (j == jdeb) i1 = ideb
        IF (j == jfin) i2 = ifin
!dir$ ivdep
!ocl novrec
!cdir nodep
        DO i = i1,i2
          y(2*i  +FFT%Ldy1*j+FFT%Ldy1*FFT%Ldy2*k) = FFT%Work(ioff+             iwork+k*nmtemp)
          y(2*i+1+FFT%Ldy1*j+FFT%Ldy1*FFT%Ldy2*k) = FFT%Work(ioff+nwork_temp/4+iwork+k*nmtemp)
          iwork = iwork+1
        END DO
      END DO
    END DO

  END DO

  ! On fait les T.F. sur la deuxieme dimension en tronconnant sur la premiere
  ! et la troisieme
  debut = .TRUE.
  fini  = .FALSE.
  DO WHILE (.NOT.fini)

    ! Tronconnage
    CALL jmdecoup3(FFT%Nx/2+1,FFT%Nz,4*FFT%Ny,FFT%WorkSize,debut,.TRUE.,ideb,ifin,kdeb,kfin,nltemp,nwork_temp,fini)
    debut = .FALSE.

    ! On copie le tableau d'entree dans le tableau de travail
    ! On prend garde a la gestion des extremites
    DO j = 0,FFT%Ny-1
      iwork = 0
      DO k = kdeb,kfin
        i1 = 0
        i2 = FFT%Nx/2
        IF (k == kdeb) i1 = ideb
        IF (k == kfin) i2 = ifin
!dir$ ivdep
!ocl novrec
!cdir nodep
        DO i = i1,i2
          FFT%Work(             iwork+j*nltemp) = &
          & y(2*i  +FFT%Ldy1*j+FFT%Ldy1*FFT%Ldy2*k)
          FFT%Work(nwork_temp/4+iwork+j*nltemp) = &
          & y(2*i+1+FFT%Ldy1*j+FFT%Ldy1*FFT%Ldy2*k)
          iwork = iwork+1
        END DO
      END DO
    END DO

    ! On fait les T.F. sur la deuxieme dimension
    ioff = 0
    CALL jmccm1d(nltemp,FFT%Ny,fact,100,nfact,FFT%Table,FFT%TableSize,100+2*FFT%Nx    ,FFT%Work,nwork_temp,ioff)

    ! On recopie dans le tableau d'arrivee
    DO j = 0,FFT%Ny-1
      iwork = 0
      DO k = kdeb,kfin
        i1 = 0
        i2 = FFT%Nx/2
        IF (k == kdeb) i1 = ideb
        IF (k == kfin) i2 = ifin
!dir$ ivdep
!ocl novrec
!cdir nodep
        DO i = i1,i2
          y(2*i  +FFT%Ldy1*j+FFT%Ldy1*FFT%Ldy2*k) = FFT%Work(ioff             +iwork+j*nltemp)
          y(2*i+1+FFT%Ldy1*j+FFT%Ldy1*FFT%Ldy2*k) = FFT%Work(ioff+nwork_temp/4+iwork+j*nltemp)
          iwork = iwork+1
        END DO
      END DO
    END DO

  END DO

  ! On fait les T.F. sur la premiere dimension en tronconnant sur la deuxieme
  ! et la troisieme
  debut = .TRUE.
  fini  = .FALSE.
  DO WHILE (.NOT.fini)

    ! Tronconnage
    CALL jmdecoup3(FFT%Ny,FFT%Nz,4*(FFT%Nx/2+1),FFT%WorkSize,debut,FFT%even_x,jdeb,jfin,kdeb,kfin,mltemp,nwork_temp,fini)
    debut = .FALSE.

    ! On copie le tableau d'entree dans le tableau de travail
    ! On prend garde a la gestion des extremites
    DO i = 0,FFT%Nx/2
      iwork = 0
      DO k = kdeb,kfin
        j1 = 0
        j2 = FFT%Ny-1
        IF (k == kdeb) j1 = jdeb
        IF (k == kfin) j2 = jfin
!dir$ ivdep
!ocl novrec
!cdir nodep
        DO j = j1,j2
          FFT%Work(             iwork+i*mltemp) = y(2*i  +FFT%Ldy1*j+FFT%Ldy1*FFT%Ldy2*k)
          FFT%Work(nwork_temp/4+iwork+i*mltemp) = y(2*i+1+FFT%Ldy1*j+FFT%Ldy1*FFT%Ldy2*k)
          iwork = iwork+1
        END DO
      END DO
    END DO

    ! On fait les T.F. sur la premiere dimension
    ioff = 0
    CALL jmcsm1d(mltemp,FFT%Nx,fact,100,0    ,FFT%Table,FFT%TableSize,100+0      ,FFT%Work,nwork_temp,ioff)

    ! On recopie dans le tableau d'arrivee
    DO i = 0,FFT%Nx-1
      iwork = 0
      DO k = kdeb,kfin
        j1 = 0
        j2 = FFT%Ny-1
        IF (k == kdeb) j1 = jdeb
        IF (k == kfin) j2 = jfin
!dir$ ivdep
!ocl novrec
!cdir nodep
        DO j = j1,j2
          y(i+FFT%Ldy1*j+FFT%Ldy1*FFT%Ldy2*k) = FFT%Work(ioff+iwork+i*mltemp)
          iwork = iwork+1
        END DO
      END DO
    END DO

  END DO

  DO k = 0, FFT%Nz-1
     DO j = 0, FFT%Ny-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        DO i = 0, FFT%Nx-1
           r_out(i,j,k) = y(i+FFT%Ldy1*j+FFT%Ldy1*FFT%Ldy2*k)
        END DO
     END DO
  END DO
  IF( PRESENT(code) ) code = 0
END SUBROUTINE GFT_do_cr_3d

SUBROUTINE GFT_do_rc_1d(FFT, isign, scale, r_in, c_out, code)
  IMPLICIT NONE

  !... Input dummy Parameters
  TYPE(GFT_RCR), INTENT(inout)                   :: FFT
  INTEGER, INTENT(in)                            :: isign
  REAL(kind=GFT_prec), INTENT(in)                :: scale
  REAL(kind=GFT_prec), INTENT(in), DIMENSION(0:) :: r_in

  !... Output dummy Parameters
  COMPLEX(kind=GFT_prec), INTENT(out), DIMENSION(0:) :: c_out
  INTEGER, OPTIONAL, INTENT(out)                     :: code

  !... Local Variables
  REAL(kind=GFT_prec), DIMENSION(0:2*(SIZE(r_in)/2)+1) :: y
  INTEGER                     :: i, nfact
  INTEGER, DIMENSION(0:99)    :: fact
  INTEGER                     :: dimx, debx, incx, jumpx
  INTEGER                     :: dimy, deby, incy, jumpy
  CHARACTER(len=*), PARAMETER :: spname = "GFT_do_rc_1d"

  !... Check conditions
  IF (.NOT.FFT%Init) CALL GFT_error(spname//": Error: GFT_init must be called prior to this routine")
  FFT%Ldx1 = SIZE(r_in)
  FFT%Ldy1 = SIZE(c_out)
  IF ( FFT%Ldx1 < FFT%Nx) CALL GFT_error(spname//":  Error: Size(R_IN) < Nx")
  IF ( FFT%Ldy1 < FFT%Nx/2+1) CALL GFT_error(spname//":  Error: Size(C_OUT) < Nx/2+1")
  IF (isign /=-1 .AND. isign /= 1) CALL GFT_error(spname//": Error: sign has invalid value")

  FFT%Isign = isign
  FFT%Scale = scale

  nfact = NINT(FFT%Table(0))
  fact(0:nfact-1) = NINT(FFT%Table(0:nfact-1))

  !... Do the job (remember that r_in is x)
  dimx = FFT%Nx         ; debx = 0 ; incx = 1 ; jumpx = 0
  dimy = 2*(FFT%Nx/2)+2 ; deby = 0 ; incy = 1 ; jumpy = 0
  CALL jmscm1dxy(1,FFT%Nx,fact,100,0,FFT%Table,FFT%TableSize,100+0, &
                 FFT%Work,FFT%WorkSize,r_in,dimx,debx,incx,jumpx,y, &
                 dimy,deby,incy,jumpy,isign,scale)
!dir$ ivdep
!ocl novrec
!cdir nodep
  DO i = 0, FFT%Nx/2
     c_out(i) = CMPLX(y(2*i), y(2*i+1),KIND=GFT_prec)
  END DO

  IF( PRESENT(code) ) code = 0
END SUBROUTINE GFT_do_rc_1d

SUBROUTINE GFT_do_rc_2d(FFT, isign, scale, r_in, c_out, code)
  IMPLICIT NONE

  !... Input dummy Parameters
  TYPE(GFT_RCR), INTENT(inout)                      :: FFT
  INTEGER, INTENT(in)                               :: isign
  REAL(kind=GFT_prec), INTENT(in)                   :: scale
  REAL(kind=GFT_prec), INTENT(in), DIMENSION(0:,0:) :: r_in

  !... Output dummy Parameters
  COMPLEX(kind=GFT_prec), INTENT(out), DIMENSION(0:,0:) :: c_out
  INTEGER, OPTIONAL, INTENT(out)                        :: code

  !... Local Variables
  REAL(kind=GFT_prec), DIMENSION(0:SIZE(r_in)-1)    :: x
  REAL(kind=GFT_prec), DIMENSION(0:2*SIZE(c_out)-1) :: y
  INTEGER                     :: i, j, ioff, nfact, mfact
  INTEGER, DIMENSION(0:99)    :: fact
  INTEGER                     :: ideb, ifin, jdeb, jfin, n_temp, m_temp, nwork_temp
  LOGICAL                     :: debut, fin
  INTEGER                     :: dimx, debx, incx, jumpx
  INTEGER                     :: dimy, deby, incy, jumpy
  INTEGER                     :: signe
  CHARACTER(len=*), PARAMETER :: spname = "GFT_do_rc_2d"

  !... Check conditions
  IF (.NOT.FFT%Init) CALL GFT_error(spname//": Error: GFT_init must be called prior to this routine")
  FFT%Ldx1 = SIZE(r_in ,dim=1)
  FFT%Ldy1 = SIZE(c_out,dim=1)
  IF (FFT%Ldx1 < FFT%Nx    ) CALL GFT_error(spname//":  Error: Size(R_IN,dim=1) < Nx")
  IF (FFT%Ldy1 < FFT%Nx/2+1) CALL GFT_error(spname//":  Error: Size(C_OUT,dim=1) < Nx/2+1")
  IF (SIZE(r_in) < FFT%Nx*FFT%Ny) CALL GFT_error(spname//":  Error: Size(R_IN) < Nx*Ny")
  IF (SIZE(c_out) < (FFT%Nx/2+1)*FFT%Ny) CALL GFT_error(spname//":  Error: Size(C_OUT) < (Nx/2+1)*Ny")

  nfact = NINT(FFT%Table(0))
  mfact = NINT(FFT%Table(nfact)) + nfact
  fact(0:mfact-1) = NINT(FFT%Table(0:mfact-1))

  DO j=0, FFT%Ny-1
!dir$ ivdep
!ocl novrec
!cdir nodep
     DO i = 0, FFT%Nx-1
        x(i+FFT%Ldx1*j) = r_in(i,j)
     END DO
  END DO

  ! On fait les T.F. sur la premiere dimension en tronconnant sur la deuxieme
  debut = .TRUE.
  DO

    ! Tronconnage
    CALL jmdecoup(FFT%Ny,2*FFT%Nx,FFT%WorkSize,debut,FFT%even_x,m_temp,jdeb,jfin,nwork_temp,fin)

    ! On fait les T.F. reelles sur la premiere dimension
    ! Tout se passe comme si on faisait une T.F. 1D multiple (d'ordre m)
    dimx = FFT%Ldx1*FFT%Ny   ; debx = jdeb*FFT%Ldx1   ; incx = 1 ; jumpx = FFT%Ldx1
    dimy = 2*FFT%Ldy1*FFT%Ny ; deby = jdeb*2*FFT%Ldy1 ; incy = 1 ; jumpy = 2*FFT%Ldy1
    signe = 1
    CALL jmscm1dxy(m_temp,FFT%Nx,fact,100,0,FFT%Table,FFT%TableSize,100+0, &
      FFT%Work,nwork_temp, x,dimx,debx,incx,jumpx,y,dimy,deby,incy,jumpy,  &
      signe,scale)

    ! A-t-on fini ?
    IF (fin) THEN
      EXIT
    ELSE
      debut = .FALSE.
      CYCLE
    END IF

  END DO

  ! On fait les T.F. sur l'autre dimension
  debut = .TRUE.
  DO

    ! Tronconnage
    CALL jmdecoup(FFT%Nx/2+1,4*FFT%Ny,FFT%WorkSize,debut,FFT%even_y,n_temp,ideb,ifin,nwork_temp,fin)

    ! On copie
    DO j = 0,FFT%Ny-1
!dir$ ivdep
!ocl novrec
!cdir nodep
      DO i = ideb,ifin
        FFT%Work(             n_temp*j+i-ideb) = y(2*i  +2*FFT%Ldy1*j) 
        FFT%Work(nwork_temp/4+n_temp*j+i-ideb) = y(2*i+1+2*FFT%Ldy1*j) 
      END DO
    END DO
    ioff = 0

    ! On fait les T.F. sur l'autre dimension (divisee par deux bien sur)
    ! On va chercher l'autre table des cosinus
    CALL jmccm1d(n_temp,FFT%Ny,fact,100,nfact,FFT%Table,FFT%TableSize,100+2*FFT%Nx,FFT%Work,nwork_temp,ioff)

    ! On recopie dans le tableau d'arrivee
    DO j = 0,FFT%Ny-1
!dir$ ivdep
!ocl novrec
!cdir nodep
      DO i = ideb,ifin
        y(2*i  +2*FFT%Ldy1*j) =       FFT%Work(ioff             +n_temp*j+i-ideb)
        y(2*i+1+2*FFT%Ldy1*j) = isign*FFT%Work(ioff+nwork_temp/4+n_temp*j+i-ideb)
      END DO
    END DO

    ! A-t-on fini ?
    IF (fin) THEN
      EXIT
    ELSE
      debut = .FALSE.
      CYCLE
    END IF

  END DO

  DO j=0, FFT%Ny-1
!dir$ ivdep
!ocl novrec
!cdir nodep
     DO i = 0, FFT%Nx/2
        c_out(i,j) = CMPLX(y(2*i+2*FFT%Ldy1*j), y(2*i+1+2*FFT%Ldy1*j),KIND=GFT_prec)
     END DO
  END DO

  IF( PRESENT(code) ) code = 0
END SUBROUTINE GFT_do_rc_2d

SUBROUTINE GFT_do_rc_3d(FFT, isign, scale, r_in, c_out, code)
  IMPLICIT NONE

  !... Input dummy Parameters
  TYPE(GFT_RCR), INTENT(inout)                         :: FFT
  INTEGER, INTENT(in)                                  :: isign
  REAL(kind=GFT_prec), INTENT(in)                      :: scale
  REAL(kind=GFT_prec), INTENT(in), DIMENSION(0:,0:,0:) :: r_in

  !... Output dummy Parameters
  COMPLEX(kind=GFT_prec), INTENT(out), DIMENSION(0:,0:,0:) :: c_out
  INTEGER, OPTIONAL, INTENT(out)                           :: code

  !... Local Variables
  REAL(kind=GFT_prec), DIMENSION(0:SIZE(r_in)-1)    :: x
  REAL(kind=GFT_prec), DIMENSION(0:2*SIZE(c_out)-1) :: y
  INTEGER                     :: i, j, k, ioff, nfact, mfact, lfact
  INTEGER, DIMENSION(0:99)    :: fact
  INTEGER                     :: ideb, ifin, jdeb, jfin, kdeb, kfin, i1, i2, j1, j2
  INTEGER                     :: nltemp, nmtemp, mltemp, nwork_temp, iwork
  LOGICAL                     :: debut, fini
  CHARACTER(len=*), PARAMETER :: spname = "GFT_do_rc_3d"

  !... Check conditions
  IF (.NOT.FFT%Init) &
         CALL GFT_error(spname//": Error: GFT_init must be called prior to this routine")
  FFT%Ldx1 = SIZE(r_in ,dim=1)
  FFT%Ldx2 = SIZE(r_in ,dim=2)
  FFT%Ldy1 = SIZE(c_out,dim=1)
  FFT%Ldy2 = SIZE(c_out,dim=2)
  IF (FFT%Ldx1 < FFT%Nx    ) CALL GFT_error(spname//":  Error: Size(R_IN,dim=1) < Nx")
  IF (FFT%Ldx2 < FFT%Ny    ) CALL GFT_error(spname//":  Error: Size(R_IN,dim=2) < Ny")
  IF (FFT%Ldy1 < FFT%Nx/2+1) CALL GFT_error(spname//":  Error: Size(C_OUT,dim=1) < Nx/2+1")
  IF (FFT%Ldy2 < FFT%Ny)     CALL GFT_error(spname//":  Error: Size(C_OUT,dim=2) < Ny")
  IF ( SIZE(r_in) < FFT%Nx*FFT%Ny*FFT%Nz ) &
                     CALL GFT_error(spname//":  Error: Size(R_IN) < Nx*Ny*Nz")
  IF ( SIZE(c_out) < (FFT%Nx/2+1)*FFT%Ny*FFT%Nz ) &
                     CALL GFT_error(spname//":  Error: Size(C_OUT) < (Nx/2+1)*Ny*Nz")
  IF (isign /=-1 .AND. isign /= 1) &
                     CALL GFT_error(spname//": Error: isign has an invalid value")

  FFT%Isign = isign
  FFT%Scale = scale

  DO k = 0, FFT%Nz-1
    DO j = 0, FFT%Ny-1
!dir$ ivdep
!ocl novrec
!cdir nodep
       DO i = 0,FFT%Nx-1
         x(i+FFT%Ldx1*j+FFT%Ldx1*FFT%Ldx2*k) = r_in(i,j,k)
       END DO
    END DO
  END DO

  nfact = NINT(FFT%Table(0))
  mfact = NINT(FFT%Table(nfact)) + nfact
  lfact = NINT(FFT%Table(mfact)) + mfact
  fact(0:lfact-1) = NINT(FFT%Table(0:lfact-1))

  ! On fait les T.F. sur la premiere dimension en tronconnant sur la deuxieme
  ! et la troisieme
  debut = .TRUE.
  fini  = .FALSE.
  DO WHILE (.NOT.fini)

    ! Tronconnage
    CALL jmdecoup3(FFT%Ny,FFT%Nz,4*(FFT%Nx/2+1),FFT%WorkSize,debut,FFT%even_x,jdeb,jfin,kdeb,kfin,mltemp,nwork_temp,fini)
    debut = .FALSE.

    ! On copie le tableau d'entree dans le tableau de travail
    ! On prend garde a la gestion des extremites
    DO i = 0,FFT%Nx-1
      iwork = 0
      DO k = kdeb,kfin
        j1 = 0
        j2 = FFT%Ny-1
        IF (k == kdeb) j1 = jdeb
        IF (k == kfin) j2 = jfin
!dir$ ivdep
!ocl novrec
!cdir nodep
        DO j = j1,j2
          FFT%Work(iwork+i*mltemp) = scale*x(i+FFT%Ldx1*j+FFT%Ldx1*FFT%Ldx2*k)
          iwork = iwork+1
        END DO
      END DO
    END DO

    ! On fait les T.F. sur la premiere dimension
    ioff = 0
    CALL jmscm1d(mltemp,FFT%Nx,fact,100,0    ,FFT%Table,FFT%TableSize,100+0      ,FFT%Work,nwork_temp,ioff)

    ! On recopie dans le tableau d'arrivee
    DO i = 0,FFT%Nx/2
      iwork = 0
      DO k = kdeb,kfin
        j1 = 0
        j2 = FFT%Ny-1
        IF (k == kdeb) j1 = jdeb
        IF (k == kfin) j2 = jfin
!dir$ ivdep
!ocl novrec
!cdir nodep
        DO j = j1,j2
           y(2*i  +2*FFT%Ldy1*j+2*FFT%Ldy1*FFT%Ldy2*k) = FFT%Work(ioff+             iwork+i*mltemp)
           y(2*i+1+2*FFT%Ldy1*j+2*FFT%Ldy1*FFT%Ldy2*k) = FFT%Work(ioff+nwork_temp/4+iwork+i*mltemp)
          iwork = iwork+1
        END DO
      END DO
    END DO

  END DO

  ! On fait les T.F. sur la troisieme dimension en tronconnant sur la premiere
  ! et la deuxieme
  debut = .TRUE.
  fini  = .FALSE.
  DO WHILE (.NOT.fini)

    ! Tronconnage
    ! Note : on met npair a .true. car il n'y a pas de restriction dans ce cas
    CALL jmdecoup3(FFT%Nx/2+1,FFT%Ny,4*FFT%Nz,FFT%WorkSize,debut,.TRUE.,ideb,ifin,jdeb,jfin,nmtemp,nwork_temp,fini)
    debut = .FALSE.

    ! On copie le tableau d'entree dans le tableau de travail
    ! On en profite pour premultiplier et pour tenir compte du signe
    ! On prend garde a la gestion des extremites
    DO k = 0,FFT%Nz-1
      iwork = 0
      DO j = jdeb,jfin
        i1 = 0
        i2 = FFT%Nx/2
        IF (j == jdeb) i1 = ideb
        IF (j == jfin) i2 = ifin
!dir$ ivdep
!ocl novrec
!cdir nodep
        DO i = i1,i2
          FFT%work(             iwork+k*nmtemp) = y(2*i  +2*FFT%Ldy1*j+2*FFT%Ldy1*FFT%Ldy2*k)
          FFT%Work(nwork_temp/4+iwork+k*nmtemp) = y(2*i+1+2*FFT%Ldy1*j+2*FFT%Ldy1*FFT%Ldy2*k)
          iwork = iwork+1
        END DO
      END DO
    END DO

    ! On fait les T.F. sur la troisieme dimension
    ioff = 0
    CALL jmccm1d(nmtemp,FFT%Nz,fact,100,mfact,FFT%Table,FFT%TableSize,100+2*(FFT%Nx+FFT%Ny),FFT%Work,nwork_temp,ioff)

    ! On recopie dans le tableau d'arrivee
    DO k = 0,FFT%Nz-1
      iwork = 0
      DO j = jdeb,jfin
        i1 = 0
        i2 = FFT%Nx/2
        IF (j == jdeb) i1 = ideb
        IF (j == jfin) i2 = ifin
!dir$ ivdep
!ocl novrec
!cdir nodep
        DO i = i1,i2
          y(2*i  +2*FFT%Ldy1*j+2*FFT%Ldy1*FFT%Ldy2*k) = FFT%work(ioff+             iwork+k*nmtemp)
          y(2*i+1+2*FFT%Ldy1*j+2*FFT%Ldy1*FFT%Ldy2*k) = FFT%work(ioff+nwork_temp/4+iwork+k*nmtemp)
          iwork = iwork+1
        END DO
      END DO
    END DO

  END DO

  ! On fait les T.F. sur la deuxieme dimension en tronconnant sur la premiere
  ! et la troisieme
  debut = .TRUE.
  fini  = .FALSE.
  DO WHILE (.NOT.fini)

    ! Tronconnage
    CALL jmdecoup3(FFT%Nx/2+1,FFT%Nz,4*FFT%Ny,FFT%WorkSize,debut,.TRUE.,ideb,ifin,kdeb,kfin,nltemp,nwork_temp,fini)
    debut = .FALSE.

    ! On copie le tableau d'entree dans le tableau de travail
    ! On prend garde a la gestion des extremites
    DO j = 0,FFT%Ny-1
      iwork = 0
      DO k = kdeb,kfin
        i1 = 0
        i2 = FFT%Nx/2
        IF (k == kdeb) i1 = ideb
        IF (k == kfin) i2 = ifin
!dir$ ivdep
!ocl novrec
!cdir nodep
        DO i = i1,i2
          FFT%work(             iwork+j*nltemp) = y(2*i  +2*FFT%Ldy1*j+2*FFT%Ldy1*FFT%Ldy2*k)
          FFT%work(nwork_temp/4+iwork+j*nltemp) = y(2*i+1+2*FFT%Ldy1*j+2*FFT%Ldy1*FFT%Ldy2*k)
          iwork = iwork+1
        END DO
      END DO
    END DO

    ! On fait les T.F. sur la deuxieme dimension
    ioff = 0
    CALL jmccm1d(nltemp,FFT%Ny,fact,100,nfact,FFT%Table,FFT%TableSize,100+2*FFT%Nx    ,FFT%Work,nwork_temp,ioff)

    ! On recopie dans le tableau d'arrivee
    DO j = 0,FFT%Ny-1
      iwork = 0
      DO k = kdeb,kfin
        i1 = 0
        i2 = FFT%Nx/2
        IF (k == kdeb) i1 = ideb
        IF (k == kfin) i2 = ifin
!dir$ ivdep
!ocl novrec
!cdir nodep
        DO i = i1,i2
          y(2*i  +2*FFT%Ldy1*j+2*FFT%Ldy1*FFT%Ldy2*k) = &
          &       FFT%Work(ioff             +iwork+j*nltemp)
          y(2*i+1+2*FFT%Ldy1*j+2*FFT%Ldy1*FFT%Ldy2*k) = &
          & isign*FFT%Work(ioff+nwork_temp/4+iwork+j*nltemp)
          iwork = iwork+1
        END DO
      END DO
    END DO

  END DO

  DO k = 0, FFT%Nz-1
     DO j = 0, FFT%Ny-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        DO i = 0, FFT%Nx/2
           c_out(i,j,k) = CMPLX(y(2*i+2*FFT%Ldy1*j+2*FFT%Ldy1*FFT%Ldy2*k), &
                          y(2*i+1+2*FFT%Ldy1*j+2*FFT%Ldy1*FFT%Ldy2*k),KIND=GFT_prec)
        END DO
     END DO
  END DO

  IF( PRESENT(code) ) code = 0
END SUBROUTINE GFT_do_rc_3d

SUBROUTINE GFT_do_mcc_1d(FFT, isign, scale, c_in, c_out, code)
  IMPLICIT NONE

  !... Input dummy Parameters
  TYPE(GFT_MCC), INTENT(inout)                         :: FFT
  INTEGER, INTENT(in)                                  :: isign
  REAL(kind=GFT_prec), INTENT(in)                      :: scale
  COMPLEX(kind=GFT_prec), INTENT(in), DIMENSION(0:,0:) :: c_in

  !... Output dummy Parameters
  COMPLEX(kind=GFT_prec), INTENT(out), DIMENSION(0:,0:) :: c_out
  INTEGER, OPTIONAL, INTENT(out)                        :: code

  !... Local variables
  REAL(kind=GFT_prec), DIMENSION(0:2*SIZE(c_in)-1)  :: x
  REAL(kind=GFT_prec), DIMENSION(0:2*SIZE(c_out)-1) :: y
  INTEGER                     :: i, j, ioff, nfact
  INTEGER, DIMENSION(0:99)    :: fact
  CHARACTER(len=*), PARAMETER :: spname = 'GFT_do_mcc_1d'

  !... Check conditions
  IF (.NOT.FFT%Init) CALL GFT_error(spname//": Error: GFT_init must be called prior to this routine")
  FFT%Ldx1  = SIZE(c_in ,dim=1)
  FFT%Ldy1  = SIZE(c_out,dim=1)
  IF ( SIZE(c_in) < FFT%Nx*FFT%Ny ) CALL GFT_error(spname//":  Error: Size(C_IN) < Nx*Ny")
  IF ( SIZE(c_out) < FFT%Nx*FFT%Ny ) CALL GFT_error(spname//":  Error: Size(C_OUT) < Nx*Ny")
  IF ( FFT%Ldx1 < FFT%Nx ) CALL GFT_error(spname//":  Error: Size(C_IN,dim=1) < Nx")
  IF ( FFT%Ldy1 < FFT%Nx ) CALL GFT_error(spname//":  Error: Size(C_OUT,dim=1) < Nx")
  IF (isign /=-1 .AND. isign /= 1) CALL GFT_error(spname//": Error: isign value is invalid")

  FFT%Isign = isign
  FFT%Scale = scale

  DO j = 0, FFT%Ny-1
!dir$ ivdep
!ocl novrec
!cdir nodep
    DO i = 0,FFT%Nx-1
      x(2*i  +j*2*FFT%Ldx1) = REAL(c_in(i,j),KIND=GFT_prec)
      x(2*i+1+j*2*FFT%Ldx1) = AIMAG(c_in(i,j))
    END DO
  END DO

  nfact = NINT(FFT%table(0))
  fact(0:nfact-1) = NINT(FFT%table(0:nfact-1))

  ! On copie le tableau d'entree dans le tableau de travail
  ! On en profite pour premultiplier et pour tenir compte du signe
  DO i = 0,FFT%Nx-1
!dir$ ivdep
!ocl novrec
!cdir nodep
    DO j = 0,FFT%Ny-1
      FFT%Work(j+FFT%Ny*i)          =       scale*x(2*i  +2*FFT%Ldx1*j)
      FFT%Work(j+FFT%Ny*(FFT%Nx+i)) = isign*scale*x(2*i+1+2*FFT%Ldx1*j)
    END DO
  END DO

  ! On appelle le sous-programme principal
  ioff = 0
  CALL jmccm1d(FFT%Ny, FFT%Nx,fact,100,0,FFT%Table,FFT%TableSize,100+0,FFT%Work,FFT%WorkSize,ioff)

  ! On recopie le tableau de travail dans le tableau de sortie
  DO i = 0,FFT%Nx-1
!dir$ ivdep
!ocl novrec
!cdir nodep
    DO j = 0,FFT%Ny-1
      y(2*i  +2*FFT%Ldy1*j) =       FFT%work(ioff+j+FFT%Ny*i)
      y(2*i+1+2*FFT%Ldy1*j) = isign*FFT%work(ioff+j+FFT%Ny*(FFT%Nx+i))
    END DO
  END DO

  DO j = 0, FFT%Ny-1
     DO i = 0, FFT%Nx-1
        c_out(i,j) = CMPLX(y(2*i+2*FFT%Ldy1*j), y(2*i+1+2*FFT%Ldy1*j),KIND=GFT_prec)
     END DO
  END DO

  IF( PRESENT(code) ) code = 0
END SUBROUTINE GFT_do_mcc_1d

SUBROUTINE GFT_do_mcr_1d(FFT, isign, scale, c_in, r_out, code)
  IMPLICIT NONE

  !... Input dummy Parameters
  TYPE(GFT_MRCR), INTENT(inout)                        :: FFT
  INTEGER, INTENT(in)                                  :: isign
  REAL(kind=GFT_prec), INTENT(in)                      :: scale
  COMPLEX(kind=GFT_prec), INTENT(in), DIMENSION(0:,0:) :: c_in

  !... Output dummy Parameters
  REAL(kind=GFT_prec), INTENT(out), DIMENSION(0:,0:) :: r_out
  INTEGER, OPTIONAL, INTENT(out)                     :: code

  !... Local variables
  REAL(kind=GFT_prec), DIMENSION(0:2*SIZE(c_in)-1) :: x
  REAL(kind=GFT_prec), DIMENSION(0:SIZE(r_out)-1)  :: y
  INTEGER                     :: i, j, nfact
  INTEGER, DIMENSION(0:99)    :: fact
  INTEGER                     :: dimx, debx, incx, jumpx
  INTEGER                     :: dimy, deby, incy, jumpy
  CHARACTER(len=*), PARAMETER :: spname = "GFT_do_mcr_1d"

  !... Check conditions
  IF (.NOT.FFT%Init) CALL GFT_error(spname//": Error: GFT_init must be called prior to this routine")
  FFT%Ldx1 = SIZE(c_in ,dim=1)
  FFT%Ldy1 = SIZE(r_out,dim=1)
  IF ( FFT%Ldx1 < FFT%Nx/2+1 ) CALL GFT_error(spname//":  Error: Size(C_IN,dim=1) < Nx/2+1")
  IF ( FFT%Ldy1 < FFT%Nx     ) CALL GFT_error(spname//":  Error: Size(R_OUT,dim=1) < Nx")
  IF ( SIZE(r_out) < FFT%Nx*FFT%Ny ) CALL GFT_error(spname//":  Error: Size(R_OUT) < Nx*Ny")
  IF ( SIZE(c_in) < (FFT%Nx/2+1)*FFT%Ny ) CALL GFT_error(spname//":  Error: Size(C_IN) < (Nx/2+1)*Ny")
  IF (isign /=-1 .AND. isign /= 1) CALL GFT_error(spname//": Error: isign has an invalid value")

  FFT%Isign = isign
  FFT%Scale = scale

  DO j = 0, FFT%Ny-1
!dir$ ivdep
!ocl novrec
!cdir nodep
     DO i = 0,FFT%Nx/2
        x(2*i  +2*FFT%Ldx1*j) = REAL(c_in(i,j),KIND=GFT_prec)
        x(2*i+1+2*FFT%Ldx1*j) = AIMAG(c_in(i,j))
     END DO
  END DO

  nfact = NINT(FFT%Table(0))
  fact(0:nfact-1) = NINT(FFT%Table(0:nfact-1))

  ! On appelle le sous-programme principal
  dimx = 2*FFT%Ldx1*FFT%Ny ; debx = 0 ; incx = 1 ; jumpx = 2*FFT%Ldx1
  dimy = FFT%Ldy1*FFT%Ny   ; deby = 0 ; incy = 1 ; jumpy = FFT%Ldy1
  CALL jmcsm1dxy(FFT%Ny,FFT%Nx,fact,100,0,FFT%Table,FFT%TableSize,100+0,  &
                 FFT%Work,FFT%WorkSize,x,dimx,debx,incx,jumpx,y,dimy,deby,&
                 incy,jumpy,isign,scale)

  DO j = 0, FFT%Ny-1
     DO i = 0, FFT%Nx-1
        r_out(i,j) = y(i + FFT%Ldy1*j)
     END DO
  END DO

  IF( PRESENT(code) ) code = 0
END SUBROUTINE GFT_do_mcr_1d

SUBROUTINE GFT_do_mrc_1d(FFT, isign, scale, r_in, c_out, code)
  IMPLICIT NONE

  !... Input dummy Parameters
  TYPE(GFT_MRCR), INTENT(inout)                     :: FFT
  INTEGER, INTENT(in)                               :: isign
  REAL(kind=GFT_prec), INTENT(in)                   :: scale
  REAL(kind=GFT_prec), INTENT(in), DIMENSION(0:,0:) :: r_in

  !... Output dummy Parameters
  COMPLEX(kind=GFT_prec), INTENT(out), DIMENSION(0:,0:) :: c_out
  INTEGER, OPTIONAL, INTENT(out)                        :: code

  !... Local Variables
  REAL(kind=GFT_prec), DIMENSION(0:SIZE(r_in)-1)    :: x
  REAL(kind=GFT_prec), DIMENSION(0:2*SIZE(c_out)-1) :: y
  INTEGER                     :: i, j, nfact
  INTEGER, DIMENSION(0:99)    :: fact
  INTEGER                     :: dimx, debx, incx, jumpx
  INTEGER                     :: dimy, deby, incy, jumpy
  CHARACTER(len=*), PARAMETER :: spname = "GFT_do_mrc_1d"

  !... Check conditions
  IF (.NOT.FFT%Init) CALL GFT_error(spname//": Error: GFT_init must be called prior to this routine")
  FFT%Ldx1 = SIZE(r_in ,dim=1)
  FFT%Ldy1 = SIZE(c_out,dim=1)
  IF (FFT%Ldx1 < FFT%Nx    ) CALL GFT_error(spname//":  Error: Size(R_IN,dim=1) < Nx")
  IF (FFT%Ldy1 < FFT%Nx/2+1) CALL GFT_error(spname//":  Error: Size(C_OUT,dim=1) < Nx/2+1")
  IF (SIZE(r_in) < FFT%Nx*FFT%Ny) CALL GFT_error(spname//":  Error: Size(R_IN) < Nx*Ny")
  IF (SIZE(c_out) < (FFT%Nx/2+1)*FFT%Ny) CALL GFT_error(spname//":  Error: Size(C_OUT) < (Nx/2+1)*Ny")

  DO j=0, FFT%Ny-1
!dir$ ivdep
!ocl novrec
!cdir nodep
     DO i = 0, FFT%Nx-1
        x(i+FFT%Ldx1*j) = r_in(i,j)
     END DO
  END DO

  nfact = NINT(FFT%Table(0))
  fact(0:nfact-1) = NINT(FFT%Table(0:nfact-1))

  ! On appelle le sous-programme principal
  dimx = FFT%Ldx1*FFT%Ny   ; debx = 0 ; incx = 1 ; jumpx = FFT%Ldx1
  dimy = 2*FFT%Ldy1*FFT%Ny ; deby = 0 ; incy = 1 ; jumpy = 2*FFT%Ldy1
  CALL jmscm1dxy(FFT%Ny,FFT%Nx,fact,100,0,FFT%Table,FFT%TableSize,100+0, &
                 FFT%Work,FFT%WorkSize,x,dimx,debx,incx,jumpx,y,dimy,deby,&
                 incy,jumpy,isign,scale)

  DO j = 0, FFT%Ny-1
     DO i = 0, FFT%Nx/2
        c_out(i,j) = CMPLX(y(2*i+2*FFT%Ldy1*j), y(2*i+1+2*FFT%Ldy1*j),KIND=GFT_prec)
     END DO
  END DO

  IF( PRESENT(code) ) code = 0
END SUBROUTINE GFT_do_mrc_1d

END MODULE GFT_do
