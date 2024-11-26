!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! -*- Mode: F90 -*- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! GFT_set.f90 --- GFT initialization routines
!!
!! Auteur          : Jalel Chergui (CNRS/IDRIS) <Jalel.Chergui@idris.fr>
!! Créé le         : Tue Feb 19 10:26:52 2002
!! Dern. mod. par  : Jalel Chergui (CNRS/IDRIS) <Jalel.Chergui@idris.fr>
!! Dern. mod. le   : Wed May 15 14:19:14 2002
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Permission is granted to copy and distribute this file or modified
! versions of this file for no fee, provided the copyright notice and
! this permission notice are preserved on all copies.
! Copyright February 2002, CNRS/IDRIS, Jalel.Chergui@idris.fr
!
MODULE GFT_set
  USE GFT_common
  USE GFT_JMFFT

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: GFT_set_bf
  INTERFACE GFT_set_bf
     MODULE PROCEDURE GFT_set_bf_cc, GFT_set_bf_rcr, GFT_set_bf_mcc, &
                      GFT_set_bf_mrcr
  END INTERFACE GFT_set_bf

  PUBLIC :: GFT_set_fft
  INTERFACE GFT_set_fft
     MODULE PROCEDURE GFT_set_cc_1d,  GFT_set_cc_2d,  GFT_set_cc_3d,  &
                      GFT_set_rcr_1d, GFT_set_rcr_2d, GFT_set_rcr_3d, &
                      GFT_set_mcc,   GFT_set_mrcr
  END INTERFACE GFT_set_fft

  CONTAINS

  SUBROUTINE GFT_set_bf_cc(BF, FFT, code)
    IMPLICIT NONE

    !... Input dummy arguments
    INTEGER, INTENT(IN) :: BF

    !... Output dummy arguments
    TYPE(GFT_CC), INTENT(OUT)      :: FFT
    INTEGER, OPTIONAL, INTENT(OUT) :: code

    FFT%BF=BF
    IF( PRESENT(code) ) code=0
  END SUBROUTINE GFT_set_bf_cc

  SUBROUTINE GFT_set_bf_rcr(BF, FFT, code)
    IMPLICIT NONE

    !... Input dummy arguments
    INTEGER, INTENT(IN) :: BF

    !... Output dummy arguments
    TYPE(GFT_RCR), INTENT(OUT)     :: FFT
    INTEGER, OPTIONAL, INTENT(OUT) :: code

    FFT%BF=BF
    IF( PRESENT(code) ) code=0
  END SUBROUTINE GFT_set_bf_rcr

  SUBROUTINE GFT_set_bf_mcc(BF, FFT, code)
    IMPLICIT NONE

    !... Input dummy arguments
    INTEGER, INTENT(IN) :: BF

    !... Output dummy arguments
    TYPE(GFT_MCC), INTENT(OUT)     :: FFT
    INTEGER, OPTIONAL, INTENT(OUT) :: code

    FFT%BF=BF
    IF( PRESENT(code) ) code=0
  END SUBROUTINE GFT_set_bf_mcc

  SUBROUTINE GFT_set_bf_mrcr(BF, FFT, code)
    IMPLICIT NONE

    !... Input dummy arguments
    INTEGER, INTENT(IN) :: BF

    !... Output dummy arguments
    TYPE(GFT_MRCR), INTENT(OUT)    :: FFT
    INTEGER, OPTIONAL, INTENT(OUT) :: code

    FFT%BF=BF
    IF( PRESENT(code) ) code=0
  END SUBROUTINE GFT_set_bf_mrcr

  SUBROUTINE GFT_set_cc_1d(Nx, FFT, code)
    IMPLICIT NONE

    !... Input dummy arguments
    INTEGER, INTENT(IN)           :: Nx

    !... Output dummy arguments
    TYPE(GFT_CC), INTENT(OUT)      :: FFT
    INTEGER, OPTIONAL, INTENT(OUT) :: code

    !... Local variables
    INTEGER, DIMENSION(0:99)    :: fact
    INTEGER                     :: nfact
    CHARACTER(LEN=*), PARAMETER :: spname="GFT_set_cc_1d"

    FFT%Nx = Nx
    IF( FFT%Nx < 1 ) CALL GFT_Error(spname//": Error: Nx < 1")

    FFT%Init      = .TRUE.
    FFT%TableSize = 100+2*FFT%Nx
    ALLOCATE(FFT%Table(0:FFT%TableSize-1))

    FFT%WorkSize = 4*FFT%BF*FFT%Nx
    ALLOCATE(FFT%Work(0:FFT%WorkSize-1))

    ! Pour la factorisation
    CALL jmfact(FFT%Nx,fact,100,    0,nfact)
    FFT%Table(0:nfact-1) = fact(0:nfact-1)

    ! Pour les sinus et cosinus
    CALL jmtable(FFT%Table,FFT%TableSize,100+0,FFT%Nx)

    IF( PRESENT(code) ) code=0
  END SUBROUTINE GFT_set_cc_1d

  SUBROUTINE GFT_set_cc_2d(Nx, Ny, FFT, code)
    IMPLICIT NONE

    !... Input dummy arguments
    INTEGER, INTENT(IN)           :: Nx, Ny

    !... Output dummy arguments
    TYPE(GFT_CC), INTENT(OUT)      :: FFT
    INTEGER, OPTIONAL, INTENT(OUT) :: code

    !... Internal variables
    INTEGER, DIMENSION(0:99)    :: fact
    INTEGER                     :: nfact, mfact
    CHARACTER(LEN=*), PARAMETER :: spname="GFT_set_cc_2d"

    FFT%Nx = Nx
    FFT%Ny = Ny
    IF( FFT%Nx < 1 ) CALL GFT_Error(spname//": Error: NX < 1")
    IF( FFT%Ny < 1 ) CALL GFT_Error(spname//": Error: NY < 1")

    FFT%Init      = .TRUE.
    FFT%TableSize = 100+2*(FFT%Nx + FFT%Ny)
    ALLOCATE(FFT%Table(0:FFT%TableSize-1))

    FFT%WorkSize = 4*FFT%BF*MAX(FFT%Nx, FFT%Ny)
    CALL JMsetnwork(FFT%WorkSize)
    ALLOCATE(FFT%Work(0:FFT%WorkSize-1))

    !... Pour la factorisation
    CALL jmfact(FFT%Nx,fact,100,    0,nfact)
    CALL jmfact(FFT%Ny,fact,100,nfact,mfact)
    FFT%Table(0:mfact-1) = fact(0:mfact-1)

    !... Pour les sinus et cosinus
    CALL jmtable(FFT%Table,FFT%TableSize,100+0       ,FFT%Nx)
    CALL jmtable(FFT%Table,FFT%TableSize,100+2*FFT%Nx,FFT%Ny)

    IF( PRESENT(code) ) code=0
  END SUBROUTINE GFT_set_cc_2d

  SUBROUTINE GFT_set_cc_3d(Nx, Ny, Nz, FFT, code)
    IMPLICIT NONE

    !... Input dummy arguments
    INTEGER, INTENT(IN)           :: Nx, Ny, Nz

    !... Output dummy arguments
    TYPE(GFT_CC), INTENT(OUT)      :: FFT
    INTEGER, OPTIONAL, INTENT(OUT) :: code

    !... Internal variables
    INTEGER, DIMENSION(0:99)    :: fact
    INTEGER                     :: nfact, mfact, lfact
    CHARACTER(LEN=*), PARAMETER :: spname="GFT_set_cc_3d"

    FFT%Nx = Nx
    FFT%Ny = Ny
    FFT%Nz = Nz
    IF( FFT%Nx < 1 ) CALL GFT_Error(spname//": Error: Nx < 1")
    IF( FFT%Ny < 1 ) CALL GFT_Error(spname//": Error: Ny < 1")
    IF( FFT%Nz < 1 ) CALL GFT_Error(spname//": Error: Nz < 1")

    FFT%Init      = .TRUE.
    FFT%TableSize = 100+2*(FFT%Nx + FFT%Ny + FFT%Nz)
    ALLOCATE(FFT%Table(0:FFT%TableSize-1))

    FFT%WorkSize = 4*FFT%BF*MAX(FFT%Nx, FFT%Ny, FFT%Nz)
    CALL JMsetnwork(FFT%WorkSize)
    ALLOCATE(FFT%Work(0:FFT%WorkSize-1))

    !... Pour la factorisation
    CALL jmfact(FFT%Nx,fact,100,    0,nfact)
    CALL jmfact(FFT%Ny,fact,100,nfact,mfact)
    CALL jmfact(FFT%Nz,fact,100,mfact,lfact)
    FFT%Table(0:lfact-1) = fact(0:lfact-1)

    !... Pour les sinus et cosinus
    CALL jmtable(FFT%Table,FFT%TableSize,100+0       ,         FFT%Nx)
    CALL jmtable(FFT%Table,FFT%TableSize,100+2*FFT%Nx,         FFT%Ny)
    CALL jmtable(FFT%Table,FFT%TableSize,100+2*(FFT%Nx+FFT%Ny),FFT%Nz)

    IF( PRESENT(code) ) code=0
  END SUBROUTINE GFT_set_cc_3d

  SUBROUTINE GFT_set_rcr_1d(Nx, FFT, code)
    IMPLICIT NONE

    !... Input dummy arguments
    INTEGER, INTENT(IN)           :: Nx

    !... Output dummy arguments
    TYPE(GFT_RCR), INTENT(OUT)     :: FFT
    INTEGER, OPTIONAL, INTENT(OUT) :: code

    !... Local variables
    INTEGER, DIMENSION(0:99)    :: fact
    INTEGER                     :: nfact
    CHARACTER(LEN=*), PARAMETER :: spname="GFT_set_rcr_1d"

    FFT%Nx = Nx
    IF( FFT%Nx < 1 ) CALL GFT_Error(spname//": Error: Nx < 1")
    FFT%even_x = (MOD(FFT%Nx,2) == 0 )
    IF( .NOT.FFT%even_x ) CALL GFT_Error(spname//": Error: Nx not even")

    FFT%Init      = .TRUE.
    FFT%TableSize = 100+2*FFT%Nx
    ALLOCATE(FFT%Table(0:FFT%TableSize-1))

    FFT%WorkSize = 2*FFT%BF*FFT%Nx
    ALLOCATE(FFT%Work(0:FFT%WorkSize-1))

    ! Pour la factorisation
    CALL jmfact(FFT%Nx,fact,100,    0,nfact)
    FFT%Table(0:nfact-1) = fact(0:nfact-1)

    ! Pour les sinus et cosinus
    CALL jmtable(FFT%Table,FFT%TableSize,100+0,FFT%Nx)

    IF( PRESENT(code) ) code=0    
  END SUBROUTINE GFT_set_rcr_1d

  SUBROUTINE GFT_set_rcr_2d(Nx, Ny, FFT, code)
    IMPLICIT NONE

    !... Input dummy arguments
    INTEGER, INTENT(IN)           :: Nx, Ny

    !... Output dummy arguments
    TYPE(GFT_RCR), INTENT(OUT)     :: FFT
    INTEGER, OPTIONAL, INTENT(OUT) :: code

    !... Internal variables
    INTEGER, DIMENSION(0:99)    :: fact
    INTEGER                     :: nfact, mfact
    CHARACTER(LEN=*), PARAMETER :: spname="GFT_set_rcr_2d"

    FFT%Nx = Nx
    FFT%Ny = Ny
    IF( FFT%Nx < 1 ) CALL GFT_Error(spname//": Error: Nx < 1")
    IF( FFT%Ny < 1 ) CALL GFT_Error(spname//": Error: Ny < 1")
    FFT%even_x = (MOD(FFT%Nx,2) == 0); FFT%even_y = (MOD(FFT%Ny,2) == 0)
    IF( .NOT.FFT%even_x .AND. .NOT.FFT%even_y) CALL GFT_Error(spname//": Error: Nx or Ny not even")

    FFT%Init      = .TRUE.
    FFT%TableSize = 100+2*(FFT%Nx + FFT%Ny)
    ALLOCATE(FFT%Table(0:FFT%TableSize-1))

    FFT%WorkSize = 4*FFT%BF*MAX(FFT%Nx, FFT%Ny)
    CALL JMsetnwork(FFT%WorkSize)
    ALLOCATE(FFT%Work(0:FFT%WorkSize-1))

    !... Pour la factorisation
    CALL jmfact(FFT%Nx,fact,100,    0,nfact)
    CALL jmfact(FFT%Ny,fact,100,nfact,mfact)
    FFT%Table(0:mfact-1) = fact(0:mfact-1)

    !... Pour les sinus et cosinus
    CALL jmtable(FFT%Table,FFT%TableSize,100+0       ,FFT%Nx)
    CALL jmtable(FFT%Table,FFT%TableSize,100+2*FFT%Nx,FFT%Ny)

    IF( PRESENT(code) ) code=0
  END SUBROUTINE GFT_set_rcr_2d

  SUBROUTINE GFT_set_rcr_3d(Nx, Ny, Nz, FFT, code)
    IMPLICIT NONE

    !... Input dummy arguments
    INTEGER, INTENT(IN)           :: Nx, Ny, Nz

    !... Output dummy arguments
    TYPE(GFT_RCR), INTENT(OUT)     :: FFT
    INTEGER, OPTIONAL, INTENT(OUT) :: code

    !... Internal variables
    INTEGER, DIMENSION(0:99)    :: fact
    INTEGER                     :: nfact, mfact, lfact
    CHARACTER(LEN=*), PARAMETER :: spname="GFT_set_rcr_3d"

    FFT%Nx = Nx
    FFT%Ny = Ny
    FFT%Nz = Nz
    IF( FFT%Nx < 1 ) CALL GFT_Error(spname//": Error: Nx < 1")
    IF( FFT%Ny < 1 ) CALL GFT_Error(spname//": Error: Ny < 1")
    IF( FFT%Nz < 1 ) CALL GFT_Error(spname//": Error: Nz < 1")
    FFT%even_x = (MOD(FFT%Nx,2) == 0)
    FFT%even_y = (MOD(FFT%Ny,2) == 0)
    FFT%even_z = (MOD(FFT%Nz,2) == 0)
    IF( .NOT.FFT%even_x .AND. .NOT.FFT%even_y .AND. .NOT.FFT%even_z) &
                  CALL GFT_Error(spname//": Error: Nx or Ny or Nz not even")

    FFT%Init      = .TRUE.
    FFT%TableSize = 100+2*(FFT%Nx + FFT%Ny + FFT%Nz)
    ALLOCATE(FFT%Table(0:FFT%TableSize-1))

    FFT%WorkSize = 4*FFT%BF*MAX(FFT%Nx, FFT%Ny, FFT%Nz)
    CALL JMsetnwork(FFT%WorkSize)
    ALLOCATE(FFT%Work(0:FFT%WorkSize-1))

    !... Pour la factorisation
    CALL jmfact(FFT%Nx,fact,100,    0,nfact)
    CALL jmfact(FFT%Ny,fact,100,nfact,mfact)
    CALL jmfact(FFT%Nz,fact,100,mfact,lfact)
    FFT%Table(0:lfact-1) = fact(0:lfact-1)

    !... Pour les sinus et cosinus
    CALL jmtable(FFT%Table,FFT%TableSize,100+0                ,FFT%Nx)
    CALL jmtable(FFT%Table,FFT%TableSize,100+2*FFT%Nx         ,FFT%Ny)
    CALL jmtable(FFT%Table,FFT%TableSize,100+2*(FFT%Nx+FFT%Ny),FFT%Nz)

    IF( PRESENT(code) ) code=0
  END SUBROUTINE GFT_set_rcr_3d

  SUBROUTINE GFT_set_mcc(Nx, Ny, FFT, code)
    IMPLICIT NONE

    !... Input dummy arguments
    INTEGER, INTENT(IN)           :: Nx, Ny

    !... Output dummy arguments
    TYPE(GFT_MCC), INTENT(OUT)     :: FFT
    INTEGER, OPTIONAL, INTENT(OUT) :: code

    !... Internal variables
    INTEGER, DIMENSION(0:99)    :: fact
    INTEGER                     :: nfact
    CHARACTER(LEN=*), PARAMETER :: spname="GFT_set_mcc"

    FFT%Nx = Nx
    FFT%Ny = Ny
    IF( FFT%Nx < 1 ) CALL GFT_Error(spname//": Error: Nx < 1")
    IF( FFT%Ny < 1 ) CALL GFT_Error(spname//": Error: Ny < 1")

    FFT%Init      = .TRUE.
    FFT%TableSize = 100+2*FFT%Nx
    ALLOCATE(FFT%Table(0:FFT%TableSize-1))

    FFT%WorkSize = 4*FFT%BF*FFT%Nx*FFT%Ny
    ALLOCATE(FFT%Work(0:FFT%WorkSize-1))

    !... Pour la factorisation
    CALL jmfact(FFT%Nx,fact,100,0,nfact)
    FFT%Table(0:nfact-1) = fact(0:nfact-1)

    !... Pour les sinus et cosinus
    CALL jmtable(FFT%Table,FFT%TableSize,100+0,FFT%Nx)

    IF( PRESENT(code) ) code=0
  END SUBROUTINE GFT_set_mcc

  SUBROUTINE GFT_set_mrcr(Nx, Ny, FFT, code)
    IMPLICIT NONE

    !... Input dummy arguments
    INTEGER, INTENT(IN)           :: Nx, Ny

    !... Output dummy arguments
    TYPE(GFT_MRCR), INTENT(OUT)    :: FFT
    INTEGER, OPTIONAL, INTENT(OUT) :: code

    !... Internal variables
    INTEGER, DIMENSION(0:99)    :: fact
    INTEGER                     :: nfact
    CHARACTER(LEN=*), PARAMETER :: spname="GFT_set_mrcr"

    FFT%Nx = Nx
    FFT%Ny = Ny
    !... Check conditions
    IF( FFT%Nx < 1 ) CALL GFT_Error(spname//": Error: Nx < 1")
    IF( FFT%Ny < 1 ) CALL GFT_Error(spname//": Error: Ny < 1")
    FFT%even_x = (MOD(FFT%Nx,2) == 0) ; FFT%even_y = (MOD(FFT%Ny,2) == 0)
    IF( .NOT.FFT%even_x .AND. .NOT.FFT%even_y ) CALL GFT_Error(spname//": Error: Nx or Ny not even")

    FFT%Init      = .TRUE.
    FFT%TableSize = 100+2*FFT%Nx
    ALLOCATE(FFT%Table(0:FFT%TableSize-1))

    FFT%WorkSize = 2*FFT%BF*FFT%Nx*FFT%Ny
    ALLOCATE(FFT%Work(0:FFT%WorkSize-1))

    ! Pour la factorisation
    CALL jmfact(FFT%Nx,fact,100,    0,nfact)
    FFT%Table(0:nfact-1) = fact(0:nfact-1)

    ! Pour les sinus et cosinus
    CALL jmtable(FFT%Table,FFT%TableSize,100+0,FFT%Nx)

    IF( PRESENT(code) ) code=0
  END SUBROUTINE GFT_set_mrcr

END MODULE GFT_set
