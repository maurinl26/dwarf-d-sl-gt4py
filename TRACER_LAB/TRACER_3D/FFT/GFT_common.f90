!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! -*- Mode: F90 -*- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! GFT_common.f90 --- GFT common objects
!! 
!! Auteur          : Jalel Chergui (CNRS/IDRIS) <Jalel.Chergui@idris.fr>
!! Créé le         : Tue Feb 19 10:02:32 2002
!! Dern. mod. par  : Jalel Chergui (CNRS/IDRIS) <Jalel.Chergui@idris.fr>
!! Dern. mod. le   : Thu Jul  4 13:13:46 2002
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Permission is granted to copy and distribute this file or modified
! versions of this file for no fee, provided the copyright notice and
! this permission notice are preserved on all copies.
! Copyright (C) Februry 2002, CNRS/IDRIS, Jalel.Chergui@idris.fr.
!
MODULE GFT_common
  IMPLICIT NONE

  !... Floating point precisions
  INTEGER, PARAMETER :: GFT_R4  = SELECTED_REAL_KIND(6)
  INTEGER, PARAMETER :: GFT_R8  = SELECTED_REAL_KIND(12)
  INTEGER, PARAMETER :: GFT_R16 = SELECTED_REAL_KIND(24)

  !... Recommended Default Floating Point Precision
  INTEGER, PARAMETER :: GFT_prec=GFT_R8

  !... Default copyright
  LOGICAL            :: copyright=.FALSE.

  !... Derived type for 1D, 2D and 3D Complex-Complex FFT
  TYPE GFT_CC
    INTEGER             :: Nx
    INTEGER             :: Ny
    INTEGER             :: Nz
    INTEGER             :: BF=1
    INTEGER             :: Isign=0
    REAL(KIND=GFT_prec) :: Scale
    INTEGER             :: Ldx1
    INTEGER             :: Ldx2
    INTEGER             :: Ldy1
    INTEGER             :: Ldy2
    INTEGER             :: TableSize
    INTEGER             :: WorkSize
    LOGICAL             :: Init=.FALSE.
    REAL(KIND=GFT_prec), POINTER, DIMENSION(:) :: Work
    REAL(KIND=GFT_prec), POINTER, DIMENSION(:) :: Table
  END TYPE GFT_CC

  !... Derived type for 1D, 2D and 3D Real-Complex and Complex-Real FFT
  TYPE GFT_RCR
    INTEGER             :: Nx
    INTEGER             :: Ny
    INTEGER             :: Nz
    INTEGER             :: BF=1
    INTEGER             :: Isign=0
    REAL(KIND=GFT_prec) :: Scale
    INTEGER             :: Ldx1
    INTEGER             :: Ldx2
    INTEGER             :: Ldy1
    INTEGER             :: Ldy2
    LOGICAL             :: even_x
    LOGICAL             :: even_y
    LOGICAL             :: even_z
    LOGICAL             :: Init=.FALSE.
    INTEGER             :: TableSize
    INTEGER             :: WorkSize
    REAL(KIND=GFT_prec), POINTER, DIMENSION(:) :: Work
    REAL(KIND=GFT_prec), POINTER, DIMENSION(:) :: Table
  END TYPE GFT_RCR


  !... Derived type for Multiple 1D Complex-Complex FFT
  TYPE GFT_MCC
    INTEGER             :: Nx
    INTEGER             :: Ny
    INTEGER             :: BF=1
    INTEGER             :: Isign=0
    REAL(KIND=GFT_prec) :: Scale
    INTEGER             :: Ldx1
    INTEGER             :: Ldy1
    LOGICAL             :: Init=.FALSE.
    INTEGER             :: TableSize
    INTEGER             :: WorkSize
    REAL(KIND=GFT_prec), POINTER, DIMENSION(:) :: Work
    REAL(KIND=GFT_prec), POINTER, DIMENSION(:) :: Table
  END TYPE GFT_MCC

  !... Derived type for Multiple 1D Real-Complex and Complex-Real FFT
  TYPE GFT_MRCR
    INTEGER             :: Nx
    INTEGER             :: Ny
    INTEGER             :: BF=1
    INTEGER             :: Isign=0
    REAL(KIND=GFT_prec) :: Scale
    INTEGER             :: Ldx1
    INTEGER             :: Ldy1
    LOGICAL             :: even_x
    LOGICAL             :: even_y
    LOGICAL             :: Init=.FALSE.
    INTEGER             :: TableSize
    INTEGER             :: WorkSize
    REAL(KIND=GFT_prec), POINTER, DIMENSION(:) :: Work
    REAL(KIND=GFT_prec), POINTER, DIMENSION(:) :: Table
  END TYPE GFT_MRCR

  CONTAINS

  SUBROUTINE GFT_Error( Message )
    IMPLICIT NONE

    !... Input dummy arguments
    CHARACTER(LEN=*) :: Message

    PRINT *,Message
    STOP
  END SUBROUTINE GFT_Error

  SUBROUTINE GFT_Warning( Message )
    IMPLICIT NONE

    !... Input dummy arguments
    CHARACTER(LEN=*) :: Message

    PRINT *,Message
  END SUBROUTINE GFT_Warning

END MODULE GFT_common
