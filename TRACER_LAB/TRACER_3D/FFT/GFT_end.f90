!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! -*- Mode: F90 -*- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! GFT_end.f90 --- Free/Reset GFT common objects
!! 
!! Auteur          : Jalel Chergui (CNRS/IDRIS) <Jalel.Chergui@idris.fr>
!! Créé le         : Tue Feb 19 10:03:47 2002
!! Dern. mod. par  : Jalel Chergui (CNRS/IDRIS) <Jalel.Chergui@idris.fr>
!! Dern. mod. le   : Mon Apr 15 14:22:48 2002
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Permission is granted to copy and distribute this file or modified
! versions of this file for no fee, provided the copyright notice and
! this permission notice are preserved on all copies.
! Copyright (C) Februry 2002, CNRS/IDRIS, Jalel.Chergui@idris.fr.
!
MODULE GFT_end
  USE GFT_common

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: GFT_end_fft
  INTERFACE GFT_end_fft
     MODULE PROCEDURE GFT_end_cc, GFT_end_rcr, GFT_end_mcc, GFT_end_mrcr
  END INTERFACE GFT_end_fft

  CONTAINS

  SUBROUTINE GFT_end_cc(FFT, code)
    IMPLICIT NONE

    !... Input dummy arguments
    TYPE(GFT_CC), INTENT(INOUT) :: FFT

    !... Output dummy arguments
    INTEGER, OPTIONAL, INTENT(OUT) :: code

    FFT%Nx    = -1
    FFT%Ny    = -1
    FFT%Nz    = -1
    FFT%BF    =  1
    FFT%Scale =  1.0_GFT_prec
    FFT%Isign =  0
    FFT%Init  = .FALSE.

    IF ( ASSOCIATED(FFT%Table) ) DEALLOCATE(FFT%Table)
    IF ( ASSOCIATED(FFT%Work) )  DEALLOCATE(FFT%Work)

    IF( PRESENT(code) ) code=0
  END SUBROUTINE GFT_end_cc

  SUBROUTINE GFT_end_rcr(FFT, code)
    IMPLICIT NONE

    !... Input dummy arguments
    TYPE(GFT_RCR), INTENT(INOUT) :: FFT

    !... Output dummy arguments
    INTEGER, OPTIONAL, INTENT(OUT) :: code

    FFT%Nx    = -1
    FFT%Ny    = -1
    FFT%Nz    = -1
    FFT%BF    =  1
    FFT%Scale =  1.0_GFT_prec
    FFT%Isign =  0
    FFT%Init  = .FALSE.

    IF ( ASSOCIATED(FFT%Table) ) DEALLOCATE(FFT%Table)
    IF ( ASSOCIATED(FFT%Work) )  DEALLOCATE(FFT%Work)

    IF( PRESENT(code) ) code=0
  END SUBROUTINE GFT_end_rcr

  SUBROUTINE GFT_end_mcc(FFT, code)
    IMPLICIT NONE

    !... Input dummy arguments
    TYPE(GFT_MCC), INTENT(INOUT) :: FFT

    !... Output dummy arguments
    INTEGER, OPTIONAL, INTENT(OUT) :: code

    FFT%Nx    = -1
    FFT%Ny    = -1
    FFT%BF    =  1
    FFT%Scale =  1.0_GFT_prec
    FFT%Isign =  0
    FFT%Init  = .FALSE.

    IF ( ASSOCIATED(FFT%Table) ) DEALLOCATE(FFT%Table)
    IF ( ASSOCIATED(FFT%Work) )  DEALLOCATE(FFT%Work)

    IF( PRESENT(code) ) code=0
  END SUBROUTINE GFT_end_mcc

  SUBROUTINE GFT_end_mrcr(FFT, code)
    IMPLICIT NONE

    !... Input dummy arguments
    TYPE(GFT_MRCR), INTENT(INOUT) :: FFT

    !... Output dummy arguments
    INTEGER, OPTIONAL, INTENT(OUT) :: code

    FFT%Nx    = -1
    FFT%Ny    = -1
    FFT%BF    =  1
    FFT%Scale =  1.0_GFT_prec
    FFT%Isign =  0
    FFT%Init  = .FALSE.

    IF ( ASSOCIATED(FFT%Table) ) DEALLOCATE(FFT%Table)
    IF ( ASSOCIATED(FFT%Work) )  DEALLOCATE(FFT%Work)

    IF( PRESENT(code) ) code=0
  END SUBROUTINE GFT_end_mrcr
END MODULE GFT_end
