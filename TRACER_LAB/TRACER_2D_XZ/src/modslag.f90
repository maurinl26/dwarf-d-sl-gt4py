!#################################  MODELE JOUET 2D  ##################################!
!#                                                                                    #!
!# auteurs : F.Voitus                                                                 #!
!# sujet   : Module servant à faire l'ensemble des calculs dynamiques                 #!
!#                                                                                    #!
!######################################################################################!

!=====================================================================!
!=========================  MODULE PRINCIPAL  ========================!
!=====================================================================!

MODULE MOD_SLAG

  USE MOD_SHARE, ONLY : NXPT,NLEV,RPI,ONE,TWO,ZERO,HALF,LPRINTLEV,RVMAX,RDT,RDX, &
                      & LSETTLS,LNESC,NITMP,NPERIO,LPERIO,LADV_PIS,LCOMAD,RFDC,&
                      & NCOMAD_OPT,C_SLTRAJ,C_SLINTP,LSLTRAJ_PREC
  USE MOD_SETUP, ONLY : SMILAG_STRUCT,RHSVAR,GEOMETRY
 
CONTAINS

!***************************************************************************************!
!***************************************************************************************!
!***************************************************************************************!
!***************************************************************************************!
  SUBROUTINE SMILAG_TRACER_TRANSPORT_SCHEME(RHS_DEP,RHS_ARR,SLSTRUC,            &
                          & PXDOT1,PXDOT0,PXDOT9,PETADOT1,PETADOT0,PETADOT9,    &
                          & PDXDOT1_DX,PDXDOT0_DX,PETA,PDELB,KSTEP,CPART)
  
  IMPLICIT NONE

  TYPE(RHSVAR),                  INTENT(OUT)    :: RHS_DEP
  TYPE(RHSVAR),                  INTENT(IN)     :: RHS_ARR
  TYPE(SMILAG_STRUCT),           INTENT(INOUT)  :: SLSTRUC
  REAL(8),                       INTENT(IN)     :: PXDOT0(NLEV,NXPT),PETADOT0(NLEV,NXPT)
  REAL(8),                       INTENT(IN)     :: PXDOT1(NLEV,NXPT),PETADOT1(NLEV,NXPT)
  REAL(8),                       INTENT(IN)     :: PXDOT9(NLEV,NXPT),PETADOT9(NLEV,NXPT)
  REAL(8),                       INTENT(IN)     :: PDXDOT1_DX(NLEV,NXPT),PDXDOT0_DX(NLEV,NXPT)
  REAL(8),                       INTENT(IN)     :: PETA(NLEV),PDELB(NLEV)
  INTEGER(8),                    INTENT(IN)     :: KSTEP
  CHARACTER(LEN=1),              INTENT(IN)     :: CPART

  REAL(8)                                       :: XDOT_F(NLEV,NXPT),ETADOT_F(NLEV,NXPT)
  REAL(8)                                       :: DXDOT_DX_F(NLEV,NXPT)
  REAL(8)                                       :: ZLIN,ZVMAX,ZWIND
  INTEGER(8)                                    :: I,J

  ZLIN = ONE
  IF (CPART == 'LI') ZLIN =ZERO
 
  IF (LSETTLS) THEN
     XDOT_F(:,:)       = ZLIN*(TWO*PXDOT0(:,:)-PXDOT9(:,:))
     ETADOT_F(:,:)     = ZLIN*(TWO*PETADOT0(:,:)-PETADOT9(:,:))
     DXDOT_DX_F(:,:)   = ZLIN*PDXDOT0_DX(:,:)
  ELSE IF (LNESC) THEN   
     XDOT_F(:,:)       = ZLIN*PXDOT0(:,:)
     ETADOT_F(:,:)     = ZLIN*PETADOT0(:,:)
     DXDOT_DX_F(:,:)   = ZLIN*PDXDOT0_DX(:,:)
  ELSE
     XDOT_F(:,:)       = ZLIN*PXDOT1(:,:)
     ETADOT_F(:,:)     = ZLIN*PETADOT1(:,:)
     DXDOT_DX_F(:,:)   = ZLIN*PDXDOT1_DX(:,:)
  ENDIF
 
  SLSTRUC%XDOT(:,:)    = ZLIN*PXDOT0(:,:)
  SLSTRUC%ETADOT(:,:)  = ZLIN*PETADOT0(:,:)

  SLSTRUC%DXDOT_DX(:,:)= ZERO
  SLSTRUC%DXDOT_DZ(:,:)= ZERO
  SLSTRUC%DZDOT_DX(:,:)= ZERO  
  SLSTRUC%DZDOT_DZ(:,:)= ZERO

  
  !*********************************************************************************!
  IF (LPRINTLEV) THEN                                                               !
     WRITE(*,*) 'NSTEP        : ',KSTEP                                             !
     WRITE(*,*) 'MAX.HOR.WIND : ',MAXVAL(ABS(PXDOT0))                               !
     WRITE(*,*) 'MAX.VER.WIND : ',MAXVAL(ABS(PETADOT0))                             !
  END IF                                                                            !
  IF (MAXVAL(PXDOT0) .GE. RVMAX) THEN                                               !
     DO I=1,NXPT                                                                    !
        DO J=1,NLEV                                                                 !
             ZWIND=PXDOT0(J,I)                                                      !
             IF (ZWIND .GE. RVMAX) THEN                                             !
               WRITE (*,*) '(I,J) = ','(',I,',',J,')',' .... ',ZWIND                !
            END IF                                                                  ! 
         END DO                                                                     !
      END DO                                                                        !
     WRITE(*,*) 'WIND TOO STRONG EXPLOSION AT TIMESTEP = ',KSTEP                    !
     STOP                                                                           !
  END IF                                                                            !
  !*********************************************************************************!
  
  
  !*---------------------------------------------------------------
  ! SL TRAJECTORIES, NEIGHBOURING GRIDPOINT AND ASSOCIATED WEIGHTS    

  IF (LSLTRAJ_PREC) THEN
     CALL COMPUTE_DEFORMATIONAL(SLSTRUC%DXDOT_DX,SLSTRUC%DXDOT_DZ,   &
       & SLSTRUC%DZDOT_DX,SLSTRUC%DZDOT_DZ,SLSTRUC%XDOT,             &
       & SLSTRUC%ETADOT,PETA,1)
  END IF
     
  CALL LARCINA(SLSTRUC,XDOT_F,ETADOT_F,PETA,C_SLTRAJ)

  IF (LCOMAD) THEN
     !Compute comad stretching coefficients
     CALL COMPUTE_COMAD_STRETCHING(SLSTRUC%ALPHA_DX,SLSTRUC%ALPHA_DZ,&
          & XDOT_F,ETADOT_F,SLSTRUC%XDOT,SLSTRUC%ETADOT,PDXDOT0_DX,  &
          & DXDOT_DX_F,SLSTRUC%NXLAG,SLSTRUC%XWEI,SLSTRUC%NZLAG,     &
          & SLSTRUC%ZWEI,PETA,1,C_SLTRAJ)
     !Compute new SL weights for Sources
     CALL ELASCAW(SLSTRUC%NXLAG,SLSTRUC%NZLAG,SLSTRUC%XWEI,          &
          & SLSTRUC%ZWEI,SLSTRUC%XTRAJ,SLSTRUC%ZTRAJ,PETA,1,         &
          & PALPHA_DX=SLSTRUC%ALPHA_DX,PALPHA_DZ=SLSTRUC%ALPHA_DZ,   &
          & LDCOMAD=.TRUE.)
  END IF
  
  CALL LARCINA_SURF(SLSTRUC,XDOT_F,PDXDOT0_DX,DXDOT_DX_F,            &
       & PDELB,LCOMAD)
  
  !*----------------------------------------------------------------
  ! SL INTERPOLATIONS USING NEIGHBOURING GRIDPOINT
  
  CALL LARCINB(RHS_DEP%Q,RHS_ARR%Q,RHS_DEP%M,RHS_ARR%M,              &
       & RHS_DEP%PIS,RHS_ARR%PIS,SLSTRUC%XWEI,SLSTRUC%NXLAG,         &
       & SLSTRUC%ZWEI,SLSTRUC%NZLAG,SLSTRUC%XWEIS,                   &
       & SLSTRUC%NXLAGS,C_SLINTP)
  
END SUBROUTINE SMILAG_TRACER_TRANSPORT_SCHEME
!#################################################################################
!#################################################################################
!#################################################################################
!#################################################################################
!#################################################################################
!#################################################################################    
SUBROUTINE LARCINA(SLSTRUC,PXDOT_ARR,PETADOT_ARR,PETA,CD_SLTRAJ)    

IMPLICIT NONE
               
TYPE(SMILAG_STRUCT),           INTENT(INOUT) :: SLSTRUC
REAL(8), DIMENSION(NLEV,NXPT), INTENT(IN)    :: PXDOT_ARR,PETADOT_ARR
REAL(8), DIMENSION(NLEV),      INTENT(IN)    :: PETA
CHARACTER(LEN=1),              INTENT(IN)    :: CD_SLTRAJ


REAL(8), DIMENSION(NLEV,NXPT)                :: XDOT_DEP,      ETADOT_DEP
REAL(8), DIMENSION(NLEV,NXPT)                :: DXDOT_DX_DEP,  DXDOT_DZ_DEP
REAL(8), DIMENSION(NLEV,NXPT)                :: DZDOT_DX_DEP,  DZDOT_DZ_DEP
INTEGER                                      :: JITER


!initialisation of iterative research
XDOT_DEP(:,:)     = SLSTRUC%XDOT(:,:)
ETADOT_DEP(:,:)   = SLSTRUC%ETADOT(:,:)

DXDOT_DX_DEP(:,:) = SLSTRUC%DXDOT_DX(:,:)
DXDOT_DZ_DEP(:,:) = SLSTRUC%DXDOT_DZ(:,:)
DZDOT_DX_DEP(:,:) = SLSTRUC%DZDOT_DX(:,:)  
DZDOT_DZ_DEP(:,:) = SLSTRUC%DZDOT_DZ(:,:)

DO JITER=1,NITMP
   !Compute SL trajectories
   CALL ELARCHE(SLSTRUC%XTRAJ,SLSTRUC%ZTRAJ,PXDOT_ARR,            &
        & PETADOT_ARR,XDOT_DEP,ETADOT_DEP,PETA,1,                 &
        & PDXDOT_DX_DEP=DXDOT_DX_DEP,PDXDOT_DZ_DEP=DXDOT_DZ_DEP,  &
        & PDZDOT_DX_DEP=DZDOT_DX_DEP,PDZDOT_DZ_DEP=DZDOT_DZ_DEP)
   !Compute SL weights
   CALL ELASCAW(SLSTRUC%NXLAG,SLSTRUC%NZLAG,SLSTRUC%XWEI,         &
        & SLSTRUC%ZWEI,SLSTRUC%XTRAJ,SLSTRUC%ZTRAJ,PETA,1)
   !Interpolations to Departure points
   CALL ELARMES(XDOT_DEP,ETADOT_DEP,DXDOT_DX_DEP,DXDOT_DZ_DEP,    &
        & DZDOT_DX_DEP,DZDOT_DZ_DEP,SLSTRUC%XDOT,SLSTRUC%ETADOT,  &
        & SLSTRUC%DXDOT_DX,SLSTRUC%DXDOT_DZ,SLSTRUC%DZDOT_DX,     &
        & SLSTRUC%DZDOT_DZ,SLSTRUC%NXLAG,SLSTRUC%NZLAG,           &
        & SLSTRUC%XWEI,SLSTRUC%ZWEI,1,CD_SLTRAJ)
END DO

END SUBROUTINE LARCINA
!*******************************************************************************!
!*******************************************************************************!
!*******************************************************************************!
!*******************************************************************************!
!*******************************************************************************!
SUBROUTINE LARCINA_SURF(SLSTRUC,PXDOT_F,PDXDOT_DX,PDXDOT_DX_F,PDELB,LD_COMAD)    

IMPLICIT NONE
               
TYPE(SMILAG_STRUCT),           INTENT(INOUT) :: SLSTRUC
REAL(8), DIMENSION(NLEV,NXPT), INTENT(IN)    :: PXDOT_F
REAL(8), DIMENSION(NLEV,NXPT), INTENT(IN)    :: PDXDOT_DX,PDXDOT_DX_F
REAL(8), DIMENSION(NLEV),      INTENT(IN)    :: PDELB
LOGICAL,                       INTENT(IN)    :: LD_COMAD

REAL(8), DIMENSION(NXPT)                     :: ZXTRAJS,ZALPHAS_DX
INTEGER(8)                                   :: II,JJ

DO II=1,NXPT
   ZXTRAJS(II)    = ZERO
   DO JJ=1,NLEV
      ZXTRAJS(II) = ZXTRAJS(II) &
           & + PDELB(JJ)*SLSTRUC%XTRAJ(JJ,II)
   END DO   
END DO

CALL ELASCAW_SURF(SLSTRUC%NXLAGS,SLSTRUC%XWEIS,ZXTRAJS)

IF (LD_COMAD) THEN
  CALL COMPUTE_COMAD_SURF(ZALPHAS_DX,PXDOT_F,SLSTRUC%XDOT,  &
       & PDXDOT_DX_F,PDXDOT_DX,SLSTRUC%NXLAGS,SLSTRUC%XWEIS,&
       & PDELB,C_SLTRAJ)
  CALL ELASCAW_SURF(SLSTRUC%NXLAGS,SLSTRUC%XWEIS,ZXTRAJS,   &
       & PALPHA=ZALPHAS_DX,LDCOMAD=LCOMAD)  
END IF   

END SUBROUTINE LARCINA_SURF
!************************************************************************************************!
!************************************************************************************************!
!************************************************************************************************!
!************************************************************************************************!
SUBROUTINE ELARMES(PXDOT_DEP,PZDOT_DEP,PDXDOT_DX_DEP,PDXDOT_DZ_DEP,PDZDOT_DX_DEP,PDZDOT_DZ_DEP, &
                 & PXDOT_ARR,PZDOT_ARR,PDXDOT_DX_ARR,PDXDOT_DZ_ARR,PDZDOT_DX_ARR,PDZDOT_DZ_ARR, &
                 & KXLAG,KZLAG,PXWEI,PZWEI,KTOP,CD_SLTRAJ)

IMPLICIT NONE

REAL(8),           INTENT(OUT) :: PXDOT_DEP(KTOP:NLEV,NXPT),     PZDOT_DEP(KTOP:NLEV,NXPT)
REAL(8),           INTENT(OUT) :: PDXDOT_DX_DEP(KTOP:NLEV,NXPT), PDXDOT_DZ_DEP(KTOP:NLEV,NXPT)
REAL(8),           INTENT(OUT) :: PDZDOT_DX_DEP(KTOP:NLEV,NXPT), PDZDOT_DZ_DEP(KTOP:NLEV,NXPT)
REAL(8),           INTENT(IN)  :: PXDOT_ARR(KTOP:NLEV,NXPT),     PZDOT_ARR(KTOP:NLEV,NXPT)
REAL(8),           INTENT(IN)  :: PDXDOT_DX_ARR(KTOP:NLEV,NXPT), PDXDOT_DZ_ARR(KTOP:NLEV,NXPT)
REAL(8),           INTENT(IN)  :: PDZDOT_DX_ARR(KTOP:NLEV,NXPT), PDZDOT_DZ_ARR(KTOP:NLEV,NXPT)
REAL(8),           INTENT(IN)  :: PXWEI(6,KTOP:NLEV,NXPT),       PZWEI(6,KTOP:NLEV,NXPT)
INTEGER(8),        INTENT(IN)  :: KXLAG(6,KTOP:NLEV,NXPT),       KZLAG(6,KTOP:NLEV,NXPT)
INTEGER,           INTENT(IN)  :: KTOP
CHARACTER(LEN=1),  INTENT(IN)  :: CD_SLTRAJ


INTEGER(8)                     :: KK  ,LL  ,II  ,JJ
INTEGER                        :: ISTA  ,IEND

IF (CD_SLTRAJ == "C") THEN
   ISTA = 1
   IEND = 4
ELSE IF (CD_SLTRAJ == "L") THEN
   ISTA = 5
   IEND = 6
ELSE
   WRITE(*,*) "UNKNOWN TRAJECTORY INTERPOLATION ORDER : ",CD_SLTRAJ
   STOP
END IF   

DO II=1,NXPT
   DO JJ=KTOP,NLEV
      PXDOT_DEP(JJ,II) = ZERO
      PZDOT_DEP(JJ,II) = ZERO
      DO KK=ISTA,IEND
         DO LL=ISTA,IEND
            PXDOT_DEP(JJ,II) = PXDOT_DEP(JJ,II)                            &
                             & + PZWEI(LL,JJ,II)*PXWEI(KK,JJ,II)*(         &
                             &   PXDOT_ARR(KZLAG(LL,JJ,II),KXLAG(KK,JJ,II)))
            PZDOT_DEP(JJ,II) = PZDOT_DEP(JJ,II)                            &
                             & + PZWEI(LL,JJ,II)*PXWEI(KK,JJ,II)*(         &
                             &   PZDOT_ARR(KZLAG(LL,JJ,II),KXLAG(KK,JJ,II)))
         ENDDO
      ENDDO
   ENDDO
ENDDO

!* set deformational term to zero
PDXDOT_DX_DEP(:,:) = ZERO
PDXDOT_DZ_DEP(:,:) = ZERO
PDZDOT_DX_DEP(:,:) = ZERO
PDZDOT_DZ_DEP(:,:) = ZERO

IF (LSLTRAJ_PREC) THEN
  DO II=1,NXPT
     DO JJ=KTOP,NLEV
        DO KK=ISTA,IEND
           DO LL=ISTA,IEND
              PDXDOT_DX_DEP(JJ,II) = PDXDOT_DZ_DEP(JJ,II)                            &
                                   & + PZWEI(LL,JJ,II)*PXWEI(KK,JJ,II)*(             &
                                   &   PDXDOT_DX_ARR(KZLAG(LL,JJ,II),KXLAG(KK,JJ,II)))
              PDXDOT_DZ_DEP(JJ,II) = PDXDOT_DZ_DEP(JJ,II)                            &
                                   & + PZWEI(LL,JJ,II)*PXWEI(KK,JJ,II)*(             &
                                   &   PDXDOT_DZ_ARR(KZLAG(LL,JJ,II),KXLAG(KK,JJ,II)))
              PDZDOT_DX_DEP(JJ,II) = PDZDOT_DZ_DEP(JJ,II)                            &
                                   & + PZWEI(LL,JJ,II)*PXWEI(KK,JJ,II)*(             &
                                   &   PDZDOT_DX_ARR(KZLAG(LL,JJ,II),KXLAG(KK,JJ,II)))
              PDZDOT_DZ_DEP(JJ,II) = PDZDOT_DZ_DEP(JJ,II)                            &
                                   & + PZWEI(LL,JJ,II)*PXWEI(KK,JJ,II)*(             &
                                   &   PDZDOT_DZ_ARR(KZLAG(LL,JJ,II),KXLAG(KK,JJ,II)))
           ENDDO
        ENDDO
     ENDDO
  ENDDO
END IF

END SUBROUTINE ELARMES
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
SUBROUTINE ELARCHE(PXTRAJ,PZTRAJ,PXDOT_ARR,PETADOT_ARR,PXDOT_DEP,PETADOT_DEP, &
                 & PETA,KTOP,PDXDOT_DX_DEP,PDXDOT_DZ_DEP,PDZDOT_DX_DEP,PDZDOT_DZ_DEP)
    
IMPLICIT NONE

REAL(8),           INTENT(INOUT) :: PXTRAJ(KTOP:NLEV,NXPT),       PZTRAJ(KTOP:NLEV,NXPT)       
REAL(8),           INTENT(IN)    :: PXDOT_ARR(KTOP:NLEV,NXPT),    PETADOT_ARR(KTOP:NLEV,NXPT)
REAL(8),           INTENT(IN)    :: PXDOT_DEP(KTOP:NLEV,NXPT),    PETADOT_DEP(KTOP:NLEV,NXPT)
REAL(8), OPTIONAL, INTENT(IN)    :: PDXDOT_DX_DEP(KTOP:NLEV,NXPT),PDZDOT_DZ_DEP(KTOP:NLEV,NXPT)
REAL(8), OPTIONAL, INTENT(IN)    :: PDXDOT_DZ_DEP(KTOP:NLEV,NXPT),PDZDOT_DX_DEP(KTOP:NLEV,NXPT)
REAL(8),           INTENT(IN)    :: PETA(KTOP:NLEV)
INTEGER,           INTENT(IN)    :: KTOP

INTEGER                          :: II         ,JJ
REAL(8)                          :: ZVETAOX  ,ZVETAON  ,ZTRAJZ_UBND    ,ZTRAJX  ,ZTRAJZ
REAL(8)                          :: ZEPS     ,ZDET     ,ZDET0_C  


ZDET0_C = ABS(1000*EPSILON(ONE))

ZVETAON=PETA(KTOP)
ZVETAOX=PETA(NLEV)

!-------------------------------
! TRAJECTORY RESEARCH
!-------------------------------

IF (LSLTRAJ_PREC) THEN
   DO II=1,NXPT
      DO JJ=KTOP,NLEV
         ZDET             = ONE+(RDT/TWO)*(PDXDOT_DX_DEP(JJ,II)+PDZDOT_DZ_DEP(JJ,II))   &
                          & + ((RDT/TWO)**2)*(PDXDOT_DX_DEP(JJ,II)*PDZDOT_DZ_DEP(JJ,II) &
                          & - PDXDOT_DZ_DEP(JJ,II)*PDZDOT_DX_DEP(JJ,II))
         IF (ABS(ZDET) > ZDET0_C) THEN
            ZTRAJX        = REAL(II-1,8)*RDX                                            &
                          & - (RDT/TWO)*( PXDOT_ARR(JJ,II)+PXDOT_DEP(JJ,II) )           &
                          & + (RDT/TWO)*( PDXDOT_DX_DEP(JJ,II)*PXTRAJ(JJ,II)            &
                          & + PDXDOT_DZ_DEP(JJ,II)*PZTRAJ(JJ,II) ) 
            ZTRAJZ        = PETA(JJ)                                                    &
                          & - (RDT/TWO)*(PETADOT_ARR(JJ,II)+PETADOT_DEP(JJ,II))         &
                          & + (RDT/TWO)*( PDZDOT_DZ_DEP(JJ,II)*PZTRAJ(JJ,II)            &
                          & + PDZDOT_DX_DEP(JJ,II)*PXTRAJ(JJ,II) )
            ZTRAJZ_UBND   = ( (ONE+(RDT/TWO)*PDXDOT_DX_DEP(JJ,II))*ZTRAJZ               &
                          & - (RDT/TWO)*PDZDOT_DX_DEP(JJ,II)*ZTRAJX )/ZDET
            PXTRAJ(JJ,II) = ( (ONE+(RDT/TWO)*PDZDOT_DZ_DEP(JJ,II))*ZTRAJX               &
                          & - (RDT/TWO)*PDXDOT_DZ_DEP(JJ,II)*ZTRAJZ )/ZDET
            PZTRAJ(JJ,II) = MAX(ZVETAON,MIN(ZVETAOX,ZTRAJZ_UBND))
         ELSE
            PXTRAJ(JJ,II) = REAL(II-1,8)*RDX                                            &
                          & - (RDT/TWO)*(PXDOT_ARR(JJ,II)+PXDOT_DEP(JJ,II))
            ZTRAJZ_UBND   = PETA(JJ)                                                    &
                          & - (RDT/TWO)*(PETADOT_ARR(JJ,II)+PETADOT_DEP(JJ,II) )
            PZTRAJ(JJ,II) = MAX(ZVETAON,MIN(ZVETAOX,ZTRAJZ_UBND))
         END IF   
      ENDDO
   ENDDO
ELSE
   DO II=1,NXPT
      DO JJ=KTOP,NLEV
         PXTRAJ(JJ,II)    = REAL(II-1,8)*RDX - (RDT/TWO)*( PXDOT_ARR(JJ,II)             &
                          & + PXDOT_DEP(JJ,II) )
         ZTRAJZ_UBND      = PETA(JJ) - (RDT/TWO)*( PETADOT_ARR(JJ,II)                   &
                          & + PETADOT_DEP(JJ,II) )
         PZTRAJ(JJ,II)    = MAX(ZVETAON,MIN(ZVETAOX,ZTRAJZ_UBND))
      ENDDO
   ENDDO
END IF   
   
END SUBROUTINE ELARCHE
!**********************************************************************************************!
!**********************************************************************************************!
!**********************************************************************************************!
!**********************************************************************************************!
!**********************************************************************************************!
!**********************************************************************************************!
!**********************************************************************************************!
!**********************************************************************************************!
SUBROUTINE ELASCAW(KXLAG,KZLAG,PXWEI,PZWEI,PXTRAJ,PZTRAJ,PETA,KTOP, &
                 & PALPHA_DX,PALPHA_DZ,LDCOMAD)

IMPLICIT NONE

REAL(8),           INTENT(OUT)   :: PXWEI(6,KTOP:NLEV,NXPT), PZWEI(6,KTOP:NLEV,NXPT)
INTEGER(8),        INTENT(OUT)   :: KXLAG(6,KTOP:NLEV,NXPT), KZLAG(6,KTOP:NLEV,NXPT)
REAL(8),           INTENT(IN)    :: PXTRAJ(KTOP:NLEV,NXPT),  PZTRAJ(KTOP:NLEV,NXPT)
REAL(8),           INTENT(IN)    :: PETA(KTOP:NLEV)
INTEGER,           INTENT(IN)    :: KTOP
REAL(8), OPTIONAL, INTENT(IN)    :: PALPHA_DX(KTOP:NLEV,NXPT),PALPHA_DZ(KTOP:NLEV,NXPT)
LOGICAL, OPTIONAL, INTENT(IN)    :: LDCOMAD

INTEGER(8)                       :: IP       ,II       ,JJ       ,JP       ,KK      ,LL    ,J0
REAL(8)                          :: ZXP1     ,ZX0      ,ZXM1     ,ZXM2     ,ZXD     ,ZALX  ,ZWI
REAL(8)                          :: ZZP1     ,ZZ0      ,ZZM1     ,ZZM2     ,ZZD     ,ZAL   ,ZWJ
REAL(8)                          :: ZEPS     ,ZDIFF    ,ZSIGN    ,ZWA      ,ZWB     ,ZWC   ,ZWD
LOGICAL                          :: LL_COMAD

ZEPS = ABS(EPSILON(ONE))

IF (PRESENT(LDCOMAD)) THEN
   LL_COMAD = LDCOMAD 
ELSE
   LL_COMAD = .FALSE.
END IF

!--------------------------------------------------------------------------------------------------
! COMPUTE WEIGHTS FOR SL HORIZONTAL INTERPOLATIONS
IF (LL_COMAD) THEN
   
   IF (LPERIO) THEN      
      DO II=1,NXPT
         DO JJ=KTOP,NLEV
         
            ZXD  = PXTRAJ(JJ,II) 
            ZALX = REAL(II-1,8)-(ZXD/RDX)
            IF (ABS(ZALX).LT.ZEPS) ZALX = ZERO
        
            IP = FLOOR(ZALX)

            ZXM1= REAL(II-IP-1-1,8)*RDX
            ZWI = (ZXD-ZXM1)/RDX
            ZWC = ZWI+(ZWI-HALF)*(PALPHA_DX(JJ,II)-ONE)

            !-------------------------------------------------------------------------------------
            ! Cubic interpolating weights between II-IP-2, II-IP-1, II-IP and II-IP+1
            DO KK=1,4
               KXLAG(KK,JJ,II) = NPERIO(II-IP+2-KK) 
            END DO
            PXWEI(1,JJ,II)=((ZWC**2)-ONE)*(ZWC/6.d0)
            PXWEI(2,JJ,II)=ZWC+((ZWC**2)/TWO)*(ONE-ZWC)
            PXWEI(3,JJ,II)=ONE-(ZWC**2)+(ZWC/TWO)*((ZWC**2)-ONE)
            PXWEI(4,JJ,II)=ONE-PXWEI(1,JJ,II)-PXWEI(2,JJ,II)-PXWEI(3,JJ,II)
            !-------------------------------------------------------------------------------------
            ! Linear interpolating weights between II-IP-1, II-IP
            DO KK=5,6
               KXLAG(KK,JJ,II) = NPERIO(II-IP+5-KK) 
            END DO
            PXWEI(5,JJ,II)=ZWC
            PXWEI(6,JJ,II)=ONE-ZWC
            !-------------------------------------------------------------------------------------
         END DO  
      END DO   
   ELSE
      DO II=1,NXPT
         DO JJ=KTOP,NLEV
        
            ZXD  = PXTRAJ(JJ,II) 
            ZALX = REAL(II-1,8)-(ZXD/RDX)
            IF (ABS(ZALX).LT.ZEPS) ZALX = ZERO
            IP   = FLOOR(ZALX)

            !*-------------------------------------------------------------------------------------
            ! cubic interpolation between II-IP-2, II-IP-1, II-IP and II-IP+1
            DO KK=1,4
               KXLAG(KK,JJ,II) = MIN(NXPT,MAX(II-IP+2-KK,1)) 
            END DO
            IF (((II-IP+1).LE.NXPT).AND.((II-IP-2).GE.1)) THEN
               ZXM1=REAL(II-IP-1-1,8)*RDX
               ZWI =(ZXD-ZXM1)/RDX
               ZWC =ZWI+(ZWI-HALF)*(PALPHA_DX(JJ,II)-ONE)
               PXWEI(1,JJ,II)=((ZWC**2)-ONE)*(ZWC/6.d0)
               PXWEI(2,JJ,II)=ZWC+((ZWC**2)/TWO)*(ONE-ZWC)
               PXWEI(3,JJ,II)=ONE-(ZWC**2)+(ZWC/TWO)*((ZWC**2)-ONE)
               PXWEI(4,JJ,II)=ONE-PXWEI(1,JJ,II)-PXWEI(2,JJ,II)-PXWEI(3,JJ,II)
            ELSE IF (II-IP.LE.1) THEN
            ! trajectory truncation at x=0
               PXWEI(1,JJ,II)=ZERO
               PXWEI(2,JJ,II)=ONE
               PXWEI(3,JJ,II)=ZERO
               PXWEI(4,JJ,II)=ZERO
            ELSE IF (II-IP-1.EQ.1) THEN
            ! linear interpolation between II-IP and II-IP+1 at x=0
               ZXM1=REAL(II-IP-1-1,8)*RDX
               ZWI =(ZXD-ZXM1)/RDX
               ZWC =ZWI+(ZWI-HALF)*(PALPHA_DX(JJ,II)-ONE)
               PXWEI(1,JJ,II)=ZERO
               PXWEI(2,JJ,II)=ZWC
               PXWEI(3,JJ,II)=ONE-ZWC
               PXWEI(4,JJ,II)=ZERO
            ELSE IF (II-IP-1.GE.NXPT) THEN
            ! trajectory truncation at x=L
               PXWEI(1,JJ,II)=ZERO
               PXWEI(2,JJ,II)=ZERO
               PXWEI(3,JJ,II)=ONE
               PXWEI(4,JJ,II)=ZERO
            ELSE IF (II-IP.EQ.NXPT) THEN
            ! linear interpolation between II-IP-2 and II-IP-1 at x=L
               ZXM1=REAL(II-IP-1-1,8)*RDX
               ZWI =(ZXD-ZXM1)/RDX
               ZWC =ZWI+(ZWI-HALF)*(PALPHA_DX(JJ,II)-ONE)
               PXWEI(1,JJ,II)=ZERO
               PXWEI(2,JJ,II)=ZWC
               PXWEI(3,JJ,II)=ONE-ZWC
               PXWEI(4,JJ,II)=ZERO
            ELSE
               WRITE(*,*) 'bug in the code in the interpolation'
               STOP
            END IF
        
            !*---------------------------------------------------------
            ! Linear interpolation between II-IP-1, and II-IP 
            DO KK=5,6
               KXLAG(KK,JJ,II) = MIN(NXPT,MAX(II-IP+5-KK,1)) 
            END DO
            IF (((II-IP).LE.NXPT).AND.((II-IP-1).GE.1)) THEN
               ZXM1=REAL(II-IP-1-1,8)*RDX
               ZWI =(ZXD-ZXM1)/RDX
               ZWC =ZWI+(ZWI-HALF)*(PALPHA_DX(JJ,II)-ONE)
               PXWEI(5,JJ,II)=ZWC
               PXWEI(6,JJ,II)=ONE-ZWC
            ELSE IF (II-IP.LE.1) THEN
               PXWEI(5,JJ,II)=ONE
               PXWEI(6,JJ,II)=ZERO
            ELSE IF (II-IP-1.GE.NXPT) THEN
               PXWEI(5,JJ,II)=ZERO
               PXWEI(6,JJ,II)=ONE
            END IF
         END DO  
      END DO
   END IF
   
ELSE  !**** NOT COMAD *****!
   
   IF (LPERIO) THEN      
      DO II=1,NXPT
         DO JJ=KTOP,NLEV
         
            ZXD  = PXTRAJ(JJ,II) 
            ZALX = REAL(II-1,8)-(ZXD/RDX)
            IF (ABS(ZALX).LT.ZEPS) ZALX = ZERO
        
            IP = FLOOR(ZALX)

            DO KK=1,4
               KXLAG(KK,JJ,II) = NPERIO(II-IP+2-KK) 
            END DO
         
            ZXP1=REAL(II-IP+1-1,8)*RDX
            ZX0 =REAL(II-IP+0-1,8)*RDX
            ZXM1=REAL(II-IP-1-1,8)*RDX
            ZXM2=REAL(II-IP-2-1,8)*RDX

            !---------------------------------------------------------------------------------------
            ! Cubic interpolating weights between II-IP-2, II-IP-1, II-IP and II-IP+1
            PXWEI(1,JJ,II)=((ZXD-ZX0)/(ZXP1-ZX0))*((ZXD-ZXM1)/(ZXP1-ZXM1))*((ZXD-ZXM2)/(ZXP1-ZXM2))
            PXWEI(2,JJ,II)=((ZXD-ZXP1)/(ZX0-ZXP1))*((ZXD-ZXM1)/(ZX0-ZXM1))*((ZXD-ZXM2)/(ZX0-ZXM2))
            PXWEI(3,JJ,II)=((ZXD-ZXP1)/(ZXM1-ZXP1))*((ZXD-ZX0)/(ZXM1-ZX0))*((ZXD-ZXM2)/(ZXM1-ZXM2))
            PXWEI(4,JJ,II)=((ZXD-ZXP1)/(ZXM2-ZXP1))*((ZXD-ZX0)/(ZXM2-ZX0))*((ZXD-ZXM1)/(ZXM2-ZXM1))

            !---------------------------------------------------------------------------------------
            ! Linear interpolating weights between II-IP-1, II-IP
            DO KK=5,6
               KXLAG(KK,JJ,II) = NPERIO(II-IP+5-KK) 
            END DO
            PXWEI(5,JJ,II)=(ZXD-ZXM1)/(ZX0-ZXM1)
            PXWEI(6,JJ,II)=(ZXD-ZX0)/(ZXM1-ZX0)
            !---------------------------------------------------------------------------------------
         END DO  
      END DO
   ELSE !****NOT LPERIO****!
      DO II=1,NXPT
         DO JJ=KTOP,NLEV
        
            ZXD   = PXTRAJ(JJ,II) 
            ZALX  = REAL(II-1,8)-(ZXD/RDX)
            IF (ABS(ZALX).LT.ZEPS) ZALX = ZERO
         
            IP = FLOOR(ZALX)
         
            DO KK=1,4
               KXLAG(KK,JJ,II) = MIN(NXPT,MAX(II-IP+2-KK,1)) 
            END DO
            IF (((II-IP+1).LE.NXPT).AND.((II-IP-2).GE.1)) THEN
              ! cubic interpolation between II-IP-2, II-IP-1, II-IP and II-IP+1
              ZXP1=REAL(II-IP+1-1,8)*RDX
              ZX0 =REAL(II-IP+0-1,8)*RDX
              ZXM1=REAL(II-IP-1-1,8)*RDX
              ZXM2=REAL(II-IP-2-1,8)*RDX
              PXWEI(1,JJ,II)=((ZXD-ZX0)/(ZXP1-ZX0))*((ZXD-ZXM1)/(ZXP1-ZXM1))*((ZXD-ZXM2)/(ZXP1-ZXM2))
              PXWEI(2,JJ,II)=((ZXD-ZXP1)/(ZX0-ZXP1))*((ZXD-ZXM1)/(ZX0-ZXM1))*((ZXD-ZXM2)/(ZX0-ZXM2))
              PXWEI(3,JJ,II)=((ZXD-ZXP1)/(ZXM1-ZXP1))*((ZXD-ZX0)/(ZXM1-ZX0))*((ZXD-ZXM2)/(ZXM1-ZXM2))
              PXWEI(4,JJ,II)=((ZXD-ZXP1)/(ZXM2-ZXP1))*((ZXD-ZX0)/(ZXM2-ZX0))*((ZXD-ZXM1)/(ZXM2-ZXM1))
            ELSE IF (II-IP.LE.1) THEN
              ! trajectory truncation at x=0
              PXWEI(1,JJ,II)=ZERO
              PXWEI(2,JJ,II)=ONE
              PXWEI(3,JJ,II)=ZERO
              PXWEI(4,JJ,II)=ZERO
            ELSE IF (II-IP-1.EQ.1) THEN
              ! linear interpolation between II-IP and II-IP+1 at x=0
              ZXM1=REAL(II-IP-1-1,8)*RDX
              ZX0 =REAL(II-IP+0-1,8)*RDX
              PXWEI(1,JJ,II)=ZERO
              PXWEI(2,JJ,II)=(ZXD-ZXM1)/(ZX0-ZXM1)
              PXWEI(3,JJ,II)=(ZXD-ZX0)/(ZXM1-ZX0)
              PXWEI(4,JJ,II)=ZERO
            ELSE IF (II-IP-1.GE.NXPT) THEN
              ! trajectory truncation at x=L
              PXWEI(1,JJ,II)=ZERO
              PXWEI(2,JJ,II)=ZERO
              PXWEI(3,JJ,II)=ONE
              PXWEI(4,JJ,II)=ZERO
            ELSE IF (II-IP.EQ.NXPT) THEN
              ! linear interpolation between II-IP-2 and II-IP-1 at x=L
              ZXM1=REAL(II-IP-1-1,8)*RDX
              ZX0=REAL(II-IP-1,8)*RDX
              PXWEI(1,JJ,II)=ZERO
              PXWEI(2,JJ,II)=(ZXD-ZXM1)/(ZX0-ZXM1)
              PXWEI(3,JJ,II)=(ZXD-ZX0)/(ZXM1-ZX0)
              PXWEI(4,JJ,II)=ZERO
            ELSE
              WRITE(*,*) 'bug in the code in the interpolation'
              STOP
            END IF

            ! Linear interpolation between II-IP-1, and II-IP 
            DO KK=5,6
                KXLAG(KK,JJ,II) = MIN(NXPT,MAX(II-IP+5-KK,1)) 
            END DO
            IF (((II-IP).LE.NXPT).AND.((II-IP-1).GE.1)) THEN
              ZX0 =REAL(II-IP+0-1,8)*RDX
              ZXM1=REAL(II-IP-1-1,8)*RDX
              PXWEI(5,JJ,II)=(ZXD-ZXM1)/(ZX0-ZXM1)
              PXWEI(6,JJ,II)=(ZXD-ZX0)/(ZXM1-ZX0)
            ELSE IF (II-IP.LE.1) THEN
              PXWEI(5,JJ,II)=ONE
              PXWEI(6,JJ,II)=ZERO
            ELSE IF (II-IP-1.GE.NXPT) THEN
              PXWEI(5,JJ,II)=ZERO
              PXWEI(6,JJ,II)=ONE
            END IF  
         END DO  
      END DO
   END IF 

END IF
   
         !***********************************!
         !* VERTICAL INTERPOLATIONS WEIGHTS *!
         !***********************************!

IF (LL_COMAD) THEN
   DO II=1,NXPT
      DO JJ=KTOP,NLEV
         ZZD   = PZTRAJ(JJ,II) 
         ZAL   = PETA(JJ)-PZTRAJ(JJ,II)
         IF (ABS(ZAL).LT.ZEPS) THEN
            ZAL = ZERO
            JP  = JJ
         ELSE
            JP = JJ
            ZSIGN = MAX(ZERO,SIGN(ONE,ZAL))
            IF (ZSIGN == ABS(ONE)) THEN
               DO J0=JJ,KTOP+1,-1
                  ZDIFF = MAX(ZERO,SIGN(ONE,(ZZD-PETA(J0))*(ZZD-PETA(J0-1))))
                  IF (ZDIFF==ZERO) THEN
                     JP=J0
                     EXIT
                  END IF
                  JP=KTOP
               END DO
            ELSE
               JP = JJ
               DO J0=JJ+1,NLEV
                  ZDIFF = MAX(ZERO,SIGN(ONE,(ZZD-PETA(J0))*(ZZD-PETA(J0-1))))
                  IF (ZDIFF==ZERO) THEN
                     JP=J0
                     EXIT   
                  END IF
                  JP = NLEV
               END DO  
            END IF
         END IF

         !--------------------------------------------------------------------
         ! cubic interpolation between JJ-JP-2, JJ-JP-1, JJ-IP and JJ-JP+1
         DO LL=1,4
            KZLAG(LL,JJ,II) = MAX(KTOP,MIN(NLEV,JP+2-LL)) 
         END DO
         IF ((JP+1.LE.NLEV).AND.(JP-2.GE.KTOP)) THEN

            ZZP1=PETA(JP+1)
            ZZ0 =PETA(JP+0)
            ZZM1=PETA(JP-1)
            ZZM2=PETA(JP-2)
            
            ZWJ=(ZZD-ZZM1)/(ZZ0-ZZM1) 
            ZWA=(ZZM1-ZZM2)/(ZZ0-ZZM1)
            ZWB=(ZZP1-ZZM1)/(ZZ0-ZZM1)
            ZWC=ZWJ+(ZWJ-HALF)*(PALPHA_DZ(JJ,II)-ONE)
            ZWD=ONE/((ONE+ZWA)*(ONE+ZWB)*(ZWA-ZWB)) 
           
            PZWEI(1,JJ,II)=ZWD*(((ONE-(ZWA**2))/ZWB)*ZWC*(ZWC-ONE)                  &
                 & + ZWC*(ONE-(ZWC**2))*((ONE+ZWA)/ZWB))
            PZWEI(2,JJ,II)=ZWC+ZWD*(ZWC*(ZWC-ONE)*((ONE-(ZWB**2))-(ONE-(ZWA**2)))   &
                 & + ZWC*((ZWC**2)-ONE)*(ZWA+ZWB))
            PZWEI(3,JJ,II)=(ONE-ZWC)+ZWD*( ZWC*(ZWC-ONE)*(                          &
                 & (ONE-(ZWA**2))*(ONE-(ONE/ZWB))-(ONE-(ZWB**2))*(ONE+(ONE/ZWA)))   &
                 & + ZWC*(ONE-(ZWC**2))*(((ONE-ZWB)/ZWA)+((ONE+ZWA)/ZWB)-(ZWA+ZWB)) )
            PZWEI(4,JJ,II)=ONE-PZWEI(1,JJ,II)-PZWEI(2,JJ,II)-PZWEI(3,JJ,II)
            
         ELSEIF (JP-1.EQ.KTOP) THEN
         ! linear interpolation between JJ-JP-1, JJ-JP and JJ-JP+1 at x=0
            ZZM1=PETA(JP-1)
            ZZ0 =PETA(JP+0)
            ZWJ =(ZZD-ZZM1)/(ZZ0-ZZM1)
            ZWC =ZWJ+(ZWJ-HALF)*(PALPHA_DZ(JJ,II)-ONE)
            PZWEI(1,JJ,II)=ZERO
            PZWEI(2,JJ,II)=ZWC
            PZWEI(3,JJ,II)=ONE-ZWC
            PZWEI(4,JJ,II)=ZERO
         ELSEIF (JP.LE.KTOP) THEN
         ! Trajectory truncation at x=L
            PZWEI(1,JJ,II)=ZERO
            PZWEI(2,JJ,II)=ONE
            PZWEI(3,JJ,II)=ZERO
            PZWEI(4,JJ,II)=ZERO    
         ELSEIF (JP.EQ.NLEV) THEN
         ! linear interpolation between JJ-JP-2, JJ-JP-1 and JJ-JP at x=L
            ZZ0 =PETA(JP+0)
            ZZM1=PETA(JP-1)
            ZWJ =(ZZD-ZZM1)/(ZZ0-ZZM1)
            ZWC =ZWJ+(ZWJ-HALF)*(PALPHA_DZ(JJ,II)-ONE)
            PZWEI(1,JJ,II)=ZERO
            PZWEI(2,JJ,II)=ZWC
            PZWEI(3,JJ,II)=ONE-ZWC
            PZWEI(4,JJ,II)=ZERO
         ELSEIF (JP-1.GE.NLEV) THEN
         ! Trajectory truncation at x=L
            PZWEI(1,JJ,II)=ZERO
            PZWEI(2,JJ,II)=ZERO
            PZWEI(3,JJ,II)=ONE
            PZWEI(4,JJ,II)=ZERO
         ENDIF

         !--------------------------------------------------------------------
         !... Linear vertical interpolating
         DO LL=5,6
            KZLAG(LL,JJ,II) = MAX(KTOP,MIN(NLEV,JP+5-LL)) 
         END DO
         IF ((JP.LE.NLEV).AND.(JP-1.GE.KTOP)) THEN
            ZZ0 =PETA(JP+0)
            ZZM1=PETA(JP-1)
            ZWJ =(ZZD-ZZM1)/(ZZ0-ZZM1)
            ZWC =ZWJ+(ZWJ-HALF)*(PALPHA_DZ(JJ,II)-ONE)
            PZWEI(5,JJ,II)=ZWC
            PZWEI(6,JJ,II)=ONE-ZWC
         ELSEIF (JP.LE.KTOP) THEN
         ! Trajectory truncation at x=1
            PZWEI(5,JJ,II)=ONE
            PZWEI(6,JJ,II)=ZERO     
         ELSEIF (JP-1.GE.NLEV) THEN
         ! Trajectory truncation at x=L
            PZWEI(5,JJ,II)=ZERO 
            PZWEI(6,JJ,II)=ONE
         END IF
         !--------------------------------------------------------------------
      END DO
   END DO
   
ELSE !****** NOT COMAD *****!
   
   DO II=1,NXPT
      DO JJ=KTOP,NLEV
         ZZD   = PZTRAJ(JJ,II) 
         ZAL   = PETA(JJ)-PZTRAJ(JJ,II)
         IF (ABS(ZAL).LT.ZEPS) THEN
            ZAL = ZERO
            JP  = JJ
         ELSE
            JP = JJ
            ZSIGN = MAX(ZERO,SIGN(ONE,ZAL))
            IF (ZSIGN == ABS(ONE)) THEN
              DO J0=JJ,KTOP+1,-1
                 ZDIFF = MAX(ZERO,SIGN(ONE,(ZZD-PETA(J0))*(ZZD-PETA(J0-1))))
                 IF (ZDIFF==ZERO) THEN
                    JP=J0
                    EXIT
                 END IF
                 JP=KTOP
              END DO
            ELSE
              JP = JJ
              DO J0=JJ+1,NLEV
                 ZDIFF = MAX(ZERO,SIGN(ONE,(ZZD-PETA(J0))*(ZZD-PETA(J0-1))))
                 IF (ZDIFF==ZERO) THEN
                    JP=J0
                    EXIT   
                 END IF
                 JP = NLEV
              END DO  
            END IF
         END IF

         !*-------------------------------------------------------------------------------------
         ! cubic interpolation between JJ-JP-2, JJ-JP-1, JJ-IP and JJ-JP+1
         DO LL=1,4
            KZLAG(LL,JJ,II) = MAX(KTOP,MIN(NLEV,JP+2-LL)) 
         END DO
         IF ((JP+1.LE.NLEV).AND.(JP-2.GE.KTOP)) THEN
           ZZP1=PETA(JP+1)
           ZZ0 =PETA(JP+0)
           ZZM1=PETA(JP-1)
           ZZM2=PETA(JP-2)
           PZWEI(1,JJ,II)=((ZZD-ZZ0)/(ZZP1-ZZ0))*((ZZD-ZZM1)/(ZZP1-ZZM1))*((ZZD-ZZM2)/(ZZP1-ZZM2))
           PZWEI(2,JJ,II)=((ZZD-ZZP1)/(ZZ0-ZZP1))*((ZZD-ZZM1)/(ZZ0-ZZM1))*((ZZD-ZZM2)/(ZZ0-ZZM2))
           PZWEI(3,JJ,II)=((ZZD-ZZP1)/(ZZM1-ZZP1))*((ZZD-ZZ0)/(ZZM1-ZZ0))*((ZZD-ZZM2)/(ZZM1-ZZM2))
           PZWEI(4,JJ,II)=((ZZD-ZZP1)/(ZZM2-ZZP1))*((ZZD-ZZ0)/(ZZM2-ZZ0))*((ZZD-ZZM1)/(ZZM2-ZZM1)) 
         ELSEIF (JP-1.EQ.KTOP) THEN
         ! linear interpolation between JJ-JP-1, JJ-JP  at eta=1
           ZZ0 =PETA(JP+0)
           ZZM1=PETA(JP-1)
           PZWEI(1,JJ,II)=ZERO
           PZWEI(2,JJ,II)=(ZZD-ZZM1)/(ZZ0-ZZM1)
           PZWEI(3,JJ,II)=(ZZD-ZZ0)/(ZZM1-ZZ0)
           PZWEI(4,JJ,II)=ZERO
         ELSEIF (JP.LE.KTOP) THEN
           ! Trajectory truncation at eta=1
           PZWEI(1,JJ,II)=ZERO
           PZWEI(2,JJ,II)=ONE
           PZWEI(3,JJ,II)=ZERO
           PZWEI(4,JJ,II)=ZERO    
         ELSEIF (JP.EQ.NLEV) THEN
           ! linear interpolation between JJ-JP-1 and JJ-JP at eta=L
           ZZ0 =PETA(JP+0)
           ZZM1=PETA(JP-1)
           PZWEI(1,JJ,II)=ZERO
           PZWEI(2,JJ,II)=(ZZD-ZZM1)/(ZZ0-ZZM1)
           PZWEI(3,JJ,II)=(ZZD-ZZ0)/(ZZM1-ZZ0)
           PZWEI(4,JJ,II)=ZERO  
         ELSEIF (JP-1.GE.NLEV) THEN
           ! Trajectory truncation at eta=L
           PZWEI(1,JJ,II)=ZERO
           PZWEI(2,JJ,II)=ZERO
           PZWEI(3,JJ,II)=ONE
           PZWEI(4,JJ,II)=ZERO
         ENDIF

         !---------------------------------------------------------------------
         !... Linear vertical interpolating wieghts
         DO LL=5,6
            KZLAG(LL,JJ,II) = MAX(KTOP,MIN(NLEV,JP+5-LL)) 
         END DO
         IF ((JP.LE.NLEV).AND.(JP-1.GE.KTOP)) THEN
           ZZ0 =PETA(JP+0)
           ZZM1=PETA(JP-1)
           PZWEI(5,JJ,II)=(ZZD-ZZM1)/(ZZ0-ZZM1)
           PZWEI(6,JJ,II)=(ZZD-ZZ0)/(ZZM1-ZZ0)
         ELSEIF (JP.LE.KTOP) THEN
           ! Trajectory truncation at eta=1
           PZWEI(5,JJ,II)=ONE
           PZWEI(6,JJ,II)=ZERO     
         ELSEIF (JP-1.GE.NLEV) THEN
           ! Trajectory truncation at eta=L
           PZWEI(5,JJ,II)=ZERO 
           PZWEI(6,JJ,II)=ONE
         END IF
         !*---------------------------------------------------------------------
      END DO
   END DO
END IF   


END SUBROUTINE ELASCAW
!*******************************************************************************!
!*******************************************************************************!
!*******************************************************************************!
!*******************************************************************************!
SUBROUTINE ELASCAW_SURF(KXLAG,PXWEI,PXTRAJ,PALPHA,LDCOMAD)

IMPLICIT NONE

REAL(8),                INTENT(OUT)   :: PXWEI(6,NXPT)
INTEGER(8),             INTENT(OUT)   :: KXLAG(6,NXPT)
REAL(8),                INTENT(IN)    :: PXTRAJ(NXPT)
REAL(8),    OPTIONAL,   INTENT(IN)    :: PALPHA(NXPT)
LOGICAL,    OPTIONAL,   INTENT(IN)    :: LDCOMAD

INTEGER(8)                     :: IP       ,II       ,KK
REAL(8)                        :: ZXP1     ,ZX0      ,ZXM1     ,ZXM2     ,ZXD
REAL(8)                        :: ZAL      ,ZEPS     ,ZWI      ,ZWC
LOGICAL                        :: LL_COMAD

ZEPS = ABS(EPSILON(ONE))

IF (PRESENT(LDCOMAD)) THEN
   LL_COMAD = LDCOMAD 
ELSE
   LL_COMAD = .FALSE.
END IF


IF (LL_COMAD) THEN

   ! COMPUTE WEIGHTS FOR SL HORIZONTAL INTERPOLATIONS
  IF (LPERIO) THEN      
     DO II=1,NXPT
        ZXD = PXTRAJ(II) 
        ZAL = REAL(II-1,8)-(ZXD/RDX)
        IF (ABS(ZAL).LT.ZEPS) ZAL = ZERO
        IP  = FLOOR(ZAL)

        ZXM1= REAL(II-IP-1-1,8)*RDX
        ZWI = (ZXD-ZXM1)/RDX
        ZWC = ZWI+(ZWI-HALF)*(PALPHA(II)-ONE)

        !---------------------------------------------------------------------------------------
        ! Linear interpolating between II-IP-1, II-IP
        DO KK=5,6
           KXLAG(KK,II) = NPERIO(II-IP+5-KK) 
        END DO
        PXWEI(5,II)=ZWC
        PXWEI(6,II)=ONE-ZWC
        
        !---------------------------------------------------------------------------------------
        ! cubic interpolation between II-IP-2, II-IP-1, II-IP and II-IP+1
        DO KK=1,4
           KXLAG(KK,II) = NPERIO(II-IP+2-KK) 
        END DO
        PXWEI(1,II)=((ZWC**2)-ONE)*(ZWC/6.d0)
        PXWEI(2,II)=ZWC+((ZWC**2)/TWO)*(ONE-ZWC)
        PXWEI(3,II)=ONE-(ZWC**2)+(ZWC/TWO)*((ZWC**2)-ONE)
        PXWEI(4,II)=ONE-PXWEI(1,II)-PXWEI(2,II)-PXWEI(3,II)
        !---------------------------------------------------------------------------------------
    END DO   
  ELSE
    DO II=1,NXPT
       ZXD = PXTRAJ(II) 
       ZAL = REAL(II-1,8)-(ZXD/RDX)
       IF (ABS(ZAL).LT.ZEPS) ZAL = ZERO
       IP  = FLOOR(ZAL)

       !*---------------------------------------------------------------------------------------
       ! Linear interpolation between II-IP-1, and II-IP
       DO KK=5,6
          KXLAG(KK,II) = MIN(NXPT,MAX(II-IP+5-KK,1)) 
       END DO
       IF (((II-IP).LE.NXPT).AND.((II-IP-1).GE.1)) THEN
           ZXM1=REAL(II-IP-1-1,8)*RDX
           ZWI = (ZXD-ZXM1)/RDX
           ZWC = ZWI+(ZWI-HALF)*(PALPHA(II)-ONE)
           PXWEI(5,II)=ZWC
           PXWEI(6,II)=ONE-ZWC
       ELSE IF (II-IP.LE.1) THEN
           PXWEI(5,II)=ONE
           PXWEI(6,II)=ZERO
       ELSE IF (II-IP-1.GE.NXPT) THEN
           PXWEI(5,II)=ZERO
           PXWEI(6,II)=ONE
       END IF     
      
       DO KK=1,4
          KXLAG(KK,II) = MIN(NXPT,MAX(II-IP+2-KK,1)) 
       END DO
       IF (((II-IP+1).LE.NXPT).AND.((II-IP-2).GE.1)) THEN
       ! cubic interpolation between II-IP-2, II-IP-1, II-IP and II-IP+1
           ZXM1= REAL(II-IP-1-1,8)*RDX
           ZWI = (ZXD-ZXM1)/RDX
           ZWC = ZWI+(ZWI-HALF)*(PALPHA(II)-ONE)
           PXWEI(1,II)=((ZWC**2)-ONE)*(ZWC/6.d0)
           PXWEI(2,II)=ZWC+((ZWC**2)/TWO)*(ONE-ZWC)
           PXWEI(3,II)=ONE-(ZWC**2)+(ZWC/TWO)*((ZWC**2)-ONE)
           PXWEI(4,II)=ONE-PXWEI(1,II)-PXWEI(2,II)-PXWEI(3,II)
       ELSE IF (II-IP.LE.1) THEN
       ! trajectory truncation at x=0
           PXWEI(1,II)=ZERO
           PXWEI(2,II)=ONE
           PXWEI(3,II)=ZERO
           PXWEI(4,II)=ZERO
       ELSE IF (II-IP-1.EQ.1) THEN
           ! linear interpolation between II-IP and II-IP-1 at x=0
           ZXM1=REAL(II-IP-1-1,8)*RDX
           ZWI = (ZXD-ZXM1)/RDX
           ZWC = ZWI+(ZWI-HALF)*(PALPHA(II)-ONE)
           PXWEI(1,II)=ZERO
           PXWEI(2,II)=ZWC
           PXWEI(3,II)=ONE-ZWC
           PXWEI(4,II)=ZERO
       ELSE IF (II-IP-1.GE.NXPT) THEN
         ! trajectory truncation at x=L
           PXWEI(1,II)=ZERO
           PXWEI(2,II)=ZERO
           PXWEI(3,II)=ONE
           PXWEI(4,II)=ZERO
       ELSE IF (II-IP.EQ.NXPT) THEN
           ! linear interpolation between II-IP-2 and II-IP-1 at x=L
           ZXM1=REAL(II-IP-1-1,8)*RDX
           ZWI = (ZXD-ZXM1)/RDX
           ZWC = ZWI+(ZWI-HALF)*(PALPHA(II)-ONE)
           PXWEI(1,II)=ZERO
           PXWEI(2,II)=ZWC
           PXWEI(3,II)=ONE-ZWC
           PXWEI(4,II)=ZERO
       ELSE
           WRITE(*,*) 'bug in the code in the interpolation'
           STOP
       END IF
    END DO
  END IF
 
ELSE !****** NOT COMAD ****!
   
  ! COMPUTE WEIGHTS FOR SL HORIZONTAL INTERPOLATIONS
  IF (LPERIO) THEN      
     DO II=1,NXPT
        ZXD   = PXTRAJ(II) 
        ZAL   = REAL(II-1,8)-(ZXD/RDX)
        IF (ABS(ZAL).LT.ZEPS) ZAL = ZERO
        
        IP = FLOOR(ZAL)
        DO KK=1,4
           KXLAG(KK,II) = NPERIO(II-IP+2-KK) 
        END DO
        ! cubic interpolation between II-IP-2, II-IP-1, II-IP and II-IP+1
        ZXP1=REAL(II-IP+1-1,8)*RDX
        ZX0 =REAL(II-IP+0-1,8)*RDX
        ZXM1=REAL(II-IP-1-1,8)*RDX
        ZXM2=REAL(II-IP-2-1,8)*RDX
        PXWEI(1,II)=((ZXD-ZX0)/(ZXP1-ZX0))*((ZXD-ZXM1)/(ZXP1-ZXM1))*((ZXD-ZXM2)/(ZXP1-ZXM2))
        PXWEI(2,II)=((ZXD-ZXP1)/(ZX0-ZXP1))*((ZXD-ZXM1)/(ZX0-ZXM1))*((ZXD-ZXM2)/(ZX0-ZXM2))
        PXWEI(3,II)=((ZXD-ZXP1)/(ZXM1-ZXP1))*((ZXD-ZX0)/(ZXM1-ZX0))*((ZXD-ZXM2)/(ZXM1-ZXM2))
        PXWEI(4,II)=((ZXD-ZXP1)/(ZXM2-ZXP1))*((ZXD-ZX0)/(ZXM2-ZX0))*((ZXD-ZXM1)/(ZXM2-ZXM1))

        !---------------------------------------------------------------------------------------
        DO KK=5,6
           KXLAG(KK,II) = NPERIO(II-IP+5-KK) 
        END DO
        
        ! Linear interpolating weights between II-IP-1, II-IP
        PXWEI(5,II)=(ZXD-ZXM1)/(ZX0-ZXM1)
        PXWEI(6,II)=(ZXD-ZX0)/(ZXM1-ZX0)
        !---------------------------------------------------------------------------------------
    END DO   
  ELSE
    DO II=1,NXPT
       ZXD   = PXTRAJ(II) 
       ZAL   = REAL(II-1,8)-(ZXD/RDX)
       IF (ABS(ZAL).LT.ZEPS) ZAL = ZERO
       IP = FLOOR(ZAL)
       
       DO KK=1,4
          KXLAG(KK,II) = MIN(NXPT,MAX(II-IP+2-KK,1)) 
       END DO
       IF (((II-IP+1).LE.NXPT).AND.((II-IP-2).GE.1)) THEN
       ! cubic interpolation between II-IP-2, II-IP-1, II-IP and II-IP+1
           ZXP1=REAL(II-IP+1-1,8)*RDX
           ZX0 =REAL(II-IP+0-1,8)*RDX
           ZXM1=REAL(II-IP-1-1,8)*RDX
           ZXM2=REAL(II-IP-2-1,8)*RDX
           PXWEI(1,II)=((ZXD-ZX0)/(ZXP1-ZX0))*((ZXD-ZXM1)/(ZXP1-ZXM1))*((ZXD-ZXM2)/(ZXP1-ZXM2))
           PXWEI(2,II)=((ZXD-ZXP1)/(ZX0-ZXP1))*((ZXD-ZXM1)/(ZX0-ZXM1))*((ZXD-ZXM2)/(ZX0-ZXM2))
           PXWEI(3,II)=((ZXD-ZXP1)/(ZXM1-ZXP1))*((ZXD-ZX0)/(ZXM1-ZX0))*((ZXD-ZXM2)/(ZXM1-ZXM2))
           PXWEI(4,II)=((ZXD-ZXP1)/(ZXM2-ZXP1))*((ZXD-ZX0)/(ZXM2-ZX0))*((ZXD-ZXM1)/(ZXM2-ZXM1))
       ELSE IF (II-IP.LE.1) THEN
       ! trajectory truncation at x=0
           PXWEI(1,II)=ZERO
           PXWEI(2,II)=ONE
           PXWEI(3,II)=ZERO
           PXWEI(4,II)=ZERO
       ELSE IF (II-IP-1.EQ.1) THEN
           ! linear interpolation between II-IP and II-IP+1 at x=0
           ZXM1=REAL(II-IP-1-1,8)*RDX
           ZX0 =REAL(II-IP+0-1,8)*RDX
           PXWEI(1,II)=ZERO
           PXWEI(2,II)=(ZXD-ZXM1)/(ZX0-ZXM1)
           PXWEI(3,II)=(ZXD-ZX0)/(ZXM1-ZX0)
           PXWEI(4,II)=ZERO
       ELSE IF (II-IP-1.GE.NXPT) THEN
         ! trajectory truncation at x=L
           PXWEI(1,II)=ZERO
           PXWEI(2,II)=ZERO
           PXWEI(3,II)=ONE
           PXWEI(4,II)=ZERO
       ELSE IF (II-IP.EQ.NXPT) THEN
           ! linear interpolation between II-IP-2 and II-IP-1 at x=L
           ZX0=REAL(II-IP+0-1,8)*RDX
           ZXM1=REAL(II-IP-1-1,8)*RDX
           PXWEI(1,II)=ZERO
           PXWEI(2,II)=(ZXD-ZXM1)/(ZX0-ZXM1)
           PXWEI(3,II)=(ZXD-ZX0)/(ZXM1-ZX0)
           PXWEI(4,II)=ZERO
       ELSE
           WRITE(*,*) 'bug in the code in the interpolation'
           STOP
       END IF

       DO KK=5,6
          KXLAG(KK,II) = MIN(NXPT,MAX(II-IP+5-KK,1)) 
       END DO
       IF (((II-IP).LE.NXPT).AND.((II-IP-1).GE.1)) THEN
       ! Linear interpolation between II-IP-1, and II-IP 
           ZX0 =REAL(II-IP+0-1,8)*RDX
           ZXM1=REAL(II-IP-1-1,8)*RDX
           PXWEI(5,II)=(ZXD-ZXM1)/(ZX0-ZXM1)
           PXWEI(6,II)=(ZXD-ZX0)/(ZXM1-ZX0)
       ELSE IF (II-IP.LE.1) THEN
           PXWEI(5,II)=ONE
           PXWEI(6,II)=ZERO
       ELSE IF (II-IP-1.GE.NXPT) THEN
           PXWEI(5,II)=ZERO
           PXWEI(6,II)=ONE
       END IF     
    END DO
  END IF
 
END IF

END SUBROUTINE ELASCAW_SURF
!*************************************************************************************!
!*************************************************************************************!
!*************************************************************************************!
!*************************************************************************************!
!*************************************************************************************!
!*************************************************************************************!
!*************************************************************************************!
SUBROUTINE COMPUTE_DEFORMATIONAL(PDXDOT_DX,PDXDOT_DZ,PDZDOT_DX,PDZDOT_DZ,&
                               & PXDOT,PZDOT,PETA,KTOP)
    
IMPLICIT NONE
  
INTEGER,           INTENT(IN)    :: KTOP
REAL(8),           INTENT(IN)    :: PXDOT(KTOP:NLEV,NXPT),      PZDOT(KTOP:NLEV,NXPT)
REAL(8),           INTENT(IN)    :: PETA(KTOP:NLEV)
REAL(8),           INTENT(OUT)   :: PDXDOT_DX(KTOP:NLEV,NXPT),  PDZDOT_DZ(KTOP:NLEV,NXPT)
REAL(8),           INTENT(OUT)   :: PDXDOT_DZ(KTOP:NLEV,NXPT),  PDZDOT_DX(KTOP:NLEV,NXPT)

INTEGER(8)                       :: II         ,JJ
REAL(8)                          :: ZETA(0:NLEV+1),        ZDETA(0:NLEV)
REAL(8)                          :: ZXDOT(0:NLEV+1,NXPT),  ZZDOT(0:NLEV+1,NXPT),   ZVP(KTOP:NLEV,NXPT)
REAL(8)                          :: ZDISCR,                ZEPS,                   ZLIPCHITZ,   ZTAU

ZEPS = ABS(EPSILON(ONE))
ZTAU = ONE-ZEPS

ZETA(0)      = ZERO
ZETA(NLEV+1) = ONE

DO JJ=KTOP,NLEV
   ZETA(JJ)  = PETA(JJ)
END DO   

DO JJ=0,NLEV
   ZDETA(JJ)    = ZETA(JJ+1)-ZETA(JJ)
END DO

ZXDOT(0,:)      = PXDOT(KTOP,:)
ZXDOT(NLEV+1,:) = PXDOT(NLEV,:)

ZZDOT(0,:)      = ZERO
ZZDOT(NLEV+1,:) = ZERO

DO JJ=KTOP,NLEV
   ZXDOT(JJ,:)  = PXDOT(JJ,:)
   ZZDOT(JJ,:)  = PZDOT(JJ,:)
END DO

!----------------------------------
! COMPUTE DEFORMATIONAL TERM
!----------------------------------

DO II=1,NXPT
   DO JJ=KTOP,NLEV
      PDXDOT_DX(JJ,II) = RFDC(II)*( PXDOT(JJ,NPERIO(II+1))                 &
                       & - PXDOT(JJ,NPERIO(II-1)) )/RDX
      PDXDOT_DZ(JJ,II) = 0.5d0*((ZXDOT(JJ,II)-ZXDOT(JJ-1,II))/ZDETA(JJ-1)) &
                       & + 0.5d0*((ZXDOT(JJ+1,II)-ZXDOT(JJ,II))/ZDETA(JJ))
      PDZDOT_DX(JJ,II) = RFDC(II)*( PZDOT(JJ,NPERIO(II+1))                 &
                       & - PZDOT(JJ,NPERIO(II-1)) )/RDX
      PDZDOT_DZ(JJ,II) = 0.5d0*((ZZDOT(JJ,II)-ZZDOT(JJ-1,II))/ZDETA(JJ-1)) &
                       & + 0.5d0*((ZZDOT(JJ+1,II)-ZZDOT(JJ,II))/ZDETA(JJ))
   ENDDO
ENDDO   

!----------------------------------
! COMPUTE LIPCHITZ NUMBERS
!----------------------------------

DO II=1,NXPT
   DO JJ=KTOP,NLEV
      ZDISCR        = ((PDXDOT_DX(JJ,II)-PDZDOT_DZ(JJ,II))**2)  &
                    & + 4.d0*PDXDOT_DZ(JJ,II)*PDZDOT_DX(JJ,II)
      IF (ZDISCR .GE. ZEPS) THEN
         ZVP(JJ,II) = 0.25d0*RDT*MAX(ABS((PDXDOT_DX(JJ,II)      &
                    & +PDZDOT_DZ(JJ,II))+SQRT(ZDISCR)), &
                    & ABS((PDXDOT_DX(JJ,II)+PDZDOT_DZ(JJ,II))   &
                    & -SQRT(ZDISCR)))
      ELSE
         ZVP(JJ,II) = 0.25d0*RDT*SQRT(((PDXDOT_DX(JJ,II)        &
                    & +PDZDOT_DZ(JJ,II))**2)+ABS(ZDISCR))
      END IF   
   ENDDO
ENDDO 

!****************************************************************************!
IF (LPRINTLEV) THEN                                                          !
   IF (MAXVAL(ZVP) .GE. ZTAU) THEN                                           !
      WRITE(*,*) 'VIOLATION OF THE LIPCHITZ CRITERIUM AT POINTS : '          !       
      !DO I=1,NXPT                                                           !
      !   DO J=1,NLEV                                                        !
      !      ZLIPCHITZ=ZVP(J,I)                                              !
      !      IF (ZLIPCHITZ .GE. ZTAU) THEN                                   !
      !         WRITE (*,*) '(I,J) = ','(',I,',',J,')',' .... ',ZLIPCHITZ    !
      !      END IF                                                          !
      !   END DO                                                             !
      !END DO                                                                !     
   END IF                                                                    !             
END IF                                                                       !
!****************************************************************************!

END SUBROUTINE COMPUTE_DEFORMATIONAL
!*************************************************************************************************!
!*************************************************************************************************!
!*************************************************************************************************!
!*************************************************************************************************!
SUBROUTINE COMPUTE_COMAD_STRETCHING(PALPHA_DX,PALPHA_DZ,PXDOT_F,PZDOT_F,PXDOT,PZDOT,PDXDOT_DX, &
                                  & PDXDOT_DX_F,KXLAG,PXWEI,KZLAG,PZWEI,PETA,KTOP,CD_SLTRAJ)
    
IMPLICIT NONE
  
INTEGER,           INTENT(IN)    :: KTOP
REAL(8),           INTENT(IN)    :: PETA(KTOP:NLEV)
REAL(8),           INTENT(IN)    :: PXDOT(KTOP:NLEV,NXPT),       PZDOT(KTOP:NLEV,NXPT)
REAL(8),           INTENT(IN)    :: PXDOT_F(KTOP:NLEV,NXPT),     PZDOT_F(KTOP:NLEV,NXPT)
REAL(8),           INTENT(IN)    :: PDXDOT_DX(KTOP:NLEV,NXPT),   PDXDOT_DX_F(KTOP:NLEV,NXPT)
REAL(8),           INTENT(IN)    :: PXWEI(6,KTOP:NLEV,NXPT),     PZWEI(6,KTOP:NLEV,NXPT)
INTEGER(8),        INTENT(IN)    :: KXLAG(6,KTOP:NLEV,NXPT),     KZLAG(6,KTOP:NLEV,NXPT)
REAL(8),           INTENT(OUT)   :: PALPHA_DX(KTOP:NLEV,NXPT),   PALPHA_DZ(KTOP:NLEV,NXPT)
CHARACTER(LEN=1),  INTENT(IN)    :: CD_SLTRAJ

INTEGER(8)                       :: II,  JJ,  KK,  LL
INTEGER                          :: ISTA,     IEND
REAL(8)                          :: ZDIVZ,            ZEPS,       ZDIVX,      ZSDT       
REAL(8)                          :: ZETA(0:NLEV+1),              ZDETA(0:NLEV)
REAL(8)                          :: ZZDOT(0:NLEV+1,NXPT),        ZZDOT_F(0:NLEV+1,NXPT)
REAL(8)                          :: ZDXDOT_DX_F(KTOP:NLEV,NXPT), ZDZDOT_DZ_F(KTOP:NLEV,NXPT)
REAL(8)                          :: ZDXDOT_DX_O(KTOP:NLEV,NXPT), ZDZDOT_DZ_O(KTOP:NLEV,NXPT)
REAL(8)                          :: ZDXDOT_DX(KTOP:NLEV,NXPT),   ZDZDOT_DZ(KTOP:NLEV,NXPT)

ZEPS = ABS(1000*EPSILON(ONE))
ZSDT = (TWO/RDT)*(ONE-ZEPS)

ZETA(0)      = ZERO
ZETA(NLEV+1) = ONE

DO JJ=KTOP,NLEV
   ZETA(JJ)  = PETA(JJ)
END DO   

DO JJ=0,NLEV
   ZDETA(JJ)    = ZETA(JJ+1)-ZETA(JJ)
END DO

ZZDOT_F(0,:)      = ZERO
ZZDOT_F(NLEV+1,:) = ZERO

DO JJ=KTOP,NLEV
   ZZDOT_F(JJ,:)  = PZDOT_F(JJ,:)
   ZZDOT(JJ,:)    = PZDOT(JJ,:)
END DO

ZZDOT(0,:)        = ZERO
ZZDOT(NLEV+1,:)   = ZERO



IF (NCOMAD_OPT == 0) THEN
!----------------------------------
! COMPUTE STRETCHING COEFFICIENTS   
  DO II=1,NXPT
     DO JJ=KTOP,NLEV
        PALPHA_DX(JJ,II)       = MAX(0.0001d0,MIN((ONE+RDT*PDXDOT_DX(JJ,II)),ONE))
        PALPHA_DZ(JJ,II)       = MAX(0.0001d0,MIN((ONE-RDT*PDXDOT_DX(JJ,II)),ONE))
     ENDDO     
  ENDDO     
ELSE   
  !----------------------------------
  ! COMPUTE DEFORMATIONAL TERMS
  IF (NCOMAD_OPT == 1) THEN
      ZDXDOT_DX(:,:)           =  PDXDOT_DX(:,:)
      ZDXDOT_DX_F(:,:)         =  PDXDOT_DX_F(:,:)
      ZDZDOT_DZ(:,:)           = -PDXDOT_DX(:,:)
      ZDZDOT_DZ_F(:,:)         = -PDXDOT_DX_F(:,:)
  ELSE IF (NCOMAD_OPT == 2) THEN  
      DO II=1,NXPT
         DO JJ=KTOP,NLEV
            ZDXDOT_DX(JJ,II)   = PDXDOT_DX(JJ,II)
            ZDZDOT_DZ(JJ,II)   = 0.5d0*((ZZDOT(JJ,II)-ZZDOT(JJ-1,II))/ZDETA(JJ-1))     &
                               & + 0.5d0*((ZZDOT(JJ+1,II)-ZZDOT(JJ,II))/ZDETA(JJ))
            ZDXDOT_DX_F(JJ,II) = PDXDOT_DX_F(JJ,II)
            ZDZDOT_DZ_F(JJ,II) = 0.5d0*((ZZDOT_F(JJ,II)-ZZDOT_F(JJ-1,II))/ZDETA(JJ-1)) &
                               & + 0.5d0*((ZZDOT_F(JJ+1,II)-ZZDOT_F(JJ,II))/ZDETA(JJ))
         ENDDO
      ENDDO
  ELSE
      DO II=1,NXPT
         DO JJ=KTOP,NLEV
            ZDXDOT_DX(JJ,II)   = RFDC(II)*( PXDOT(JJ,NPERIO(II+1))                     &
                               & - PXDOT(JJ,NPERIO(II-1)) )/RDX
            ZDZDOT_DZ(JJ,II)   = 0.5d0*((ZZDOT(JJ,II)-ZZDOT(JJ-1,II))/ZDETA(JJ-1))     &
                               & + 0.5d0*((ZZDOT(JJ+1,II)-ZZDOT(JJ,II))/ZDETA(JJ))
            ZDXDOT_DX_F(JJ,II) = RFDC(II)*( PXDOT_F(JJ,NPERIO(II+1))                   &
                               & - PXDOT_F(JJ,NPERIO(II-1)) )/RDX
            ZDZDOT_DZ_F(JJ,II) = 0.5d0*((ZZDOT_F(JJ,II)-ZZDOT_F(JJ-1,II))/ZDETA(JJ-1)) &
                               & + 0.5d0*((ZZDOT_F(JJ+1,II)-ZZDOT_F(JJ,II))/ZDETA(JJ))
         ENDDO
      ENDDO
  END IF   
  !----------------------------------
  ! INTERPOLATED DEFORMATIONAL
  IF (C_SLTRAJ == "C") THEN
     ISTA = 1
     IEND = 4
  ELSE IF (C_SLTRAJ == "L") THEN
     ISTA = 5
     IEND = 6
  END IF

  DO II=1,NXPT
     DO JJ=KTOP,NLEV
        ZDXDOT_DX_O(JJ,II) = ZERO
        ZDZDOT_DZ_O(JJ,II) = ZERO
        DO KK=ISTA,IEND
           DO LL=ISTA,IEND
              ZDXDOT_DX_O(JJ,II) = ZDXDOT_DX_O(JJ,II)                          &
                                 & + PZWEI(LL,JJ,II)*PXWEI(KK,JJ,II)*(         &
                                 &   ZDXDOT_DX(KZLAG(LL,JJ,II),KXLAG(KK,JJ,II)))
              ZDZDOT_DZ_O(JJ,II) = ZDZDOT_DZ_O(JJ,II)                          &
                                 & + PZWEI(LL,JJ,II)*PXWEI(KK,JJ,II)*(         &
                                 &   ZDZDOT_DZ(KZLAG(LL,JJ,II),KXLAG(KK,JJ,II)))
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  !----------------------------------
  ! COMPUTE STRETCHING COEFFICIENTS   
  DO II=1,NXPT
     DO JJ=KTOP,NLEV
        ZDIVX = ZDXDOT_DX_O(JJ,II)
        IF (ZDIVX .LE. ZSDT) THEN
           PALPHA_DX(JJ,II) = MAX(0.0001d0,MIN((ONE+(RDT/TWO)*ZDXDOT_DX_F(JJ,II))  &
                            & /(ONE-(RDT/TWO)*ZDXDOT_DX_O(JJ,II)),ONE))
        ELSE
           PALPHA_DX(JJ,II) = ONE
        END IF    
        ZDIVZ = ZDZDOT_DZ_O(JJ,II)
        IF (ZDIVZ .LE. ZSDT) THEN
           PALPHA_DZ(JJ,II) = MAX(0.0001d0,MIN((ONE+(RDT/TWO)*ZDZDOT_DZ_F(JJ,II)) &
                            & /(ONE-(RDT/TWO)*ZDZDOT_DZ_O(JJ,II)),ONE))
        ELSE
           PALPHA_DZ(JJ,II) = ONE
        END IF
     ENDDO     
  ENDDO       
END IF

END SUBROUTINE COMPUTE_COMAD_STRETCHING
!*************************************************************************************************!
!*************************************************************************************************!
!*************************************************************************************************!
!*************************************************************************************************!
!*************************************************************************************************!
!*************************************************************************************************!
!*************************************************************************************************!
SUBROUTINE COMPUTE_COMAD_SURF(PALPHA_DX,PXDOT_F,PXDOT,PDXDOT_DX_F, &
                            & PDXDOT_DX,KXLAGS,PXWEIS,PDELB,CD_SLTRAJ)
    
IMPLICIT NONE
  
REAL(8),           INTENT(IN)    :: PDELB(NLEV)
REAL(8),           INTENT(IN)    :: PXDOT(NLEV,NXPT),     PXDOT_F(NLEV,NXPT)
REAL(8),           INTENT(IN)    :: PDXDOT_DX(NLEV,NXPT), PDXDOT_DX_F(NLEV,NXPT)
REAL(8),           INTENT(IN)    :: PXWEIS(6,NXPT)
INTEGER(8),        INTENT(IN)    :: KXLAGS(6,NXPT)
REAL(8),           INTENT(OUT)   :: PALPHA_DX(NXPT)
CHARACTER(LEN=1),  INTENT(IN)    :: CD_SLTRAJ

INTEGER(8)                       :: I,     J,    K
INTEGER                          :: ISTA,  IEND
REAL(8)                          :: ZEPS,  ZDIVX,  ZSDT    
REAL(8)                          :: ZXDOTS(NXPT),ZXDOTS_F(NXPT)
REAL(8)                          :: ZDXDOT_DX_F(NXPT),ZDXDOT_DX_O(NXPT),ZDXDOT_DX(NXPT)

ZEPS = ABS(1000*EPSILON(ONE))
ZSDT = (TWO/RDT)*(ONE-ZEPS)

IF (NCOMAD_OPT == 0) THEN
!----------------------------------
! COMPUTE STRETCHING COEFFICIENTS   
   DO I=1,NXPT
      ZDXDOT_DX(I)    = ZERO
      DO J=1,NLEV
         ZDXDOT_DX(I) = ZDXDOT_DX(I)+PDELB(J)*PDXDOT_DX(J,I)
      END DO   
      PALPHA_DX(I)    = MAX(0.0001d0,MIN((ONE+RDT*ZDXDOT_DX(I)),ONE)) 
   ENDDO     
ELSE   
  !----------------------------------
  ! COMPUTE DEFORMATIONAL TERMS
  IF ((NCOMAD_OPT == 1).OR.(NCOMAD_OPT == 2)) THEN
     DO I=1,NXPT
        ZDXDOT_DX(I)      = ZERO
        ZDXDOT_DX_F(I)    = ZERO
        DO J=1,NLEV
           ZDXDOT_DX(I)   = ZDXDOT_DX(I)+PDELB(J)*PDXDOT_DX(J,I)
           ZDXDOT_DX_F(I) = ZDXDOT_DX_F(I)+PDELB(J)*PDXDOT_DX_F(J,I)
        END DO   
     END DO   
  ELSE
     DO I=1,NXPT
        ZXDOTS(I)         = ZERO
        ZXDOTS_F(I)       = ZERO
        DO J=1,NLEV
           ZXDOTS(I)      = ZXDOTS(I)+PDELB(J)*PXDOT(J,I)
           ZXDOTS_F(I)    = ZXDOTS_F(I)+PDELB(J)*PXDOT_F(J,I)
        END DO   
     END DO
     DO I=1,NXPT
        ZDXDOT_DX(I)      = RFDC(I)*( ZXDOTS(NPERIO(I+1))   &
                          & - ZXDOTS(NPERIO(I-1)) )/RDX
        ZDXDOT_DX_F(I)    = RFDC(I)*( ZXDOTS_F(NPERIO(I+1)) &
                          & - ZXDOTS_F(NPERIO(I-1)) )/RDX
     ENDDO
  END IF   
  !----------------------------------
  ! INTERPOLATED DEFORMATIONAL
  IF (C_SLTRAJ == "C") THEN
     ISTA = 1
     IEND = 4
  ELSE IF (C_SLTRAJ == "L") THEN
     ISTA = 5
     IEND = 6
  END IF
  DO I=1,NXPT
     ZDXDOT_DX_O(I)    = ZERO
     DO K=ISTA,IEND
        ZDXDOT_DX_O(I) = ZDXDOT_DX_O(I)           &
              & + PXWEIS(K,I)*ZDXDOT_DX(KXLAGS(K,I))
    ENDDO  
  ENDDO
  !----------------------------------
  ! COMPUTE STRETCHING COEFFICIENTS   
  DO I=1,NXPT
     ZDIVX = ZDXDOT_DX_O(I)
     IF (ZDIVX .LE. ZSDT) THEN
        PALPHA_DX(I) = MAX(0.0001d0,MIN((ONE+(RDT/TWO)*ZDXDOT_DX_F(I))  &
                     & /(ONE-(RDT/TWO)*ZDXDOT_DX_O(I)),ONE))
     ELSE
        PALPHA_DX(I) = ONE
     END IF    
  ENDDO       
END IF

END SUBROUTINE COMPUTE_COMAD_SURF
!*************************************************************************************************!
!*************************************************************************************************!
!*************************************************************************************************!
!*************************************************************************************************!
!*************************************************************************************************!
SUBROUTINE LARCINB(RQ_DEP,RQ_ARR,RM_DEP,RM_ARR,RPIS_DEP,RPIS_ARR, &
                 & PXWEI,KXLAG,PZWEI,KZLAG,PXWEIS,KXLAGS,CD_SLINTP)
  
IMPLICIT NONE

REAL(8),          INTENT(OUT) :: RQ_DEP(NLEV,NXPT),RM_DEP(NLEV,NXPT),RPIS_DEP(NXPT)
REAL(8),          INTENT(IN)  :: RQ_ARR(NLEV,NXPT),RM_ARR(NLEV,NXPT),RPIS_ARR(NXPT)
REAL(8),          INTENT(IN)  :: PXWEI(6,NLEV,NXPT),PZWEI(6,NLEV,NXPT),PXWEIS(6,NXPT)
INTEGER(8),       INTENT(IN)  :: KXLAG(6,NLEV,NXPT),KZLAG(6,NLEV,NXPT),KXLAGS(6,NXPT)
CHARACTER(LEN=1), INTENT(IN)  :: CD_SLINTP

INTEGER(8)                    :: KK ,LL ,II   ,JJ
INTEGER                       :: ISTA  ,IEND


IF (CD_SLINTP == "C") THEN
   ISTA = 1
   IEND = 4
ELSE IF (CD_SLINTP == "L") THEN
   ISTA = 5
   IEND = 6
ELSE
   WRITE(*,*) "UNKNOWN TRAJECTORY INTERPOLATION ORDER : ",CD_SLINTP
   STOP
END IF

DO II=1,NXPT
   DO JJ=1,NLEV
      RQ_DEP(JJ,II) = ZERO
      RM_DEP(JJ,II) = ZERO
      DO KK=ISTA,IEND
         DO LL=ISTA,IEND
            RQ_DEP(JJ,II) = RQ_DEP(JJ,II)                 &
              & + PZWEI(LL,JJ,II)*PXWEI(KK,JJ,II)*(       &
              &   RQ_ARR(KZLAG(LL,JJ,II),KXLAG(KK,JJ,II)) )
            RM_DEP(JJ,II) = RM_DEP(JJ,II)                 &
              & + PZWEI(LL,JJ,II)*PXWEI(KK,JJ,II)*(       &
              &   RQ_ARR(KZLAG(LL,JJ,II),KXLAG(KK,JJ,II)) )
       ENDDO
     ENDDO
   ENDDO
   RPIS_DEP(II) = ZERO
   DO KK=ISTA,IEND
      RPIS_DEP(II) = RPIS_DEP(II)                 &
          & + PXWEIS(KK,II)*RPIS_ARR(KXLAGS(KK,II))
   ENDDO  
 ENDDO   

END SUBROUTINE LARCINB
!*************************************************************************************************!
!*************************************************************************************************!
!*************************************************************************************************!
!*************************************************************************************************!
!*************************************************************************************************!
!*************************************************************************************************!
!*************************************************************************************************!
!*************************************************************************************************!
SUBROUTINE ILMC_FILTER(PDEP,PARR,PM,KXLAG,KYLAG)
  
IMPLICIT NONE

REAL(8),          INTENT(INOUT) :: PDEP(NLEV,NXPT)
REAL(8),          INTENT(IN)    :: PARR(NLEV,NXPT),PDELPI(NLEV,NXPT)
INTEGER(8),       INTENT(IN)    :: KXLAG(6,NLEV,NXPT),KZLAG(6,NLEV,NXPT)

INTEGER(8)                      :: I  ,J
REAL(8)                         :: ZM(24)
REAL(8)                         :: ZMAX_LOC(NLEV,NXPT), ZMIN_LOC(NLEV,NXPT)


DO I=1,NLEV
   DO J=1,NXPT
      ZMAX_LOC(I,J) =  MAX(                   &
           & PARR(KXLAG(2,I,J),KYLAG(2,I,J)), &
           & PARR(KXLAG(2,I,J),KYLAG(3,I,J)), &
           & PARR(KXLAG(3,I,J),KYLAG(2,I,J)), &
           & PARR(KXLAG(3,I,J),KYLAG(3,I,J))  )
      ZMIN_LOC(I,J) =  MIN(                   &
           & PARR(KXLAG(2,I,J),KYLAG(2,I,J)), &
           & PARR(KXLAG(2,I,J),KYLAG(3,I,J)), &
           & PARR(KXLAG(3,I,J),KYLAG(2,I,J)), &
           & PARR(KXLAG(3,I,J),KYLAG(3,I,J))  )
   ENDDO
ENDDO


!_____._____._____._____._____._____._____.
!     !     !     !     !     !     !     !
! 42  !  43 !  44 !  45 !  46 !  47 !  48 ! 
!_____!_____!_____!_____!_____!_____!_____!
!     !     !     !     !     !     !     !
! 40  !  20 !  21 !  22 !  23 !  24 !  41 !
!_____!_____!_____!_____!_____!_____!_____!
!     !     !     !     !     !     !     !
! 38  !  18 !  6  !  7  !  8  !  19 !  39 !     
!_____!_____!_____!_____!_____!_____!_____!
!     !     !     !     !     !     !     !
! 36  !  16 !  4  !  X  !  5  !  17 !  37 !
!_____!_____!_____!_____!_____!_____!_____!
!     !     !     !     !     !     !     !
! 34  !  14 !  1  !  2  !  3  !  15 !  35 !
!_____!_____!_____!_____!_____!_____!_____!
!     !     !     !     !     !     !     !
! 32  !  9  !  10 !  11 !  12 !  13 !  33 !
!_____!_____!_____!_____!_____!_____!_____!
!     !     !     !     !     !     !     !
! 25  !  26 !  27 !  28 !  29 !  30 !  31 !
!_____!_____!_____!_____!_____!_____!_____!


DO I=1,NXPT
   DO J=2,NLEV-1
      !----------------------------!
      !  DEALING WITH OVERSHOOTS   !
      !----------------------------!
      IF (PDEP(I,J) > ZMAX_LOC(I,J)) THEN
         ZM(0) = (PDEP(J,I)-ZMAX_LOC(J,I))*PDELPI(J,I)
         ZM(1) = MAX(ZMAX_LOC(J-1,NPERIO(I-1))-PDEP(J-1,NPERIO(I-1)),ZERO)*PDELPI(J-1,NPERIO(I-1))
         ZM(2) = MAX(ZMAX_LOC(J-1,I)-PDEP(J-1,I),ZERO)*PDELPI(J-1,I)
         ZM(3) = MAX(ZMAX_LOC(J-1,NPERIO(I+1))-PDEP(J-1,NPERIO(I+1)),ZERO)*PDELPI(J-1,NPERIO(I+1))
         ZM(4) = MAX(ZMAX_LOC(J,NPERIO(I-1))-PDEP(J,NPERIO(I-1)),ZERO)*PDELPI(J,NPERIO(I-1))
         ZM(5) = MAX(ZMAX_LOC(J,NPERIO(I+1))-PDEP(J,NPERIO(I+1)),ZERO)*PDELPI(J,NPERIO(I+1))
         ZM(6) = MAX(ZMAX_LOC(J+1,NPERIO(I-1))-PDEP(J+1,NPERIO(I-1)),ZERO)*PDELPI(J+1,NPERIO(I-1))
         ZM(7) = MAX(ZMAX_LOC(J+1,I)-PDEP(J+1,I),ZERO)*PDELPI(J+1,I)
         ZM(8) = MAX(ZMAX_LOC(J+1,NPERIO(I+1))-PDEP(J+1,NPERIO(I+1)),ZERO)*PDELPI(J+1,NPERIO(I+1))
         ZM_A  = ZERO
         DO K=1,8
            ZM_A = ZM_A + ZM(K)
         END DO
         IF (ZM_A > ZM(0)) THEN
            PDEP(J,I)     = ZMAX_LOC(J,I)
            PDEP(J-1,NPERIO(I-1)) = PDEP(J-1,NPERIO(I-1)) + (ZM(0)/PDELPI(J-1,NPERIO(I-1)))*(ZM(1)/ZM_A)
            PDEP(J-1,I)           = PDEP(J-1,I)           + (ZM(0)/PDELPI(J-1,I))*(ZM(2)/ZM_A)
            PDEP(J-1,NPERIO(I+1)  = PDEP(J-1,NPERIO(I+1)) + (ZM(0)/PDELPI(J-1,NPERIO(I+1)))*(ZM(3)/ZM_A)
            PDEP(J,NPERIO(I-1))   = PDEP(J,NPERIO(I-1))   + (ZM(0)/PDELPI(J,NPERIO(I-1)))*(ZM(4)/ZM_A)
            PDEP(J,NPERIO(I+1))   = PDEP(J,NPERIO(I+1))   + (ZM(0)/PDELPI(J,NPERIO(I+1)))*(ZM(5)/ZM_A)
            PDEP(J+1,NPERIO(I-1)) = PDEP(J+1,NPERIO(I-1)) + (ZM(0)/PDELPI(J+1,NPERIO(I-1)))*(ZM(6)/ZM_A)
            PDEP(J+1,I)           = PDEP(J+1,I)           + (ZM(0)/PDELPI(J+1,I))*(ZM(7)/ZM_A)
            PDEP(J+1,NPERIO(I+1)) = PDEP(J+1,NPERIO(I+1)) + (ZM(0)/PDELPI(J+1,NPERIO(I+1)))*(ZM(8)/ZM_A)
         ELSE IF ((ZM_A .LE. ZM(0)).AND.(J.GE.3).AND.(J.LE.NLEV-2)) THEN
            ZM(9)  = MAX(ZMAX_LOC(J-2,NPERIO(I-2))-PDEP(J-2,NPERIO(I-2)),ZERO)*PDELPI(J-2,NPERIO(I-2))
            ZM(10) = MAX(ZMAX_LOC(J-2,NPERIO(I-1))-PDEP(J-2,NPERIO(I-1)),ZERO)*PDELPI(J-2,NPERIO(I-1))
            ZM(11) = MAX(ZMAX_LOC(J-2,I)-PDEP(J-2,I),ZERO)*PDELPI(J-2,I)
            ZM(12) = MAX(ZMAX_LOC(J-2,NPERIO(I+1))-PDEP(J-2,NPERIO(I+1)),ZERO)*PDELPI(J-2,NPERIO(I+1))
            ZM(13) = MAX(ZMAX_LOC(J-2,NPERIO(I+2))-PDEP(J-2,NPERIO(I+2)),ZERO)*PDELPI(J-2,NPERIO(I+2))
            ZM(14) = MAX(ZMAX_LOC(J-1,NPERIO(I-2))-PDEP(J-1,NPERIO(I-2)),ZERO)*PDELPI(J-1,NPERIO(I-2))
            ZM(15) = MAX(ZMAX_LOC(J-1,NPERIO(I+2))-PDEP(J-1,NPERIO(I+2)),ZERO)*PDELPI(J-1,NPERIO(I+2))
            ZM(16) = MAX(ZMAX_LOC(J,NPERIO(I-2))-PDEP(J,NPERIO(I-2)),ZERO)*PDELPI(J,NPERIO(I-2))
            ZM(17) = MAX(ZMAX_LOC(J,NPERIO(I+2))-PDEP(J,NPERIO(I+2)),ZERO)*PDELPI(J,NPERIO(I+2))
            ZM(18) = MAX(ZMAX_LOC(J+1,NPERIO(I-2))-PDEP(J+1,NPERIO(I-2)),ZERO)*PDELPI(J+1,NPERIO(I-2))
            ZM(19) = MAX(ZMAX_LOC(J+1,NPERIO(I+2))-PDEP(J+1,NPERIO(I+2)),ZERO)*PDELPI(J+1,NPERIO(I+2))
            ZM(20) = MAX(ZMAX_LOC(J+2,NPERIO(I-2))-PDEP(J+2,NPERIO(I-2)),ZERO)*PDELPI(J+2,NPERIO(I-2))
            ZM(21) = MAX(ZMAX_LOC(J+2,NPERIO(I-1))-PDEP(J+2,NPERIO(I-1)),ZERO)*PDELPI(J+2,NPERIO(I-1))
            ZM(22) = MAX(ZMAX_LOC(J+2,I)-PDEP(J+2,I),ZERO)*PDELPI(J+2,I)
            ZM(23) = MAX(ZMAX_LOC(J+2,NPERIO(I+1))-PDEP(J+2,NPERIO(I+1)),ZERO)*PDELPI(J+2,NPERIO(I+1))
            ZM(24) = MAX(ZMAX_LOC(J+2,NPERIO(I+2))-PDEP(J+2,NPERIO(I+2)),ZERO)*PDELPI(J+2,NPERIO(I+2))
            DO K=9,24
               ZM_A = ZM_A + ZM(K)
            END DO
            IF (ZM_A > ZM(0)) THEN
               PDEP(I,J)             = ZMAX_LOC(I,J)
               PDEP(J-1,NPERIO(I-1)) = PDEP(J-1,NPERIO(I-1)) + (ZM(0)/PDELPI(J-1,NPERIO(I-1)))*(ZM(1)/ZM_A)
               PDEP(J-1,I)           = PDEP(J-1,I)           + (ZM(0)/PDELPI(J-1,I))*(ZM(2)/ZM_A)
               PDEP(J-1,NPERIO(I+1)  = PDEP(J-1,NPERIO(I+1)) + (ZM(0)/PDELPI(J-1,NPERIO(I+1)))*(ZM(3)/ZM_A)
               PDEP(J,NPERIO(I-1))   = PDEP(J,NPERIO(I-1))   + (ZM(0)/PDELPI(J,NPERIO(I-1)))*(ZM(4)/ZM_A)
               PDEP(J,NPERIO(I+1))   = PDEP(J,NPERIO(I+1))   + (ZM(0)/PDELPI(J,NPERIO(I+1)))*(ZM(5)/ZM_A)
               PDEP(J+1,NPERIO(I-1)) = PDEP(J+1,NPERIO(I-1)) + (ZM(0)/PDELPI(J+1,NPERIO(I-1)))*(ZM(6)/ZM_A)
               PDEP(J+1,I)           = PDEP(J+1,I)           + (ZM(0)/PDELPI(J+1,I))*(ZM(7)/ZM_A)
               PDEP(J+1,NPERIO(I+1)) = PDEP(J+1,NPERIO(I+1)) + (ZM(0)/PDELPI(J+1,NPERIO(I+1)))*(ZM(8)/ZM_A)
               PDEP(J-2,NPERIO(I-2)) = PDEP(J-2,NPERIO(I-2)) + ZM(0)*(ZM(9)/ZM_A)/PDELPI(J-2,NPERIO(I-2))
               PDEP(J-2,NPERIO(I-1)) = PDEP(J-2,NPERIO(I-1)) + ZM(0)*(ZM(10)/ZM_A)/PDELPI(J-2,NPERIO(I-1))
               PDEP(J-2,I)           = PDEP(J-2,I)           + ZM(0)*(ZM(11)/ZM_A)/PDELPI(J-2,I)
               PDEP(J-2,NPERIO(I+1)) = PDEP(J-2,NPERIO(I+1)) + ZM(0)*(ZM(12)/ZM_A)/PDELPI(J-2,NPERIO(I+1))
               PDEP(J-2,NPERIO(I+2)) = PDEP(J-2,NPERIO(I+2)) + ZM(0)*(ZM(13)/ZM_A)/PDELPI(J-2,NPERIO(I+2))
               PDEP(J-1,NPERIO(I-2)) = PDEP(J-1,NPERIO(I-2)) + ZM(0)*(ZM(14)/ZM_A)/PDELPI(J-1,NPERIO(I-2))
               PDEP(J-1,NPERIO(I+2)) = PDEP(J-1,NPERIO(I+2)) + ZM(0)*(ZM(15)/ZM_A)/PDELPI(J-1,NPERIO(I+2))
               PDEP(J,NPERIO(I-2))   = PDEP(J,NPERIO(I-2))   + ZM(0)*(ZM(16)/ZM_A)/PDELPI(J,NPERIO(I-2))
               PDEP(J,NPERIO(I+2))   = PDEP(J,NPERIO(I+2))   + ZM(0)*(ZM(17)/ZM_A)/PDELPI(J,NPERIO(I+2))
               PDEP(J+1,NPERIO(I-2)) = PDEP(J+1,NPERIO(I-2)) + ZM(0)*(ZM(18)/ZM_A)/PDELPI(J+1,NPERIO(I-2))
               PDEP(J+1,NPERIO(I+2)) = PDEP(J+1,NPERIO(I+2)) + ZM(0)*(ZM(19)/ZM_A)/PDELPI(J+1,NPERIO(I+2))
               PDEP(J+2,NPERIO(I-2)) = PDEP(J+2,NPERIO(I-2)) + ZM(0)*(ZM(20)/ZM_A)/PDELPI(J+2,NPERIO(I-2))
               PDEP(J+2,NPERIO(I-1)) = PDEP(J+2,NPERIO(I-1)) + ZM(0)*(ZM(21)/ZM_A)/PDELPI(J+2,NPERIO(I-1))
               PDEP(J+2,I)           = PDEP(J+2,I)           + ZM(0)*(ZM(22)/ZM_A)/PDELPI(J+2,I)
               PDEP(J+2,NPERIO(I+1)) = PDEP(J+2,NPERIO(I+1)) + ZM(0)*(ZM(23)/ZM_A)/PDELPI(J+2,NPERIO(I+1))
               PDEP(J+2,NPERIO(I+2)) = PDEP(J+2,NPERIO(I+2)) + ZM(0)*(ZM(24)/ZM_A)/PDELPI(J+2,NPERIO(I+2))
            END IF   
         END IF
      !----------------------------!
      !  DEALING WITH UNDERSHOOTS  !
      !----------------------------!   
      ELSE IF (PDEP(I,J) < ZMIN_LOC(I,J)) THEN
         ZM(0) = ZMIN_LOC(I,J)-PDEP(I,J)
         ZM(1) = MAX(PDEP(J-1,NPERIO(I-1))-ZMIN_LOC(J-1,NPERIO(I-1)),ZERO)*PDELPI(J-1,NPERIO(I-1))
         ZM(2) = MAX(PDEP(J-1,I)-ZMIN_LOC(J-1,I),ZERO)*PDELPI(J-1,I)
         ZM(3) = MAX(PDEP(J-1,NPERIO(I+1))-ZMIN_LOC(J-1,NPERIO(I+1)),ZERO)*PDELPI(J-1,NPERIO(I+1))
         ZM(4) = MAX(PDEP(J,NPERIO(I-1))-ZMIN_LOC(J,NPERIO(I-1)),ZERO)*PDELPI(J,NPERIO(I-1))
         ZM(5) = MAX(PDEP(J,NPERIO(I+1))-ZMIN_LOC(J,NPERIO(I+1)),ZERO)*PDELPI(J,NPERIO(I+1))
         ZM(6) = MAX(PDEP(J+1,NPERIO(I-1))-ZMIN_LOC(J+1,NPERIO(I-1)),ZERO)*PDELPI(J+1,NPERIO(I-1))
         ZM(7) = MAX(PDEP(J+1,I)-ZMIN_LOC(J+1,I),ZERO)*PDELPI(J+1,I)
         ZM(8) = MAX(PDEP(J+1,NPERIO(I+1))-ZMIN_LOC(J+1,NPERIO(I+1)),ZERO)*PDELPI(J+1,NPERIO(I+1))
         ZM_A  = ZERO 
         DO K=1,8
            ZM_A = ZM_A + ZM(K)
         END DO
         IF (ZM_A > ZM(0)) THEN
            PDEP(I,J)     = ZMIN_LOC(I,J)
            PDEP(J-1,NPERIO(I-1)) = PDEP(J-1,NPERIO(I-1)) - (ZM(0)/PDELPI(J-1,NPERIO(I-1)))*(ZM(1)/ZM_A)
            PDEP(J-1,I)           = PDEP(J-1,I)           - (ZM(0)/PDELPI(J-1,I))*(ZM(2)/ZM_A)
            PDEP(J-1,NPERIO(I+1)  = PDEP(J-1,NPERIO(I+1)) - (ZM(0)/PDELPI(J-1,NPERIO(I+1)))*(ZM(3)/ZM_A)
            PDEP(J,NPERIO(I-1))   = PDEP(J,NPERIO(I-1))   - (ZM(0)/PDELPI(J,NPERIO(I-1)))*(ZM(4)/ZM_A)
            PDEP(J,NPERIO(I+1))   = PDEP(J,NPERIO(I+1))   - (ZM(0)/PDELPI(J,NPERIO(I+1)))*(ZM(5)/ZM_A)
            PDEP(J+1,NPERIO(I-1)) = PDEP(J+1,NPERIO(I-1)) - (ZM(0)/PDELPI(J+1,NPERIO(I-1)))*(ZM(6)/ZM_A)
            PDEP(J+1,I)           = PDEP(J+1,I)           - (ZM(0)/PDELPI(J+1,I))*(ZM(7)/ZM_A)
            PDEP(J+1,NPERIO(I+1)) = PDEP(J+1,NPERIO(I+1)) - (ZM(0)/PDELPI(J+1,NPERIO(I+1)))*(ZM(8)/ZM_A)
         ELSE IF ((ZM_A .LE. ZM(0)).AND.(J.GE.3).AND.(J.LE.NLEV-2)) THEN
            ZM(9)  = MAX(PDEP(J-2,NPERIO(I-2))-ZMIN_LOC(J-2,NPERIO(I-2)),ZERO)*PDELPI(J-2,NPERIO(I-2))
            ZM(10) = MAX(PDEP(J-2,NPERIO(I-1))-ZMIN_LOC(J-2,NPERIO(I-1)),ZERO)*PDELPI(J-2,NPERIO(I-1))
            ZM(11) = MAX(PDEP(J-2,I)-ZMIN_LOC(J-2,I),ZERO)*PDELPI(J-2,I)
            ZM(12) = MAX(PDEP(J-2,NPERIO(I+1))-ZMIN_LOC(J-2,NPERIO(I+1)),ZERO)*PDELPI(J-2,NPERIO(I+1))
            ZM(13) = MAX(PDEP(J-2,NPERIO(I+2))-ZMIN_LOC(J-2,NPERIO(I+2)),ZERO)*PDELPI(J-2,NPERIO(I+2))
            ZM(14) = MAX(PDEP(J-1,NPERIO(I-2))-ZMIN_LOC(J-1,NPERIO(I-2)),ZERO)*PDELPI(J-1,NPERIO(I-2))
            ZM(15) = MAX(PDEP(J-1,NPERIO(I+2))-ZMIN_LOC(J-1,NPERIO(I+2)),ZERO)*PDELPI(J-1,NPERIO(I+2))
            ZM(16) = MAX(PDEP(J,NPERIO(I-2))-ZMIN_LOC(J,NPERIO(I-2)),ZERO)*PDELPI(J,NPERIO(I-2))
            ZM(17) = MAX(PDEP(J,NPERIO(I+2))-ZMIN_LOC(J,NPERIO(I+2)),ZERO)*PDELPI(J,NPERIO(I+2))
            ZM(18) = MAX(PDEP((J+1,NPERIO(I-2))-ZMIN_LOC((J+1,NPERIO(I-2)),ZERO)*PDELPI(J+1,NPERIO(I-2))
            ZM(19) = MAX(PDEP((J+1,NPERIO(I+2))-ZMIN_LOC((J+1,NPERIO(I+2)),ZERO)*PDELPI(J+1,NPERIO(I+2))
            ZM(20) = MAX(PDEP(J+2,NPERIO(I-2))-ZMIN_LOC(J+2,NPERIO(I-2)),ZERO)*PDELPI(J+2,NPERIO(I-2))
            ZM(21) = MAX(PDEP(J+2,NPERIO(I-1))-ZMIN_LOC(J+2,NPERIO(I-1)),ZERO)*PDELPI(J+2,NPERIO(I-1))
            ZM(22) = MAX(PDEP(J+2,I)-ZMIN_LOC(J+2,I),ZERO)*PDELPI(J+2,I)
            ZM(23) = MAX(PDEP(J+2,NPERIO(I+1))-ZMIN_LOC(J+2,NPERIO(I+1)),ZERO)*PDELPI(J+2,NPERIO(I+1))
            ZM(24) = MAX(PDEP(J+2,NPERIO(I+2))-ZMIN_LOC(J+2,NPERIO(I+2)),ZERO)*PDELPI(J+2,NPERIO(I+2))
            DO K=9,24
               ZM_A = ZM_A + ZM(K)
            END DO
            IF (ZM_A > ZM(0)) THEN
               PDEP(I,J)     = ZMIN_LOC(I,J)
               PDEP(J-1,NPERIO(I-1)) = PDEP(J-1,NPERIO(I-1)) - (ZM(0)/PDELPI(J-1,NPERIO(I-1)))*(ZM(1)/ZM_A)
               PDEP(J-1,I)           = PDEP(J-1,I)           - (ZM(0)/PDELPI(J-1,I))*(ZM(2)/ZM_A)
               PDEP(J-1,NPERIO(I+1)  = PDEP(J-1,NPERIO(I+1)) - (ZM(0)/PDELPI(J-1,NPERIO(I+1)))*(ZM(3)/ZM_A)
               PDEP(J,NPERIO(I-1))   = PDEP(J,NPERIO(I-1))   - (ZM(0)/PDELPI(J,NPERIO(I-1)))*(ZM(4)/ZM_A)
               PDEP(J,NPERIO(I+1))   = PDEP(J,NPERIO(I+1))   - (ZM(0)/PDELPI(J,NPERIO(I+1)))*(ZM(5)/ZM_A)
               PDEP(J+1,NPERIO(I-1)) = PDEP(J+1,NPERIO(I-1)) - (ZM(0)/PDELPI(J+1,NPERIO(I-1)))*(ZM(6)/ZM_A)
               PDEP(J+1,I)           = PDEP(J+1,I)           - (ZM(0)/PDELPI(J+1,I))*(ZM(7)/ZM_A)
               PDEP(J+1,NPERIO(I+1)) = PDEP(J+1,NPERIO(I+1)) - (ZM(0)/PDELPI(J+1,NPERIO(I+1)))*(ZM(8)/ZM_A)
               PDEP(J-2,NPERIO(I-2)) = PDEP(J-2,NPERIO(I-2)) - ZM(0)*(ZM(9)/ZM_A)/PDELPI(J-2,NPERIO(I-2))
               PDEP(J-2,NPERIO(I-1)) = PDEP(J-2,NPERIO(I-1)) - ZM(0)*(ZM(10)/ZM_A)/PDELPI(J-2,NPERIO(I-1))
               PDEP(J-2,I)           = PDEP(J-2,I)           - ZM(0)*(ZM(11)/ZM_A)/PDELPI(J-2,I)
               PDEP(J-2,NPERIO(I+1)) = PDEP(J-2,NPERIO(I+1)) - ZM(0)*(ZM(12)/ZM_A)/PDELPI(J-2,NPERIO(I+1))
               PDEP(J-2,NPERIO(I+2)) = PDEP(J-2,NPERIO(I+2)) - ZM(0)*(ZM(13)/ZM_A)/PDELPI(J-2,NPERIO(I+2))
               PDEP(J-1,NPERIO(I-2)) = PDEP(J-1,NPERIO(I-2)) - ZM(0)*(ZM(14)/ZM_A)/PDELPI(J-1,NPERIO(I-2))
               PDEP(J-1,NPERIO(I+2)) = PDEP(J-1,NPERIO(I+2)) - ZM(0)*(ZM(15)/ZM_A)/PDELPI(J-1,NPERIO(I+2))
               PDEP(J,NPERIO(I-2))   = PDEP(J,NPERIO(I-2))   - ZM(0)*(ZM(16)/ZM_A)/PDELPI(J,NPERIO(I-2))
               PDEP(J,NPERIO(I+2))   = PDEP(J,NPERIO(I+2))   - ZM(0)*(ZM(17)/ZM_A)/PDELPI(J,NPERIO(I+2))
               PDEP(J+1,NPERIO(I-2)) = PDEP(J+1,NPERIO(I-2)) - ZM(0)*(ZM(18)/ZM_A)/PDELPI(J+1,NPERIO(I-2))
               PDEP(J+1,NPERIO(I+2)) = PDEP(J+1,NPERIO(I+2)) - ZM(0)*(ZM(19)/ZM_A)/PDELPI(J+1,NPERIO(I+2))
               PDEP(J+2,NPERIO(I-2)) = PDEP(J+2,NPERIO(I-2)) - ZM(0)*(ZM(20)/ZM_A)/PDELPI(J+2,NPERIO(I-2))
               PDEP(J+2,NPERIO(I-1)) = PDEP(J+2,NPERIO(I-1)) - ZM(0)*(ZM(21)/ZM_A)/PDELPI(J+2,NPERIO(I-1))
               PDEP(J+2,I)           = PDEP(J+2,I)           - ZM(0)*(ZM(22)/ZM_A)/PDELPI(J+2,I)
               PDEP(J+2,NPERIO(I+1)) = PDEP(J+2,NPERIO(I+1)) - ZM(0)*(ZM(23)/ZM_A)/PDELPI(J+2,NPERIO(I+1))
               PDEP(J+2,NPERIO(I+2)) = PDEP(J+2,NPERIO(I+2)) - ZM(0)*(ZM(24)/ZM_A)/PDELPI(J+2,NPERIO(I+2))
            END IF   
         END IF
      END IF   
   ENDDO
ENDDO


END SUBROUTINE ILMC_FILTER 
!*************************************************************************************!
!*************************************************************************************!
!*************************************************************************************!
!*************************************************************************************!
!*************************************************************************************!
!*************************************************************************************!
SUBROUTINE BC_MASS_FIXER(PDEP,PARR,PXWEI,KXLAG,PYWEI,KYLAG,PMASS0)
  
IMPLICIT NONE

REAL(8),          INTENT(INOUT) :: PDEP(GXPT,GYPT)
REAL(8),          INTENT(IN)    :: PARR(GXPT,GYPT)
REAL(8),          INTENT(IN)    :: PMASS0
REAL(8),          INTENT(IN)    :: PXWEI(6,GXPT,GYPT),PYWEI(6,GXPT,GYPT)
INTEGER(8),       INTENT(IN)    :: KXLAG(6,GXPT,GYPT),KYLAG(6,GXPT,GYPT)

INTEGER(8)                      :: I,  J,  K,  L
REAL(8)                         :: ZMASS,ZDM,ZDEP_L,ZEPS
REAL(8)                         :: ZMIN_LOC,ZMAX_LOC,ZWS
REAL(8)                         :: ZW(GXPT,GYPT),ZDEP(GXPT,GYPT)

ZEPS = ABS(100*EPSILON(ONE))

!* Compute Mass error
ZMASS = ZERO
DO I=1,GXPT
   DO J=1,GYPT
      ZMASS = ZMASS + PDEP(I,J)
   ENDDO
ENDDO
ZDM = ZMASS - PMASS0


IF (ABS(ZDM) < ZEPS) THEN
   !* mass conserving coefficient
   DO I=1,GXPT
      DO J=1,GYPT
         ZDEP_L = ZERO
         DO K=5,6
            DO L=5,6
               ZDEP_L = ZDEP_L + PXWEI(K,I,J)*PYWEI(L,I,J)*(   &
                      &   PARR(KXLAG(K,I,J),KYLAG(L,I,J)) )
            ENDDO
         ENDDO
         ZW(I,J) = MAX(ZERO,SIGN(ZDM)*SIGN(PDEP(I,J)-ZDEP_L)*( &
                 & EXP(RQ*LOG(ABS(PDEP(I,J)-ZDEP_L))) ) )
      ENDDO
   ENDDO
   !* compute sum of weights
   ZWS = ZERO
   DO I=1,GXPT
      DO J=1,GYPT
         ZWS = ZWS + ZW(I,J)
      ENDDO
   ENDDO
   !* mass correction
   IF (ZWS > ZEPS) THEN
      DO I=1,GXPT
         DO J=1,GYPT
            PDEP(I,J) = PDEP(I,J) - (ZDM/ZWS)*ZW(I,J)
         ENDDO
      ENDDO
   END IF
   ! shape preservation
   DO I=1,GXPT
      DO J=1,GYPT
         ZMAX_LOC  =  MAX(                            &
                   & PARR(KXLAG(2,I,J),KYLAG(2,I,J)), &
                   & PARR(KXLAG(2,I,J),KYLAG(3,I,J)), &
                   & PARR(KXLAG(3,I,J),KYLAG(2,I,J)), &
                   & PARR(KXLAG(3,I,J),KYLAG(3,I,J))  )
         ZMIN_LOC  =  MIN(                            &
                   & PARR(KXLAG(2,I,J),KYLAG(2,I,J)), &
                   & PARR(KXLAG(2,I,J),KYLAG(3,I,J)), &
                   & PARR(KXLAG(3,I,J),KYLAG(2,I,J)), &
                   & PARR(KXLAG(3,I,J),KYLAG(3,I,J))  )
         ZDEP(I,J) = MAX(ZERO,SIGN(ZDM))*MAX(ZMIN_LOC,PDEP(I,J))  &
                   & + MAX(ZERO,SIGN(-ZDM))*MIN(ZMAX_LOC,PDEP(I,J))
      ENDDO
   ENDDO
   !* Compute Mass error and correction
   ZMASS = ZERO
   DO I=1,GXPT
      DO J=1,GYPT
         ZMASS = ZMASS + ZDEP(I,J)
      ENDDO
   ENDDO
   DO I=1,GXPT
      DO J=1,GYPT
         PDEP(I,J) = (PMASS0/ZMASS)*ZDEP(I,J)
      ENDDO
   ENDDO
END IF

END SUBROUTINE BC_MASS_FIXER 
!************************************************************************************************!
!************************************************************************************************!
!************************************************************************************************!
!************************************************************************************************!
!################################################################################################!
!************************************************************************************************!
!################################################################################################!
!************************************************************************************************!
!************************************************************************************************!

END MODULE MOD_SLAG
