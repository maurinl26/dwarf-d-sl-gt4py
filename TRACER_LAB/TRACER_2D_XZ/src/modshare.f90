!#################################  MODELE JOUET 2D  ##################################!
!#                                                                                    #!
!# auteurs : F,Voitus                                                                 #!
!# sujet   : Liste des fonctions servant a initialiser le programme                   #!
!#                                                                                    #!
!######################################################################################!

!=====================================================================!
!=========================  MODULE PRINCIPAL  ========================!
!=====================================================================!

MODULE MOD_SHARE

    INTEGER(8), PARAMETER              :: NLEV=50
    INTEGER(8), PARAMETER              :: NXPT=500
    INTEGER(8), PARAMETER              :: NORD1=2
    INTEGER(8), PARAMETER              :: NORD2=3
    REAL(8),    PARAMETER              :: RPI=3.14159265358979323846264338327950288419716939937510d0
    REAL(8),    PARAMETER              :: RCP=1004.709d0
    REAL(8),    PARAMETER              :: RG=9.80665d0
    REAL(8),    PARAMETER              :: RCV=(5.d0/7.d0)*RCP
    REAL(8),    PARAMETER              :: RR=RCP-RCV
    REAL(8),    PARAMETER              :: RKAPPA=(RR/RCP)
    REAL(8),    PARAMETER              :: RP00=100000.d0
    REAL(8),    PARAMETER              :: RT00=288.15d0
    REAL(8),    PARAMETER              :: RP=2.6d0
    REAL(8),    PARAMETER              :: TWO=2.d0
    REAL(8),    PARAMETER              :: ZERO=0.d0
    REAL(8),    PARAMETER              :: ONE=1.d0
    REAL(8),    PARAMETER              :: HALF=0.5d0
    REAL(8),    PARAMETER              :: REPS=1.d-10
    REAL(8),    PARAMETER              :: RVMAX=300.d0
    
    INTEGER(8)                         :: NTIMESTEP,NITMP,NSMAX,HALO,NCOMAD_OPT,NLBC
    INTEGER(8)                         :: NFREQOUTPUT,NTOP,NLEFT,NRIGHT,NCASE,NILEV
    INTEGER                            :: NFD_ORD

    REAL(8)                            :: RDX,RDZ,RDT,RLX,RDELTR,RA,RB,RC,RLD,RN,RS,RFD_COEF
    REAL(8)                            :: RXC,RZC,RAX,RAZ,RX0,RZ0,RHMAX,RMILIEU,RZ1,RPSTA,RZTOP
    REAL(8)                            :: RU00,RTSUR,RQPP,RPSUR,RTTROPO,RFCT,RZ2,RSLAG,RLZ
    REAL(8)                            :: RDZMAX,RDZMIN,RZLIM,RDTDZ_STD,RDZTOP,RLAMBDA,RZLOW
    REAL(8)                            :: RVFE_ALPHA,RVFE_BETA

    REAL(8), DIMENSION(4)              :: RFD_WEIGHT
    REAL(8), DIMENSION(NXPT)           :: RFDC
    
    LOGICAL                            :: LPERIO,LALIAZ,LADV_SPLIT,LSLAG,LSETTLS,LNESC,LREGETA
    LOGICAL                            :: LN000,LN001,LN002,LCOMAD,LPSTA,LREFINE_SLOROG,LPOST_Z
    LOGICAL                            :: LFCT_MONO,LADV_PIS,LPRINTLEV,LSPEC,LSLTRAJ_PREC
 
    CHARACTER(LEN=1)                   :: CTYPE_OROG,C_ADVECT,C_TRUNCS
    CHARACTER(LEN=3)                   :: CTYPE_AB_LEVEL
    CHARACTER(LEN=6)                   :: CTYPE_ETAH
    CHARACTER(LEN=1)                   :: C_SLTRAJ,C_SLINTP
    
CONTAINS

  !*********************************************************************!
  !*********************************************************************!
  !************************  LISTE DES ROUTINES  ***********************!
  !*********************************************************************!
  !*********************************************************************!
  
  SUBROUTINE INITIAL_SETTING

    IMPLICIT NONE

    INTEGER    :: IL,IXPT
    REAL(8)    :: ZR0,ZRI,ZL

    NAMELIST/NAMDIM/  NTIMESTEP,RDT,RDX,RDZ,NITMP,LSLAG,LSETTLS,LNESC,LCOMAD,  &
                    & LADV_SPLIT,LADV_PIS,LFCT_MONO,LSLTRAJ_PREC,NCOMAD_OPT
    NAMELIST/NAMPHY/  NCASE,CTYPE_OROG,RHMAX,LN000,LN001,LN002,LPSTA,LREGETA,  &
                    & RDZMAX,RZLIM
    NAMELIST/NAMDYN/  C_ADVECT,CTYPE_AB_LEVEL,CTYPE_ETAH,LPERIO,LALIAZ,RDELTR, &
                    & C_SLTRAJ,C_SLINTP,C_TRUNCS,RVFE_ALPHA,RVFE_BETA
    NAMELIST/NAMVIEW/ NFREQOUTPUT,LPRINTLEV,LPOST_Z
    
    OPEN(UNIT=4,FILE='init',STATUS='OLD')
    READ(4,NAMDIM)
    READ(4,NAMPHY)
    READ(4,NAMDYN)
    READ(4,NAMVIEW)
    CLOSE(4)
      
    NTOP    = 1
    NLEFT   = 1
    NRIGHT  = NXPT
    
    RLX     = REAL(NXPT,8)*RDX

    RTTROPO = 120.d0
    RA      = 25000.d0
    RLAMBDA = 8000.d0

    RAZ     = 3000.d0
    RZ1     = 4000.d0
    RZ2     = 5000.d0
   
    RX0     = -50000.d0
    RZ0     = 9000.d0
    RZC     = 0.d0
    RXC     = 0.d0

    RXC     = (RLX/TWO) + RXC
    RA      = RPI/RLAMBDA
    RB      = RPI/(TWO*RAX)
    RC      = ONE
    
    RFCT    = ZERO
    RS      = ZERO
    RSLAG   = ZERO
    RPSTA   = ZERO

    IF (LSLAG)      RSLAG = ONE
    IF (LFCT_MONO)  RFCT  = ONE
    IF (NCASE == 2) RS    = ONE
    IF (LPSTA)      RPSTA = ONE

    IF (.NOT.LSLAG) RDELTR = ZERO
    
    IF (LN002) THEN
       RN        = 0.02d0
       RTSUR     = ((RG**2)/(RN**2))/RCP
       RPSUR     = 100000.d0
    ELSE IF (LN001) THEN
       RN        = 0.01d0
       RTSUR     = 288.d0
       RPSUR     = 100000.d0
    ELSE IF (LN000) THEN
       RN        = 0.d0
       RTSUR     = 300.d0
       RPSUR     = 100000.d0
    ELSE
       RTSUR     = 288.15d0
       RPSUR     = 100000.d0
    END IF
        
    IF (NFREQOUTPUT .LT. 0) THEN
       NFREQOUTPUT = (-NFREQOUTPUT)*NTIMESTEP/100
    END IF
    IF (NFREQOUTPUT .EQ. 0) THEN
       NFREQOUTPUT = NTIMESTEP+1
    END IF


    SELECT CASE(NFD_ORD)
    CASE (2)
       HALO          = 1
       RFD_COEF      = 1.d0
       RFD_WEIGHT    = (/ 0.5d0/RDX, 0.d0, 0.d0, 0.d0 /)
    CASE (4)
       HALO          = 2
       RFD_COEF      = 1.37222d0
       RFD_WEIGHT    = (/ 2.d0/(3.d0*RDX), -1.d0/(12.d0*RDX), 0.d0, 0.d0 /)
    CASE (6)
       HALO          = 3
       RFD_COEF      = 1.58598d0
       RFD_WEIGHT    = (/ 0.75d0/RDX, -0.15d0/RDX, 1.d0/(60.d0*RDX), 0.d0 /)
    CASE (8)
       HALO          = 4
       RFD_COEF      = 1.7306d0
       RFD_WEIGHT    = (/ 0.8d0/RDX, -0.2d0/RDX, 4.d0/(105.d0*RDX), -1.d0/(280.d0*RDX) /)
    CASE DEFAULT
       HALO          = 0
       RFD_COEF      = RPI
       RFD_WEIGHT(:) = ZERO
    END SELECT


    IF (C_TRUNCS == "C") THEN
      NSMAX   = INT(NXPT/4)
    ELSE IF (C_TRUNCS == "Q") THEN
      NSMAX   = INT((NXPT-1)/3)
    ELSE 
      NSMAX   = INT(NXPT/2)
    ENDIF

    RMILIEU   = (RLX/TWO) + RXC

    RFDC(:) = 0.5d0
    IF (.NOT.LPERIO) THEN
       RFDC(NXPT) = 1.d0
       RFDC(1)    = 1.d0
    END IF   
    
    !----------------------------------------------------------------------------!
    ! vertical stretched height-coordinate

    RZTOP   = RDZ*REAL(NLEV-1,8)
    IF (.NOT.LREGETA) RZLOW = (RDZ/TWO)
    
    RDZMIN = RZLOW*(ONE+(((RR*RTSUR)/RG)/(((RR*RTSUR)/RG)-RZLOW)))
    RDZTOP = ((RR*RTTROPO)/RG)*(ONE+(RCP/RR))
    NILEV  = 1_8 + FLOOR(REAL(NLEV-1,8)*(RZLIM/TWO)*(RZLIM+3.d0)/(RZLIM+ONE))

    IF (CTYPE_AB_LEVEL == "SIG") THEN
       IF (LREGETA) THEN
          RLZ    = ZERO
          DO IL = 1,NLEV-1
             RLZ = RLZ + RDZ
          END DO
       ELSE
          ZRI    = RZLIM*RDZMAX+(ONE-RZLIM)*RDZMIN
          ZR0    = (REAL(NLEV-2,8)*ZRI-REAL(NILEV-1,8)*RDZMAX) &
                 & /REAL(NLEV-1-NILEV,8)
          RLZ    = ZERO 
          DO IL = 1,NILEV
             ZL  = REAL(IL-1,8)/REAL(NILEV-1,8)
             RLZ = RLZ + ZR0 + ZL*(ZRI-ZR0)            &
                 & + ((ONE-ZL)**2)*(RDZMIN-ZR0) 
          END DO
          DO IL=NILEV+1,NLEV-1
             ZL  = REAL(IL-1,8)/REAL(NLEV-2,8)
             RLZ = RLZ + ZR0 + ZL*(RDZMAX-ZR0) 
          END DO   
       END IF
       IF (RLZ.LT.RZTOP) THEN
          WRITE(*,*) 'INSUFFICIENT NUMBER OF NLEV :', RLZ, RZTOP 
          STOP
       END IF
    END IF   
   
    
  END SUBROUTINE INITIAL_SETTING
  !*********************************************************************!
  !**************  Cree une montagne au milieu du domaine  *************!
  FUNCTION RHEIGHT(PIN) RESULT(POUT)
    REAL(8)            :: PIN, POUT

    IF (CTYPE_OROG .EQ. "S") THEN
       POUT = RHMAX*(COS(RPI*PIN/RLAMBDA)**2) &
            & *(COS((RPI/(TWO*RA))*MAX(-RA,MIN(PIN,RA)))**2)   
    ELSE 
       POUT = 0.d0
    END IF

  END FUNCTION RHEIGHT
  !**********************************************************************!
  !**********************************************************************!
  FUNCTION DRHEIGHTDX(PIN) RESULT(POUT)
    REAL(8)            :: PIN, POUT

    IF (CTYPE_OROG .EQ. "S") THEN
       POUT = -RPI*RHMAX*( (SIN(TWO*RPI*PIN/RLAMBDA)&
            & *(COS((RPI/(TWO*RA))*MAX(-RA,MIN(PIN,RA)))**2)/RLAMBDA) &
            & +(SIN((RPI*PIN/RA)*MAX(-RA,MIN(PIN,RA)))&
            & *(COS(RPI*PIN/RLAMBDA)**2)/(TWO*RA)) ) 
    ELSE 
      POUT = 0.d0
    END IF

  END FUNCTION DRHEIGHTDX
  !*********************************************************************!
  !*********************************************************************!
  INTEGER FUNCTION NPERIO(IK)

  IMPLICIT NONE

  INTEGER(8)      :: IK
  LOGICAL         :: LLOK

  IF (LPERIO) THEN
     NPERIO=IK-1
     IF(NPERIO.GE.NXPT)THEN
       LLOK=.FALSE.
       DO WHILE (.NOT.LLOK)
         NPERIO=NPERIO-NXPT
         IF (NPERIO.LT.NXPT) LLOK=.TRUE.
       ENDDO
     ELSEIF (NPERIO.LT.0) THEN
       LLOK=.FALSE.
       DO WHILE (.NOT.LLOK) 
         NPERIO=NPERIO+NXPT
         IF(NPERIO.GE.0) LLOK=.TRUE.
       ENDDO
     ENDIF
     NPERIO=NPERIO+1
  ELSE
     NPERIO=IK
     IF (NPERIO.GE.NXPT) THEN
        NPERIO=NXPT
     ELSEIF (NPERIO.LE.1) THEN
        NPERIO=1
     ENDIF
  END IF

  END FUNCTION NPERIO
  !*********************************************************************!
  !*********************************************************************!
  INTEGER FUNCTION NBC(IK)

  IMPLICIT NONE

  INTEGER(8)      :: IK
  
  NBC=IK
  IF (NBC.GE.NLEV) THEN
     NBC=NLEV
  ELSEIF (NBC.LE.1) THEN 
     NBC=1
  ENDIF

  END FUNCTION NBC
  !*********************************************************************!
  !*********************************************************************!
  INTEGER FUNCTION NBL(IK)

  IMPLICIT NONE

  INTEGER(8)      :: IK
  
  NBL=IK
  IF (NBL.GE.NLEV+1) THEN
     NBL=NLEV+1
  ELSEIF (NBL.LE.1) THEN 
     NBL=1
  ENDIF

  END FUNCTION NBL
  !*********************************************************************!
  !*********************************************************************!
  !*********************************************************************!
  FUNCTION RBOYD(PIN) RESULT(POUT)
    REAL(8)            :: PIN
    REAL(8)            :: POUT
    REAL(8), PARAMETER :: B=2.5d0

    IF (PIN.GE.ONE) THEN
      POUT = ONE
    ELSE IF (PIN.LE.-ONE) THEN
      POUT = ZERO
    ELSE   
      POUT= (ONE+ERF(B*PIN/SQRT(ONE-PIN**2)))/TWO
    END IF

  END FUNCTION RBOYD
  !************************************************************************!
  !************************************************************************!
  !************************************************************************!
  !************************************************************************!
  FUNCTION MEANVAL(PVAL,KX) RESULT(POUT)

  IMPLICIT NONE

  REAL(8)    :: POUT
  REAL(8)    :: PVAL
  INTEGER(8) :: KX

  INTEGER(8) :: JI
  REAL(8)    :: ZSUM

  POUT=ZERO
  DO JI=1,KX
     POUT = POUT + (PVAL(JI)/REAL(KX,8))
  END DO   
  
  END FUNCTION MEANVAL
  !--------------------------------------------------
  
!************************************************************************!
!************************************************************************!
!************************************************************************!
!************************************************************************!
!************************************************************************!
!************************************************************************!

END MODULE MOD_SHARE

!=====================================================================!
!=====================================================================!
!=====================================================================!
