!##################   PARAMETRES   ###################!
!#                                                   #!
!# auteur : F.Voitus                                 #!
!# sujet  : Définit les parametres du modèle         #!
!#                                                   #!
!#####################################################!

!=====================================================!
!================  MODULE PARAMETRES  ================!
!=====================================================!

MODULE MOD_PARAM

  IMPLICIT NONE

  !**************  Liste des constantes  ***************!
  ! Dimension du domaine global
  INTEGER, PARAMETER   :: NLEV  = 21
  INTEGER, PARAMETER   :: GXPT  = 600
  INTEGER, PARAMETER   :: GYPT  = 10
  INTEGER, PARAMETER   :: GEXT  = 0
  INTEGER, PARAMETER   :: NLBC  = 0
  INTEGER, PARAMETER   :: NORD1 = 2
  INTEGER, PARAMETER   :: NORD2 = 3

  ! Constantes algébrique
  REAL(8), PARAMETER :: ZERO = 0.d0
  REAL(8), PARAMETER :: ONE  = 1.d0
  REAL(8), PARAMETER :: TWO  = 2.d0
  REAL(8), PARAMETER :: HALF = 0.5d0
  REAL(8), PARAMETER :: RPI  = 3.14159265358979323846264338327950288419716939937510d0

  ! Constantes physique
  ! termodynamie des gaz parfait
  REAL(8), PARAMETER :: RCP    = 1004.709d0       ! 1004.5d0 baroclinic test value
  REAL(8), PARAMETER :: RCV    = (5.d0/7.d0)*RCP
  REAL(8), PARAMETER :: RR     = RCP-RCV
  REAL(8), PARAMETER :: RKAPPA = RR/RCP
  REAL(8), PARAMETER :: RGAMMA = RCP/RCV          ! 0.005d0 baroclinic test value
  REAL(8), PARAMETER :: RPATM  = 101325.d0
  REAL(8), PARAMETER :: RZ00   = 500.d0
  REAL(8), PARAMETER :: RP00   = 100000.d0
  REAL(8), PARAMETER :: RT00   = 288.15d0
  ! pesanteur
  REAL(8), PARAMETER :: RG     = 9.80665d0        ! 9.80616d0 baroclinic test value
  REAL(8), PARAMETER :: RAT    = 6371229.d0
  REAL(8), PARAMETER :: ROMEGA = 7.292d-5
  REAL(8), PARAMETER :: RGAMA0 = 0.005d0

  ! Constantes numériques
  REAL(8), PARAMETER :: RP   = 2.16d0
  REAL(8), PARAMETER :: RNU  = 0.05d0
  REAL(8), PARAMETER :: REPS = ABS(1000*EPSILON(ONE)) 
  REAL(8), PARAMETER :: RTAU = ONE + ABS(10*EPSILON(ONE))
  !*****************************************************!
  

  INTEGER(8)         :: NSMAX,MSMAX
  INTEGER(8)         :: NTIMESTEP,NFREQOUTPUT
  INTEGER            :: NITMP,NORDER,NCOMAD_OPT,NCASE
  
  REAL(8)            :: RECHEANCE,RFCT
  REAL(8)            :: RU00,RV00,RQ00,RN,RS,RT,RAX,RAY
  REAL(8)            :: RX0,RY0,RXC,RYC,RXA,RYA,RXT,RYT
  REAL(8)            :: RDT,RDX,RDY,RLX,RLY,RXS,RYS,RLD
  REAL(8)            :: RPSUR,RHMAX,RDELTR,RPSI0,RX1,RY1
  
  CHARACTER(LEN=1)   :: C_ADVECT,C_TRUNCS
  CHARACTER(LEN=1)   :: C_SLTRAJ,C_SLINTP
  
  LOGICAL            :: LPERIO,LALIAZ,LSLAG,LFCTM,LNESC,LSETTLS
  LOGICAL            :: LPRINTLEV,LSPEC,LCOMAD,LSLTRAJ_PREC

CONTAINS

  !*****************************************************!
  !***************  LISTE DES ROUTINES  ****************!
  !*****************************************************!


  !*******************************************************************!
  !**   Initialise les variables du domaine, physique et dynamique  **!
  !*******************************************************************!
  SUBROUTINE DEF_CONFIG_PARAMETER

    IMPLICIT NONE

    INTEGER :: NBCOEUR,NINTER
    REAL(8) :: ZXC,ZYC

    NAMELIST/NAMDIM/ RECHEANCE,RDT,RDX,RDY,LSPEC,LALIAZ,C_TRUNCS
    NAMELIST/NAMPHY/ NCASE
    NAMELIST/NAMDYN/ C_ADVECT,LSLAG,NITMP,C_SLTRAJ,C_SLINTP,LNESC, &
                   & LSETTLS,LCOMAD,LSLTRAJ_PREC,RDELTR,LFCTM
    NAMELIST/NAMCT0/ NFREQOUTPUT,LPRINTLEV
    
    OPEN(UNIT=10,FILE='init',STATUS='OLD')
    READ(10,NAMDIM)
    READ(10,NAMPHY)
    READ(10,NAMDYN)
    READ(10,NAMCT0)
    CLOSE(10)

    NTIMESTEP  = CEILING(REAL(RECHEANCE*3600,8)/RDT)

    RLY        = RDY*REAL(GYPT,8)
    RLX        = RDX*REAL(GXPT,8)
    RX0        = REAL(GXPT,8)*RDX/TWO
    RY0        = REAL(GYPT,8)*RDY/TWO
    RAX        = 25000.d0
    RAY        = 3000.d0
    RLD        = 8000.d0
    RX0        = -50000.d0
    RY0        = 9000.d0
    RXC        = 0.d0
    RYC        = 0.d0
    !RA         = RPI/RLD
    !RB         = RPI/(TWO*RAX)
    RXC        = RXC + RX0
    RYC        = RYC + RY0
    RYS        = RY0
    RFCT       = ZERO
    RS         = ZERO 

    IF (.NOT.LSLAG) RDELTR = ZERO 
    IF (LFCTM)      RFCT  = ONE
    IF (NCASE==2)   RS  = ONE
    
    IF (NFREQOUTPUT .LT. 0) THEN
       NFREQOUTPUT = (-NFREQOUTPUT)*NTIMESTEP/100
    END IF
    IF (NFREQOUTPUT .EQ. 0) THEN
       NFREQOUTPUT = NTIMESTEP+1
    END IF

    IF (C_TRUNCS == "Q") THEN
       NSMAX = (GXPT-1)/3
       MSMAX = (GYPT-1)/3
       IF (MOD(NSMAX,2) /= 0) NSMAX = NSMAX - 1_8
       IF (MOD(MSMAX,2) /= 0) MSMAX = MSMAX - 1_8
    ELSE IF (C_TRUNCS == "L") THEN
       NSMAX = GXPT/2
       MSMAX = GYPT/2
       IF (MOD(NSMAX,2) /= 0) NSMAX = NSMAX - 1_8
       IF (MOD(MSMAX,2) /= 0) MSMAX = MSMAX - 1_8
    ENDIF
    
  END SUBROUTINE DEF_CONFIG_PARAMETER
  !****************************************************************************************!
  !****************************************************************************************!
  !****************************************************************************************!
  !****************************************************************************************!

  !#########################################################
  FUNCTION RBOYD(PIN) RESULT(POUT)
    REAL(8)            :: PIN
    REAL(8)            :: POUT
    REAL(8), PARAMETER :: B=2.5

    IF (PIN.GE.ONE) THEN
      POUT = ONE
    ELSE IF (PIN.LE.-ONE) THEN
      POUT = ZERO
    ELSE   
      POUT= (ONE+ERF(B*PIN/SQRT(ONE-PIN**2)))/TWO
    END IF

  END FUNCTION RBOYD
  !########################################################
  !########################################################
  
  !**************  Cree une montagne au milieu du domaine  *************!
  FUNCTION STREAM_FUNCTION(PX,PY,PT) RESULT(POUT)
    REAL(8)            :: PX,PY,PT, POUT
    
    IF ((NCASE == 1).OR.(NCASE == 2)) THEN
       POUT = (4.d0*RPI/RT)*( ((((PX-RXC)**2)+((PY-RYC)**2))/TWO)          &
            & +COS(TWO*RPI*PT/RT)*( (((((PX-RXC)**2)+((PY-RYC)**2))/TWO)   &
            & +(1.d0/96.d0)*LOG(ONE-16.d0*((((PX-RXC)**2)+((PY-RYC)**2)))  &
            & + 25.d0*((((PX-RXC)**2)+((PY-RYC)**2))**2))                  &
            & -(1.d0/48.d0)*LOG(ONE+16.d0*((((PX-RXC)**2)+((PY-RYC)**2)))) &
            & -(SQRT(3.d0)/48.d0)*ATAN((32.d0*((((PX-RXC)**2)              &
            & +((PY-RYC)**2)))-ONE)/SQRT(3.d0)))) )
    ELSE IF (NCASE == 3) THEN
       POUT = (RPSI0/RT)*((RLX/(TWO*RPI))**2)*(                            &
            & SIN(TWO*RPI*(((PX-RXC)/RLX)-(PT/RT)))**2)*(                  &
            & COS(RPI*(PT/RT))*(COS(RPI*((PY-RYC)/RLY))**2))               &
            & -(RLX/RT)*(PY-RYC)
    ELSE 
       POUT = 0.d0
    END IF
  END FUNCTION STREAM_FUNCTION
  !**********************************************************************!
  !********************************************************************!
  !####################################################################!
  !####################################################################!
  !********************************************************************!
  INTEGER FUNCTION NPERIOX(IK)

  IMPLICIT NONE

  INTEGER(8)      :: IK
  LOGICAL         :: LLOK

  IF (LPERIO) THEN
     NPERIOX=IK-1
     IF(NPERIOX.GE.GXPT)THEN
       LLOK=.FALSE.
       DO WHILE (.NOT.LLOK)
         NPERIOX=NPERIOX-GXPT
         IF (NPERIOX.LT.GXPT) LLOK=.TRUE.
       ENDDO
     ELSEIF (NPERIOX.LT.0) THEN
       LLOK=.FALSE.
       DO WHILE (.NOT.LLOK) 
         NPERIOX=NPERIOX+GXPT
         IF(NPERIOX.GE.0) LLOK=.TRUE.
       ENDDO
     ENDIF
     NPERIOX=NPERIOX+1
  ELSE
     NPERIOX=IK
     IF (NPERIOX.GE.GXPT) THEN
        NPERIOX=GXPT
     ELSEIF (NPERIOX.LE.1) THEN
        NPERIOX=1
     ENDIF
  END IF
     
  END FUNCTION NPERIOX
  !**********************************************
  INTEGER FUNCTION NPERIOY(IK)

  IMPLICIT NONE

  INTEGER(8)      :: IK
  LOGICAL         :: LLOK

  IF (LPERIO) THEN
     NPERIOY=IK-1
     IF(NPERIOY.GE.GYPT)THEN
       LLOK=.FALSE.
       DO WHILE (.NOT.LLOK)
         NPERIOY=NPERIOY-GYPT
         IF (NPERIOY.LT.GYPT) LLOK=.TRUE.
       ENDDO
     ELSEIF (NPERIOY.LT.0) THEN
       LLOK=.FALSE.
       DO WHILE (.NOT.LLOK) 
         NPERIOY=NPERIOY+GYPT
         IF(NPERIOY.GE.0) LLOK=.TRUE.
       ENDDO
     ENDIF
     NPERIOY=NPERIOY+1
  ELSE
     NPERIOY=IK
     IF (NPERIOY.GE.GYPT) THEN
        NPERIOY=GYPT
     ELSEIF (NPERIOY.LE.1) THEN
        NPERIOY=1
     ENDIF
  END IF   
     
  END FUNCTION NPERIOY
  !#####################################################################!
  !#####################################################################!

END MODULE MOD_PARAM

!=====================================================!
!=====================================================!
!=====================================================!
