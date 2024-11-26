!##################   PARAMETRES   ###################!
!#                                                   #!
!# auteur : F.Voitus                                 #!
!# sujet  : Définit les parametres du modèle         #!
!#                                                   #!
!#####################################################!

!=====================================================!
!================  MODULE PARAMETRES  ================!
!=====================================================!

MODULE MOD_PARAMETER

  IMPLICIT NONE

  !**************  Liste des constantes  ***************!
  ! Dimension du domaine global
  INTEGER, PARAMETER   :: NLEV  = 21
  INTEGER, PARAMETER   :: GXPT  = 600
  INTEGER, PARAMETER   :: GYPT  = 10
  INTEGER, PARAMETER   :: GEXT  = 0
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
  REAL(8), PARAMETER :: RA     = 6371229.d0
  REAL(8), PARAMETER :: RB     = 2.d0
  REAL(8), PARAMETER :: ROMEGA = 7.292d-5
  REAL(8), PARAMETER :: RCORI0 = TWO*ROMEGA*DSIN((RPI/TWO)/TWO)
  REAL(8), PARAMETER :: RBETA0 = TWO*(ROMEGA/RA)*DCOS((RPI/TWO)/TWO)
  REAL(8), PARAMETER :: RGAMA0 = 0.005d0

  ! Constantes numériques
  REAL(8), PARAMETER :: RP   = 2.16d0
  REAL(8), PARAMETER :: RNU  = 0.05d0
  REAL(8), PARAMETER :: RTAU = ONE + ABS(10*EPSILON(ONE))
  !*****************************************************!
  
 
  INTEGER            :: NB_COEUR
  INTEGER            :: NB_PROCS
  INTEGER            :: LXPT,LYPT,HALO
  INTEGER            :: TRANCHEX,TRANCHEY
  INTEGER            :: NPARTXY,NPARTX,NPARTY
  INTEGER(8)         :: NDIM,NSMAX,MSMAX,NSMAX_ORO
  INTEGER(8)         :: NTIMESTEP,NFREQOUTPUT
  INTEGER(8)         :: NCASE,NSITER,NSLTER,MSMAX_ORO
  INTEGER            :: NORDER
  
  REAL(8)            :: RECHEANCE,RPIS2D,RZ1,RZ2,RFCT
  REAL(8)            :: RU00,RV00,RTPP,RQPP,RN,RP00,RT00,RB,RS
  REAL(8)            :: RX0,RY0,RZ0,RXC,RYC,RZC,RXA,RYA,RZT,RXT
  REAL(8)            :: RDT,RDX,RDY,RDZ,RLX,RLY,RXS,RYS,RZS,RLD
  REAL(8)            :: RPSUR,RTSUR,RHMAX,RSLAG,RDELTR,RPSTA
  
  CHARACTER(LEN=1)   :: C_ADVECT,C_TRUNCS,CTRUNCO
  CHARACTER(LEN=1)   :: C_OPTION,C_OROTYP
  
  LOGICAL            :: LPERIO,LRALIAZ,LSLAG,LADV_PIS2D,LFCTM
  LOGICAL            :: LSLICE_XY,LSLICE_XZ,LPRINTLEV,LPSTA

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

    NAMELIST/NAMDIM/ RECHEANCE,RDT,RDX,RDY,RDZ,NSITER,NSLTER, &
                   & NVERINT,LSLAG,LADV_PIS2D,LFCTM
    NAMELIST/NAMPHY/ NCASE,C_OPTION,RHMAX,LPSTA
    NAMELIST/NAMDYN/ C_ADVECT,C_TRUNCS,CTRUNCO,RDELTR
    NAMELIST/NAMCT0/ NFREQOUTPUT,LSLICE_XY,LSLICE_XZ,LPRINTLEV
    
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
    RAZ        = 3000.d0
    RZ1        = 4000.d0
    RZ2        = 5000.d0
    RLD        = 8000.d0
    RX0        = -50000.d0
    RZ0        = 9000.d0
    RZC        = 0.d0
    RXC        = 0.d0
    RA         = RPI/RLD
    RB         = RPI/(TWO*RAX)
    RXC        = RXC + RX0
    RYC        = RYC + RY0
    RYS        = RY0
    RFCT       = ZERO
    RS         = ZERO
    RSLAG      = 0.d0
    RPIS2D     = 1.d0
    RPSTA      = 0.d0 

    IF (.NOT.LADV_PIS2D) THEN
       RDELTR = ZERO
       IF (.NOT.LSLAG) RPIS2D = ZERO 
    END IF
    
    IF (LPSTA)  RPSTA = ONE
    IF (LSLAG)  RSLAG = ONE
    IF (LFCTM)  RFCT  = ONE
    IF (NCASE==2) RS  = ONE
    
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
    ELSE
       NSMAX = GXPT
       MSMAX = GYPT
    ENDIF

    IF (C_TRUNCO == "Q") THEN
       NSMAX_ORO = (GXPT-1)/3
       MSMAX_ORO = (GYPT-1)/3
       IF (MOD(NSMAX_ORO,2) /= 0) NSMAX_ORO = NSMAX_ORO - 1_8
       IF (MOD(MSMAX_ORO,2) /= 0) MSMAX_ORO = MSMAX_ORO - 1_8
    ELSE IF (C_TRUNCO == "L") THEN
       NSMAX_ORO = GXPT/2
       MSMAX_ORO = GYPT/2
       IF (MOD(NSMAX_ORO,2) /= 0) NSMAX_ORO = NSMAX_ORO - 1_8
       IF (MOD(MSMAX_ORO,2) /= 0) MSMAX_ORO = MSMAX_ORO - 1_8
    ELSE
       NSMAX_ORO = GXPT
       MSMAX_ORO = GYPT
    ENDIF
    
  END SUBROUTINE DEF_CONFIG_PARAMETER
  !***************************************************************!
  !*******  Défintion des dimension du domaine  ******************!
  !***************************************************************!
  SUBROUTINE DEF_MPI_DIMENSIONING(ABS_GLOBALOUT,ORD_GLOBALOUT,KEIGHOUR,NB_PROCSIN,KRANG)

    USE OMP_LIB

    IMPLICIT NONE

    INTEGER(8), DIMENSION(:),ALLOCATABLE,   INTENT(OUT) :: ABS_GLOBALOUT,ORD_GLOBALOUT
    INTEGER, DIMENSION(4),                  INTENT(OUT) :: KEIGHOUR
    INTEGER,                                INTENT(IN)  :: KRANG
    INTEGER,                                INTENT(IN)  :: NB_PROCSIN
    
    INTEGER                                             :: NX,NY,IQ

    !$OMP PARALLEL

    NB_COEUR = OMP_GET_NUM_THREADS()

    !$OMP END PARALLEL

    NB_PROCS = NB_PROCSIN

    ! ------------------------------ !
    ! DEFINITION DES DOMAINES LOCAUX !
    LYPT = GYPT/2
    LXPT = GXPT/2
    ! ------------------------------ !

    !*********************************************
    !* SL HALO DIMENSIONING & MPI-TASK DIMENSION
    !*********************************************
    
    HALO     = 0

    TRANCHEX = LXPT+2*HALO
    TRANCHEY = LYPT+2*HALO
    
    NDIM     = TRANCHEX*TRANCHEY

    
    ! Partition du domaine horizontal
    NPARTX = LXPT/NB_COEUR
    IF (NPARTX*NB_COEUR .NE. LXPT) THEN
       NPARTX = NPARTX + 1
    END IF
    NPARTY = LYPT/NB_COEUR
    IF (NPARTY*NB_COEUR .NE. LYPT) THEN
       NPARTY = NPARTY + 1
    END IF
    NPARTXY = LXPT*LYPT/NB_COEUR
    IF (NPARTXY*NB_COEUR .NE. LXPT*LYPT) THEN
       NPARTXY = NPARTXY + 1
    END IF

    ! Test de cohérence entre le nbe de sous-domaines et de procs
    ! RGGOBALOUT : Placement des sous-domaines sur le domaine

    ALLOCATE(ABS_GLOBALOUT(0:NB_PROCSIN-1))
    ALLOCATE(ORD_GLOBALOUT(0:NB_PROCSIN-1))
    
    IF ( (GYPT/LYPT)*(GXPT/LXPT) .NE. NB_PROCSIN ) THEN
       IF (KRANG .EQ. 0) THEN
          WRITE(*,*) 'ERREUR DE PARTITIONNEMENT DU DOMMAINE DANS DEF_MPI_DIMENSIONING'
          STOP
       END IF
    ELSE
       DO NY = 0, (GYPT/LYPT)-1
          DO NX = 0, (GXPT/LXPT)-1
             ABS_GLOBALOUT(NX+NY*(GXPT/LXPT)) = NX*LXPT
             ORD_GLOBALOUT(NX+NY*(GXPT/LXPT)) = NY*LYPT
          END DO
       END DO
    END IF

     ! Identification des voisins du proc X dans la variable KEIGHOUR :
    !   -KEIGHOUR(1) est son voisin du Nord
    !   -KEIGHOUR(2) ''  ''    ''   de l'Ouest
    !   -KEIGHOUR(3) ''  ''    ''   de l'Est
    !   -KEIGHOUR(4) ''  ''    ''   du Sud
    !        -------
    !        |     |
    !        |  1  |
    !        |     |
    !  -------     -------
    !  |                 |
    !  |  2     X     3  |
    !  |                 |
    !  -------     -------
    !        |     |
    !        |  4  |
    !        |     |
    !        -------

    IQ = KRANG/(GXPT/LXPT)

    KEIGHOUR(1) = MOD(KRANG+GXPT/LXPT,NB_PROCS)
    KEIGHOUR(4) = MOD(KRANG-GXPT/LXPT+NB_PROCS,NB_PROCS)
    KEIGHOUR(2) = MOD(KRANG-1+GXPT/LXPT,GXPT/LXPT)+IQ*(GXPT/LXPT)
    KEIGHOUR(3) = MOD(KRANG+1,GXPT/LXPT)+IQ*(GXPT/LXPT)
    

  END SUBROUTINE DEF_MPI_DIMENSIONING
  !****************************************************************************************!
  !****************************************************************************************!
  !****************************************************************************************!
  !****************************************************************************************!

  !#########################################################
  INTEGER FUNCTION NTAG(IFLD,IGPT,IPROC)

    IMPLICIT NONE

    INTEGER         :: IFLD,IPROC
    INTEGER(8)      :: IGPT

    NTAG = 10000*INT(IPROC,8) + 1000*INT(IFLD,8) + IGPT 

  END FUNCTION NTAG
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
  FUNCTION ROROG(PX,PY) RESULT(POUT)
    REAL(8)            :: PX,PY, POUT
    
    IF (C_OROTYP .EQ. "A") THEN
       POUT = RHMAX/( ONE + ((PX/RXA)**2) + ((PY/RYA)**2) )
    ELSE IF (C_OROTYP .EQ. "M") THEN
       POUT = RHMAX*(COS(RA*MAX(-RXA,MIN(PIN,RXA)))**2) &
            & *(COS(RB*MAX(-RXA,MIN(PIN,RXA)))**2)     
    ELSE 
       POUT = 0.d0
    END IF
  END FUNCTION ROROG
   !**********************************************************************!
  FUNCTION DOROGDX(PX,PY) RESULT(POUT)
    REAL(8)            :: PX,PY,POUT

    IF (C_OROTYP .EQ. "A") THEN
       POUT = -2.d0*RHMAX*(PX/(RXA**2)) &
            & /((1.d0+((PX/RXA)**2)+((PY/RYA)**2))**2)
    ELSE IF (C_OROTYP .EQ. "M") THEN
       POUT = - RHMAX*(                               &
            &   RB*SIN(TWO*RB*MAX(-RXA,MIN(PIN,RXA))) &
            &   *(COS(RA*MAX(-RXA,MIN(PIN,RXA)))**2)  &
            & + RA*SIN(TWO*RA*MAX(-RXA,MIN(PIN,RXA))) &
            &   *(COS(RB*MAX(-RXA,MIN(PIN,RXA)))**2) )
    ELSE 
      POUT = 0.d0
    END IF

  END FUNCTION DOROGDX
  !*************************************************************************
  FUNCTION DOROGDY(PX,PY) RESULT(POUT)
    REAL(8)            :: PX,PY,POUT

    IF (C_OROTYP .EQ. "A") THEN
       POUT = -2.d0*RHMAX*(PY/(RYA**2)) &
            & /((1.d0+((PX/RXA)**2)+((PY/RYA)**2))**2)
    ELSE IF (C_OROTYP .EQ. "M") THEN
       POUT = ZERO   
    ELSE 
      POUT = 0.d0
    END IF

  END FUNCTION DOROGDY
  !********************************************************************!
  !####################################################################!
  !####################################################################!
  !********************************************************************!
  INTEGER FUNCTION NPERIOX(IK)

  IMPLICIT NONE

  INTEGER(8)      :: IK
  LOGICAL         :: LLOK

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

  END FUNCTION NPERIOX
  !**********************************************
  INTEGER FUNCTION NPERIOY(IK)

  IMPLICIT NONE

  INTEGER(8)      :: IK
  LOGICAL         :: LLOK

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

  END FUNCTION NPERIOY
  !#####################################################################
  !#####################################################################
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

END MODULE MOD_PARAMETER

!=====================================================!
!=====================================================!
!=====================================================!
