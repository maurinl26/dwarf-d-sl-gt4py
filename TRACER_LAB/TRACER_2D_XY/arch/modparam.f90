!##################   PARAMETRES   ###################!
!#                                                   #!
!# auteur : C.Colavolpe & F.Voitus                   #!
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
  INTEGER, PARAMETER   :: GXPT  = 96
  INTEGER, PARAMETER   :: GYPT  = 96
  INTEGER, PARAMETER   :: GEXT  = 0
  INTEGER, PARAMETER   :: NLBCX = 8
  INTEGER, PARAMETER   :: NLBCY = 8

  ! Constantes algébrique
  REAL(8), PARAMETER :: ZERO = 0.d0
  REAL(8), PARAMETER :: ONE  = 1.d0
  REAL(8), PARAMETER :: TWO  = 2.d0
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
  INTEGER            :: TRANCHEX,TRANCHEY,NDIM
  INTEGER            :: NPARTXY,NPARTX,NPARTY
  INTEGER            :: NHALO,NSLICEX,NSLICEY,NDIMLOC

  INTEGER(8)         :: NTIMESTEP
  INTEGER(8)         :: NFREQOUTPUT,NFPLEV,NFIRST,NDLNPR
  INTEGER(8)         :: NHD_ORDER,NZBC,NLBC,NSPG1,NSPG2
  INTEGER(8)         :: NSMAX,MSMAX,NSMAX_ORO,MSMAX_ORO,NBDIAG
  INTEGER(8)         :: NSITER,NSLTER
  INTEGER            :: NORDER
  
  REAL(8)            :: RKDF,RECHEANCE
  REAL(8)            :: RSITR,RSIPR,RSITRA
  REAL(8)            :: RHST,RCST2,RNST2
  REAL(8)            :: RU00,RV00,RUPP,RTPP,RTHETAS,RN,RP00,RT00
  REAL(8)            :: RX0,RY0,RZ0,RXC,RYC,RZC,RXA,RYA,RZT,RXT
  REAL(8)            :: RDT,RDX,RDY,RDZ,RLX,RLY,RXS,RYS,RZS
  REAL(8)            :: RDXTAU,RDAMPDIV,RDAMPVOR,RDAMPT,RDAMPDW
  REAL(8)            :: REGVINT,RPSUR,RTSUR,RTTROPO,RSMAX
  REAL(8)            :: RSPONX,RSPONT,RMUV,RMUW,RDAMPZ,RCSH
  REAL(8)            :: RHMAX,REXPT,RCC,RBBC,RNOFL,RWW,RDELTR

  REAL(8), DIMENSION(4) :: RWEIGHT
  
  CHARACTER(LEN=2)   :: C_ADJUST,C_ADVECT
  CHARACTER(LEN=1)   :: C_OPTION,C_TRUNCS,C_TRUNCO,C_OROTYP
 
  LOGICAL            :: LOROG,LBUBBLE,LNZERO,LN001,LN002 
  LOGICAL            :: LPERIO,LNUDGING,LSPONGE,LHDIFF,LVDIFF
  LOGICAL            :: LBIGW,LNHGRAD,LRALIAZ,LALGOR
  LOGICAL            :: LVEREG,LPENTADIAG,LNOFLOW,LSHEAR
  LOGICAL            :: LBUBBLE_HOT,LBUBBLE_WARM
  LOGICAL            :: LOROG_SCHAR,LOROG_GAUSS,LOROG_AGNESI
  LOGICAL            :: LSLICE_XY,LSLICE_XZ

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

    NAMELIST/NAMDIM/ RECHEANCE,RDT,RDX,RDY,RDZ,NDLNPR,NSITER,NSLTER,LRALIAZ,LALGOR,LBIGW,LVEREG
    NAMELIST/NAMPHY/ LOROG_AGNESI,LOROG_GAUSS,LOROG_SCHAR,LBUBBLE_HOT,LBUBBLE_WARM,LNOFLOW,LSHEAR
    NAMELIST/NAMDYN/ C_ADJUST,C_ADVECT,C_OPTION,C_TRUNCS,C_TRUNCO,RSITR,RSIPR,RSITRA, &
                   & LHDIFF,LVDIFF,LSPONGE,LNUDGING,RDELTR,RSMAX
    NAMELIST/NAMCT0/ NFREQOUTPUT,NFIRST,NFPLEV,LSLICE_XY,LSLICE_XZ
    
    OPEN(UNIT=10,FILE='init',STATUS='OLD')
    READ(10,NAMDIM)
    READ(10,NAMPHY)
    READ(10,NAMDYN)
    READ(10,NAMCT0)
    CLOSE(10)

    NTIMESTEP  = CEILING(REAL(RECHEANCE*3600,8)/RDT)

    RLY        = RDY*REAL(GYPT,8)
    RLX        = RDX*REAL(GXPT,8)
    RX0        = REAL(GXPT-GEXT,8)*RDX/TWO
    RY0        = REAL(GYPT-GEXT,8)*RDY/TWO

    ! Parametrisation des coefficients de la dérivée horziontale
    SELECT CASE(ORDRE)
    CASE (2)
       NHALO   = 1
       RKDF    = 1.d0
       RWEIGHT = (/ 0.5d0, 0.d0, 0.d0, 0.d0 /)
    CASE (4)
       NHALO   = 2
       RKDF    = 1.37222d0
       RWEIGHT = (/ 2.d0/3.d0, -1.d0/12.d0, 0.d0, 0.d0 /)
    CASE (6)
       NHALO   = 3
       RKDF    = 1.58598d0
       RWEIGHT = (/ 0.75d0, -0.15d0, 1.d0/60.d0, 0.d0 /)
    CASE (8)
       NHALO   = 4
       RKDF    = 1.7306d0
       RWEIGHT = (/ 0.8d0, -0.2d0, 4.d0/105.d0, -1.d0/280.d0 /)
    CASE DEFAULT
       WRITE(*,*) "Erreur de lecture de l'ordre de discrétisation"
    END SELECT
    
    RCST2      = RR*RSITR*(RCP/RCV)
    RHST       = (RR*RSITR)/RG
    RNST2      = (RG/RSITR)*(RG/RCP)

    RSPONX     = 1.d0 - 10d-2
    RSPONT     = 5.d0
    REXPT      = 1.d0
    RCC        = 1.d0
    RWW        = 0.d0
    NLBC       = 8_8
    NZBC       = 8_8
    RTTROPO    = 120.d0
    RP00       = 100000.d0
    RT00       = 288.15d0
   
    RDAMPDIV   = 1.d0
    RDAMPVOR   = 5.d0
    RDAMPT     = 5.d0
    RDAMPDW    = 1.d0
    REGVINT    = 0.d0
    RBBC       = 0.d0
    RNOFL      = 0.d0
    RCSH       = 0.d0
    NHD_ORDER  = 4_8
    NBDIAG     = 3
    LPENTADIAG = .FALSE.
    LNHGRAD    = .FALSE.

    IF (LBIGW) THEN
       LPENTADIAG = .TRUE.
       LNHGRAD    = .TRUE.
       NBDIAG     = 5
       RBBC       = ONE
       RWW        = ONE
    END IF   

    IF (LVEREG)  REGVINT = ONE
    IF (LNOFLOW) RNOFL   = ONE
    IF (LSHEAR)  RCSH    = 0.00025d0
    
    LBUBBLE = (LBUBBLE_WARM).OR.(LBUBBLE_HOT)
    LOROG   = (LOROG_AGNESI).OR.(LOROG_SCHAR).OR.(LOROG_GAUSS)
    
    LNZERO  = (LBUBBLE_WARM)
    LN001   = (LBUBBLE_HOT).OR.(LOROG_SCHAR)
    LN002   = (LOROG_AGNESI).OR.(LOROG_GUASS)

    IF (LOROG_AGNESI) THEN
       !NLEV      = 100
       !NXPT      = 360
       !RDT       = 2.d0
       !NTIMESTEP = CEILING(REAL(8000,8)/RDT)
       !RDX       = 200.d0
       !RDZ       = 300.d0
       RXA       = 2500.d0
       RYA       = 1000.d0
       RHMAX     = RSMAX         !500.d0
       C_OROTYP  = 'A'
       RU00      = 10.d0
       RV00      = 0.d0
       RN        = 0.02d0
       RTSUR     = ((RG**2)/(RN**2))/RCP
       RPSUR     = 101325.d0
       RXC       = 0.d0
       RYC       = 0.d0
       RZC       = 0.d0
       RXT       = 0.d0
       RZT       = 0.d0
       NLBC      = 16_8
       NZBC      = 35
       NSPG1     = 5
       NSPG2     = 35
       RSPONX    = 1.d0 - 75d-2
       RSPONT    = 5.d0
       LPERIO    = .FALSE.
       LHDIFF    = .TRUE.
       LSPONGE   = .TRUE.
       LNUDGING  = .FALSE.
       RDXTAU    = 100.d0
       RDAMPZ    = 0.d0
       RTTROPO   = 140.d0
       RMUW      = 1.2d0 
       RMUV      = 0.d0
    END IF
    ! Schar orography (Schär (2002))
    ! NXPT=1000, NLEV=65, avec RDX=500m et RDZ=312m
    IF (LOROG_SCHAR) THEN
       !NLEV      = 65
       !NXPT      = 1000
       !RDX       = 500
       !RDZ       = 312
       !RDT       = 4
       !NTIMESTEP = CEILING(REAL(5*3600,8)/RDT)
       RXA       = 5000.d0
       RYA       = 10000.d0
       RHMAX     = 250.d0
       C_OROTYP  = 'S'
       RU00      = 10.d0
       RV00      = 0.d0
       RTSUR     = 288.15d0
       RPSUR     = 100000.d0
       RN        = 0.01d0
       RZC       = 0.d0
       RXC       = 0.d0
       RYC       = 0.d0
       RXT       = 0.d0
       RZT       = 0.d0
       NSPG1     = 5
       NSPG2     = 30
       NZBC      = 25
       NLBC      = 16_8
       LPERIO    = .FALSE.
       LHDIFF    = .FALSE.
       LSPONGE   = .TRUE.
       LNUDGING  = .FALSE.
       RSPONX    = 1.d0 - 25d-2
       RSPONT    = 5.d0
       RDAMPZ    = 0.d0
       RDXTAU    = 0.d0
       RMUW      = 1.2d0
       RMUV      = 0.d0 
    END IF
    ! no flow orography test (Klemp (2011))
    ! NXPT=100, NLEV=40, avec RDX=500m et RDZ=500m
    IF (LOROG_GAUSS) THEN
       !NLEV      = 40
       !NXPT      = 100
       !RDX       = 500
       !RDZ       = 500
       !RDT       = 4
       RXA       = 2000.d0
       RYA       = 2000.d0
       RHMAX     = RSMAX          !6000 
       C_OROTYP  = 'G'
       RU00      = (ONE-RNOFL)*20.d0
       RV00      = 0.d0
       RN        = 0.02d0
       RTSUR     = ((RG**2)/(RN**2))/RCP
       RPSUR     = 100000.d0
       RZC       = 2000.d0
       RXT       = 0.d0
       RZT       = 3000.d0
       RXC       = 0.d0
       NSPG1     = 5
       NSPG2     = 20
       NZBC      = 40
       NLBC      = 16_8
       LPERIO    = .FALSE.
       LHDIFF    = .TRUE.
       LVDIFF    = .FALSE.
       LSPONGE   = .FALSE.
       LNUDGING  = .FALSE.
       RSPONX    = 1.d0 - 25d-2
       RSPONT    = 5.d0
       RDXTAU    = 100.d0
       RDAMPZ    = 0.d0
       RMUW      = 0.d0 
       RMUV      = 0.d0 
    END IF
    ! Hot bubble (Klemp (1994))
    ! NLEV=20, NXPT=600, RDX=500m,et RDZ=500m
    IF (LBUBBLE_HOT) THEN
       !NLEV      = 21
       !NXPT      = 600
       !RDX       = 500
       !RDZ       = 500
       !RDT       = 6
       !NTIMESTEP = 500
       C_OROTYP  = 'O'
       RHMAX     = 0.d0
       RU00      = 20.d0
       RV00      = 0.d0
       RPSUR     = 100000.d0
       RTHETAS   = 300.d0
       RTSUR     = RTHETAS*((RPSUR/RP00)**(RKAPPA))
       RN        = 0.01d0
       RZC       = 0.d0
       RXT       = 0.d0
       RZT       = 10000.d0
       RXC       = -60000.d0
       RYC       = 0.d0
       RP00      = 100000.d0
       RT00      = (RG**2)/(RN**2)/RCP
       RTPP      = 0.01d0
       RUPP      = 0.d0
       RZT       = 10000.d0
       RYC       = 0.d0
       RXC       = -60000.d0
       RXA       = 5000.d0
       RYA       = 1000.d0
       RXC       = RXC + RX0
       RYC       = RYC + RY0
       NSPG1     = 0
       NSPG2     = 1
       LPERIO    = .TRUE.
       LHDIFF    = .FALSE.
       LVDIFF    = .FALSE.
       LSPONGE   = .FALSE.
       LNUDGING  = .FALSE.
       RDXTAU    = 0.d0
       RDAMPZ    = 0.d0
    END IF
    ! Warm bubble (Janjic (2001))
    IF (LBUBBLE_WARM) THEN
       ! NLEV      = 80
       ! NXPT      = 160
       ! RDZ       = 125
       ! RDX       = 125
       ! RDT       = 2
       ! RDAMPX    = 300.d0
       C_OROTYP  = 'O'
       RHMAX     = 0.d0
       RU00      = 0.d0
       RTSUR     = 300.d0
       RPSUR     = 86100.d0
       RTHETAS   = RTSUR*((RP00/RPSUR)**(RKAPPA))
       RTPP      = 2.d0
       RN        = 0.d0
       RZC       = 2000.d0
       RXT       = 2000.d0
       RZT       = 2000.d0
       RXC       = 0.d0
       RYC       = 0.d0
       RXA       = 0.d0
       RYA       = 0.d0
       NSPG1     = 5
       NSPG2     = 15
       LPERIO    = .TRUE.
       LHDIFF    = .TRUE.
       LVDIFF    = .TRUE.
       LSPONGE   = .FALSE.
       LNUDGING  = .FALSE.
       RDXTAU    = 0.d0
       RDAMPZ    = 0.025d0
       NHD_ORDER = 2_8 
    END IF
    
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
  SUBROUTINE DEF_MPI_DIMENSIONING(ABS_GLOBALOUT,ORD_GLOBALOUT,NB_PROCSIN,RANGIN)

    USE OMP_LIB

    IMPLICIT NONE

    INTEGER(8), DIMENSION(:),ALLOCATABLE,   INTENT(OUT) :: ABS_GLOBALOUT,ORD_GLOBALOUT
    INTEGER,                                INTENT(IN)  :: RANGIN
    INTEGER,                                INTENT(IN)  :: NB_PROCSIN
    INTEGER                                             :: NX,NY

    !$OMP PARALLEL

    NB_COEUR = OMP_GET_NUM_THREADS()

    !$OMP END PARALLEL

    NB_PROCS = NB_PROCSIN

    ! ------------------------------ !
    ! DEFINITION DES DOMAINES LOCAUX !
    LYPT = GYPT/2
    LXPT = GXPT/2
    ! ------------------------------ !

    HALO = 0

    TRANCHEX = LXPT+2*HALO
    TRANCHEY = LYPT+2*HALO

    NDIM = TRANCHEX*TRANCHEY

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
       IF (RANGIN .EQ. 0) THEN
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

  END SUBROUTINE DEF_MPI_DIMENSIONING
  !****************************************************************************************!
  !****************************************************************************************!
  !**************************  Défintion des dimension du domaine *************************!
  !****************************************************************************************!
  SUBROUTINE DEF_MPI_DIMENSION(KGABS,KGORD,KEIGHOUR,KNBPROCS,KRANG)

    IMPLICIT NONE

    INTEGER, DIMENSION(:),ALLOCATABLE, INTENT(OUT) :: KGABS,KGORD
    INTEGER, DIMENSION(4),             INTENT(OUT) :: KEIGHOUR
    INTEGER,                           INTENT(IN)  :: KNBPROCS,KRANG
    
    INTEGER                                        :: IX,IY,IQ

    ! Mise en mémoire globale du nombre de processeur
    NBPROCS = KNBPROCS

    ! ------------------------------ !
    ! DEFINITION DES DOMAINES LOCAUX !
    LYPT = GYPT/2
    LXPT = GXPT/2
    ! ------------------------------ !

    ! Test sur les dimensions du domaine pour que chaque
    ! sous-domaine des processeurs aient la même taille 
    IF ( (GYPT/LYPT)*(GXPT/LXPT) .NE. NB_PROCSIN ) THEN
       IF (RANGIN .EQ. 0) THEN
          WRITE(*,*) 'ERREUR DE PARTITIONNEMENT DU DOMMAINE DANS DEF_MPI_DIMENSIONING'
          STOP
       END IF
    END IF

    ! Identifie la position de chaque processeur dans
    ! le domaine global
    ALLOCATE(KGABS(0:NB_PROCSIN-1))
    ALLOCATE(KGORD(0:NB_PROCSIN-1))
    DO IY = 0, (GYPT/LYPT)-1
       DO IX = 0, (GXPT/LXPT)-1
          KGABS(IX+IY*(GXPT/LXPT)) = IX*LXPT
          KGORD(IX+IY*(GXPT/LXPT)) = IY*LYPT
       END DO
    END DO

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
    IQ = RANGIN/(GXPT/LXPT)
    KEIGHOUR(1) = MOD(KRANG+GXPT/LXPT,NBPROCS)
    KEIGHOUR(4) = MOD(KRANG-GXPT/LXPT+NBPROCS,NBPROCS)
    KEIGHOUR(2) = MOD(KRANG-1+GXPT/LXPT,GXPT/LXPT)+IQ*(GXPT/LXPT)
    KEIGHOUR(3) = MOD(KRANG+1,GXPT/LXPT)+IQ*(GXPT/LXPT)

  END SUBROUTINE DEF_MPI_DIMENSION
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
    REAL(8), PARAMETER :: B=4000.d0

    IF (C_OROTYP .EQ. "A") THEN
       POUT = RHMAX/( ONE + ((PX/RXA)**2) + ((PY/RYA)**2) )
    ELSE IF (C_OROTYP .EQ. "S") THEN
       POUT = RHMAX*DEXP(-(PX/RXA)**2)*DCOS(RPI*PX/B)**2
    ELSE IF (C_OROTYP .EQ. "G") THEN
       POUT = RHMAX*DEXP(-(PX/RXA)**2)*DEXP(-(PY/RYA)**2)   
    ELSE 
       POUT = 0.d0
    END IF
  END FUNCTION ROROG
   !**********************************************************************!
  FUNCTION DOROGDX(PX,PY) RESULT(POUT)
    REAL(8)            :: PX,PY, POUT
    REAL(8), PARAMETER :: B=4000.d0

    IF (C_OROTYP .EQ. "A") THEN
       POUT = -2.d0*RHMAX*(PX/(RXA**2)) &
            & /((1.d0+((PX/RXA)**2)+((PY/RYA)**2))**2)
    ELSE IF (C_OROTYP .EQ. "S") THEN
       POUT = -2.d0*RHMAX*DEXP(-(PX/RXA)**2)*DCOS(RPI*PX/B)*((RPI/B)*DSIN(RPI*PX/B) & 
            & +(PX/RXA**2)*DCOS(RPI*PX/B))
    ELSE IF (C_OROTYP .EQ. "G") THEN
       POUT = -2.d0*PX*(RHMAX/(RXA**2))*DEXP(-(PX/RXA)**2)*DEXP(-(PY/RYA)**2)
    ELSE 
      POUT = 0.d0
    END IF

  END FUNCTION DOROGDX
  !*************************************************************************
  FUNCTION DOROGDY(PX,PY) RESULT(POUT)
    REAL(8)            :: PX,PY, POUT
    REAL(8), PARAMETER :: B=4000.d0

    IF (C_OROTYP .EQ. "A") THEN
       POUT = -2.d0*RHMAX*(PY/(RYA**2)) &
            & /((1.d0+((PX/RXA)**2)+((PY/RYA)**2))**2)
    ELSE IF (C_OROTYP .EQ. "S") THEN
       POUT = ZERO
    ELSE IF (C_OROTYP .EQ. "G") THEN
       POUT = -2.d0*PY*(RHMAX/(RYA**2))*DEXP(-(PX/RXA)**2)*DEXP(-(PY/RYA)**2)   
    ELSE 
      POUT = 0.d0
    END IF

  END FUNCTION DOROGDY
  !************************************************************************!
  FUNCTION D2OROGDX2(PX,PY) RESULT(POUT)
    REAL(8)            :: PX,PY, POUT
    REAL(8), PARAMETER :: B=4000.d0

    IF (C_OROTYP .EQ. "A") THEN
       POUT = -(2.d0*RHMAX/(RXA**2))*(1.d0-3.d0*((PX/RXA)**2)+((PY/RYA)**2)) &
            & /((1.d0+((PX/RXA)**2)+((PY/RYA)**2))**3)
    ELSE IF (C_OROTYP .EQ. "S") THEN
       POUT = 2.d0*RHMAX*DEXP(-(PX/RXA)**2)*(((2*PX/RXA**2)*DCOS(RPI*PX/B) & 
            & + (RPI/B)*DSIN(RPI*PX/B))*((RPI/B)*DSIN(RPI*PX/B) & 
            & +(PX/RXA**2)*DCOS(RPI*PX/B)) &
            & + DCOS(RPI*PX/B)*((PX*RPI/(B*RXA**2))*DSIN(RPI*PX/B) & 
            & -((RPI/B)**2+1.d0/RXA**2)*DCOS(RPI*PX/B)))
    ELSE IF (C_OROTYP .EQ. "G") THEN
       POUT = -2.d0*(RHMAX/(RXA**2))*(1d0-2.d0*((PX/RXA)**2)) &
            & *DEXP(-(PX/RXA)**2)*DEXP(-(PY/RYA)**2)      
    ELSE 
      POUT = 0.d0
    END IF

  END FUNCTION D2OROGDX2
  !*******************************************************
  FUNCTION D2OROGDY2(PX,PY) RESULT(POUT)
    REAL(8)            :: PX,PY, POUT
    REAL(8), PARAMETER :: B=4000.d0

    IF (C_OROTYP .EQ. "A") THEN
       POUT = -(2.d0*RHMAX/(RYA**2))*(1.d0-3.d0*((PY/RYA)**2)+((PX/RXA)**2)) &
            & /((1.d0+((PX/RXA)**2)+((PY/RYA)**2))**3)
    ELSE IF (C_OROTYP .EQ. "S") THEN
       POUT = ZERO
    ELSE IF (C_OROTYP .EQ. "G") THEN
       POUT = -2.d0*(RHMAX/(RYA**2))*(1d0-2.d0*((PY/RYA)**2)) &
            & *DEXP(-(PX/RXA)**2)*DEXP(-(PY/RYA)**2) 
    ELSE 
      POUT = 0.d0
    END IF

  END FUNCTION D2OROGDY2
  !***************************************
  FUNCTION D2OROGDXY(PX,PY) RESULT(POUT)
    REAL(8)            :: PX,PY, POUT
    REAL(8), PARAMETER :: B=4000.d0

    IF (C_OROTYP .EQ. "A") THEN
       POUT = 8.d0*RHMAX*(PX/(RXA**2))*(PY/(RYA**2)) &
            & /((1.d0+((PX/RXA)**2)+((PY/RYA)**2))**3)
    ELSE IF (C_OROTYP .EQ. "S") THEN
       POUT = ZERO
    ELSE IF (C_OROTYP .EQ. "G") THEN
       POUT = 4.d0*(PY/(RYA**2))*(PX/(RXA**2))*RHMAX*DEXP(-(PX/RXA)**2)*DEXP(-(PY/RYA)**2)  
    ELSE 
      POUT = 0.d0
    END IF

  END FUNCTION D2OROGDXY
  !***************************************************************************************!
  
  !#####################################################################
  !#####################################################################
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
  
  !*********************************************************************!
  !*********************************************************************!
  !*********************************************************************!

END MODULE MOD_PARAMETER

!=====================================================!
!=====================================================!
!=====================================================!
