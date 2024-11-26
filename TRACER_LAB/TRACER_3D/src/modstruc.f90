!##################   PARAMETRES   ###################!
!#                                                   #!
!# auteur : F. Voitus                                #!
!# sujet  : Définit les parametres du modèle         #!
!#                                                   #!
!#####################################################!

!=====================================================!
!================  MODULE PARAMETRES  ================!
!=====================================================!

MODULE MOD_STRUCTURE


  USE MOD_PARAMETER, ONLY :  NLEV,GXPT,GYPT,GEXT,ZERO,ONE,TWO,RDX,RDY,RDZ,RLX,RLY,     & 
                          &  RN,RG,RCP,RR,RKAPPA,RTSUR,RPSUR,RPI,RT00,RP00,RU00,RV00,  &
                          &  RY0,RX0,RZ0,RXS,RYS,RZS,RTPP,RQPP,RXC,RYC,RZC,RXA,RYA,RZA
                          
  IMPLICIT NONE
  
  ! type contenant l'ensemble des variables
  TYPE VARPROG
     REAL(8), DIMENSION(:,:), ALLOCATABLE   :: U,V,ETADOTF
     REAL(8), DIMENSION(:,:), ALLOCATABLE   :: VH,UH,ETADOTH
     REAL(8), DIMENSION(:,:), ALLOCATABLE   :: T,MF
     REAL(8), DIMENSION(:),   ALLOCATABLE   :: PIS
     REAL(8), DIMENSION(:,:), ALLOCATABLE   :: Q,MH     
     REAL(8), DIMENSION(:,:), ALLOCATABLE   :: DIV
  END TYPE VARPROG
  
  ! type contenant l'ensemble des variables dynamique
  TYPE VARGLOB
     REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: T,U,V
     REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: DIV
     REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: Q
     REAL(8), DIMENSION(:,:),   ALLOCATABLE :: PIS
  END TYPE VARGLOB

  ! type 
  TYPE RHS_LOC
     REAL(8), DIMENSION(:,:), ALLOCATABLE   :: T
     REAL(8), DIMENSION(:,:), ALLOCATABLE   :: Q
     REAL(8), DIMENSION(:),   ALLOCATABLE   :: PIS
  END TYPE RHS_LOC

  ! type  
  TYPE RHS_GLB
     REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: T
     REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: Q
     REAL(8), DIMENSION(:,:),   ALLOCATABLE :: PIS
  END TYPE RHS_GLB
  
  ! type contenant l'ensemble des variables des couplages
  TYPE BACKGROUND
     REAL(8), DIMENSION(NLEV)               :: U,V,T
     REAL(8)                                :: PIS
  END TYPE BACKGROUND

  ! type contenant l'ensemble des vecteurs servant aux calculs des operateurs verticaux
  TYPE GEOMETRY
     REAL(8), DIMENSION(0:NLEV)             :: AH,BH
     REAL(8), DIMENSION(0:NLEV)             :: ETAH,DETAH        
     REAL(8), DIMENSION(NLEV)               :: ETAF,DETAF
     REAL(8), DIMENSION(NLEV)               :: DELTB
     REAL(8), DIMENSION(:),     ALLOCATABLE :: PIS               ! \pi_{S}*
     REAL(8), DIMENSION(:),     ALLOCATABLE :: ZS,DZSDX,DZSDY
     REAL(8), DIMENSION(:,:),   ALLOCATABLE :: DELTB
     REAL(8), DIMENSION(:,:),   ALLOCATABLE :: PIH               ! \pi sur demi-niveaux 
     REAL(8), DIMENSION(:,:),   ALLOCATABLE :: PIF               ! \pi sur niveaux pleins 
     REAL(8), DIMENSION(:,:),   ALLOCATABLE :: DELTAPI           ! m
     REAL(8), DIMENSION(:,:),   ALLOCATABLE :: ALPHA,DELTA       ! \alpha et \delta
     REAL(8), DIMENSION(:,:),   ALLOCATABLE :: EPS               ! vecteur sur les demi-niveau
  END TYPE GEOMETRY

   ! type contenant l'ensemble des variables des couplages
  TYPE DOMAIN
     REAL(8), DIMENSION(NLEV)               :: ZF,ETAF
     REAL(8), DIMENSION(0:NLEV)             :: ZH,ETAH
     REAL(8), DIMENSION(0:NLEV)             :: AH,BH
     REAL(8), DIMENSION(:,:),   ALLOCATABLE :: YX,XY
     REAL(8), DIMENSION(:,:),   ALLOCATABLE :: XZ,ZX
     REAL(8), DIMENSION(:,:),   ALLOCATABLE :: ZS,DZSDX,DZSDY
     REAL(8), DIMENSION(:,:),   ALLOCATABLE :: PIS_REF
     REAL(8), DIMENSION(:,:),   ALLOCATABLE :: DPISDX_REF,DPISDY_REF
  END TYPE DOMAIN

  TYPE SMILAG_STRUCT
     REAL(8),    DIMENSION(:,:,:),   ALLOCATABLE  :: XDOT,    YDOT,   ZDOT
     REAL(8),    DIMENSION(:,:,:),   ALLOCATABLE  :: XDOTF,   YDOTF,  ZDOTF
     REAL(8),    DIMENSION(:,:,:),   ALLOCATABLE  :: XTRAJ,   YTRAJ,  ZTRAJ
     INTEGER(8), DIMENSION(:,:,:,:), ALLOCATABLE  :: NXLAG,   NYLAG,  NZLAG
     REAL(8),    DIMENSION(:,:,:,:), ALLOCATABLE  :: XWEI,    YWEI,   ZWEI
     INTEGER(8), DIMENSION(:,:,:,:), ALLOCATABLE  :: NXLAGH,  NYLAGH, NZLAGH
     REAL(8),    DIMENSION(:,:,:,:), ALLOCATABLE  :: XWEIH,   YWEIH,  ZWEIH
     INTEGER(8), DIMENSION(:,:,:,:), ALLOCATABLE  :: NXLAGH,  NYLAGH, NZLAGH
     REAL(8),    DIMENSION(:,:,:,:), ALLOCATABLE  :: XWEIH,   YWEIH,  ZWEIH
     INTEGER(8), DIMENSION(:,:,:,:), ALLOCATABLE  :: NXLAGS,  NYLAGS
     REAL(8),    DIMENSION(:,:,:,:), ALLOCATABLE  :: XWEIS,   YWEIS
  END TYPE SMILAG_STRUCT

  ! type mpdata 
  TYPE MPDATA_STRUCT
     REAL(8),    DIMENSION(:,:,:,:), ALLOCATABLE  :: MF,    MH
     REAL(8),    DIMENSION(:,:,:),   ALLOCATABLE  :: XDOT,  XDOTF,  XDOTH
     REAL(8),    DIMENSION(:,:,:),   ALLOCATABLE  :: YDOT,  YDOTF,  YDOTH
     REAL(8),    DIMENSION(:,:,:,:), ALLOCATABLE  :: ZDOT,  ZDOTF,  ZDOTH
     REAL(8),    DIMENSION(:,:),     ALLOCATABLE  :: XDOTS, YDOTS
  END TYPE MPDATA_STRUCT
  
CONTAINS

  !*****************************************************!
  !***************  LISTE DES ROUTINES  ****************!
  !*****************************************************!

  !***************************************************************************************
  !***************************************************************************************
  !***************************************************************************************
  SUBROUTINE DEF_DOMAIN_GEOMETRY(DOM,PAH,PBH)

    USE MOD_PARAMETER, ONLY : ROROG,DOROGDX,DOROGDY
    USE MOD_FFT,       ONLY : FFTDIR_2D,FFTORO_2D
    
    IMPLICIT NONE

    TYPE(GEODOM),               INTENT(INOUT) :: DOM
    REAL(8), DIMENSION(0:NLEV), INTENT(IN)    :: PAH,PBH
    
    INTEGER(8)                                :: I,J,L
    REAL(8)                                   :: ZZ,ZX,ZY,ZRELAX,ZPERIO,ZUNIT
    REAL(8)                                   :: ZORODX(GXPT,GYPT),ZORODY(GXPT,GYPT)

    ZUNIT = 1000.d0
    
    DOM%AH(:) = PAH(:)
    DOM%BH(:) = PBH(:)
    
    ! Définition de la géométrie du problème 
    DOM%ZH(0)      = REAL(NLEV,8)*RDZ
    DO L = 1, NLEV
       DOM%ZF(L)   = (REAL(NLEV-L,8)+0.5d0)*RDZ
       DOM%ZH(L)   = REAL(NLEV-L,8)*RDZ
    END DO
    
    DO I = 1, GXPT
       DOM%XY(I,:) = (REAL(I-1,8)*RDX - RXC)/ZUNIT     
       DOM%XZ(:,I) = (REAL(I-1,8)*RDX - RXC)/ZUNIT 
    END DO
    DO J = 1, GYPT 
       DOM%YX(:,J) = (REAL(J-1,8)*RDY - RYC)/ZUNIT       
    END DO
    
    ! orography definition
    DO I = 1, GXPT
       ZX=REAL(I-1,8)*RDX - RXC
       DO J = 1,GYPT
          ZY=REAL(J-1,8)*RDY - RYC
          DOM%ZS(I,J) = ROROG(ZX,ZY)
          ZORODX(I,J) = DOROGDX(ZX,ZY)
          ZORODY(I,J) = DOROGDY(ZX,ZY)
       END DO   
    END DO
    
    CALL FFTDIR_2D(DOM%ZS)
    CALL FFTORO_2D(DOM%ZS,DOM%DZSDX,DOM%DZSDY,LALIAZ=LRALIAZ)

    DO I = 1, GXPT
       DO J = 1,GYPT
          DOM%PIS_REF(I,J)    =  RP00*DEXP(-(RG/(RR*RT00))*DOM%ZS(I,J))
          DOM%DPISDX_REF(I,J) = -(RG/(RR*RT00))*DOM%PIS_REF(I,J)*DOM%DZSDX(I,J)
          DOM%DPISDY_REF(I,J) = -(RG/(RR*RT00))*DOM%PIS_REF(I,J)*DOM%DZSDY(I,J) 
       END DO   
    END DO
    
    DO I=1,GXPT
       DOM%ZX(NLEV,I) = DOM%ZS(I,INT(RYS/RDY)+1)
       DO L = 0, NLEV-1
          DOM%ZX(L,:) = (REAL(NLEV-L,8))*RDZ
       END DO
    END DO

    DOM%ETAH(0)       = DOM%BH(0)+(ONE/RP00)*DOM%AH(0)                      
    DO J=1,NLEV
       DOM%ETAH(J)    = DOM%BH(J)+(ONE/RP00)*DOM%AH(J)
       DOM%ETAF(J)    = HALF*((DOM%BH(J)+DOM%BH(J-1))  &
                      & +(ONE/RP00)*(DOM%AH(J)+DOM%AH(J-1)))
    END DO

    WRITE(*,*) " "
    WRITE(*,*) 'MAX OROGRAPHY SLOPE (%) : ', &
    & 100*MAX(MAXVAL(DOM%DZSDX),MAXVAL(DOM%DZSDY))
    WRITE(*,*) 'MAX OROGRAPHY SLOPE (°) : ', &
    & (180/RPI)*ATAN(MAX(MAXVAL(DOM%DZSDX),MAXVAL(DOM%DZSDY)))
    
  END SUBROUTINE DEF_DOMAIN_GEOMETRY
  !***************************************************************************************
  !***************************************************************************************
  !***************************************************************************************
  !***************************************************************************************
  SUBROUTINE DEF_BACKGROUND_STATE(STB)

    USE MOD_PARAMETER, ONLY : RU00,RV00,RTSUR,RPSUR,RZ1,RZ2
    
    IMPLICIT NONE

    TYPE(BACKGROUND),          INTENT(OUT)   :: STB
    INTEGER                                  :: J
    REAL(8)                                  :: ZZF,ZIJ
    
    STB%T(:)       = RTSUR
    STB%PIS        = RPSUR

    IF (NCASE == 1) THEN
       DO J = 1,NLEV
          ZZF      = (REAL(NLEV-J,8)+0.5d0)*RDZ
          ZIJ      = MAX(ZERO,MIN((ZZF-RZ1)/(RZ2-RZ1),ONE))   
          STB%U(J) = RU00*(DSIN((RPI/TWO)*ZIJ)**2)
          STB%V(J) = RV00
       END DO
    ELSE
       STB%U(:)    = RU00
       STB%V(:)    = RV00
    END IF
    
  END SUBROUTINE DEF_BACKGROUND_STATE
  !***************************************************************************************
  !***************************************************************************************
  !***************************************************************************************
  !***************************************************************************************
  SUBROUTINE DEF_INITIAL_PERTURBATION(STG,STB,DOM)
    
    USE MOD_FFT,       ONLY : FFTDIR_UV,FFTDIV
    
    IMPLICIT NONE

    TYPE(VARGLOB),                        INTENT(INOUT) :: STG
    TYPE(BACKGROUND),                     INTENT(IN)    :: STB
    TYPE(GEODOM),                         INTENT(IN)    :: DOM
        
    INTEGER                                             :: JL,NX,NY
    REAL(8)                                             :: ZDISC,ZX,ZY,ZZ
         
    DO NX = 1,GXPT
       ZX = REAL(NX-1,8)*RDX
       DO NY = 1,GYPT
          ZY = REAL(NY-1,8)*RDY
          STG%PIS(NX,NY)  = RPSUR*DEXP(-RPSTA*(RG/(RR*RTSUR))*DOM%ZS(NX,NY))
          DO JL=1,NLEV
             ZZ           = (REAL(NLEV-J,8)+0.5d0)*RDZ
             ZDISC        = MIN(DSQRT((((ZX-RX0)/RAX)**2)                    &
                          & +(((ZY-RY0)/RAX)**2)+(((ZZ-RZ0)/RAZ)**2)),ONE)
             STG%U(JL,NX,NY) = STB%U(JL)*(DOM%ETAH(JL)-DOM%ETAH(JL-1))       &
                          & /( ((DOM%AH(JL)-DOM%AH(JL-1))/RP00)              &
                          & +(RS*(STG%PIS(NX,NY)/RP00)+(ONE-RS))*(DOM%BH(JL) &
                          & - DOM%BH(JL-1)) )
             STG%V(JL,NX,NY) = STB%V(JL)
             STG%T(JL,NX,NY) = STB%T(JL)+RTPP*((DCOS((RPI/TWO)*ZDISC))**2)
          END DO
          DO JL=0,NLEV
             ZZ           = REAL(NLEV-J,8)*RDZ
             ZDISC        = MIN(DSQRT((((ZX-RX0)/RAX)**2)                    &
                          & +(((ZY-RY0)/RAX)**2)+(((ZZ-RZ0)/RAZ)**2)),ONE)
             STG%Q(JL,NX,NY) = RQPP*((DCOS((RPI/TWO)*ZDISC))**2)
          END DO
       END DO   
    END DO

    DO JL=1,NLEV
       CALL FFTDIR_UV(STG%U(JL,:,:),STG%V(JL,:,:))
       CALL FFTDIV(STG%U(JL,:,:),STG%V(JL,:,:),STG%DIV(JL,:,:),LALIAZ=LRALIAZ) 
    END DO 
        
  END SUBROUTINE DEF_INITIAL_PERTURBATION
  !*******************************************************************************!
  !*******************************************************************************!
  !*******************************************************************************!
  SUBROUTINE DEF_AB_FUNCTIONS(PAH,PBH,PETAF,PETAH,PDETAF,PDETAH,PDELB,PEPS,C_OPT)

    IMPLICIT NONE

    REAL(8), DIMENSION(0:NLEV),         INTENT(OUT) :: PAH,PBH
    REAL(8), DIMENSION(0:NLEV),         INTENT(OUT) :: PETAH,PDETAH,PEPS
    REAL(8), DIMENSION(NLEV),           INTENT(OUT) :: PETAF,PDETAF,PDELB
    CHARACTER(LEN=*),                   INTENT(IN)  :: C_OPT

    INTEGER                                         :: I,J
    REAL(8)                                         :: ZEPS,ZDZDH
    REAL(8), DIMENSION(0:90)                        :: ZA1,ZB1
    REAL(8), DIMENSION(0:137)                       :: ZA2,ZB2

    ZEPS = 100*ABS(EPSILON(ONE))
    
    PAH(:)    = ZERO
    PBH(:)    = ZERO

    IF (C_OPT == 'Z') THEN
       PBH(NLEV) = ONE
       DO J = 1,NLEV-1
          ZDZDH       = MAX((RG/RR)*(RDZ/TWO)/RTSUR,ZEPS)
          PBH(NLEV-J) = PBH(NLEV-J+1)*((DSQRT((ZDZDH**2)+ONE)-ZDZDH)**2)
       END DO
    ELSE IF ((C_OPT == 'S').AND.(NLEV == 90)) THEN
       ZA1(0) = 0.d0
       ZA1(1) = 2718.2817
       ZA1(2) = 3747.0893
       ZA1(3) = 4915.5871
       ZA1(4) = 6184.6290
       ZA1(5) = 7511.8256
       ZA1(6) = 8869.0999
       ZA1(7) =10030.3026
       ZA1(8) =10981.6510
       ZA1(9) =11785.4396
       ZA1(10)=12467.1137
       ZA1(11)=13036.9950
       ZA1(12)=13504.5335
       ZA1(13)=13878.6376
       ZA1(14)=14167.4665
       ZA1(15)=14378.1637
       ZA1(16)=14516.6238
       ZA1(17)=14587.3575
       ZA1(18)=14595.7546
       ZA1(19)=14546.6959
       ZA1(20)=14444.6128
       ZA1(21)=14293.6727
       ZA1(22)=14097.9871
       ZA1(23)=13861.8303
       ZA1(24)=13589.8532
       ZA1(25)=13287.2816
       ZA1(26)=12960.0840
       ZA1(27)=12612.2519
       ZA1(28)=12247.2879
       ZA1(29)=11868.2461
       ZA1(30)=11477.7792
       ZA1(31)=11078.1891
       ZA1(32)=10671.4786
       ZA1(33)=10259.4024
       ZA1(34)= 9843.5158
       ZA1(35)= 9425.2210
       ZA1(36)= 9005.8078
       ZA1(37)= 8586.4919
       ZA1(38)= 8168.6961
       ZA1(39)= 7753.7612
       ZA1(40)= 7342.9551
       ZA1(41)= 6937.4811
       ZA1(42)= 6538.4835
       ZA1(43)= 6147.0529
       ZA1(44)= 5764.2288
       ZA1(45)= 5391.0012
       ZA1(46)= 5028.3104
       ZA1(47)= 4677.0451
       ZA1(48)= 4338.0394
       ZA1(49)= 4012.0686
       ZA1(50)= 3699.8436
       ZA1(51)= 3402.0051
       ZA1(52)= 3119.0079
       ZA1(53)= 2851.1415
       ZA1(54)= 2598.5379
       ZA1(55)= 2361.1824
       ZA1(56)= 2138.9263
       ZA1(57)= 1931.5000
       ZA1(58)= 1738.5284
       ZA1(59)= 1559.5458
       ZA1(60)= 1394.0114
       ZA1(61)= 1241.3250   
       ZA1(62)= 1100.8420
       ZA1(63)=  971.8880
       ZA1(64)=  853.7728
       ZA1(65)=  746.1066
       ZA1(66)=  648.4464
       ZA1(67)=  560.3053  
       ZA1(68)=  481.1623
       ZA1(69)=  410.4714
       ZA1(70)=  347.6705
       ZA1(71)=  292.1898
       ZA1(72)=  243.4601
       ZA1(73)=  200.9194
       ZA1(74)=  164.0199
       ZA1(75)=  132.2333
       ZA1(76)=  105.0561
       ZA1(77)=   82.0135
       ZA1(78)=   62.7448
       ZA1(79)=   46.8857
       ZA1(80)=   34.0715
       ZA1(81)=   23.9407
       ZA1(82)=   16.1385
       ZA1(83)=   10.3203
       ZA1(84)=    6.1553
       ZA1(85)=    3.3306
       ZA1(86)=    1.5543
       ZA1(87)=    0.5601
       ZA1(88)=    0.1103
       ZA1(89)=    0.d0
       ZA1(90)=    0.d0

       ZB1(0) =0.d0
       ZB1(1) =0.d0
       ZB1(2) =0.d0
       ZB1(3) =0.d0
       ZB1(4) =0.d0
       ZB1(5) =0.d0
       ZB1(6) =0.d0
       ZB1(7) =0.00208087
       ZB1(8) =0.00633423
       ZB1(9) =0.01217321
       ZB1(10)=0.01939861
       ZB1(11)=0.02783712
       ZB1(12)=0.03735295
       ZB1(13)=0.04784293
       ZB1(14)=0.05923520
       ZB1(15)=0.07148733
       ZB1(16)=0.08457464
       ZB1(17)=0.09837415
       ZB1(18)=0.11283833
       ZB1(19)=0.12793151
       ZB1(20)=0.14361851
       ZB1(21)=0.15985947
       ZB1(22)=0.17660445
       ZB1(23)=0.19378795
       ZB1(24)=0.21132357
       ZB1(25)=0.22909875
       ZB1(26)=0.24697007
       ZB1(27)=0.26490217
       ZB1(28)=0.28286756
       ZB1(29)=0.30084563
       ZB1(30)=0.31882157
       ZB1(31)=0.33678539
       ZB1(32)=0.35473076
       ZB1(33)=0.37265407
       ZB1(34)=0.39055337
       ZB1(35)=0.40842741
       ZB1(36)=0.42627480
       ZB1(37)=0.44409309
       ZB1(38)=0.46186741
       ZB1(39)=0.47958299
       ZB1(40)=0.49722495
       ZB1(41)=0.51477798
       ZB1(42)=0.53222617
       ZB1(43)=0.54955275
       ZB1(44)=0.56673991
       ZB1(45)=0.58376861
       ZB1(46)=0.60061846
       ZB1(47)=0.61726751
       ZB1(48)=0.63369222
       ZB1(49)=0.64986725
       ZB1(50)=0.66576547
       ZB1(51)=0.68135784
       ZB1(52)=0.69661932
       ZB1(53)=0.71152830
       ZB1(54)=0.72606668
       ZB1(55)=0.74021997
       ZB1(56)=0.75397731
       ZB1(57)=0.76733145
       ZB1(58)=0.78027882
       ZB1(59)=0.79281944
       ZB1(60)=0.80495693
       ZB1(61)=0.81669843
       ZB1(62)=0.82805459
       ZB1(63)=0.83903946
       ZB1(64)=0.84967049
       ZB1(65)=0.85993862
       ZB1(66)=0.86983693
       ZB1(67)=0.87936065
       ZB1(68)=0.88850710
       ZB1(69)=0.89727571
       ZB1(70)=0.90566797
       ZB1(71)=0.91368743
       ZB1(72)=0.92133963
       ZB1(73)=0.92863212
       ZB1(74)=0.93557441
       ZB1(75)=0.94217792
       ZB1(76)=0.94845599
       ZB1(77)=0.95442386
       ZB1(78)=0.96007297
       ZB1(79)=0.96539491
       ZB1(80)=0.97038132
       ZB1(81)=0.97502394
       ZB1(82)=0.97931457
       ZB1(83)=0.98324507
       ZB1(84)=0.98680735
       ZB1(85)=0.98999338
       ZB1(86)=0.99279512
       ZB1(87)=0.99520458
       ZB1(88)=0.99721377
       ZB1(89)=0.99881471
       ZB1(90)=1.d0

       DO J=1,NLEV
          PAH(J) = ZA1(J)
          PBH(J) = ZB1(J)
       END DO
       
    ELSE IF ((C_OPT == 'S').AND.(NLEV == 137)) THEN
       ZA2(0)   = 0.000000
       ZA2(1)   = 2.000365
       ZA2(2)   = 3.102241
       ZA2(3)   = 4.666084
       ZA2(4)   = 6.827977
       ZA2(5)   = 9.746966
       ZA2(6)   = 13.605424
       ZA2(7)   = 18.608931
       ZA2(8)   = 24.985718
       ZA2(9)   = 32.985710
       ZA2(10)  = 42.879242
       ZA2(11)  = 54.955463
       ZA2(12)  = 69.520576
       ZA2(13)  = 86.895882
       ZA2(14)  = 107.415741
       ZA2(15)  = 131.425507
       ZA2(16)  = 159.279404
       ZA2(17)  = 191.338562
       ZA2(18)  = 227.968948
       ZA2(19)  = 269.539581
       ZA2(20)  = 316.420746
       ZA2(21)  = 368.982361
       ZA2(22)  = 427.592499
       ZA2(23)  = 492.616028
       ZA2(24)  = 564.413452
       ZA2(25)  = 643.339905
       ZA2(26)  = 729.744141
       ZA2(27)  = 823.967834
       ZA2(28)  = 926.344910
       ZA2(29)  = 1037.201172
       ZA2(30)  = 1156.853638
       ZA2(31)  = 1285.610352
       ZA2(32)  = 1423.770142
       ZA2(33)  = 1571.622925
       ZA2(34)  = 1729.448975
       ZA2(35)  = 1897.519287
       ZA2(36)  = 2076.095947
       ZA2(37)  = 2265.431641
       ZA2(38)  = 2465.770508
       ZA2(39)  = 2677.348145
       ZA2(40)  = 2900.391357
       ZA2(41)  = 3135.119385
       ZA2(42)  = 3381.743652
       ZA2(43)  = 3640.468262
       ZA2(44)  = 3911.490479
       ZA2(45)  = 4194.930664
       ZA2(46)  = 4490.817383
       ZA2(47)  = 4799.149414
       ZA2(48)  = 5119.895020
       ZA2(49)  = 5452.990723
       ZA2(50)  = 5798.344727
       ZA2(51)  = 6156.074219
       ZA2(52)  = 6526.946777
       ZA2(53)  = 6911.870605
       ZA2(54)  = 7311.869141
       ZA2(55)  = 7727.412109
       ZA2(56)  = 8159.354004
       ZA2(57)  = 8608.525391
       ZA2(58)  = 9076.400391
       ZA2(59)  = 9562.682617
       ZA2(60)  = 10065.978516
       ZA2(61)  = 10584.631836
       ZA2(62)  = 11116.662109
       ZA2(63)  = 11660.067383
       ZA2(64)  = 12211.547852
       ZA2(65)  = 12766.873047
       ZA2(66)  = 13324.668945
       ZA2(67)  = 13881.331055
       ZA2(68)  = 14432.139648
       ZA2(69)  = 14975.615234
       ZA2(70)  = 15508.256836
       ZA2(71)  = 16026.115234
       ZA2(72)  = 16527.322266
       ZA2(73)  = 17008.789063
       ZA2(74)  = 17467.613281
       ZA2(75)  = 17901.621094
       ZA2(76)  = 18308.433594
       ZA2(77)  = 18685.718750
       ZA2(78)  = 19031.289063
       ZA2(79)  = 19343.511719
       ZA2(80)  = 19620.042969
       ZA2(81)  = 19859.390625
       ZA2(82)  = 20059.931641
       ZA2(83)  = 20219.664063
       ZA2(84)  = 20337.863281
       ZA2(85)  = 20412.308594
       ZA2(86)  = 20442.078125
       ZA2(87)  = 20425.718750
       ZA2(88)  = 20361.816406
       ZA2(89)  = 20249.511719
       ZA2(90)  = 20087.085938
       ZA2(91)  = 19874.025391
       ZA2(92)  = 19608.572266
       ZA2(93)  = 19290.226563
       ZA2(94)  = 18917.460938
       ZA2(95)  = 18489.707031
       ZA2(96)  = 18006.925781
       ZA2(97)  = 17471.839844
       ZA2(98)  = 16888.687500
       ZA2(99)  = 16262.046875
       ZA2(100) = 15596.695313
       ZA2(101) = 14898.453125
       ZA2(102) = 14173.324219
       ZA2(103) = 13427.769531
       ZA2(104) = 12668.257813
       ZA2(105) = 11901.339844
       ZA2(106) = 11133.304688
       ZA2(107) = 10370.175781
       ZA2(108) = 9617.515625
       ZA2(109) = 8880.453125
       ZA2(110) = 8163.375000
       ZA2(111) = 7470.343750
       ZA2(112) = 6804.421875
       ZA2(113) = 6168.531250
       ZA2(114) = 5564.382813
       ZA2(115) = 4993.796875
       ZA2(116) = 4457.375000
       ZA2(117) = 3955.960938
       ZA2(118) = 3489.234375
       ZA2(119) = 3057.265625
       ZA2(120) = 2659.140625
       ZA2(121) = 2294.242188
       ZA2(122) = 1961.500000
       ZA2(123) = 1659.476563
       ZA2(124) = 1387.546875
       ZA2(125) = 1143.250000
       ZA2(126) = 926.507813
       ZA2(127) = 734.992188
       ZA2(128) = 568.062500
       ZA2(129) = 424.414063
       ZA2(130) = 302.476563
       ZA2(131) = 202.484375
       ZA2(132) = 122.101563
       ZA2(133) = 62.781250
       ZA2(134) = 22.835938
       ZA2(135) = 3.757813
       ZA2(136) = 0.000000
       ZA2(137) = 0.000000

       ZB2(0)   = 0.0000000000
       ZB2(1)   = 0.0000000000
       ZB2(2)   = 0.0000000000
       ZB2(3)   = 0.0000000000
       ZB2(4)   = 0.0000000000
       ZB2(5)   = 0.0000000000
       ZB2(6)   = 0.0000000000
       ZB2(7)   = 0.0000000000
       ZB2(8)   = 0.0000000000
       ZB2(9)   = 0.0000000000
       ZB2(10)  = 0.0000000000
       ZB2(11)  = 0.0000000000
       ZB2(12)  = 0.0000000000
       ZB2(13)  = 0.0000000000
       ZB2(14)  = 0.0000000000
       ZB2(15)  = 0.0000000000
       ZB2(16)  = 0.0000000000
       ZB2(17)  = 0.0000000000
       ZB2(18)  = 0.0000000000
       ZB2(19)  = 0.0000000000
       ZB2(20)  = 0.0000000000
       ZB2(21)  = 0.0000000000
       ZB2(22)  = 0.0000000000
       ZB2(23)  = 0.0000000000
       ZB2(24)  = 0.0000000000
       ZB2(25)  = 0.0000000000
       ZB2(26)  = 0.0000000000
       ZB2(27)  = 0.0000000000
       ZB2(28)  = 0.0000000000
       ZB2(29)  = 0.0000000000
       ZB2(30)  = 0.0000000000
       ZB2(31)  = 0.0000000000
       ZB2(32)  = 0.0000000000
       ZB2(33)  = 0.0000000000
       ZB2(34)  = 0.0000000000
       ZB2(35)  = 0.0000000000
       ZB2(36)  = 0.0000000000
       ZB2(37)  = 0.0000000000
       ZB2(38)  = 0.0000000000
       ZB2(39)  = 0.0000000000
       ZB2(40)  = 0.0000000000
       ZB2(41)  = 0.0000000000
       ZB2(42)  = 0.0000000000
       ZB2(43)  = 0.0000000000
       ZB2(44)  = 0.0000000000
       ZB2(45)  = 0.0000000000
       ZB2(46)  = 0.0000000000
       ZB2(47)  = 0.0000000000
       ZB2(48)  = 0.0000000000
       ZB2(49)  = 0.0000000000
       ZB2(50)  = 0.0000000000
       ZB2(51)  = 0.0000000000
       ZB2(52)  = 0.0000000000
       ZB2(53)  = 0.0000000000
       ZB2(54)  = 0.0000000382
       ZB2(55)  = 0.0000067607
       ZB2(56)  = 0.0000243480
       ZB2(57)  = 0.0000589220
       ZB2(58)  = 0.0001119143
       ZB2(59)  = 0.0001985774
       ZB2(60)  = 0.0003403797
       ZB2(61)  = 0.0005615553
       ZB2(62)  = 0.0008896979
       ZB2(63)  = 0.0013528055
       ZB2(64)  = 0.0019918380
       ZB2(65)  = 0.0028571242
       ZB2(66)  = 0.0039709536
       ZB2(67)  = 0.0053778146
       ZB2(68)  = 0.0071333768
       ZB2(69)  = 0.0092614600
       ZB2(70)  = 0.0118060224
       ZB2(71)  = 0.0148156285
       ZB2(72)  = 0.0183184519
       ZB2(73)  = 0.0223548450
       ZB2(74)  = 0.0269635208
       ZB2(75)  = 0.0321760960
       ZB2(76)  = 0.0380263999
       ZB2(77)  = 0.0445479602
       ZB2(78)  = 0.0517730154
       ZB2(79)  = 0.0597284138
       ZB2(80)  = 0.0684482530
       ZB2(81)  = 0.0779583082
       ZB2(82)  = 0.0882857367
       ZB2(83)  = 0.0994616672
       ZB2(84)  = 0.1115046516
       ZB2(85)  = 0.1244481280
       ZB2(86)  = 0.1383128911
       ZB2(87)  = 0.1531250328
       ZB2(88)  = 0.1689104140
       ZB2(89)  = 0.1856894493
       ZB2(90)  = 0.2034912109
       ZB2(91)  = 0.2223328650
       ZB2(92)  = 0.2422440052
       ZB2(93)  = 0.2632418871
       ZB2(94)  = 0.2853540182
       ZB2(95)  = 0.3085984588
       ZB2(96)  = 0.3329390883
       ZB2(97)  = 0.3582541943
       ZB2(98)  = 0.3843633235
       ZB2(99)  = 0.4111247659
       ZB2(100) = 0.4383912086
       ZB2(101) = 0.4660032988
       ZB2(102) = 0.4938003123
       ZB2(103) = 0.5216192007
       ZB2(104) = 0.5493011475
       ZB2(105) = 0.5766921639
       ZB2(106) = 0.6036480665
       ZB2(107) = 0.6300358176
       ZB2(108) = 0.6557359695
       ZB2(109) = 0.6806430221
       ZB2(110) = 0.7046689987
       ZB2(111) = 0.7277387381
       ZB2(112) = 0.7497965693
       ZB2(113) = 0.7707975507
       ZB2(114) = 0.7907167673
       ZB2(115) = 0.8095360398
       ZB2(116) = 0.8272560835
       ZB2(117) = 0.8438811302
       ZB2(118) = 0.8594318032
       ZB2(119) = 0.8739292622
       ZB2(120) = 0.8874075413
       ZB2(121) = 0.8999004960
       ZB2(122) = 0.9114481807
       ZB2(123) = 0.9220956564
       ZB2(124) = 0.9318807721
       ZB2(125) = 0.9408595562
       ZB2(126) = 0.9490644336
       ZB2(127) = 0.9565495253
       ZB2(128) = 0.9633517265
       ZB2(129) = 0.9695134163
       ZB2(130) = 0.9750784039
       ZB2(131) = 0.9800716043
       ZB2(132) = 0.9845418930
       ZB2(133) = 0.9884995222
       ZB2(134) = 0.9919840097
       ZB2(135) = 0.9950025082
       ZB2(136) = 0.9976301193
       ZB2(137) = 1.0000000000

       DO J=1,NLEV
          PAH(J) = ZA2(J)
          PBH(J) = ZB2(J)
       END DO
       
    ENDIF

    PETAH(0)      = PBH(0)+(ONE/RP00)*PAH(0)                      
    DO J=1,NLEV
       PETAH(J)   = PBH(J)+(ONE/RP00)*PAH(J)
       PETAF(J)   = HALF*((PBH(J)+PBH(J-1))  &
                     & +(ONE/RP00)*(PAH(J)+PAH(J-1)))
       PDETAF(J)  = (PBH(J)-PBH(J-1))        &
                     & + (ONE/RP00)*(PAH(J)-PAH(J-1))
       PDELTB(J)  = (PBH(J) - PBH(J-1))
    END DO
    
    DO J=1,NLEV-1
       PDETAH(J)  = PETAF(J+1)-PETAF(J)
       PEPS(J)    = ((PAH(J+1)-PAH(J))+RP00*(PBH(J+1)-PBH(J))) &
                  & /((PAH(J+1)-PAH(J-1))+RP00*(PBH(J+1)-PBH(J-1)))
    END DO

    PDETAH(0)     = PETAF(1)-PETAH(0)
    PDETAH(NLEV)  = PETAH(NLEV)-PETAF(NLEV)

    PEPS(0)       = ZERO
    PEPS(NLEV)    = ONE

    
  END SUBROUTINE DEF_AB_FUNCTIONS
  !*******************************************************************************!
  !*******************************************************************************!
  !*******************************************************************************!
  !*******************************************************************************!
  SUBROUTINE DEF_VERTICAL_METRICS(PPIH,PPIF,PDELTP,PDELTA,PALPHA,PAH,PBH,PPIS)

    USE MOD_PARAMETER,  ONLY : RP00

    IMPLICIT NONE

    REAL(8), DIMENSION(NLEV),   INTENT(OUT)  :: PPIF,PDELTP,PALPHA
    REAL(8), DIMENSION(NLEV),   INTENT(OUT)  :: PDELTA
    REAL(8), DIMENSION(0:NLEV), INTENT(OUT)  :: PPIH
    REAL(8), DIMENSION(0:NLEV), INTENT(IN)   :: PAH,PBH
    REAL(8)                   , INTENT(IN)   :: PPIS
   
    INTEGER                                  :: J
    
    PPIH(0) = PAH(0)+PBH(0)*ZPIS
    
    DO J=1,NLEV
       PPIH(J)    = PAH(J)+PBH(J)*ZPIS
       PDELTP(J)  = (PAH(J)-PAH(J-1))+(PBH(J)-PBH(J-1))*PPIS
    END DO

    PALPHA(1)    = ONE
    PDELTA(1)    = ONE + ONE/RKAPPA
    PPIF(1)      = PDELTP(1)/PDELTA(1)
    
    DO J = 2,NLEV
       PPIF(J)   = DSQRT(PPIH(J)*PPIH(J-1))
       PALPHA(J) = ONE-(PPIF(J)/PPIH(J))
       PDELTA(J) = PDELTP(J)/PPIF(J)
    END DO
   

  END SUBROUTINE DEF_VERTICAL_METRICS
  !*******************************************************************************!
  !*******************************************************************************!
  !*******************************************************************************!
  !*******************************************************************************! 
  !*******************************************************************************!
  SUBROUTINE SUALLOC_MPI(X0,XC,M2,M1,M0,L0,LC,GEO,DIM)

    IMPLICIT NONE

    TYPE(VARPROG),                 INTENT(INOUT) :: X0
    TYPE(VARLBC),                  INTENT(INOUT) :: XC
    TYPE(RHS_LOC),                 INTENT(INOUT) :: M2,M1,M0,L0
    TYPE(RHSLBC),                  INTENT(INOUT) :: LC
    TYPE(GEO_NHNL),                INTENT(INOUT) :: GEO
    INTEGER(8),                    INTENT(IN)    :: DIM

    !DIM = TRANCHEX*TRANCHEY  
 
    ALLOCATE(X0%U(NLEV,DIM),X0%V(NLEV,DIM),X0%T(NLEV,DIM),X0%DW(NLEV,DIM),X0%Q(NLEV,DIM),X0%PIS(DIM))
    ALLOCATE(X0%GW(0:NLEV,DIM),X0%DIV(NLEV,DIM),X0%X(NLEV,DIM),X0%Y(0:NLEV,DIM),X0%ETADOT(NLEV,DIM))
    ALLOCATE(X0%DTDX(NLEV,DIM),X0%DQDX(NLEV,DIM),X0%DPISDX(DIM))
    ALLOCATE(X0%DTDY(NLEV,DIM),X0%DQDY(NLEV,DIM),X0%DPISDY(DIM))
   
    ALLOCATE(XC%U(NLEV,DIM),XC%V(NLEV,DIM),XC%T(NLEV,DIM),XC%DW(NLEV,DIM),XC%Q(NLEV,DIM))
    ALLOCATE(XC%DIV(NLEV,DIM),XC%DTDX(NLEV,DIM),XC%DTDY(NLEV,DIM),XC%DQDX(NLEV,DIM),XC%DQDY(NLEV,DIM))
    ALLOCATE(XC%PIS(DIM),XC%DPISDX(DIM),XC%DPISDY(DIM))
    
    ALLOCATE(M0%U(NLEV,DIM),M0%V(NLEV,DIM),M0%GW(0:NLEV,DIM),M0%T(NLEV,DIM))
    ALLOCATE(M0%DW(NLEV,DIM),M0%Q(NLEV,DIM),M0%PI(NLEV,DIM),M0%PIS(DIM))
    ALLOCATE(L0%U(NLEV,DIM),L0%V(NLEV,DIM),L0%GW(0:NLEV,DIM),L0%T(NLEV,DIM))
    ALLOCATE(L0%DW(NLEV,DIM),L0%Q(NLEV,DIM),L0%PI(NLEV,DIM),L0%PIS(DIM))
    ALLOCATE(M1%U(NLEV,DIM),M1%V(NLEV,DIM),M1%GW(0:NLEV,DIM),M1%T(NLEV,DIM))
    ALLOCATE(M1%DW(NLEV,DIM),M1%Q(NLEV,DIM),M1%PI(NLEV,DIM),M1%PIS(DIM))
    ALLOCATE(M2%U(NLEV,DIM),M2%V(NLEV,DIM),M2%GW(0:NLEV,DIM),M2%T(NLEV,DIM))
    ALLOCATE(M2%DW(NLEV,DIM),M2%Q(NLEV,DIM),M2%PI(NLEV,DIM),M2%PIS(DIM))
    ALLOCATE(LC%U(NLEV,DIM),LC%V(NLEV,DIM),LC%T(NLEV,DIM),LC%DW(NLEV,DIM))
    ALLOCATE(LC%Q(NLEV,DIM),LC%PIS(DIM),LC%RELAX(DIM))

    ALLOCATE(GEO%DELTA(NLEV,DIM),GEO%DELTB(NLEV,DIM),GEO%DELTC(NLEV,DIM),GEO%DELTAPI(NLEV,DIM),GEO%ALPHA(NLEV,DIM))
    ALLOCATE(GEO%PIF(NLEV,DIM),GEO%PIH(0:NLEV,DIM),GEO%PIS(DIM),GEO%EPS(0:NLEV,DIM),GEO%XIM(NLEV,DIM),GEO%XIP(NLEV,DIM))
    ALLOCATE(GEO%DDELTADX(NLEV,DIM),GEO%DDELTADY(NLEV,DIM),GEO%DALPHADX(NLEV,DIM),GEO%DALPHADY(NLEV,DIM))
    ALLOCATE(GEO%ZS(DIM),GEO%DZSDX(DIM),GEO%DZSDY(DIM),GEO%D2ZSDX2(DIM),GEO%D2ZSDY2(DIM),GEO%D2ZSDXY(DIM))

    
  END SUBROUTINE SUALLOC_MPI
  !**********************************************************************************
  !**********************************************************************************
  !**********************************************************************************
  SUBROUTINE SUALLOC_GBL(XG,YSL,DOM)

    IMPLICIT NONE

    TYPE(RHS_GLB),                 INTENT(INOUT) :: XG
    TYPE(SL_STRUCT),               INTENT(INOUT) :: YSL
    TYPE(GEODOM),                  INTENT(INOUT) :: DOM
   
    ALLOCATE(DOM%XY(GXPT,GYPT),DOM%YX(GXPT,GYPT),DOM%XZ(NLEV,GXPT),DOM%ZX(0:NLEV,GXPT))
    ALLOCATE(DOM%ZS(GXPT,GYPT),DOM%DZSDX(GXPT,GYPT),DOM%DZSDY(GXPT,GXPT))
    ALLOCATE(DOM%D2ZSDX2(GXPT,GYPT),DOM%D2ZSDY2(GXPT,GYPT),DOM%D2ZSDXY(GXPT,GYPT))
    ALLOCATE(DOM%RELAX(GXPT,GYPT))
    
    ALLOCATE(XG%U(NLEV,GXPT,GYPT),XG%V(NLEV,GXPT,GYPT),XG%DW(NLEV,GXPT,GYPT),XG%GW(NLEV,GXPT,GYPT))
    ALLOCATE(XG%DIV(NLEV,GXPT,GYPT),XG%VOR(NLEV,GXPT,GYPT),XG%PI(NLEV,GXPT,GYPT))
    ALLOCATE(XG%T(NLEV,GXPT,GYPT),XG%Q(NLEV,GXPT,GYPT),XG%PIS(GXPT,GYPT))

    ALLOCATE(YSL%XDOT(NLEV,GXPT,GYPT),YSL%YDOT(NLEV,GXPT,GYPT),YSL%ZDOT(NLEV,GXPT,GYPT))
    ALLOCATE(YSL%XDOTF(NLEV,GXPT,GYPT),YSL%YDOTF(NLEV,GXPT,GYPT),YSL%ZDOTF(NLEV,GXPT,GYPT))
    ALLOCATE(YSL%XTRAJ(NLEV,GXPT,GYPT),YSL%YTRAJ(NLEV,GXPT,GYPT),YSL%ZTRAJ(NLEV,GXPT,GYPT))
    ALLOCATE(YSL%NXLAG(4,NLEV,GXPT,GYPT),YSL%NYLAG(4,NLEV,GXPT,GYPT),YSL%NZLAG(4,NLEV,GXPT,GYPT))
    ALLOCATE(YSL%NXLAGH(4,NLEV,GXPT,GYPT),YSL%NYLAGH(4,NLEV,GXPT,GYPT),YSL%NZLAGH(4,NLEV,GXPT,GYPT))
    ALLOCATE(YSL%XWEI(4,NLEV,GXPT,GYPT),YSL%YWEI(4,NLEV,GXPT,GYPT),YSL%ZWEI(4,NLEV,GXPT,GYPT))
    ALLOCATE(YSL%XWEIH(4,NLEV,GXPT,GYPT),YSL%YWEIH(4,NLEV,GXPT,GYPT),YSL%ZWEIH(4,NLEV,GXPT,GYPT))
    
  END SUBROUTINE SUALLOC_GBL
  !**********************************************************************************
  !**********************************************************************************
  !**********************************************************************************
  SUBROUTINE FREEMEM_GBL(XG,YSL,DOM)

    IMPLICIT NONE

    TYPE(RHS_GLB),                 INTENT(INOUT) :: XG
    TYPE(GEODOM),                  INTENT(INOUT) :: DOM
    TYPE(SL_STRUCT),               INTENT(INOUT) :: YSL

    DEALLOCATE(XG%U,XG%V,XG%T,XG%PIS,XG%DW,XG%Q,XG%DIV,XG%VOR,XG%PI)
    DEALLOCATE(DOM%YX,DOM%XY,DOM%XZ,DOM%ZX,DOM%RELAX)
    DEALLOCATE(DOM%ZS,DOM%DZSDX,DOM%DZSDY,DOM%D2ZSDX2,DOM%D2ZSDY2,DOM%D2ZSDXY)

    DEALLOCATE(YSL%XDOT,YSL%YDOT,YSL%ZDOT)
    DEALLOCATE(YSL%XDOTF,YSL%YDOTF,YSL%ZDOTF)
    DEALLOCATE(YSL%XTRAJ,YSL%YTRAJ,YSL%ZTRAJ)
    DEALLOCATE(YSL%NXLAG,YSL%NYLAG,YSL%NZLAG)
    DEALLOCATE(YSL%NXLAGH,YSL%NYLAGH,YSL%NZLAGH)
    DEALLOCATE(YSL%XWEI,YSL%YWEI,YSL%ZWEI)
    DEALLOCATE(YSL%XWEIH,YSL%YWEIH,YSL%ZWEIH)
    
  END SUBROUTINE FREEMEM_GBL
  !**********************************************************************************
  !**********************************************************************************
  !**********************************************************************************
  SUBROUTINE FREEMEM_MPI(X0,XC,M2,M1,M0,L0,LC,GEO)

    IMPLICIT NONE
    
    TYPE(VARPROG),                 INTENT(INOUT) :: X0
    TYPE(VARLBC),                  INTENT(INOUT) :: XC
    TYPE(RHS_LOC),                 INTENT(INOUT) :: M2,M1,M0,L0
    TYPE(RHSLBC),                  INTENT(INOUT) :: LC
    TYPE(GEO_NHNL),                INTENT(INOUT) :: GEO

    DEALLOCATE(XC%U,XC%V,XC%T,XC%DW,XC%Q,XC%PIS)
    DEALLOCATE(XC%DTDX,XC%DQDX,XC%DIV,XC%DPISDX)
    DEALLOCATE(XC%DTDY,XC%DQDY,XC%DPISDY)
    
    DEALLOCATE(X0%U,X0%V,X0%T,X0%GW,X0%Q,X0%DW,X0%DIV,X0%X,X0%Y,X0%PIS)
    DEALLOCATE(X0%DTDX,X0%DQDX,X0%DPISDX,X0%DTDY,X0%DQDY,X0%DPISDY,X0%ETADOT)
    
    DEALLOCATE(M0%U,M0%V,M0%T,M0%DW,M0%Q,M0%GW,M0%PI,M0%PIS)
    DEALLOCATE(M1%U,M1%V,M1%T,M1%DW,M1%Q,M1%GW,M1%PI,M1%PIS)
    DEALLOCATE(M2%U,M2%V,M2%T,M2%DW,M2%Q,M2%GW,M2%PI,M2%PIS)
    DEALLOCATE(L0%U,L0%V,L0%T,L0%DW,L0%Q,L0%GW,L0%PI,L0%PIS)
    DEALLOCATE(LC%U,LC%V,LC%T,LC%DW,LC%Q,LC%PIS,LC%RELAX)

    DEALLOCATE(GEO%PIS,GEO%PIF,GEO%PIH,GEO%DELTA,GEO%ALPHA,GEO%EPS,GEO%XIM,GEO%XIP)
    DEALLOCATE(GEO%DELTAPI,GEO%DELTB,GEO%DELTC,GEO%DDELTADX,GEO%DDELTADY,GEO%DALPHADX)
    DEALLOCATE(GEO%ZS,GEO%DZSDX,GEO%DZSDY,GEO%D2ZSDX2,GEO%D2ZSDY2,GEO%D2ZSDXY,GEO%DALPHADY)
    
  END SUBROUTINE FREEMEM_MPI

  !*****************************************************!
  !*****************************************************!
  !*****************************************************!

END MODULE MOD_STRUCTURE

!=====================================================!
!=====================================================!
!=====================================================!
