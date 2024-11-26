!##################   PARAMETRES   ###################!
!#                                                   #!
!# auteur : F. Voitus                                #!
!# sujet  : Définit les parametres du modèle         #!
!#                                                   #!
!#####################################################!

!=====================================================!
!================  MODULE PARAMETRES  ================!
!=====================================================!

MODULE MOD_SETUP

  USE MOD_SHARE, ONLY : NLEV,NXPT,ZERO,ONE,TWO,HALF,RDX,RDZ,RLX,RMILIEU,RZ1,RP,    &
                      & RX0,RZ0,RZC,RXC,RN,RG,RCP,RR,RKAPPA,RTSUR,RS,RAX,RAZ,RZ2,  &
                      & RPSUR,RPI,RP00,RU00,RT00,RTSUR,RDT,RC,LPERIO,LSPEC,LALIAZ, &
                      & LALIAZ,LN000,LN001,LN002,RTTROPO,RDELTR,NCASE,RPSTA,RQPP,  &
                      & RVFE_BETA,RVFE_ALPHA,NLBC
                          
  IMPLICIT NONE
  
  !* type variables 
  TYPE STATEVAR
     REAL(8), DIMENSION(NLEV,NXPT)        :: U,ETADOT
     REAL(8), DIMENSION(0:NLEV,NXPT)      :: ETADOTH
     REAL(8), DIMENSION(NLEV:NLEV,NXPT)   :: PIS  
     REAL(8), DIMENSION(NLEV,NXPT)        :: M,Q
     REAL(8), DIMENSION(NLEV,NXPT)        :: DIV
  END TYPE STATEVAR

  !* type RHS
  TYPE RHSVAR
     REAL(8), DIMENSION(NLEV,NXPT)        :: Q
     REAL(8), DIMENSION(NLEV,NXPT)        :: M
     REAL(8), DIMENSION(NXPT)             :: PIS  
  END TYPE RHSVAR

  !* type background state
  TYPE BACKGROUND
     REAL(8), DIMENSION(NLEV,NXPT)        :: U,UA   
     REAL(8), DIMENSION(NLEV,NXPT)        :: T
     REAL(8), DIMENSION(NLEV:NLEV,NXPT)   :: PIS
  END TYPE BACKGROUND
  
  !* type geometry
  TYPE GEOMETRY
     REAL(8), DIMENSION(0:NLEV)           :: AH,BH            
     REAL(8), DIMENSION(0:NLEV)           :: ETAH
     REAL(8), DIMENSION(0:NLEV+1)         :: ETA
     REAL(8), DIMENSION(0:NLEV)           :: DETAH,RDETAH        
     REAL(8), DIMENSION(NLEV)             :: DETA,RDETA
     REAL(8), DIMENSION(NLEV)             :: DELB,DELA
     REAL(8), DIMENSION(NLEV:NLEV,NXPT)   :: ZS,DZSDX          
     REAL(8), DIMENSION(NLEV,NXPT)        :: X                 
     REAL(8), DIMENSION(0:NLEV,NXPT)      :: Z                      
     REAL(8), DIMENSION(NLEV:NLEV,NXPT)   :: PIS_REF,DPISDX_REF
     REAL(8), DIMENSION(NXPT)             :: RELAX
  END TYPE GEOMETRY

  !* type semi-Lag
  TYPE SMILAG_STRUCT
     REAL(8),    DIMENSION(NLEV,NXPT)     :: XDOT,    ETADOT
     INTEGER(8), DIMENSION(6,NLEV,NXPT)   :: NXLAG,   NZLAG
     REAL(8),    DIMENSION(6,NLEV,NXPT)   :: XWEI,    ZWEI
     REAL(8),    DIMENSION(NLEV,NXPT)     :: XTRAJ,   ZTRAJ
     REAL(8),    DIMENSION(NXPT)          :: XDOTS
     INTEGER(8), DIMENSION(6,NXPT)        :: NXLAGS
     REAL(8),    DIMENSION(6,NXPT)        :: XWEIS
     REAL(8),    DIMENSION(NLEV,NXPT)     :: DXDOT_DX,  DXDOT_DZ
     REAL(8),    DIMENSION(NLEV,NXPT)     :: DZDOT_DX,  DZDOT_DZ
     REAL(8),    DIMENSION(NLEV,NXPT)     :: ALPHA_DX,  ALPHA_DZ
     REAL(8),    DIMENSION(NXPT)          :: ALPHAS_DX
  END TYPE SMILAG_STRUCT

  !* type mpdata 
  TYPE MPDATA_STRUCT
     REAL(8), DIMENSION(0:3,NLEV,NXPT)    :: MF
     REAL(8), DIMENSION(NLEV,NXPT)        :: XDOTF
     REAL(8), DIMENSION(3,0:NLEV,NXPT)    :: ZDOTH
     REAL(8), DIMENSION(NXPT)             :: XDOTS
  END TYPE MPDATA_STRUCT
  
  
CONTAINS

  !*****************************************************!
  !***************  LISTE DES ROUTINES  ****************!
  !*****************************************************!

  !**********************************************************************************!  
  !**********************************************************************************!
  !**********************************************************************************!
  !**********************************************************************************!
  SUBROUTINE SET_DOMAIN_GEOMETRY(GEO)

    USE MOD_SHARE, ONLY : RHEIGHT,DRHEIGHTDX,LSPEC,LALIAZ,RDTDZ_STD,&
                        & CTYPE_AB_LEVEL,CTYPE_ETAH,LREFINE_SLOROG
    USE MOD_GRAD,  ONLY : GRADP
    
    IMPLICIT NONE

    TYPE(GEOMETRY),          INTENT(INOUT) :: GEO
    
    INTEGER(8)                             :: J,I
    REAL(8)                                :: ZX,ZRELAX
    REAL(8)                                :: ZS,ZA,ZH
    REAL(8), DIMENSION(NXPT)               :: ZSLOPE
    LOGICAL                                :: LL_ORO



    GEO%RELAX(:)=ZERO

    DO I=1,NLBC
       ZX     = REAL(NLBC-I,8)/REAL(NLBC-1,8)
       GEO%RELAX(I)         = ((RP+ONE)*(ZX**RP)-RP*(ZX**(RP+ONE)))
       GEO%RELAX(NXPT-I+1)  = ((RP+ONE)*(ZX**RP)-RP*(ZX**(RP+ONE)))
    END DO
    
    !----------------------------------------------------------
    !.... orography-related definition
    DO I = 1, NXPT
       GEO%ZS(NLEV,I) = (ONE-GEO%RELAX(I))*RHEIGHT(REAL(I-1,8)*RDX-RXC)
       ZSLOPE(I)      = DRHEIGHTDX(REAL(I-1,8)*RDX-RXC)
    END DO

    CALL GRADP(GEO%DZSDX,GEO%ZS,NLEV,LSPEC,LALIAZ)

    WRITE(*,*) " "                                              
    WRITE(*,*) 'MAX ANA. OROGRAPHY SLOPE                : '  &   
         & ,100*MAXVAL(ZSLOPE),   ' (%), '                   &   
         & ,(180/RPI)*ATAN(MAXVAL(ZSLOPE)),   ' (°)'            
    WRITE(*,*) 'MAX NUM. OROGRAPHY SLOPE                : '  &   
         & ,100*MAXVAL(GEO%DZSDX),' (%), '                   &   
         & ,(180/RPI)*ATAN(MAXVAL(GEO%DZSDX)),' (°)'

    
    !*-----------------------------------------------------------
    !------ hydrostatic surface pressure of reference 
    ZA = (RR/RG)*RDTDZ_STD
    ZH = (RR/RG)*RT00
    
    IF (LREFINE_SLOROG) THEN
       DO I=1,NXPT
          ZS = GEO%ZS(NLEV,I)
          GEO%PIS_REF(NLEV,I) = RP00*DEXP(-(ZS/ZH))*(        &
                              & ONE+(ZA/2.d0)*(ZS/ZH)        &
                              & -((ZA**2)/3.d0)*((ZS/ZH)**2) )
       END DO 
    ELSE
       DO I=1,NXPT
          GEO%PIS_REF(NLEV,I) = RP00*DEXP(-GEO%ZS(NLEV,I)/ZH) 
       END DO
    END IF   
    CALL GRADP(GEO%DPISDX_REF,GEO%PIS_REF,NLEV,LSPEC,LALIAZ)

    !----------------------------------------------------------
    !*Définition de la géométrie du problème 
    DO I = 1,NXPT
       GEO%X(:,I)    = REAL(I-1,8)*RDX
       GEO%Z(NLEV,I) = GEO%ZS(NLEV,I)
       DO J = 0, NLEV-1
          GEO%Z(J,I) = (REAL(NLEV-J,8))*RDZ
       END DO
    END DO

    !*--------------------------------------------------------
    !* define A(eta) and B(eta) vertical functions
    CALL SET_AB_FUNCTIONS(GEO%AH,GEO%BH,CTYPE_AB_LEVEL)

    !*--------------------------------------------------------
    !* Explicit definition for eta-coordinate
    CALL SET_ETA_LEVELS(GEO,CTYPE_ETAH)

    
  END SUBROUTINE SET_DOMAIN_GEOMETRY
  !*******************************************************************************!
  !*******************************************************************************!
  !*******************************************************************************!
  !*******************************************************************************!
  !*******************************************************************************!
  SUBROUTINE SET_P_FROM_Z(PP)
    
    USE MOD_SHARE, ONLY : RDZMAX,RDZMIN,RZLIM,NILEV,RZLOW,LREGETA
    USE MOD_DIAG,  ONLY : DEF_DZ
    
    IMPLICIT NONE

    REAL(8), DIMENSION(0:NLEV), INTENT(OUT)  :: PP
    
    INTEGER                                  :: J
    REAL(8)                                  :: ZDZDH,ZEPS
    REAL(8), DIMENSION(NLEV)                 :: ZZF,ZDZ,ZT,ZDH

    ZEPS = 1000*ABS(EPSILON(ONE))

    !-----------------------------------------------------------------
    ! Basic-state height
    CALL DEF_DZ(ZDZ,RDZMIN,RDZMAX,RZLIM,NILEV,LREGETA)

    ZZF(NLEV)      = RZLOW
    DO J =1,NLEV-1
       ZZF(NLEV-J) = ZZF(NLEV-J+1) + ZDZ(NLEV-J+1)
    END DO
    
    !-----------------------------------------------------------------
    ! Basic-state Temperature profil
    IF (LN000) THEN
       DO J=1,NLEV
          ZT(J) = MAX(RTSUR-(RG/RCP)*ZZF(J),RTTROPO)
       END DO
    ELSE IF (LN001) THEN
       DO J=1,NLEV
          ZT(J) = MAX((((RG**2)/(RN**2))/RCP)*  &
                & (ONE-DEXP((RN**2/RG)*ZZF(J))) &
                & +RTSUR*DEXP((RN**2/RG)*ZZF(J)),RTTROPO)
       END DO   
    ELSE IF (LN002) THEN
       ZT(:) = ((RG**2)/(RN**2))/RCP
    ELSE
       ZT(:) = RTSUR
    END IF

    !-----------------------------------------------------------------
    ! Basic-state pressure profil
    PP(0)    = ZERO
    PP(NLEV) = RPSUR
    DO J = 1, NLEV-1
       ZDZDH      = MAX((RG/RR)*(ZDZ(NLEV-J+1)/TWO)/ZT(NLEV-J+1),ZEPS)
       PP(NLEV-J) = PP(NLEV-J+1)*((DSQRT((ZDZDH**2)+ONE)-ZDZDH)**2)
    END DO

    !PP(0)    = ZERO
    !ZDH(1)   = MAX((RG/RR)*(ZDZ(1)/ZT(1)),ZEPS)
    !PP(NLEV) = RPSUR
    !DO J = 1, NLEV-1
    !   ZDH(NLEV-J+1) = MAX((RG/RR)*(ZDZ(NLEV-J+1)/ZT(NLEV-J+1)),ZEPS)
    !   PP(NLEV-J)    = PP(NLEV-J+1)*DEXP(-ZDH(NLEV-J+1))
    !END DO

  END SUBROUTINE SET_P_FROM_Z
  !***********************************************************************************!
  !***********************************************************************************!
  !***********************************************************************************!
  !***********************************************************************************!
  !***********************************************************************************!
  !***********************************************************************************!
  SUBROUTINE SET_AB_FUNCTIONS(PAH,PBH,C_OPT)

    IMPLICIT NONE

    REAL(8), DIMENSION(0:NLEV),         INTENT(OUT) :: PAH,PBH
    CHARACTER(LEN=*),                   INTENT(IN)  :: C_OPT

    INTEGER                                         :: I,J
    REAL(8), DIMENSION(0:90)                        :: ZA1,ZB1
    REAL(8), DIMENSION(0:137)                       :: ZA2,ZB2
    REAL(8), DIMENSION(0:NLEV)                      :: ZPP
    
    PAH(:) = ZERO
    PBH(:) = ZERO

    IF (C_OPT == 'SIG') THEN
       CALL SET_P_FROM_Z(ZPP)
       DO J = 1,NLEV
          PBH(J) = ZPP(J)/ZPP(NLEV)
       END DO
    ELSE IF ((C_OPT == 'ARO').AND.(NLEV == 90)) THEN
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
       
    ELSE IF ((C_OPT == 'IFS').AND.(NLEV == 137)) THEN
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

  END SUBROUTINE SET_AB_FUNCTIONS
  !*******************************************************************************!
  !*******************************************************************************!
  !*******************************************************************************!
  !*******************************************************************************!
  SUBROUTINE SET_ETA_LEVELS(GEO,CTYPE_ETAH)
  
    IMPLICIT NONE

    TYPE(GEOMETRY),          INTENT(INOUT)   :: GEO
    CHARACTER(LEN=*),        INTENT(IN)      :: CTYPE_ETAH
    
    INTEGER                          :: J,I,JLEV,ITER
    REAL(KIND=8)                     :: ZP1,ZP2,ZETA_REG,ZETA_COS,ZETA
    REAL(KIND=8)                     :: ZA,ZB,ZC,ZG,ZGD,ZZ
    REAL(KIND=8)                     :: ZS(0:NLEV+1)


    !*-------------------------------------------------------------------*!
    !* various definitions of eta-levels
    IF (CTYPE_ETAH=='REGETA') THEN

      DO JLEV = 0,NLEV
        ZETA_REG        = REAL(JLEV,8)/REAL(NLEV,8)
        GEO%ETAH(JLEV) = ZETA_REG
      ENDDO
     
    ELSEIF (CTYPE_ETAH=='MCONST') THEN

      ! centripetal definition
      WRITE(*,*) "eta-levels definition = centripetal - regeta"
      ZS(0) = 0.d0
      DO JLEV = 1,NLEV
        ZP1      = (ONE/RP00)*GEO%AH(JLEV-1)+GEO%BH(JLEV-1)
        ZP2      = (ONE/RP00)*GEO%AH(JLEV)+GEO%BH(JLEV)
        ZS(JLEV) = ZS(JLEV-1)+(ZP2-ZP1)**RVFE_ALPHA
      ENDDO
      WRITE(*,*) "eta-levels definition: density increased on boundaries acc.to RVFE_BETA"
      DO JLEV = 0,NLEV
        ZETA_REG        = ZS(JLEV)/ZS(NLEV)
        ZETA_COS        = 0.5d0*(1.d0-COS(RPI*ZETA_REG))
        GEO%ETAH(JLEV) = (1.d0-RVFE_BETA)*ZETA_REG + RVFE_BETA*ZETA_COS
      ENDDO
      GEO%ETAH(0)    = 0.d0
      GEO%ETAH(NLEV) = 1.d0

    ELSEIF (CTYPE_ETAH=='MCNTRL') THEN
      
      WRITE(*,*) "VFE eta definition = dpi/deta = a + b eta + c eta^2"
      ! RVFE_ALPHA (dpi/deta_top) and RVFE_BETA (dpi/deta_surf) represent here 
      ! relative density of layers in eta space:
      ! ALPHA/BETA = 0, regular distribution
      ! ALPHA/BETA > 0, denser close to top/bottom boundary
      ! ALPHA/BETA < 0, denser towards inner domain 

      ZA = RVFE_ALPHA + 1.d0
      ZB = -2.d0*(2.d0*RVFE_ALPHA + RVFE_BETA)
      ZC =  3.d0*(RVFE_ALPHA + RVFE_BETA)

      GEO%ETAH(0) = 0.d0
      DO JLEV = 1,NLEV
        ZP2  = REAL(JLEV,8)/REAL(NLEV,8)
        ZETA = GEO%ETAH(JLEV-1)
        DO  ITER = 1, 10
          ZG   = ZETA * (ZA + ZETA * (ZB / 2.d0 + ZETA * ZC / 3.d0)) - ZP2
          ZGD  = ZA + ZETA * (ZB + ZETA * ZC)
          ZETA = ZETA - ZG / ZGD
        ENDDO
        GEO%ETAH(JLEV) = ZETA
      ENDDO
      GEO%ETAH(NLEV) = 1.d0
      
    ELSE
      
      ! Standard definition (so called chordal);
      ! 1/pref * dpi/deta == 1 by definition since A/pref + B = pi/pref = eta
      WRITE(*,*) "VFE eta definition = A/p_r + B (dpi/deta =~ 1)"
      DO JLEV=0,NLEV
        GEO%ETAH(JLEV) = (ONE/RP00)*GEO%AH(JLEV)+GEO%BH(JLEV)
      ENDDO
     
    ENDIF
   
   !*-------------------------------------------------------------------*!
   !*useful definitions *!   
    DO J=1,NLEV
       GEO%ETA(J)    = (GEO%ETAH(J)+GEO%ETAH(J-1))/TWO
    END DO
    GEO%ETA(0)       = GEO%ETAH(0)
    GEO%ETA(NLEV+1)  = GEO%ETAH(NLEV)

    DO J=1,NLEV
       GEO%DETA(J)   = GEO%ETAH(J)-GEO%ETAH(J-1)
       GEO%RDETA(J)  = ONE/GEO%DETA(J)
    END DO
    DO J=0,NLEV
       GEO%DETAH(J)  = (GEO%ETA(J+1)-GEO%ETA(J))
       GEO%RDETAH(J) = ONE/(GEO%ETA(J+1)-GEO%ETA(J))
    END DO
    DO J=1,NLEV
       GEO%DELA(J)   = GEO%AH(J)-GEO%AH(J-1)
       GEO%DELB(J)   = GEO%BH(J)-GEO%BH(J-1)
    END DO
    !*-------------------------------------------------------------------*!

  END SUBROUTINE SET_ETA_LEVELS
  !*************************************************************************************!
  !*************************************************************************************!
  !*************************************************************************************!
  !*************************************************************************************!
  !*************************************************************************************!   
  SUBROUTINE SET_BACKGROUND_STATE(STB,GEO)

    USE MOD_DIAG, ONLY : DEF_Z,DEF_P
    
    IMPLICIT NONE

    TYPE(BACKGROUND),          INTENT(INOUT) :: STB
    TYPE(GEOMETRY),            INTENT(IN)    :: GEO
    
    INTEGER                                  :: I,J
    REAL(8)                                  :: ZTS,ZIJ,ZHS,ZAK
    REAL(8), DIMENSION(NLEV,NXPT)            :: ZZ,ZP,ZU
    REAL(8), DIMENSION(0:NLEV,NXPT)          :: ZH,PSI
    
    
    !----------------------------------------------------------------------
    !* Background surface pressure
    ZHS = RPSTA*(RG/(RR*RTSUR))
    DO I=1,NXPT
       STB%PIS(NLEV,I) = RPSUR*DEXP(-GEO%ZS(NLEV,I)/ZHS)
    END DO
    
    !---------------------------------------------------------------------
    !* Background temperature   
    IF ((LN000).OR.(LN001).OR.(LN002)) THEN
       ZTS = RCP*((RN/RG)**2)*RTSUR
       DO I=1,NXPT
          CALL DEF_P(ZP(:,I),STB%PIS(NLEV,I),GEO%AH,GEO%BH,1)
          DO J=1,NLEV
             ZAK        = DEXP(RKAPPA*DLOG(ZP(J,I)/STB%PIS(NLEV,I)))
             STB%T(J,I) = MAX(RTSUR*(ZAK/(ONE+ZTS*(ZAK-ONE))),RTTROPO)
          END DO
       END DO 
    ELSE   
       STB%T(:,:) = RTSUR
    END IF      

    !----------------------------------------------------------------------
    !* Background Horizontal wind shear 
    IF (NCASE == 1) THEN
       DO I = 1,NXPT
          !* compute analytic wind shear  
          CALL DEF_Z(1,ZZ(:,I),STB%T(:,I),STB%PIS(NLEV,I), &
               & GEO%AH,GEO%BH,GEO%ZS(NLEV,I))
          DO J = 1,NLEV
             ZIJ = MAX(ZERO,MIN((ZZ(J,I)-RZ1)/(RZ2-RZ1),ONE))   
             STB%UA(J,I) = RU00*(DSIN((RPI/TWO)*ZIJ)**2)  
          END DO
    
          !* compute from stream-function defined at half-level
          CALL DEF_Z(0,ZH(:,I),STB%T(:,I),STB%PIS(NLEV,I), &
               & GEO%AH,GEO%BH,GEO%ZS(NLEV,I))
          DO J=0,NLEV
             ZIJ = ZH(J,I)
             IF (ZIJ .GT. RZ2) THEN
                PSI(J,I) = (RU00/TWO)*(RZ2+RZ1-TwO*ZIJ)
             ELSE IF (ZIJ .LE. RZ1) THEN
                PSI(J,I) = ZERO
             ELSE
                PSI(J,I) = (RU00/TWO)*( (RZ1-ZIJ) &
                         & + ((RZ2-RZ1)/RPI)*SIN(RPI*((ZIJ-RZ1)/(RZ2-RZ1))) )
             END IF   
          END DO
          DO J=1,NLEV
             STB%U(J,I)  = -(PSI(J,I)-PSI(J-1,I))*GEO%RDETA(J) 
          END DO
       END DO    
    ELSE
       STB%UA(:,:) = RU00
       STB%U(:,:)  = RU00     
    END IF
    
  END SUBROUTINE SET_BACKGROUND_STATE
  !************************************************************************************!
  !************************************************************************************!
  !************************************************************************************!
  !************************************************************************************!
  SUBROUTINE SET_INITIAL_PERTURBATION(ST,STB,GEO)

    USE MOD_DIAG, ONLY : DEF_M,DEF_Z
    
    IMPLICIT NONE

    TYPE(STATEVAR),            INTENT(INOUT) :: ST
    TYPE(BACKGROUND),          INTENT(IN)    :: STB
    TYPE(GEOMETRY),            INTENT(IN)    :: GEO
    
    INTEGER                                  :: I,J
    REAL(8)                                  :: ZR
    REAL(8), DIMENSION(NLEV,NXPT)            :: ZZ
    

    DO I = 1, NXPT
       !* compute reference hydrostatic PS
       ST%PIS(NLEV,I) = STB%PIS(NLEV,I)    

       !* compute full-level height
       CALL DEF_Z(1,ZZ(:,I),STB%T(:,I),STB%PIS(NLEV,I),     &
            & GEO%AH,GEO%BH,GEO%ZS(NLEV,I))
       
       !* define initial tracer
       DO J = 1,NLEV
          ZR          = MIN(DSQRT((((GEO%X(J,I)-RX0)/RAX)**2) &
                      & +(((ZZ(J,I)-RZ0)/RAZ)**2)),ONE)   
          ST%Q(J,I)   = RQPP*((DCOS((RPI/TWO)*ZR))**2)
       END DO
       
       !* compute pseudo-metric 
       CALL DEF_M(ST%M(:,I),STB%PIS(NLEV,I),GEO%DELA,         &
            & GEO%DELB,GEO%RDETA)
    END DO

    CALL SET_WIND_STATE(ST,STB,GEO)
    
  END SUBROUTINE SET_INITIAL_PERTURBATION
  !************************************************************************************!
  !************************************************************************************!
  !************************************************************************************!
  !************************************************************************************!
  SUBROUTINE SET_WIND_STATE(ST,STB,GEO)

    USE MOD_GRAD, ONLY : GRADP
    USE MOD_DIAG, ONLY : DEF_Z,DEF_ETADOT
    
    IMPLICIT NONE

    TYPE(STATEVAR),            INTENT(INOUT) :: ST
    TYPE(BACKGROUND),          INTENT(IN)    :: STB
    TYPE(GEOMETRY),            INTENT(IN)    :: GEO
    
    INTEGER                                  :: I,J
    REAL(8)                                  :: ZJ
    REAL(8), DIMENSION(0:NLEV,NXPT)          :: ZH,PSI
    REAL(8), DIMENSION(NLEV:NLEV,NXPT)       :: ZDPISDX

    
    DO I=1,NXPT
       
       !* compute height at half-levels
       CALL DEF_Z(0,ZH(:,I),STB%T(:,I),ST%PIS(NLEV,I), &
            & GEO%AH,GEO%BH,GEO%ZS(NLEV,I))

       IF (NCASE == 1) THEN
          !* compute stream-function
          DO J=0,NLEV
             ZJ = ZH(J,I)
             IF (ZJ .GT. RZ2) THEN
                PSI(J,I) = (RU00/TWO)*(RZ2+RZ1-TwO*ZJ)
             ELSE IF (ZJ .LE. RZ1) THEN
                PSI(J,I) = ZERO
             ELSE
                PSI(J,I) = (RU00/TWO)*( (RZ1-ZJ) &
                         & + ((RZ2-RZ1)/RPI)*SIN(RPI*((ZJ-RZ1)/(RZ2-RZ1))) )
             END IF   
          END DO
       
          !* horizontal wind component
          DO J=1,NLEV
             ST%U(J,I)   = -(PSI(J,I)-PSI(J-1,I))*GEO%RDETA(J)
          END DO
       ELSE
          ST%U(:,I) = STB%U(:,I)
       END IF
       
    END DO

    !* compute horizontal divergence
    CALL GRADP(ST%DIV,ST%U,1_8,LSPEC,LALIAZ)

    IF (NCASE == 1) THEN
       !* vertical velocity
       CALL GRADP(ST%ETADOTH,PSI,0_8,LSPEC,LALIAZ)
       DO J=1,NLEV
          ST%ETADOT(J,:) = HALF*(ST%ETADOTH(J,:)          &
                         & + ST%ETADOTH(J-1,:))
       END DO   
    ELSE
       !* compute horizontal divergence
       CALL GRADP(ZDPISDX,ST%PIS,NLEV,LSPEC,LALIAZ)
       DO I=1,NXPT
          CALL DEF_ETADOT(ST%ETADOT(:,I),ST%ETADOTH(:,I), &
               & ST%U(:,I),ST%DIV(:,I),ST%PIS(NLEV,I),    &
               & ZDPISDX(NLEV,I),GEO%DELA,GEO%DELB,       &
               & GEO%RDETA,GEO%AH,GEO%BH,GEO%RDETAH)
       END DO
    END IF
    
  END SUBROUTINE SET_WIND_STATE
  !************************************************************************************!
  !************************************************************************************!
  !************************************************************************************!
  !************************************************************************************!
  
END MODULE MOD_SETUP

!=====================================================!
!=====================================================!
!=====================================================!
