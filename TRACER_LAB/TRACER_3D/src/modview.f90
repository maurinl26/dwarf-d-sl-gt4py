!####################   SORTIE   #####################!
!#                                                   #!
!# auteur : F.Voitus                                 #!
!# sujet  : Gère les différentes sorties du modèle   #!
!#                                                   #!
!#####################################################!

!=====================================================!
!==================  MODULE SORTIE  ==================!
!=====================================================!

MODULE MOD_VIEW

  USE MOD_PARAMETER, ONLY : NLEV,GXPT,GYPT,RLY,RLX,RDX,RDY,RDZ,RR,RCP,RG,RKAPPA,RPI,  &
                          & NTIMESTEP,NDLNPR,NFPLEV,ZERO,ONE,TWO,RXC,RYC,RX0,RY0,RZ0, &
                          & RXS,RYS,RZS,LSLICE_XY,LSLICE_XZ,LBIGW,RDT
  USE MOD_STRUCTURE, ONLY : GEODOM,RHS_GLB,BACKGROUND
  
CONTAINS

  !*****************************************************!
  !***************  LISTE DES ROUTINES  ****************!
  !*****************************************************!

  !*********  Affichage en-tête et nettoyage   *********!
  SUBROUTINE SCREEN_START(TEMPSOUT,JOUROUT,DOM,STB)

    USE MOD_DIAGNOSTIC,  ONLY: DIAG_PHI,DIAG_ALPHADELTA

    IMPLICIT NONE

    INTEGER,           INTENT(OUT) :: TEMPSOUT,JOUROUT
    TYPE(BACKGROUND),  INTENT(IN)  :: STB
    TYPE(GEODOM),      INTENT(IN)  :: DOM
    
    INTEGER                        :: L
    INTEGER, DIMENSION(3)          :: JJOUR,JMAINTENANT
    REAL(8), DIMENSION(NLEV)       :: ZPIF,ZALPHA,ZDELTA,ZDELTAPI
    REAL(8), DIMENSION(NLEV)       :: ZPHIF,ZNUL
    REAL(8), DIMENSION(0:NLEV)     :: ZPHIH

    CALL IDATE(JJOUR)
    JOUROUT = JJOUR(2)
    
    ZNUL(:)=ZERO
    
    CALL DIAG_ALPHADELTA(ZPIF,ZALPHA,ZDELTA,ZDELTAPI,STB%PIS,DOM%AH,DOM%BH)
    CALL DIAG_PHI(0_8,ZPHIH,STB%T,ZNUL,ZALPHA,ZDELTA,ZERO)
    CALL DIAG_PHI(1_8,ZPHIF,STB%T,ZNUL,ZALPHA,ZDELTA,ZERO)
    
     WRITE(*,*) "     **************************"
    WRITE(*,*) "     *       3D EE MODEL      *"
    WRITE(*,*) "     **************************"
    WRITE(*,*) " "
    IF (JJOUR(2) .LT. 10) THEN
       IF (JJOUR(1) .LT. 10) THEN
          WRITE(*,'(A12,I1,A2,I1,A1,I4)') 'date     : ',JJOUR(1),'/0',JJOUR(2),'/',JJOUR(3)
       ELSE
          WRITE(*,'(A12,I2,A2,I1,A1,I4)') 'date     : ',JJOUR(1),'/0',JJOUR(2),'/',JJOUR(3)
       END IF
    ELSE
       IF (JJOUR(1) .LT. 10) THEN
          WRITE(*,'(A12,I1,A1,I2,A1,I4)') 'date     : ',JJOUR(1),'/',JJOUR(2),'/',JJOUR(3)
       ELSE
          WRITE(*,'(A12,I2,A1,I2,A1,I4)') 'date     : ',JJOUR(1),'/',JJOUR(2),'/',JJOUR(3)
       END IF
    END IF
    WRITE(*,*)
    WRITE(*,*)
    WRITE(*,'(A11,1X,L2)') ' LBIGW    : ',LBIGW
    WRITE(*,'(A11,1X,I2)') ' NDLNPR   : ',NDLNPR
    WRITE(*,*)
    WRITE(*,*)
    WRITE(*,'(A15)') ' BASIC STATE '
    WRITE(*,*)
    WRITE(*,'(6X,A6,5X,A6,5X,A6,5X,A6,5X,A14,5X,A12,5X,A9,5X,A14)') &
         & 'niveau', 'A(eta)', 'B(eta)','S(eta)','pression(hPa)', &
         & 'altitude(km)','DeltaZ(m)','temperature(K)'
    DO L = 1, NLEV
       WRITE(*,'(6X,I3,7X,F6.3,5X,F6.3,5X,F6.3,7X,F8.2,12X,F6.2,6X,F12.6,6X,F8.2)') L, &
            & DOM%AH(L),DOM%BH(L),DOM%SH(L),ZPIF(L)*0.01d0,(ZPHIF(L)/RG)*0.001d0,      &
            & (ZPHIH(L-1)-ZPHIH(L))/RG,STB%T(L)
    END DO
    WRITE(*,*)
    WRITE(*,'(A23,6X,A20)') 'PARAMETRES DYNAMIQUES :' , 'DOMAINE NUMÉRIQUE :'
    WRITE(*,*)
    WRITE(*,'(1X,A7,I6,16X,A11,I3,A2,I2,A4,I2,A1)')'NSTEP = ',NTIMESTEP,'Echeance : ',                     &
             & INT(NTIMESTEP*RDT)/3600,'h ',(INT(NTIMESTEP*RDT)-(INT(NTIMESTEP*RDT)/3600)*3600)/60,'min ', &
             & INT(NTIMESTEP*RDT)-INT(NTIMESTEP*RDT)/(3600**2) -(INT(NTIMESTEP*RDT)                        &
             & -INT(NTIMESTEP*RDT)/(3600**2))/(60**2),'s'
    WRITE(*,*)
    WRITE(*,'(1X,A7,I6,16X,A8,I6,A3)')'GXPT = ', GXPT,'Lon.    : ', INT(RLX/1000),' km'
    WRITE(*,'(1X,A7,I6,16X,A8,I6,A3)')'GYPT = ', GYPT,'Lat.    : ', INT(RLY/1000),' km'
    WRITE(*,'(1X,A7,I6,16X,A8,I4,A3)')'NLEV = ', NLEV,'Levels. : ', INT(ZPHIH(0)/1000/RG),' km'
    WRITE(*,'(1X,A7,F6.2,A2)')'dt   = ', RDT, ' s'
    WRITE(*,'(1X,A7,I6,A2)')  'dx   = ', INT(RDX),' m'
    WRITE(*,'(1X,A7,I6,A2)')  'dy   = ', INT(RDY),' m'
    WRITE(*,'(1X,A7,I6,A2)')  'dz   = ', INT(RDZ),' m'
    WRITE(*,*)
   
   
    CALL ITIME(JMAINTENANT)
    WRITE(*,'(A20,I2,A1,I2,A1,I2)') '  début du calcul : ',&
         & JMAINTENANT(1),':',JMAINTENANT(2),':',JMAINTENANT(3)
    TEMPSOUT = JMAINTENANT(1)*3600 + JMAINTENANT(2)*60 + JMAINTENANT(3)

    WRITE(*,*)
    WRITE(*,'(A54)') " 0%  10%  20%  30%  40%  50%  60%  70%  80%  90%  100%"
    WRITE(*,'(A52)') " |----|----|----|----|----|----|----|----|----|----|"
    WRITE(*,'(A2)',ADVANCE='NO') ' *'

  END SUBROUTINE SCREEN_START
  !*****************************************************!
  !**********  Affichage durant la dynamique  **********!
  SUBROUTINE SCREEN(TEMPSINOUT,RHS,DOM,BOUCLEIN)

    IMPLICIT NONE

    INTEGER,                            INTENT(INOUT) :: TEMPSINOUT
    INTEGER(8),                         INTENT(IN)    :: BOUCLEIN
    TYPE(RHS_GLB),                      INTENT(IN)    :: RHS
    TYPE(GEODOM),                       INTENT(IN)    :: DOM
    
    INTEGER, DIMENSION(3)                             :: IMAINTENANT,IJOUR
    CHARACTER(len=3)                                  :: STARS    
    REAL(8)                                           :: ZMAXU,ZMAXV,ZSEUIL
    LOGICAL                                           :: LL_UNSTABLE

    ZMAXU=MAXVAL(ABS(RHS%U))
    ZMAXV=MAXVAL(ABS(RHS%V))
    ZSEUIL = 300.d0

    LL_UNSTABLE = (ZMAXU .GT. ZSEUIL).OR.(ZMAXV .GT. ZSEUIL)
    
    IF (LL_UNSTABLE) THEN
       WRITE(*,*) "\n  UNSTABLE SCHEME : WIND EXPLOSION AT STEP : ",BOUCLEIN
       CALL ITIME(IMAINTENANT)
       WRITE(*,*)
       WRITE(*,'(A20,I2,A1,I2,A1,I2)') "    fin du calcul : ", &
            & IMAINTENANT(1),":",IMAINTENANT(2),':',IMAINTENANT(3)

       TEMPSINOUT     = IMAINTENANT(1)*3600+IMAINTENANT(2)*60+IMAINTENANT(3)-TEMPSINOUT
       IMAINTENANT(1) = TEMPSINOUT/3600
       IMAINTENANT(2) = (TEMPSINOUT-IMAINTENANT(1)*3600)/60
       IMAINTENANT(3) = TEMPSINOUT-IMAINTENANT(1)*3600-IMAINTENANT(2)*60
       WRITE(*,'(A20,I2,A2,I2,A4,I2,A1)') "  temps de calcul : ",IMAINTENANT(1),"h ",&
            & IMAINTENANT(2),"min ",IMAINTENANT(3),'s'
       WRITE(*,*) "\n Fin du programme \n"
       STOP
    END IF

    IF ((MOD(50*BOUCLEIN,NTIMESTEP) .EQ. 0) .AND. (BOUCLEIN .NE. 0)) THEN
       WRITE(*,'(A1)',ADVANCE='NO') "*"
    END IF

    !WRITE(*,'(F20.16,2X,F20.16,2X,F10.6)') ZMAXU, ZMAXV, ZMINT
    
  END SUBROUTINE SCREEN
  !*********************************************************************************************************!
  !*********************************************************************************************************!
  !*********************************************************************************************************!
  SUBROUTINE SCREEN_END(TEMPSINOUT,JOURIN)

    IMPLICIT NONE

    INTEGER, INTENT(INOUT) :: TEMPSINOUT
    INTEGER, INTENT(IN)    :: JOURIN
    INTEGER, DIMENSION(3)  :: ITEMPS,IJOUR

    WRITE(*,*)
    CALL ITIME(ITEMPS)
    WRITE(*,*)
    WRITE(*,'(A20,I2,A1,I2,A1,I2)') "    fin du calcul : ",ITEMPS(1),":",ITEMPS(2),':',ITEMPS(3)
    CALL IDATE(IJOUR)
    IF (JOURIN .EQ. IJOUR(2)) THEN
       TEMPSINOUT = ITEMPS(1)*3600+ITEMPS(2)*60+ITEMPS(3) - TEMPSINOUT
       ITEMPS(1) = TEMPSINOUT/3600
       ITEMPS(2) = (TEMPSINOUT-ITEMPS(1)*3600)/60
       ITEMPS(3) = TEMPSINOUT-ITEMPS(1)*3600-ITEMPS(2)*60
       WRITE(*,'(A20,I2,A2,I2,A4,I2,A1)') "  temps de calcul : ",ITEMPS(1),"h ",ITEMPS(2),"min ",ITEMPS(3),"s"
    END IF

    WRITE(*,*) '\n Fin du programme \n '

  END SUBROUTINE SCREEN_END
  !*********************************************************************************************************!
  !*********************************************************************************************************!
  !*********************************************************************************************************!
  SUBROUTINE OUTPUT_INITIAL_FIELDS(PT,PDTDX,PDTDY,PPIS,PDPISDX,PDPISDY,PDW,STB,DOM)

    USE MOD_DIAGNOSTIC,  ONLY : DIAG_BACKGROUND_THERMO,DIAG_THERMO,DIAG_PHI,&
                              & DIAG_WPOS,DIAG_ALPHADELTA
    
    IMPLICIT NONE

    TYPE(GEODOM),                               INTENT(IN)    :: DOM
    TYPE(BACKGROUND),                           INTENT(IN)    :: STB
    REAL(8), DIMENSION(NLEV,GXPT,GYPT),         INTENT(IN)    :: PDW
    REAL(8), DIMENSION(NLEV,GXPT,GYPT),         INTENT(IN)    :: PT,PDTDX,PDTDY
    REAL(8), DIMENSION(GXPT,GYPT),              INTENT(IN)    :: PPIS,PDPISDX,PDPISDY
    
    INTEGER                                                   :: NY,NX
    INTEGER                                                   :: K, I ,J
    INTEGER                                                   :: K0,I0,J0
    INTEGER                                                   :: NYS
    INTEGER(8), DIMENSION(NLEV)                               :: JLOC

    REAL(8)                                                   :: ZGAMMA, ZSECUR, ZDIF
    REAL(8)                                                   :: PIJ, PLAST, ZIJ
    
    REAL(8), DIMENSION(NLEV)                                  :: ZPIF,ZALPHA,ZDELTA,ZDELTAPI
    REAL(8), DIMENSION(NLEV)                                  :: ZPHIF,ZTH,ZTH0,ZNUL
    REAL(8), DIMENSION(0:NLEV)                                :: ZPHIH,ZW
     
    REAL(8), DIMENSION(GXPT,GYPT)                             :: ZZXY
    REAL(8), DIMENSION(NLEV,GXPT)                             :: ZTXZ, ZWXZ
   
    LOGICAL, DIMENSION(NLEV)                                  :: LLDONE, LLEXTRA

    ZNUL(:)   = ZERO
    ZSECUR    = ABS(1000*EPSILON(1.d0))
    NYS       = 1+INT(RYS/RDY)
    
    
    ! X-Z CROSS-SECTIONS
    IF (LSLICE_XZ) THEN
        
       DO NX = 1, GXPT
   
          CALL DIAG_ALPHADELTA(ZPIF,ZALPHA,ZDELTA,ZDELTAPI,PPIS(NX,NYS),DOM%AH,DOM%BH)
           
          ! INTERPOLATION POTENTIAL TEMPERATURE
          CALL DIAG_THERMO(PT(:,NX,NYS),ZNUL,ZPIF,ZTH,'TH')
          CALL DIAG_BACKGROUND_THERMO(STB%T,STB%PIS,DOM%AH,DOM%BH,ZTH0,'TH')
          CALL DIAG_PHI(1_8,ZPHIF,PT(:,NX,NYS),ZNUL,ZALPHA,ZDELTA,DOM%ZS(NX,NYS))
           
          PLAST = ZPHIF(NLEV)
          DO J= 1, NLEV-1
             LLDONE(J) = .FALSE.
             PIJ       = RG*RDZ*REAL(NLEV-J,8) 
             JLOC(J)   = J
             DO K = 1, NLEV-1
                ZDIF = (ZPHIF(K)-PIJ)*(PIJ-ZPHIF(K+1))        
                IF (ZDIF > -ZSECUR) THEN
                  JLOC(J)   = K
                  LLDONE(J) = .TRUE.
                  EXIT
                END IF 
             END DO   
             IF (LLDONE(J)) THEN
                LLEXTRA(J) = .FALSE.
             ELSE
                LLEXTRA(J) = (PLAST > PIJ)
             END IF  
             IF (.NOT.LLEXTRA(J)) THEN  
                ZGAMMA      = (ZPHIF(JLOC(J))-PIJ)/(ZPHIF(JLOC(J))-ZPHIF(JLOC(J)+1))
                ZTXZ(J,NX)  = ZGAMMA*(ZTH(JLOC(J)+1)-ZTH0(JLOC(J)+1)) &
                            & +(ONE-ZGAMMA)*(ZTH(JLOC(J))-ZTH0(JLOC(J)))
             ELSE
                ZTXZ(J,NX)  = ZTH(NLEV)-ZTH0(NLEV)
             END IF
          END DO
          ZTXZ(NLEV,NX)     = ZTH(NLEV)-ZTH0(NLEV)

          ! INTERPOLATION VERTICAL WIND
          CALL DIAG_WPOS(ZW,PDW(:,NX,NYS),STB%U,STB%V,PT(:,NX,NYS),PDTDX(:,NX,NYS),PDTDY(:,NX,NYS),      &
               & ZNUL,ZNUL,ZNUL,PPIS(NX,NYS),PDPISDX(NX,NYS),PDPISDY(NX,NYS),ZPIF,ZALPHA,ZDELTA,ZDELTAPI, &
               & DOM%AH,DOM%BH,DOM%SH,DOM%DZSDX(NX,NYS),DOM%DZSDY(NX,NYS),0_8)
          CALL DIAG_PHI(0_8,ZPHIH,PT(:,NX,NYS),ZNUL,ZALPHA,ZDELTA,DOM%ZS(NX,NYS))

          PLAST = ZPHIH(NLEV)
          DO J = 1,NLEV-1
            LLDONE(J) = .FALSE.
            PIJ       = RG*RDZ*REAL(NLEV-J,8)
            JLOC(J)   = J
            DO K = 1, NLEV
              ZDIF = (ZPHIH(K-1)-PIJ)*(PIJ-ZPHIH(K))
              IF (ZDIF > -ZSECUR) THEN
                JLOC(J)   = K-1
                LLDONE(J) = .TRUE.
                EXIT
              END IF 
            END DO   
            IF (LLDONE(J)) THEN
               LLEXTRA(J) = .FALSE.
            ELSE
               LLEXTRA(J) = (PLAST > PIJ)
            END IF
            IF (.NOT.LLEXTRA(J)) THEN
               ZGAMMA     = (ZPHIH(JLOC(J))-PIJ)/(ZPHIH(JLOC(J))-ZPHIH(JLOC(J)+1))
               ZWXZ(J,NX) = ZGAMMA*ZW(JLOC(J)+1)+(ONE-ZGAMMA)*ZW(JLOC(J)) 
            ELSE
               ZWXZ(J,NX) = ZW(NLEV)
            END IF
          END DO      
           ZWXZ(NLEV,NX)  = ZW(NLEV)
        END DO
        
        CALL XZ_OUT('ZX',DOM%ZX,0_8)
        CALL XZ_OUT('XZ',DOM%XZ,0_8)
        CALL XZ_OUT('T',ZTXZ,0_8)
        CALL XZ_OUT('W',ZWXZ,0_8)
         
    END IF
     ! X-Y CROSS-SECTIONS
    IF (LSLICE_XY) THEN

       ZZXY(:,:) = ONE
       
       DO NX=1,GXPT
          DO NY=1,GYPT
             ZIJ   = DOM%ZS(NX,NY)
             IF (ZIJ < RZS) THEN
                ZZXY(NX,NY) = ZERO
             END IF
          END DO
       END DO

       CALL XY_OUT('YX',DOM%YX,0_8)
       CALL XY_OUT('XY',DOM%XY,0_8)
       CALL XY_OUT('ZXY',ZZXY,0_8)
          
    END IF
    
  END SUBROUTINE OUTPUT_INITIAL_FIELDS
  !************************************************************************************************************************!
  !************************************************************************************************************************!
  !************************************************************************************************************************!
  SUBROUTINE OUTPUT_CURRENT_FIELDS(RHS,PDTDX,PDQDX,PDPISDX,PDTDY,PDQDY,PDPISDY,STB,DOM,KSTEP)

    USE MOD_DIAGNOSTIC,  ONLY : DIAG_BACKGROUND_THERMO,DIAG_THERMO,DIAG_PHI, &
                              & DIAG_WPOS,DIAG_ALPHADELTA
    
    IMPLICIT NONE

    TYPE(GEODOM),                               INTENT(IN)    :: DOM
    TYPE(BACKGROUND),                           INTENT(IN)    :: STB
    TYPE(RHS_GLB),                              INTENT(IN)    :: RHS
    REAL(8), DIMENSION(NLEV,GXPT,GYPT),         INTENT(IN)    :: PDTDX,PDTDY
    REAL(8), DIMENSION(NLEV,GXPT,GYPT),         INTENT(IN)    :: PDQDX,PDQDY
    REAL(8), DIMENSION(GXPT,GYPT),              INTENT(IN)    :: PDPISDX,PDPISDY
    INTEGER(8),                                 INTENT(IN)    :: KSTEP
    
    INTEGER                                                   :: NX,NY,NYS
    INTEGER                                                   :: J0,J,K,K0,I0
    INTEGER(8), DIMENSION(NLEV)                               :: JLOC

    REAL(8)                                                   :: ZGAMMA, ZSECUR, ZDIF, PIJ, PLAST
    REAL(8), DIMENSION(NLEV)                                  :: ZPIF,ZALPHA,ZDELTA,ZDELTAPI
    REAL(8), DIMENSION(NLEV)                                  :: ZPHIF,ZTH,ZTH0,ZNUL,ZU,ZV
    REAL(8), DIMENSION(0:NLEV)                                :: ZPHIH,ZW
    
    REAL(8), DIMENSION(GXPT,GYPT)                             :: ZUXY, ZVXY, ZWXY
    REAL(8), DIMENSION(NLEV,GXPT)                             :: ZTXZ, ZWXZ

    LOGICAL, DIMENSION(NLEV)                                  :: LLDONE, LLEXTRA

    ZNUL(:)    = ZERO
    ZSECUR     = ABS(1000*EPSILON(ONE))
    NYS        = INT(RYS/RDY)+1
    
    ! X-Z CROSS-SECTIONS
    IF (LSLICE_XZ) THEN
        
       DO NX = 1, GXPT
   
          CALL DIAG_ALPHADELTA(ZPIF,ZALPHA,ZDELTA,ZDELTAPI,RHS%PIS(NX,NYS),DOM%AH,DOM%BH)
           
          ! INTERPOLATION POTENTIAL TEMPERATURE
          CALL DIAG_THERMO(RHS%T(:,NX,NYS),RHS%Q(:,NX,NYS),ZPIF,ZTH,'TH')
          CALL DIAG_BACKGROUND_THERMO(STB%T,STB%PIS,DOM%AH,DOM%BH,ZTH0,'TH')
          CALL DIAG_PHI(1_8,ZPHIF,RHS%T(:,NX,NYS),RHS%Q(:,NX,NYS),ZALPHA,ZDELTA,DOM%ZS(NX,NYS))
           
          PLAST = ZPHIF(NLEV)
          DO J= 1, NLEV-1
             LLDONE(J) = .FALSE.
             PIJ       = RG*RDZ*REAL(NLEV-J,8) 
             JLOC(J)   = J
             DO K = 1, NLEV-1
                ZDIF = (ZPHIF(K)-PIJ)*(PIJ-ZPHIF(K+1))        
                IF (ZDIF > -ZSECUR) THEN
                  JLOC(J)   = K
                  LLDONE(J) = .TRUE.
                  EXIT
                END IF 
             END DO   
             IF (LLDONE(J)) THEN
                LLEXTRA(J) = .FALSE.
             ELSE
                LLEXTRA(J) = (PLAST > PIJ)
             END IF  
             IF (.NOT.LLEXTRA(J)) THEN  
                ZGAMMA      = (ZPHIF(JLOC(J))-PIJ)/(ZPHIF(JLOC(J))-ZPHIF(JLOC(J)+1))
                ZTXZ(J,NX)  = ZGAMMA*(ZTH(JLOC(J)+1)-ZTH0(JLOC(J)+1)) &
                            & +(ONE-ZGAMMA)*(ZTH(JLOC(J))-ZTH0(JLOC(J)))
             ELSE
                ZTXZ(J,NX)  = ZTH(NLEV)-ZTH0(NLEV)
             END IF
          END DO
          ZTXZ(NLEV,NX)     = ZTH(NLEV)-ZTH0(NLEV)

          ! INTERPOLATION VERTICAL WIND
          CALL DIAG_WPOS(ZW,RHS%DW(:,NX,NYS),RHS%U(:,NX,NYS),RHS%V(:,NX,NYS),RHS%T(:,NX,NYS),&
               & PDTDX(:,NX,NYS),PDTDY(:,NX,NYS),RHS%Q(:,NX,NYS),PDQDX(:,NX,NYS),   &
               & PDQDY(:,NX,NYS),RHS%PIS(NX,NYS),PDPISDX(NX,NYS),PDPISDY(NX,NYS),   &
               & ZPIF,ZALPHA,ZDELTA,ZDELTAPI,DOM%AH,DOM%BH,DOM%SH,DOM%DZSDX(NX,NYS),         &
               & DOM%DZSDY(NX,NYS),KSTEP)
          CALL DIAG_PHI(0_8,ZPHIH,RHS%T(:,NX,NYS),RHS%PIS(NX,NYS),ZALPHA,ZDELTA,DOM%ZS(NX,NYS))

          PLAST = ZPHIH(NLEV)
          DO J = 1,NLEV-1
            LLDONE(J) = .FALSE.
            PIJ       = RG*RDZ*REAL(NLEV-J,8)
            JLOC(J)   = J
            DO K = 1, NLEV
              ZDIF = (ZPHIH(K-1)-PIJ)*(PIJ-ZPHIH(K))
              IF (ZDIF > -ZSECUR) THEN
                JLOC(J)   = K-1
                LLDONE(J) = .TRUE.
                EXIT
              END IF 
            END DO   
            IF (LLDONE(J)) THEN
               LLEXTRA(J) = .FALSE.
            ELSE
               LLEXTRA(J) = (PLAST > PIJ)
            END IF
            IF (.NOT.LLEXTRA(J)) THEN
               ZGAMMA     = (ZPHIH(JLOC(J))-PIJ)/(ZPHIH(JLOC(J))-ZPHIH(JLOC(J)+1))
               ZWXZ(J,NX) = ZGAMMA*ZW(JLOC(J)+1)+(ONE-ZGAMMA)*ZW(JLOC(J)) 
            ELSE
               ZWXZ(J,NX) = ZW(NLEV)
            END IF
          END DO      
           ZWXZ(NLEV,NX)  = ZW(NLEV)
        END DO
        
        CALL XZ_OUT('T',ZTXZ,KSTEP)
        CALL XZ_OUT('W',ZWXZ,KSTEP)
         
    END IF
     ! X-Y CROSS-SECTIONS
    IF (LSLICE_XY) THEN

         DO NX=1,GXPT
            DO NY=1,GYPT

               CALL DIAG_ALPHADELTA(ZPIF,ZALPHA,ZDELTA,ZDELTAPI,RHS%PIS(NX,NY),DOM%AH,DOM%BH)
               CALL DIAG_WPOS(ZW,RHS%DW(:,NX,NY),RHS%U(:,NX,NY),RHS%V(:,NX,NY),                   &
                    & RHS%T(:,NX,NY),PDTDX(:,NX,NY),PDTDY(:,NX,NY),RHS%Q(:,NX,NY),PDQDX(:,NX,NY), &
                    & PDQDY(:,NX,NY),RHS%PIS(NX,NY),PDPISDX(NX,NY),PDPISDY(NX,NY),                &
                    & ZPIF,ZALPHA,ZDELTA,ZDELTAPI,DOM%AH,DOM%BH,DOM%SH,DOM%DZSDX(NX,NY),          &
                    & DOM%DZSDY(NX,NY),KSTEP)
               CALL DIAG_PHI(0_8,ZPHIH,RHS%T(:,NX,NY),RHS%Q(:,NX,NY),ZALPHA,ZDELTA,DOM%ZS(NX,NY))
             
               LLDONE(1) = .FALSE.
               PLAST     = ZPHIH(NLEV)
               PIJ       = RG*RZS
               K0        = NLEV
               DO K = NLEV,1,-1
                  ZDIF = (ZPHIH(K-1)-PIJ)*(PIJ-ZPHIH(K))        
                  IF (ZDIF > -ZSECUR) THEN
                     K0 = K
                     LLDONE(1) = .TRUE.
                     EXIT
                  END IF 
               END DO   
               IF (LLDONE(1)) THEN
                  LLEXTRA(1) = .FALSE.
               ELSE
                  LLEXTRA(1) = (PLAST > PIJ)
               END IF  
               IF (.NOT.LLEXTRA(1)) THEN  
                  ZGAMMA       = (ZPHIH(K0-1)-PIJ)/(ZPHIH(K0-1)-ZPHIH(K0))
                  ZUXY(NX,NY)  = ZGAMMA*RHS%U(K0,NX,NY) + (ONE-ZGAMMA)*RHS%U(K0-1,NX,NY)
                  ZVXY(NX,NY)  = ZGAMMA*RHS%V(K0,NX,NY) + (ONE-ZGAMMA)*RHS%V(K0-1,NX,NY)
                  ZWXY(NX,NY)  = ZGAMMA*ZW(K0)          + (ONE-ZGAMMA)*ZW(K0-1)
               ELSE
                  ZUXY(NX,NY)  = RHS%U(NLEV,NX,NY)
                  ZVXY(NX,NY)  = RHS%V(NLEV,NX,NY)
                  ZWXY(NX,NY)  = ZW(NLEV) 
               END IF
            END DO
         END DO

         CALL XY_OUT('U',ZUXY,KSTEP)
         CALL XY_OUT('V',ZVXY,KSTEP)
         CALL XY_OUT('W',ZWXY,KSTEP)
         
    END IF   
    
  END SUBROUTINE OUTPUT_CURRENT_FIELDS
  !************************************************************************************************************************!
  !************************************************************************************************************************!
  !************************************************************************************************************************!
  SUBROUTINE XZ_OUT(NOMIN,XIN,ITERATIONIN)

    IMPLICIT NONE

    CHARACTER (LEN=*),               INTENT(IN) :: NOMIN
    INTEGER(8),                      INTENT(IN) :: ITERATIONIN
    REAL(8), DIMENSION(NLEV,GXPT),   INTENT(IN) :: XIN

    INTEGER(8)                                  :: N,L
    CHARACTER(LEN=11)                           :: CLOOP
    CHARACTER(LEN=20)                           :: CNAME

    WRITE(CLOOP,'(I11.11)') ITERATIONIN
    OPEN(UNIT=40,FILE='./res/XZ/'//NOMIN//CLOOP)

    DO L = 1, NLEV
       DO N = 1, GXPT-1
          WRITE(40,'(E18.8E4,1X)',ADVANCE='NO') XIN(L,N)
       END DO
       WRITE(40,'(E18.8E4)') XIN(L,GXPT)
    END DO

    CLOSE(40)

  END SUBROUTINE XZ_OUT
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE XY_OUT(NOMIN,XIN,ITERATIONIN)

    IMPLICIT NONE

    CHARACTER (LEN=*),             INTENT(IN) :: NOMIN
    INTEGER(8),                    INTENT(IN) :: ITERATIONIN
    REAL(8), DIMENSION(GXPT,GYPT), INTENT(IN) :: XIN
    INTEGER(8)                                :: NX,NY
    CHARACTER(LEN=11)                         :: CLOOP
    CHARACTER(LEN=20)                         :: CNAME

    WRITE(CLOOP,'(I11.11)') ITERATIONIN
    OPEN(UNIT=60,FILE='./res/XY/'//NOMIN//CLOOP)

    DO NY = 1, GYPT
       DO NX = 1, GXPT-1
          WRITE(60,'(E18.8E4,1X)',ADVANCE='NO') XIN(NX,NY)
       END DO
       WRITE(60,'(E18.8E4)') XIN(GXPT,NY)
    END DO

    CLOSE(60)

  END SUBROUTINE XY_OUT
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE MOD_VIEW

!=================================================================================!
!=================================================================================!
!=================================================================================!
