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

  USE MOD_PARAM, ONLY : GXPT,GYPT,RLY,RLX,RDX,RDY,RG,RPI,NTIMESTEP, &
                      & ZERO,ONE,TWO,HALF,RDT
  USE MOD_SETUP, ONLY : DOMAIN,PROGVAR
  
CONTAINS

  !*****************************************************!
  !***************  LISTE DES ROUTINES  ****************!
  !*****************************************************!

  !*********  Affichage en-tête et nettoyage   *********!
  SUBROUTINE SCREEN_START(TEMPSOUT,JOUROUT,DOM,ST0)

    IMPLICIT NONE

    INTEGER,           INTENT(OUT) :: TEMPSOUT,JOUROUT
    TYPE(PROGVAR),     INTENT(IN)  :: ST0
    TYPE(DOMAIN),      INTENT(IN)  :: DOM
    
    INTEGER                        :: L
    INTEGER, DIMENSION(3)          :: JJOUR,JMAINTENANT

    CALL IDATE(JJOUR)
    JOUROUT = JJOUR(2)
    
    
    WRITE(*,*) "     ***********************"
    WRITE(*,*) "     *   TRACER XY MODEL   *"
    WRITE(*,*) "     ***********************"
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
    WRITE(*,'(1X,A7,F6.2,A2)')'dt   = ', RDT, ' s'
    WRITE(*,'(1X,A7,I6,A2)')  'dx   = ', INT(RDX),' m'
    WRITE(*,'(1X,A7,I6,A2)')  'dy   = ', INT(RDY),' m'
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
  SUBROUTINE SCREEN(TEMPSINOUT,ST,DOM,BOUCLEIN)

    IMPLICIT NONE

    INTEGER,                            INTENT(INOUT) :: TEMPSINOUT
    INTEGER(8),                         INTENT(IN)    :: BOUCLEIN
    TYPE(PROGVAR),                      INTENT(IN)    :: ST
    TYPE(DOMAIN),                       INTENT(IN)    :: DOM
    
    INTEGER, DIMENSION(3)                             :: IMAINTENANT,IJOUR
    CHARACTER(len=3)                                  :: STARS    
    REAL(8)                                           :: ZMAXU,ZMAXV,ZSEUIL
    LOGICAL                                           :: LL_UNSTABLE

    ZMAXU=MAXVAL(ABS(ST%U))
    ZMAXV=MAXVAL(ABS(ST%V))
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

    !WRITE(*,'(F20.16,2X,F20.16,2X,F10.6)') ZMAXU, ZMAXV
    
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
  SUBROUTINE OUTPUT_INITIAL_FIELDS(ST,DOM)
    
    IMPLICIT NONE

    TYPE(DOMAIN),                               INTENT(IN)    :: DOM
    TYPE(PROGVAR),                              INTENT(IN)    :: ST
    
    ! X-Y CROSS-SECTION

    CALL X_OUT('X',DOM%X,0_8)
    CALL Y_OUT('Y',DOM%Y,0_8)
    CALL XY_OUT('Q',ST%Q,0_8)
    CALL XY_OUT('U',ST%Q,0_8)
    CALL XY_OUT('V',ST%Q,0_8)
    
    
  END SUBROUTINE OUTPUT_INITIAL_FIELDS
  !**********************************************************************************************!
  !**********************************************************************************************!
  !**********************************************************************************************!
  SUBROUTINE OUTPUT_CURRENT_FIELDS(ST,DOM,KSTEP)

    USE MOD_DIAG,  ONLY : SP_NORMS
    
    IMPLICIT NONE

    TYPE(DOMAIN),                               INTENT(IN)    :: DOM
    TYPE(PROGVAR),                              INTENT(IN)    :: ST
    INTEGER(8),                                 INTENT(IN)    :: KSTEP
    
    CALL XY_OUT('U',ST%U,KSTEP)
    CALL XY_OUT('V',ST%V,KSTEP)
    CALL XY_OUT('Q',ST%Q,KSTEP)    
    
  END SUBROUTINE OUTPUT_CURRENT_FIELDS
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
       DO NX = 1, GXPT
          WRITE(60,'(E18.8E4,1X)',ADVANCE='NO') XIN(NX,NY)
       END DO
       !WRITE(60,'(E18.8E4)') XIN(GXPT,NY)
    END DO

    CLOSE(60)

  END SUBROUTINE XY_OUT
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE X_OUT(NOMIN,XIN,ITERATIONIN)

    IMPLICIT NONE

    CHARACTER (LEN=*),             INTENT(IN) :: NOMIN
    INTEGER(8),                    INTENT(IN) :: ITERATIONIN
    REAL(8), DIMENSION(GXPT),      INTENT(IN) :: XIN
    INTEGER(8)                                :: NX
    CHARACTER(LEN=11)                         :: CLOOP
    CHARACTER(LEN=20)                         :: CNAME

    WRITE(CLOOP,'(I11.11)') ITERATIONIN
    OPEN(UNIT=60,FILE='./res/XY/'//NOMIN//CLOOP)

    DO NX=1,GXPT
       WRITE(60,'(E18.8E4,1X)',ADVANCE='NO') XIN(NX)
    END DO

    CLOSE(60)

  END SUBROUTINE X_OUT
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE Y_OUT(NOMIN,XIN,ITERATIONIN)

    IMPLICIT NONE

    CHARACTER (LEN=*),             INTENT(IN) :: NOMIN
    INTEGER(8),                    INTENT(IN) :: ITERATIONIN
    REAL(8), DIMENSION(GYPT),      INTENT(IN) :: XIN
    INTEGER(8)                                :: NY
    CHARACTER(LEN=11)                         :: CLOOP
    CHARACTER(LEN=20)                         :: CNAME

    WRITE(CLOOP,'(I11.11)') ITERATIONIN
    OPEN(UNIT=60,FILE='./res/XY/'//NOMIN//CLOOP)

    DO NY = 1,GYPT
       WRITE(60,'(E18.8E4,1X)',ADVANCE='NO') XIN(NY)
    END DO

    CLOSE(60)

  END SUBROUTINE Y_OUT
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
END MODULE MOD_VIEW

!=================================================================================!
!=================================================================================!
!=================================================================================!
