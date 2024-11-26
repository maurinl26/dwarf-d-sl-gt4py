!###############################  MODELE JOUET 2D  ###############################!
!#                                                                               #!
!# auteur : F.Voitus                                                             #!
!# sujet  : Gere l'ensemble des fonctions d'affichage et nettoyage               #!
!#                                                                               #!
!#################################################################################!

!=================================================================================!
!=================================  MODULE VIEW  =================================!
!=================================================================================!

MODULE MOD_VIEW

  USE MOD_SHARE, ONLY : NXPT,NLEV,NTIMESTEP,RDT,RDX,RDZ,RHMAX,RA,C_ADVECT,RG,LSLAG,NORD2,&
                      & LFCT_MONO,NCASE,NITMP,ZERO,ONE,TWO,LSETTLS,LADV_PIS,NTIMESTEP,   &
                      & NTOP,RKAPPA,RG,RP00,RPSUR,NLEFT,NRIGHT,RXC,LPOST_Z,LNESC,NORD1,  &
                      & LCOMAD,NCOMAD_OPT,LSLTRAJ_PREC,LADV_SPLIT,C_SLTRAJ,C_SLINTP,LSPEC
  USE MOD_SETUP, ONLY : BACKGROUND,GEOMETRY,STATEVAR
  
CONTAINS

  !*******************************************************************************!
  !*****************************  LISTE DES ROUTINES  ****************************!
  !*******************************************************************************!

  !**********************  Affichage en-tête et nettoyage   **********************!
  SUBROUTINE SCREENFIRST(PTIME,PDAY,GEO,STB)

    USE MOD_DIAG,  ONLY : DEF_Z,DEF_P

    IMPLICIT NONE

    INTEGER(4),        INTENT(OUT)  :: PTIME,PDAY
    TYPE(BACKGROUND),  INTENT(IN)   :: STB
    TYPE(GEOMETRY),    INTENT(IN)   :: GEO
    
    INTEGER(8)                      :: I,L,NX
    INTEGER(4), DIMENSION(3)        :: JDAYSTART,JNOW
    REAL(8), DIMENSION(0:NLEV,NXPT) :: ZH
    REAL(8), DIMENSION(NLEV,NXPT)   :: ZZ,ZP
    REAL(8), DIMENSION(NLEV)        :: ZNUL
    
    CALL IDATE(JDAYSTART)
    PDAY = JDAYSTART(2)

    WRITE(*,*) '     **************************'
    WRITE(*,*) "     *  2D SLICE TRACER TEST  *"      
    WRITE(*,*) "     *    OVER A MOUNTAIN     *"
    WRITE(*,*) '     **************************'
    WRITE(*,*) ' '
    IF (JDAYSTART(2) .LT. 10) THEN
       IF (JDAYSTART(1) .LT. 10) THEN
          WRITE(*,'(A11,I1,A2,I1,A1,I4)') 'date    : ',   &
          & JDAYSTART(1),'/0',JDAYSTART(2),'/',JDAYSTART(3)
       ELSE
          WRITE(*,'(A11,I2,A2,I1,A1,I4)') 'date    : ',   &
          & JDAYSTART(1),'/0',JDAYSTART(2),'/',JDAYSTART(3)
       END IF
    ELSE
       IF (JDAYSTART(1) .LT. 10) THEN
          WRITE(*,'(A11,I1,A1,I2,A1,I4)') 'date    : ',  &
          & JDAYSTART(1),'/',JDAYSTART(2),'/',JDAYSTART(3)
       ELSE
          WRITE(*,'(A11,I2,A1,I2,A1,I4)') 'date    : ',  &
          & JDAYSTART(1),'/',JDAYSTART(2),'/',JDAYSTART(3)
       END IF
    END IF

    WRITE(*,*)
    WRITE(*,*)
    WRITE(*,'(A11,1X,I2)') 'NCASE      : ',NCASE
    WRITE(*,'(A11,1X,L2)') 'LSPEC      : ',LSPEC
    WRITE(*,*)
    WRITE(*,*)
    WRITE(*,*) '* ------------------------- *'
    WRITE(*,*) '* SMILAG TRANSPORT SCHEME :  '
    WRITE(*,*) '* ------------------------- *'
    WRITE(*,*)
    WRITE(*,'(A11,1X,L2)') 'LSLAG      : ',LSLAG
    WRITE(*,'(A11,1X,I2)') 'NITMP      : ',NITMP
    WRITE(*,'(A11,1X,L2)') 'LSETTLS    : ',LSETTLS
    WRITE(*,'(A11,1X,L2)') 'LNESC      : ',LNESC
    WRITE(*,'(A11,1X,L2)') 'LPIS2D     : ',LADV_PIS
    WRITE(*,'(A11,1X,A2)') 'C_SLTRAJ   : ',C_SLTRAJ
    WRITE(*,'(A11,1X,A2)') 'C_SLINTP   : ',C_SLINTP
    WRITE(*,'(A11,1X,L2)') 'LSLTR_PREC : ',LSLTRAJ_PREC
    WRITE(*,'(A11,1X,L2)') 'LCOMAD     : ',LCOMAD
    WRITE(*,'(A11,1X,I2)') 'NCOMAD_OPT : ',NCOMAD_OPT
    WRITE(*,*)
    WRITE(*,*) '* ------------------------- *'
    WRITE(*,*) '* MPDATA TRANSPORT SCHEME :  '
    WRITE(*,*) '* ------------------------- *'
    WRITE(*,*)
    WRITE(*,'(A11,1X,L2)') 'LMPDATA_SLP: ',LADV_SPLIT
    WRITE(*,'(A11,1X,I2)') 'NORD1      : ',NORD1
    WRITE(*,'(A11,1X,I2)') 'NORD2      : ',NORD2
    WRITE(*,'(A11,1X,L2)') 'LFCTM      : ',LFCT_MONO
    WRITE(*,*)
    WRITE(*,*)
    

    NX = INT(NXPT/2)
    
    DO I=1,NXPT
       CALL DEF_Z(1,ZZ(:,I),STB%T,STB%PIS(NLEV,I), &
            & GEO%AH,GEO%BH,GEO%ZS(NLEV,I))
       CALL DEF_Z(0,ZH(:,I),STB%T,STB%PIS(NLEV,I), &
            & GEO%AH,GEO%BH,GEO%ZS(NLEV,I))
       CALL DEF_P(ZP(:,I),STB%PIS(NLEV,I),GEO%AH,GEO%BH,1)
    END DO
       
    WRITE(*,'(8X,A5,7X,A7,13X,A7,7X,A9,4X,A15,2X,A12,2X,A12,2X,A12)') 'level',       &     
         & 'A(eta)','B(eta)','temp (K)',' pression (hPa) ','U-prof (m/s)','Z-full (km)','DeltaZ (m)'
    DO L = 1, NLEV
       WRITE(*,'(8X,I3,2X,F18.12,2X,F18.12,3X,F8.2,3X,F18.12,4X,F6.3,4X,F6.3)') &
            & L,GEO%AH(L),GEO%BH(L),STB%T(L,NX),ZP(L,NX)*0.01d0,&
            & STB%U(L,NX),ZZ(L,NX)*0.001d0,ZH(L-1,NX)-ZH(L,NX)
    END DO
    WRITE(*,*)
    WRITE(*,'(A24,5X,A24)') ' DIMENSIONS DU DOMAINE :', ' PARAMETRES DYNAMIQUES :' 
    WRITE(*,*)
    WRITE(*,'(A7,I6,16X,A5,F6.2,A2)')'STEP = ',NTIMESTEP,'dt = ', RDT,' s'
    WRITE(*,'(A7,I6,16X,A5,F6.0,A2)')'NXPT = ', NXPT,'dx = ', RDX,' m'
    WRITE(*,'(A7,I6,16X,A5,F6.0,A2)')'NLEV = ', NLEV,'dz = ', RDZ,' m'
    WRITE(*,'(A7,F6.0,A2)')'HMAX = ', RHMAX, ' m'
    WRITE(*,'(A7,F6.0,A2)')'RA   = ', RA,' m'
    WRITE(*,'(A7,I6,A3)')  'Lx   = ', INT(RDX*NXPT/1000),' km'
    WRITE(*,'(A7,I6,A3)')  'Lz   = ', INT(0.001d0*ZH(0,INT(NXPT/2))),' km'
    WRITE(*,*)  
    WRITE(*,'(A7,I3,A2,I2,A4,I2,A1)') 'TIME : ',INT(NTIMESTEP*RDT)/3600,'h ',  &
         & (INT(NTIMESTEP*RDT)-(INT(NTIMESTEP*RDT)/3600)*3600)/60,'min ',      & 
         & INT(NTIMESTEP*RDT)-INT(NTIMESTEP*RDT)/3600*3600-(INT(NTIMESTEP*RDT) &
         & -INT(NTIMESTEP*RDT)/3600*3600)/60*60,'s'

    WRITE(*,*) 
 
    CALL ITIME(JNOW)
    WRITE(*,'(A20,I2,A1,I2,A1,I2)') '  début du calcul : ', &
         & JNOW(1),':',JNOW(2),':',JNOW(3)
    PTIME = JNOW(1)*3600+JNOW(2)*60+JNOW(3)

    WRITE(*,*)
    WRITE(*,'(A54)') " 0%  10%  20%  30%  40%  50%  60%  70%  80%  90%  100%"
    WRITE(*,'(A52)') " |----|----|----|----|----|----|----|----|----|----|"
    WRITE(*,'(A2)',ADVANCE='NO') ' *'

  END SUBROUTINE SCREENFIRST
  !*******************************************************************************!
  !***********************  Affichage durant la dynamique   **********************!
  !*******************************************************************************!
  SUBROUTINE SCREENDYN(TLOOP,ST,PTIME)

    IMPLICIT NONE

    INTEGER(4),     INTENT(INOUT) :: PTIME
    INTEGER(8),     INTENT(IN)    :: TLOOP
    TYPE(STATEVAR), INTENT(IN)    :: ST

    REAL(8)                       :: PU
    INTEGER(4), DIMENSION(3)      :: INOW,ITODAY
    CHARACTER(len=3)              :: STARS

    IF (PU .GT. 1.d4) THEN
       WRITE(*,*) "\n          SCHEMA INSTABLE" 
       WRITE(*,*) "     Calculs stoppes au rang : ",TLOOP
       CALL ITIME(INOW)
       WRITE(*,*)
       WRITE(*,'(A20,I2,A1,I2,A1,I2)') "    fin du calcul : ", &
            & INOW(1),":",INOW(2),':',INOW(3)
       CALL IDATE(ITODAY)
       PTIME = INOW(1)*3600+INOW(2)*60+INOW(3)-PTIME
       INOW(1) = PTIME/3600
       INOW(2) = (PTIME-INOW(1)*3600)/60
       INOW(3) = PTIME-INOW(1)*3600-INOW(2)*60
       WRITE(*,'(A20,I2,A2,I2,A4,I2,A1)') "  temps de calcul : ", &
            & INOW(1),"h ",INOW(2),"min ",INOW(3),'s'
       WRITE(*,*) "\n Fin du programme \n"
       STOP
    END IF

    IF ((MOD(50*TLOOP,NTIMESTEP) == 0) .AND. (TLOOP /= 0)) THEN
       WRITE(*,'(A1)',ADVANCE='NO') "*"
    END IF
    
  END SUBROUTINE SCREENDYN
  !*******************************************************************************!
  !*****************************   Affichage finale  *****************************!
  !*******************************************************************************!
  SUBROUTINE SCREENFINAL(PTIME,PDAY)

    USE MOD_SHARE, ONLY : NTIMESTEP

    IMPLICIT NONE

    INTEGER(4), INTENT(INOUT) :: PTIME
    INTEGER(4), INTENT(IN)    :: PDAY
    INTEGER(4), DIMENSION(3)  :: INOW,ITODAY

    WRITE(*,*)
    CALL ITIME(INOW)
    WRITE(*,*)
    WRITE(*,'(A20,I5)')             "    nbr timesteps : ",NTIMESTEP
    WRITE(*,'(A20,I2,A1,I2,A1,I2)') "    fin du calcul : ", &
         & INOW(1),":",INOW(2),':',INOW(3)
    CALL IDATE(ITODAY)
    IF (PDAY == ITODAY(2)) THEN
       PTIME = INOW(1)*3600+INOW(2)*60+INOW(3)-PTIME
       INOW(1) = PTIME/3600
       INOW(2) = (PTIME-INOW(1)*3600)/60
       INOW(3) = PTIME-INOW(1)*3600-INOW(2)*60
       WRITE(*,'(A20,I2,A2,I2,A4,I2,A1)') "  temps de calcul : ", &
            & INOW(1),"h ",INOW(2),"min ",INOW(3),"s"
    END IF
    WRITE(*,*) '\n Fin du programme \n '

  END SUBROUTINE SCREENFINAL
  !*******************************************************************************!
  !***********************  Imprime l'etat des variables  ************************!
  !*******************************************************************************!
  SUBROUTINE PRINT_STATE(ST,STB,GEO,NLOOP)

    USE MOD_DIAG,  ONLY : DEF_Z
    
    IMPLICIT NONE

    INTEGER(8),               INTENT(IN) :: NLOOP
    TYPE(STATEVAR),           INTENT(IN) :: ST
    TYPE(BACKGROUND),         INTENT(IN) :: STB
    TYPE(GEOMETRY),           INTENT(IN) :: GEO
    
    INTEGER(8)                           :: I,J,K
    INTEGER(8), DIMENSION(NLEV)          :: ILOC
    REAL(8),    DIMENSION(NLEV,NXPT)     :: ZZ,ZQ
    REAL(8)                              :: ZGAM,ZSECUR,ZDIFF
    REAL(8)                              :: ZIJ, ZLAST, ZMIN, ZTOP
    LOGICAL,    DIMENSION(NLEV)          :: LLDONE, LLEXTRA
    
    
    !************************************************************************************!   
    IF (LPOST_Z) THEN

       DO I=1,NXPT
          CALL DEF_Z(1,ZZ(:,I),STB%T(:,I), &
               & ST%PIS(NLEV,I),GEO%AH,GEO%BH,GEO%ZS(NLEV,I))
       END DO
    
       ZSECUR  =  ABS(1000.d0*EPSILON(1.d0))
       ZMIN    =  MINVAL(GEO%ZS(NLEV,:))

       DO I = 1, NXPT
          ZLAST = ZZ(NLEV,I)
          ZTOP  = ZZ(1,I)
          DO J= 1, NLEV-1
             LLDONE(J) = .FALSE.
             ZIJ       = ZMIN + RDZ*REAL(NLEV-J,8)
             ILOC(J)   = J
             DO K = 1, NLEV-1
                ZDIFF  = (ZZ(K,I)-ZIJ)*(ZIJ-ZZ(K+1,I))        
                IF (ZDIFF > -ZSECUR) THEN
                   ILOC(J)   = K
                   LLDONE(J) = .TRUE.
                   EXIT
                END IF 
             END DO
             IF (LLDONE(J)) THEN
                LLEXTRA(J) = .FALSE.
             ELSE
                LLEXTRA(J) = (ZLAST > ZIJ)
             END IF  
             IF (.NOT.LLEXTRA(J)) THEN  
                ZGAM    = (ZZ(ILOC(J),I)-ZIJ)/(ZZ(ILOC(J),I)-ZZ(ILOC(J)+1,I))
                ZQ(J,I) = ZGAM*ST%Q(ILOC(J)+1,I)+(ONE-ZGAM)*ST%Q(ILOC(J),I)
             ELSE
                ZQ(J,I) = ST%Q(NLEV,I)
             END IF
          END DO
          ZQ(NLEV,I)    = ST%Q(NLEV,I)
       END DO
    ELSE
       ZQ(:,:) = ST%Q(:,:)
    END IF
    
    CALL PRINT_FIELD_UNIT('Q',ZQ,NLOOP)

  END SUBROUTINE PRINT_STATE
  !*******************************************************************************!
  !******************  Fonction qui cree les fichiers de sortie  *****************!
  !*******************************************************************************!
  SUBROUTINE PRINT_FIELD_UNIT(CPREFI,V,NLOOP)

    IMPLICIT NONE

    CHARACTER (LEN=*),                    INTENT(IN) :: CPREFI
    INTEGER(8),                           INTENT(IN) :: NLOOP
    REAL(8), DIMENSION(NTOP:NLEV,NXPT),   INTENT(IN) :: V
    INTEGER(8)                                       :: I,J
    CHARACTER(LEN=11)                                :: CLOOP
    CHARACTER(LEN=20)                                :: CNAME

    WRITE(CLOOP,'(I11.11)') NLOOP
    CNAME = CPREFI//CLOOP
    OPEN(UNIT=40,FILE='./res/'//CNAME)

    DO J = 1,NLEV
       DO I = NLEFT, NRIGHT-1
          WRITE(40,'(E18.8E4,1X)',ADVANCE='NO') V(J,I)
       END DO
       WRITE(40,'(E18.8E4)') V(J,NRIGHT)
    END DO

    CLOSE(40)

  END SUBROUTINE PRINT_FIELD_UNIT
  !*******************************************************************************!

END MODULE MOD_VIEW

!=================================================================================!
!=================================================================================!
!=================================================================================!
