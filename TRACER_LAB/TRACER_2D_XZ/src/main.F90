!#################################  MODELE JOUET 2D  ##################################!
!#                                                                                    #!
!# auteurs : F,Voitus                                                                 #!
!#                                                                                    #!
!######################################################################################!

!=====================================================================!
!=======================  PROGRAMME PRINCIPAL  =======================!
!=====================================================================!

PROGRAM MAIN

  USE MOD_SHARE
  USE MOD_SETUP
  USE MOD_STATE
  USE MOD_SLAG
  USE MOD_MPDTA
  USE MOD_VIEW

  IMPLICIT NONE

  INTEGER(4)                    :: TIME,DAY
  INTEGER(8)                    :: NSTEP
  TYPE(GEOMETRY)                :: GEO
  TYPE(SMILAG_STRUCT)           :: YSLAG
  TYPE(MPDATA_STRUCT)           :: YMDTA
  TYPE(BACKGROUND)              :: STB
  TYPE(RHSVAR)                  :: RHS_A ,RHS_D ,RHS_M   
  TYPE(STATEVAR)                :: ST9   ,ST0   ,ST1
  

  WRITE(*,*) "   Initialisation des parameteres......"
  CALL INITIAL_SETTING
  WRITE(*,*) "...Création des condition initiales...."
  CALL INITIAL_STATE(ST0,ST9,STB,GEO)
  
  !WRITE(*,*) "...Ecriture des condition initiales...."
  CALL PRINT_FIELD_UNIT('Z',GEO%Z(NTOP:NLEV,:)/1000.d0,0_8)
  CALL PRINT_FIELD_UNIT('X',(GEO%X(NTOP:NLEV,:)-RXC)/1000.d0,0_8)
  CALL PRINT_STATE(ST0,STB,GEO,0_8)
  !WRITE(*,*)
  WRITE(*,*) " "
  WRITE(*,*) "   Integration........................."
  WRITE(*,*) "  "
  CALL SCREENFIRST(TIME,DAY,GEO,STB)
  ! WRITE(*,*) " "
  
  DO NSTEP = 1, NTIMESTEP
     
     CALL RHS_MDL(RHS_M,ST0,GEO)
     
     CALL RHS_EXP(RHS_A,ST0,RDT/TWO,RHS_M,GEO)
     CALL SWAPP_VAR(ST1,ST0,LWIND=.TRUE.)

     IF (LSLAG) THEN
        CALL SMILAG_TRACER_TRANSPORT_SCHEME(RHS_D,RHS_A,YSLAG,  &
             & ST1%U,ST0%U,ST9%U,ST1%ETADOT,ST0%ETADOT,         &
             & ST9%ETADOT,ST1%DIV,ST0%DIV,GEO%ETA,GEO%DELB,     &
             & NSTEP,C_ADVECT)        
     ELSE
        CALL MPDATA_TRACER_TRANSPORT_SCHEME(RHS_D,RHS_A,YMDTA,  &
             & ST0%M,ST1%U,ST0%U,ST9%U,ST1%ETADOTH,ST0%ETADOTH,&
             & ST9%ETADOTH,GEO%DETA,GEO%DELB,NSTEP,C_ADVECT)
     ENDIF

     CALL SWAPP_VAR(ST9,ST0)
     CALL RHS_IMP(ST0,ST1,STB,RHS_D,GEO)
          
     ! ecriture des fichiers de sorties
     IF (MOD(NSTEP,NFREQOUTPUT) == 0) THEN
        CALL PRINT_STATE(ST0,STB,GEO,NSTEP)
     END IF
     CALL SCREENDYN(NSTEP,ST0,TIME)
     
  END DO
  
  IF (MOD(NTIMESTEP,NFREQOUTPUT) .NE. 0) THEN
    CALL PRINT_STATE(ST0,STB,GEO,NTIMESTEP)
  END IF  
  CALL SCREENFINAL(TIME,DAY)

END PROGRAM

!=====================================================================!
!=====================================================================!
!=====================================================================!
