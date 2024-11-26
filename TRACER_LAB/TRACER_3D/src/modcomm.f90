!##################   PARAMETRES   ###################!
!#                                                   #!
!# auteur : F.Voitus                                 #!
!# sujet  : Définit les parametres du modèle         #!
!#                                                   #!
!#####################################################!

!=====================================================!
!================  MODULE PARAMETRES  ================!
!=====================================================!

MODULE MOD_COMM

  USE MOD_PARAMETER, ONLY : GXPT,GYPT,NLEV,RPI,RDX,RDY,RLX,RLY,ZERO,ONE,TWO,HALF, &
                          & HALO,LXPT,LYPT,NDIM,TRANCHEX,TRANCHEY,NTAG,NB_PROCS,  &
                          & RR,RG,RT00,RP00
  USE MOD_STRUCTURE, ONLY : VARPROG,VARGLOB,GEOMETRY,DOMAIN
  
CONTAINS

  !*****************************************************!
  !***************  LISTE DES ROUTINES  ****************!
  !*****************************************************!

  !****************************************************************************************!
  !****************************************************************************************!
  !****************************************************************************************!
  SUBROUTINE TRANSGLB(GLOBALOUT,LOCALIN,TOPIN,ABS_GLB,ORD_GLB,MYPROC,CODE,NUM_FIELD)

    USE MPI

    IMPLICIT NONE
  
    REAL(8), DIMENSION(TOPIN:NLEV,GXPT,GYPT),         INTENT(OUT)   :: GLOBALOUT
    REAL(8), DIMENSION(TOPIN:NLEV,NDIM),              INTENT(IN)    :: LOCALIN
    INTEGER(8), DIMENSION(0:NB_PROCS-1),              INTENT(IN)    :: ABS_GLB,ORD_GLB
    INTEGER,                                          INTENT(IN)    :: TOPIN
    INTEGER,                                          INTENT(IN)    :: MYPROC
    INTEGER,                                          INTENT(IN)    :: NUM_FIELD
    INTEGER,                                          INTENT(INOUT) :: CODE
    
    INTEGER(8)                                                      :: JL,NY,NL
    INTEGER                                                         :: NPROC
    INTEGER, DIMENSION(MPI_STATUS_SIZE)                             :: STATUS

    DO JL = TOPIN, NLEV
       IF (MYPROC .NE. 0) THEN
          DO NY = 1, LYPT
             NL = (NY-1)*TRANCHEX+HALO
             CALL MPI_SEND(LOCALIN(JL,NL+1:NL+LXPT),LXPT,MPI_REAL8,0, &
                  & NTAG(NUM_FIELD,NY,MYPROC),MPI_COMM_WORLD,CODE)
          END DO
       ELSE
          DO NY = 1, LYPT
             NL = (NY-1)*TRANCHEX+HALO
             GLOBALOUT(JL,1:LXPT,NY) = LOCALIN(JL,NL+1:NL+LXPT)
             DO NPROC = 1, NB_PROCS-1
                CALL MPI_RECV(GLOBALOUT(JL,ABS_GLB(NPROC)+1:ABS_GLB(NPROC)+LXPT,ORD_GLB(NPROC)+NY),&
                     & LXPT,MPI_REAL8,NPROC,NTAG(NUM_FIELD,NY,NPROC),MPI_COMM_WORLD,STATUS,CODE)
             END DO
          END DO
       END IF
       CALL MPI_BARRIER(MPI_COMM_WORLD,CODE)
    END DO

  END SUBROUTINE TRANSGLB
  !****************************************************************************************!
  !****************************************************************************************!
  !****************************************************************************************!
  SUBROUTINE TRANSLOC(LOCALOUT,GLOBALIN,TOPIN,ABS_GLB,ORD_GLB,MYPROC,CODE,NUM_FIELD)

    USE MPI

    IMPLICIT NONE

    REAL(8),    DIMENSION(TOPIN:NLEV,NDIM),           INTENT(OUT)   :: LOCALOUT
    REAL(8),    DIMENSION(TOPIN:NLEV,GXPT,GYPT),      INTENT(IN)    :: GLOBALIN
    INTEGER(8), DIMENSION(0:NB_PROCS-1),              INTENT(IN)    :: ABS_GLB,ORD_GLB
    INTEGER,                                          INTENT(IN)    :: TOPIN
    INTEGER,                                          INTENT(IN)    :: MYPROC
    INTEGER,                                          INTENT(IN)    :: NUM_FIELD
    INTEGER,                                          INTENT(INOUT) :: CODE
    
    INTEGER(8)                                                      :: JL,NY,NL
    INTEGER                                                         :: NPROC
    INTEGER, DIMENSION(MPI_STATUS_SIZE)                             :: STATUS

    DO JL = TOPIN, NLEV
       IF (MYPROC .EQ. 0) THEN
          DO NY = 1, LYPT
             NL = (NY-1)*TRANCHEX+HALO
             LOCALOUT(JL,NL+1:NL+LXPT) = GLOBALIN(JL,1:LXPT,NY)
             DO NPROC = 1, NB_PROCS-1          
                CALL MPI_SEND(GLOBALIN(JL,ABS_GLB(NPROC)+1:ABS_GLB(NPROC)+LXPT,ORD_GLB(NPROC)+NY),&
                     & LXPT,MPI_REAL8,NPROC,NTAG(NUM_FIELD,NY,NPROC),MPI_COMM_WORLD,CODE)
             END DO
          END DO
       ELSE
          DO NY = 1, LYPT
             NL = (NY-1)*TRANCHEX+HALO
             CALL MPI_RECV(LOCALOUT(JL,NL+1:NL+LXPT),LXPT,MPI_REAL8,0, &
                  & NTAG(NUM_FIELD,NY,MYPROC),MPI_COMM_WORLD,STATUS,CODE)
          END DO
       END IF
       CALL MPI_BARRIER(MPI_COMM_WORLD,CODE)
    END DO

  END SUBROUTINE TRANSLOC
  !****************************************************************************************!
  !****************************************************************************************!
  !****************************************************************************************!
  SUBROUTINE TRANSLOC_SURF(PSI_L,DPSIDX_L,DPSIDY_L,PSI_G,DPSIDX_G,DPSIDY_G, &
                     & ABS_GLB,ORD_GLB,MYPROC,CODE,NUM)

    USE MPI

    IMPLICIT NONE

    REAL(8), DIMENSION(NDIM),              INTENT(OUT)   :: PSI_L, DPSIDX_L, DPSIDY_L 
    REAL(8), DIMENSION(GXPT,GYPT),         INTENT(IN)    :: PSI_G, DPSIDX_G, DPSIDY_G
    INTEGER(8), DIMENSION(0:NB_PROCS-1),   INTENT(IN)    :: ABS_GLB,ORD_GLB
    INTEGER,                               INTENT(IN)    :: MYPROC
    INTEGER,                               INTENT(IN)    :: NUM
    INTEGER,                               INTENT(INOUT) :: CODE
    
    INTEGER(8)                                           :: JL,NY,NL
    INTEGER                                              :: NPROC
    INTEGER, DIMENSION(MPI_STATUS_SIZE)                  :: STATUS

    IF (MYPROC .EQ. 0) THEN
       DO NY = 1, LYPT
          NL = (NY-1)*TRANCHEX+HALO
          PSI_L(NL+1:NL+LXPT)      = PSI_G(1:LXPT,NY)
          DPSIDX_L(NL+1:NL+LXPT)   = DPSIDX_G(1:LXPT,NY)
          DPSIDY_L(NL+1:NL+LXPT)   = DPSIDY_G(1:LXPT,NY)
          DO NPROC = 1, NB_PROCS-1
             CALL MPI_SEND(PSI_G(ABS_GLB(NPROC)+1:ABS_GLB(NPROC)+LXPT,ORD_GLB(NPROC)+NY),    &
                  & LXPT,MPI_REAL8,NPROC,NTAG(NUM,NY,NPROC),MPI_COMM_WORLD,CODE)
             CALL MPI_SEND(DPSIDX_G(ABS_GLB(NPROC)+1:ABS_GLB(NPROC)+LXPT,ORD_GLB(NPROC)+NY), &
                  & LXPT,MPI_REAL8,NPROC,NTAG(NUM+2,NY,NPROC),MPI_COMM_WORLD,CODE)
             CALL MPI_SEND(DPSIDY_G(ABS_GLB(NPROC)+1:ABS_GLB(NPROC)+LXPT,ORD_GLB(NPROC)+NY), &
                  & LXPT,MPI_REAL8,NPROC,NTAG(NUM+4,NY,NPROC),MPI_COMM_WORLD,CODE)
          END DO
       END DO   
    ELSE
       DO NY = 1, LYPT
          NL = (NY-1)*TRANCHEX+HALO
          CALL MPI_RECV(PSI_L(NL+1:NL+LXPT),LXPT,MPI_REAL8,0,NTAG(NUM,NY,MYPROC),      &
               & MPI_COMM_WORLD,STATUS,CODE)
          CALL MPI_RECV(DPSIDX_L(NL+1:NL+LXPT),LXPT,MPI_REAL8,0,NTAG(NUM+2,NY,MYPROC), &
               & MPI_COMM_WORLD,STATUS,CODE)
          CALL MPI_RECV(DPSIDY_L(NL+1:NL+LXPT),LXPT,MPI_REAL8,0,NTAG(NUM+4,NY,MYPROC), &
               & MPI_COMM_WORLD,STATUS,CODE)
       END DO
    END IF
    CALL MPI_BARRIER(MPI_COMM_WORLD,CODE)

  END SUBROUTINE TRANSLOC_SURF
  !******************************************************************************************!
  !******************************************************************************************!
  !******************************************************************************************!
  SUBROUTINE TRANSGLB_RHS(RHSG,RHSL,ABS_GLB,ORD_GLB,MYPROC,CODE,NUM_FIELD)

    USE MPI

    IMPLICIT NONE
  
    TYPE(RHS_GLB),                                    INTENT(INOUT) :: RHSG
    TYPE(RHS_LOC),                                    INTENT(IN)    :: RHSL
    INTEGER(8), DIMENSION(0:NB_PROCS-1),              INTENT(IN)    :: ABS_GLB,ORD_GLB
    INTEGER,                                          INTENT(IN)    :: MYPROC
    INTEGER,                                          INTENT(IN)    :: NUM_FIELD
    INTEGER,                                          INTENT(INOUT) :: CODE
    
    INTEGER(8)                                                      :: JL,NY,NL
    INTEGER                                                         :: NPROC
    INTEGER, DIMENSION(MPI_STATUS_SIZE)                             :: STATUS

    DO JL = 1,NLEV
       IF (MYPROC .NE. 0) THEN
          DO NY = 1, LYPT
             NL = (NY-1)*TRANCHEX+HALO
             CALL MPI_SEND(RHSL%T(JL,NL+1:NL+LXPT),LXPT,MPI_REAL8,0,   &
                  & NTAG(NUM_FIELD+1,NY,MYPROC),MPI_COMM_WORLD,CODE)
             CALL MPI_SEND(RHSL%Q(JL-1,NL+1:NL+LXPT),LXPT,MPI_REAL8,0, &
                  & NTAG(NUM_FIELD+2,NY,MYPROC),MPI_COMM_WORLD,CODE)
          END DO
       ELSE
          DO NY = 1, LYPT
             NL = (NY-1)*TRANCHEX+HALO
             RHSG%T(JL,1:LXPT,NY)   = RHSL%T(JL,NL+1:NL+LXPT)
             RHSG%Q(JL-1,1:LXPT,NY) = RHSL%Q(JL-1,NL+1:NL+LXPT)
             DO NPROC = 1, NB_PROCS-1
                CALL MPI_RECV(RHSB%T(JL,ABS_GLB(NPROC)+1:ABS_GLB(NPROC)+LXPT,ORD_GLB(NPROC)+NY),   & 
                     & LXPT,MPI_REAL8,NPROC,NTAG(NUM_FIELD+1,NY,NPROC),MPI_COMM_WORLD,STATUS,CODE)
                CALL MPI_RECV(RHSG%Q(JL-1,ABS_GLB(NPROC)+1:ABS_GLB(NPROC)+LXPT,ORD_GLB(NPROC)+NY), &
                     & LXPT,MPI_REAL8,NPROC,NTAG(NUM_FIELD+2,NY,NPROC),MPI_COMM_WORLD,STATUS,CODE)
             END DO
          END DO
       END IF
       CALL MPI_BARRIER(MPI_COMM_WORLD,CODE)
    END DO

    IF (MYPROC .NE. 0) THEN
       DO NY = 1, LYPT
          NL = (NY-1)*TRANCHEX+HALO
          CALL MPI_SEND(RHSL%Q(NLEV,NL+1:NL+LXPT),LXPT,MPI_REAL8,0, &
               & NTAG(NUM_FIELD+8,NY,MYPROC),MPI_COMM_WORLD,CODE)
          CALL MPI_SEND(RHSL%PIS(NL+1:NL+LXPT),LXPT,MPI_REAL8,0,    &
               & NTAG(NUM_FIELD+9,NY,MYPROC),MPI_COMM_WORLD,CODE)
       END DO
    ELSE
       DO NY = 1, LYPT
          NL = (NY-1)*TRANCHEX+HALO
          RHSG%Q(NLEV,1:LXPT,NY) = RHSL%Q(NLEV,NL+1:NL+LXPT)
          RHSG%PIS(1:LXPT,NY)    = RHSL%PIS(NL+1:NL+LXPT)
          DO NPROC = 1, NB_PROCS-1
             CALL MPI_RECV(RHSG%Q(NLEV,ABS_GLOBAL(NPROC)+1:ABS_GLOBAL(NPROC)+LXPT,ORD_GLOBAL(NPROC)+NY),  &
                  & LXPT,MPI_REAL8,NPROC,NTAG(NUM_FIELD+8,NY,NPROC),MPI_COMM_WORLD,STATUS,CODE)
             CALL MPI_RECV(RHSG%PIS(ABS_GLOBAL(NPROC)+1:ABS_GLOBAL(NPROC)+LXPT,ORD_GLOBAL(NPROC)+NY),     &
                  & LXPT,MPI_REAL8,NPROC,NTAG(NUM_FIELD+9,NY,NPROC),MPI_COMM_WORLD,STATUS,CODE)
          END DO
       END DO
    END IF
    CALL MPI_BARRIER(MPI_COMM_WORLD,CODE)
    
  END SUBROUTINE TRANSGLB_RHS
   !*********************************************************************************************!
  !**********************************************************************************************!
  !**********************************************************************************************!
  SUBROUTINE TRANSLOC_RHS(RHSL,RHSG,ABS_GLB,ORD_GLB,MYPROC,CODE,NUM_FIELD)

    USE MPI

    IMPLICIT NONE

    TYPE(RHS_LOC),                                    INTENT(INOUT) :: RHSL
    TYPE(RHS_GLB),                                    INTENT(IN)    :: RHSG
    INTEGER(8), DIMENSION(0:NB_PROCS-1),              INTENT(IN)    :: ABS_GLB,ORD_GLB
    INTEGER,                                          INTENT(IN)    :: MYPROC
    INTEGER,                                          INTENT(IN)    :: NUM_FIELD
    INTEGER,                                          INTENT(INOUT) :: CODE
    
    INTEGER(8)                                                      :: JL,NY,NL
    INTEGER                                                         :: NPROC
    INTEGER, DIMENSION(MPI_STATUS_SIZE)                             :: STATUS

    DO JL = 1, NLEV
       IF (MYPROC .EQ. 0) THEN
          DO NY = 1, LYPT
             NL = (NY-1)*TRANCHEX+HALO
             RHSL%T(JL,NL+1:NL+LXPT)   = RHSG%T(JL,1:LXPT,NY)
             RHSL%Q(JL-1,NL+1:NL+LXPT) = RHSG%Q(JL-1,1:LXPT,NY)
             DO NPROC = 1, NB_PROCS-1          
                CALL MPI_SEND(RHSG%T(JL,ABS_GLB(NPROC)+1:ABS_GLB(NPROC)+LXPT,ORD_GLB(NPROC)+NY),   &
                     & LXPT,MPI_REAL8,NPROC,NTAG(NUM_FIELD+1,NY,NPROC),MPI_COMM_WORLD,CODE)
                CALL MPI_SEND(RHSG%Q(JL-1,ABS_GLB(NPROC)+1:ABS_GLB(NPROC)+LXPT,ORD_GLB(NPROC)+NY), &
                     & LXPT,MPI_REAL8,NPROC,NTAG(NUM_FIELD+2,NY,NPROC),MPI_COMM_WORLD,CODE)
             END DO
          END DO
       ELSE
          DO NY = 1, LYPT
             NL = (NY-1)*TRANCHEX+HALO
             CALL MPI_RECV(RHSL%T(JL,NL+1:NL+LXPT),LXPT,MPI_REAL8,0, &
                  & NTAG(NUM_FIELD+1,NY,MYPROC),MPI_COMM_WORLD,STATUS,CODE)
             CALL MPI_RECV(RHSL%Q(JL-1,NL+1:NL+LXPT),LXPT,MPI_REAL8,0, &
                  & NTAG(NUM_FIELD+2,NY,MYPROC),MPI_COMM_WORLD,STATUS,CODE)
          END DO
       END IF
       CALL MPI_BARRIER(MPI_COMM_WORLD,CODE)
    END DO

    IF (MYPROC .EQ. 0) THEN
       DO NY = 1, LYPT
          NL = (NY-1)*TRANCHEX+HALO
          RHSL%Q(NLEV,NL+1:NL+LXPT) = RHSG%Q(NLEV,1:LXPT,NY)
          RHSL%PIS(NL+1:NL+LXPT)    = RHSG%PIS(1:LXPT,NY)
          DO NPROC = 1, NB_PROCS-1
             CALL MPI_SEND(RHSG%Q(NLEV,ABS_GLB(NPROC)+1:ABS_GLB(NPROC)+LXPT,ORD_GLB(NPROC)+NY), &
                  & LXPT,MPI_REAL8,NPROC,NTAG(NUM_FIELD+3,NY,NPROC),MPI_COMM_WORLD,CODE)
             CALL MPI_SEND(RHSG%PIS(ABS_GLB(NPROC)+1:ABS_GLB(NPROC)+LXPT,ORD_GLB(NPROC)+NY),    &
                  & LXPT,MPI_REAL8,NPROC,NTAG(NUM_FIELD+4,NY,NPROC),MPI_COMM_WORLD,CODE)
          END DO
      END DO   
    ELSE
      DO NY = 1, LYPT
         NL = (NY-1)*TRANCHEX+HALO
         CALL MPI_RECV(RHSL%Q(NLEV,NL+1:NL+LXPT),LXPT,MPI_REAL8,0, &
              & NTAG(NUM_FIELD+3,NY,MYPROC),MPI_COMM_WORLD,STATUS,CODE)
         CALL MPI_RECV(RHSL%PIS(NL+1:NL+LXPT),LXPT,MPI_REAL8,0,    &
              & NTAG(NUM_FIELD+4,NY,MYPROC),MPI_COMM_WORLD,STATUS,CODE)
      END DO
    END IF
   
    CALL MPI_BARRIER(MPI_COMM_WORLD,CODE)
   
  END SUBROUTINE TRANSLOC_RHS
  !****************************************************************************************!
  !****************************************************************************************!
  !****************************************************************************************!
  SUBROUTINE WIND_GLB_TRANSFER(GLB1,GLB2,GLB3,LOC1,LOC2,LOC3, &
             & TOPIN,ABS_GLB,ORD_GLB,MYPROC,CODE,NUM_FIELD)

    USE MPI

    IMPLICIT NONE
  
    REAL(8), DIMENSION(TOPIN:NLEV,GXPT,GYPT),         INTENT(OUT)   :: GLB1,GLB2,GLB3
    REAL(8), DIMENSION(TOPIN:NLEV,NDIM),              INTENT(IN)    :: LOC1,LOC2,LOC3
    INTEGER(8), DIMENSION(0:NB_PROCS-1),              INTENT(IN)    :: ABS_GLB,ORD_GLB
    INTEGER,                                          INTENT(IN)    :: TOPIN
    INTEGER,                                          INTENT(IN)    :: MYPROC
    INTEGER,                                          INTENT(IN)    :: NUM_FIELD
    INTEGER,                                          INTENT(INOUT) :: CODE
    
    INTEGER(8)                                                      :: JL,NY,NL
    INTEGER                                                         :: NPROC
    INTEGER, DIMENSION(MPI_STATUS_SIZE)                             :: STATUS

    DO JL = TOPIN, NLEV
       IF (MYPROC .NE. 0) THEN
          DO NY = 1, LYPT
             NL = (NY-1)*TRANCHEX+HALO
             CALL MPI_SEND(LOC1(JL,NL+1:NL+LXPT),LXPT,MPI_REAL8,0, &
                  & NTAG(NUM_FIELD+1,NY,MYPROC),MPI_COMM_WORLD,CODE)
             CALL MPI_SEND(LOC2(JL,NL+1:NL+LXPT),LXPT,MPI_REAL8,0, &
                  & NTAG(NUM_FIELD+2,NY,MYPROC),MPI_COMM_WORLD,CODE)
             CALL MPI_SEND(LOC3(JL,NL+1:NL+LXPT),LXPT,MPI_REAL8,0, &
                  & NTAG(NUM_FIELD+3,NY,MYPROC),MPI_COMM_WORLD,CODE)
          END DO
       ELSE
          DO NY = 1, LYPT
             NL = (NY-1)*TRANCHEX+HALO
             GLB1(JL,1:LXPT,NY) = LOC1(JL,NL+1:NL+LXPT)
             GLB2(JL,1:LXPT,NY) = LOC2(JL,NL+1:NL+LXPT)
             GLB3(JL,1:LXPT,NY) = LOC3(JL,NL+1:NL+LXPT)
             DO NPROC = 1, NB_PROCS-1
                CALL MPI_RECV(GLB1(JL,ABS_GLB(NPROC)+1:ABS_GLB(NPROC)+LXPT,ORD_GLB(NPROC)+NY),   &
                     & LXPT,MPI_REAL8,NPROC,NTAG(NUM_FIELD+1,NY,NPROC),MPI_COMM_WORLD,STATUS,CODE)
                CALL MPI_RECV(GLB2(JL,ABS_GLB(NPROC)+1:ABS_GLB(NPROC)+LXPT,ORD_GLB(NPROC)+NY),   &
                     & LXPT,MPI_REAL8,NPROC,NTAG(NUM_FIELD+2,NY,NPROC),MPI_COMM_WORLD,STATUS,CODE)
                CALL MPI_RECV(GLB3(JL,ABS_GLB(NPROC)+1:ABS_GLB(NPROC)+LXPT,ORD_GLB(NPROC)+NY),   &
                     & LXPT,MPI_REAL8,NPROC,NTAG(NUM_FIELD+3,NY,NPROC),MPI_COMM_WORLD,STATUS,CODE)
             END DO
          END DO
       END IF
       CALL MPI_BARRIER(MPI_COMM_WORLD,CODE)
    END DO

  END SUBROUTINE WIND_GLB_TRANSFER
  !***********************************************************************************************!
  !***********************************************************************************************!
  !***********************************************************************************************!
  !***********************************************************************************************!
  !***********************************************************************************************!
  !***********************************************************************************************!
  SUBROUTINE TRANSLOC_VAR(STL,STG,ABS_GLB,ORD_GLB,MYPROC,CODE)

    USE MPI

    IMPLICIT NONE

    TYPE(VARPROG),                                    INTENT(INOUT) :: STL
    TYPE(VARGLOB),                                    INTENT(INOUT) :: STG
    INTEGER(8), DIMENSION(0:NB_PROCS-1),              INTENT(IN)    :: ABS_GLB,ORD_GLB
    INTEGER,                                          INTENT(IN)    :: MYPROC
    INTEGER,                                          INTENT(INOUT) :: CODE
    
    INTEGER(8)                                                      :: JL,NY,NL
    INTEGER                                                         :: NPROC
    INTEGER, DIMENSION(MPI_STATUS_SIZE)                             :: STATUS
    
    DO JL=1,NLEV
       IF (MYPROC .EQ. 0) THEN
         DO NY = 1, LYPT
            NL = (NY-1)*TRANCHEX+HALO
            STL%U(JL,NL+1:NL+LXPT)   = STG%U(JL,1:LXPT,NY)
            STL%V(JL,NL+1:NL+LXPT)   = STG%V(JL,1:LXPT,NY)
            STL%T(JL,NL+1:NL+LXPT)   = STG%T(JL,1:LXPT,NY)
            STL%Q(JL-1,NL+1:NL+LXPT) = STG%Q(JL-1,1:LXPT,NY)
            STL%DIV(JL,NL+1:NL+LXPT) = STG%DIV(JL,1:LXPT,NY)
            DO NPROC = 1, NB_PROCS-1
               CALL MPI_SEND(STG%U(JL,ABS_GLB(NPROC)+1:ABS_GLB(NPROC)+LXPT,ORD_GLB(NPROC)+NY),   &
                    & LXPT,MPI_REAL8,NPROC,NTAG(1,NY,NPROC),MPI_COMM_WORLD,CODE)
               CALL MPI_SEND(STG%V(JL,ABS_GLB(NPROC)+1:ABS_GLB(NPROC)+LXPT,ORD_GLB(NPROC)+NY),   &
                    & LXPT,MPI_REAL8,NPROC,NTAG(2,NY,NPROC),MPI_COMM_WORLD,CODE)
               CALL MPI_SEND(STG%T(JL,ABS_GLB(NPROC)+1:ABS_GLB(NPROC)+LXPT,ORD_GLB(NPROC)+NY),   &
                    & LXPT,MPI_REAL8,NPROC,NTAG(3,NY,NPROC),MPI_COMM_WORLD,CODE)
               CALL MPI_SEND(STG%Q(JL-1,ABS_GLB(NPROC)+1:ABS_GLB(NPROC)+LXPT,ORD_GLB(NPROC)+NY), &
                    & LXPT,MPI_REAL8,NPROC,NTAG(4,NY,NPROC),MPI_COMM_WORLD,CODE)
               CALL MPI_SEND(STG%DIV(JL,ABS_GLB(NPROC)+1:ABS_GLB(NPROC)+LXPT,ORD_GLB(NPROC)+NY), &
                    & LXPT,MPI_REAL8,NPROC,NTAG(5,NY,NPROC),MPI_COMM_WORLD,CODE)
            END DO
         END DO   
       ELSE
         DO NY = 1, LYPT
            NL = (NY-1)*TRANCHEX+HALO
            CALL MPI_RECV(STL%U(JL,NL+1:NL+LXPT),LXPT,MPI_REAL8,0,NTAG(1,NY,MYPROC),   &
                 & MPI_COMM_WORLD,STATUS,CODE)
            CALL MPI_RECV(STL%V(JL,NL+1:NL+LXPT),LXPT,MPI_REAL8,0,NTAG(2,NY,MYPROC),   &
                 & MPI_COMM_WORLD,STATUS,CODE)
            CALL MPI_RECV(STL%T(JL,NL+1:NL+LXPT),LXPT,MPI_REAL8,0,NTAG(3,NY,MYPROC),   &
                 & MPI_COMM_WORLD,STATUS,CODE)
            CALL MPI_RECV(STL%Q(JL-1,NL+1:NL+LXPT),LXPT,MPI_REAL8,0,NTAG(4,NY,MYPROC), &
                 & MPI_COMM_WORLD,STATUS,CODE)
            CALL MPI_RECV(STL%DIV(JL,NL+1:NL+LXPT),LXPT,MPI_REAL8,0,NTAG(5,NY,MYPROC), &
                 & MPI_COMM_WORLD,STATUS,CODE)
         END DO
       END IF
       CALL MPI_BARRIER(MPI_COMM_WORLD,CODE)
    END DO

    IF (MYPROC .EQ. 0) THEN
       DO NY = 1, LYPT
          NL = (NY-1)*TRANCHEX+HALO
          STL%Q(NLEV,NL+1:NL+LXPT) = STG%Q(NLEV,1:LXPT,NY)
          STL%PIS(NL+1:NL+LXPT)    = STG%PIS(1:LXPT,NY)
          DO NPROC = 1, NB_PROCS-1
             CALL MPI_SEND(STG%Q(NLEV,ABS_GLB(NPROC)+1:ABS_GLB(NPROC)+LXPT,ORD_GLB(NPROC)+NY), &
                  & LXPT,MPI_REAL8,NPROC,NTAG(6,NY,NPROC),MPI_COMM_WORLD,CODE)
             CALL MPI_SEND(STG%PIS(ABS_GLB(NPROC)+1:ABS_GLB(NPROC)+LXPT,ORD_GLB(NPROC)+NY),    &
                  & LXPT,MPI_REAL8,NPROC,NTAG(7,NY,NPROC),MPI_COMM_WORLD,CODE)
          END DO
      END DO   
    ELSE
      DO NY = 1, LYPT
         NL = (NY-1)*TRANCHEX+HALO
         CALL MPI_RECV(STL%Q(NLEV,NL+1:NL+LXPT),LXPT,MPI_REAL8,0,NTAG(6,NY,MYPROC), &
              & MPI_COMM_WORLD,STATUS,CODE)
         CALL MPI_RECV(STL%PIS(NL+1:NL+LXPT),LXPT,MPI_REAL8,0,NTAG(7,NY,MYPROC),    &
              & MPI_COMM_WORLD,STATUS,CODE)
      END DO
   END IF
   CALL MPI_BARRIER(MPI_COMM_WORLD,CODE)
       

  END SUBROUTINE TRANSLOC_VAR
  !******************************************************************************************************!
  !******************************************************************************************************!
  !******************************************************************************************************!
  SUBROUTINE TRANSGLB_VAR(STG,STL,ABS_GLB,ORD_GLB,MYPROC,CODE)

    USE MPI

    IMPLICIT NONE

    TYPE(VARGLOB),                                    INTENT(INOUT) :: STG
    TYPE(VARPROG),                                    INTENT(IN)    :: STL
    INTEGER(8), DIMENSION(0:NB_PROCS-1),              INTENT(IN)    :: ABS_GLB,ORD_GLB
    INTEGER,                                          INTENT(IN)    :: MYPROC
    INTEGER,                                          INTENT(INOUT) :: CODE
    
    INTEGER(8)                                                      :: JL,NY,NL
    INTEGER                                                         :: NPROC
    INTEGER, DIMENSION(MPI_STATUS_SIZE)                             :: STATUS

     DO JL = 1, NLEV
       IF (MYPROC .NE. 0) THEN
          DO NY = 1, LYPT
             NL = (NY-1)*TRANCHEX+HALO
             CALL MPI_SEND(STL%U(JL,NL+1:NL+LXPT),LXPT,MPI_REAL8,0,   &
                  & NTAG(3,NY,MYPROC),MPI_COMM_WORLD,CODE)
             CALL MPI_SEND(STL%V(JL,NL+1:NL+LXPT),LXPT,MPI_REAL8,0,   &
                  & NTAG(4,NY,MYPROC),MPI_COMM_WORLD,CODE)
             CALL MPI_SEND(STL%T(JL,NL+1:NL+LXPT),LXPT,MPI_REAL8,0,   &
                  & NTAG(5,NY,MYPROC),MPI_COMM_WORLD,CODE)
             CALL MPI_SEND(STL%Q(JL-1,NL+1:NL+LXPT),LXPT,MPI_REAL8,0, &
                  & NTAG(6,NY,MYPROC),MPI_COMM_WORLD,CODE)
          END DO
       ELSE
          DO NY = 1, LYPT
             NL = (NY-1)*TRANCHEX+HALO
             STG%U(JL,1:LXPT,NY)   = STL%U(JL,NL+1:NL+LXPT)
             STG%V(JL,1:LXPT,NY)   = STL%V(JL,NL+1:NL+LXPT)
             STG%T(JL,1:LXPT,NY)   = STL%T(JL,NL+1:NL+LXPT)
             STG%Q(JL-1,1:LXPT,NY) = STL%Q(JL-1,NL+1:NL+LXPT)
             DO NPROC = 1, NB_PROCS-1
                CALL MPI_RECV(STG%U(JL,ABS_GLOBAL(NPROC)+1:ABS_GLOBAL(NPROC)+LXPT,ORD_GLOBAL(NPROC)+NY),   &
                     & LXPT,MPI_REAL8,NPROC,NTAG(3,NY,NPROC),MPI_COMM_WORLD,STATUS,CODE)
                CALL MPI_RECV(STG%V(JL,ABS_GLOBAL(NPROC)+1:ABS_GLOBAL(NPROC)+LXPT,ORD_GLOBAL(NPROC)+NY),   &
                     & LXPT,MPI_REAL8,NPROC,NTAG(4,NY,NPROC),MPI_COMM_WORLD,STATUS,CODE)
                CALL MPI_RECV(STG%T(JL,ABS_GLOBAL(NPROC)+1:ABS_GLOBAL(NPROC)+LXPT,ORD_GLOBAL(NPROC)+NY),   &
                     & LXPT,MPI_REAL8,NPROC,NTAG(5,NY,NPROC),MPI_COMM_WORLD,STATUS,CODE)
                CALL MPI_RECV(STG%Q(JL-1,ABS_GLOBAL(NPROC)+1:ABS_GLOBAL(NPROC)+LXPT,ORD_GLOBAL(NPROC)+NY), &
                     & LXPT,MPI_REAL8,NPROC,NTAG(6,NY,NPROC),MPI_COMM_WORLD,STATUS,CODE)
             END DO
          END DO
       END IF
       CALL MPI_BARRIER(MPI_COMM_WORLD,CODE)
    END DO  
    
    IF (MYPROC .NE. 0) THEN
       DO NY = 1, LYPT
          NL = (NY-1)*TRANCHEX+HALO
          CALL MPI_SEND(STL%Q(NLEV,NL+1:NL+LXPT),LXPT,MPI_REAL8,0, &
               & NTAG(8,NY,MYPROC),MPI_COMM_WORLD,CODE)
          CALL MPI_SEND(STL%PIS(NL+1:NL+LXPT),LXPT,MPI_REAL8,0, &
               & NTAG(9,NY,MYPROC),MPI_COMM_WORLD,CODE)
       END DO
    ELSE
       DO NY = 1, LYPT
          NL = (NY-1)*TRANCHEX+HALO
          STG%Q(NLEV,1:LXPT,NY) = STL%Q(NLEV,NL+1:NL+LXPT)
          STG%PIS(1:LXPT,NY)    = STL%PIS(NL+1:NL+LXPT)
          DO NPROC = 1, NB_PROCS-1
             CALL MPI_RECV(STG%Q(NLEV,ABS_GLOBAL(NPROC)+1:ABS_GLOBAL(NPROC)+LXPT,ORD_GLOBAL(NPROC)+NY),  &
                  & LXPT,MPI_REAL8,NPROC,NTAG(8,NY,NPROC),MPI_COMM_WORLD,STATUS,CODE)
             CALL MPI_RECV(STG%PIS(ABS_GLOBAL(NPROC)+1:ABS_GLOBAL(NPROC)+LXPT,ORD_GLOBAL(NPROC)+NY),     &
                  & LXPT,MPI_REAL8,NPROC,NTAG(9,NY,NPROC),MPI_COMM_WORLD,STATUS,CODE)
          END DO
       END DO
    END IF
    CALL MPI_BARRIER(MPI_COMM_WORLD,CODE)

  END SUBROUTINE TRANSGLB_VAR  
  !******************************************************************************************************!
  !******************************************************************************************************!
  !******************************************************************************************************!
  SUBROUTINE TRANSLOC_GEO(GEO,DOM,ABS_GLB,ORD_GLB,MYPROC,CODE)

    USE MOD_PARAMETER,  ONLY : NB_PROCS,TRANCHEX,TRANCHEY,HALO,LYPT,LXPT,NTAG
  
    IMPLICIT NONE 
  
    TYPE(GEOMETRY),                      INTENT(INOUT) :: GEO
    TYPE(GEODOM),                        INTENT(IN)    :: DOM
    INTEGER(8), DIMENSION(0:NB_PROCS-1), INTENT(IN)    :: ABS_GLB,ORD_GLB
    INTEGER,                             INTENT(IN)    :: MYPROC
    INTEGER,                             INTENT(INOUT) :: CODE

    INTEGER(8)                                         :: NY,NX,JL,NL
    INTEGER                                            :: NPROC
    INTEGER, DIMENSION(MPI_STATUS_SIZE)                :: STATUS  
    
    IF (MYPROC .EQ. 0) THEN
       DO NY = 1, LYPT
          NL = (NY-1)*TRANCHEX+HALO
          GEO%ZS(NL+1:NL+LXPT)           = DOM%ZS(1:LXPT,NY)
          GEO%DZSDX(NL+1:NL+LXPT)        = DOM%DZSDX(1:LXPT,NY)
          GEO%DZSDY(NL+1:NL+LXPT)        = DOM%DZSDY(1:LXPT,NY)
          !GEO%PIS_REF(NL+1:NL+LXPT)      = DOM%PIS_REF(1:LXPT,NY)
          !GEO%DPISDX_REF(NL+1:NL+LXPT)   = DOM%DPISDX_REF(1:LXPT,NY)
          !GEO%DPISDY_REF(NL+1:NL+LXPT)   = DOM%DPISDY_REF(1:LXPT,NY)
          DO NPROC = 1, NB_PROCS-1
             CALL MPI_SEND(DOM%ZS(ABS_GLB(NPROC)+1:ABS_GLB(NPROC)+LXPT,ORD_GLB(NPROC)+NY),      &
                  & LXPT,MPI_REAL8,NPROC,NTAG(3,NY,NPROC),MPI_COMM_WORLD,CODE)
             CALL MPI_SEND(DOM%DZSDX(ABS_GLB(NPROC)+1:ABS_GLB(NPROC)+LXPT,ORD_GLB(NPROC)+NY),   &
                  & LXPT,MPI_REAL8,NPROC,NTAG(4,NY,NPROC),MPI_COMM_WORLD,CODE)
             CALL MPI_SEND(DOM%DZSDY(ABS_GLB(NPROC)+1:ABS_GLB(NPROC)+LXPT,ORD_GLB(NPROC)+NY),   &
                  & LXPT,MPI_REAL8,NPROC,NTAG(5,NY,NPROC),MPI_COMM_WORLD,CODE)
             !CALL MPI_SEND(DOM%PIS_REF(ABS_GLB(NPROC)+1:ABS_GLB(NPROC)+LXPT,ORD_GLB(NPROC)+NY), &
             !     & LXPT,MPI_REAL8,NPROC,NTAG(6,NY,NPROC),MPI_COMM_WORLD,CODE)
             !CALL MPI_SEND(DOM%DPISDX_REF(ABS_GLB(NPROC)+1:ABS_GLB(NPROC)+LXPT,ORD_GLB(NPROC)+NY), &
             !     & LXPT,MPI_REAL8,NPROC,NTAG(7,NY,NPROC),MPI_COMM_WORLD,CODE)
             !CALL MPI_SEND(DOM%DPISDY_REF(ABS_GLB(NPROC)+1:ABS_GLB(NPROC)+LXPT,ORD_GLB(NPROC)+NY), &
             !     & LXPT,MPI_REAL8,NPROC,NTAG(8,NY,NPROC),MPI_COMM_WORLD,CODE)
          END DO
       END DO   
    ELSE
       DO NY = 1, LYPT
          NL = (NY-1)*TRANCHEX+HALO
          CALL MPI_RECV(CPL%RELAX(NL+1:NL+LXPT),LXPT,MPI_REAL8,0,NTAG(2,NY,MYPROC),   &
               & MPI_COMM_WORLD,STATUS,CODE)
          CALL MPI_RECV(GEO%ZS(NL+1:NL+LXPT),LXPT,MPI_REAL8,0,NTAG(3,NY,MYPROC),      &
               & MPI_COMM_WORLD,STATUS,CODE)
          CALL MPI_RECV(GEO%DZSDX(NL+1:NL+LXPT),LXPT,MPI_REAL8,0,NTAG(4,NY,MYPROC),   &
               & MPI_COMM_WORLD,STATUS,CODE)
          CALL MPI_RECV(GEO%DZSDY(NL+1:NL+LXPT),LXPT,MPI_REAL8,0,NTAG(5,NY,MYPROC),   &
               & MPI_COMM_WORLD,STATUS,CODE)
          !CALL MPI_RECV(GEO%PIS_REF(NL+1:NL+LXPT),LXPT,MPI_REAL8,0,NTAG(6,NY,MYPROC), &
          !     & MPI_COMM_WORLD,STATUS,CODE)
          !CALL MPI_RECV(GEO%DPISDX_REF(NL+1:NL+LXPT),LXPT,MPI_REAL8,0,NTAG(7,NY,MYPROC), &
          !     & MPI_COMM_WORLD,STATUS,CODE)
          !CALL MPI_RECV(GEO%DPISDY_REF(NL+1:NL+LXPT),LXPT,MPI_REAL8,0,NTAG(8,NY,MYPROC), &
          !     & MPI_COMM_WORLD,STATUS,CODE)
       END DO
    END IF
    CALL MPI_BARRIER(MPI_COMM_WORLD,CODE)

    DO NY = 1,LYPT
       DO NX = 1,LXPT
          NL = (HALO+NY-1)*TRANCHEX + HALO + NX
          GEO%PIS_REF(NL)    =  RP00*DEXP(-(RG/(RR*RT00))*GEO%ZS(NL))
          GEO%DPISDX_REF(NL) = -(RG/(RR*RT00))*GEO%PIS_REF(NL)*GEO%DZSDX(NL)
          GEO%DPISDY_REF(NL) = -(RG/(RR*RT00))*GEO%PIS_REF(NL)*GEO%DZSDY(NL)
       END DO   
    END DO
    
  END SUBROUTINE TRANSLOC_GEO
  !*********************************************************************************************!
  !*********************************************************************************************!
  !*********************************************************************************************!
  !*********************************************************************************************!
  !*********************************************************************************************!
  !*********************************************************************************************!
  !*********************************************************************************************!
  !*********************************************************************************************!
  !*********************************************************************************************!
  !*********************************************************************************************!
  !*********************************************************************************************!
  SUBROUTINE COMM_INTER_(PX,KDIM,KLEV,KPROC,KEIGHOUR,KERROR)

    USE MPI
    USE MOD_PARAMETER, ONLY : NLEV,LXPT,LYPT,NDIM,HALO

    IMPLICIT NONE

    REAL(8), DIMENSION(KLEV:NLEV,NDIM),            INTENT(OUT)   :: PX
    INTEGER,                                       INTENT(IN)    :: KDIM,KLEV,KPROC
    INTEGER, DIMENSION(4),                         INTENT(IN)    :: KEIGHOUR
    INTEGER,                                       INTENT(INOUT) :: KERROR
    INTEGER                                                      :: IDIM,ISUP,IINF
    INTEGER, DIMENSION(2)                                        :: IREQUEST
    INTEGER, DIMENSION(MPI_STATUS_SIZE)                          :: ISTATUS
    REAL(8), DIMENSION(2,2,HALO*(NLEV+1-KLEV)*KDIM)              :: ZBOX

    IDIM = HALO*(NLEV+1-KLEV)*KDIM
    IF (KDIM .EQ. LXPT) THEN
       ISUP = 3
       IINF = 2
    ELSE
       ISUP = 4
       IINF = 1
    END IF
    
    CALL PACKING_(ZBOX(1,:,:),PVAR,KDIM,KLEV)

    CALL MPI_IRECV(ZBOX(2,2,1:IDIM),IDIM,MPI_REAL8,KEIGHOUR(ISUP), &
         & KEIGHOUR(ISUP),MPI_COMM_WORLD,IREQUEST(1),KERROR)
    CALL MPI_SEND(ZBOX(1,1,1:IDIM),IDIM,MPI_REAL8,KEIGHOUR(IINF),  &
         & KPROC,MPI_COMM_WORLD,KERROR)
    CALL MPI_WAIT(IREQUEST(1),ISTATUS,CODEINOUT)

    CALL MPI_IRECV(ZBOX(2,1,1:IDIM),IDIM,MPI_REAL8,KEIGHOUR(IINF), &
         & KEIGHOUR(IINF),MPI_COMM_WORLD,IREQUEST(2),KERROR)
    CALL MPI_SEND(ZBOX(1,2,1:IDIM),IDIM,MPI_REAL8,KEIGHOUR(ISUP),  &
         & KPROC,MPI_COMM_WORLD,KERROR)
    CALL MPI_WAIT(IREQUEST(2),ISTATUS,KERROR)

    CALL UNPACKING_(PVAR,BOX(2,:,:),KDIM,KLEV)

    CALL MPI_BARRIER(MPI_COMM_WORLD,KERROR)

  END SUBROUTINE COMM_INTER_
  !************************************************************************************************!
  !************************************************************************************************!
  !************************************************************************************************!
  !************************************************************************************************!
  SUBROUTINE PACKING_(PBOX,PX,KDIM,KLEV)
    ! Construction des boite de messages pour l'envoi
    ! le premier indice montre que la boite lit
    ! le deuxième montre quelle partie du domaine il lit
    ! le troisième est l'indice de lecture
    
    USE MOD_PARAMETER, ONLY : NLEV,LXPT,LYPT,HALO,TRANCHEX,TRANCHEY

    IMPLICIT NONE
    
    REAL(8), DIMENSION(2,HALO*(NLEV+1-KLEV)*KDIM),  INTENT(OUT) :: PBOX
    REAL(8), DIMENSION(KLEV:NLEV,NDIM),             INTENT(IN)  :: PX
    INTEGER,                                        INTENT(IN)  :: KDIM,KLEV
    INTEGER                                                     :: IX,IY,IL,INB,INV

    IF (KDIM .EQ. LYPT) THEN

       ! Remplissage de la boite d'envoie
       DO IL = 1, NLEV
          DO IY = 1, LYPT
             DO IX = 1, HALO
                ! indice des variables pronostiques à l'Ouest
                INV = HALO + IX + (IY-1)*TRANCHEX
                INB = IX + (IY-1)*HALO + (IL-1)*HALO*LYPT
                PBOX(1,INB) = PX(IL,INV)

                ! indice des variables pronostiques à l'Est
                INV = LXPT + IX + (IY-1)*TRANCHEX
                INB = IX + (IY-1)*HALO + (IL-1)*HALO*LYPT
                PBOX(2,INB) = PVAR%UX(IL,INV)
             END DO
          END DO
       END DO

    ELSE

       ! Remplissage de la boite d'envoie
       DO IL = 1, NLEV
          DO IY = 1, NHALO
             DO IX = 1, LXPT
                ! indice des variables pronostiques au Nord
                INV = HALO + IX + (LYPT+IY-1)*TRANCHEX
                INB = IX + (IY-1)*HALO + (IL-1)*HALO*LXPT
                PBOX(1,INB) = PX(IL,INV)

                ! indice des variables pronostiques au Sud
                INV = HALO + IX + (HALO+IY-1)*TRANCHEX
                INB = IX + (IY-1)*HALO + (IL-1)*HALO*LXPT
                PBOX(2,INB) = PX(IL,INV)

             END DO
          END DO
       END DO
       
    END IF

  END SUBROUTINE PACKING_
  !********************************************************************************************!
  !********************************************************************************************!
  !********************************************************************************************!
  !********************************************************************************************!
  SUBROUTINE UNPACKING_(PBOX,PX,KDIM,KLEV)
    ! Construction des boite de messages pour l'envoi
    ! le premier indice montre que la boite est lu
    ! le deuxième montre quelle partie du domaine est lue
    ! le troisième est l'indice de lecture
    
    USE MOD_PARAMETER, ONLY : NLEV,LXPT,LYPT,HALO,TRANCHEX,TRANCHEY

    IMPLICIT NONE
    
    REAL(8), DIMENSION(KLEV:NLEV,NDIM),             INTENT(INOUT) :: PX
    REAL(8), DIMENSION(2,HALO*(NLEV+1-KLEV)*KDIM),  INTENT(IN)    :: PBOX
    INTEGER,                                        INTENT(IN)    :: KDIM,KLEV
    INTEGER                                                       :: IX,IY,IL,INB,INV

    IF (KDIM .EQ. LYPT) THEN

       ! Lecture de la boite
       DO IL = 1, NLEV
          DO IY = 1, LYPT
             DO IX = 1, HALO
                ! indice de la boite
                INB = IX + (IY-1)*HALO + (IL-1)*HALO*LYPT
                ! indice de la variable pronostique à l'Ouest
                INV = IX + (IY-1)*TRANCHEX
                PX(IL,INV) = PBOX(1,INB)

                ! indice de la variable pronostique à l'Est
                INV = HALO + LXPT + IX + (IY-1)*TRANCHEX
                PX(IL,INV) = PBOX(2,INB)
             END DO
          END DO
       END DO

    ELSE

       ! Lecture de la boite
       DO IL = 1, NLEV
          DO IY = 1, HALO
             DO IX = 1, LXPT
                ! indice de la boite
                INB = IX + (IY-1)*HALO + (IL-1)*HALO*LXPT
                ! indice de la variable pronostique au Nord
                INV = HALO + IX + (HALO+LYPT+IY-1)*TRANCHEX
                PX(IL,INV) = PBOX(1,INB)

                ! indice de la variable pronostique au Sud
                INV = HALO + IX + (IY-1)*TRANCHEX
                PX(IL,INV) = BOX(2,INB)
             END DO
          END DO
       END DO

    END IF

  END SUBROUTINE UNPACKING_
  !*********************************************************************************************!
  !*********************************************************************************************!
  !*********************************************************************************************!
  !*********************************************************************************************!
  !*********************************************************************************************!

   
END MODULE MOD_COMM

!=====================================================!
!=====================================================!
!=====================================================!
