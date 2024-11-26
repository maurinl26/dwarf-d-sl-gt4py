!##################   PARAMETRES   ###################!
!#                                                   #!
!# auteur : C.Colavolpe & F.Voitus                   #!
!# sujet  : Définit les parametres du modèle         #!
!#                                                   #!
!#####################################################!

!=====================================================!
!================  MODULE PARAMETRES  ================!
!=====================================================!

MODULE MOD_COMM

  USE MOD_PARAMETER, ONLY : GXPT,GYPT,NSMAX,MSMAX,RPI,RDX,RDY,RLX,RLY,&
                          & ZERO,ONE,TWO,HALO,LXPT,LYPT,NDIM

CONTAINS

  !*****************************************************!
  !***************  LISTE DES ROUTINES  ****************!
  !*****************************************************!


  !****************************************************************************************
  !****************************************************************************************
  !****************************************************************************************
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
  !****************************************************************************************
  !****************************************************************************************
  !****************************************************************************************
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
  !****************************************************************************************
  !****************************************************************************************
  !****************************************************************************************
  SUBROUTINE TRANSLOC_MPI(PSI_L,DPSIDX_L,DPSIDY_L,PSI_G,DPSIDX_G,DPSIDY_G, &
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
             CALL MPI_SEND(PSI_G(ABS_GLB(NPROC)+1:ABS_GLB(NPROC)+LXPT,ORD_GLB(NPROC)+NY),        &
                  & LXPT,MPI_REAL8,NPROC,NTAG(NUM,NY,NPROC),MPI_COMM_WORLD,CODE)
             CALL MPI_SEND(DPSIDX_G(ABS_GLB(NPROC)+1:ABS_GLB(NPROC)+LXPT,ORD_GLB(NPROC)+NY),     &
                  & LXPT,MPI_REAL8,NPROC,NTAG(NUM+2,NY,NPROC),MPI_COMM_WORLD,CODE)
             CALL MPI_SEND(DPSIDY_G(ABS_GLB(NPROC)+1:ABS_GLB(NPROC)+LXPT,ORD_GLB(NPROC)+NY),     &
                  & LXPT,MPI_REAL8,NPROC,NTAG(NUM+4,NY,NPROC),MPI_COMM_WORLD,CODE)
          END DO
       END DO   
    ELSE
       DO NY = 1, LYPT
          NL = (NY-1)*TRANCHEX+HALO
          CALL MPI_RECV(PSI_L(NL+1:NL+LXPT),LXPT,MPI_REAL8,0,NTAG(NUM,NY,MYPROC),        &
               & MPI_COMM_WORLD,STATUS,CODE)
          CALL MPI_RECV(DPSIDX_L(NL+1:NL+LXPT),LXPT,MPI_REAL8,0,NTAG(NUM+2,NY,MYPROC),   &
               & MPI_COMM_WORLD,STATUS,CODE)
          CALL MPI_RECV(DPSIDY_L(NL+1:NL+LXPT),LXPT,MPI_REAL8,0,NTAG(NUM+4,NY,MYPROC),   &
               & MPI_COMM_WORLD,STATUS,CODE)
       END DO
    END IF
    CALL MPI_BARRIER(MPI_COMM_WORLD,CODE)

  END SUBROUTINE TRANSLOC_MPI
  !****************************************************************************************
  !****************************************************************************************
  !****************************************************************************************
  
  !********  Comunications inter-processorale  *********!
   SUBROUTINE COMM_INTER_(PX,KDIM,KLEV,KPROC,KEIGHOUR,KERROR)

    USE MPI
    USE MOD_PARAMETER, ONLY : NLEV,LXPT,LYPT,NSLICEX,NSLICEY,RWEIGHT

    IMPLICIT NONE

    REAL(8), DIMENSION(KLEV:NLEV,NSLICEX*NSLICEY), INTENT(OUT)   :: PX
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
    CALL MPI_SEND( ZBOX(1,1,1:IDIM),IDIM,MPI_REAL8,KEIGHOUR(IINF), &
            & KPROC,MPI_COMM_WORLD,KERROR)
    CALL MPI_WAIT(IREQUEST(1),ISTATUS,CODEINOUT)

    CALL MPI_IRECV(ZBOX(2,1,1:IDIM),IDIM,MPI_REAL8,KEIGHOUR(IINF), &
            & KEIGHOUR(IINF),MPI_COMM_WORLD,IREQUEST(2),KERROR)
    CALL MPI_SEND( ZBOX(1,2,1:IDIM),IDIM,MPI_REAL8,KEIGHOUR(ISUP), &
            & KPROC,MPI_COMM_WORLD,KERROR)
    CALL MPI_WAIT(IREQUEST(2),ISTATUS,KERROR)

    CALL UNPACKING_(PVAR,BOX(2,:,:),KDIM,KLEV)

    CALL MPI_BARRIER(MPI_COMM_WORLD,KERROR)

  END SUBROUTINE COMM_INTER_
  !*************************************************************************
  !*************************************************************************
  !*************************************************************************
  SUBROUTINE PACKING_(PBOX,PX,KDIM,KLEV)
    ! Construction des boite de messages pour l'envoi
    ! le premier indice montre que la boite lit
    ! le deuxième montre quelle partie du domaine il lit
    ! le troisième est l'indice de lecture
    
    USE MOD_PARAMETER, ONLY : NLEV,LXPT,LYPT,NHALO,NSLICEX,NSLICEY

    IMPLICIT NONE
    
    REAL(8), DIMENSION(2,NHALO*(NLEV+1-KLEV)*KDIM), INTENT(OUT) :: PBOX
    REAL(8), DIMENSION(KLEV:NLEV,NSLICEX*NSLICEY),  INTENT(IN)  :: PX
    INTEGER,                                        INTENT(IN)  :: KDIM,KLEV
    INTEGER                                                     :: IX,IY,IL,INB,INV

    IF (KDIM .EQ. LYPT) THEN

       ! Remplissage de la boite d'envoie
       DO IL = 1, NLEV
          DO IY = 1, LYPT
             DO IX = 1, NHALO
                ! indice des variables pronostiques à l'Ouest
                INV = NHALO + IX + (IY-1)*NSLICEX
                INB = IX + (IY-1)*NHALO + (IL-1)*NHALO*LYPT
                PBOX(1,INB) = PX(IL,INV)

                ! indice des variables pronostiques à l'Est
                INV = LXPT + IX + (IY-1)*NSLICEX
                INB = IX + (IY-1)*NHALO + (IL-1)*NHALO*LYPT
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
                INV = NHALO + IX + (LYPT+IY-1)*NSLICEX
                INB = IX + (IY-1)*NHALO + (IL-1)*NHALO*LXPT
                PBOX(1,INB) = PX(IL,INV)

                ! indice des variables pronostiques au Sud
                INV = HALO + IX + (NHALO+IY-1)*NSLICEX
                INB = IX + (IY-1)*NHALO + (IL-1)*NHALO*LXPT
                PBOX(2,INB) = PX(IL,INV)

             END DO
          END DO
       END DO
       
    END IF

  END SUBROUTINE PACKING_
  !**********************************************************************
  !**********************************************************************
  !**********************************************************************
  SUBROUTINE UNPACKING_(PBOX,PX,KDIM,KLEV)
    ! Construction des boite de messages pour l'envoi
    ! le premier indice montre que la boite est lu
    ! le deuxième montre quelle partie du domaine est lue
    ! le troisième est l'indice de lecture
    
    USE MOD_PARAMETER, ONLY : NLEV,LXPT,LYPT,NHALO,NSLICEX,NSLICEY

    IMPLICIT NONE
    
    REAL(8), DIMENSION(KLEV:NLEV,NSLICEX*NSLICEY),  INTENT(INOUT) :: PX
    REAL(8), DIMENSION(2,NHALO*(NLEV+1-KLEV)*KDIM), INTENT(IN)    :: PBOX
    INTEGER,                                        INTENT(IN)    :: KDIM,KLEV
    INTEGER                                                       :: IX,IY,IL,INB,INV

    IF (KDIM .EQ. LYPT) THEN

       ! Lecture de la boite
       DO IL = 1, NLEV
          DO IY = 1, LYPT
             DO IX = 1, NHALO
                ! indice de la boite
                INB = IX + (IY-1)*NHALO + (IL-1)*NHALO*LYPT
                ! indice de la variable pronostique à l'Ouest
                INV = IX + (IY-1)*NSLICEX
                PX(IL,INV) = PBOX(1,INB)

                ! indice de la variable pronostique à l'Est
                INV = HALO + LXPT + IX + (IY-1)*NSLICEX
                PX(IL,INV) = PBOX(2,INB)
             END DO
          END DO
       END DO

    ELSE

       ! Lecture de la boite
       DO IL = 1, NLEV
          DO IY = 1, NHALO
             DO IX = 1, LXPT
                ! indice de la boite
                INB = IX + (IY-1)*NHALO + (IL-1)*NHALO*LXPT
                ! indice de la variable pronostique au Nord
                INV = NHALO + IX + (HALO+LYPT+IY-1)*NSLICEX
                PX(IL,INV) = PBOX(1,INB)

                ! indice de la variable pronostique au Sud
                INV = HALO + IX + (IY-1)*NSLICEX
                PX(IL,INV) = BOX(2,INB)
             END DO
          END DO
       END DO

    END IF

  END SUBROUTINE UNPACKING_
  !*******************************************************************************************!
  !*******************************************************************************************!
  !*******************************************************************************************!
  !*******************************************************************************************!
  SUBROUTINE GRAD_DX(PDXDX,PX,KLEV,KPROC,KEIGHOUR,KERROR)

    USE MOD_PARAMETER, ONLY : NLEV,LXPT,LYPT,NSLICEX,NSLICEY,RWEIGHT

    IMPLICIT NONE

    REAL(8), DIMENSION(KLEV:NLEV,NSLICEX*NSLICEY), INTENT(OUT)   :: PDXDX
    REAL(8), DIMENSION(KLEV:NLEV,NSLICEX*NSLICEY), INTENT(INOUT) :: PX
    INTEGER,                                       INTENT(IN)    :: KLEV,KPROC
    INTEGER, DIMENSION(4),                         INTENT(IN)    :: KEIGHOUR
    INTEGER,                                       INTENT(INOUT) :: KERROR
    INTEGER                                                      :: IX,IY,IN,IK

    CALL COMM_INTER_(PX,LXPT,KLEV,KPROC,KEIGHOUR,KERROR)

    DO IY = 1, LYPT
       DO IX = 1, LXPT
          IN = (NHALO+IY-1)*NSLICEX + NHALO + IX
          
          PDXDX(:,IN) = 0.d0
          DO IK = 1, NHALO
             PDXDX(:,IN) = PDXDX(:,IN) + RWEIGHT(IK)*(PX(:,IN+IK) - PX(:,IN-IK))
          END DO
          
       END DO
    END DO

  END SUBROUTINE GRAD_DX
  !***************************************************************************************!
  !***************************************************************************************!
  SUBROUTINE GRAD_DY(PDXDY,PX,KLEV,KEIGHOUR,KPROC,KERROR)

    USE MOD_PARAMETER, ONLY : NLEV,LXPT,LYPT,NSLICEX,NSLICEY,RWEIGHT

    IMPLICIT NONE

    REAL(8), DIMENSION(KLEV:NLEV,NSLICEX*NSLICEY), INTENT(OUT)   :: PDXDY
    REAL(8), DIMENSION(KLEV:NLEV,NSLICEX*NSLICEY), INTENT(INOUT) :: PX
    INTEGER,                                       INTENT(IN)    :: KLEV,KPROC
    INTEGER, DIMENSION(4),                         INTENT(IN)    :: KEIGHOUR
    INTEGER,                                       INTENT(INOUT) :: KERROR
    INTEGER                                                      :: IX,IY,IN,IK

    CALL COMM_INTER_(PX,LYPT,KLEV,KPROC,KEIGHOUR,KERROR)

    DO IX = 1, LXPT
       DO IY = 1, LYPT
          IN = (NHALO+IY-1)*NSLICEX + NHALO + IX
          
          PDXDY(:,IN) = 0.d0
          DO IK = 1, NHALO
             PDXDY(:,IN) = PDXDY(:,IN) + RWEIGHT(IK)*(PX(:,IN+IK*NSLICEX) - PX(:,IN-IK*NSLICEX))
          END DO
          
       END DO
    END DO

  END SUBROUTINE GRAD_DY
  !*********************************************************************************************!
  !*********************************************************************************************!
  !*********************************************************************************************!
  !*********************************************************************************************!

END MODULE MOD_PARAMETER

!=====================================================!
!=====================================================!
!=====================================================!
