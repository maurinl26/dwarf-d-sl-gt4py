!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! -*- Mode: F90 -*- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! GFT_jmfft.f90 --- J.-M. FFT routines
!! 
!! Auteur          : Jean-Marie Teuler, CNRS-IDRIS
!! Créé le         : Tue Feb 19 10:22:47 2002
!! Dern. mod. par  : Jalel Chergui (CNRS/IDRIS) <Jalel.Chergui@idris.fr>
!! Dern. mod. le   : Thu Jun 13 12:01:48 2002
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Permission is granted to copy and distribute this file or modified
! versions of this file for no fee, provided the copyright notice and
! this permission notice are preserved on all copies.
! Copyright (C) Februry 2002, CNRS/IDRIS, Jean-Marie Teuler.
!
MODULE GFT_JMFFT
USE GFT_common

CONTAINS

SUBROUTINE jmtable(table,ntable,itable,n)

  IMPLICIT NONE

  ! Arguments
  INTEGER, INTENT(in) :: ntable, itable, n
  REAL(kind=GFT_prec), INTENT(out), DIMENSION(0:ntable-1) :: table

  ! Variables locales
  REAL(kind=GFT_prec), SAVE :: twopi
  REAL(kind=GFT_prec) :: temp, temp1, temp2

  INTEGER :: i

  twopi = 2.0_GFT_prec * ACOS(-1.0_GFT_prec)

  ! Calcul des sinus et cosinus

  ! Si n est multiple de 4, astuces en serie
  IF (MOD(n,4) == 0) THEN
    ! On se debarrasse des cas limite
    table(itable+      0) =  1
    table(itable+n+    0) =  0
    table(itable+    n/4) =  0
    table(itable+n+  n/4) =  1
    table(itable+    n/2) = -1
    table(itable+n+  n/2) =  0
    table(itable+  3*n/4) =  0
    table(itable+n+3*n/4) = -1
    ! Cas general
!dir$ ivdep
!ocl novrec
!cdir nodep
    DO i = 1,n/4-1
      temp = COS(twopi*REAL(i,kind=GFT_prec)/REAL(n,kind=GFT_prec))
      table(itable+    i)     =  temp
      table(itable+    n/2-i) = -temp
      table(itable+    n/2+i) = -temp
      table(itable+    n-i)   =  temp
      table(itable+n+  n/4+i) =  temp
      table(itable+n+  n/4-i) =  temp
      table(itable+n+3*n/4+i) = -temp
      table(itable+n+3*n/4-i) = -temp
    END DO

  ! Si n est simplement multiple de 2 (sans etre multiple de 4)
  ELSE IF (MOD(n,2) == 0) THEN
    ! On se debarrasse des cas limite
    table(itable+    0) =  1
    table(itable+n+  0) =  0
    table(itable+  n/2) = -1
    table(itable+n+n/2) =  0
    ! Cas general
!dir$ ivdep
!ocl novrec
!cdir nodep
    DO i = 1,n/2-1
      temp1 = COS(twopi*REAL(i,kind=GFT_prec)/REAL(n,kind=GFT_prec))
      table(itable+      i) =  temp1
      table(itable+  n/2+i) = -temp1
      temp2 = SIN(twopi*REAL(i,kind=GFT_prec)/REAL(n,kind=GFT_prec))
      table(itable+n+    i) =  temp2
      table(itable+n+n/2+i) = -temp2
    END DO

  ! Si n est impair
  ELSE
    ! On se debarrasse des cas limite
    table(itable+  0) =  1
    table(itable+n+0) =  0
!dir$ ivdep
!ocl novrec
!cdir nodep
    DO i = 1,n/2
      temp1 = COS(twopi*REAL(i,kind=GFT_prec)/REAL(n,kind=GFT_prec))
      table(itable+    i) =  temp1
      table(itable+  n-i) =  temp1
      temp2 = SIN(twopi*REAL(i,kind=GFT_prec)/REAL(n,kind=GFT_prec))
      table(itable+n+  i) =  temp2
      table(itable+n+n-i) = -temp2
    END DO

  END IF

END SUBROUTINE jmtable

SUBROUTINE jmfact(n,fact,nfact,ideb,ifin)

  IMPLICIT NONE

  ! Arguments
  INTEGER, INTENT(in) :: n, nfact, ideb
  INTEGER, INTENT(out) :: ifin
  INTEGER, INTENT(inout), DIMENSION(0:nfact-1) :: fact

  ! Variables locales
  INTEGER :: m
  INTEGER :: n2, p2, n3, p3, n5, p5
  CHARACTER(len=*), PARAMETER :: nomsp = 'JMFACT'
  ! Nombres premiers
  INTEGER, PARAMETER :: npremiers = 7
  INTEGER, DIMENSION(0:npremiers-1) :: premiers = (/7,11,13,17,19,23,29/)
  INTEGER :: ip, premier, pp, np

  m = n

  ! Etude des puissances de deux
  p2 = 0
  n2 = 1
  DO
    IF (MOD(m,2) == 0) THEN
      p2 = p2+1
      n2 = n2*2
      m  = m/2
    ELSE
      EXIT
    END IF
  END DO
  ifin = ideb+3
  IF (ifin > nfact) &
  & CALL jmerreur2(nomsp,7,nfact,ifin)
  fact(ideb+1) = n2
  fact(ideb+2) = p2

  ! Etude des puissances de trois
  p3 = 0
  n3 = 1
  DO
    IF (MOD(m,3) == 0) THEN
      p3 = p3+1
      n3 = n3*3
      m  = m/3
    ELSE
      EXIT
    END IF
  END DO
  ifin = ifin+2
  IF (ifin > nfact) &
  & CALL jmerreur2(nomsp,7,nfact,ifin)
  fact(ideb+3) = n3
  fact(ideb+4) = p3

  ! Etude des puissances de cinq
  p5 = 0
  n5 = 1
  DO
    IF (MOD(m,5) == 0) THEN
      p5 = p5+1
      n5 = n5*5
      m  = m/5
    ELSE
      EXIT
    END IF
  END DO
  ifin = ifin+2
  IF (ifin > nfact) &
  & CALL jmerreur2(nomsp,7,nfact,ifin)
  fact(ideb+5) = n5
  fact(ideb+6) = p5

  ! On met a jour le nombre de termes
  fact(ideb) = 7

  ! Si on a fini
  IF (n2*n3*n5 == n) RETURN

  ! Il reste maintenant des facteurs premiers bizarres
  ! On va boucler tant qu'on n'a pas fini ou tant qu'on n'a pas epuise la liste

  DO ip = 0,npremiers-1

    premier = premiers(ip)

    pp = 0
    np = 1
    DO
      IF (MOD(m,premier) == 0) THEN
        pp = pp+1
        np = np*premier
        m  = m/premier
      ELSE
        EXIT
      END IF
    END DO
    ifin = ifin+2
    IF (ifin > nfact) &
    & CALL jmerreur2(nomsp,7,nfact,ifin)
    fact(ifin-2) = pp
    fact(ifin-1) = premier
    fact(ideb) = fact(ideb) + 2

    ! Si le nombre est completement factorise, inutile de continuer
    IF (m == 1) EXIT

  END DO

  ! On regarde si la factorisation est terminee
  IF (m == 1) THEN
    RETURN
  ELSE
    CALL jmerreur1(nomsp,3,n)
  END IF
  
END SUBROUTINE jmfact

! Pour tronconner de facon a tenir dans le nwork disponible
SUBROUTINE jmdecoup(n,nr,nwork,debut,mpair,n_temp,ideb,ifin,nwork_temp,fin)

  IMPLICIT NONE

  ! Arguments
  INTEGER, INTENT(in) :: n, nr, nwork
  LOGICAL, INTENT(in) :: debut, mpair
  INTEGER, INTENT(out) :: n_temp, ideb, nwork_temp
  INTEGER, INTENT(inout) :: ifin
  LOGICAL, INTENT(out) :: fin

  ! Variables locales
  CHARACTER(len=*), PARAMETER :: nomsp = 'JMDECOUP'

  ! n*nr est l'espace total qu'il faudrait pour work.
  ! Malheureusement, on n'a que nwork au plus
  ! On va donc decouper n en morceaux pour tenir

  ! Gestion de debut
  IF (debut) THEN
    ideb = 0
  ELSE
    ideb = ifin+1
  END IF

  ! Gestion de n_temp et ifin
  n_temp = nwork/nr
  ! Si m impair, on doit eviter que n_temp soit impair (routine cs et sc)
  IF (.NOT.mpair .AND. MOD(n_temp,2) /= 0) n_temp = n_temp-1
  ifin = MIN(ideb+n_temp-1,n-1)
  n_temp = ifin-ideb+1
  ! On verifie que n_temp n'est pas nul
  IF (n_temp <= 0) THEN
    CALL jmerreur3(nomsp,6,n,nr,nwork)
  END IF
  nwork_temp = n_temp*nr

  ! Gestion de fin
  IF (ifin == n-1) THEN
    fin = .TRUE.
  ELSE
    fin = .FALSE.
  END IF

END SUBROUTINE jmdecoup

! Pour tronconner en dimension 3 de facon a tenir dans le nwork disponible
SUBROUTINE jmdecoup3(n,m,nmr,nwork,debut,lpair,ideb,ifin,jdeb,jfin,nmtemp,nwork_temp,fini)

  IMPLICIT NONE

  ! Arguments
  INTEGER, INTENT(in) :: n, m, nmr, nwork
  LOGICAL, INTENT(in) :: debut, lpair
  INTEGER, INTENT(out) :: nmtemp, nwork_temp
  INTEGER, INTENT(out)   :: ideb, jdeb
  INTEGER, INTENT(inout) :: ifin, jfin
  LOGICAL, INTENT(out) :: fini

  ! Variables locales
  INTEGER :: ijdeb, ijfin
  CHARACTER(len=*), PARAMETER :: nomsp = 'JMDECOUP3'

  ! n*m*nr est l'espace total qu'il faudrait pour work.
  ! Malheureusement, on n'a que nwork au plus
  ! On va donc decouper n et m en morceaux pour tenir

  ! Gestion de debut
  IF (debut) THEN
    ideb = 0
    jdeb = 0
  ELSE
    IF (ifin < n-1) THEN
      ideb = ifin+1
      jdeb = jfin
    ELSE
      ideb = 0
      jdeb = jfin+1
    END IF
  END IF

  ! Gestion de nmtemp
  nmtemp = nwork/nmr
  ! Si l impair, on doit eviter que nmtemp soit impair (routine cs et sc)
  IF (.NOT.lpair .AND. MOD(nmtemp,2) /= 0) nmtemp = nmtemp-1
  ! Pour simplifier, on passe par des indices 2d
  ijdeb = ideb+jdeb*n
  ijfin = MIN(ijdeb+nmtemp-1,n*m-1)
  nmtemp = ijfin-ijdeb+1
  ! On verifie que nmtemp n'est pas nul
  IF (nmtemp <= 0) THEN
    CALL jmerreur4(nomsp,6,n,m,nmr,nwork)
  END IF
  nwork_temp = nmtemp*nmr

  ! On deduit ifin et jfin de ijfin
  jfin = ijfin/n
  ifin = ijfin-n*jfin

  ! Gestion de fin
  IF (ifin == n-1 .AND. jfin == m-1) THEN
    fini = .TRUE.
  ELSE
    fini = .FALSE.
  END IF

END SUBROUTINE jmdecoup3

SUBROUTINE jmccm1d(m,n,fact,nfact,ifact,table,ntable,itable,work,nwork,ioff)

  IMPLICIT NONE

  ! Arguments
  INTEGER, INTENT(in) :: m, n
  INTEGER, INTENT(in) :: nfact, ifact
  INTEGER, INTENT(in), DIMENSION(0:nfact-1) :: fact
  INTEGER, INTENT(in) :: ntable,itable
  REAL(kind=GFT_prec), INTENT(in), DIMENSION(0:ntable-1) :: table
  INTEGER, INTENT(in) :: nwork
  REAL(kind=GFT_prec), INTENT(inout), DIMENSION(0:nwork-1) :: work
  INTEGER, INTENT(inout) :: ioff

  ! Variables locales
  INTEGER :: nterms
  INTEGER :: np, pp, lastnp, premier
  INTEGER :: nprod, nprod1, nprod2
  INTEGER :: n2, p2, n3, p3, n5, p5
  INTEGER :: i

  ! On recupere les facteurs
  nterms = fact(ifact)
  n2 = fact(ifact+1)
  p2 = fact(ifact+2)
  n3 = fact(ifact+3)
  p3 = fact(ifact+4)
  n5 = fact(ifact+5)
  p5 = fact(ifact+6)
  nprod = n2*n3*n5
  DO i = 7,nterms-1,2
    nprod = nprod*fact(ifact+i+1)**fact(ifact+i)
  END DO

  ! On fait n3*n5 T.F. de n2 (qui est en puissances de 2)
  IF (n2 /= 1) THEN
    CALL jmccm1d2(p2,n2,m*(nprod/n2),table,ntable,itable,n,n/n2,work,nwork,ioff)
  END IF

  ! On transpose (on tient compte de ioff) en permutant les deux parties
  ! On en profite pour multiplier par le bon wij
  IF (n2 /= 1 .AND. nprod /= n2) THEN
    CALL jmcctranspcs(m,n,n2,nprod/n2,table,ntable,itable,work,nwork,ioff)
  END IF

  ! On fait n5*n2 T.F. de n3 (en puissances de 3)
  IF (n3 /= 1) THEN
    CALL jmccm1d3(p3,n3,m*(nprod/n3),table,ntable,itable,n,n/n3,work,nwork,ioff)
  END IF

  ! On transpose (on tient compte de ioff) en permutant les deux parties
  ! On en profite pour multiplier par le bon wij
  IF (n3 /= 1 .AND. nprod /= n3) THEN
    CALL jmcctranspcs(m*n2,n,n3,nprod/(n2*n3), &
    & table,ntable,itable,work,nwork,ioff)
  END IF

  ! On fait n2*n3 T.F. de n5 (en puissances de 5)
  IF (n5 /= 1) THEN
    CALL jmccm1d5(p5,n5,m*(nprod/n5),table,ntable,itable,n,n/n5,work,nwork,ioff)
  END IF

  ! On transpose s'il y a lieu (si on a fait quelque chose et s'il reste des
  ! termes a traiter
  IF (n5 /= 1 .AND. nprod /= n5 .AND. nterms > 7) THEN
    CALL jmcctranspcs(m*n2*n3,n,n5,nprod/(n2*n3*n5), &
    & table,ntable,itable,work,nwork,ioff)
  END IF
  nprod1 = m*n2*n3
  nprod2 = n2*n3*n5
  lastnp = n5

  ! On passe aux nombres premiers autres que 2, 3 et 5
  DO i = 7,nterms-1,2

    pp = fact(ifact+i)
    premier = fact(ifact+i+1)
    np = premier**pp

    CALL jmccm1dp(premier,pp,m*(nprod/np), &
    & table,ntable,itable,n,n/np,work,nwork,ioff)

    nprod1 = nprod1 * lastnp
    nprod2 = nprod2 * np
    IF (np /= 1 .AND. nprod /= np .AND. nterms > i+1) THEN
      CALL jmcctranspcs(nprod1,n,np,nprod/nprod2, &
      & table,ntable,itable,work,nwork,ioff)
    END IF
    lastnp = np

  END DO

END SUBROUTINE jmccm1d

SUBROUTINE jmccm1d2(p,n,m,table,ntable,itable,ntable2,mtable,work,nwork,ioff)

  IMPLICIT NONE

  ! Arguments
  INTEGER, INTENT(in) :: p, n, m
  INTEGER, INTENT(in) :: ntable,itable,ntable2,mtable
  REAL(kind=GFT_prec), INTENT(in), DIMENSION(0:ntable-1) :: table
  INTEGER, INTENT(in) :: nwork
  REAL(kind=GFT_prec), INTENT(inout), DIMENSION(0:nwork-1) :: work
  INTEGER, INTENT(inout) :: ioff

  ! Variables locales
  INTEGER :: k, jl
  INTEGER :: it1,iu1,it2,iu2
  INTEGER :: jt1,ju1,jt2,ju2
  REAL(kind=GFT_prec) :: x1, x2, y1, y2
  REAL(kind=GFT_prec) :: c, s
  INTEGER :: ioff1, ioff2

  ! On joue sur ioff pour alterner entre le haut et le bas de work
  ioff1 = ioff
  ioff2 = nwork/2-ioff1

  IF (MOD(p,2)==0) THEN

    ! Si p est pair, on peut travailler entierement en base 4
    CALL jmccm1d4(p/2,n,m,table,ntable,itable,ntable2,mtable,work,nwork,ioff1)
    ioff = ioff1

  ELSE

    ! On fait les premieres etapes en base 4
    CALL jmccm1d4(p/2,n,2*m,table,ntable,itable,ntable2,mtable*2,work,nwork,ioff1)
    ioff2 = nwork/2-ioff1
    ! On fait la derniere etape en base 2
    IF (m >= 16 .OR. 2**(p-1) < 8) THEN
      DO k = 0,2**(p-1)-1

        ! Les sinus et cosinus
        c = table(itable+        mtable*k)
        s = table(itable+ntable2+mtable*k)

        ! Les indices
        it1 = ioff1        +m*(k*2  )
        iu1 = ioff1+nwork/4+m*(k*2  )
        it2 = ioff1        +m*(k*2+1)
        iu2 = ioff1+nwork/4+m*(k*2+1)
        jt1 = ioff2        +m*( k          )
        ju1 = ioff2+nwork/4+m*( k          )
        jt2 = ioff2        +m*((k+2**(p-1)))
        ju2 = ioff2+nwork/4+m*((k+2**(p-1)))

!dir$ ivdep
!ocl novrec
!cdir nodep
        DO jl = 0,m-1
          x1 = work(it1+jl)
          y1 = work(iu1+jl)
          x2 = work(it2+jl)
          y2 = work(iu2+jl)
          work(jt1+jl) = x1 + ( x2*c - y2*s )
          work(ju1+jl) = y1 + ( x2*s + y2*c )
          work(jt2+jl) = x1 - ( x2*c - y2*s )
          work(ju2+jl) = y1 - ( x2*s + y2*c )
        END DO
      END DO
    ELSE
      DO jl = 0,m-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        DO k = 0,2**(p-1)-1
          x1 = work(ioff1+jl        +m*(k*2  ))
          y1 = work(ioff1+jl+nwork/4+m*(k*2  ))
          x2 = work(ioff1+jl        +m*(k*2+1))
          y2 = work(ioff1+jl+nwork/4+m*(k*2+1))
          ! Les sinus et cosinus
          c = table(itable+        mtable*k)
          s = table(itable+ntable2+mtable*k)
          work(ioff2+jl        +m*( k          )) = x1 + ( x2*c - y2*s )
          work(ioff2+jl+nwork/4+m*( k          )) = y1 + ( x2*s + y2*c )
          work(ioff2+jl        +m*((k+2**(p-1)))) = x1 - ( x2*c - y2*s )
          work(ioff2+jl+nwork/4+m*((k+2**(p-1)))) = y1 - ( x2*s + y2*c )
        END DO
      END DO
    END IF

    ioff = ioff2

  END IF

END SUBROUTINE jmccm1d2

SUBROUTINE jmccm1d3(p,n,m,table,ntable,itable,ntable2,mtable,work,nwork,ioff)

  IMPLICIT NONE

  ! Arguments
  INTEGER, INTENT(in) :: p, n, m
  INTEGER, INTENT(in) :: ntable,itable,ntable2, mtable
  REAL(kind=GFT_prec), INTENT(in), DIMENSION(0:ntable-1) :: table
  INTEGER, INTENT(in) :: nwork
  REAL(kind=GFT_prec), INTENT(inout), DIMENSION(0:nwork-1) :: work
  INTEGER, INTENT(inout) :: ioff

  ! Variables locales
  INTEGER :: i, k, jl
  REAL(kind=GFT_prec) :: x1, x2, x3, y1, y2, y3, t1, t2, t3, u1, u2, u3
  REAL(kind=GFT_prec) :: c2, s2, c3, s3
  INTEGER :: it1,iu1,it2,iu2,it3,iu3
  INTEGER :: jt1,ju1,jt2,ju2,jt3,ju3
  REAL(kind=GFT_prec) :: r,s,t,u

  ! Gestion des constantes cosinus
  REAL(kind=GFT_prec), SAVE :: ctwopi3, stwopi3
  LOGICAL, SAVE :: first = .TRUE.

  INTEGER :: ioff1, ioff2

  ! On recupere cos et sin de 2*pi/3
  IF (first) THEN
    first = .FALSE.
    ctwopi3 = -1.0_GFT_prec/2.0_GFT_prec
    stwopi3 = SQRT(3.0_GFT_prec)/2.0_GFT_prec
  END IF

  ! On joue sur ioff pour alterner entre le haut et le bas de work
  ioff1 = ioff
  ioff2 = nwork/2-ioff1

  ! Boucle sur les etapes
  DO i = 0, p-1

    IF (m*3**(p-i-1) >= 16 .OR. 3**i < 8) THEN

      DO k = 0,3**i-1

        ! Les sinus et cosinus
        c2 = table(itable+        mtable*  3**(p-i-1)*k)
        s2 = table(itable+ntable2+mtable*  3**(p-i-1)*k)
        c3 = table(itable+        mtable*2*3**(p-i-1)*k)
        s3 = table(itable+ntable2+mtable*2*3**(p-i-1)*k)

        ! Les indices
        it1 = ioff1        +m*(k*3**(p-i)             )
        iu1 = ioff1+nwork/4+m*(k*3**(p-i)             )
        it2 = ioff1        +m*(k*3**(p-i)+  3**(p-i-1))
        iu2 = ioff1+nwork/4+m*(k*3**(p-i)+  3**(p-i-1))
        it3 = ioff1        +m*(k*3**(p-i)+2*3**(p-i-1))
        iu3 = ioff1+nwork/4+m*(k*3**(p-i)+2*3**(p-i-1))
        jt1 = ioff2        +m*( k        *3**(p-i-1))
        ju1 = ioff2+nwork/4+m*( k        *3**(p-i-1))
        jt2 = ioff2        +m*((k+  3**i)*3**(p-i-1))
        ju2 = ioff2+nwork/4+m*((k+  3**i)*3**(p-i-1))
        jt3 = ioff2        +m*((k+2*3**i)*3**(p-i-1))
        ju3 = ioff2+nwork/4+m*((k+2*3**i)*3**(p-i-1))

!dir$ ivdep
!ocl novrec
!cdir nodep
        DO jl = 0,m*3**(p-i-1)-1

          r = (c2*work(it2+jl))-(s2*work(iu2+jl))
          s = (c2*work(iu2+jl))+(s2*work(it2+jl))
          t = (c3*work(it3+jl))-(s3*work(iu3+jl))
          u = (c3*work(iu3+jl))+(s3*work(it3+jl))
          x1 = work(it1+jl)
          y1 = work(iu1+jl)
          work(jt1+jl) = x1 + r + t
          work(ju1+jl) = y1 + s + u
          work(jt2+jl) = x1 + ctwopi3*(r+t) - stwopi3*(s-u)
          work(ju2+jl) = y1 + ctwopi3*(s+u) + stwopi3*(r-t)
          work(jt3+jl) = x1 + ctwopi3*(r+t) + stwopi3*(s-u)
          work(ju3+jl) = y1 + ctwopi3*(s+u) - stwopi3*(r-t)

        END DO

      END DO

    ELSE

      DO jl = 0,m*3**(p-i-1)-1

!dir$ ivdep
!ocl novrec
!cdir nodep
        DO k = 0,3**i-1

          t1 = work(ioff1+jl        +m*(k*3**(p-i)             ))
          u1 = work(ioff1+jl+nwork/4+m*(k*3**(p-i)             ))
          t2 = work(ioff1+jl        +m*(k*3**(p-i)+  3**(p-i-1)))
          u2 = work(ioff1+jl+nwork/4+m*(k*3**(p-i)+  3**(p-i-1)))
          t3 = work(ioff1+jl        +m*(k*3**(p-i)+2*3**(p-i-1)))
          u3 = work(ioff1+jl+nwork/4+m*(k*3**(p-i)+2*3**(p-i-1)))

          ! Les sinus et cosinus
          c2 = table(itable+        mtable*  3**(p-i-1)*k)
          s2 = table(itable+ntable2+mtable*  3**(p-i-1)*k)
          c3 = table(itable+        mtable*2*3**(p-i-1)*k)
          s3 = table(itable+ntable2+mtable*2*3**(p-i-1)*k)

          ! On premultiplie
          x1 = t1
          y1 = u1
          x2 = c2*t2-s2*u2
          y2 = c2*u2+s2*t2
          x3 = c3*t3-s3*u3
          y3 = c3*u3+s3*t3

          ! Il reste a multiplier par les twopi3
          work(ioff2+jl        +m*( k        *3**(p-i-1))) = &
          & x1 + x2                    + x3
          work(ioff2+jl+nwork/4+m*( k        *3**(p-i-1))) = &
          & y1 + y2                    + y3
          work(ioff2+jl        +m*((k+  3**i)*3**(p-i-1))) = &
          & x1 + ctwopi3*x2-stwopi3*y2 + ctwopi3*x3+stwopi3*y3
          work(ioff2+jl+nwork/4+m*((k+  3**i)*3**(p-i-1))) = &
          & y1 + ctwopi3*y2+stwopi3*x2 + ctwopi3*y3-stwopi3*x3
          work(ioff2+jl        +m*((k+2*3**i)*3**(p-i-1))) = &
          & x1 + ctwopi3*x2+stwopi3*y2 + ctwopi3*x3-stwopi3*y3
          work(ioff2+jl+nwork/4+m*((k+2*3**i)*3**(p-i-1))) = &
          & y1 + ctwopi3*y2-stwopi3*x2 + ctwopi3*y3+stwopi3*x3

        END DO

      END DO

    END IF

    ! On alterne les offsets
    ioff1 = nwork/2-ioff1
    ioff2 = nwork/2-ioff2

  ! Fin boucle sur le nombre d'etapes
  END DO

  ioff = ioff1

END SUBROUTINE jmccm1d3

SUBROUTINE jmccm1d4(p,n,m,table,ntable,itable,ntable2,mtable,work,nwork,ioff)

  IMPLICIT NONE

  ! Arguments
  INTEGER, INTENT(in) :: p, n, m
  INTEGER, INTENT(in) :: ntable,itable,ntable2, mtable
  REAL(kind=GFT_prec), INTENT(in), DIMENSION(0:ntable-1) :: table
  INTEGER, INTENT(in) :: nwork
  REAL(kind=GFT_prec), INTENT(inout), DIMENSION(0:nwork-1) :: work
  INTEGER, INTENT(inout) :: ioff

  ! Variables locales
  INTEGER :: i, k, jl
  REAL(kind=GFT_prec) :: x0,x1,x2,x3,y0,y1,y2,y3,t0,t1,t2,t3,u0,u1,u2,u3
  REAL(kind=GFT_prec) :: x0px2,x0mx2,x1px3,x1mx3
  REAL(kind=GFT_prec) :: y0py2,y0my2,y1py3,y1my3
  REAL(kind=GFT_prec) :: c1, s1, c2, s2, c3, s3
  INTEGER :: ioff1, ioff2
  INTEGER :: it0,iu0,it1,iu1,it2,iu2,it3,iu3
  INTEGER :: jt0,ju0,jt1,ju1,jt2,ju2,jt3,ju3

  ! On joue sur ioff pour alterner entre le haut et le bas de work
  ioff1 = ioff
  ioff2 = nwork/2-ioff1

  ! Boucle sur les etapes
  DO i = 0, p-1

    IF (m*4**(p-i-1) >= 16 .OR. 4**i < 8) THEN

      DO k = 0,4**i-1

        ! Les sinus et cosinus
        c1 = table(itable+        mtable*  4**(p-i-1)*k)
        s1 = table(itable+ntable2+mtable*  4**(p-i-1)*k)
        c2 = table(itable+        mtable*2*4**(p-i-1)*k)
        s2 = table(itable+ntable2+mtable*2*4**(p-i-1)*k)
        c3 = table(itable+        mtable*3*4**(p-i-1)*k)
        s3 = table(itable+ntable2+mtable*3*4**(p-i-1)*k)

        ! Les indices
        it0 = ioff1        +m*(k*4**(p-i)             )
        iu0 = ioff1+nwork/4+m*(k*4**(p-i)             )
        it1 = ioff1        +m*(k*4**(p-i)+  4**(p-i-1))
        iu1 = ioff1+nwork/4+m*(k*4**(p-i)+  4**(p-i-1))
        it2 = ioff1        +m*(k*4**(p-i)+2*4**(p-i-1))
        iu2 = ioff1+nwork/4+m*(k*4**(p-i)+2*4**(p-i-1))
        it3 = ioff1        +m*(k*4**(p-i)+3*4**(p-i-1))
        iu3 = ioff1+nwork/4+m*(k*4**(p-i)+3*4**(p-i-1))
        jt0 = ioff2        +m*( k        *4**(p-i-1))
        ju0 = ioff2+nwork/4+m*( k        *4**(p-i-1))
        jt1 = ioff2        +m*((k+  4**i)*4**(p-i-1))
        ju1 = ioff2+nwork/4+m*((k+  4**i)*4**(p-i-1))
        jt2 = ioff2        +m*((k+2*4**i)*4**(p-i-1))
        ju2 = ioff2+nwork/4+m*((k+2*4**i)*4**(p-i-1))
        jt3 = ioff2        +m*((k+3*4**i)*4**(p-i-1))
        ju3 = ioff2+nwork/4+m*((k+3*4**i)*4**(p-i-1))

!dir$ ivdep
!ocl novrec
!cdir nodep
        DO jl = 0,m*4**(p-i-1)-1

          x0px2 = work(it0+jl) + (c2*work(it2+jl)-s2*work(iu2+jl))
          x0mx2 = work(it0+jl) - (c2*work(it2+jl)-s2*work(iu2+jl))
          y0py2 = work(iu0+jl) + (c2*work(iu2+jl)+s2*work(it2+jl))
          y0my2 = work(iu0+jl) - (c2*work(iu2+jl)+s2*work(it2+jl))
          x1px3 = (c1*work(it1+jl)-s1*work(iu1+jl))+(c3*work(it3+jl)-s3*work(iu3+jl))
          x1mx3 = (c1*work(it1+jl)-s1*work(iu1+jl))-(c3*work(it3+jl)-s3*work(iu3+jl))
          y1py3 = (c1*work(iu1+jl)+s1*work(it1+jl))+(c3*work(iu3+jl)+s3*work(it3+jl))
          y1my3 = (c1*work(iu1+jl)+s1*work(it1+jl))-(c3*work(iu3+jl)+s3*work(it3+jl))

          ! Il reste a multiplier par les twopi4
          work(jt0+jl) = (x0px2)+(x1px3)
          work(jt2+jl) = (x0px2)-(x1px3)
          work(ju0+jl) = (y0py2)+(y1py3)
          work(ju2+jl) = (y0py2)-(y1py3)
          work(jt1+jl) = (x0mx2)-(y1my3)
          work(jt3+jl) = (x0mx2)+(y1my3)
          work(ju1+jl) = (y0my2)+(x1mx3)
          work(ju3+jl) = (y0my2)-(x1mx3)

        END DO

      END DO

    ELSE

      DO jl = 0,m*4**(p-i-1)-1

!dir$ ivdep
!ocl novrec
!cdir nodep
        DO k = 0,4**i-1

          t0 = work(ioff1+jl        +m*(k*4**(p-i)             ))
          u0 = work(ioff1+jl+nwork/4+m*(k*4**(p-i)             ))
          t1 = work(ioff1+jl        +m*(k*4**(p-i)+  4**(p-i-1)))
          u1 = work(ioff1+jl+nwork/4+m*(k*4**(p-i)+  4**(p-i-1)))
          t2 = work(ioff1+jl        +m*(k*4**(p-i)+2*4**(p-i-1)))
          u2 = work(ioff1+jl+nwork/4+m*(k*4**(p-i)+2*4**(p-i-1)))
          t3 = work(ioff1+jl        +m*(k*4**(p-i)+3*4**(p-i-1)))
          u3 = work(ioff1+jl+nwork/4+m*(k*4**(p-i)+3*4**(p-i-1)))

          ! Les sinus et cosinus
          c1 = table(itable+        mtable*  4**(p-i-1)*k)
          s1 = table(itable+ntable2+mtable*  4**(p-i-1)*k)
          c2 = table(itable+        mtable*2*4**(p-i-1)*k)
          s2 = table(itable+ntable2+mtable*2*4**(p-i-1)*k)
          c3 = table(itable+        mtable*3*4**(p-i-1)*k)
          s3 = table(itable+ntable2+mtable*3*4**(p-i-1)*k)

          ! On premultiplie
          x0 = t0
          y0 = u0
          x1 = c1*t1-s1*u1
          y1 = c1*u1+s1*t1
          x2 = c2*t2-s2*u2
          y2 = c2*u2+s2*t2
          x3 = c3*t3-s3*u3
          y3 = c3*u3+s3*t3

          ! Il reste a multiplier par les twopi4
          work(ioff2+jl        +m*( k        *4**(p-i-1))) = x0+x1+x2+x3
          work(ioff2+jl+nwork/4+m*( k        *4**(p-i-1))) = y0+y1+y2+y3
          work(ioff2+jl        +m*((k+  4**i)*4**(p-i-1))) = x0-y1-x2+y3
          work(ioff2+jl+nwork/4+m*((k+  4**i)*4**(p-i-1))) = y0+x1-y2-x3
          work(ioff2+jl        +m*((k+2*4**i)*4**(p-i-1))) = x0-x1+x2-x3
          work(ioff2+jl+nwork/4+m*((k+2*4**i)*4**(p-i-1))) = y0-y1+y2-y3
          work(ioff2+jl        +m*((k+3*4**i)*4**(p-i-1))) = x0+y1-x2-y3
          work(ioff2+jl+nwork/4+m*((k+3*4**i)*4**(p-i-1))) = y0-x1-y2+x3

        END DO

      END DO

    END IF

    ! On alterne les offsets
    ioff1 = nwork/2-ioff1
    ioff2 = nwork/2-ioff2

  ! Fin boucle sur le nombre d'etapes
  END DO

  ioff = ioff1

END SUBROUTINE jmccm1d4

SUBROUTINE jmccm1d5(p,n,m,table,ntable,itable,ntable2,mtable,work,nwork,ioff)

  IMPLICIT NONE

  ! Arguments
  INTEGER, INTENT(in) :: p, n, m
  INTEGER, INTENT(in) :: ntable,itable,ntable2, mtable
  REAL(kind=GFT_prec), INTENT(in), DIMENSION(0:ntable-1) :: table
  INTEGER, INTENT(in) :: nwork
  REAL(kind=GFT_prec), INTENT(inout), DIMENSION(0:nwork-1) :: work
  INTEGER, INTENT(inout) :: ioff

  ! Variables locales
  INTEGER :: i, k, jl
  REAL(kind=GFT_prec) :: x0,x1,x2,x3,x4,y0,y1,y2,y3,y4,t0,t1,t2,t3,t4,u0,u1,u2,u3,u4
  REAL(kind=GFT_prec) :: c1, s1, c2, s2, c3, s3, c4, s4
  INTEGER :: ioff1, ioff2

  ! Gestion des constantes cosinus
  REAL(kind=GFT_prec), SAVE :: twopi5
  REAL(kind=GFT_prec), SAVE :: ctwopi51, ctwopi52, ctwopi53, ctwopi54
  REAL(kind=GFT_prec), SAVE :: stwopi51, stwopi52, stwopi53, stwopi54
  LOGICAL, SAVE :: first = .TRUE.

  ! On recupere cos et sin de 2*pi/5
  IF (first) THEN
    first = .FALSE.
    twopi5   = 2.0_GFT_prec*ACOS(-1.0_GFT_prec)/5.0_GFT_prec
    ctwopi51 = COS(  twopi5)
    stwopi51 = SIN(  twopi5)
    ctwopi52 = COS(2.0_GFT_prec*twopi5)
    stwopi52 = SIN(2.0_GFT_prec*twopi5)
    ctwopi53 = COS(3.0_GFT_prec*twopi5)
    stwopi53 = SIN(3.0_GFT_prec*twopi5)
    ctwopi54 = COS(4.0_GFT_prec*twopi5)
    stwopi54 = SIN(4.0_GFT_prec*twopi5)
  END IF

  ! On joue sur ioff pour alterner entre le haut et le bas de work
  ioff1 = ioff
  ioff2 = nwork/2-ioff1

  ! Boucle sur les etapes
  DO i = 0, p-1

    IF (m*5**(p-i-1) >= 16 .OR. 5**i < 8) THEN

      DO k = 0,5**i-1

!dir$ ivdep
!ocl novrec
!cdir nodep
        DO jl = 0,m*5**(p-i-1)-1

          t0 = work(ioff1+jl        +m*(k*5**(p-i)             ))
          u0 = work(ioff1+jl+nwork/4+m*(k*5**(p-i)             ))
          t1 = work(ioff1+jl        +m*(k*5**(p-i)+  5**(p-i-1)))
          u1 = work(ioff1+jl+nwork/4+m*(k*5**(p-i)+  5**(p-i-1)))
          t2 = work(ioff1+jl        +m*(k*5**(p-i)+2*5**(p-i-1)))
          u2 = work(ioff1+jl+nwork/4+m*(k*5**(p-i)+2*5**(p-i-1)))
          t3 = work(ioff1+jl        +m*(k*5**(p-i)+3*5**(p-i-1)))
          u3 = work(ioff1+jl+nwork/4+m*(k*5**(p-i)+3*5**(p-i-1)))
          t4 = work(ioff1+jl        +m*(k*5**(p-i)+4*5**(p-i-1)))
          u4 = work(ioff1+jl+nwork/4+m*(k*5**(p-i)+4*5**(p-i-1)))

          ! Les sinus et cosinus
          c1 = table(itable+        mtable*  5**(p-i-1)*k)
          s1 = table(itable+ntable2+mtable*  5**(p-i-1)*k)
          c2 = table(itable+        mtable*2*5**(p-i-1)*k)
          s2 = table(itable+ntable2+mtable*2*5**(p-i-1)*k)
          c3 = table(itable+        mtable*3*5**(p-i-1)*k)
          s3 = table(itable+ntable2+mtable*3*5**(p-i-1)*k)
          c4 = table(itable+        mtable*4*5**(p-i-1)*k)
          s4 = table(itable+ntable2+mtable*4*5**(p-i-1)*k)

          ! On premultiplie
          x0 = t0
          y0 = u0
          x1 = c1*t1-s1*u1
          y1 = c1*u1+s1*t1
          x2 = c2*t2-s2*u2
          y2 = c2*u2+s2*t2
          x3 = c3*t3-s3*u3
          y3 = c3*u3+s3*t3
          x4 = c4*t4-s4*u4
          y4 = c4*u4+s4*t4

          ! Il reste a multiplier par les twopi5
          work(ioff2+jl        +m*( k        *5**(p-i-1))) =   &
          & x0 + x1                    + x2                    &
          &    + x3                    + x4
          work(ioff2+jl+nwork/4+m*( k        *5**(p-i-1))) =   &
          & y0 + y1                    + y2                    &
          &    + y3                    + y4
          work(ioff2+jl        +m*((k+  5**i)*5**(p-i-1))) =   &
          & x0 + ctwopi51*x1 - stwopi51*y1 &
          &    + ctwopi52*x2 - stwopi52*y2 &
          &    + ctwopi53*x3 - stwopi53*y3 &
          &    + ctwopi54*x4 - stwopi54*y4
          work(ioff2+jl+nwork/4+m*((k+  5**i)*5**(p-i-1))) =   &
          & y0 + ctwopi51*y1 + stwopi51*x1 &
          &    + ctwopi52*y2 + stwopi52*x2 &
          &    + ctwopi53*y3 + stwopi53*x3 &
          &    + ctwopi54*y4 + stwopi54*x4
          work(ioff2+jl        +m*((k+2*5**i)*5**(p-i-1))) =   &
          & x0 + ctwopi52*x1 - stwopi52*y1 &
          &    + ctwopi54*x2 - stwopi54*y2 &
          &    + ctwopi51*x3 - stwopi51*y3 &
          &    + ctwopi53*x4 - stwopi53*y4
          work(ioff2+jl+nwork/4+m*((k+2*5**i)*5**(p-i-1))) =   &
          & y0 + ctwopi52*y1 + stwopi52*x1 &
          &    + ctwopi54*y2 + stwopi54*x2 &
          &    + ctwopi51*y3 + stwopi51*x3 &
          &    + ctwopi53*y4 + stwopi53*x4
          work(ioff2+jl        +m*((k+3*5**i)*5**(p-i-1))) =   &
          & x0 + ctwopi53*x1 - stwopi53*y1 &
          &    + ctwopi51*x2 - stwopi51*y2 &
          &    + ctwopi54*x3 - stwopi54*y3 &
          &    + ctwopi52*x4 - stwopi52*y4
          work(ioff2+jl+nwork/4+m*((k+3*5**i)*5**(p-i-1))) =   &
          & y0 + ctwopi53*y1 + stwopi53*x1 &
          &    + ctwopi51*y2 + stwopi51*x2 &
          &    + ctwopi54*y3 + stwopi54*x3 &
          &    + ctwopi52*y4 + stwopi52*x4
          work(ioff2+jl        +m*((k+4*5**i)*5**(p-i-1))) =   &
          & x0 + ctwopi54*x1 - stwopi54*y1 &
          &    + ctwopi53*x2 - stwopi53*y2 &
          &    + ctwopi52*x3 - stwopi52*y3 &
          &    + ctwopi51*x4 - stwopi51*y4
          work(ioff2+jl+nwork/4+m*((k+4*5**i)*5**(p-i-1))) =   &
          & y0 + ctwopi54*y1 + stwopi54*x1 &
          &    + ctwopi53*y2 + stwopi53*x2 &
          &    + ctwopi52*y3 + stwopi52*x3 &
          &    + ctwopi51*y4 + stwopi51*x4

        END DO

      END DO

    ELSE

      DO jl = 0,m*5**(p-i-1)-1

!dir$ ivdep
!ocl novrec
!cdir nodep
        DO k = 0,5**i-1

          t0 = work(ioff1+jl        +m*(k*5**(p-i)             ))
          u0 = work(ioff1+jl+nwork/4+m*(k*5**(p-i)             ))
          t1 = work(ioff1+jl        +m*(k*5**(p-i)+  5**(p-i-1)))
          u1 = work(ioff1+jl+nwork/4+m*(k*5**(p-i)+  5**(p-i-1)))
          t2 = work(ioff1+jl        +m*(k*5**(p-i)+2*5**(p-i-1)))
          u2 = work(ioff1+jl+nwork/4+m*(k*5**(p-i)+2*5**(p-i-1)))
          t3 = work(ioff1+jl        +m*(k*5**(p-i)+3*5**(p-i-1)))
          u3 = work(ioff1+jl+nwork/4+m*(k*5**(p-i)+3*5**(p-i-1)))
          t4 = work(ioff1+jl        +m*(k*5**(p-i)+4*5**(p-i-1)))
          u4 = work(ioff1+jl+nwork/4+m*(k*5**(p-i)+4*5**(p-i-1)))

          ! Les sinus et cosinus
          c1 = table(itable+        mtable*  5**(p-i-1)*k)
          s1 = table(itable+ntable2+mtable*  5**(p-i-1)*k)
          c2 = table(itable+        mtable*2*5**(p-i-1)*k)
          s2 = table(itable+ntable2+mtable*2*5**(p-i-1)*k)
          c3 = table(itable+        mtable*3*5**(p-i-1)*k)
          s3 = table(itable+ntable2+mtable*3*5**(p-i-1)*k)
          c4 = table(itable+        mtable*4*5**(p-i-1)*k)
          s4 = table(itable+ntable2+mtable*4*5**(p-i-1)*k)

          ! On premultiplie
          x0 = t0
          y0 = u0
          x1 = c1*t1-s1*u1
          y1 = c1*u1+s1*t1
          x2 = c2*t2-s2*u2
          y2 = c2*u2+s2*t2
          x3 = c3*t3-s3*u3
          y3 = c3*u3+s3*t3
          x4 = c4*t4-s4*u4
          y4 = c4*u4+s4*t4

          ! Il reste a multiplier par les twopi5
          work(ioff2+jl        +m*( k        *5**(p-i-1))) =   &
          & x0 + x1                    + x2                    &
          &    + x3                    + x4
          work(ioff2+jl+nwork/4+m*( k        *5**(p-i-1))) =   &
          & y0 + y1                    + y2                    &
          &    + y3                    + y4
          work(ioff2+jl        +m*((k+  5**i)*5**(p-i-1))) =   &
          & x0 + ctwopi51*x1 - stwopi51*y1 &
          &    + ctwopi52*x2 - stwopi52*y2 &
          &    + ctwopi53*x3 - stwopi53*y3 &
          &    + ctwopi54*x4 - stwopi54*y4
          work(ioff2+jl+nwork/4+m*((k+  5**i)*5**(p-i-1))) =   &
          & y0 + ctwopi51*y1 + stwopi51*x1 &
          &    + ctwopi52*y2 + stwopi52*x2 &
          &    + ctwopi53*y3 + stwopi53*x3 &
          &    + ctwopi54*y4 + stwopi54*x4
          work(ioff2+jl        +m*((k+2*5**i)*5**(p-i-1))) =   &
          & x0 + ctwopi52*x1 - stwopi52*y1 &
          &    + ctwopi54*x2 - stwopi54*y2 &
          &    + ctwopi51*x3 - stwopi51*y3 &
          &    + ctwopi53*x4 - stwopi53*y4
          work(ioff2+jl+nwork/4+m*((k+2*5**i)*5**(p-i-1))) =   &
          & y0 + ctwopi52*y1 + stwopi52*x1 &
          &    + ctwopi54*y2 + stwopi54*x2 &
          &    + ctwopi51*y3 + stwopi51*x3 &
          &    + ctwopi53*y4 + stwopi53*x4
          work(ioff2+jl        +m*((k+3*5**i)*5**(p-i-1))) =   &
          & x0 + ctwopi53*x1 - stwopi53*y1 &
          &    + ctwopi51*x2 - stwopi51*y2 &
          &    + ctwopi54*x3 - stwopi54*y3 &
          &    + ctwopi52*x4 - stwopi52*y4
          work(ioff2+jl+nwork/4+m*((k+3*5**i)*5**(p-i-1))) =   &
          & y0 + ctwopi53*y1 + stwopi53*x1 &
          &    + ctwopi51*y2 + stwopi51*x2 &
          &    + ctwopi54*y3 + stwopi54*x3 &
          &    + ctwopi52*y4 + stwopi52*x4
          work(ioff2+jl        +m*((k+4*5**i)*5**(p-i-1))) =   &
          & x0 + ctwopi54*x1 - stwopi54*y1 &
          &    + ctwopi53*x2 - stwopi53*y2 &
          &    + ctwopi52*x3 - stwopi52*y3 &
          &    + ctwopi51*x4 - stwopi51*y4
          work(ioff2+jl+nwork/4+m*((k+4*5**i)*5**(p-i-1))) =   &
          & y0 + ctwopi54*y1 + stwopi54*x1 &
          &    + ctwopi53*y2 + stwopi53*x2 &
          &    + ctwopi52*y3 + stwopi52*x3 &
          &    + ctwopi51*y4 + stwopi51*x4

        END DO

      END DO

    END IF

    ! On alterne les offsets
    ioff1 = nwork/2-ioff1
    ioff2 = nwork/2-ioff2

  ! Fin boucle sur le nombre d'etapes
  END DO

  ioff = ioff1

END SUBROUTINE jmccm1d5

SUBROUTINE jmccm1dp(p,q,m,table,ntable,itable,ntable2,mtable,work,nwork,ioff)

  ! On fait m t.f. d'ordre q en base p (p**q)
  ! Note : n n'est pas utilise ici

  IMPLICIT NONE

  ! Arguments
  INTEGER, INTENT(in) :: p, q, m
  INTEGER, INTENT(in) :: ntable,itable,ntable2, mtable
  REAL(kind=GFT_prec), INTENT(in), DIMENSION(0:ntable-1) :: table
  INTEGER, INTENT(in) :: nwork
  REAL(kind=GFT_prec), INTENT(inout), DIMENSION(0:nwork-1) :: work
  INTEGER, INTENT(inout) :: ioff

  ! Variables locales
  INTEGER :: i, k, jl, jp, kp
  REAL(kind=GFT_prec) :: ck, sk, tk, uk, cpjk, spjk
  INTEGER :: pqq, pi, pqi, pqii
  INTEGER :: ikpr, ikpi, ijpr, ijpi
  INTEGER :: itr, iti, jtr, jti
  INTEGER :: ioff1, ioff2
  REAL(kind=GFT_prec) :: c11, c12, c21, c22

  ! On joue sur ioff pour alterner entre le haut et le bas de work
  ioff1 = ioff
  ioff2 = nwork/2-ioff1

  ! Pour le calcul des cos(2*pi/p)
  pqq = p**(q-1)

  ! Boucle sur les etapes
  DO i = 0, q-1

    pi   = p**i
    pqi  = p**(q-i)
    pqii = p**(q-i-1)

    DO k = 0,pi-1

      DO jp = 0,p-1

        DO jl = 0,m*pqii-1

          ijpr = ioff2 + jl + m*((k+jp*pi)*pqii)
          ijpi = ijpr + nwork/4

          work(ijpr) = 0
          work(ijpi) = 0

        END DO

      END DO

      DO kp = 0,p-1

        itr = itable+mtable*kp*pqii*k
        iti = itr + ntable2
        ck = table(itr)
        sk = table(iti)

        DO jp = 0,p-1

          ! Gymanstique infernale pour recuperer cos(2*pi/p) etc
          jtr = itable+mtable*pqq*MOD(jp*kp,p)
          jti = jtr + ntable2
          cpjk = table(jtr)
          spjk = table(jti)
          c11 = (cpjk*ck-spjk*sk)
          c12 = (cpjk*sk+spjk*ck)
          c21 = (cpjk*sk+spjk*ck)
          c22 = (cpjk*ck-spjk*sk)

!dir$ ivdep
!ocl novrec
!cdir nodep
          DO jl = 0,m*pqii-1

            ikpr = ioff1+jl+m*(k*pqi+kp*pqii)
            ikpi = ikpr + nwork/4
            tk = work(ikpr)
            uk = work(ikpi)

            ijpr = ioff2+jl+m*((k+jp*pi)*pqii)
            ijpi = ijpr + nwork/4

            work(ijpr) = work(ijpr) + tk*c11-uk*c12
            work(ijpi) = work(ijpi) + tk*c21+uk*c22

          END DO

        END DO

      END DO

    END DO

    ! On alterne les offsets
    ioff1 = nwork/2 - ioff1
    ioff2 = nwork/2 - ioff2

  ! Fin boucle sur le nombre d'etapes
  END DO

  ioff = ioff1

END SUBROUTINE jmccm1dp

SUBROUTINE jmcsm1d(m,n,fact,nfact,ifact,table,ntable,itable,work,nwork,ioff)

  IMPLICIT NONE

  ! Arguments
  INTEGER, INTENT(in) :: m, n
  INTEGER, INTENT(in) :: nfact, ifact
  INTEGER, INTENT(inout), DIMENSION(0:nfact-1) :: fact
  INTEGER, INTENT(in) :: ntable,itable
  REAL(kind=GFT_prec), INTENT(inout), DIMENSION(0:ntable-1) :: table
  INTEGER, INTENT(in) :: nwork
  REAL(kind=GFT_prec), INTENT(inout), DIMENSION(0:nwork-1) :: work
  INTEGER, INTENT(inout) :: ioff

  ! Variables locales
  INTEGER :: ioff1, ioff2
  INTEGER :: i, j
  REAL(kind=GFT_prec) :: t, u, v, w, tt, uu, vv, ww
  REAL(kind=GFT_prec) :: c, s
  INTEGER :: it

  ! Gestion de work
  ioff1 = ioff
  ioff2 = nwork/2 - ioff1

  ! On doit faire m T.F. complexes -> reelles de longueur n
  ! Si m est pair
  IF (MOD(m,2) == 0) THEN

    ! On distribue

    IF (m/2 >= 16 .OR. n/2 < 8) THEN

      DO i = 0,n/2
!dir$ ivdep
!ocl novrec
!cdir nodep
        DO j = 0,m/2-1
          it = n-i
          IF (i == 0) it = 0
          t = work(ioff1        +i*m+j    )
          u = work(ioff1+nwork/4+i*m+j    )
          v = work(ioff1        +i*m+j+m/2)
          w = work(ioff1+nwork/4+i*m+j+m/2)
          work(ioff2        + i*m/2+j) = (t-w)
          work(ioff2+nwork/4+ i*m/2+j) = (u+v)
          work(ioff2        +it*m/2+j) = (t+w)
          work(ioff2+nwork/4+it*m/2+j) = (v-u)
        END DO
      END DO

    ELSE

      DO j = 0,m/2-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        DO i = 0,n/2
          it = n-i
          IF (i == 0) it = 0
          t = work(ioff1        +i*m+j    )
          u = work(ioff1+nwork/4+i*m+j    )
          v = work(ioff1        +i*m+j+m/2)
          w = work(ioff1+nwork/4+i*m+j+m/2)
          work(ioff2        + i*m/2+j) = (t-w)
          work(ioff2+nwork/4+ i*m/2+j) = (u+v)
          work(ioff2        +it*m/2+j) = (t+w)
          work(ioff2+nwork/4+it*m/2+j) = (v-u)
        END DO
      END DO

    END IF

    ! On fait m/2 t.f. complexes -> complexes de longueur n
    CALL jmccm1d(m/2,n,fact,nfact,ifact,table,ntable,itable,work,nwork,ioff2)
    ioff1 = nwork/2 - ioff2

    ! On reconstitue

    IF (m/2 >= 16 .OR. n < 8) THEN

      DO i = 0,n-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        DO j = 0,m/2-1
          work(ioff1+i*m+j    ) = work(ioff2        +i*m/2+j)
          work(ioff1+i*m+j+m/2) = work(ioff2+nwork/4+i*m/2+j)
        END DO
      END DO

    ELSE

      DO j = 0,m/2-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        DO i = 0,n-1
          work(ioff1+i*m+j    ) = work(ioff2        +i*m/2+j)
          work(ioff1+i*m+j+m/2) = work(ioff2+nwork/4+i*m/2+j)
        END DO
      END DO

    END IF

  ! Si m n'est pas pair mais que n l'est
  ELSE IF (MOD(n,2) == 0) THEN

    ! On distribue

    IF (m >= 16 .OR. n/2 < 8) THEN

      DO i = 0,n/2-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        DO j = 0,m-1
          ! Note : Signe - sur les parties imaginaires pour inversion
          t =  work(ioff1        +      i*m+j)
          u = -work(ioff1+nwork/4+      i*m+j)
          v =  work(ioff1        +(n/2-i)*m+j)
          w = -work(ioff1+nwork/4+(n/2-i)*m+j)
          c = table(itable+i)
          s = table(itable+i+n)
          tt = (t+v)/2.0_GFT_prec
          uu = (u-w)/2.0_GFT_prec
          vv = (c*(t-v)+s*(u+w))/2.0_GFT_prec
          ww = (c*(u+w)-s*(t-v))/2.0_GFT_prec
          ! Note : le facteur 2 et le signe - viennent de l'inversion Fourier
          work(ioff2        +m*i+j) =  2.0_GFT_prec*(tt-ww)
          work(ioff2+nwork/4+m*i+j) = -2.0_GFT_prec*(uu+vv)
        END DO
      END DO

    ELSE

      DO j = 0,m-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        DO i = 0,n/2-1
          ! Note : Signe - sur les parties imaginaires pour inversion
          t =  work(ioff1        +      i*m+j)
          u = -work(ioff1+nwork/4+      i*m+j)
          v =  work(ioff1        +(n/2-i)*m+j)
          w = -work(ioff1+nwork/4+(n/2-i)*m+j)
          c = table(itable+i)
          s = table(itable+i+n)
          tt = (t+v)/2.0_GFT_prec
          uu = (u-w)/2.0_GFT_prec
          vv = (c*(t-v)+s*(u+w))/2.0_GFT_prec
          ww = (c*(u+w)-s*(t-v))/2.0_GFT_prec
          ! Note : le facteur 2 et le signe - viennent de l'inversion Fourier
          work(ioff2        +m*i+j) =  2.0_GFT_prec*(tt-ww)
          work(ioff2+nwork/4+m*i+j) = -2.0_GFT_prec*(uu+vv)
        END DO
      END DO

    END IF

    ! On fait m t.f. complexes de taille n/2
    fact(ifact+1) = fact(ifact+1)/2 ! Revient a remplacer n2 par n2/2
    fact(ifact+2) = fact(ifact+2)-1 ! Revient a remplacer p2 par p2-1
    CALL jmccm1d(m,n,fact,nfact,ifact,table,ntable,itable,work,nwork,ioff2)
    fact(ifact+1) = fact(ifact+1)*2 ! On retablit les valeurs initiales
    fact(ifact+2) = fact(ifact+2)+1
    ioff1 = nwork/2 - ioff2

    ! On reconstitue

    IF (m >= 16 .OR. n/2 < 8) THEN

      DO i = 0, n/2-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        DO j = 0, m-1
          ! Note : le signe - vient de l'inversion
          work(ioff1+m*(2*i  )+j) =  work(ioff2        +m*i+j)
          work(ioff1+m*(2*i+1)+j) = -work(ioff2+nwork/4+m*i+j)
        END DO
      END DO

    ELSE

      DO j = 0, m-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        DO i = 0, n/2-1
          ! Note : le signe - vient de l'inversion
          work(ioff1+m*(2*i  )+j) =  work(ioff2        +m*i+j)
          work(ioff1+m*(2*i+1)+j) = -work(ioff2+nwork/4+m*i+j)
        END DO
      END DO

    END IF

  END IF

  ioff = ioff1

END SUBROUTINE jmcsm1d

! Variante de jmcsm1d ou on fournit x en entree et en sortie
SUBROUTINE jmcsm1dxy(m,n,fact,nfact,ifact,table,ntable,itable,work,nwork,x,dimx,debx,incx,jumpx,y,dimy,deby,incy,jumpy,isign,scale)

  IMPLICIT NONE

  ! Arguments
  INTEGER, INTENT(in) :: m, n
  INTEGER, INTENT(in) :: nfact, ifact
  INTEGER, INTENT(inout), DIMENSION(0:nfact-1) :: fact
  INTEGER, INTENT(in) :: ntable,itable
  REAL(kind=GFT_prec), INTENT(in), DIMENSION(0:ntable-1) :: table
  INTEGER, INTENT(in) :: nwork
  REAL(kind=GFT_prec), INTENT(inout), DIMENSION(0:nwork-1) :: work
  INTEGER, INTENT(in) :: dimx, debx, incx, jumpx
  INTEGER, INTENT(in) :: dimy, deby, incy, jumpy
  REAL(kind=GFT_prec), INTENT(in),  DIMENSION(0:dimx-1) :: x
  REAL(kind=GFT_prec), INTENT(out), DIMENSION(0:dimy-1) :: y
  INTEGER, INTENT(in) :: isign
  REAL(kind=GFT_prec), INTENT(in) :: scale

  ! Variables locales
  INTEGER :: i, j
  REAL(kind=GFT_prec) :: t, u, v, w, tt, uu, vv, ww
  REAL(kind=GFT_prec) :: c, s
  INTEGER :: it
  INTEGER :: ioff

  ! On doit faire m T.F. complexes -> reelles de longueur n
  ! Si m est pair
  IF (MOD(m,2) == 0) THEN

    ! On distribue

    IF (m/2 >= 16 .OR. n/2 < 8) THEN

      DO i = 0,n/2
!dir$ ivdep
!ocl novrec
!cdir nodep
        DO j = 0,m/2-1
          it = n-i
          IF (i == 0) it = 0
          t =       scale*x(debx+incx*(2*i  )+jumpx*(j    ))
          u = isign*scale*x(debx+incx*(2*i+1)+jumpx*(j    ))
          v =       scale*x(debx+incx*(2*i  )+jumpx*(j+m/2))
          w = isign*scale*x(debx+incx*(2*i+1)+jumpx*(j+m/2))
          work(         i*m/2+j) = (t-w)
          work(nwork/4+ i*m/2+j) = (u+v)
          work(        it*m/2+j) = (t+w)
          work(nwork/4+it*m/2+j) = (v-u)
        END DO
      END DO

    ELSE

      DO j = 0,m/2-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        DO i = 0,n/2
          it = n-i
          IF (i == 0) it = 0
          t =       scale*x(debx+incx*(2*i  )+jumpx*(j    ))
          u = isign*scale*x(debx+incx*(2*i+1)+jumpx*(j    ))
          v =       scale*x(debx+incx*(2*i  )+jumpx*(j+m/2))
          w = isign*scale*x(debx+incx*(2*i+1)+jumpx*(j+m/2))
          work(         i*m/2+j) = (t-w)
          work(nwork/4+ i*m/2+j) = (u+v)
          work(        it*m/2+j) = (t+w)
          work(nwork/4+it*m/2+j) = (v-u)
        END DO
      END DO

    END IF

    ! On fait m/2 t.f. complexes -> complexes de longueur n
    ioff = 0
    CALL jmccm1d(m/2,n,fact,nfact,ifact,table,ntable,itable,work,nwork,ioff)

    ! On reconstitue

    IF (m/2 >= 16 .OR. n < 8) THEN

      DO i = 0,n-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        DO j = 0,m/2-1
          y(deby+jumpy*(j    )+incy*i) = work(ioff        +i*m/2+j)
          y(deby+jumpy*(j+m/2)+incy*i) = work(ioff+nwork/4+i*m/2+j)
        END DO
      END DO

    ELSE

      DO j = 0,m/2-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        DO i = 0,n-1
          y(deby+jumpy*(j    )+incy*i) = work(ioff        +i*m/2+j)
          y(deby+jumpy*(j+m/2)+incy*i) = work(ioff+nwork/4+i*m/2+j)
        END DO
      END DO

    END IF

  ! Si m n'est pas pair mais que n l'est
  ELSE IF (MOD(n,2) == 0) THEN

    ! On distribue

    IF (m >= 16 .OR. n/2 < 8) THEN

      DO i = 0,n/2-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        DO j = 0,m-1
          ! Note : Signe - sur les parties imaginaires pour inversion
          t =        scale*x(debx+(2*i)  *incx+j*jumpx)
          u = -isign*scale*x(debx+(2*i+1)*incx+j*jumpx)
          v =        scale*x(debx+(2*(n/2-i)  )*incx+j*jumpx)
          w = -isign*scale*x(debx+(2*(n/2-i)+1)*incx+j*jumpx)
          c = table(itable+i)
          s = table(itable+i+n)
          tt = (t+v)/2.0_GFT_prec
          uu = (u-w)/2.0_GFT_prec
          vv = (c*(t-v)+s*(u+w))/2.0_GFT_prec
          ww = (c*(u+w)-s*(t-v))/2.0_GFT_prec
          ! Note : le facteur 2 et le signe - viennent de l'inversion Fourier
          work(        m*i+j) =  2.0_GFT_prec*(tt-ww)
          work(nwork/4+m*i+j) = -2.0_GFT_prec*(uu+vv)
        END DO
      END DO

    ELSE

      DO j = 0,m-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        DO i = 0,n/2-1
          ! Note : Signe - sur les parties imaginaires pour inversion
          t =        scale*x(debx+(2*i)  *incx+j*jumpx)
          u = -isign*scale*x(debx+(2*i+1)*incx+j*jumpx)
          v =        scale*x(debx+(2*(n/2-i)  )*incx+j*jumpx)
          w = -isign*scale*x(debx+(2*(n/2-i)+1)*incx+j*jumpx)
          c = table(itable+i)
          s = table(itable+i+n)
          tt = (t+v)/2.0_GFT_prec
          uu = (u-w)/2.0_GFT_prec
          vv = (c*(t-v)+s*(u+w))/2.0_GFT_prec
          ww = (c*(u+w)-s*(t-v))/2.0_GFT_prec
          ! Note : le facteur 2 et le signe - viennent de l'inversion Fourier
          work(        m*i+j) =  2.0_GFT_prec*(tt-ww)
          work(nwork/4+m*i+j) = -2.0_GFT_prec*(uu+vv)
        END DO
      END DO

    END IF

    ! On fait m t.f. complexes de taille n/2
    ioff = 0
    fact(ifact+1) = fact(ifact+1)/2 ! Revient a remplacer n2 par n2/2
    fact(ifact+2) = fact(ifact+2)-1 ! Revient a remplacer p2 par p2-1
    CALL jmccm1d(m,n,fact,nfact,ifact,table,ntable,itable,work,nwork,ioff)
    fact(ifact+1) = fact(ifact+1)*2 ! On retablit les valeurs initiales
    fact(ifact+2) = fact(ifact+2)+1

    ! On reconstitue

    IF (m >= 16 .OR. n/2 < 8) THEN

      DO i = 0, n/2-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        DO j = 0, m-1
          ! Note : le signe - vient de l'inversion
          y(deby+incy*(2*i  )+jumpy*j) =  work(ioff        +m*i+j)
          y(deby+incy*(2*i+1)+jumpy*j) = -work(ioff+nwork/4+m*i+j)
        END DO
      END DO

    ELSE

!dir$ ivdep
      DO j = 0, m-1
!ocl novrec
!cdir nodep
        DO i = 0, n/2-1
          ! Note : le signe - vient de l'inversion
          y(deby+incy*(2*i  )+jumpy*j) =  work(ioff        +m*i+j)
          y(deby+incy*(2*i+1)+jumpy*j) = -work(ioff+nwork/4+m*i+j)
        END DO
      END DO

    END IF

  END IF

END SUBROUTINE jmcsm1dxy

SUBROUTINE jmscm1d(m,n,fact,nfact,ifact,table,ntable,itable,work,nwork,ioff)

  IMPLICIT NONE

  ! Arguments
  INTEGER, INTENT(in) :: m, n
  INTEGER, INTENT(in) :: nfact, ifact
  INTEGER, INTENT(inout), DIMENSION(0:nfact-1) :: fact
  INTEGER, INTENT(in) :: ntable,itable
  REAL(kind=GFT_prec), INTENT(in), DIMENSION(0:ntable-1) :: table
  INTEGER, INTENT(in) :: nwork
  REAL(kind=GFT_prec), INTENT(inout), DIMENSION(0:nwork-1) :: work
  INTEGER, INTENT(inout) :: ioff

  ! Variables locales
  INTEGER :: ioff1, ioff2
  INTEGER :: i, j
  REAL(kind=GFT_prec) :: t, u, v, w
  REAL(kind=GFT_prec) :: c, s
  INTEGER :: is, it

  ! Gestion de work
  ioff1 = ioff
  ioff2 = nwork/2 - ioff1

  ! On doit faire m T.F. reelles de longueur n
  ! Si m est pair
  IF (MOD(m,2) == 0) THEN

    ! On distribue
    IF (m/2 >= 16 .OR. n < 8) THEN

      DO i = 0,n-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        DO j = 0,m/2-1
          work(ioff2        +i*m/2+j) = work(ioff1+i*m+j    )
          work(ioff2+nwork/4+i*m/2+j) = work(ioff1+i*m+j+m/2)
        END DO
      END DO

    ELSE

      DO j = 0,m/2-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        DO i = 0,n-1
          work(ioff2        +i*m/2+j) = work(ioff1+i*m+j    )
          work(ioff2+nwork/4+i*m/2+j) = work(ioff1+i*m+j+m/2)
        END DO
      END DO

    END IF
        
    ! On fait m/2 t.f. complexes de longueur n
    CALL jmccm1d(m/2,n,fact,nfact,ifact,table,ntable,itable,work,nwork,ioff2)
    ioff1 = nwork/2 - ioff2

    ! On regenere le resultat
    IF (m/2 >= 16 .OR. n/2 < 8) THEN

      DO i = 0,n/2
!dir$ ivdep
!ocl novrec
!cdir nodep
        DO j = 0,m/2-1
          it = n-i
          IF (i == 0) it = 0
          t = work(ioff2        + i*m/2+j)
          u = work(ioff2+nwork/4+ i*m/2+j)
          v = work(ioff2        +it*m/2+j)
          w = work(ioff2+nwork/4+it*m/2+j)
          work(ioff1        +i*m+j    ) = (t+v)/2.0_GFT_prec
          work(ioff1+nwork/4+i*m+j    ) = (u-w)/2.0_GFT_prec
          work(ioff1        +i*m+j+m/2) = (u+w)/2.0_GFT_prec
          work(ioff1+nwork/4+i*m+j+m/2) = (v-t)/2.0_GFT_prec
        END DO
      END DO

    ELSE

      DO j = 0,m/2-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        DO i = 0,n/2
          it = n-i
          IF (i == 0) it = 0
          t = work(ioff2        + i*m/2+j)
          u = work(ioff2+nwork/4+ i*m/2+j)
          v = work(ioff2        +it*m/2+j)
          w = work(ioff2+nwork/4+it*m/2+j)
          work(ioff1        +i*m+j    ) = (t+v)/2.0_GFT_prec
          work(ioff1+nwork/4+i*m+j    ) = (u-w)/2.0_GFT_prec
          work(ioff1        +i*m+j+m/2) = (u+w)/2.0_GFT_prec
          work(ioff1+nwork/4+i*m+j+m/2) = (v-t)/2.0_GFT_prec
        END DO
      END DO

    END IF

  ! Si m n'est pas pair mais que n l'est
  ELSE IF (MOD(n,2) == 0) THEN

    ! On distribue les indices pairs et impairs selon n

    IF (m >= 16 .OR. n/2 < 8) THEN

      DO i = 0, n/2-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        DO j = 0, m-1
          work(ioff2+        m*i+j) = work(ioff1+m*(2*i  )+j)
          work(ioff2+nwork/4+m*i+j) = work(ioff1+m*(2*i+1)+j)
        END DO
      END DO

    ELSE

      DO j = 0, m-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        DO i = 0, n/2-1
          work(ioff2+        m*i+j) = work(ioff1+m*(2*i  )+j)
          work(ioff2+nwork/4+m*i+j) = work(ioff1+m*(2*i+1)+j)
        END DO
      END DO

    END IF

    ! On fait m t.f. complexes de taille n/2
    fact(ifact+1) = fact(ifact+1)/2 ! Revient a remplacer n2 par n2/2
    fact(ifact+2) = fact(ifact+2)-1 ! Revient a remplacer p2 par p2-1
    CALL jmccm1d(m,n,fact,nfact,ifact,table,ntable,itable,work,nwork,ioff2)
    fact(ifact+1) = fact(ifact+1)*2 ! On retablit les valeurs initiales
    fact(ifact+2) = fact(ifact+2)+1
    ioff1 = nwork/2 - ioff2

    ! Maintenant, il faut reconstituer la t.f. reelle

    IF (m >= 16 .OR. n/2 < 8) THEN

      DO i = 0,n/2
!dir$ ivdep
!ocl novrec
!cdir nodep
        DO j = 0,m-1
          is = i
          it = n/2-i
          IF (i == 0 .OR. i == n/2) THEN
            is = 0
            it = 0
          END IF
          t = work(ioff2        +is*m+j)
          u = work(ioff2+nwork/4+is*m+j)
          v = work(ioff2        +it*m+j)
          w = work(ioff2+nwork/4+it*m+j)
          c = table(itable+i)
          s = table(itable+i+n)
          work(ioff1        +i*m+j) = (t+v)/2.0_GFT_prec + c*(u+w)/2.0_GFT_prec - s*(v-t)/2.0_GFT_prec
          work(ioff1+nwork/4+i*m+j) = (u-w)/2.0_GFT_prec + c*(v-t)/2.0_GFT_prec + s*(u+w)/2.0_GFT_prec
        END DO
      END DO

    ELSE

      DO j = 0,m-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        DO i = 0,n/2
          is = i
          it = n/2-i
          IF (i == 0 .OR. i == n/2) THEN
            is = 0
            it = 0
          END IF
          t = work(ioff2        +is*m+j)
          u = work(ioff2+nwork/4+is*m+j)
          v = work(ioff2        +it*m+j)
          w = work(ioff2+nwork/4+it*m+j)
          c = table(itable+i)
          s = table(itable+i+n)
          work(ioff1        +i*m+j) = (t+v)/2.0_GFT_prec + c*(u+w)/2.0_GFT_prec - s*(v-t)/2.0_GFT_prec
          work(ioff1+nwork/4+i*m+j) = (u-w)/2.0_GFT_prec + c*(v-t)/2.0_GFT_prec + s*(u+w)/2.0_GFT_prec
        END DO
      END DO

    END IF

  END IF

  ioff = ioff1

END SUBROUTINE jmscm1d

! Variante de jmscm1d ou on fournit x en entree et en sortie
SUBROUTINE jmscm1dxy(m,n,fact,nfact,ifact,table,ntable,itable,work,nwork,x,dimx,debx,incx,jumpx,y,dimy,deby,incy,jumpy,isign,scale)

  IMPLICIT NONE

  ! Arguments
  INTEGER, INTENT(in) :: m, n
  INTEGER, INTENT(in) :: nfact, ifact
  INTEGER, INTENT(inout), DIMENSION(0:nfact-1) :: fact
  INTEGER, INTENT(in) :: ntable,itable
  REAL(kind=GFT_prec), INTENT(in), DIMENSION(0:ntable-1) :: table
  INTEGER, INTENT(in) :: nwork
  REAL(kind=GFT_prec), INTENT(inout), DIMENSION(0:nwork-1) :: work
  INTEGER, INTENT(in) :: dimx, debx, incx, jumpx
  INTEGER, INTENT(in) :: dimy, deby, incy, jumpy
  REAL(kind=GFT_prec), INTENT(in) :: scale
  REAL(kind=GFT_prec), INTENT(in),  DIMENSION(0:dimx-1) :: x
  REAL(kind=GFT_prec), INTENT(out), DIMENSION(0:dimy-1) :: y
  INTEGER, INTENT(in) :: isign

  ! Variables locales
  INTEGER :: i, j
  REAL(kind=GFT_prec) :: t, u, v, w
  REAL(kind=GFT_prec) :: c, s
  INTEGER :: is, it
  INTEGER :: ioff

  ! On doit faire m T.F. reelles de longueur n
  ! Si m est pair
  IF (MOD(m,2) == 0) THEN

    ! On distribue
    IF (m/2 >= 16 .OR. n < 8) THEN

      DO i = 0,n-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        DO j = 0,m/2-1
          work(        i*m/2+j) = x(debx+i*incx+(j)    *jumpx)
          work(nwork/4+i*m/2+j) = x(debx+i*incx+(j+m/2)*jumpx)
        END DO
      END DO

    ELSE

      DO j = 0,m/2-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        DO i = 0,n-1
          work(        i*m/2+j) = x(debx+i*incx+(j)    *jumpx)
          work(nwork/4+i*m/2+j) = x(debx+i*incx+(j+m/2)*jumpx)
        END DO
      END DO

    END IF

    ! On fait m/2 t.f. complexes de longueur n
    ioff = 0
    CALL jmccm1d(m/2,n,fact,nfact,ifact,table,ntable,itable,work,nwork,ioff)

    ! On regenere le resultat
    IF (m/2 >= 16 .OR. n/2 < 8) THEN

      DO i = 0,n/2
!dir$ ivdep
!ocl novrec
!cdir nodep
        DO j = 0,m/2-1
          it = n-i
          IF (i == 0) it = 0
          t = work(ioff        + i*m/2+j)
          u = work(ioff+nwork/4+ i*m/2+j)
          v = work(ioff        +it*m/2+j)
          w = work(ioff+nwork/4+it*m/2+j)
          y(deby+(2*i)  *incy+(j)    *jumpy) =       scale*(t+v)/2.0_GFT_prec
          y(deby+(2*i+1)*incy+(j)    *jumpy) = isign*scale*(u-w)/2.0_GFT_prec
          y(deby+(2*i)  *incy+(j+m/2)*jumpy) =       scale*(u+w)/2.0_GFT_prec
          y(deby+(2*i+1)*incy+(j+m/2)*jumpy) = isign*scale*(v-t)/2.0_GFT_prec
        END DO
      END DO

    ELSE

      DO j = 0,m/2-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        DO i = 0,n/2
          it = n-i
          IF (i == 0) it = 0
          t = work(ioff        + i*m/2+j)
          u = work(ioff+nwork/4+ i*m/2+j)
          v = work(ioff        +it*m/2+j)
          w = work(ioff+nwork/4+it*m/2+j)
          y(deby+(2*i)  *incy+(j)    *jumpy) =       scale*(t+v)/2.0_GFT_prec
          y(deby+(2*i+1)*incy+(j)    *jumpy) = isign*scale*(u-w)/2.0_GFT_prec
          y(deby+(2*i)  *incy+(j+m/2)*jumpy) =       scale*(u+w)/2.0_GFT_prec
          y(deby+(2*i+1)*incy+(j+m/2)*jumpy) = isign*scale*(v-t)/2.0_GFT_prec
        END DO
      END DO

    END IF

  ! Si m n'est pas pair mais que n l'est
  ELSE IF (MOD(n,2) == 0) THEN

    ! On distribue les indices pairs et impairs selon n

    IF (m >= 16 .OR. n/2 < 8) THEN

      DO i = 0, n/2-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        DO j = 0, m-1
          work(        m*i+j) = x(debx+incx*(2*i  )+jumpx*j)
          work(nwork/4+m*i+j) = x(debx+incx*(2*i+1)+jumpx*j)
        END DO
      END DO

    ELSE

      DO j = 0, m-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        DO i = 0, n/2-1
          work(        m*i+j) = x(debx+incx*(2*i  )+jumpx*j)
          work(nwork/4+m*i+j) = x(debx+incx*(2*i+1)+jumpx*j)
        END DO
      END DO

    END IF

    ! On fait m t.f. complexes de taille n/2
    ioff = 0
    fact(ifact+1) = fact(ifact+1)/2 ! Revient a remplacer n2 par n2/2
    fact(ifact+2) = fact(ifact+2)-1 ! Revient a remplacer p2 par p2-1
    CALL jmccm1d(m,n,fact,nfact,ifact,table,ntable,itable,work,nwork,ioff)
    fact(ifact+1) = fact(ifact+1)*2 ! On retablit les valeurs initiales
    fact(ifact+2) = fact(ifact+2)+1

    ! Maintenant, il faut reconstituer la t.f. reelle

    IF (m >= 16 .OR. n/2 < 8) THEN

      DO i = 0,n/2
!dir$ ivdep
!ocl novrec
!cdir nodep
        DO j = 0,m-1
          is = i
          it = n/2-i
          IF (i == 0 .OR. i == n/2) THEN
            is = 0
            it = 0
          END IF
          t = work(ioff        +is*m+j)
          u = work(ioff+nwork/4+is*m+j)
          v = work(ioff        +it*m+j)
          w = work(ioff+nwork/4+it*m+j)
          c = table(itable+i)
          s = table(itable+i+n)
          y(deby+(2*i  )*incy+j*jumpy) = &
          &       scale*((t+v)/2.0_GFT_prec + c*(u+w)/2.0_GFT_prec - s*(v-t)/2.0_GFT_prec)
          y(deby+(2*i+1)*incy+j*jumpy) = &
          & isign*scale*((u-w)/2.0_GFT_prec + c*(v-t)/2.0_GFT_prec + s*(u+w)/2.0_GFT_prec)
        END DO
      END DO

    ELSE

      DO j = 0,m-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        DO i = 0,n/2
          is = i
          it = n/2-i
          IF (i == 0 .OR. i == n/2) THEN
            is = 0
            it = 0
          END IF
          t = work(ioff        +is*m+j)
          u = work(ioff+nwork/4+is*m+j)
          v = work(ioff        +it*m+j)
          w = work(ioff+nwork/4+it*m+j)
          c = table(itable+i)
          s = table(itable+i+n)
          y(deby+(2*i  )*incy+j*jumpy) = &
          &       scale*((t+v)/2.0_GFT_prec + c*(u+w)/2.0_GFT_prec - s*(v-t)/2.0_GFT_prec)
          y(deby+(2*i+1)*incy+j*jumpy) = &
          & isign*scale*((u-w)/2.0_GFT_prec + c*(v-t)/2.0_GFT_prec + s*(u+w)/2.0_GFT_prec)
        END DO
      END DO

    END IF

  END IF

END SUBROUTINE jmscm1dxy

SUBROUTINE jmtransp(n,m,l,work,nwork,ioff)

  IMPLICIT NONE

  ! Arguments
  INTEGER, INTENT(in) :: n, m, l
  INTEGER, INTENT(in) :: nwork
  REAL(kind=GFT_prec), INTENT(inout), DIMENSION(0:nwork-1) :: work
  INTEGER, INTENT(inout) :: ioff

  ! Variables locales
  INTEGER :: ioff1, ioff2
  INTEGER :: ij, k

  ioff1 = ioff
  ioff2 = nwork/2-ioff1

  ! On transpose (nm)(l) en (l)(nm) en distinguant les parties reelles et im.
  IF (m*n >= 16 .OR. l < 8) THEN

    DO k = 0,l-1
!dir$ ivdep
!ocl novrec
!cdir nodep
      DO ij = 0,m*n-1
        work(ioff2+      ij*l+k) = work(ioff1+      k*n*m+ij)
        work(ioff2+n*m*l+ij*l+k) = work(ioff1+n*m*l+k*n*m+ij)
      END DO
    END DO

  ELSE

    DO ij = 0,m*n-1
!dir$ ivdep
!ocl novrec
!cdir nodep
      DO k = 0,l-1
        work(ioff2+      ij*l+k) = work(ioff1+      k*n*m+ij)
        work(ioff2+n*m*l+ij*l+k) = work(ioff1+n*m*l+k*n*m+ij)
      END DO
    END DO

  END IF

  ioff = ioff2

END SUBROUTINE jmtransp

SUBROUTINE jmcctranspcs(m,n,n2,n3,table,ntable,itable,work,nwork,ioff)

  ! Cette subroutine permute le contenu du tableau work de la facon suivante
  ! On considere qu'a l'origine ce tableau est en (m,n3,n2)
  ! On doit transposer le terme d'ordre (k,j,i) en (k,i,j)
  ! On en profite pour faire les multiplications par wij
  ! Le role de n est seulement d'attaquer les bonnes valeurs du tableau table
  ! (il y a un ecart de n entre les cos et les sin, et le stride entre
  !  les cos est de n/(n2*n3)
  ! Note : le sens de n2 et n3 ici n'a rien a voir avec celui de jmccm1d

  IMPLICIT NONE

  ! Arguments
  INTEGER, INTENT(in) :: m, n
  INTEGER, INTENT(in) :: n2, n3
  INTEGER, INTENT(in) :: ntable,itable
  REAL(kind=GFT_prec), INTENT(in),  DIMENSION(0:ntable-1) :: table
  INTEGER, INTENT(in) :: nwork
  REAL(kind=GFT_prec), INTENT(inout),  DIMENSION(0:nwork-1) :: work
  INTEGER, INTENT(inout) :: ioff

  ! Variables locales
  INTEGER :: i, j, k
  REAL(kind=GFT_prec) :: t, u, c, s
  INTEGER :: ioff1, ioff2
  INTEGER :: is

  ! Gestion des offsets
  ioff1 = ioff
  ioff2 = nwork/2-ioff

  ! Gestion du stride
  is = n/(n2*n3)

  IF ( m >= 16 .OR. (n2 < 8 .AND. n3 < 8) ) THEN

    DO i = 0,n2-1
      DO j = 0,n3-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        DO k = 0,m-1
          t = work(ioff1        +k+m*(j+n3*i))
          u = work(ioff1+nwork/4+k+m*(j+n3*i))
          c = table(itable+  is*i*j)
          s = table(itable+n+is*i*j)
          work(ioff2        +k+m*(i+n2*j)) = c*t-s*u
          work(ioff2+nwork/4+k+m*(i+n2*j)) = c*u+s*t
        END DO
      END DO
    END DO

  ELSE IF ( n2 >= 16 .OR. n3 < 8 ) THEN

    DO j = 0,n3-1
      DO k = 0,m-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        DO i = 0,n2-1
          t = work(ioff1        +k+m*(j+n3*i))
          u = work(ioff1+nwork/4+k+m*(j+n3*i))
          c = table(itable+  is*i*j)
          s = table(itable+n+is*i*j)
          work(ioff2        +k+m*(i+n2*j)) = c*t-s*u
          work(ioff2+nwork/4+k+m*(i+n2*j)) = c*u+s*t
        END DO
      END DO
    END DO

  ELSE

    DO i = 0,n2-1
      DO k = 0,m-1
!dir$ ivdep
!ocl novrec
!cdir nodep
        DO j = 0,n3-1
          t = work(ioff1        +k+m*(j+n3*i))
          u = work(ioff1+nwork/4+k+m*(j+n3*i))
          c = table(itable+  is*i*j)
          s = table(itable+n+is*i*j)
          work(ioff2        +k+m*(i+n2*j)) = c*t-s*u
          work(ioff2+nwork/4+k+m*(i+n2*j)) = c*u+s*t
        END DO
      END DO
    END DO

  END IF

  ioff = ioff2

END SUBROUTINE jmcctranspcs

SUBROUTINE jmsetnwork(nwork)

  ! Subroutine appelee par l'utilisateur pour augmenter le nwork
  ! des routines 2d et 3d

  IMPLICIT NONE

  ! Arguments
  INTEGER, INTENT(in) :: nwork

  ! Variables locales
  CHARACTER(len=*), PARAMETER :: nomsp = 'JMSETNWORK'
  INTEGER :: nwork2

  IF (nwork <= 0) THEN
    CALL jmerreur1(nomsp,4,nwork)
  END IF

  nwork2 = nwork
  CALL jmgetsetnwork(nwork2,'s')

END SUBROUTINE jmsetnwork

SUBROUTINE jmgetnwork(nwork,nwork_def,nwork_min)

  ! On recupere la valeur de nwork si elle a ete augmentee par l'utilisateur
  ! Sinon on prend la valeur par defaut
  ! Il s'agit du nwork des routines 2d et 3d

  IMPLICIT NONE

  ! Arguments
  INTEGER, INTENT(out) :: nwork
  INTEGER, INTENT(in)  :: nwork_def, nwork_min

  ! Variables locales
  INTEGER :: nwork_loc
  CHARACTER(len=*), PARAMETER :: nomsp = 'JMGETNWORK'

  CALL jmgetsetnwork(nwork_loc,'g')

  ! Valeur par defaut
  IF (nwork_loc == -1) THEN
    nwork = nwork_def
  ! Valeur invalide (trop petite)
  ELSE IF (nwork_loc < nwork_min) THEN
    CALL jmerreur2(nomsp,5,nwork_loc,nwork_min)
  ! Valeur correcte
  ELSE
    nwork = nwork_loc
  END IF

END SUBROUTINE jmgetnwork


SUBROUTINE jmgetsetnwork(nwork,TYPE)

  ! Subroutine qui permet de stocker une valeur statique
  ! Ceci evite de recourir a un common ...

  IMPLICIT NONE

  ! Arguments
  INTEGER, INTENT(inout) :: nwork
  CHARACTER(len=1), INTENT(in) :: TYPE

  ! Variables locales

  ! Variable statique
  INTEGER, SAVE :: nwork_last = -1

  IF (TYPE == 's') THEN
    nwork_last = nwork
  ELSE IF (TYPE == 'g') THEN 
    nwork = nwork_last
  END IF

END SUBROUTINE jmgetsetnwork

SUBROUTINE jmgeterreur(arret)

  IMPLICIT NONE

  ! Arguments
  LOGICAL, INTENT(out) :: arret

  ! Variables locales

  CALL jmgetseterreur(arret,'g')

END SUBROUTINE jmgeterreur

SUBROUTINE jmerreur1(nomsp,code,var1)

  IMPLICIT NONE

  ! Arguments
  CHARACTER(len=*), INTENT(in) :: nomsp
  INTEGER, INTENT(in) :: code
  INTEGER, INTENT(in) :: var1

  ! Variables locales
  INTEGER :: arret
  CHARACTER(len=80) :: message

  CALL jmgetstop(arret)
  IF (arret == 1) THEN
    CALL jmgetmessage(code,message)
    PRINT *,'JMFFT error in ',TRIM(nomsp),' : ',TRIM(message), &
    & ' (',var1,')'
    STOP 1
  ELSE
    CALL jmsetcode(code)
  END IF

END SUBROUTINE jmerreur1


SUBROUTINE jmerreur2(nomsp,code,var1,var2)

  IMPLICIT NONE

  ! Arguments
  CHARACTER(len=*), INTENT(in) :: nomsp
  INTEGER, INTENT(in) :: code
  INTEGER, INTENT(in) :: var1, var2

  ! Variables locales
  INTEGER :: arret
  CHARACTER(len=80) :: message

  CALL jmgetstop(arret)
  IF (arret == 1) THEN
    CALL jmgetmessage(code,message)
    PRINT *,'JMFFT error in ',TRIM(nomsp),' : ',TRIM(message), &
    & ' (',var1,var2,')'
    STOP 1
  ELSE
    CALL jmsetcode(code)
  END IF

END SUBROUTINE jmerreur2

SUBROUTINE jmerreur3(nomsp,code,var1,var2,var3)

  IMPLICIT NONE

  ! Arguments
  CHARACTER(len=*), INTENT(in) :: nomsp
  INTEGER, INTENT(in) :: code
  INTEGER, INTENT(in) :: var1, var2, var3

  ! Variables locales
  INTEGER :: arret
  CHARACTER(len=80) :: message

  CALL jmgetstop(arret)
  IF (arret == 1) THEN
    CALL jmgetmessage(code,message)
    PRINT *,'JMFFT error in ',TRIM(nomsp),' : ',TRIM(message), &
    & ' (',var1,var2,var3,')'
    STOP 1
  ELSE
    CALL jmsetcode(code)
  END IF

END SUBROUTINE jmerreur3

SUBROUTINE jmerreur4(nomsp,code,var1,var2,var3,var4)

  IMPLICIT NONE

  ! Arguments
  CHARACTER(len=*), INTENT(in) :: nomsp
  INTEGER, INTENT(in) :: code
  INTEGER, INTENT(in) :: var1, var2, var3,var4

  ! Variables locales
  INTEGER :: arret
  CHARACTER(len=80) :: message

  CALL jmgetstop(arret)
  IF (arret == 1) THEN
    CALL jmgetmessage(code,message)
    PRINT *,'JMFFT error in ',TRIM(nomsp),' : ',TRIM(message), &
    & ' (',var1,var2,var3,var4,')'
    STOP 1
  ELSE
    CALL jmsetcode(code)
  END IF

END SUBROUTINE jmerreur4

SUBROUTINE jmgetstop(arret)

  IMPLICIT NONE

  ! Arguments
  INTEGER, INTENT(out) :: arret

  ! Variables locales

  CALL jmgetsetstop(arret,'g')

END SUBROUTINE jmgetstop

SUBROUTINE jmgetsetstop(arret,TYPE)

  ! Subroutine qui permet de stocker une valeur statique
  ! Ceci evite de recourir a un common ...

  IMPLICIT NONE

  ! Arguments
  INTEGER, INTENT(inout) :: arret
  CHARACTER(len=1), INTENT(in) :: TYPE

  ! Variables locales

  ! Variable statique
  INTEGER, SAVE :: arret_last = 1

  IF (TYPE == 's') THEN
    arret_last = arret
  ELSE IF (TYPE == 'g') THEN 
    arret = arret_last
  END IF

END SUBROUTINE jmgetsetstop

SUBROUTINE jmgetmessage(code,message)

  IMPLICIT NONE

  ! Arguments
  INTEGER, INTENT(in) :: code
  CHARACTER(len=*), INTENT(out) :: message

  ! Variables locales
  INTEGER, PARAMETER :: mm = 26
  CHARACTER(len=34), DIMENSION(0:mm-1) :: messages = (/ &
  & "No error                         ",                &
  & "Isign must be equal -1 or 1      ",                &
  & "Isign must be equal 0, -1 or 1   ",                &
  & "Prime numbers to large           ",                &
  & "Size of work is less or equal 0  ",                &
  & "Size of work is to small         ",                &
  & "Impossible to decompose          ",                &
  & "Too many prime factors           ",                &
  & "Nz must be >= 1                  ",                &
  & "ldx doit etre >= n               ",                &
  & "ldx doit etre >= n/2+1           ",                &
  & "ldx1 doit etre >= n              ",                &
  & "ldx1 doit etre >= n/2+1          ",                &
  & "ldx2 doit etre >= m              ",                &
  & "ldy doit etre >= n               ",                &
  & "ldy doit etre >= n+2             ",                &
  & "ldy doit etre >= n/2+1           ",                &
  & "ldy1 doit etre >= n              ",                &
  & "ldy1 doit etre >= n+2            ",                &
  & "ldy1 doit etre >= n/2+1          ",                &
  & "ldy2 doit etre >= m              ",                &
  & "Nz must be >= 1                  ",                &
  & "Ny or Nx must be even            ",                &
  & "Nx must be >= 1                  ",                &
  & "Nx must be even                  ",                &
  & "Nx or Ny or Nz must be even      "                 &
  & /)

  IF (code < 0 .OR. code >= mm) THEN
    PRINT *,'JMFFT GETMESSAGE invalid code : ',code
  END IF

  message = messages(code)

END SUBROUTINE jmgetmessage

SUBROUTINE jmsetcode(code)

  IMPLICIT NONE

  ! Arguments
  INTEGER, INTENT(in) :: code

  ! Variables locales
  INTEGER :: errcode

  errcode = code
  CALL jmgetsetcode(errcode,'s')

END SUBROUTINE jmsetcode

SUBROUTINE jmgetsetcode(code,TYPE)

  ! Subroutine qui permet de stocker le dernier code de retour obtenu
  ! Ceci evite de recourir a un common ...

  IMPLICIT NONE

  ! Arguments
  INTEGER, INTENT(inout) :: code
  CHARACTER(len=1), INTENT(in) :: TYPE

  ! Variables locales

  ! Variable statique
  INTEGER, SAVE :: code_last = 0

  IF (TYPE == 's') THEN
    code_last = code
  ELSE IF (TYPE == 'g') THEN 
    code = code_last
  END IF

END SUBROUTINE jmgetsetcode

SUBROUTINE jmgetseterreur(arret,TYPE)

  ! Subroutine qui permet de stocker une valeur statique
  ! Ceci evite de recourir a un common ...

  IMPLICIT NONE

  ! Arguments
  LOGICAL, INTENT(inout) :: arret
  CHARACTER(len=1), INTENT(in) :: TYPE

  ! Variables locales

  ! Variable statique
  LOGICAL, SAVE :: arret_last = .TRUE.

  IF (TYPE == 's') THEN
    arret_last = arret
  ELSE IF (TYPE == 'g') THEN 
    arret = arret_last
  END IF

END SUBROUTINE jmgetseterreur

END MODULE GFT_JMFFT
