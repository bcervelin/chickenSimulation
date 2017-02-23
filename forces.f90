!Copyright 2017 Bruno Henrique Cervelin
!This file is part of chickenSimulation
!
!chickenSimulation is free software: you can redistribute it and/or modify
!it under the terms of the GNU General Public License as published by
!the Free Software Foundation, either version 3 of the License, or
!(at your option) any later version.
!
!This program is distributed in the hope that it will be useful,
!but WITHOUT ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!GNU General Public License for more details.
!
!You should have received a copy of the GNU General Public License
!along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
MODULE forces
IMPLICIT NONE

PUBLIC :: force_fome, decide_fome_sede_temp,force_temp
PRIVATE :: force_sede, temp_subroutine, maior, menor

CONTAINS

FUNCTION force_fome(x,fome) RESULT(f)
  USE update_parameters
  IMPLICIT NONE
  REAL :: x(2), fome, f(2)
  !local
  REAL cx,cy,norm,influ
  !comeca rotina
  f=0.0
  !IF (fome .GT. fmin) THEN
    CALL nearest_eater(x,cx,cy)
    cx = (cx - x(1))
    cy = (cy - x(2))
    norm = SQRT(cx**2+cy**2)
    IF (norm .GT. force_max) THEN
      cx = force_max*cx/norm
      cy = force_max*cy/norm
      influ = sqrt(fome)
      f(1) = influ*cx
      f(2) = influ*cy
    END IF
  !END IF
END FUNCTION

FUNCTION force_sede(x,sede,calor) RESULT(f)
  USE update_parameters
  IMPLICIT NONE
  INTEGER :: n
  REAL :: x(2), sede, f(2), calor
  !local
  REAL cx,cy,norm
  REAL influ
  !comeca rotina
  f =0.0
  CALL nearest_drinker(x,cx,cy)
  cx = (cx - x(1))
  cy = (cy - x(2))
  norm = SQRT(cx**2+cy**2)
  IF (norm .GT. force_max) THEN
    cx = force_max*cx/norm
    cy = force_max*cy/norm
    influ = SQRT(sede)*(1+min(max(calor-calormin,0.0)/DeltaCalor,1.0))
    f(1) = influ*cx
    f(2) = influ*cy
  END IF
END FUNCTION

function dist_eater(x) result(r)
  use update_parameters
  implicit NONE
  real :: x(2), r
  !local
  real :: cx,cy
  CALL nearest_eater(x,cx,cy)
  r = SQRT((x(1)-cx)**2+(x(2)-cy)**2)
end function

function dist_drinker(x) result(r)
  use update_parameters
  implicit NONE
  real :: x(2), r
  !local
  real :: cx,cy
  CALL nearest_drinker(x,cx,cy)
  r = SQRT((x(1)-cx)**2+(x(2)-cy)**2)
end function

SUBROUTINE decide_fome_sede_temp(n,x,fome,sede,calor,frio,f,vivos,n_vivos)
  USE ff
  use update_parameters
  IMPLICIT NONE
  INTEGER :: n,n_vivos,vivos(n)
  REAL :: x(n,2), fome(n),sede(n),f(n,2),calor(n),frio(n)
  !local
  INTEGER :: i,j
  REAL ::maximo
  !comeca rotina
  DO j = 1,n_vivos
    i = vivos(j)
    !se frango esta bebendo fica comendo até matar a sede.
    if ((sede(i) .gt. smin/2.0) .and. (dist_drinker(x(i,:)) .lt. force_max)) THEN
      f(i,:) = 0.0
    else
      maximo = MAX(fome(i),sede(i),calor(i),frio(i))
      IF (maximo .EQ. fome(i)) THEN
        f(i,:) = force_fome(x(i,:),fome(i))
      ELSEIF (maximo .EQ. sede(i)) THEN
        f(i,:) = force_sede(x(i,:),sede(i),calor(i))
      ELSE
        f(i,:) = force_temp(x(i,:),calor(i),frio(i),Tindex(i))
      END IF
    ENDIF
  END DO
END SUBROUTINE

FUNCTION force_temp(x,calor,frio,influ) RESULT(f)
  USE ff
  USE update_parameters
  IMPLICIT NONE
  REAL :: x(2), f(2),fi(2),calor,frio
  !local
  INTEGER :: px,py, ki,kj,ii,jj
  INTEGER :: poss(8), n_poss,aux
  REAL :: influ
  !comeca rotina
  f = 0.0
  IF (frio .GT. friomin) THEN
    CALL find_part(x,px,py)
    CALL temp_subroutine(x,px,py,maior,fi)
    f = influ*fi
  ELSEIF (calor .GT. calormin) THEN
    CALL find_part(x,px,py)
    CALL temp_subroutine(x,px,py,menor,fi)
    f = influ*fi
  END IF
END FUNCTION

SUBROUTINE temp_subroutine(x,px,py,compara,f)
  USE ff
  use update_parameters
  USE aleatory, ONLY : rand_int
  IMPLICIT NONE
  INTEGER :: px,py
  LOGICAL :: compara
  REAL x(2),f(2)
  !local
  INTEGER :: n_poss,ind,i,j
  REAL :: r,poss(8,2)
  !comeca rotina
  f = 0.0
  !verifica possibilidade para se mover
  n_poss = 0
  poss = 0.0
  DO i = MAX(px-1,1),MIN(px+1,npx)
    DO j = MAX(py-1,1),MIN(py+1,npy)
      IF (compara(Temp(i,j),Temp(px,py))) THEN
        n_poss = n_poss + 1
        poss(n_poss,1) = i
        poss(n_poss,2) = j
      END IF
    END DO
  END DO
  !casa contrario frango se move na direcao da particao
  IF (n_poss > 0) THEN
    !escolhe aleatorio das entre as possibilidades
    ind = rand_int(1,n_poss)
    i = poss(ind,1)
    j = poss(ind,2)
    f(1) = (i*szp-0.5*szp - sxD2) - x(1)
    f(2) = (j*szp-0.5*szp - syD2) - x(2)
    r = SQRT(f(1)**2+f(2)**2)
    IF (r .GT. force_max) THEN
      f = force_max*f/r
    END IF
  END IF
END SUBROUTINE

LOGICAL FUNCTION maior(a,b)
  IMPLICIT NONE
  REAL :: a,b
  maior = (a .GT. b)
END FUNCTION

LOGICAL FUNCTION menor(a,b)
  IMPLICIT NONE
  REAL :: a,b
  menor = (a .LT. b)
END FUNCTION

FUNCTION force_parede(n,x,v) RESULT(f)
  USE ff
  use update_parameters
  IMPLICIT NONE
  INTEGER :: n
  REAL :: x(n,2),v(n,2),f(n,2)
  !local
  INTEGER :: i,j
  REAL :: tol = 1.0 ,rx,ry
  !comeca rotina
  f = 0.0
  DO i = 1,n
    rx = (ABS(x(i,1) +sxD2))
    ry = (ABS(x(i,2) +syD2))
    IF (rx .LT. tol) THEN
      if (v(i,1) .LT. 0.0 ) f(i,1) = -v(i,1)
      f(i,1) = f(i,1) + force_max*sidex/800.0
    ENDIF
    IF (ry .LT. tol) THEN
      if (v(i,2) .LT. 0.0 ) f(i,2) = -v(i,2)
      f(i,2) = f(i,2) + force_max*sidey/200.0
    ENDIF
    rx = (ABS(sxD2 - x(i,1)))
    ry = (ABS(syD2 - x(i,2)))
    IF (rx .LT. tol) THEN
      if (v(i,1) .GT. 0.0 ) f(i,1) = -v(i,1)
      f(i,1) = f(i,1) - force_max*sidex/800.0
    ENDIF
    IF (ry .LT. tol) THEN
      if (v(i,2) .GT. 0.0 ) f(i,2) = -v(i,2)
      f(i,2) = f(i,2) - force_max*sidey/200.0
    ENDIF
  END DO
END FUNCTION

END MODULE
