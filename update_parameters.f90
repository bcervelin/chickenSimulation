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
MODULE update_parameters
  IMPLICIT NONE
  !thirst parameters
  REAL, PARAMETER :: smin = 1.0/15.0
  REAL, PARAMETER :: smax = 1.0
  REAL, PARAMETER :: DeltaS = smax-smin
  !hunger parameters
  REAL, PARAMETER :: fmin = 1/15.0
  REAL, PARAMETER :: fmax = 1.0
  REAL, PARAMETER :: DeltaF = Smax-Smin
  !temperature parameters
  REAL, PARAMETER :: calormin = 2.0 /15.0
  REAL, PARAMETER :: calormax = 5*calormin
  REAL, PARAMETER :: DeltaCalor = Calormax-Calormin
  REAL, PARAMETER :: friomin = 1.0 /1.0d1
  REAL, PARAMETER :: friomax = 7*friomin
  REAL, PARAMETER :: DeltaFrio = Friomax-Friomin
  REAL, PARAMETER :: Tmax = 23.0
  REAL, PARAMETER :: Tmin = 20.0
  REAL, PARAMETER :: Tchicken = 37.0
  !sleep parameters
  REAL, PARAMETER :: sonomin = 0.8
  !comfort parameters
  REAL, PARAMETER :: confmin = 0.8
  !radius parameters
  REAL, PARAMETER :: rmax = 1.4
  REAL, PARAMETER :: rmin = 0.3
  REAL, PARAMETER :: DeltaR = rmax-rmin
  !increase and decrease parameters
  REAL,PARAMETER :: upcalor   = calormin/1.2d2,        &
                    downcalor = 2*calormin/6.0d1,      &
                    upfrio    = friomin/1.8d2,         &
                    downfrio  = 2*friomin/1.2d2  ,     &
                    dormirfome = 9.0 *fmin/(6*3.6d3 ),&
                    downfome =  5.0 *fmin/1.2d2,       &
                    upfome =   fmin/6.0d2,             &
                    dormirsede = 9.0 *smin/(6*3.6d3 ),&
                    downsede = 6.0 *smin/2.0d1,       &
                    upsede =   smin/1.3d2,             &
                    upsono = sonomin/(19*3600.0),      &
                    downsono = 1.0/(5*3600.0),         &
                    vida = 42*24*3600.0, &
                    upmorto = 1.0/60./60./2. , &
                    downmorto = upmorto*3.
  !maximum force
  REAL, PARAMETER :: force_max = 3.0

  private :: upcalor ,downcalor,upfrio,downfrio,dormirfome, &
             downfome,upfome,dormirsede,downsede,upsede,&
             upsono,downsono,vida



CONTAINS

!funcao calcula o indice TEIbc - disponicel em Math Model for Thermal
!entradas :: m,n       = dimensoes das particoes no espaco
!         :: Temp(m,n) = temperatura sentida em cada particao
!                        (considerando apenas concentracao de frangos)
!         :: RH        = umidade relativa do ar (constante em todo aviario)
!         :: V         = velocidade do ar (constante em todo aviario)
!saida    :: Tindex    = indice TEIbc
real function TEIbc(Temp,RH,V)
  implicit none
  real :: Temp,RH,V
  !comeca rotina
  TEIbc = 45.6026 - 2.3107*Temp - 0.3683*RH + 9.7092*V + 0.05492*Temp**2 + &
           0.00121*RH**2 + 0.66329*V**2 + 0.0128968*Temp*RH - 0.300928*Temp*V - &
           0.05952*RH*V
  !write(*,*) Temp,RH,V,TEIbc
  !normaliza o indice - menor que 27 - conforto
  !                     maior que 35 - desconforto extremo
  TEIbc = max(TEIbc-27.,0.)/3.
end function
!find partition the chicken is
SUBROUTINE find_part(x,px,py)
  USE ff
  IMPLICIT NONE
  INTEGER :: px,py
  REAL :: x(2)
  px = CEILING((x(1)+sidex/2)/szp)
  py = CEILING((x(2)+sidey/2)/szp)
  IF (py .EQ. 0) py = 1
  IF (px .EQ. 0) px = 1
END SUBROUTINE

SUBROUTINE nearest_eater(x,pi,pj)
  USE ff
  IMPLICIT NONE
  REAL :: x(2),pi,pj
  !local
  REAL :: xi,pf,pc
  !comeca rotina
  xi =  x(1) + sxD2
  IF (xi .LE. espc) THEN
    pi = espc-sxD2
  ELSEIF (xi .GE. nc*espc) THEN
    pi = nc*espc-sxD2
  ELSE
    pf = FLOOR(xi/espc)*espc
    pc = CEILING(xi/espc)*espc
    IF (ABS(xi - pf) .GT. ABS(pc-xi)) THEN
      pi = pc - sxD2
    ELSE
      pi = pf - sxD2
    END IF
  END IF
  !encontra posicao y do comedouro mais proximo
  IF (x(2) .GT. 0.0 ) THEN
    pj = 4*syD5 - syD2
  ELSE
    pj = syD5 - syD2
  END IF
END SUBROUTINE

SUBROUTINE nearest_drinker(x,pi,pj)
  USE ff
  IMPLICIT NONE
  REAL :: x(2),pi,pj
  !local
  REAL :: xi,pf,pc
  !comeca rotina
  xi =  x(1) + sxD2
  IF (xi .LE. espb) THEN
    pi = espb-sxD2
  ELSEIF (xi .GE. nb*espb) THEN
    pi = nb*espb-sxD2
  ELSE
    pf = FLOOR(xi/espb)*espb
    pc = CEILING(xi/espb)*espb
    IF (ABS(xi - pf) .GT. ABS(pc-xi)) THEN
      pi = pc - sxD2
    ELSE
      pi = pf - sxD2
    END IF
  END IF
  !encontra posicao y do comedouro mais proximo
  IF (x(2) .GT. 0.0 ) THEN
    pj = 3*syD5 - syD2
  ELSE
    pj = 2*syD5 - syD2
  END IF
END SUBROUTINE

subroutine vec_vivos(n,morto,n_vivos,vivos)
  implicit NONE
  integer :: n, n_vivos
  integer :: vivos(n)
  real :: morto(n)
  !local
  integer :: i,j,count,aux(n)
  !comeca rotina
  count = 0
  do j  = 1,n_vivos
    i = vivos(j)
    if(morto(i) .lt. 1.) then
      count = count +1
      aux(count) = i
    end if
  end do
  n_vivos = count
  vivos(1:n_vivos) = aux(1:n_vivos)
end subroutine


SUBROUTINE update_all(n,x,dt,hora,fome,sede,calor,frio,sono,conf,r,dormindo,vivos,n_vivos,morto)
  USE ff
  IMPLICIT NONE
  INTEGER :: n,dormindo,n_vivos,vivos(n)
  REAL :: x(n,2),dt,hora,fome(n),sede(n),calor(n),frio(n),sono(n),morto(n)
  REAL :: conf,r
  !local
  INTEGER :: i,ii,jj,j
  REAL :: xi,yi,influ,sz, confi(n)
  !comeca rotina
  DO j = 1,n_vivos
      i = vivos(j)
      if (morto(i) .le. 1) then
      IF (dormindo .eq. 0) THEN
        !fome
        CALL nearest_eater(x(i,:),xi,yi)
        sz = SQRT((xi-x(i,1))**2+(yi-x(i,2))**2)
        IF (sz .LE. force_max ) THEN
          fome(i) = MAX(fome(i)-downfome*dt,0.0 )
        ELSE
          fome(i) = MIN(fome(i)+upfome*dt,1.0 )
        END IF
        !sede
        !calor influencia na sede
        CALL nearest_drinker(x(i,:),xi,yi)
        sz = SQRT((xi-x(i,1))**2+(yi-x(i,2))**2)
        IF (sz .LE. force_max) THEN
          sede(i) = MAX(sede(i)-downsede*dt,0.0 )
        ELSE
          sede(i) = MIN(sede(i)+upsede*(1+calor(i))*dt,1.0 )
        END IF
      ELSE
        fome(i) = MIN(fome(i) + dormirfome*dt,1.0 )
        sede(i) = MIN(sede(i) + dormirsede*(1+calor(i))*dt,1.0 )
      END IF
      !termico
      CALL find_part(x(i,:),ii,jj)
      !atualiza o Indice termico
      Tindex(i) = TEIbc(Temp(ii,jj),RH,AS)
      !desconforto termico
      if (Tindex(i) .gt. 0) THEN
        influ = min(Tindex(i),1.0)
        !frio - desconforto eh maior do que se a temperatura fosse maior
        if (Tindex(i) .gt. TEIbc(Temp(ii,jj)+1,RH,AS)) THEN
          IF (calor(i) .LT. calormin) THEN
            frio(i) = MIN(frio(i) + influ*upfrio*dt,1.0 )
          ENDIF
          calor(i) = MAX(calor(i) - (1+influ)*downcalor*dt,0.0 )
        !calor
        ELSE
          IF (frio(i) .LT. friomin) THEN
            calor(i) = MIN(calor(i) + influ*upcalor*dt,1.0 )
          END IF
          frio(i) = MAX(frio(i) - (1+influ)*downfrio*dt,0.0 )
        end if
      !conforto termico
      ELSE
        frio(i) = MAX(frio(i) - downfrio*dt,0.0 )
        calor(i) = MAX(calor(i) - downcalor*dt,0.0 )
      END IF
    end if
  END DO
  !conforto
  confi = (1.0  - MAX(fome-3*fmin,0.0 )/(1.0-3*fmin))* &
          (1.0  - MAX(sede-3*smin,0.0 )/(1.0-3*smin))* &
          (1.0  - MAX(frio-1.5*friomin,0.0 )/(1.0-1.5*friomin))* &
          (1.0  - MAX(calor-1.5*calormin,0.0 )/(1.0-1.5*calormin))
  call update_morto(n,confi,vivos,n_vivos,morto,dt)
  !sono
  IF (dormindo .eq. 0) THEN
    sono = MIN(sono + upsono*dt,1.0 )
  ELSE
    sono = MAX(sono - (1 - (MAX(confmin - confi,0.0 ))/confmin)*downsono*dt,0.0 )
  END IF
  !atualiza conforto com sono
  confi = confi/(1+MAX(sono-sonomin,0.0 )/(1-sonomin))
  conf = SUM(confi(vivos(1:n_vivos)))/n_vivos
  !tamanho
  !atualiza tamanho pela versao exponencial
  call update_radius(r,dt*MIN(confmin,conf)/confmin)
  !r = r + dt*(MIN(confmin,conf)/confmin)*DeltaR/vida
END SUBROUTINE

subroutine update_radius(r,dt)
  implicit NONE
  real :: r,dt
  !local
  integer :: i,j
  real :: fatorial
  real :: k = 1-exp(-1.0)
  real :: soma = 0
  !comeca rotina
  r = r + (DeltaR/K + rmin - r)*(1-exp(-dt/vida))
end subroutine

SUBROUTINE update_morto(n,confi,vivos,n_vivos,morto,dt)
  implicit NONE
  integer:: n, n_vivos,vivos(n)
  real :: confi(n),morto(n),dt
  !local
  integer :: i,j
  !comeca rotina
  do j = 1,n_vivos
    i = vivos(j)
    if (confi(i) .eq. 0) then
      morto(i) = min(morto(i) + upmorto*dt,1.0)
    else
      morto(i) = max(morto(i) - downmorto*dt,0.0)
    end if
  end do
end SUBROUTINE

FUNCTION average_filter(m,n,A) RESULT(FA)
  IMPLICIT NONE
  INTEGER :: m,n,A(m,n)
  REAL FA(m,n)
  !local
  INTEGER :: i,j,ii,jj
  !comeca rotina
  DO i =1,m
    DO j = 1,n
      FA(i,j) = 0.0
      DO ii = MAX(i-1,1),MIN(i+1,m)
        DO jj = MAX(j-1,1),MIN(j+1,m)
          IF ((ii .NE. i) .AND. (jj.NE.j)) THEN
            FA(i,j) = FA(i,j) + A(ii,jj)
          ELSEIF ((ii .NE. i) .OR. (jj.NE.j)) THEN
            FA(i,j) = FA(i,j) + 2*A(ii,jj)
          ELSE
            FA(i,j) = FA(i,j) + 4*A(ii,jj)
          END IF
        END DO
      END DO
    END DO
  END DO
  FA = FA/16.0
END FUNCTION

SUBROUTINE update_part(n,x,r,Tamb,vivos,n_vivos)
  USE ff
  IMPLICIT NONE
  INTEGER :: n,vivos(n),n_vivos
  REAL :: x(n,2),Tamb(npx),r
  REAL, PARAMETER :: pi = 4. *ATAN(1. )
  !local
  INTEGER :: i,j,px,py,np
  !comeca rotina
  num_part = 0
  parti = 0
  DO j = 1,n_vivos
    i = vivos(j)
    CALL find_part(x(i,:),px,py)
    num_part(px,py) = num_part(px,py) + 1
    np = num_part(px,py)
    parti(px,py,np) = i
  END DO
  !calcula a temperatura em cada particao
  Temp = average_filter(npx,npy,num_part)
  do i = 1,npx
    Temp(i,:) = Tamb(i) + (Tchicken-Tamb(i))*Temp(i,:)*(pi*r**2)/(szp**2)
  end do
END SUBROUTINE

END MODULE
