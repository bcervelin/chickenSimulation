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
MODULE ff
  IMPLICIT NONE
  !dias de simulacao
  integer, parameter :: dias = 3
  !numero de frangos
  INTEGER, PARAMETER :: n_frangos = 1000
  !nuymero de comedouros e beberos
  INTEGER, PARAMETER :: nc=450!INT(n_frangos/2/45)
  INTEGER, PARAMETER :: nb=2000!INT(n_frangos/2/10)
  !tamanho do aviario
  REAL, PARAMETER :: sidex = 1500.  !1pixel = 0.1 m
  REAL, PARAMETER :: sidey = 200.
  !constantes usadas repetidamente
  real :: sxD2 = sidex/2.0
  real :: syD2 = sidey/2.0
  real :: syD5 = sidey/5.0
  real :: segundos = 24*3600.0
  !espacao entre comedouros
  REAL, PARAMETER :: espc = sidex/(nc+1)
  !espacao entre beberoudos
  REAL, PARAMETER :: espb = sidex/(nb+1)

  !tamanho das particoes
  REAL, PARAMETER :: szp = 10.0

  !numero de particoes
  INTEGER, PARAMETER :: npx = CEILING(sidex/szp)
  INTEGER, PARAMETER :: npy = CEILING(sidey/szp)
  !numero de francos em cada particao
  INTEGER :: num_part(npx,npy)
  !temperatura em cada particao
  REAL :: Temp(npx,npy)
  !indice de temperatura para cada frango
  real :: Tindex(n_frangos)
  !temperatura na entrada e na saida de ar (ºC)
  real, parameter :: Tin = 20, Tout = 22
  !umidade relativa (%) e velocidade do ar (m/s)
  real :: RH = 50, AS = 1

CONTAINS
  subroutine image_vector(n,x)
    IMPLICIT NONE
    INTEGER :: n
    REAL :: x(n,2)
    !comeca rotina
    x(:,1) = MIN(MAX(x(:,1),-sxD2),sxD2)
    x(:,2) = MIN(MAX(x(:,2),-syD2),syD2)
  END subroutine image_vector
END MODULE ff
