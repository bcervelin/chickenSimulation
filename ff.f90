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
  !number of days to be simulated
  integer, parameter :: dias = 3
  !number of chickens
  INTEGER, PARAMETER :: n_frangos = 1000
  !number of eates and number of water fountains
  INTEGER, PARAMETER :: nc=450!INT(n_frangos/2/45)
  INTEGER, PARAMETER :: nb=2000!INT(n_frangos/2/10)
  !size of the broiler house in dm
  REAL, PARAMETER :: sidex = 1500.  !1pixel = 0.1 m
  REAL, PARAMETER :: sidey = 200.
  !some constants
  real, parameter :: sxD2 = sidex/2.0
  real, parameter :: syD2 = sidey/2.0
  real, parameter :: syD5 = sidey/5.0
  real, parameter :: segundos = 24*3600.0
  !distance between eates
  REAL, PARAMETER :: espc = sidex/(nc+1)
  !distance between water fountains
  REAL, PARAMETER :: espb = sidex/(nb+1)

  !size of partition in dm
  REAL, PARAMETER :: szp = 10.0

  !number of partitions
  INTEGER, PARAMETER :: npx = CEILING(sidex/szp)
  INTEGER, PARAMETER :: npy = CEILING(sidey/szp)
  !number of chickens in each partition
  INTEGER :: num_part(npx,npy)
  !temperature felt in each partition
  REAL :: Temp(npx,npy)
  !Normalized thermal index for each chicken
  real :: Tindex(n_frangos)
  !Temperatures at the entrance and at the exit of the house in graus C
  real, parameter :: Tin = 20, Tout = 22
  !Relative Humidity (%) e air velocity (m/s)
  real :: RH = 50, AS = 1

CONTAINS
  subroutine image_vector(n,x)
  !subroutine that prevents the chicken from exiting the house limits  
    IMPLICIT NONE
    INTEGER :: n
    REAL :: x(n,2)
    !comeca rotina
    x(:,1) = MIN(MAX(x(:,1),-sxD2),sxD2)
    x(:,2) = MIN(MAX(x(:,2),-syD2),syD2)
  END subroutine image_vector
END MODULE ff
