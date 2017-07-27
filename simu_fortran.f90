!Copyright 2017 Bruno Henrique Cervelin
!This file is part of chickenSimulation
!
!chickenSimulation is free software: you can redistribute it and/or modify
!it under the terms of the GNU General Public License as published by
!the Free Software Foundation version 3.
!
!This program is distributed in the hope that it will be useful,
!but WITHOUT ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!GNU General Public License for more details.
!
!You should have received a copy of the GNU General Public License
!along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
PROGRAM simu_frangos
  !module with the simulation internal parameters - modify it to use with the
  !especifities of your problem
  USE ff
  !module with the forces applied over the chickens
  USE forces
  !module with the random variables subroutines
  USE aleatory
  !module with the parameter update subroutines
  USE update_parameters
  IMPLICIT NONE
  !number of chickens in broiler house
  INTEGER, PARAMETER :: n = n_frangos
  INTEGER :: i , seed
  !time discretization, total simulation time, in seconds and number of steps needed to perform the simulation
  REAL, parameter ::  dt = 1.0,  tprod = dias*24.0*3600.0
  integer, parameter :: nsteps = INT( tprod / dt )
  !friction parameter
  real :: lambda = 0.5
  !time counting variables
  REAL :: comeco,fim
  !chicken parameters
  REAL :: fome(n),sede(n),calor(n) = 0,frio(n) = 0,sono(n) = 0,conf = 1,morto(n) = 0
  !matrices with position, velocities and forces
  REAL :: x(n,2), v(n,2) = 0, f(n,2) = 0, f_ant(n,2) = 0,vec_aux(n,2) = 0,f_parede(n,2) = 0
  !environmental temperature on horizontal discretization
  real :: Tamb(npx) = (/(Tin + (Tout-Tin)*(i-1.0)/(npx-1.0), i = 1,npx) /)
  !hour and hour update parameter
  !simulation stats at 0:00 a.m. of day 0
  REAL :: hora = 0.0,up_hora = dt/3600.0
  !day and index used to save chicken posistions
  INTEGER :: dia = 0,ind = 0
  !chicken radius
  REAL :: r = rmin
  !time discretization divided by 20
  real :: dtDiv20 = dt/20.0
  integer :: dorme
  !number of alive chickens, and index of alive chickens
  integer :: n_vivos = n, vivos(n) = (/(i, i = 1,n)/)
  !mede tempo
  CALL CPU_TIME(comeco)
  ! Initialize random number generator
  CALL seed_from_time(seed)
  CALL init_random_number(seed)
  ! Initial coordinates
  CALL initial(n,x)
  CALL update_part(n,x,r,Tamb,vivos,n_vivos)
  !initialize hunger and thirst
  CALL RANDOM_NUMBER(fome)
  CALL RANDOM_NUMBER(sede)
  fome = fome*fmin
  sede = sede*smin
  ! outputs
  !chicken positions
  CALL printx(n,x,ind,vivos,n_vivos)
  !parameters as .m files
  call printMatLabinit(n_vivos,vivos,n,conf, r, sede,fome,frio,calor)
  !print information on screen
  call imprimetela(n,fome,sede,Tamb,frio,calor,sono,conf,r,dia,vivos,n_vivos)
  !star simlation
  DO i = 1, nsteps
    !update hour
    hora = hora + up_hora
    !save parameters every hour
    if (MOD(hora,1.0) .lt. up_hora) call printMatLab(n_vivos,vivos,n,conf, r, sede,fome,frio,calor)
    !verify if chickens are sleeping or not
    dorme = sleep(hora)
    !Updating positions
    IF ((dorme .eq. 0).or. (dorme.eq. 1)) THEN
      x(vivos(1:n_vivos),:) = x(vivos(1:n_vivos),:) + v(vivos(1:n_vivos),:)*dt + &
                              0.5 *f(vivos(1:n_vivos),:)*dt**2
      call image_vector(n,x)
      !update partition information
      CALL update_part(n,x,r,Tamb,vivos,n_vivos)
    END IF
    !update parameters
    CALL update_all(n,x,dt,fome,sede,calor,frio,sono,conf,r,dorme,vivos,n_vivos,morto)
    call vec_vivos(n,morto,n_vivos,vivos)
    !if chickens are awaken
    IF (dorme .eq. 0) THEN
      !save force from last step
      f_ant(vivos(1:n_vivos),:) = f(vivos(1:n_vivos),:)
      !evalute force from parameter
      CALL decide_fome_sede_temp(n,x,fome,sede,calor,frio,f,vivos,n_vivos)
      ! Add friction force
      f(vivos(1:n_vivos),:) = f(vivos(1:n_vivos),:) - lambda *v(vivos(1:n_vivos),:)
      !wall force
      f_parede(vivos(1:n_vivos),:) = force_parede(n_vivos,x(vivos(1:n_vivos),:),v(vivos(1:n_vivos),:))
      !random force
      vec_aux(:,1) = (1-MIN(fome/fmax,1.0 ))*(1-MIN(sede/smax,1.0 ))*&
                     (1-MIN(calor/calormax,1.0 ))*(1-MIN(frio/friomax,1.0 ))
      vec_aux(:,2) = vec_aux(:,1)
      vec_aux = vec_aux*mat_normal(n,2,   0.0 , 0.2*dt)
      ! Updating velocities
      v(vivos(1:n_vivos),:) = v(vivos(1:n_vivos),:) + &
        dt*(f(vivos(1:n_vivos),:)-f_ant(vivos(1:n_vivos),:)+f_parede(vivos(1:n_vivos),:) + &
        vec_aux(vivos(1:n_vivos),:))
    !if chickens are sleeping, or going to sleep
    ELSE
      !wall force and paremeters force is null
      f = 0.0
      f_parede = 0.0
      f_ant = 0.0
      !if chickens are going to sleep, their velocity is decreasing
      IF (dorme .eq. 1) THEN
        v(vivos(1:n_vivos),:) = v(vivos(1:n_vivos),:)/(1+dtDiv20 )
      !if chickens are sleeping, their velocity is 0
      else
        v = 0.0
      end if
    END IF

    !print dailly information on screen
    IF (hora .GE. 24.0 ) THEN
      hora = 0.0
      dia = dia+1
      call imprimetela(n,fome,sede,Tamb,frio,calor,sono,conf,r,dia,vivos,n_vivos)
    END IF
    !print chicken position at 9:00 a.m.
    if (abs(hora - 9.0) .lt. 1e-16) then
      ind = ind+1
      CALL printx(n,x,ind,vivos,n_vivos)
    end if
    !if all chickens are dead, end simulation
    if (n_vivos .eq. 0) goto 321
  END DO
  !close .m files
321 call printMatLabfinal()
  !time required to perform simulation
  CALL CPU_TIME(fim)
  WRITE(*,*) 'time spent (seconds) :: ', fim - comeco

contains
  !verify if chickens are sleeping
  integer function sleep(hora)
    !reurns 0 if they are awaken
    !reurns 1 if they are going to sleep
    !reurns 2 if they are sleeping
    implicit none
    real :: hora
    if (((hora .gt. 6.0) .and. (hora .lt. 22.0)) .or. &
        ((hora .gt. 0.0) .and. (hora .lt. 1.0))  .or. &
        ((hora .gt. 3.0) .and. (hora .lt. 4.0))) then
      sleep = 0
    elseif (((hora .ge. 22.0) .and. (hora .lt. 22.03)) .or. &
            ((hora .ge. 1.0)  .and. (hora .lt. 1.03))  .or. &
            ((hora .ge. 4.0)  .and. (hora .lt. 4.03))) then
      sleep = 1
    else
      sleep = 2
    end if
  end function
  !inital posistion of chickens
  SUBROUTINE initial(n,x)
    USE ff
    USE update_parameters
    USE aleatory
    IMPLICIT NONE
    INTEGER :: n
    REAL :: x(n,2)
    ! Creating random initial coordinates
    CALL RANDOM_NUMBER(x)
    x(:,1) = -sidex/2.0  + sidex*x(:,1)
    x(:,2) = -sidey/2.0  + sidey*x(:,2)
  END SUBROUTINE initial
  !print information on screen
  subroutine imprimetela(n,fome,sede,Tamb,frio,calor,sono,conf,r,dia,vivos,n_vivos)
    use ff
    implicit none
    integer :: n,n_vivos,vivos(n),dia
    real :: fome(n),sede(n),Tamb(npx),frio(n),calor(n),sono(n),conf,r
    WRITE(*,'(" Day= ", I9, " hunger ", f9.2, " thirst ", f9.2, " temp  ", f9.2, " Tidx ",f9.2)') &
            dia, SUM(fome(vivos(1:n_vivos)))/n_vivos, SUM(sede(vivos(1:n_vivos)))/n_vivos, &
            sum(Tamb)/npx, sum(Tindex(vivos(1:n_vivos)))/n_vivos
    WRITE(*,'("                cold   ", f9.2, " hot    ", f9.2, " sleep ", f9.2, " comf ", f9.2)') &
            SUM(frio(vivos(1:n_vivos)))/n_vivos, SUM(calor(vivos(1:n_vivos)))/n_vivos,SUM(sono(vivos(1:n_vivos)))/n_vivos, conf
    WRITE(*,'("                radius ", f9.2, " alive  ", I9)') r, n_vivos
  end subroutine
  !save chicken position in file
  SUBROUTINE printx(n,x,ind,vivo,n_vivo)
    USE ff
    IMPLICIT NONE
    INTEGER :: n, i,ind,n_vivo,j
    REAL :: x(n,2)
    integer :: vivo(n)
    Character(len=11) :: nome
    !find closest pixel to chicken position
    call escrevenome(nome,ind)
    open(unit=76,file=nome)
    do j = 1,n_vivo
      i = vivo(j)
      write(76,*) x(i,1), '   ', x(i,2)
    end do
  END SUBROUTINE printx
  !chicken output position file name
  subroutine escrevenome(nome,ind)
    implicit none
    integer :: ind
    Character(len=11) :: nome
    !comeca rotina
    if (ind < 10) then
      write(nome,'("out000",I1,".out")') ind
    elseif (ind < 100) then
      write(nome,'("out00",I2,".out")') ind
    elseif (ind < 1000) then
      write(nome,'("out0",I3,".out")') ind
    else
      write(nome,'("out",I4,".out")') ind
    end if
  end subroutine
  !write parameters in .m files
  subroutine printMatLab(n_vivos,vivos,n,conf, r, sede,fome,frio,calor)
    implicit none
    integer :: n_vivos,n,vivos(n)
    real :: conf, r, sede(n),fome(n),frio(n),calor(n)
    write(32,*) conf, ";"
    write(33,*) n_vivos, ";"
    write(34,*) r, ";"
    write(35,*) sum(sede(vivos(1:n_vivos)))/n_vivos,';'
    write(37,*) sum(fome(vivos(1:n_vivos)))/n_vivos,';'
    write(38,*) sum(frio(vivos(1:n_vivos)))/n_vivos,';'
    write(39,*) sum(calor(vivos(1:n_vivos)))/n_vivos,';'
  end subroutine
  !open .m files
  subroutine printMatLabinit(n_vivos,vivos,n,conf, r, sede,fome,frio,calor)
    implicit none
    integer :: n_vivos,n,vivos(n)
    real :: conf, r, sede(n),fome(n),frio(n),calor(n)
    !open files
    open(32,file='conf.m')
    open(33,file='vivo.m')
    open(34,file='raio.m')
    open(35,file='sede.m')
    open(37,file='fome.m')
    open(38,file='frio.m')
    open(39,file='calor.m')
    !save files
    write(32,*) 'c =[',conf,';'
    write(33,*) 'v =[',n_vivos,';'
    write(34,*) 'r =[',r,';'
    write(35,*) 's =[',sum(sede)/n,';'
    write(37,*) 'f =[',sum(fome)/n,';'
    write(38,*) 'fr =[',sum(frio)/n,';'
    write(39,*) 'ca =[',sum(calor)/n,';'
  end subroutine
  !close .m files
  subroutine printMatLabfinal()
    implicit none
    write(32,*) '];'
    CLOSE(32)
    write(33,*) '];'
    close(33)
    write(34,*) '];'
    close(34)
    write(35,*) '];'
    close(35)
    write(37,*) '];'
    close(37)
    write(38,*) '];'
    close(38)
    write(39,*) '];'
    close(39)
   end subroutine
END PROGRAM
