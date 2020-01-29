module simulacao
implicit none
!temporary variable used to scale and simplify simulation (to be removed)
  real :: escala = 500.0
!experiment variables
  !broiler size (0.1m)
  REAL :: sidex, sidey
  !number of chickens, eaters and water fountains
  INTEGER :: nf, nE, nW
  !Temperatura at air input and output (ºC) Relative humidity (%) and air speed (ºC)
  real :: Tin, Tout, RH, AS
!some global variables
  !partition size
  REAL :: szp
  !number of partition
  INTEGER :: npx, npy
  !number of chickens in each partition
  INTEGER, allocatable :: nPart(:,:)
  !temperatura in each partition
  REAL, allocatable :: Temp(:,:)
  !emviromental temperature
  real, allocatable :: Tamb(:)
!some reused values
  real :: syD2,  syD5
  REAL :: espc, espb
  real :: pi, piDszp2, umme
  real :: vida
  REAL :: Deltaho !will be initialized later
!change the variables values according to your data
  !Duration of simulation in days
  integer, parameter :: dias = 42
!chicken radius parameters
  REAL, PARAMETER :: rmax = 1.4, rmin = 0.3, DeltaR = rmax-rmin
!temperature parameters (ºC)
  REAL, PARAMETER :: Tmax = 23.0, Tmin = 20.0, Tchicken = 37.0
!internal parameters used as variables for the optimization problem
  REAL :: homin, homax,  comin,  comax,  thmin, thmax,   humin, humax,  &
          slmin, uho,    dho,    uco,    dco,   slhu,    dhu,   uhu,    &
          slth,  usl,    dth,    uth,    dsl,   udd,     ddd,   fmax,   &
          cfmin, lambda, fwallx, fwally, sigma, lambdasl
contains

  !function used as inital input in optmization problem
  function initInput() result(x)
    real :: x(30)
    x = (/2./15.,   1./3.,    0.1,      0.7,     1./15.,   1.,       1./15.,   &
          1.,       0.8,      1./9.e3,  1./255., 1./1.8e3, 1/6.e2,   1./3.6e4, &
          1./90.,   1./9.e2,  1./3.6e4, 0.02,    1./1950., 1/8.55e4, 1./1.8e4, &
          1./7.2e3, 1./2.4e4, 3.,       0.8,     0.5,      1./8.e2,  1./2.e2,  &
          0.2,      2./2.1 /)
  end function

  !initiate global variables paramters (may be changed)
  subroutine initGlobal()
    sidex = 1500./sqrt(escala); sidey = 200./sqrt(escala)
    nf = ceiling(20000/escala); nE=ceiling(450/escala); nW=ceiling(2000/escala)
    Tin = 20.0; Tout = 22.0; RH = 50.0; AS = 1.0
    szp = min(10.0,sidex,sidey)
    npx = CEILING(sidex/szp); npy = CEILING(sidey/szp)
    allocate(nPart(npx,npy),Temp(npx,npy),Tamb(npx))
    syD2 = sidey/2.0;  syD5 = sidey/5.0
    espc = sidex/(nE+1); espb = sidex/(nW+1)
    pi =4.0*ATAN(1.0); piDszp2=pi/(szp**2); umme=1.0-exp(-1.0)
  end subroutine

  !dealocate global variables
  subroutine dealocGlobal()
    deallocate(nPart,Temp,Tamb)
  endsubroutine

  !test function used as objective function in optimization process
  function func(x) result(f)
    real :: x(30), f
    integer, parameter :: ndt = 3
    real :: fvec(ndt), data(ndt)
    !data with original produced mass
    data = (/2.48,0.0,0.0/)
    call initGlobal()
    call initParam(x)
    !run the simulation for optimal conditions
    Tin = 20.0; Tout = 22.0; RH = 50.0; AS = 1.0
    fvec(1) =  simula()
    !run the simulation for cold conditions
    Tin = 0.0; Tout = 2.0; RH = 50.0; AS = 1.0
    fvec(2) =  simula()
    !run the simulation for optimal conditions
    Tin = 30.0; Tout = 40.0; RH = 50.0; AS = 1.0
    fvec(3) =  simula()
    !f is sum of the squarte diference between the predicted average mass and
    !the experimental one
    write(*,*) fvec
    write(*,*) data
    f = sum((fvec-data)**2)
    call dealocGlobal()
  endfunction

  !allocate internal parameters as the vector x of the optimization problem
  subroutine initParam(x)
    real, intent(in) :: x(30)
    homin   = x(1);   homax   = x(2);   comin   = x(3);   comax   = x(4)
    thmin   = x(5);   thmax   = x(6);   humin   = x(7);   humax   = x(8)
    slmin   = x(9);   uho     = x(10);  dho     = x(11);  uco     = x(12)
    dco     = x(13);  slhu    = x(14);  dhu     = x(15);  uhu     = x(16)
    slth    = x(17);  dth     = x(18);  uth     = x(19);  usl     = x(20)
    dsl     = x(21);  udd     = x(22);  ddd     = x(23);  fmax    = x(24)
    cfmin   = x(25);  lambda  = x(26);  fwallx  = x(27);  fwally  = x(28)
    sigma   = x(29);  lambdasl= x(30)
    Deltaho = homax - homin
  end subroutine

  !simulation returns the mass produced by the simulation
  function simula() result(massa)
    real :: massa
    integer :: n
    !internal simulation parameters
    real, parameter :: dt = 1.0 !time discretitzation parameter (sec)
    real, parameter :: dt2 = dt**2 !save some evaluations
    integer, parameter :: nsteps = int(24.*3600.*dias/dt) !number of steps
    !chicken parameters
    REAL, allocatable :: hu(:), th(:), x1(:), x2(:), tind(:), Dind(:), Tc(:)
    real, allocatable :: ho(:), co(:), sl(:), dd(:), v1(:), v2(:), f1(:), f2(:),&
                        fl1(:),fl2(:),fw1(:),fw2(:),fr1(:),fr2(:)
    logical, allocatable :: ind(:)
    !auxiliary parameters with postion and distance of nearest eater and
    !water fountain
    real, allocatable :: neE1(:),neE2(:),diE(:),neW1(:),neW2(:),diW(:)
    real :: cf=1.0, r
    integer, allocatable:: px(:),py(:)
    !simulation stats at 8:00 a.m. of day 0
    REAL :: hr = 8.0, uhr = dt/3600.0
    !local variables
    integer :: i,dorme
    !initializa some variables
    n = nf
    call initLoc(nf, hu,  th,  x1, x2, tind, Dind, Tc, ho, co, sl,  dd,&
                  v1, v2,  f1,  f2, fl1,fl2,  fw1, fw2, fr1,fr2,neE1,neE2,&
                  diE,neW1,neW2,diW,px,py,ind)
    r = rmin
    vida=dt/(42.*24.*3600.)
    Tamb=(/(Tin + (Tout-Tin)*(i-1.0)/(npx-1.0), i = 1,npx) /)
    call random_number(hu); call random_number(th)
    hu = hu*humin;          th = th*thmin
    call random_number(x1); call random_number(x2)
    x1 = x1*sidex;          x2 = x2*sidey
    !set the temperatura in each area partition
    !(considering the amount of chickens)
    call upPart(n,x1(1:n),x2(1:n),r,px(1:n),py(1:n),tc(1:n))
    !start the iterative process
    !write(*,*) "T_amb  T_chkn  T_ind  D_ind  hunge  thirs  hotne  coldn  "//&
    !            "conf   alive   PMax  PX    PY   PMed   radiu  front"
    do i = 1,nsteps
      dorme = sleep(hr) !verify if chickens are sleeping
      if (dorme.ne.2) then !if chickens ar not sleeping, move then
        call upPos(n,   x1(1:n),x2(1:n),v1(1:n),v2(1:n),f1(1:n),f2(1:n),dt, &
                   dt2, sidex,  sidey)
        !update partition and temperature informations
        call upPart(n,x1(1:n),x2(1:n),r,px(1:n),py(1:n),tc(1:n))
        !update nearest eater position
        ind(1:n) = x2(1:n) .gt. syD2
        call nearest(x1(1:n),x2(1:n),nE,espc,ind(1:n),4*syD5,syD5,neE1(1:n),&
                     neE2(1:n),diE(1:n))
        !update nearest water fountain position
        call nearest(x1(1:n),x2(1:n),nW,espb,ind(1:n),3*syD5,2*syD5,neE1(1:n),&
                     neE2(1:n),diE(1:n))
        !update temperature index
        call TEIbc(tc(1:n),RH,AS,tind(1:n),dind(1:n))
      end if
      !update chickens internal parameters
      call upParameters(n,   x1(1:n),   x2(1:n),   v1(1:n),   v2(1:n),  &
                  f1(1:n),   f2(1:n),   hu(1:n),   th(1:n),   ho(1:n),  &
                  co(1:n),   sl(1:n),   dd(1:n),   tc(1:n),   cf,       &
                  neE1(1:n), neE2(1:n), diE(1:n),  neW1(1:n), neW2(1:n),&
                  diW(1:n),  tind(1:n), dind(1:n), r,     dt, dorme)
      !if all chickens are dead, finish simulation
      if (n.eq. 0) goto 321
      !update velocity and force variables
      if (dorme.eq.0) then !check if the chickens are awaken
        !save force from last step
        fl1(1:n) = f1(1:n);fl2(1:n) = f2(1:n)
        !evalute force from hunger, thirst, hot or cold
        call fParamer(n,         x1(1:n),   x2(1:n),  f1(1:n),   f2(1:n),  &
                      hu(1:n),   th(1:n),   ho(1:n),  co(1:n),   tind(1:n),&
                      neE1(1:n), neE2(1:n), diE(1:n), neW1(1:n), neW2(1:n),&
                      diW(1:n),  px(1:n),   py(1:n))
        !add friction force
        f1(1:n)=f1(1:n)-lambda*v1(1:n);f2(1:n)=f2(1:n)-lambda*v2(1:n)
        !evaluate the wall force
        call fWall(n,x1(1:n),x2(1:n),v1(1:n),v2(1:n),fw1(1:n),fw2(1:n),dt)
        !evaluate the random force
        call fRand(n,fr1(1:n),fr2(1:n),sigma,dt,hu(1:n)/humax,th(1:n)/thmax,&
                   ho(1:n)/homax,co(1:n)/comax)
        !update velocities
        v1(1:n)=v1(1:n)+dt*(f1(1:n)-fl1(1:n)+fw1(1:n)+fr1(1:n))
        v2(1:n)=v2(1:n)+dt*(f2(1:n)-fl2(1:n)+fw2(1:n)+fr2(1:n))
      else !if the chickens are sleeping or goint to
        !set forces as null
        f1(1:n) = 0.0; fw1(1:n)=0.0;fl1(1:n)=0.0
        f2(1:n) = 0.0; fw2(1:n)=0.0;fl2(1:n)=0.0
        !if chickens are going to sleep, apply friction force
        if (dorme.eq.1) then
          v1(1:n)=v1(1:n)*lambdasl*dt; v2(1:n)=v2(1:n)*lambdasl*dt
        else
          v1(1:n) = 0.0; v2(1:n) = 0.0
        endif
      endif
      hr = hr + uhr !update hour parameter
      if (hr.ge.24.0) hr = 0.0!update day
    end do
321 massa = (r**3)*n/nf
    !dealocate internal parameters
    call dealoc_loc(hu,  th,  x1, x2, tind, Dind, Tc, ho, co, sl,  dd,&
                    v1, v2,  f1,  f2, fl1,fl2,  fw1, fw2, fr1,fr2,neE1,neE2,&
                    diE,neW1,neW2,diW,px,py,ind)
  end function

  !update chicken position
  subroutine upPos(n,x1,x2,v1,v2,f1,f2,dt,dt2,sidex,sidey)
    integer :: n
    real :: x1(n),x2(n),v1(n),v2(n),f1(n),f2(n),dt,dt2,sidex,sidey
    x1=x1+v1*dt+0.5*f1*dt2; x1=MIN(MAX(x1,0.0),sidex)
    x2=x2+v2*dt+0.5*f2*dt2; x2=MIN(MAX(x2,0.0),sidey)
  end subroutine

  !evaluate random force parameter
  subroutine fRand(n,f1,f2,sigma,dt,hu,th,ho,co)
    integer :: n
    real :: f1(n),f2(n),sigma,dt,hu(n),th(n),ho(n),co(n)
    !aleatory force
    f2=(1-MIN(hu,1.0))*(1-MIN(th,1.0))*(1-MIN(ho,1.0))*(1-MIN(co,1.0))
    f1=f2*normalVec(n,0.0,sigma*dt)
    f2=f2*normalVec(n,0.0,sigma*dt)
  end subroutine

  !random normal vector
  fUNCTION normalVec(n,mean,sigma) RESULT(v)
    IMPLICIT NONE
    INTEGER :: n
    REAL :: mean,sigma,v(n)
    !local
    REAL :: u1,u2,mp3s,mm3s
    integer :: i
    !comeca rotina
    mp3s = mean+3.0*sigma; mm3s = mean-3.0*sigma
    do i = 1,n
  11  CALL RANDOM_NUMBER(u1); CALL RANDOM_NUMBER(u2)
      v(i) = sigma*SQRT(-2. *LOG(u1))*COS(2. *pi*u2) + mean
      if ((v(i) .gt. mp3s) .or. (v(i) .lt. mm3s)) goto 11
    end do
  END FUNCTION

  !random integer
  INTEGER FUNCTION  randInt(i,f)
    IMPLICIT NONE
    INTEGER :: i,f
    REAL :: u
    CALL RANDOM_NUMBER(u)
    randInt = i + FLOOR((f+1-i)*u)  ! We want to choose one from m-n+1 integers
  END FUNCTION

  !evaluate force parameter
  subroutine fParamer(n,   x1,  x2, f1,  f2,  hu,  th, ho,   co, tind,&
                      neE1,neE2,diE,neW1,neW2,diW, px, py)
    integer :: n
    integer :: px(n), py(n)
    real :: x1(n),  x2(n),  f1(n),  f2(n), hu(n),  th(n),  ho(n), co(n),  &
            tind(n),neE1(n),neE2(n),diE(n),neW1(n),neW2(n),diW(n)
    !auxiliary variables
    integer :: i
    real :: M(n)
    !initiate force with 0
    f1 = 0.0; f2 = 0.0
    !check if chicken is not drinking water
    where((th .le. thmin/2.0) .or. (diW .ge. fmax))
      M = max(hu,th,ho,co) !check the biggest parameter
      where (igualv(M,ho) .and. M.gt.co)
        where(M.le.hu .and. diE.gt.fmax)!use hunger force
          f2=sqrt(hu)*fmax/diE;
          f1=f2*neE1;f2=f2*neE2;
        elsewhere(igualv(M,th) .and. diW.gt.fmax)!use thirst parameter
          f2=SQRT(th)*(1+min(max(ho-homin,0.0)/deltaho,1.0))*fmax/diW
          f1=f2*neW1;f2=f2*neW2;
        endwhere
      endwhere
    endwhere
    !use temperature force
    do i = 1,n
      !check if chicken is feeling enough hot
      if(igual(M(i),ho(i)) .and. ho(i).gt. homin ) then
        call fTemp(x1(i),x2(i),f1(i),f2(i),tind(i),px(i),py(i),menor)
      !check if chicken is feeling enough cold
      elseif(igual(M(i),co(i)) .and. co(i).gt. comin ) then
        call fTemp(x1(i),x2(i),f1(i),f2(i),tind(i),px(i),py(i),maior)
      end if
    end do
  endsubroutine

  !check if two float vector components are equal
  function igualv(a,b) result(sim)
    implicit none; real :: a(:),b(:); logical :: sim(size(a));
    sim = abs(a-b).lt. 1e-7
  end function

  !check if two floats are equal
  logical function igual(a,b)
    implicit none; real :: a,b; igual = abs(a-b).lt. 1e-7
  end function

  !check if two a float is greater than another float
  LOGICAL FUNCTION maior(a,b)
    IMPLICIT NONE; REAL :: a,b; maior = (a .GT. b)
  END FUNCTION

  !check if two a float is lesser than another float
  LOGICAL FUNCTION menor(a,b)
    IMPLICIT NONE; REAL :: a,b; menor = (a .LT. b)
  END FUNCTION

  !evaluate temperature force
  subroutine fTemp(x1,x2,f1,f2,tind,px,py,compare)
    real :: x1,x2,f1,f2,tind
    integer :: px,py
    logical :: compare
    !auxiliary variables
    integer :: i,j,nP,ind, Pos1(8),Pos2(8)
    real :: r
    !initiate force with 0
    f1=0.0;f2=0.0
    !intialize number of possibilities with 0
    np = 0; Pos1 = 0; Pos2 = 0
    !search for colder partitions
    DO i = MAX(px-1,1),MIN(px+1,npx); DO j = MAX(py-1,1),MIN(py+1,npy)
      IF (compare(Temp(i,j),Temp(px,py))) THEN
        nP = nP + 1; Pos1(nP) = i; Pos2(nP) = j
      END IF
    END DO; END DO
    !if there are colder/hotter partitions, move to one of them
    if (nP .gt. 0) then
      ind = randInt(1,nP); i = Pos1(ind); j = Pos2(ind)
      f1 = (i-0.5)*szp - x1; f2 = (j-0.5)*szp - x2;
      r = sqrt(f1**2+f2**2)
      if (r .gt. fmax) then
        f1 = fmax*f1/r; f2 = fmax*f2/r
      end if
      f1 = f1*tind; f2 = f2*tind
    end if
  end subroutine

  !evaluate wall force
  subroutine fWall(n,x1,x2,v1,v2,fw1,fw2,dt)
    integer :: n
    real :: dt,x1(n),x2(n),v1(n),v2(n),fw1(n),fw2(n)
    !initiate force with 0
    fw1=0.0;fw2=0.0
    !apply the force in x-axis
    !check if chicken is close to the left wall
    where(x1 .lt. 1.0)
      !if chicken is moving towards the wall, apply an oposite direction force
      where(v1.lt.0.0) fw1=-v1
      fw1 = fw1 + fmax*sidex*fwallx*dt
    !check if chickens are close to the right wall
    elsewhere(x1.gt. sidex-1.0)
      !if chicken is moving towards the wall, apply an oposite direction force
      where(v1.gt.0.0) fw1=-v1
      fw1 = fw1 - fmax*sidex*fwallx*dt
    endwhere
    !apply the force in y-axis
    !check if chicken is close to the botton wall
    where(x2 .lt. 1.0)
      !if chicken is moving towards the wall, apply an oposite direction force
      where(v2.lt.0.0) fw2=-v2
      fw2 = fw2 + fmax*sidey*fwally*dt
    !check if chickens are close to the upper wall
    elsewhere(x2.gt. sidey-1.0)
      !if chicken is moving towards the wall, apply an oposite direction force
      where(v2.gt.0.0) fw2=-v2
      fw2 = fw2 - fmax*sidey*fwally*dt
    endwhere
  end subroutine

  !update partition information
  subroutine upPart(n,x1,x2,r,px,py,tc)
    integer :: n
    real, intent(in)     :: x1(n), x2(n), r
    real, intent(out)    :: tc(n)
    integer, intent(out) :: px(n),py(n)
    !auxiliary variables
    integer :: i
    !starts routine
    nPart = 0
    !find partition of chicken
    px = max(CEILING(x1/szp),1); py = max(CEILING(x2/szp),1)
    !count the number of chickens in each partition
    do i = 1,n; nPart(px(i),py(i)) = nPart(px(i),py(i)) + 1; end do
    !uses gaussian filter to find the temperature in each partition
    Temp = gaussFilter2(npx,npy,nPart*piDszp2*r**2)
    !temp = npart*pidszp2*r**2
    do i = 1,npx; Temp(i,:) = Tamb(i) + (Tchicken-Tamb(i))*Temp(i,:); end do
    !find the felt temperature for each chicken
    tc = (/(temp(px(i),py(i)),i=1,n)/)
  end subroutine upPart

  !update chicken's internal parameters
  subroutine upParameters(n,     x1,    x2,   v1,    v2,    f1,   f2,        &
                          hu,    th,    ho,   co,    sl,    dd,   tc,   cf,  &
                          neE1,  neE2,  diE,  neW1,  neW2,  diW,  tind, dind,&
                          r,     dt,    dorme)
    integer :: n,dorme
    real :: x1(n),  x2(n),  v1(n), v2(n),  f1(n),  f2(n),&
            hu(n),  th(n),  ho(n), co(n),  sl(n),  dd(n), tc(n), &
            neE1(n),neE2(n),diE(n),neW1(n),neW2(n),diW(n), &
            tind(n),dind(n)
    real :: cf,r,dt
    !auxiliary variables
    integer::nn
    real :: cfi(n)
    logical :: alive(n)
    !if chickens are waken, check if they are eating or drinking
    if (dorme.eq. 0) THEN
      !update hunger parameter
      where(diE.le.fmax);hu=max(hu-dhu*dt,0.0);
      ELSEwhere; hu=min(hu+uhu*dt,1.0); ENDwhere
      !update thist parameter
      where(diW.le.fmax);th=max(th-dth*dt,0.0);
      ELSEwhere; th=min(th+uth*dt*(1.0+ho),1.0); ENDwhere
    else !if chickens are sleeping, increase the hunger and thist parameter
      hu = MIN(hu + slhu*dt,1.0 ); th = MIN(th + slth*(1.0+ho)*dt,1.0 )
    end if
    !check for thermal discomfort
    where(tind.gt.0.0)
      !if derivative is postive => thermal index increase with temperature
      !Dind<0 => cold problem
      where(Dind.lt.0.0)
        !if chicken isn't feeling hot, cold is increased
        where(ho.lt.homin) co=min(co+tind*uco*dt,1.0)
        !chicken stops feeling hot, since it is in a cold region
        where(ho.gt.0.0) ho=max(ho-(1+tind)*dho*dt,0.0)
      !Dind>0 => hot problem
      elsewhere(Dind.gt.0.0)
        !if chicken isn't feeling cold, hot is increased
        where(co.lt.comin) ho=min(ho+tind*uho*dt,1.0)
        !chicken stops feeling cold, since it is in a hot region
        where(co.gt.0.0) co=max(co-(1+tind)*dco*dt,0.0)
      endwhere
    !if chicken is in a comfortable area, cold and hot decreases
    elsewhere
      co = MAX(co-dco*dt,0.0); ho = MAX(ho-dho*dt,0.0)
    endwhere
    !evaluate individual chicken comfort
    cfi = (1.0-MAX(hu-3.0*humin,0.0)/(1.0-3.0*humin))* &
          (1.0-MAX(th-3.0*thmin,0.0)/(1.0-3.0*thmin))* &
          (1.0-MAX(co-1.5*comin,0.0)/(1.0-1.5*comin))* &
          (1.0-MAX(ho-1.5*homin,0.0)/(1.0-1.5*homin))
    !update sleep parameter
    IF (dorme .eq. 0) THEN
      sl = MIN(sl + usl*dt,1.0 )
    ELSE
      sl = MAX(sl-(1-(MAX(cfmin-cfi,0.0))/cfmin)*dsl*dt,0.0 )
    END IF
    !update comfort parameter using sleep parameter
    cfi = cfi/(1+MAX(sl-slmin,0.0 )/(1-slmin))
    !evaluate avertare comfort
    cf =  sum(cfi)/n
    !update chicken size
    r = r + (DeltaR/umme + rmin - r)*(1-exp(-vida))
    !check for dead chickens
    where(cfi.lt.1.0e-7);dd=min(dd+udd*dt,1.0);
    elsewhere;dd=max(dd-ddd*dt,0.0);endwhere
    alive=dd.lt.1.0; nn = count(alive)
    !if the number of alive chicken has decrease, remove dead ones
    if ((nn .gt. 0).and.(nn.ne.n)) then
      x1(1:nn)   = pack(x1(1:n),alive);   x2(1:nn)   = pack(x2(1:n),alive)
      v1(1:nn)   = pack(v1(1:n),alive);   v2(1:nn)   = pack(v2(1:n),alive)
      f1(1:nn)   = pack(f1(1:n),alive);   f2(1:nn)   = pack(f2(1:n),alive)
      hu(1:nn)   = pack(hu(1:n),alive);   th(1:nn)   = pack(th(1:n),alive)
      ho(1:nn)   = pack(ho(1:n),alive);   co(1:nn)   = pack(co(1:n),alive)
      sl(1:nn)   = pack(sl(1:n),alive);   dd(1:nn)   = pack(dd(1:n),alive)
      tc(1:nn)   = pack(tc(1:n),alive);   tind(1:nn) = pack(tind(1:n),alive)
      dind(1:nn) = pack(dind(1:n),alive)
      neE1(1:nn) = pack(neE1(1:n),alive); neE2(1:nn) = pack(neE2(1:n),alive)
      diE(1:nn)  = pack(diE(1:n),alive)
      neW1(1:nn) = pack(neW1(1:n),alive); neW2(1:nn) = pack(neW2(1:n),alive)
      diW(1:nn)  = pack(diW(1:n),alive)
    end if
    if (nn.ne.n) n = nn
  end subroutine upParameters

  !apply a 2-D gaussian filter
  FUNCTION gaussFilter2(m,n,A) RESULT(FA)
    IMPLICIT NONE
    INTEGER :: m,n
    REAL A(m,n),FA(m,n)
    !local
    INTEGER :: i,j,ii,jj
    !comeca rotina
    FA = 0.0
    DO i =1,m; DO j = 1,n
      DO ii = MAX(i-1,1),MIN(i+1,m); DO jj = MAX(j-1,1),MIN(j+1,n)
        IF (ii .NE. i .AND. jj.NE.j) THEN
          FA(i,j) = FA(i,j) + A(ii,jj)
        ELSEIF (ii .NE. i .OR. jj.NE.j) THEN
          FA(i,j) = FA(i,j) + 2.0*A(ii,jj)
        ELSE
          FA(i,j) = FA(i,j) + 4.0*A(ii,jj)
        END IF
      END DO; end do
    END DO; END DO
    FA = FA/16.0
  END FUNCTION

  !verify if chickens are sleeping
  !reurns 0 if they are awaken
  !reurns 1 if they are going to sleep
  !reurns 2 if they are sleeping
  integer function sleep(hr)
    implicit none
    real :: hr
    if (inInt(hr,6.,22.) .or.inInt(hr,0.0,1.0).or.inInt(hr,3.0,4.0)) then
      sleep = 0
    elseif (inInt(hr,22.,22.03).or.inInt(hr,1.0,1.03).or.inInt(hr,4.0,4.03)) then
      sleep = 1
    else
      sleep = 2
    end if
  end function

  logical function inInt(x,a,b) !x \in (a,b]
    real :: x,a,b
    inInt = (x .le. b) .and. (x.gt.a)
  end function

  !TEIbc index - available at Math Model for Thermal
  !Dind = TEIbc derivative
  subroutine TEIbc(Temp,RH,V,Tind,Dind)
    real :: Temp(:),RH,V
    real :: Tind(:),Dind(:)
    Tind=45.6026 - 2.3107*Temp - 0.3683*RH + 9.7092*V + 0.05492*Temp**2 + &
         0.00121*RH**2 + 0.66329*V**2 + 0.0128968*Temp*RH - 0.300928*Temp*V - &
         0.05952*RH*V
    Tind = min(max(Tind-27.,0.)/3.,1.)
    Dind =  -2.3107 + 0.10984*Temp + 0.0128968*RH - 0.300928*V
  end subroutine

  !inital chicken location
  subroutine initLoc(nf, hu,  th,  x1, x2, tind, Dind, Tc, ho, co, sl,  dd,&
                v1, v2,  f1,  f2, fl1,fl2,  fw1, fw2, fr1,fr2,neE1,neE2,&
                diE,neW1,neW2,diW,px,py,ind)
    integer :: nf
    REAL, allocatable ::   hu(:), th(:), x1(:), x2(:), tind(:), Dind(:), Tc(:),&
      ho(:), co(:), sl(:), dd(:), v1(:), v2(:), f1(:), f2(:),   fl1(:), fl2(:),&
      fw1(:),fw2(:),fr1(:),fr2(:),neE1(:),neE2(:),diE(:),neW1(:),neW2(:),diW(:)
    logical, allocatable :: ind(:)
    integer, allocatable :: px(:),py(:)
    allocate(hu(nf), th(nf), x1(nf), x2(nf))
    allocate(tind(nf), Dind(nf), Tc(nf))
    allocate(ho(nf), co(nf), sl(nf), dd(nf), v1(nf),   v2(nf),   f1(nf))
    allocate(f2(nf), fl1(nf),fl2(nf),fw1(nf),fw2(nf),  fr1(nf),  fr2(nf))
    allocate(neE1(nf),neE2(nf),diE(nf),neW1(nf),neW2(nf),diW(nf))
    allocate(px(nf),py(nf))
    allocate(ind(nf))
    ho=0.0; co=0.0; sl=0.0; dd=0.0; v1=0.0; v2=0.0; f1=0.0; f2=0.0
    fl1=0.0; fl2=0.0; fw1=0.0; fw2 =0.0
  end subroutine

  !deallocate memory of varibles initiated at initLoc routine
  subroutine dealoc_loc(hu,  th,  x1, x2, tind, Dind, Tc, ho, co, sl,  dd,&
                v1, v2,  f1,  f2, fl1,fl2,  fw1, fw2, fr1,fr2,neE1,neE2,&
                diE,neW1,neW2,diW,px,py,ind)
    REAL, allocatable :: hu(:), th(:), x1(:), x2(:), tind(:), Dind(:), Tc(:)
    real, allocatable :: ho(:) , co(:) , sl(:) , dd(:)  ,&
            v1(:) , v2(:) , f1(:) , f2(:)  ,&
            fl1(:), fl2(:), fw1(:), fw2(:) ,&
            fr1(:),fr2(:)
    integer , allocatable :: px(:),py(:)
    logical, allocatable :: ind(:)
    real, allocatable :: neE1(:),neE2(:),diE(:),neW1(:),neW2(:),diW(:)
    deallocate(hu, th, x1, x2, tind, Dind, Tc,&
             ho, co, sl, dd, v1,   v2,   f1,&
             f2, fl1,fl2,fw1,fw2,  fr1,  fr2,&
             neE1,neE2,diE,neW1,neW2,diW,px,py,ind)
  end subroutine

  !find nearest eater or water fountain
  subroutine nearest(x1,x2,n,esp,ind,y1,y2,ne1,ne2,di)
    real,intent(in) :: x1(:),x2(:), y1,y2,esp
    integer,intent(in) :: n
    logical,intent(in) :: ind(:)
    real, intent(out) :: ne1(:),ne2(:),di(:)
    ne1 = max(min(int(x1/esp),n),1)*esp-x1
    where(ind); ne2=y1-x2; ELSEwhere; ne2=y2-x2; ENDwhere
    di=sqrt(ne1**2+ne2**2)
  end subroutine
end module

!test program used to evalute objective function of the optimization problem
program teste
  use simulacao
  implicit NONE
  real :: x(30), f
  x = initInput()
  f = func(x)
  write(*,*) 'resultado :: ', f
end program
