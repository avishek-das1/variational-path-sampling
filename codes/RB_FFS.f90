!Forward Flux Sampling with the Rosenbluth algorithm
program ffs
use mpi
implicit none

integer,parameter:: Nw=20000 !number of initial conditions at R_A
double precision,parameter:: vo=9.0d0
integer, parameter:: N=80
double precision, parameter:: tau=1d-5 !timestep
double precision:: x0(N,2), theta0(N),xcont(N,2),thetacont(N)

!MPI
integer :: ierr,tag,cores,id,my_id,status(MPI_STATUS_SIZE),iexit !mpi stuff
double precision :: start,finish !keeps time

!initialization
integer, parameter:: frames=10000
integer:: iframe
double precision:: init(frames,N,3)
double precision:: rtemp !random number

!counters
integer:: i,j,lamcount,kcount,tempcounter

!ffs
integer,parameter:: nlam=12 !number of interfaces
double precision,parameter,dimension(nlam):: lams=(/1.25d0,1.29d0,1.33d0,&
        & 1.38d0,1.43d0,1.50d0,1.57d0,1.65d0,1.72d0,1.77d0,1.81d0,1.85d0/)
integer,parameter:: kvals=100 !each interface starts kvals trajs
integer:: Sarray(nlam-1) !number of successful trajs
integer:: indicator(kvals) !within one branching step lists which trajectories make it
double precision:: ttotA !total time spent in A during initialization
integer:: trajinit,trajlens(2) !time spent per individual trajectory doing initiatlization and propagation
integer,parameter:: dtmax=2d4 !maximum length of each segment
double precision:: trajx(kvals,N,2),trajtheta(kvals,N)
integer:: trajlen(kvals) !length of the trajectory segments
integer:: iinit !random number for which endpoint to start the next trajectories from
double precision:: xlam1(N,2),thetalam1(N) !coordinates at lam1
integer,parameter:: maxlam1=1d7 !maximum time spent in A
integer:: ind1,ind2 !indicator for crossing lam1 or reaching lamend
integer:: breakind,sumind !temporary variables

!random number stuff
integer :: sizerand
integer, allocatable :: putrand(:)

!mpi stuff
call MPI_Init(ierr)
call MPI_Comm_rank(MPI_COMM_WORLD,my_id,ierr)
call MPI_Comm_size(MPI_COMM_WORLD,cores,ierr)

if (cores/=kvals) then
    print *, 'number of threads must be same as kvals'
    stop
end if

!seed random numbers
call random_seed(size=sizerand)
allocate(putrand(sizerand))
putrand(:)=my_id+10
call random_seed(put=putrand)
do i=1,50
    call random_number(rtemp)
end do

if (my_id==0) then
    !read initial conditions
    open(unit=140,file='v9_frames.txt')
    do i=1,frames
      do j=1,N
          read(140,*)init(i,j,:)
      end do
    end do

    !output files
    open(unit=10,file='N0_ttot.txt')
    open(unit=20,file='Stot.txt')
    open(unit=30,file='dttot.txt')
end if

ttotA=0.0d0
tempcounter=0
ind1=2
do j=1,Nw

start=MPI_Wtime()
if (my_id==0) then !initialize only in main node

!start the initial segment loop here
ind2=0
trajinit=0
do while(ind1==2 .or. ind2==0) !denotes it ended in B well or that old traj has to be continued
if (ind1==2) then!choose a new frame
    call random_number(rtemp)
    iframe=floor(rtemp*frames)+1 !an integer between [1,frames]
    x0=init(iframe,:,1:2)
    theta0=init(iframe,:,3)
end if !else x0,theta0 already has previous info
!call a trajectory
call inittraj(N,tau,x0,theta0,vo,my_id,j,lams(1),lams(nlam),maxlam1,&
&       xlam1,thetalam1,trajlen(1),ind1,ind2,tempcounter)
!initial flux calculation, outputs total time instants spent, the initial traj and its length
if (ind1==2) write (*,*) 'ended in B for my_id,Nwcount',my_id,j

ttotA=ttotA+(trajlen(1)-1)*tau
trajinit=trajinit+trajlen(1)-1
tempcounter=tempcounter+1

end do !end initialization loop

if (j<Nw) then
    x0=xlam1
    theta0=thetalam1
end if

!AD
write(*,*) 'Initialization done for Nwcount,tempcounter,iframe',j,tempcounter,iframe
!call flush()
trajlens(:)=0

end if !other cores don't need to do initialization

breakind=0
do lamcount=2,nlam
if (breakind==1) exit
!trajs start from interface at lams(lamcount-1)

if (my_id==0) then !choose the starting point at main node
!choose a starting point randomly
if (lamcount==2) then
    xcont=xlam1
    thetacont=thetalam1
else
    do while(.true.)
        call random_number(rtemp)
        iinit=floor(rtemp*kvals)+1 !an integer between [1,kvals]
        if (indicator(iinit)==1) then
            xcont=trajx(iinit,:,:)
            thetacont=trajtheta(iinit,:)
            trajlens(2)=trajlens(2)+trajlen(iinit)-1
            exit
        end if
    end do
end if

end if !main node condition ends
!broadcast the init cond
call MPI_Barrier(MPI_COMM_WORLD,ierr)
call MPI_Bcast(xcont,N*2,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(thetacont,N,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
call MPI_Barrier(MPI_COMM_WORLD,ierr)

indicator=0
call proptraj(N,tau,xcont,thetacont,&
        &       vo,my_id,j,lams(lamcount),lams(1),dtmax,&
        &       trajx(my_id+1,:,:),trajtheta(my_id+1,:),&
        &       trajlen(my_id+1),indicator(my_id+1),&
        &       lamcount,my_id+1)
if (trajlen(my_id+1)==dtmax .and. indicator(my_id+1)==0) then
    write(*,*) 'Short segment at my_id,Nwcount,lamcount,kcount',j,lamcount,my_id+1
end if
!print *, 'propagation done for Nw,lamcount,kcount',j,lamcount,my_id+1
!receive in main node trajx,trajtheta,trajlen,indicator from other nodes
if (my_id>0) then
    id=0
    call MPI_Send(trajx(my_id+1,:,:),N*2,MPI_DOUBLE,id,100*my_id+4*(lamcount-2)+1,MPI_COMM_WORLD,ierr)
    call MPI_Send(trajtheta(my_id+1,:),N,MPI_DOUBLE,id,100*my_id+4*(lamcount-2)+2,MPI_COMM_WORLD,ierr)
    call MPI_Send(indicator(my_id+1),1,MPI_INTEGER,id,100*my_id+4*(lamcount-2)+3,MPI_COMM_WORLD,ierr)
    call MPI_Send(trajlen(my_id+1),1,MPI_INTEGER,id,100*my_id+4*(lamcount-2)+4,MPI_COMM_WORLD,ierr)
else
    do i=1,cores-1
        id=i
        call MPI_Recv(trajx(id+1,:,:),N*2,MPI_DOUBLE,id,100*id+4*(lamcount-2)+1,MPI_COMM_WORLD,status,ierr)
        call MPI_Recv(trajtheta(id+1,:),N,MPI_DOUBLE,id,100*id+4*(lamcount-2)+2,MPI_COMM_WORLD,status,ierr)
        call MPI_Recv(indicator(id+1),1,MPI_INTEGER,id,100*id+4*(lamcount-2)+3,MPI_COMM_WORLD,status,ierr)
        call MPI_Recv(trajlen(id+1),1,MPI_INTEGER,id,100*id+4*(lamcount-2)+4,MPI_COMM_WORLD,status,ierr)
        !print *, 'received for core', i
    end do
    trajlens(1)=trajlens(1)+sum(trajlen)-kvals
    !write the results of this step
    sumind=sum(indicator)
    write(20,'(I5.1)',advance='no') sumind
    !print *, j,lamcount,sumind
    !if no trajectories make it, break the lamcount loop and start fresh
    if (sumind==0) then
        do i=lamcount+1,nlam
            write(20,'(I5.1)',advance='no') 0
        end do
        breakind=1
    end if

end if !printed results of this lamcount loop
!broadcast breakind
call MPI_Barrier(MPI_COMM_WORLD,ierr)
call MPI_Bcast(breakind,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Barrier(MPI_COMM_WORLD,ierr)

!print *, 'kloop done for lamcount,my_id is', lamcount,my_id
!call flush()
end do !end lamcount loop

finish=MPI_Wtime()
if (my_id==0) then
    print *, 'Nwcount, time taken is',j,finish-start
    write(20,*) !add newline
    write(30,*) trajinit,trajlens(:)
    !call flush()
end if

end do !end Nw loop
if (my_id==0) write(10,*) Nw,ttotA
!call flush()

call MPI_Finalize(ierr)

end program ffs
!*************************************************************************!
!*************************************************************************!
!*************************************************************************!
!subroutine for running a single trajectory
subroutine inittraj(N,tau,x0,theta0,vo,my_id,Nwcount,lam1,lamend,maxlam1,&
        &       xlam1,thetalam1,trajlen,ind1,ind2,tempcounter)

implicit none

!arguments
integer, intent(in):: N
double precision, intent(in):: x0(N,2),theta0(N)
double precision, intent(in):: lam1,lamend,vo,tau
integer, intent(in):: my_id,Nwcount,maxlam1,tempcounter
double precision:: xlam1(N,2),thetalam1(N)
integer:: trajlen,ind1,ind2

!trajectory parameters
!integer, parameter:: N=80
double precision, parameter:: sigma=1.0d0
double precision, parameter:: sigmap6=sigma**6, sigmap12=sigmap6**2
double precision, parameter:: cutwca=2**(1./6) *sigma
double precision, parameter:: rthick=0.5d0*sigma !rskin-rcut
double precision, parameter:: rskin=cutwca+rthick
double precision, parameter:: Dd=0.25d0,Db=0.5d0,D_theta=1.0d0 !diffusion constant for dimer and bath
double precision, parameter:: fd=2.0d0,fb=1.0d0 !friction coefficient for dimer and bath
double precision:: etavard,etavarb,eta_theta_var!prefactor for gaussian noise
double precision,parameter:: pi=4.0d0*atan(1.0d0)
double precision, parameter:: rho=0.6 !A density
!double precision, parameter:: L=(NA*pi* (sigmaA)**2/(4.0d0*rho))**(1.0d0/2.0d0)
double precision:: L!=(N/rho)**(1.0d0/2.0d0)

!force calc
double precision,dimension(N,2):: delx !holder for jump in every timestep
double precision,dimension(N):: deltheta
double precision,dimension(N,2):: F0,F0wca,F0dw,Fa !force
double precision:: potdw,rint,rintold
double precision:: frc

!random number stuff
integer :: sizerand
integer, allocatable :: putrand(:)
!double precision:: rand1(N,2),rand2(N,2),r(N,2),theta(N,2),eta(N,2)
double precision::rand1(N,2),rand2(N,2),rand3(N),rand4(N),r(N,2),thetaBM1(N,2),thetaBM2(N),eta(N,2),eta_theta(N)

!neighbor list
double precision, dimension(N,2):: dxnei
double precision:: drnei,drneimax,drneimax2
integer:: nnflag
integer,dimension(N,N):: nlist

!counters and temporary variables
integer:: i,j,k,p
double precision:: xr(2),absr
integer:: indicator,counter,boxx

!determine length of box
L=(N/rho)**(1.0d0/2.0d0)
etavard=(2.0d0*Dd*tau)**0.5d0
etavarb=(2.0d0*Db*tau)**0.5d0 
eta_theta_var=(2.0d0*D_theta*tau)**0.5d0

!seed random numbers
call random_seed(size=sizerand)
allocate(putrand(sizerand))
putrand(:)=1d6*my_id+1d3*Nwcount+tempcounter
call random_seed(put=putrand)
do i=1,50
    call random_number(rand1(1,1))
end do
!x0(1,:)=0.0d0!0.5*L !first particle in the middle of the box
!x0(2,1)=0.0d0+0.25*2.0d0
!x0(2,2)=0.0d0
!do i=3,N
!!    print *, 'starting particle', i
!    indicator=1
!    do while (indicator==1)
!        call random_number(rand1(i,:)) !three random coordinates for ith particle
!        rand1(i,:)=(rand1(i,:)-0.5d0)*L !rescales the coordinates to the whole box
!        indicator=0
!        do j=1,i-1  !Check overlap with previous particles
!            xr=rand1(i,:)-x0(j,:)
!            !!periodic boundaries
!            xr(1)=xr(1)-L*nint(xr(1)/L)
!            xr(2)=xr(2)-L*nint(xr(2)/L)
!            absr=dot_product(xr,xr)**0.5
!            if (absr<0.95*sigma) indicator=1  !overlap!
!        end do
!    end do
!    x0(i,:)=rand1(i,:)
!    if (my_id==0) print *, 'ending particle', i
!end do

!boxx=floor(L/sigma) !side length of box
!counter=0
!do i=1,boxx
!    do j=1,boxx
!        counter=counter+1
!        x0(counter,1)=-0.5*L+sigma*(i-1)
!        x0(counter,2)=-0.5*L+sigma*(j-1)
!        if (my_id==0) print *, counter, x0(counter,:)
!        if (counter==N) exit
!    end do
!    if (counter==N) exit
!end do

frc=24.0d0*(2.0d0*sigmap12/cutwca**13-sigmap6/cutwca**7)
dxnei=1000.0d0
xlam1=x0
thetalam1=theta0
ind1=0
counter=0

!begin actual transition path trajectory
do while (.true.)
       
    counter=counter+1

    !generate random noise using Box Mueller algo
    if (mod(counter,2)==1) then
        call random_number(rand1)
        call random_number(rand2)
        r=(-2.0d0*log(rand2))**0.5d0
        thetaBM1=2.0d0*pi*rand1
        rand1=r*cos(thetaBM1)
        rand2=r*sin(thetaBM1)
        eta=rand1

        call random_number(rand3)
        call random_number(rand4)
        r(:,1)=(-2.0d0*log(rand4))**0.5d0
        thetaBM2=2.0d0*pi*rand3
        rand3=r(:,1)*cos(thetaBM2)
        rand4=r(:,1)*sin(thetaBM2)
        eta_theta=rand3
    else
        eta=rand2
        eta_theta=rand4
    end if
    
    !Updating neighbor list if required by calculating two maximum displacements
    drneimax=0.0
    drneimax2=0.0
    do j=1,N
        drnei=sqrt(dot_product(dxnei(j,:),dxnei(j,:)))
        if (drnei > drneimax) then
            drneimax2=drneimax
            drneimax=drnei
        else
            if (drnei > drneimax2) then
                drneimax2=drnei
            endif
        endif
    enddo
    if (drneimax+drneimax2 > rthick) then
        nnflag=1 !update neighbor list in the next call for force
    !print *,'Updating neighbor list at step', i
    end if

    if (nnflag==1) dxnei=0

    if (counter>1) rintold=rint

    !force calc
    call force(N,sigma,sigmap6,sigmap12,cutwca,&
    &   xlam1,F0,F0wca,F0dw,rint,L,frc,&
    &   rskin,nlist,nnflag)

    if (counter>1) then
        if (rint>lam1 .and. rintold<lam1) then
            ind1=1 !crossed the first interface left to right
            ind2=1
            exit
        end if
    end if

    if (rint>lamend) then
        ind1=2 !reached lam2 in previous loop
        exit
    end if

    if (counter==maxlam1+1) then
        write(*,*) 'my_id,Nwcount',my_id,Nwcount,'initial segment unsuccessful'
        ind1=2 !need to use a different initial condition !!this biases the trajectory ensemble
        exit
    end if

    Fa(:,1)=vo*cos(thetalam1)
    Fa(:,2)=vo*sin(thetalam1)
    Fa(1,:)=0.0d0
    Fa(2,:)=0.0d0
    F0=F0+Fa

    !propagate position
    delx(3:N,:)=tau*F0(3:N,:)/fb+etavarb*eta(3:N,:)
    delx(1:2,:)=tau*F0(1:2,:)/fd+etavard*eta(1:2,:)
    xlam1=xlam1+delx
    dxnei=dxnei+delx
   
    deltheta=eta_theta_var*eta_theta
    thetalam1=thetalam1+deltheta
 
    do j=1,N
        xlam1(j,1)=xlam1(j,1)-L*nint(xlam1(j,1)/L)
        xlam1(j,2)=xlam1(j,2)-L*nint(xlam1(j,2)/L)
    end do

end do !trajectory has reached lam1
trajlen=counter

end subroutine inittraj
!*************************************************************************!
!*************************************************************************!
!*************************************************************************!
!subroutine for running a single trajectory
subroutine proptraj(N,tau,x0,theta0,vo,my_id,Nwcount,lam,lam1,dtmax,&
        &       trajx,trajtheta,trajlen,indicator,lamcount,kcount)

implicit none

!arguments
integer, intent(in):: N
double precision, intent(in):: x0(N,2),theta0(N)
double precision, intent(in):: lam,lam1,vo,tau
integer, intent(in):: my_id,Nwcount,dtmax,lamcount,kcount
double precision:: trajx(N,2),trajtheta(N)
integer:: trajlen,indicator

!trajectory parameters
!integer, parameter:: N=80
double precision, parameter:: sigma=1.0d0
double precision, parameter:: sigmap6=sigma**6, sigmap12=sigmap6**2
double precision, parameter:: cutwca=2**(1./6) *sigma
double precision, parameter:: rthick=0.5d0*sigma !rskin-rcut
double precision, parameter:: rskin=cutwca+rthick
double precision, parameter:: Dd=0.25d0,Db=0.5d0,D_theta=1.0d0 !diffusion constant for dimer and bath
double precision, parameter:: fd=2.0d0,fb=1.0d0 !friction coefficient for dimer and bath
double precision :: etavard,etavarb,eta_theta_var!prefactor for gaussian noise
double precision,parameter:: pi=4.0d0*atan(1.0d0)
double precision, parameter:: rho=0.6 !A density
!double precision, parameter:: L=(NA*pi* (sigmaA)**2/(4.0d0*rho))**(1.0d0/2.0d0)
double precision:: L!=(N/rho)**(1.0d0/2.0d0)

!force calc
double precision,dimension(N,2):: delx !holder for jump in every timestep
double precision,dimension(N):: deltheta
double precision,dimension(N,2):: F0,F0wca,F0dw,Fa !force
double precision:: potdw,rint
double precision:: frc

!random number stuff
integer :: sizerand
integer, allocatable :: putrand(:)
!double precision:: rand1(N,2),rand2(N,2),r(N,2),theta(N,2),eta(N,2)
double precision::rand1(N,2),rand2(N,2),rand3(N),rand4(N),r(N,2),thetaBM1(N,2),thetaBM2(N),eta(N,2),eta_theta(N)

!neighbor list
double precision, dimension(N,2):: dxnei
double precision:: drnei,drneimax,drneimax2
integer:: nnflag
integer,dimension(N,N):: nlist

!counters and temporary variables
integer:: i,j,k,p
double precision:: xr(2),absr
integer:: ind1,ind2,counter,boxx

!determine length of box
L=(N/rho)**(1.0d0/2.0d0)
etavard=(2.0d0*Dd*tau)**0.5d0
etavarb=(2.0d0*Db*tau)**0.5d0 
eta_theta_var=(2.0d0*D_theta*tau)**0.5d0

!seed random numbers
call random_seed(size=sizerand)
allocate(putrand(sizerand))
putrand(:)=1d6*my_id+1d4*Nwcount+1d2*lamcount+kcount !AD makes sense when kvals=100
call random_seed(put=putrand)
do i=1,50
    call random_number(rand1(1,1))
end do

!x0(1,:)=0.0d0!0.5*L !first particle in the middle of the box
!x0(2,1)=0.0d0+0.25*2.0d0
!x0(2,2)=0.0d0
!do i=3,N
!!    print *, 'starting particle', i
!    indicator=1
!    do while (indicator==1)
!        call random_number(rand1(i,:)) !three random coordinates for ith particle
!        rand1(i,:)=(rand1(i,:)-0.5d0)*L !rescales the coordinates to the whole box
!        indicator=0
!        do j=1,i-1  !Check overlap with previous particles
!            xr=rand1(i,:)-x0(j,:)
!            !!periodic boundaries
!            xr(1)=xr(1)-L*nint(xr(1)/L)
!            xr(2)=xr(2)-L*nint(xr(2)/L)
!            absr=dot_product(xr,xr)**0.5
!            if (absr<0.95*sigma) indicator=1  !overlap!
!        end do
!    end do
!    x0(i,:)=rand1(i,:)
!    if (my_id==0) print *, 'ending particle', i
!end do

!boxx=floor(L/sigma) !side length of box
!counter=0
!do i=1,boxx
!    do j=1,boxx
!        counter=counter+1
!        x0(counter,1)=-0.5*L+sigma*(i-1)
!        x0(counter,2)=-0.5*L+sigma*(j-1)
!        if (my_id==0) print *, counter, x0(counter,:)
!        if (counter==N) exit
!    end do
!    if (counter==N) exit
!end do

frc=24.0d0*(2.0d0*sigmap12/cutwca**13-sigmap6/cutwca**7)
dxnei=1000.0d0
trajx=x0
trajtheta=theta0
ind1=0 !this notes whether lam1 has been reached
indicator=0 !this notes whether lam has been reached
counter=0

!begin actual transition path trajectory
do while (.true.) !hasn't reached lam1 or lam
        
    counter=counter+1

    !generate random noise using Box Mueller algo
    if (mod(counter,2)==1) then
        call random_number(rand1)
        call random_number(rand2)
        r=(-2.0d0*log(rand2))**0.5d0
        thetaBM1=2.0d0*pi*rand1
        rand1=r*cos(thetaBM1)
        rand2=r*sin(thetaBM1)
        eta=rand1

        call random_number(rand3)
        call random_number(rand4)
        r(:,1)=(-2.0d0*log(rand4))**0.5d0
        thetaBM2=2.0d0*pi*rand3
        rand3=r(:,1)*cos(thetaBM2)
        rand4=r(:,1)*sin(thetaBM2)
        eta_theta=rand3
    else
        eta=rand2
        eta_theta=rand4
    end if
    
    !Updating neighbor list if required by calculating two maximum displacements
    drneimax=0.0
    drneimax2=0.0
    do j=1,N
        drnei=sqrt(dot_product(dxnei(j,:),dxnei(j,:)))
        if (drnei > drneimax) then
            drneimax2=drneimax
            drneimax=drnei
        else
            if (drnei > drneimax2) then
                drneimax2=drnei
            endif
        endif
    enddo
    if (drneimax+drneimax2 > rthick) then
        nnflag=1 !update neighbor list in the next call for force
    !print *,'Updating neighbor list at step', i
    end if

    if (nnflag==1) dxnei=0

    !force calc
    call force(N,sigma,sigmap6,sigmap12,cutwca,&
    &   trajx,F0,F0wca,F0dw,rint,L,frc,&
    &   rskin,nlist,nnflag)

    if (rint<lam1) then 
        ind1=1 !reached lam1 in previous loop
        exit
    else if (rint>lam) then
        indicator=1 !reached lam in previous loop
        exit
    end if

    if (counter==dtmax+1) exit

    Fa(:,1)=vo*cos(trajtheta)
    Fa(:,2)=vo*sin(trajtheta)
    Fa(1,:)=0.0d0
    Fa(2,:)=0.0d0
    F0=F0+Fa

    !propagate position
    delx(3:N,:)=tau*F0(3:N,:)/fb+etavarb*eta(3:N,:)
    delx(1:2,:)=tau*F0(1:2,:)/fd+etavard*eta(1:2,:)
    trajx=trajx+delx
    dxnei=dxnei+delx
   
    deltheta=eta_theta_var*eta_theta
    trajtheta=trajtheta+deltheta
 
    do j=1,N
        trajx(j,1)=trajx(j,1)-L*nint(trajx(j,1)/L)
        trajx(j,2)=trajx(j,2)-L*nint(trajx(j,2)/L)
    end do

end do !trajectory has reached lam or lam1
trajlen=counter

end subroutine proptraj
!*************************************************************************!
!*************************************************************************!
!*************************************************************************!
!*************************************************************************!
!subroutine for force calculation
subroutine force(N,sigma,sigmap6,sigmap12,cutwca,&
    &   x,F0,F0wca,F0dw,rint,L,frc,&
    &   rskin,nlist,nnflag)
implicit none

integer, intent(in):: N
double precision, intent(in):: x(N,2),L,frc,rskin
integer:: nnflag, nlist(N,N)
double precision, intent(in):: sigma,sigmap6,sigmap12,cutwca
double precision:: F0(N,2),F0wca(N,2),F0dw(N,2)
double precision:: potdw,rint
integer:: i,j,k

double precision:: xr(2),rcap(2),absr,rprod1,rprod2,Fint(2),fac1,fac2
double precision, parameter:: ewca=1.0d0,hdw=3.5d0,wdw=0.45d0

F0=0.0d0
F0wca=0.d0
F0dw=0.0d0

!generic forces common for all particles
if (nnflag==0) then !do calculation according to neighbor list

do i=1,N
!put 1body force here

!do j=i+1,N
do k=1,nlist(i,1)
  j=nlist(i,k+1)
  if (j<=i) cycle

!!define absr and xr right away and work with those
xr=x(i,:)-x(j,:)
!!periodic boundaries
xr(1)=xr(1)-nint(xr(1)/L)*L
xr(2)=xr(2)-nint(xr(2)/L)*L
absr=dot_product(xr,xr)**0.5
rcap=xr/absr

!!WCA evaluation
if (absr<cutwca) then
rprod1=absr*absr
rprod1=rprod1*rprod1
rprod1=rprod1*rprod1/absr !7th power
rprod2=rprod1*rprod1/absr !13th power
Fint=ewca*(24.0d0*(2.0d0*sigmap12/rprod2-sigmap6/rprod1)-frc)*rcap
F0wca(i,:)=F0wca(i,:)+Fint
F0wca(j,:)=F0wca(j,:)-Fint
end if
end do

end do

else !make new neighbor list
nlist(:,1)=0

do i=1,N
!put 1body force here

do j=i+1,N

!!define absr and xr right away and work with those
xr=x(i,:)-x(j,:)
!!periodic boundaries
xr(1)=xr(1)-nint(xr(1)/L)*L
xr(2)=xr(2)-nint(xr(2)/L)*L
absr=dot_product(xr,xr)**0.5
rcap=xr/absr

if (absr<=rskin) then !i and j are neighbors
nlist(i,nlist(i,1)+2)=j
nlist(j,nlist(j,1)+2)=i
nlist(i,1)=nlist(i,1)+1
nlist(j,1)=nlist(j,1)+1

!!WCA evaluation
if (absr<cutwca) then
rprod1=absr*absr
rprod1=rprod1*rprod1
rprod1=rprod1*rprod1/absr !7th power
rprod2=rprod1*rprod1/absr !13th power
Fint=ewca*(24.0d0*(2.0d0*sigmap12/rprod2-sigmap6/rprod1)-frc)*rcap
F0wca(i,:)=F0wca(i,:)+Fint
F0wca(j,:)=F0wca(j,:)-Fint
end if
end if !end calculation for i and j being neighbors

end do !loop over 2body forces

end do !loop over 1body forces
nnflag=0
end if

!double well forces for first two particles
xr=x(1,:)-x(2,:)
!!periodic boundaries
xr(1)=xr(1)-nint(xr(1)/L)*L
xr(2)=xr(2)-nint(xr(2)/L)*L
absr=dot_product(xr,xr)**0.5
rcap=xr/absr
fac1=absr-cutwca-wdw
fac2=1.0d0-(fac1)**2/wdw/wdw
rint=absr
potdw=hdw*fac2*fac2
Fint=4.0d0*hdw/wdw/wdw *fac1*fac2*rcap
F0dw(1,:)=Fint
F0dw(2,:)=-Fint

F0=F0wca+F0dw

end subroutine force
!*************************************************************************!
!*************************************************************************!
!*************************************************************************!
!*************************************************************************!
!subroutine for calculating total potential
subroutine potential(N,sigma,sigmap6,sigmap12,cutwca,x0,L,frc,pot0)
implicit none
integer, intent(in):: N
double precision, intent(in)::sigma,sigmap6,sigmap12,cutwca,x0(N,2),L,frc
double precision:: urc,pot0
integer:: i,j,k

double precision:: xr(2),absr,rprod1,rprod2,fac1,fac2
double precision, parameter:: ewca=1.0d0,hdw=3.5d0,wdw=0.45d0

urc=-1.0d0 !for wca potential
pot0=0.0d0

do i=1,N
do j=i+1,N

!!define absr and xr right away and work with those
xr=x0(i,:)-x0(j,:)
!!periodic boundaries
xr(1)=xr(1)-nint(xr(1)/L)*L
xr(2)=xr(2)-nint(xr(2)/L)*L
absr=dot_product(xr,xr)**0.5

!!WCA evaluation
if (absr<cutwca) then
rprod1=absr**6
rprod2=rprod1*rprod1
pot0=pot0+ewca*(4.0d0*(sigmap12/rprod2-sigmap6/rprod1)+(absr-cutwca)*frc-urc)
end if
end do
end do

!double well forces for first two particles
xr=x0(1,:)-x0(2,:)
!!periodic boundaries
xr(1)=xr(1)-nint(xr(1)/L)*L
xr(2)=xr(2)-nint(xr(2)/L)*L
absr=dot_product(xr,xr)**0.5
fac1=absr-cutwca-wdw
fac2=1.0d0-(fac1)**2/wdw/wdw
pot0=pot0+hdw*fac2*fac2

end subroutine potential
!*************************************************************************!
!*************************************************************************!
!*************************************************************************!
!*************************************************************************!

