!generating initial configurations for Variational Path Sampling, minimally modified from mcvb.f90
!"Direct evaluation of rare events in active matter from variational path sampling", Avishek Das, Benjamin Kuznets-Speck, David T. Limmer
!Passive dimer in 2d active bath
!main program, manages the mpi side of things, calls subroutine "traj"
program mcvb
    use mpi
    implicit none

    double precision,parameter:: s=-100 !bias
    double precision,parameter :: vo = 9.0d0 !magnitude for self-propulsion
    integer, parameter:: N=80,Mr=1,Mt=1 ! Total number of basis is Mr*Mt
    double precision:: chi(Mr,Mt) !coefficients for force as a function of x,t
    double precision:: omega(4),omegan(4) !4 elements are the full estimator, average of indicator, ordinary average of \Delta U, conditioned average of \Delta U
    double precision:: delomega(Mr,Mt),delomegan(Mr,Mt) !gradient of estimator
    double precision:: r0(N,2), theta0(N) !initial positions and active director angles

    !value function data structures
    integer, parameter:: Mvr=Mr,Mvt=Mt !Number of basis is Mvr*Mvt for the value function !Downstream code assumes Mvr=Mr and Mvt=Mt
    double precision:: delL(Mvr,Mvt),delLn(Mvr,Mvt),psi(Mvr,Mvt) !gradients of lossfunction for value, and parameters for value function

    !MPI
    integer:: Nw !number of trajectories
    integer :: ierr,tag,cores,id,my_id,status(MPI_STATUS_SIZE),npercore,iexit !mpi parameters
    double precision :: start,finish !keeps time

    !Learning
    double precision:: alphachi,alphapsi !learning rates of force and value function

    !initial frames sampled from steady state
    integer, parameter:: frames=1 !total number of initial bath configurations
    integer:: iframe !counter for random selection
    double precision:: init(frames,N,3) !reads in the coords and the theta from all frames
    double precision:: rtemp !random number for selecting frame

    !counters
    integer:: i,j,npercorecount,loopcount

    !mpi initialization
    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD,my_id,ierr)
    call MPI_Comm_size(MPI_COMM_WORLD,cores,ierr)
    Nw=1
    npercore=floor(Nw*1.0d0/cores)
    Nw=npercore*cores

    !initialize coefficients
    if (my_id==0) then
        !open(unit=10,file='chi0.txt')
        !open(unit=12,file='psi0.txt')
        !  do i=1,Mr
        !    read(10,*) chi(i,:)
        !    read(12,*) psi(i,:)
        !  end do
        !open(unit=20,file='omega.txt')
        !open(unit=120,file='chirunning.txt')
        !open(unit=130,file='psirunning.txt')
        !open(unit=140,file='v9_frames.txt')
        !do i=1,frames
        !    do j=1,N
        !        read(140,*)init(i,j,:)
        !    end do
        !end do
    end if
    chi=0.0d0
    psi=0.0d0
    omega=0.0d0
    delomega=0.0d0
    delL=0.0d0

    alphachi=0.0d0
    alphapsi=200.0d0

    !optimization loop
    do loopcount=0,0
        start=MPI_Wtime()

        if (loopcount==201) then
            alphachi=40.0d0
            alphapsi=200.0d0
        end if

        omegan=0.0d0
        delomegan=0.0d0
        delLn=0.0d0

        !grad descent step
        chi=chi-alphachi*delomega
        psi=psi+alphapsi*delL

        !Broadcast the new coefficients
        call MPI_Barrier(MPI_COMM_WORLD,ierr)
        call MPI_Bcast(chi,Mr*Mt,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(psi,Mvr*Mvt,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)

        !run trajectories and sum gradients
        do npercorecount=1,npercore
            !randomly choose and broadcast the initial coords
            if (my_id==0) then
                do i=1,cores-1
                    call random_number(rtemp)
                    iframe=floor(rtemp*frames)+1 !an integer between [1,frames]
                    r0=init(iframe,:,1:2)
                    theta0=init(iframe,:,3)
                    id=i
                    call MPI_Send(r0,N*2,MPI_DOUBLE,id,1000*i+1,MPI_COMM_WORLD,ierr)
                    call MPI_Send(theta0,N,MPI_DOUBLE,id,1000*i+2,MPI_COMM_WORLD,ierr)
                end do

                call random_number(rtemp)
                iframe=floor(rtemp*frames)+1 !an integer between [1,frames]
                r0=init(iframe,:,1:2)
                theta0=init(iframe,:,3)
            else
                id=0
                call MPI_Recv(r0,N*2,MPI_DOUBLE,id,1000*my_id+1,MPI_COMM_WORLD,status,ierr)
                call MPI_Recv(theta0,N,MPI_DOUBLE,id,1000*my_id+2,MPI_COMM_WORLD,status,ierr)
            end if

            call traj(N,r0,theta0,vo,Mr,Mt,Mvr,Mvt,chi,psi,s,omega,delomega,delL,my_id,npercorecount,loopcount)
            !print *, my_id,npercorecount,omega
            omegan=omegan+omega
            delomegan=delomegan+delomega
            delLn=delLn+delL
        end do !end npercore loop

        !!Send the estimator and gradients to the primary node
        if (my_id==0) then
            do i=1,cores-1
                id=i
                call MPI_Recv(omega,4,MPI_DOUBLE,id,100*i+1,MPI_COMM_WORLD,status,ierr)
                omegan=omegan+omega
                call MPI_Recv(delomega,Mr*Mt,MPI_DOUBLE,id,100*i+2,MPI_COMM_WORLD,status,ierr)
                delomegan=delomegan+delomega
                call MPI_Recv(delL,Mvr*Mvt,MPI_DOUBLE,id,100*i+11,MPI_COMM_WORLD,status,ierr)
                delLn=delLn+delL
            end do


            omega(4)=omegan(4)/omegan(2) !conditioned average of action
            omega(1)=omegan(1)/Nw !average of estimator
            omega(2)=omegan(2)/Nw !average of indicator
            omega(3)=omegan(3)/Nw !normal average of action
            delomega=delomegan/Nw
            delL=delLn/Nw

        else
            id=0
            call MPI_Send(omegan,4,MPI_DOUBLE,id,100*my_id+1,MPI_COMM_WORLD,ierr)
            call MPI_Send(delomegan,Mr*Mt,MPI_DOUBLE,id,100*my_id+2,MPI_COMM_WORLD,ierr)
            call MPI_Send(delLn,Mvr*Mvt,MPI_DOUBLE,id,100*my_id+11,MPI_COMM_WORLD,ierr)
        end if

        if (isnan(delomega(1,1)) .or. isnan(delL(1,1))) then
            delomega=0.0d0
            delL=0.0d0
            if (my_id==0) write(*,*) 'NaN gradients'
        end if

        !Broadcast the averaged quantities from primary node
        call MPI_Barrier(MPI_COMM_WORLD,ierr)
        call MPI_Bcast(omega,4,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(delomega,Mr*Mt,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(delL,Mvr*Mvt,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
        call MPI_Barrier(MPI_COMM_WORLD,ierr)

        !Output quantities in primary node
        if (my_id==0) then

            !Output the estimator
            !write(20,'(5(2x,E15.4E3))') omega(1), omega(2), omega(3), omega(4),log(omega(2))-omega(4)
            !!write(*,*) omega(1), omega(2), omega(3), omega(4),log(omega(2))+omega(4),omega(5)

            !Output the variational parameters
            !if (modulo(loopcount,50)==0) then
            !        do i=1,Mr
            !            write(120,'(41(2x,E15.4E3))') chi(i,:)
            !            write(130,'(41(2x,E15.4E3))') psi(i,:)
            !        end do
            !end if
            !call flush
        end if

        finish=MPI_Wtime()
        if (my_id==0) print *, 'Step, Nw, time taken is',loopcount,Nw,finish-start
    end do !end learning loop

    call MPI_Finalize(ierr)
end program mcvb
!*************************************************************************!
!*************************************************************************!
!*************************************************************************!
!subroutine for running a single trajectory, calls subroutine "force"
subroutine traj(N,r0,theta0,vo,Mr,Mt,Mvr,Mvt,chi,psi,s,omega,delomega,delL,my_id,npercorecount,loopcount)
    implicit none

    !arguments
    integer, intent(in):: N,Mr,Mt,Mvr,Mvt
    double precision,intent(in):: chi(Mr,Mt),r0(N,2), theta0(N)
    double precision, intent(in):: s,vo,psi(Mvr,Mvt)
    integer, intent(in):: my_id,npercorecount,loopcount
    double precision:: omega(4),delomega(Mr,Mt),delL(Mvr,Mvt)

    !trajectory parameters
    !integer, parameter:: N=80
    double precision, parameter:: sigma=1.0d0
    double precision, parameter:: sigmap6=sigma**6, sigmap12=sigmap6**2
    double precision, parameter:: rwca=2**(1./6) *sigma
    double precision, parameter:: rthick=0.5d0*sigma !for defining neighbor list
    double precision, parameter:: rskin=rwca+rthick
    integer, parameter:: steps=20000 !total number of timesteps
    double precision, parameter:: dt=1d-5 !timestep
    double precision, parameter:: temperature=0.5d0
    double precision, parameter:: gammad=2.0d0,gammab=1.0d0 !friction coefficient for dimer and bath
    double precision, parameter:: Dd=temperature/gammad,Db=temperature/gammab,D_theta=1.d0 !diffusion constant for dimer and bath and angular diffusion constant of directors
    double precision, parameter :: etavard=(2.0d0*Dd*dt)**0.5d0,etavarb=(2.0d0*Db*dt)**0.5d0 !prefactor for gaussian noise
    double precision, parameter :: eta_theta_var=(2.0d0*D_theta*dt)**0.5d0 !prefactor for rotational diffusion of director
    double precision, parameter:: pi=4.0d0*atan(1.0d0)
    double precision, parameter:: rho=0.6 !A density
    !double precision, parameter:: L=(NA*pi* (sigmaA)**2/(4.0d0*rho))**(1.0d0/2.0d0)
    double precision:: L!=(N/rho)**(1.0d0/2.0d0) !periodic box is [-L/2,L/2]

    !force calc, gradient structures
    double precision, dimension(N,2):: r1,r2,delr !holders for current coordinates and jump in every timestep
    double precision, dimension(N):: theta1,theta2,deltheta !holders for current director and jump in every timestep
    double precision, dimension(N,2):: F1,F !F1=F+lam where lam is the additional force
    double precision, dimension(N,2):: Fwca,Fdw,Fa !WCA, double-well and active forces
    double precision:: dellam(Mr,Mt,2) !contains partial derivative of additional force w.r.t. parameters, last rank is for component
    double precision:: gaussianr(Mr,2),gaussiant(Mt,2) !contains centers and variances of the gaussians for force basis
    double precision:: y(Mr,Mt),dely(Mr,Mt) !malliavin weight and its jump in every timestep
    double precision:: potdw,r12 !Double well potential and intra-dimer distance
    double precision:: frc !shifted forces cutoff for force, trivially 0 for WCA

    !random number generator data structures
    integer :: sizerand !size of seed
    integer, allocatable :: putrand(:) !seed
    double precision:: rand1(N,2),rand2(N,2),rand3(N),rand4(N) !uniform random number holders
    double precision:: rBM1(N,2),thetaBM1(N,2),rBM2(N),thetaBM2(N) !Box Mueller algorithm variables
    double precision:: eta(N,2),eta_theta(N) !gaussian noise in position and director angle

    !value function structures
    double precision:: V !value function
    double precision:: gradV(Mvr,Mvt) !gradient of value function
    !value function uses the same gaussian basis as the force
    double precision:: z(Mvr,Mvt) !time integrated gradV
    double precision:: xi !increment in esimator in every timestep

    !data structures for constructing and updating a neighbor list
    double precision, dimension(N,2):: dxnei !keeps track of displacements
    double precision:: drnei,drneimax,drneimax2 !keeps track of total displacements
    integer:: nnflag !flags on for updating neighbor list
    integer,dimension(N,N):: nlist !neighobor list

    !counters and temporary variables
    integer:: i,j,k,p,p2,counter,boxx
    double precision:: xr(2),absr

    !open file for trajectory
    open(unit=39,file='v9_frames.txt')

    !determine length of box
    L=(N/rho)**(1.0d0/2.0d0)

    !seed random numbers
    call random_seed(size=sizerand)
    allocate(putrand(sizerand))
    putrand(:)=1d6*my_id+npercorecount+1d3*loopcount
    call random_seed(put=putrand)

    !initialize all time integrals/averages
    omega=0.0d0
    delomega=0.0d0
    delL=0.0d0
    y=0.0d0
    z=0.0d0

    !Define gaussian basis centers and variance
    do j=1,Mt
      gaussiant(j,1)=1+(j-1)*(steps-1.0d0)/(Mt-1)
    end do
    gaussiant(:,2)=((steps-1.0d0)/(Mt-1)/2.0d0)**2 !variance

    do j=1,Mr
      gaussianr(j,1)=0.9+(j-1)*(2.3-0.9)/(Mr-1)
    end do
    gaussianr(:,2)=((2.3-0.9)/(Mr-1)/2.0d0)**2 !variance

    !shifted forces cutoff
    frc=24.0d0*(2.0d0*sigmap12/rwca**13-sigmap6/rwca**7)

    !initializing positions on square lattice
    boxx=floor(L/sigma) !side length of box
    counter=0
    do i=1,boxx
        do j=1,boxx
            counter=counter+1
            r0(counter,1)=-0.5*L+sigma*(i-1)
            r0(counter,2)=-0.5*L+sigma*(j-1)
            !if (my_id==0) print *, counter, r0(counter,:)
            if (counter==N) exit
        end do
        if (counter==N) exit
    end do

    !initializing directors as random
    call random_number(theta0)
    theta0=theta0*2.0d0*pi

    !initialize coordinates
    r2=r0
    theta2=theta0
    dxnei=1000.0d0

    !begin transition path trajectory
    p=0
    p2=0
    do while(p<10000)

        p2=p2+1

        !transfer new coordinates to old holde
        r1=r2
        theta1=theta2

        !generate random noise using Box Mueller algo
        if (mod(i,2)==1) then
 
            !for position
            call random_number(rand1)
            call random_number(rand2)
            rBM1=(-2.0d0*log(rand2))**0.5d0
            thetaBM1=2.0d0*pi*rand1
            rand1=rBM1*cos(thetaBM1)
            rand2=rBM1*sin(thetaBM1)
            eta=rand1

            !for angle
            call random_number(rand3)
            call random_number(rand4)
            rBM2=(-2.0d0*log(rand4))**0.5d0
            thetaBM2=2.0d0*pi*rand3
            rand3=rBM2*cos(thetaBM2)
            rand4=rBM2*sin(thetaBM2)
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
        end if

        if (nnflag==1) dxnei=0 !reset displacement counter

        !calculate original forces and additional driving forces
        call force(N,sigma,sigmap6,sigmap12,rwca,Mr,Mt,Mvr,Mvt,&
            &   r1,i,chi,psi,F1,F,Fwca,Fdw,r12,V,dellam,gradV,&
            &   gaussianr,gaussiant,steps,dt,potdw,L,frc,&
            &   rskin,nlist,nnflag)

        !active forces
        Fa(:,1) = vo*cos(theta1)
        Fa(:,2) = vo*sin(theta1)
        
        !keep the dimer particles, 1 and 2, passive
        Fa(1,:) = 0.0
        Fa(2,:) = 0.0

        !update total force
        F1=F1+Fa
        F=F+Fa

        !propagate position
        delr(3:N,:)=dt*F1(3:N,:)/gammab+etavarb*eta(3:N,:)
        delr(1:2,:)=dt*F1(1:2,:)/gammad+etavard*eta(1:2,:)
        r2=r1+delr
        dxnei=dxnei+delr

        !propagate angle
        deltheta =  eta_theta_var*eta_theta
        theta2 = theta1 + deltheta

        !periodic boundaries
        do j=1,N
            r2(j,1)=r2(j,1)-L*nint(r2(j,1)/L)
            r2(j,2)=r2(j,2)-L*nint(r2(j,2)/L)
        end do

        !print configuration
        if (p2>2d6 .and. mod(p2,10000)==0 .and. r12<=1.25d0) then
            do j=1,N
                write(39,*) r1(j,:), theta1(j)
            end do
            p=p+1
        end if

        !print trajectories
        !if (mod(i,10)==0) then
        !    do j=1,N
        !        write(39,*)r2(j,:)
        !    end do
        !    write(40,*)r12,potdw
        !    !write(39,*)r2,F1
        !end if

        !propagate integrals (malliavin weight and the value function integral)
        !dely(:,:)=etavard*((eta(1,1)-eta(2,1))*dellam(:,:,1)+(eta(1,2)-eta(2,2))*dellam(:,:,2))/(2.0d0*Dd*gammad)
        !y=y+dely
        !z=z+dt*gradV

        !propagate estimator and gradient integral
        !xi=-sum(((delr(1:2,:)/dt-F1(1:2,:)/gammad)**2-(delr(1:2,:)/dt-F(1:2,:)/gammad)**2)/(4.0d0*Dd)) &
        !&  -sum(((delr(3:N,:)/dt-F1(3:N,:)/gammab)**2-(delr(3:N,:)/dt-F(3:N,:)/gammab)**2)/(4.0d0*Db))
        !omega(3)=omega(3)+xi
        !delomega=delomega+xi *y -dely/dt *V
        !delL=delL+xi*z-gradV*V

    end do !trajectory propagation ends

    !making estimator and gradients dimensionally correct
    omega(3)=omega(3)*dt
    delomega=delomega*dt
    delL=delL*dt

    !recalculate r12 for evaluating indicator function conditioning
    xr=r2(1,:)-r2(2,:)
    !!periodic boundaries
    xr(1)=xr(1)-nint(xr(1)/L)*L
    xr(2)=xr(2)-nint(xr(2)/L)*L
    r12=dot_product(xr,xr)**0.5

    !impose endpoint biasing
    if (r12>=1.85d0) then !reached target region
        omega(2)=1 !indicator function is 1 here
        delomega=delomega+s*omega(2)*y
        delL=delL+s*omega(2)*z
    end if

    omega(1)=s*omega(2)+omega(3)-s !estimator with full Lagrange multiplier term
    omega(4)=omega(3)*omega(2) !joint expectation of action and indicator

end subroutine traj
!*************************************************************************!
!*************************************************************************!
!*************************************************************************!
!*************************************************************************!
!subroutine for force calculation using gaussian basis
subroutine force(N,sigma,sigmap6,sigmap12,rwca,Mr,Mt,Mvr,Mvt,&
    &   r,stepcount,chi,psi,F1,F,Fwca,Fdw,r12,V,dellam,gradV,&
    &   gaussianr,gaussiant,steps,dt,potdw,L,frc,&
    &   rskin,nlist,nnflag)

    implicit none

    integer, intent(in):: N,Mr,Mt,Mvr,Mvt,stepcount,steps !!assumes Mvr,Mvt=Mr,Mt
    double precision, intent(in):: chi(Mr,Mt), psi(Mvr,Mvt),r(N,2),dt,L,frc,rskin
    integer:: nnflag, nlist(N,N)
    double precision, intent(in):: gaussianr(Mr,2),gaussiant(Mt,2)
    double precision, intent(in):: sigma,sigmap6,sigmap12,rwca
    double precision:: F1(N,2),F(N,2),Fwca(N,2),Fdw(N,2),V, dellam(Mr,Mt,2),gradV(Mvr,Mvt)
    double precision:: potdw,r12
    integer:: i,j,k

    double precision:: xr(2),rcap(2),absr,rprod1,rprod2,Fint(2),fac1,fac2 !temporary variables for force calculation
    double precision, parameter:: ewca=1.0d0,hdw=3.5d0,wdw=0.45d0 !potential energy parameters

    dellam=0.0d0
    gradV=0.d0
    F1=0.0d0
    F=0.0d0
    Fwca=0.d0
    Fdw=0.0d0
    V=0.0d0

    !generic forces common for all particles
    if (nnflag==0) then !do calculation according to neighbor list

        do i=1,N
            !put 1body force here

            !do j=i+1,N
            do k=1,nlist(i,1)
                j=nlist(i,k+1)
                if (j<=i) cycle

                    !!define absr and xr right away and work with those
                    xr=r(i,:)-r(j,:)
                    !!periodic boundaries
                    xr(1)=xr(1)-nint(xr(1)/L)*L
                    xr(2)=xr(2)-nint(xr(2)/L)*L
                    absr=dot_product(xr,xr)**0.5
                    rcap=xr/absr

                    !!WCA evaluation
                    if (absr<rwca) then
                        rprod1=absr*absr
                        rprod1=rprod1*rprod1
                        rprod1=rprod1*rprod1/absr !7th power
                        rprod2=rprod1*rprod1/absr !13th power
                        Fint=ewca*(24.0d0*(2.0d0*sigmap12/rprod2-sigmap6/rprod1)-frc)*rcap
                        Fwca(i,:)=Fwca(i,:)+Fint
                        Fwca(j,:)=Fwca(j,:)-Fint
                    end if
            end do

        end do

    else !make new neighbor list
        nlist(:,1)=0

        do i=1,N
            !put 1body force here

            do j=i+1,N

                !!define absr and xr right away and work with those
                xr=r(i,:)-r(j,:)
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
                    if (absr<rwca) then
                        rprod1=absr*absr
                        rprod1=rprod1*rprod1
                        rprod1=rprod1*rprod1/absr !7th power
                        rprod2=rprod1*rprod1/absr !13th power
                        Fint=ewca*(24.0d0*(2.0d0*sigmap12/rprod2-sigmap6/rprod1)-frc)*rcap
                        Fwca(i,:)=Fwca(i,:)+Fint
                        Fwca(j,:)=Fwca(j,:)-Fint
                    end if
                end if !end calculation for i and j being neighbors

            end do !loop over 2body forces

        end do !loop over 1body forces
        nnflag=0
    end if

    !double well forces for first two particles
    xr=r(1,:)-r(2,:)
    !!periodic boundaries
    xr(1)=xr(1)-nint(xr(1)/L)*L
    xr(2)=xr(2)-nint(xr(2)/L)*L
    absr=dot_product(xr,xr)**0.5
    rcap=xr/absr
    fac1=absr-rwca-wdw
    fac2=1.0d0-(fac1)**2/wdw/wdw
    r12=absr
    potdw=hdw*fac2*fac2
    Fint=4.0d0*hdw/wdw/wdw *fac1*fac2*rcap
    Fdw(1,:)=Fint
    Fdw(2,:)=-Fint

    !account for all original forces
    F=Fwca+Fdw
    F1=F1+F

    !additional driving force and value function in terms of Gaussians
    !do i=1,Mr
    !    do k=1,Mt
    !        fac1=exp(-(absr-gaussianr(i,1))**2/(2.0d0*gaussianr(i,2))&
    !        &   -(stepcount-gaussiant(k,1))**2/(2.0d0*gaussiant(k,2)))
    !        dellam(i,k,:)=fac1*rcap
    !        gradV(i,k)=fac1
    !        F1(1,:)=F1(1,:)+chi(i,k)*dellam(i,k,:)
    !        F1(2,:)=F1(2,:)-chi(i,k)*dellam(i,k,:)
    !        V=V+psi(i,k)*gradV(i,k)
    !    end do
    !end do

end subroutine force
!*************************************************************************!
!*************************************************************************!
!*************************************************************************!
