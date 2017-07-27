program main3c

  use nonlinear3c ; use watering

  implicit none

  character(20) 	:: fnamef

  logical       	:: modeltype
  logical       	:: water

  integer			  :: Nspr	, aktif
  integer			  :: Npg
  integer			  :: i
  integer			  :: elif 
  integer			  :: j
  integer			  :: i1 
  integer			  :: totsoil
  integer			  :: numsoil 	
  real 			      :: gammaref
  real 			      :: Gmax
  real 			      :: totime
  real 			      :: intime
  real 			      :: rho
  real 			      :: Vs
  real 			      :: Vp
  real 			      :: Kmax
  real 			      :: Emax
  real 			      :: ni
  real 			      :: elo
  real 			      :: eps1(6)
  real 			      :: eps2(6)
  real 			      :: deps(6)
  real 			      :: sigmaC1
  real 			      :: sigold(6)
  real 			      :: k0
  real 			      :: epsold(6)
  real 			      :: sidold(6)
  real 			      :: Told
  real 			      :: Peff0
  real 			      :: Gm0E
  real 			      :: Wn
  real 			      :: Wsold
  real 			      :: S0old
  real 			      :: Sold
  real 			      :: sigmefold
  real 			      :: sigefold(6)
  real 			      :: sigmax
  real 			      :: grefEA
  real 			      :: G0EA
  real 			      :: E0EA
  real 			      :: KbEA
  real 			      :: sigef(6)
  real 			      :: sigmap(3)
  real 			      :: G0EAold
  real 			      :: S
  real 			      :: S0
  real 			      :: T
  real 			      :: Q01
  real 			      :: Q02
  real 			      :: r1
  real 			      :: sigmto
  real 			      :: Ws
  real 			      :: pore
  real 			      :: poreold
  real 			      :: sigmef
  real 			      :: SigOct
  real 			      :: EpsOct
  real 			      :: ello
  real 			      :: mm1
  real 			      :: mm2
  real 			      :: pp1
  real 			      :: pp2
  real 			      :: ww1
  real 			      :: ss1
  real 			      :: cc1
  real 			      :: realn
  real 			      :: realp
  real 			      :: sigpold
  real 			      :: Sigma(6)
  real 			      :: dSigma(6)
  real 			      :: S1(6)
  real 			      :: S2(6)
  real 			      :: dS(6)
  real 			      :: E (6,6)
  real*8 		      :: Ed(6,6)
  real 			      :: dplastic(6)
  real,  allocatable, dimension(:)    ::  straing   
  real,  allocatable, dimension(:)	  ::  GoverG
  real,  allocatable, dimension(:)	  ::  H
  real,  allocatable, dimension(:)	  ::  R
  real,  allocatable, dimension(:)	  ::  gamma
  real,  allocatable, dimension(:)	  ::  CNinv
  real,  allocatable, dimension(:)	  ::  F1
  real,  allocatable, dimension(:)	  ::  F2
  real,  allocatable, dimension(:)	  ::  dF1
  real,  allocatable, dimension(:)	  ::  dF2
  real,  allocatable, dimension(:,:)	  ::  eps
  real,  allocatable, dimension(:,:)	  ::  Sa1
  real,  allocatable, dimension(:,:)	  ::  Sa2
  real,  allocatable, dimension(:,:)	  ::  properties



    ! Reading soil curve ppts
    call get_G_info (modeltype,water,Nspr,Npg,gammaref, straing, GoverG)


    ! Reading soil& saturation parameters  
    open (11, file = 'soiltype', status = 'old')
    read (11,*)   totsoil
    read (11,*)   numsoil
    allocate(properties(12,totsoil))
    do i = 1, totsoil
      read(11,*)
      do j = 1,12
          read(11,*)    properties(j,i)
      end do
    end do

    rho       =   properties(1, numsoil)
    Vs        =   properties(2, numsoil)
    Gmax      =   properties(3, numsoil)
    ni        =   properties(4, numsoil)
    mm1       =   properties(5, numsoil)
    mm2       =   properties(6, numsoil)
    pp1       =   properties(7, numsoil)
    pp2       =   properties(8, numsoil)
    ss1       =   properties(9, numsoil)
    ww1       =   properties(10,numsoil)
    cc1       =   properties(11,numsoil)
    gammaref  =   properties(12,numsoil)

    Emax      =   2.0* Gmax* (1.0+ ni)                      
    Kmax      =   Emax/ (3.0* (1.0- 2.0* ni))               

    print *, rho
    print *, Vs
    print *, Gmax
    print *, gammaref
    print *, ni
    print *, Emax
    print *, Emax* ni/(1.0+ ni)/(1.0- 2.0*ni)
    print *, Kmax



!!! DIEGO- IWAN GROUP MODIFICATION
Nspr = 51
!Nspr  = 201


    allocate( R     (Nspr   ))
    allocate( gamma (Nspr   ))
    allocate( CNinv (Nspr- 1))
    R 			=	0.0
    CNinv		=	0.0

  
  if (.NOT. water)    then
      if (.NOT. modeltype) then
          call nonline_curve(Gmax,Nspr,R,gamma,CNinv )
          print*, "EXPERIMENTAL MODEL"
          open(11, file="outputfiles/SOILEXP",status="unknown")
          do i=1, Nspr
              write(11,*)   gamma(i)*100.0, R(i)/1.e3
            enddo
        else
          call nonline_curve_hyper(gammaref,Gmax,Npg,Nspr,straing,GoverG,H,R,&
                              gamma,CNinv )  
            print*, "HYPERBOLIC MODEL"
            open(11, file="outputfiles/SOILHYPER",status="unknown")
            do i=1, Nspr
              write(11,*)   gamma(i)*100.0, R(i)/1.e3
            enddo
        endif
        close(11)
   endif


	! SATURATED SOIL - Initial conditions 1 ! [sproperties] 
	!if (water) then
	sigmap   	= 0.0 
  	Ws       	= 0.0
  	sigmto   	= 0.0
	T        	= 0.0 
	S0       	= 0.0
  	S        	= 0.0 
	sigmef   	= 0.0
    sigef    	= 0.0
  	pore     	= 0.0
  	poreold  	= 0.0
  	SigOct   	= 0.0
  	EpsOct   	= 0.0
  	epsold   	= 0.0
  	sigold   	= 0.0
  	S0old    	= 1.0
  	Sold     	= 1.0
  	sigmaC1  	= 0.0
  	k0       	= 1.0
	!endif

	call openfile(water)
	call Initialcdts (water,Gm0E,Peff0,sigold,eps1,Wn,Wsold,Told,sidold,S0old,k0,mm1,mm2,rho,Gmax,sigmaC1,pp1,pp2,ww1,ss1)

	! ALL SOILS - Iwan initial conditions ! 
	allocate(F1  (Nspr)	) 
	allocate(F2  (Nspr)	)
	allocate(dF1 (Nspr)	)
	allocate(dF2 (Nspr)	)
	allocate(Sa1 (Nspr,6)	)
	allocate(Sa2 (Nspr,6)	)
	F1     = 0.0
	F2     = 0.0
	dF1    = 0.0 
	dF2    = 0.0 
	dsigma = 0.0 
	S2     = 0.0 
	dS     = 0.0
	SA1    = 0.0
	Sa2    = 0.0
	E      = 0.0
	Ed     = 0.0
	eps2   = 0.0
	eps1   = 0.0
	deps   = 0.0

    dplastic = 0.0

	! Additional conditions 
	G0EAold    	= Gm0E
	sigmefold  	= Peff0
	sigefold   	= sigold 
	S1         	= max(S1, 0.001)




!!!
aktif = 0



! ! DIEGO- IWAN GROUP MODIFICATION
!    open(11, file= "input/strain.dat",status="old")


   totime = 3.0
   intime = 0.001


	allocate( eps(int(totime/intime),3) )
	eps    	= 0.0
	eps1   	= 0.0 
	eps2   	= 0.0
	deps   	= 0.0

	sigma 	= 0.0

   
	! STARTING UP
	do i = 1,int(totime/intime)  
   	i1 = i   
   	
    ! Updating the parameters
    if (water)   &
    call update(water,cc1,ni,eps1,S0old,Gm0E,sigmax,Peff0,mm1,mm2,Sold,gammaref,grefEA,G0EA,E0EA,KbEA,nspr,gamma,R,CNinv,i1)

   
    ! Input strain definition
    eps(i,1) = sin((i- 1)*intime*4.0)*0.075*0.015* exp((i- 1)*intime/2.0)*5.0
!     read(11,*)  eps(i,1)


   	eps2(4) = eps(i,1)* 1.0  		!xz					
   	eps2(5) = eps(i,1)* 0.0			!yz		
   	eps2(6) = eps(i,1)* 0.0 		!zz	


!    	eps2(4) = eps(i,1)* 1.0  		!xz					
!    	eps2(5) = eps(i,1)* 2.			!yz		
!    	eps2(6) = eps(i,1)* 1.5 		!zz		
   	deps    = eps2- eps1


 		E      	= 0.0
 		dsigma 	= 0.0


    ! Stress increment computation
		call iwan3c(i,Nspr,Emax,Gmax,Kmax,E,deps,dsigma,dS,F1,F2,S1,Sa1,R,Ed,CNinv,Sa2,dF1,dF2,S2,elif,ni,dplastic)

 		S2  = S1+dS
 		S1  = S2

 		Sa1 = Sa2
 		dF1 = dF2

 		sigma = sigma + dsigma



    ! Writing into files

    ! Strain vs Time
    write(21,*)  (i-1)*intime, eps2(4)  
    write(22,*)  (i-1)*intime, eps2(5)  
    write(23,*)  (i-1)*intime, eps2(6)

    ! Stress vs Strain
    write(31,*) eps2(4)*100.0 ,  (sigma(4))/1.e3 
    write(32,*) eps2(5)*100.0 ,  (sigma(5))/1.e3
    write(33,*) eps2(6)*100.0 ,  (sigma(6))/1.e3      


61 continue


    eps1	 = eps2
  
    dsigma = 0.0
    dS     = 0.0
    deps   = 0.0

enddo   ! TIME'S UP



print*, "PROGRAM BITTI"
end program main3c
