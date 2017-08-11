module nonlinear3c
  
  contains

! ----------------------------------------------------------

  subroutine get_G_info(water,Nspr,gammaref)

  implicit none

  logical, intent(OUT) ::  water
  integer, intent(OUT) :: Nspr
  real, intent(OUT) :: gammaref

  integer :: i

  ! Reading the soil G/Gmax data
  open(61, file="GoverGmax.dat", status ="old")
  read(61,*)
  read(61,*)  water
  read(61,*)  Nspr
  read(61,*)
  read(61,*)  gammaref
  !
  end subroutine get_G_info
! 
! ----------------------------------------------------------

   subroutine nonline_curve_hyper(gammaref,Gmax,Nspr,H,R,gamma,CNinv)

   implicit none

   real, intent(IN) :: Gmax, gammaref
   integer, intent(IN) :: Nspr
   real, dimension(:), intent(INOUT) :: H,R,gamma, Cninv

   integer :: i,j,n     !!!
   real :: x0,xu,dx
   real :: G(Nspr)


   x0 = -6.0    
   xu = log10(0.1)
   dx = (xu-x0)/(Nspr-1)


   do i=1,Nspr              
      gamma(i) =  10.0**(x0 + dx*(i-1))
      G(i) = (1.0/(1.0+ abs(gamma(i)/gammaref)))
!      Elasticity test
      !G(i) = 1.0
      R(i) = Gmax* G(i)* gamma(i)
   enddo

   call CNinverse(Gmax,Nspr,gamma,R,CNinv)

   return

   end subroutine nonline_curve_hyper
!
! ----------------------------------------------------------

   subroutine CNinverse(Gmax,nspr,gamma,R,CNinv)

   implicit none

   real, intent(IN)     :: Gmax
   integer, intent(IN)  :: nspr
   real, intent(IN)   :: gamma(nspr)
   real, intent(IN)     :: R(nspr)

   real :: CNinv(nspr-1), sumdummy
   integer :: i

   sumdummy = 0.0
   do i= 1,nspr-1
      CNinv(i) = ((gamma(i+1)/2.0-gamma(i)/2.0)/(R(i+1)-R(i))) - 0.5/Gmax - sumdummy
      sumdummy = sumdummy+ CNinv(i)
   enddo

   end subroutine Cninverse
!
! ----------------------------------------------------------
!

  subroutine iwan3c(it,Nspr,Emax,Gmax,Kmax,E,deps,dsigma,dS,F1,F2,S1,Sa1,R,Ed,CNinv,Sa2,dF1,dF2,S2,elif,ni,dplastic)

  implicit none

  integer,  intent(IN)                    ::  it, Nspr
  real,     intent(INOUT)                 ::  Emax, Gmax, Kmax
  real,     dimension(:),   intent(INOUT) ::  deps, dsigma, dS,F1,F2,S1,R,CNinv,dF1,dF2,S2
  real,     dimension(:,:), intent(INOUT) ::  E,Sa1,Sa2
  real*8,   dimension(:,:), intent(INOUT) ::  Ed
  integer,  intent(OUT)                   ::  elif
  real,     intent(IN)      ::  ni
  real,     intent(INOUT)       ::  dplastic(6)

  real*8      :: dsigm, ss(6), Ln(Nspr), depsm, de(6), Kbm(6,6), Sm(6,6),Af(6,6),m1(6,6),m2(6,6),m3(6,6)
  real*8    :: Esd(6,6)
  integer   :: INDX(6)
  integer   :: D, j, m, k, errorflag, aktif
  real      :: lambda



  lambda  = Emax* ni/(1.0+ ni)/(1.0- 2.0*ni)

  dsigm = (dsigma(1)+ dsigma(2)+ dsigma(6))/3.0
  depsm = dsigm/(Kmax)    


  dS    = 0.0
  dS(1) = dsigma(1)-dsigm
  dS(2) = dsigma(2)-dsigm
  dS(3) = dsigma(3)
  dS(4) = dsigma(4)
  dS(5) = dsigma(5)
  dS(6) = dsigma(6)-dsigm

  Ed    = 0.0

  ss(1) = 1.0
  ss(2) = 1.0
  ss(3) = 2.0
  ss(4) = 2.0
  ss(5) = 2.0
  ss(6) = 1.0


  F2    = 0.0
  dF2   = 0.0
  aktif = 0


  do j = 1,Nspr

    do k=1,6
        F2(j)  = F2(j) + 0.5*ss(k)*(S1(k) - Sa1(j, k))**2
        dF2(j) = dF2(j) + ss(k)*(S1(k) - Sa1(j,k))*dS(k)     
    enddo

    if ((dF2(j) .GE. 0.0) .and. (F2(j) .ge. R(j)**2.0)) then
      aktif = aktif + 1
      do m=1,6
        do k=1,6

          !!!
          F2(j) = max(1.e-5, F2(j))

          Ed(m,k) = Ed(m,k)+ CNinv(j) * ss(k)*(S1(m)-Sa1(j,m)) &
                   *(S1(k)-Sa1(j,k))/ (2.0* F2(j))
      
        enddo
        Sa2(j,m) = S1(m)-R(j)/sqrt(F2(j))*(S1(m)- Sa1(j,m))
      enddo
    else
        exit
    endif
  enddo


  ! print*, 'Number of active surfaces   ', aktif
  Ed(1,1) = Ed(1,1)+ 0.5/Gmax
  Ed(2,2) = Ed(2,2)+ 0.5/Gmax
  Ed(3,3) = Ed(3,3)+ 0.5/Gmax
  Ed(4,4) = Ed(4,4)+ 0.5/Gmax
  Ed(5,5) = Ed(5,5)+ 0.5/Gmax
  Ed(6,6) = Ed(6,6)+ 0.5/Gmax


  ! print*, "ITERATION STEP",it, "3ND CDT"
  ! No inversion
  de = MATMUL(Ed, dS)


  S2 = S1+ dS

  ! Total strains
  deps(1) = de(1)+ depsm
  deps(2) = de(2)+ depsm
  deps(3) = de(3)*2.0
  deps(4) = de(4)*2.0
  deps(5) = de(5)*2.0
  deps(6) = de(6)+ depsm

  ! Plastic strain increment
  dplastic    = 0.0
  dplastic(4) = deps(4)- dsigma(4)/Gmax 
  dplastic(5) = deps(5)- dsigma(5)/Gmax 
  dplastic(6) = deps(6)- dsigma(6)/(lambda+ 2.0*Gmax) 

  end subroutine iwan3c


end module nonlinear3c