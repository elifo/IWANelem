module nonlinear3c
  use lu

contains

! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  subroutine get_G_info(modeltype,water,Nspr,Npg,gammaref,straing, GoverG)

  implicit none

  logical, intent(OUT) :: modeltype, water
  integer, intent(OUT) :: Nspr,Npg
  real,    intent(OUT) :: gammaref
  real, dimension(:), allocatable, intent(INOUT) :: straing, GoverG

  integer :: i

   ! Reading the soil G/Gmax data
   open(61, file="GoverGmax.dat", status ="old")
   read(61,*)
   read(61,*)  modeltype
   read(61,*)  water
   read(61,*)  Nspr
   read(61,*)
   if ( modeltype ) then			! If hyperbolic   NL model :
      read(61,*)  gammaref
  else 					              ! If experimental NL model :
      read(61,*) ; read(61,*) ;read(61,*) ; read(61,*) 	
      read(61,*)  Npg

   	  allocate(straing(Npg), GoverG(Npg))
   	  straing = 0.0;      GoverG = 0.0

    	do i=1,Npg
	      read(61,*)  straing(i), GoverG(i)
	      straing(i) = straing(i)/100.0
	    enddo
  endif

   end subroutine get_G_info
! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

   subroutine nonline_curve(Gmax,Nspr,R,gamma,Cninv)

   implicit none

   !--
   real, intent(IN) :: Gmax
   integer, intent(IN) ::  Nspr
   real, dimension(:), intent(INOUT) :: R,gamma, Cninv


   !-Local variables

   integer :: i,j,n			!!!



!   do j=1,Npg
!      x(j) = alog10(straing(j))
!      strain(j) = straing(j)
!   enddo

!   n=Npg

!   x0 = dlog10(dble(strain(1)))
!   xu = dlog10(dble(strain(n)))
!   dx = (xu-x0)/(Nspr-1)

!   do i=1,Nspr
!      gamma(i) = 10.0D0**(x0 + dx*(i-1))
!      call hyper_interp(strain,GoverG,n,sngl(gamma(i)),Gout)
!      G(i) = dble(Gout)
!      R(i) = Gmax* G(i)* gamma(i)
!   enddo

     open(12, file="input/BackboneCurve50Elements.txt", status="old")
!    open(12, file="input/BackboneCurve200Elements.txt", status="old")

   do i=1,Nspr
      read(12,*) gamma(i), R(i)
   enddo
   gamma = gamma/100.0
   R     = R* 1000.0


   call CNinverse(Gmax,Nspr,gamma,R,CNinv)

   return


  end subroutine nonline_curve
! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

   subroutine nonline_curve_hyper(gammaref,Gmax,Npg,Nspr,straing,GoverG,H,R,gamma,CNinv)

   implicit none

   !--
   real, intent(IN) :: Gmax, gammaref
   integer, intent(IN) :: Npg, Nspr
   real, dimension(:), intent(IN) :: straing, GoverG
   real, dimension(:), intent(INOUT) :: H,R,gamma, Cninv



   !-Local variables
   integer :: i,j,n			!!!
   real :: x0,xu,dx
   real :: G(Nspr)


   x0 = -6.0		
   xu = log10(0.1)
   dx = (xu-x0)/(Nspr-1)


   do i=1,Nspr							
      gamma(i) =  10.0**(x0 + dx*(i-1))
      G(i) = (1.0/(1.0+ abs(gamma(i)/gammaref)))

      ! ! Elasticity test
      !G(i) = 1.0

      R(i) = Gmax* G(i)* gamma(i)
   enddo


   call CNinverse(Gmax,Nspr,gamma,R,CNinv)

   return

   end subroutine nonline_curve_hyper
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine hyper_interp(x,y,n,xi,yi)
!	A PROGRAM FOR INTERPOLATION BY ARC HYPERBOLIC
!	SINE METHOD BETWEEN TWO POINTS
!	x--ABSCISSA VALUE OF INPUT POINTS;  y--ORDINATE
! 	VALUE OF INPUT POINTS;  xi--INPUT ABSCISSA VALUE;
!	yi--OUTPUT ORDINATE VALUE

   dimension x(n), y(n)

	 if (x(1)-1.e-10 < 1.e-10) x(1)=1.e-10
   if (xi <= x(1)) then
      yi=y(1)
      return
   else
      j=2
   end if
30 if (xi < x(j)) then
	    ars  = alog(xi*1.e+4+((xi*1.e+4)**2+1.)**0.5)
	    ars1 = alog(x(j-1)*1.e+4+((x(j-1)*1.e+4)**2+1.)**0.5)
	    ars2 = alog(x(j)*1.e+4+((x(j)*1.e+4)**2+1.)**0.5)
	    a    = (y(j)*x(j)*1.e+4-y(j-1)*x(j-1)*1.e+4)/(ars2-ars1)
	    b    = (y(j)*x(j)*1.e+4*ars1-y(j-1)*x(j-1)*1.e+4*ars2) / (ars1-ars2)
      yi   = a/xi/1.e+4*ars + b/xi/1.e+4
      return
   else if (xi == x(j)) then
      yi = y(j)
      return
   else
      j=j+1
   end if
   if (j <= n) then
      goto 30
   else
      yi = y(n)
   end if

end subroutine hyper_interp
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

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
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


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
depsm = (deps(1) + deps(2)+deps(6))/3.0
dsigm = 3.0 * Kmax* depsm


de       = 0.0

!!!
de(1) =  deps(1)- depsm
de(2) =  deps(2)- depsm
de(3) =  deps(3) /2.0
de(4) =  deps(4) /2.0
de(5) =  deps(5) /2.0
de(6) =  deps(6)- depsm


Ed = 0.0
ss(1) = 1.0
ss(2) = 1.0
ss(3) = 2.0
ss(4) = 2.0
ss(5) = 2.0
ss(6) = 1.0

aktif = 0
do j = 1,Nspr
  F2(j) = 0
  dF2(j) = 0
  do k=1,6
      F2(j) = F2(j) + 0.5*ss(k)*(S1(k) - Sa1(j, k))**2
      dF2(j) = dF2(j) + ss(k)*(S1(k) - Sa1(j,k))*de(k)
  enddo

  if ((dF2(j) .GE. 0.0) .and. (F2(j) .ge. R(j)**2.0)) then
    aktif = aktif + 1
    do m=1,6
      do k=1,6
        Ed(m,k) = Ed(m,k)+ CNinv(j) * ss(k)*(S1(m)-Sa1(j,m)) &
                   *(S1(k)-Sa1(j,k))/ (2.0* F2(j))
      enddo
      Sa2(j,m) = S1(m)-R(j)/sqrt(F2(j))*(S1(m)- Sa1(j,m))
    enddo
  else
      exit
  endif
enddo

print*, 'Number of active surfaces   ', aktif

Ed(1,1) = Ed(1,1)+ 0.5/Gmax
Ed(2,2) = Ed(2,2)+ 0.5/Gmax
Ed(3,3) = Ed(3,3)+ 0.5/Gmax
Ed(4,4) = Ed(4,4)+ 0.5/Gmax
Ed(5,5) = Ed(5,5)+ 0.5/Gmax
Ed(6,6) = Ed(6,6)+ 0.5/Gmax

print*, "ITERATION STEP",it, "3ND CDT"

call LUDCMP(Ed, 6, INDX, D, errorflag)
if (errorflag .ne. 0) stop 'not invertible matrix'
call LUBKSB(Ed, 6, INDX, de) ! solve Ed·x = de (de is used as input/ouput)
dS = de

S2 = S1+ dS

! Total Stresses
dsigma(1) = dS(1)+ dsigm
dsigma(2) = dS(2)+ dsigm
dsigma(3) = dS(3)
dsigma(4) = dS(4)
dsigma(5) = dS(5)
dsigma(6) = dS(6)+ dsigm

end subroutine iwan3c

end module nonlinear3c
