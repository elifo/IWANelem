module watering

contains

! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine openfile(water) 

implicit none

logical, intent(IN), target :: water

if (water) then
open(11, file = "outputfiles/LINEM1" , status = "unknown")
open(12, file = "outputfiles/LINEM2" , status = "unknown")
open(13, file = "outputfiles/LINEM3" , status = "unknown")
open(14, file = "outputfiles/LINEM1-", status = "unknown")
open(15, file = "outputfiles/LINEM2-", status = "unknown")

open(51, file = "outputfiles/poretime",            status="unknown")
open(52, file = "outputfiles/deviatoriceffective", status="unknown")
open(53, file = "outputfiles/deviatoriceffective2", status="unknown")
open(54, file = "outputfiles/principals",status = "unknown")

open(94, file = "outputfiles/Gmodulus",      status = "unknown")
open(95, file = "outputfiles/gref",  status = "unknown")
open(96, file = "outputfiles/Sw", status = "unknown")
open(97, file = "outputfiles/S0",     status = "unknown")
open(98, file = "outputfiles/S",      status = "unknown")

open(99, file = "outputfiles/KONTROL",status = "unknown")

endif


open(21, file = "outputfiles/epsxz", status = "unknown")
open(22, file = "outputfiles/epsyz", status = "unknown")
open(23, file = "outputfiles/epszz", status = "unknown")

open(31, file = "outputfiles/stressstrainxz", status = "unknown")
open(32, file = "outputfiles/stressstrainyz", status = "unknown")
open(33, file = "outputfiles/stressstrainzz", status = "unknown")

open(34, file = "outputfiles/stresstimexz", status = "unknown")
open(35, file = "outputfiles/stresstimeyz", status = "unknown")
open(36, file = "outputfiles/stresstimezz", status = "unknown")

open(41, file = "outputfiles/OCTOstressstrain", status = "unknown")
open(42, file = "outputfiles/OCTOstresstime"  , status = "unknown")



end subroutine openfile
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine Initialcdts (water,Gm0E,Peff0,sigold,epsold,Wn,Wsold,Told,sidold,S0old,k0, mm1,mm2, rho, Gmax,sigmaC1,pp1,pp2,ww1,ss1)

implicit none

logical,  intent(IN)    :: water
real,     intent(INOUT) :: Gm0E,Peff0,sigold(6),epsold(6),Wn,Wsold,Told,sidold(6),S0old,sigmaC1
real,     intent(IN)    :: k0, mm1,mm2, rho, Gmax,pp1,pp2,ww1,ss1

real                     :: r0,a,b,c,delta, x1,x2, min,max, w0, TG


if (water) then
  sigold(1) =  k0* 9.8* (rho- 1000.0)* 10.0           ! Depth = 10m assumed
  sigold(2) =  k0* 9.8* (rho- 1000.0)* 10.0   
  sigold(6) =  9.8* (rho- 1000.0)* 10.0 

! Anisotropy to have the same Peff0
!  sigold(1) =  k0* 147000.0           
!  sigold(2) =  k0* 147000.0     
!  sigold(6) =  147000.0   

else
  sigold(1) =  k0* 9.8* rho*10.0             
  sigold(2) =  k0* 9.8* rho*10.0    
  sigold(6) =  9.8* rho* 10.0 
       
  Gm0E      =  Gmax
  Peff0     =  (sigold(1)+sigold(2)+sigold(6))/3.0
       
  sidold(1) =  sigold(1)- Peff0
  sidold(2) =  sigold(2)- Peff0
  sidold(6) =  sigold(6)- Peff0

  ! Normalized initial stress
  sigmaC1 = Peff0*(1.0+ 2.0*k0)/3.0

  go to 12
endif

epsold = 0.0
sidold = 0.0

! Normalized initial stress
Told      = abs(sigold(6)-sigold(2))/2.0
Peff0     = (sigold(1)+sigold(2)+sigold(6))/3.0

sigmaC1   = Peff0* (1.0+ 2.0*k0)/3.0
Gm0E      = Gmax* ((abs(Peff0/sigmaC1))**0.5)      

Wn        = ((abs(Peff0*mm1))**2.0)/(2.0*Gm0E)

sidold(1) = sigold(1)- Peff0
sidold(2) = sigold(2)- Peff0
sidold(6) = sigold(6)- Peff0

! Initial shear work
r0        = Told/(Peff0)
a         = (1.0- 0.33*mm2/mm1)**2.0-(0.67*mm2/mm1)**2.0-(0.33*mm2/mm1)**2.0 
b         = -2.0* (1.0- 0.33*mm2/mm1- 0.67*mm2*r0/(mm1*mm1))
c         = 1.0- (r0/mm1)**2.0
delta     = b*b - 4.0*a*c

if (delta .lt. 0.0) then
  S0old = 1.0
else 
  if (a .eq. 0.0) then
    write(*,*)  "ERROR : Initial shear work - parameter a is zero"
    S0old = 1.0
    go to 11
  endif

  x1    = (-b - sqrt(delta))/(2.0*a)
  x2    = (-b + sqrt(delta))/(2.0*a)
  
  if(x1 < x2) then                      
    min = x1
    max = x2
  else 
    min = x2
    max = x1
  end if                             

  if (0.0 < min) then
    if (r0 < (min*0.67*mm2) ) then
      S0old = 1.0
    else
      if (r0 < (max*0.67*mm2) ) then
        S0old = min
      else
        S0old = min
      endif
    endif
  else 
    if ( 0.0 < max) then
      if ( r0 < (max*0.67*mm2) ) then
        S0old = 1.0
      else
        S0old = max
      endif
    else
      S0old = 1.0
    endif
  endif

  if (S0old > 1.0)      S0old = 1.0
endif 
   
11 continue

if ( S0old < 0.4 ) then
  TG    =   (abs((S0old-ss1)/ (0.4-ss1)))**(1.0/pp2)
  w0    =   ww1/TG
else
  TG = (abs((1.0- S0old)/0.6))**(1.0/pp1)
  w0 = ww1* TG
endif

Wsold = w0 * Wn

12 continue

end subroutine Initialcdts
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine update(water,cc1,ni,epsold,S0old,Gm0E,sigmax,Peff0,mm1,mm2,Sold,gammaref,grefEA,G0EA,E0EA,KbEA,nspr,strEA,&
bbcEA,CNinv,i1)

implicit none

logical, intent(IN)       :: water
real,    intent(IN)       :: cc1
integer, intent(IN)       :: nspr
integer, intent(IN)       :: i1
real,    intent(IN)       :: S0old
real,    intent(IN)       :: Sold
real,    intent(IN)       :: Gm0E
real,    intent(IN)       :: mm1
real,    intent(IN)       :: mm2
real,    intent(IN)       :: gammaref
real,    intent(IN)       :: Peff0
real,    intent(IN)       :: ni
real,    intent(IN)       :: epsold(6)
real,    intent(INOUT)    :: sigmax
real,    intent(INOUT)    :: grefEA
real,    intent(INOUT)    :: G0EA
real,    intent(INOUT)    :: KbEA
real,    intent(INOUT)    :: E0EA
real,    intent(INOUT)    :: strEA(nspr)
real,    intent(INOUT)    :: bbcEA(nspr)
real,    intent(INOUT)    :: CNinv(nspr- 1)

real                      :: delta
real                      :: maxeps
real                      :: sumdummy
real                      :: x0
real                      :: xu
real                      :: dx
integer                   :: k

! [Iai et al. 1990 - Eqns 86-91]
delta     =  0.0
if  (S0old > 0.4) then
  ! Shear strength
  sigmax = Peff0* mm1* Sold         
  ! Actual strain
  grefEA = gammaref
  ! Max shear modulus
  !G0EA   = Gm0E*Sold
  G0EA   = sigmax/grefEA
else
  delta  = (mm1-mm2)*(0.4 - S0old)*Peff0      
  sigmax = Peff0*mm1*Sold + delta             
  grefEA = gammaref/ (S0old/0.4)
  G0EA   = sigmax/grefEA  
endif

! New definitions
E0EA    =   2.0 * G0EA* (1.0+ ni)
KbEA    =   E0EA/3.0/(1.0- 2.0*ni)

! Cn coefficients
x0      = -6.0
xu      = log10(0.10/2.0)         ! Elif's approach
!xu      = log10(grefEA)       ! Viet's approach
dx      = (xu-x0)/(nspr-1)

strEA(1) = 0.0
bbcEA(1) = 0.0 

do k = 2, nspr
  strEA(k) = 10.0**(x0+(k- 1)*dx)      
  bbcEA(k) = G0EA* 2.0* strEA(k)/(1.0+ abs(2.0* strEA(k)/grefEA))       ! here gamma not epsilon
enddo

sumdummy =  0.0
do k = 1,nspr- 1
  CNinv(k) = ((strEA(k+ 1)- strEA(k))/ (bbcEA(k+ 1)- bbcEA(k)))- 0.5/ G0EA- sumdummy     ! here epsilon needed
  sumdummy = sumdummy+ CNinv(k)
enddo

62 continue

end subroutine update
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine Saturation (gammaref,sigma,sigmaini,sigmap,eps,epsold,poreold,Wsold,Sold,S0old,G0EAold,G0EA,Gm0E,Peff0,sigmefold,&
Told,T,r,sigmto,Ws,mm1,mm2,cc1,Wn,ww1,pp1,pp2,ss1,S0,S,sigmef,pore,sigef, realn, realp,i1,dplastic)

implicit none

integer,   intent(IN)    :: i1
real,      intent(INOUT)    :: sigma(6), sigmaini(6)
real,      intent(IN)    :: gammaref
real,      intent(IN)    :: eps(6)
real,      intent(IN)    :: epsold(6)
real,      intent(IN)    :: poreold
real,      intent(IN)    :: mm1
real,      intent(IN)    :: mm2
real,      intent(IN)    :: cc1
real,      intent(IN)    :: ww1
real,      intent(IN)    :: Wn
real,      intent(IN)    :: pp1
real,      intent(IN)    :: pp2
real,      intent(IN)    :: ss1
real,      intent(INOUT) :: sigmap(3)
real,      intent(INOUT) :: realn
real,      intent(INOUT) :: realp
real,      intent(INOUT) :: Wsold
real,      intent(INOUT) :: S0
real,      intent(INOUT) :: Sold
real,      intent(INOUT) :: S0old
real,      intent(INOUT) :: G0EAold
real,      intent(INOUT) :: G0EA
real,      intent(INOUT) :: Gm0E
real,      intent(INOUT) :: Peff0
real,      intent(INOUT) :: sigmefold
real,      intent(INOUT) :: Told
real,      intent(INOUT) :: T
real,      intent(INOUT) :: sigmto
real,      intent(INOUT) :: Ws
real,      intent(INOUT) :: S
real,      intent(INOUT) :: r
real,      intent(INOUT) :: sigmef
real,      intent(INOUT) :: pore 
real,      intent(INOUT) :: sigef(6)
real,      intent(IN)    :: dplastic(6)


real       :: y(3)
real       :: a
real       :: b
real       :: c
real       :: d
real       :: delta
real       :: gamma
real       :: TG
real       :: x1
real       :: x2
real       :: x3
real       :: x4
real       :: x5
real       :: x
real*8     :: e(3,3)
real*8     :: q(3,3)
real*8     :: f(5)
real*8     :: ch
integer    :: i
integer    :: j
integer    :: l

real       :: yy1(6)
real       :: depsv
real       :: CoR
real       :: dWse
real       :: dWst
real       :: dWs

real       :: w
real       :: rr2
real       :: mm3
real       :: rr3
real       :: S2
real       :: yy(6) 
logical    :: eliflogic

! ---------------------------------------------- !
! PRINCIPAL STRESS !
! ---------------------------------------------- !

!sigma = sigma- sigmaini                        


if ( sigma(3)== 0.0 .AND. sigma(4)== 0.0 .AND. sigma(5)== 0.0 ) then
  y(1) = sigma(1)
  y(2) = sigma(2)
  y(3) = sigma(6)
  go to 20
endif

! Coefficients of the eqn of 3rd power
a   = -1.0
b   = sigma(1)+ sigma(2)+ sigma(6)
 c  = sigma(3)**2.0 + sigma(4)**2.0+ sigma(5)**2.0 &
    - sigma(1)*sigma(2)&
    - sigma(2)*sigma(6)&
    - sigma(6)*sigma(1)
d   = sigma(1)*sigma(2)*sigma(6) &
  + 2.0*sigma(3)*sigma(4)*sigma(5) &
  - sigma(1)*sigma(4)**2.0 &
  - sigma(2)*sigma(5)**2.0 &
  - sigma(6)*sigma(3)**2.0

delta = b**2.0- 3.0*a*c

if (delta .ne. 0.0) then
  gamma = (9.0*a*b*c- 2.0*b**3.0- 27.0*(a**2.0)*d)/(2.0*sqrt(abs(delta**3.0)))
endif

if (delta > 0.0) then
  if (abs(gamma) .le. 1.0) then
    x1 = (2.0* sqrt(delta)*(cos(acos(gamma)/3.0))-b)/(3.0*a)
    x2 = (2.0* sqrt(delta)* cos(acos(gamma)/3.0 - 2.0* 4.0* atan(1.0)/3.0)-b)/(3.0*a)
    x3 = (2.0* sqrt(delta)* cos(acos(gamma)/3.0 + 2.0* 4.0* atan(1.0)/3.0)-b)/(3.0*a)
  else
    !go to 10
    go to 61
    eliflogic = .False. 
  endif
else if (delta .ge. 0.0) then
!  x4 = b**3.0 - 27.0*d*a**2.0
!  x5 = abs(b**3.0- 27.0*d*a**2.0)**(1.0/3.0)
!  x  = (-b+ SIGN(x5,x4))/(3.0*a) 
!  x1 = x
!  x2 = x
!  x3 = x
!else
!   go to 10
   go to 61
   eliflogic = .False.
endif

y(1) = x1
y(2) = x2
y(3) = x3

go to 20

!10   continue

!e(1,1)=sigma(1)
!e(2,2)=sigma(2)
!e(1,2)=sigma(3)
!e(2,1)=sigma(3)
!e(1,3)=sigma(5)
!e(3,1)=sigma(5)
!e(2,3)=sigma(4)
!e(3,2)=sigma(4)
!e(3,3)=sigma(6)

!call Rjacob(e,f,q,ch)
!if(ch .gt. 5.d-15) then
!  print *, "CH", ch
!end if
                          
!y(1)  = f(1)
!y(2)  = f(2)
!y(3)  = f(3)

20    continue

write(99,*)  delta, sigma(:)

do i = 1,2     
  l = i
  do j = i,3
    if (y(l) .gt. y(j) )   l = j
  enddo

  if(l.NE.i) then
    TG    = y(i)
    y(i)  = y(l)
    y(l)  = TG
  end if
end do 

do i = 1,3
  sigmap(i) = y(4-i)
end do

61 continue

! ---------------------------------------------- !
! SHEAR WORK !
! ---------------------------------------------- !
!sigma = sigma + sigmaini                        ! Total stress matrix including initial stress


yy1(1)  = 1.0
yy1(2)  = 1.0
yy1(6)  = 1.0
yy1(3)  = 0.0
yy1(4)  = 0.0
yy1(5)  = 0.0

! Volumetric strain increment
!depsv   = eps(1)- epsold(1) &
!            + eps(2)- epsold(2) &
!            + eps(6)- epsold(6) 

! Hydrostatic  stress
!sigmto  = (sigmap(1)+ sigmap(2)+ sigmap(3))/ 3.0              ! Principal axes rotation                  


! Shear work increment [Viet's PhD Eqn. 2.3.21]
!dWst    = 0.0
!do i = 1,6
!  !dWst  = dWst + sigma(i)* (eps(i)- epsold(i))* z(i) 
!  dWst  = dWst + (sigma(i)- poreold*yy1(i))* (eps(i)- epsold(i)) !* z(i)    
!enddo

!write(70,*) (i1- 1)*0.0001, depsv
!write(77,*) (i1- 1)*0.0001, (sigmto)* depsv/ 3.0   !, sigmto, sigmap


!dWst = abs(dWst - ((sigmto- poreold)* depsv/ 3.0))      
!dWst = abs(dWst - ((sigmto- poreold)* depsv))      

!write(71,*) (i1- 1)*0.0001, dWst


!T = (sigmap(1)-sigmap(3))/2.0


! Deviatoric stress
if (eliflogic) then
  T = (sigmap(1)-sigmap(3))/2.0
else
  T = Told
endif


!if (realp .lt. 1.0E-15 .or. realn .lt. 1.0E-15) then
!  dWse  = 0.0
!else
!  !dWse  =  abs(T* (T/G0EA - Told/G0EAOld))
!  dWse  = abs(T* (T/realn - (Told/realp) )) 
!endif

!dWs  =  dWst- cc1* dWse

!if (dWs .le. 0.0)  dWs = 0.0



! After Iai's mail
dWs = 0.0

do j = 1,6
  dWs = dWs+ (sigma(j)- yy1(j)*(sigma(1)+sigma(2)+sigma(6))/3.0)* dplastic(j)
enddo



r = Told/Peff0                                         

if (Sold  .ge. 0.4) then                    
    CoR = 1.0                   
  if (r/S0old .le.  0.67*mm2) then         
    CoR = 1.0
  else
    CoR = (mm1- r/Sold)/(mm1- 0.67*mm2)
  endif
else
  if (r .le. 0.4*0.67*mm2) then                         
    CoR = 1.0
  else
    CoR = (0.4*mm1- r)/(0.4*(mm1- 0.67*mm2))
  endif
endif

if(dWs .gt. 0.0)        dWs     = CoR* dWs
Ws      = Wsold+ dWs


! ---------------------------------------------- !
! S0  !
! ---------------------------------------------- ! 
w       = Ws/Wn

write(96,*)   w/ww1, S0old

write(74,*) (i1- 1)*0.0001, w



if (w .le. 0.0) then
  S0 = 1.0
elseif (w .le. ww1) then
  S0 = 1.0- 0.6* ((w/ww1)**pp1)
elseif (w .gt. ww1) then
  S0 = (0.4- ss1)* ((ww1/w)**pp2)+ ss1
endif


! ---------------------------------------------- !
! S   !
! ---------------------------------------------- !
r     = T/Peff0                                  

rr2   = mm2 * S0
mm3   = mm2 * 0.67 
rr3   = mm3 * S0
S2    = S0- (rr2-rr3)/mm1

if ( r .le. rr3) then
  S = S0
else
  S = S2 + sqrt((S0-S2)**2.0 + ((r-rr3)/mm1)**2.0)
endif




! ---------------------------------------------- !
! EFFECTIVE STRESS   !
! ---------------------------------------------- !

   yy(1)=1.0
   yy(2)=1.0
   yy(6)=1.0

   yy(3)=0.0
   yy(4)=0.0
   yy(5)=0.0

pore   = Peff0 * (1.0- S) 
sigmef = S* Peff0

do i=1,6
  sigef(i) = sigma(i) - pore*yy(i)
enddo  

! Shear work differential parameters
!realp     = Gm0E* (Sold**0.5)
!realn     = Gm0E* (S**0.5)

end subroutine Saturation
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine OctoParam(sigef,eps,SigOct,EpsOct)

implicit none

real, intent(IN) :: sigef(6), eps(6)
real, intent(INOUT) :: SigOct, EpsOct


SigOct  = (1.0/3.0)* sqrt((sigef(1)-sigef(2))**2.0  &
                       + (sigef(2)-sigef(6))**2.0  &
                       + (sigef(6)-sigef(1))**2.0  &
                       + 6.0* ((sigef(3))**2.0 + (sigef(4))**2.0 +(sigef(5))**2.0))


EpsOct  = (2.0/3.0)* sqrt(2.0*(eps(6))**2.0 + 6.0*((eps(4))**2.0+ (eps(5))**2.0))



end subroutine OctoParam
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      Subroutine Rjacob(a,d,q,ch)

!     Calculer les valeurs principales par la méthode Jacobi.

      implicit none

      integer  maxstep
      real*8 a(3,3), q(3,3), d(5)
      
      real*8 ome,s,c,r,t,ch,x,y,del,u,v
!       parameter (maxstep=30)
      integer count,i,j,id,jc,n

      n       = 3
      maxstep = 10000
      ch      = 0.d0

      do i=1,3
         d(i)=a(i,i)
         q(i,i)=1.d0
         do j=1,3
            if(j.LT.i) ch=ch+a(i,j)*a(i,j)
            if(j.NE.i) q(i,j)=0.d0
         end do 
      end do

      count = 0.d0
      del = dsqrt(ch)/2.D0/3.D0

2     continue

      do i=1,maxstep
         del = del/(10.d0**i)
         do id=1,3
            do jc=1,id-1
               if(abs(a(id,jc)).LE.1.d-20) then
                  a(id,jc)=0.d0
               else 
                  x=d(jc)-d(id)
                  if(x==0.d0) then
                     r=1.D0
                  else
                     if(abs(a(id,jc)).LT.1.d-12) then
                        r=a(id,jc)/x
                     else
                        ome=0.5D0*x/a(id,jc)
                        r=1.d0/(abs(ome)+dsqrt(1.d0+ome*ome))
                        if(ome.LT.0.d0) r=-r
                     endif
                  endif
                  c=1.d0/dsqrt(1.d0+r*r)
                  s=c*r
                  t=s/(1.d0+c)
                  y=r*a(id,jc)
                  d(id)=d(id)-y
                  d(jc)=d(jc)+y
                  ch=ch - a(id,jc)*a(id,jc)
                  a(id,jc)=0.d0

                  do j=1,3
                     u=q(j,id)
                     v=q(j,jc)
                     q(j,id)=u-s*(v+t*u)
                     q(j,jc)=v+s*(u-t*v)
                  end do

                  if(ch.gt.1.d-15)then
                     do j=1,jc-1
                        u=a(id,j)
                        v=a(jc,j)
                        a(id,j)=u-s*(v+t*u)
                        a(jc,j)=v+s*(u-t*v)
                     end do
 
                     do j=jc+1,id
                        u=a(id,j)
                        v=a(j,jc)
                        a(id,j)=u-s*(v+t*u)
                        a(j,jc)=v+s*(u-t*v)
                     end do

                     do j=id+1,n
                        u=a(j,id)
                        a(j,id)=u-s*(v+t*u)
                        a(j,jc)=v+s*(u-t*v)
                     end do
 
                     count=count+1

                  else

                     d(4)= dble(i)
                     d(5)= dble(count)
                     return

                  endif
               endif
            end do
         end do
         if(i==maxstep) then
            if(ch.gt.1.d-15) then
               goto 1
            end if
          end if
      enddo

1     Continue

      maxstep = maxstep + 500
        
      if(maxstep>25000) return

      goto 2
!      pause 'not converged, choose other init. guess'
      
      return

      end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
end module watering
