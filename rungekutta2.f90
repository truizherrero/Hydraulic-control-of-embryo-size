module constants
implicit none
save
real(8),parameter::mu=6e-9, R0=20d0,h0=3d0,b0=2.5d0,hc=10d0  
real(8)::Jv,Jc,sigma1,sigma2,b,a2,a3,dt,h22,EE,Kcell,Kjunct
real(8),parameter::pi=3.14159274101257
real(8),parameter::tmax=50000
end module constants

program rungekuttaSWEEP  
use constants
implicit none

real(8),allocatable::Xinit(:)

real(8)::q1,q2,q3,q4,w1,w2,w3,w4,e1,e2,e3,e4,z1,z2,z3,z4,eps,deps
real(8)::x1,x2,x3,x4,t,sigmat,sigmac,pmax,hmax,rmax,sigmax,trupt,epsmax
integer::i,nwrite



open(1,file='parameter2.dat')

read(1,*)Jv,Jc,EE,sigma2

i=0
dt=1d-7
nwrite=1

Kcell=0.1
Kjunct=0


sigma1=sigma2-0.5*sigma2  

 
a2=1d0
a3=a2

allocate(Xinit(6))

!initial conditions
Xinit(1)=0.5   !p
Xinit(2)=h0  !h 
Xinit(3)=R0  !R
Xinit(4)=0   !rpore

x1=Xinit(1)  !p
x2=Xinit(2)  !h
x3=Xinit(3)  !R
x4=Xinit(4)  !rpore

sigmat=x1*x3/(2*x2)
sigmac=0


t=0

do  
   i=i+1
   t=t+dt


	q1  = dt*dx1(x1,x2,x3,x4)
	w1  = dt*dx2(x1,x2,x3,x4) 
	e1  = dt*dx3(x1,x2,x3,x4) 
        z1  = dt*dx4(x1,x2,x3,x4)
    	
	q2 =  dt*dx1(x1+0.5*q1,x2+0.5*w1,x3+0.5*e1,x4+0.5*z1)
	w2 =  dt*dx2(x1+0.5*q1,x2+0.5*w1,x3+0.5*e1,x4+0.5*z1)
	e2 =  dt*dx3(x1+0.5*q1,x2+0.5*w1,x3+0.5*e1,x4+0.5*z1)
        z2 =  dt*dx4(x1+0.5*q1,x2+0.5*w1,x3+0.5*e1,x4+0.5*z1)
    
	q3 =  dt*dx1(x1+0.5*q2,x2+0.5*w2,x3+0.5*e2,x4+0.5*z2)
	w3 =  dt*dx2(x1+0.5*q2,x2+0.5*w2,x3+0.5*e2,x4+0.5*z2)
	e3 =  dt*dx3(x1+0.5*q2,x2+0.5*w2,x3+0.5*e2,x4+0.5*z2)
        z3 =  dt*dx4(x1+0.5*q2,x2+0.5*w2,x3+0.5*e2,x4+0.5*z2)

	q4 =  dt*dx1(x1+0.5*q3,x2+0.5*w3,x3+0.5*e3,x4+0.5*z3)
	w4 =  dt*dx2(x1+0.5*q3,x2+0.5*w3,x3+0.5*e3,x4+0.5*z3)
	e4 =  dt*dx3(x1+0.5*q3,x2+0.5*w3,x3+0.5*e3,x4+0.5*z3)
	z4 =  dt*dx4(x1+0.5*q3,x2+0.5*w3,x3+0.5*e3,x4+0.5*z3)
        

	x1 = x1+(q1+2*q2+2*q3+q4)/6  !p	
	x2 = x2+(w1+2*w2+2*w3+w4)/6   !h   
        x3 = x3+(e1+2*e2+2*e3+e4)/6  !R
        x4 = x4+(z1+2*z2+2*z3+z4)/6  !rpore
        
!print*,'h2',x5        
        sigmat=sigmat+(EE/x3)*(e1+2*e2+2*e3+e4)/6
        
       

        if ((abs(x4).lt.0.001).and.(sigmat>sigma2)) then   !pore opens

           x4=b0
           dt=1d-13
           rmax=x3
           sigmax=sigmat
           pmax=x1
           hmax=x2
          
           write(400,*)t,t-trupt
           trupt=t
 
           write(100,'(f18.4,8f12.6)')t,x1,x2,x3,x4,sigmat,sigmac
        end if


        if ((abs(x4-b0).lt.0.0001).and.((sigmat<sigma1).or.(abs(x1)<0.01))) then     !pore closes

           x4=0d0
           dt=1d-7
           write(300,'(10f18.5)')t,rmax,rmax-x3,hmax,hmax-x2,sigmax,sigmax-sigmat,pmax,pmax-x1

           write(100,'(f18.4,8f12.6)')t,x1,x2,x3,x4,sigmat,sigmac
  
        end if

    

        if (i==nwrite) then
           write(100,'(f18.4,8f12.6)')t,x1,x2,x3,x4,sigmat,sigmac   !t,p,h,R,rporo,sigmat,sigmac
           nwrite=nwrite+100000
        end if

        if ((t>tmax).or.(x3<=x2)) then
           print*,'t,x2,x3',t,x2,x3
           exit
        end if
end do


contains!---EXTERNAL FUNCTIONS--------------------------------------------------------------


real(8) function dx1(p,h,R,rporo) !dp
use constants
implicit none

real(8)::p,h,R,rporo,ff,Ra
real(8)::A,Q,sigma,htot
 
Q=(rporo**4/mu)*p/h

A=-Kcell*Jv-(Kcell+Kjunct)*p-Q

dx1 = (1/R**2)*(2*A*h*EE-3*p*A*R+ p*Jc/(4*pi*h))


end function

real(8) function dx2(p,h,R,rporo) !dh
use constants
implicit none
real(8)::p,h,R,rporo
real(8)::A,Q

Q=(rporo**4/mu)*p/h

A=-Kcell*Jv-(Kcell+Kjunct)*p-Q


dx2 = Jc/(4*pi*R**2) - 2*A*h/R
end function

real(8) function dx3(p,h,R,rporo) !dr
use constants
implicit none
real(8)::p,h,R,rporo
real(8)::A,Q

Q=(rporo**4/mu)*p/h
A=-Kcell*Jv-(Kcell+Kjunct)*p-Q


dx3 = A
end function


real(8) function dx4(p,h,R,rporo) !drporo
use constants
implicit none
real(8)::p,h,R,rporo
real(8)::A,Q,sigma,a1

dx4=0d0
end function





end program
 
