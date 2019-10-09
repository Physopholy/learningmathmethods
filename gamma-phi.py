import numpy as np

#just comment out inputs/routines not needed
print('Enter guesses for unknown variables.')
t=float(input("Enter T(degC):"))
t=273+t
p=float(input("Enter P(bar):"))
x1=float(input("Enter x1:"))
x2=1-x1
y1=float(input("Enter y1:"))
y2=1-y1
psat1=float(input("Enter Psat1(bar):"))
psat2=float(input("Enter Psat2(bar):"))
b1=float(input("Enter b1 from binaryinteractions.py:"))
b2=float(input("Enter b2 from binaryinteractions.py:"))
b12=float(input("Enter b12 from binaryinteractions.py:"))
del12=2*b12-b1-b2


a=float(input("Enter A_prime for margules eqn.:"))

def margules(x,a_prime): #returns gamma for the given x
	g=2.718281828459045**(a_prime*(1-x)**2)
	return(g)

g1=margules(x1,a)
g2=margules(x2,a)

def gammaphi(y1,y2,p,t,x1,x2,g1,g2,psat1,psat2):

	def phi_solver(b,p,psat,y,del12,t): #just a subroutine for phi
		phi=2.718281828459045**((b*(p-psat)+p*y**2*del12)/(83.1451*t))
		return(phi)

	def obj(y,p,phi,x,gamma,psat): #function to be rooted
		lhs=y*p*phi
		rhs=x*gamma*psat
		f=lhs-rhs
		return(f)

	def deriv(y,dy,p,dp,phi,x,gamma,psat): #centered-difference approximation
		f_0=obj(y-dy,p-dp,phi,x,gamma,psat)
		f_1=obj(y+dy,p+dp,phi,x,gamma,psat)
		dy_dx=(f_1-f_0)/(2*(dy+dp))
		return(dy_dx)

	

	phi1=phi_solver(b1,p,psat1,y2,del12,t) #these lines just initialize error
	phi2=phi_solver(b2,p,psat2,y1,del12,t)
	a=obj(y1,p,phi1,x1,g1,psat1)
	b=obj(y2,p,phi2,x2,g2,psat2)
	Fyp=np.array([[a],[b]])
	error=np.linalg.norm(Fyp)

	while error>0.000001:
		phi1=phi_solver(b1,p,psat1,y2,del12,t)
		phi2=phi_solver(b2,p,psat2,y1,del12,t)

		dobj1_dy1=deriv(y1,.00000001,p,0,phi1,x1,g1,psat1) #check these
		dobj1_dp=deriv(y1,0,p,.00000001,phi1,x1,g1,psat1)
		dobj2_dy1=deriv(y2,.00000001,p,0,phi2,x2,g2,psat2)
		dobj2_dp=deriv(y2,0,p,.00000001,phi2,x2,g2,psat2)

		x=np.array([[y1],[p]])
		a=obj(y1,p,phi1,x1,g1,psat1)
		b=obj(y2,p,phi2,x2,g2,psat2)
		Fyp=np.array([[obj(y1,p,phi1,x1,g1,psat1)],[obj(y2,p,phi2,x2,g2,psat2)]])
		jacob=np.array([[dobj1_dy1,dobj1_dp],[dobj2_dy1,dobj2_dp]])
		jinv=np.linalg.inv(jacob)

		x=x-np.dot(jinv,Fyp)

		#start here, need to figure out how to slice numpy array x

		y1=x[0]
		y2=1-y1
		p=x[1]
		error=np.linalg.norm(Fyp)

	return(p,t,y1,y2,phi1,phi2,x1,x2,g1,g2,psat1,psat2)

answer=gammaphi(y1,y2,p,t,x1,x2,g1,g2,psat1,psat2)
p=answer[0]
t=answer[1]
y1=answer[2]
y2=answer[3]
phi1=answer[4]
phi2=answer[5]
x1=answer[6]
x2=answer[7]
g1=answer[8]
g2=answer[9]
psat1=answer[10]
psat2=answer[11]

print('P(bar)=',p)
print('T(K)=',t)
print('y1=',y1)	
print('y2=',y2)
print('phi1=',phi1)
print('phi2=',phi2)
print('x1=',x1)
print('x2=',x2)
print('gamma_1= ',g1)
print('gamma_2= ',g2)
print('Psat1(bar)= ',psat1)
print('Psat2(bar)= ',psat2)

