import numpy as np

#just comment out inputs/routines not needed
print('Enter guesses for unknown variables.')
#t=float(input("Enter T(degC):"))
t=45
t=273+t
#p=float(input("Enter P(bar):"))
p=1.35
#x1=float(input("Enter x1:"))
x1=.259
x2=1-x1
#y1=float(input("Enter y1:"))
y1=.735
y2=1-y1
#psat1=float(input("Enter Psat1(bar):"))
psat1=2.34218
#psat2=float(input("Enter Psat2(bar):"))
psat2=.44265
#b1=float(input("Enter b1 from binaryinteractions.py:"))
b1=-609.83
#b2=float(input("Enter b2 from binaryinteractions.py:"))
b2=-765.94
#b12=float(input("Enter b12 from binaryinteractions.py:"))
b12=-680.65
del12=2*b12-b1-b2


#a=float(input("Enter A_prime for margules eqn.:"))
a=0.93

def margules(x,a_prime): #returns gamma for the given x
	g=2.718281828459045**(a_prime*(1-x)**2)
	return(g)

g1=margules(x1,a)
g2=margules(x2,a)

def gammaphi(y1,y2,p,t,x1,x2,g1,g2,psat1,psat2):

	def phi_solver(b,p,psat,y,del12,t): #just a subroutine for phi
		phi=2.718281828459045**((b*(p-psat)+p*(1-y)**2*del12)/(83.1451*t))
		return(phi)

	def obj(y,p,b,del12,t,x,gamma,psat): #function to be rooted
		lhs=y*p*2.718281828459045**((b*(p-psat)+p*(1-y)**2*del12)/(83.1451*t))
		rhs=x*gamma*psat
		f=lhs-rhs
		return(f)

	def deriv(y,dy,p,dp,b,del12,t,x,gamma,psat): #centered-difference approximation
		f_0=obj(y-dy,p-dp,b,del12,t,x,gamma,psat)
		#print('f_0=',f_0)
		f_1=obj(y+dy,p+dp,b,del12,t,x,gamma,psat)
		#print('f_1=',f_1)
		dy_dx=(f_1-f_0)/(2*(dy+dp))
		#print('dydx=',dy_dx)
		return(dy_dx)

	

	phi1=phi_solver(b1,p,psat1,y1,del12,t) #these lines just initialize error
	phi2=phi_solver(b2,p,psat2,y2,del12,t)
	Fyp=np.array([[obj(y1,p,b1,del12,t,x1,g1,psat1)],[obj(y2,p,b2,del12,t,x2,g2,psat2)]])
	error=np.linalg.norm(Fyp)
	index=0
	#x=np.array([[y1],[p]])

	while error>0.000001 and index<100:
		index=index+1
		print(index)

		dobj1_dy1=deriv(y1,.00000001,p,0,b1,del12,t,x1,g1,psat1) #check these
		dobj1_dp=deriv(y1,0,p,.00000001,b1,del12,t,x1,g1,psat1)
		dobj2_dy1=deriv(y2,-.00000001,p,0,b2,del12,t,x2,g2,psat2)
		dobj2_dp=deriv(y2,0,p,.00000001,b2,del12,t,x2,g2,psat2)

		x=np.array([[y1],[p]])
		print('x',x)
		Fyp=np.array([[obj(y1,p,b1,del12,t,x1,g1,psat1)],[obj(y2,p,b2,del12,t,x2,g2,psat2)]])
		print('fyp',Fyp)
		jacob=np.array([[dobj1_dy1,dobj1_dp],[dobj2_dy1,dobj2_dp]])
		print('jacob',jacob)
		jinv=np.linalg.inv(jacob)
		print('jinv',jinv)

		x=x-np.dot(jinv,Fyp)
		print('x',x)

		#start here, need to figure out how to slice numpy array x

		y1=float(x[0])
		print('y1',y1)
		y2=1-y1
		print('y2',y2)
		p=float(x[1])
		print('p',p)
		Fyp=np.array([[obj(y1,p,b1,del12,t,x1,g1,psat1)],[obj(y2,p,b2,del12,t,x2,g2,psat2)]])
		error=np.linalg.norm(Fyp)
		print('error',error)

	phi1=phi_solver(b1,p,psat1,y1,del12,t)
	phi2=phi_solver(b2,p,psat2,y2,del12,t)

	return(p,t,y1,y2,phi1,phi2,x1,x2,g1,g2,psat1,psat2,index)

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
print('Steps taken:',answer[12])

