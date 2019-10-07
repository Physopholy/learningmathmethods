#x1=x0-f(x0)/f'(x0), approximate a root

#initialize parameters
x=float(input("Enter initial guess:"))
dx=float(input("Enter step-size:"))
conv=float(input("Enter convergence tolerance:"))

def newton(x,dx,conv):
	def obj(x): #function to be rooted (hah)
		f=2-x**2
		return(f)
	def deriv(x,dx): #centered-difference approximation
		f_0=obj(x-dx)
		f_1=obj(x+dx)
		dy_dx=(f_1-f_0)/(2*dx)
		return(dy_dx)
	def e(x): #calculate error
		e=obj(x)
		return(e)

	e=e(x) #initialize error
	abse=abs(e) 

	t=0 #just a step counter


	while abse>conv:
		x=x-obj(x)/deriv(x,dx)
		t=t+1
		e=e(x)
		abse=abs(e)

	print('x=',x)
	print('steps taken:',t)
	print('final error:',e)

	return(x)

newton(x,dx,conv)


