import numpy 

def idft(inp_arr,x):
	inp_arr=numpy.zero(x)
	basis=numpy.exp(2*numpy.pi*i*inp_arr)
	y=funct(inp_arr,y)

	four_t +=basis*y	
	print 'the invesre fourier transform is' + repr(four_t)

def funct(inp_arr,y):
 	for i in range (1,inp_arr.size):	
		return y[:,i]==y[0::2]
	
if '_name_' == "_main_":
	inp_arr=numpy.arange(-5,5,0.1)
	y=numpy.exp(x**2)
	y1=funct(inp_arr,y)
	y2=fft.fft1(y)	
	print y1
	
	import numpy
from matplotlib import pylab as plt

import numpy
from matplotlib import pyplot as plt
t=numpy.arange(-5,5,0.1)
x_true=t**3-0.5*t**2
x=x_true+10*numpy.random.randn(t.size)

npoly=20  #let's fit 4th order polynomial
ndata=t.size
A=numpy.zeros([ndata,npoly])
A[:,0]=1.0
for i in range(1,npoly):
    A[:,i]=A[:,i-1]*t
#Let's ignore noise for now.  New equations are:
#m=(A^TA)^{-1}*(A^Td)
A=numpy.matrix(A)
d=numpy.matrix(x).transpose()
lhs=A.transpose()*A
rhs=A.transpose()*d
fitp=numpy.linalg.inv(lhs)*rhs
pred=A*fitp
plt.clf();plt.plot(t,x,'*');plt.plot(t,pred,'r');
plt.draw()
plt.savefig('polyfit_example_high.png')

def process_chains(t,sig=0.5,amp=1):
    return (simulate_gaussian(t,sig=0.5,amp=1))
    
def simulate_gaussian(t,sig=0.5,amp=1):
    dat=numpy.exp(-0.5*(t)**2/sig**2)*amp
    dat+=numpy.random.randn(t.size)
    return dat
	print y2
