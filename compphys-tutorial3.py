from matplotlib import pyplot as plt
from numpy.fft import fft,ifft
import numpy

#Question 1

def makeshift(x,n=0):
	zerovec=0*x 
	zerovec[n]=1
	vecft=fft(zerovec)
	xft=fft(x)
	return numpy.real(ifft(xft*vecft))
if __name__=='__main__':
	x=numpy.arange(-30,30,0.1)
	sigma=2
	y=numpy.exp(-0.5*x**2/sigma**2)
	yshift=makeshift(y,y.size/2)
	plt.ion()
	plt.plot(x,y)
	plt.plot(x,yshift)

#Question 2
from matplotlib import pyplot as plt
from numpy.fft import fft,ifft
import numpy
def takecorr(x,y):
	assert(x.size==y.size)
	xft=fft(x)
	yft=fft(y)
	yftconj=numpy.conj(yft)
	return numpy.real(ifft(xft*yftconj))
if __name__=='__main__':
	x=numpy.arange(-30,30,0.1)
	sigma=2
	y=numpy.exp(-0.5*x**2/sigma**2)
	ycorrela=takecorrela(y,y)
	plt.plot(x,ycorrela)
	plt.show()

#Question 3
from numpy.fft import fft,ifft
import numpy
from matplotlib import pyplot as plt
def myshift(x,n=0):
	zerovec=0*x 
	zerovec[n]=1
	vecft=fft(zerovec)
	xft=fft(x)
	return numpy.real(ifft(xft*vecft))
def takecorr(x,y):
	assert(x.size==y.size) 
	xft=fft(x)
	yft=fft(y)
	yftconj=numpy.conj(yft)
	return numpy.real(ifft(xft*yftconj))
if __name__=='__main__':
	x=numpy.arange(-20,20,0.1)
	sigma=2
	y=numpy.exp(-0.5*x**2/sigma**2)
	ycorr=takecorr(y,y)
	yshift=myshift(y,y.size/4)
	yshiftcorr=takecorr(yshift,yshift)
	meanerr=numpy.mean(numpy.abs(ycorr-yshiftcorr))
	print 'mean difference between the two correlation functions is ' + repr(meanerr)
	plt.plot(x,ycorr)
	plt.plot(x,yshiftcorr)
	plt.show()

#Question 4
from numpy.fft import fft,ifft
import numpy
from matplotlib import pyplot as plt
def conv_nowrap(x,y):
	assert(x.size==y.size) 
	xx=numpy.zeros(2*x.size)
	xx[0:x.size]=x
	yy=numpy.zeros(2*y.size)
	yy[0:y.size]=y
	xxft=fft(xx)
	yyft=fft(yy)
	vec=numpy.real(ifft(xxft*yyft))
	return vec[0:x.size]
if __name__=='__main__':
	x=numpy.arange(-30,30,0.1)
	sigma=2
	y=numpy.exp(-0.5*x**2/sigma**2)
	y=y/y.sum()
	yconv=conv_nowrap(y,y)
	plt.plot(x,y)
	plt.plot(x,yconv)
	plt.show()

#Bonus
import numpy
from numpy import concatenate,exp,pi,arange,complex
def newfft(vec):
	n=vec.size

	if n==1:
		return vec

	myeven=vec[0::2]
	myodd=vec[1::2]
	nn=n/2;
	j=complex(0,1)

	twid=exp(-2*pi*j*arange(0,nn)/n)

	eft=newfft(myeven)
	oft=newfft(myodd)

	result=concatenate((eft+twid*oft,eft-twid*oft))
	return result
def fft3(vec):
	n=vec.size
	if n==1:
		return vec
	a=vec[0::3]
	b=vec[1::3]
	c=vec[2::3]
	j=complex(0,1)
	nn=n/3
	twid1=exp(-2*pi*j*arange(0,nn)/n)
	twid2=exp(-4*pi*j*arange(0,nn)/n)
	f1=exp(-2*pi*j/3) 
	f2=exp(-4*pi*j/3)
	f1b=f2; 
	f2b=f1; 
	aft=fft3(a)
	bft=fft3(b)*twid1
	cft=fft3(c)*twid2
	ft1=aft+bft+cft
	ft2=aft+bft*f1+cft*f2
	ft3=aft+bft*f1b+cft*f2b
	ft=concatenate((ft1,ft2,ft3))
	return ft
