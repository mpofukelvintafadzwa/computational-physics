#question 1c: error does go down with the square of the number of points  for the Gaussian initial conditions

#Question 1a

import numpy
from matplotlib import pyplot as plt

class advect:
    def __init__(self,n=3000,u=1.0,xmax=1.0,C=0.1):
        self.rho=numpy.zeros(n)
        self.rho[n/3:2*n/3]=1.0
        self.u=u
        self.n=n
        self.dx=(0.0+xmax)/(0.0+n)
        self.C=C
    def get_ts(self):
        return numpy.abs(self.dx/self.u)*self.C

    def get_bc(self):
        self.rho[0:2]=self.rho[-4:-2]
        self.rho[-2:]=self.rho[2:4]
    
    def ic_gaussian(self,sig=0.1):
        xvec=numpy.arange(0,self.rho.size)*self.dx
        xvec=xvec-xvec.mean()
        self.rho=numpy.exp(-0.5*xvec**2/sig**2)
 
#Question 1b 
class advect1(advect):
    def update(self):
        self.get_bc()
        df=self.rho[0:-1]-self.rho[1:]
        self.rho[1:]+=df*self.C

    def test_update(self):
        rho_org=self.rho.copy()
        rhoft=numpy.fft.fft(rho_org)

        kvec=numpy.arange(0,self.rho.size)
        kvec=kvec/(self.rho.size+0.0)*2*numpy.pi
        i=numpy.complex(0,1.0)
        fac=1+self.C*numpy.exp(-i*kvec)-self.C
        newft=rhoft*fac
        dat_check=numpy.real(numpy.fft.ifft(newft)) 
        self.update()
        mean_err=numpy.mean(numpy.abs(dat_check-self.rho))
        if mean_err/numpy.mean(numpy.abs(self.rho))<1e-13:
            print 'test successful with mean error ' + repr(mean_err)
        else:
            print 'test failed with mean error ' + repr(mean_err)

        self.rho=rho_org

class advect2(advect):
    def get_interfaces(self):        
        self.myderiv=0.5*(self.rho[2:]-self.rho[0:-2])/self.dx
        self.right=self.rho[1:-1]+0.5*self.dx*(1.0-self.C)*self.myderiv
        #self.left=self.rho[1:-1]-0.5*self.dx*(1.0+self.C)*self.myderiv       
    def update(self):
        #dt=self.get_ts()
        self.get_bc()
        self.get_interfaces()
        if self.u>0:
            self.rho[2:-2]=self.rho[2:-2]+(self.right[0:-2]-self.right[1:-1])*self.C
        else:
            self.rho[2:-2]=self.rho[2:-2]+(self.left[1:]-self.left[0:-2])*self.C
        
if __name__=="__main__":
    fac=10;
    npt=200
    C=1.0/fac
    fwee1=advect1(n=npt,C=C)
    fwee2=advect2(n=npt,C=C)
    fwee_org=advect2(n=npt,C=C)
    if False:  
        sig=0.1
        fwee1.ic_gaussian(sig=sig)
        fwee2.ic_gaussian(sig=sig)
        fwee_org.ic_gaussian(sig=sig)
        mytag='gaussian'
    else:
        mytag='rectangle'
    plt.ion()
    fwee1.test_update()
    for i in range(0,npt-4):
        for ii in range(0,fac):
            fwee1.update()
            fwee2.update()
        if i%10==-1:  
            plt.clf()
            plt.plot(fwee1.rho)
            plt.draw()
    print numpy.mean(numpy.abs(fwee1.rho-fwee_org.rho))
    print numpy.mean(numpy.abs(fwee2.rho-fwee_org.rho))
plt.clf()
plt.plot(fwee1.rho)
plt.plot(fwee2.rho)
plt.plot(fwee_org.rho)
plt.draw()
plt.savefig('tut1_' + mytag + '.png')


import numpy
from matplotlib import pyplot as plt
class advect:
    def __init__(self,n=3000,u=1.0,xmax=1.0,C=0.1):
        self.rho=numpy.zeros(n)
        self.rho[n/3:2*n/3]=1.0
        self.u=u
        self.n=n
        self.dx=(0.0+xmax)/(0.0+n)
        self.C=C
    def get_ts(self):
        return numpy.abs(self.dx/self.u)*self.C
    def get_bc(self):
        self.rho[0:2]=self.rho[-4:-2]
        self.rho[-2:]=self.rho[2:4]
    def ic_gaussian(self,sig=0.1):
        xvec=numpy.arange(0,self.rho.size)*self.dx
        xvec=xvec-xvec.mean()
        self.rho=numpy.exp(-0.5*xvec**2/sig**2)
        
class advect1(advect):
    def update(self):
        self.get_bc()
        df=self.rho[0:-1]-self.rho[1:]
        self.rho[1:]+=df*self.C

class advect2(advect):
    def get_interfaces(self):              
        self.myderiv=0.5*(self.rho[2:]-self.rho[0:-2])/self.dx       
        self.right=self.rho[1:-1]+0.5*self.dx*(1.0-self.C)*self.myderiv
        #self.left=self.rho[1:-1]-0.5*self.dx*(1.0+self.C)*self.myderiv
        
    def update(self):
        #dt=self.get_ts()
        self.get_bc()
        self.get_interfaces()
        if self.u>0:
            self.rho[2:-2]=self.rho[2:-2]+(self.right[0:-2]-self.right[1:-1])*self.C
        else:
            self.rho[2:-2]=self.rho[2:-2]+(self.left[1:]-self.left[0:-2])*self.C

#question 2b: Yes the minimod limiter  gets rid of the overshooting. It is not very accurate, halfway between 1st and 2nd order.

#Question 2a 
class advect_minimod(advect2):
    def get_interfaces(self):         
        #self.myderiv=0.5*(self.rho[2:]-self.rho[0:-2])/self.dx
        deriv1=(self.rho[2:]-self.rho[1:-1])/self.dx
        deriv2=(self.rho[1:-1]-self.rho[0:-2])/self.dx
        mysign=deriv1>0
        deriv1=numpy.abs(deriv1)
        deriv2=numpy.abs(deriv2)
        ii=(deriv1>deriv2)+0.0
        deriv=deriv1.copy()
        deriv[deriv1>deriv2]=deriv2[deriv1>deriv2]
        deriv=deriv*(2.0*mysign-1.0)       
        deriv[deriv1*deriv2<0]=0.0
        self.myderiv=deriv
        self.right=self.rho[1:-1]+0.5*self.dx*(1.0-self.C)*self.myderiv
        #self.left=self.rho[1:-1]-0.5*self.dx*(1.0+self.C)*self.myderiv
       
if __name__=="__main__":

    fac=10;
    npt=500
    C=1.0/fac

    fwee1=advect1(n=npt,C=C)
    fwee2=advect2(n=npt,C=C)
    fwee_minimod=advect_minimod(n=npt,C=C)

    fwee_org=advect2(n=npt,C=C)
    if True:  
        sig=0.1
        fwee1.ic_gaussian(sig=sig)
        fwee2.ic_gaussian(sig=sig)
        fwee_minimod.ic_gaussian(sig=sig)
        fwee_org.ic_gaussian(sig=sig)
        mytag='gaussian'
    else:
        mytag='rectangle'
    plt.ion()
    for i in range(0,npt-4):
        for ii in range(0,fac):
            fwee1.update()
            fwee2.update()
            fwee_minimod.update()
        if i%10==-1:  
            plt.clf()
            plt.plot(fwee1.rho)
            plt.draw()
    print numpy.mean(numpy.abs(fwee1.rho-fwee_org.rho))
    print numpy.mean(numpy.abs(fwee2.rho-fwee_org.rho))
    print numpy.mean(numpy.abs(fwee_minimod.rho-fwee_org.rho))
plt.clf()
plt.plot(fwee1.rho)
plt.plot(fwee2.rho)
plt.plot(fwee_minimod.rho)
plt.plot(fwee_org.rho)
plt.draw()
plt.savefig('tut2_' + mytag + '.png')
