import numpy
from matplotlib import pyplot as plt
class advection:
    def __init__(self,npart=300,u=-1.0,dx=1.0):
        x=numpy.zeros(npart)
        x[npart/3:2*npart/3]=1.0;
        self.x=x
        self.u=u
        self.dx=dx
    def get_bc_periodic(self):
        self.x[0]=self.x[-2]
        self.x[-1]=self.x[1]
    def update(self,dt=1.0):
        self.get_bc_periodic()
        delt=self.x[1:]-self.x[0:-1]
        self.x[1:-1]+=self.u*dt/self.dx*delt[1:]

if __name__=='__main__':
    material=advection()
    plt.ion()
    plt.plot(material.x)
    plt.show()
    for i in range(0,300):
        material.update()
        plt.clf()
        plt.plot(material.x)
        plt.draw()


#question 2
#increasing the grid resolution by a factor of 10 decreases the time step by a factor of 10.According to the CFL
#conditions


#question 3
import numpy
from matplotlib import pyplot as plt

class Particles:
    def __init__(self,npart=300,xmax=1.0,u=1.0):
        self.x=numpy.arange(npart)/(0.0+npart)*xmax
        self.u=u*(xmax-self.x)/xmax
    def update(self,dt=0.01):
        self.x+=self.u*dt
    def get_density(self,dx=0.01):
        xmin=numpy.min(self.x)
        xmax=numpy.max(self.x)
        nbin=numpy.round(1+(xmax-xmin)/dx)
        myind=numpy.round( (self.x-xmin)/dx)
        rho=numpy.zeros(nbin)

        assert(myind.max()<nbin) 
        for i in numpy.arange(0,myind.size):
            rho[myind[i]]+=1.0
        xvec=numpy.arange(0,nbin)*dx+xmin
        return rho,xvec
if __name__=='__main__':
    part=Particles(npart=30000)
    plt.ion()
    plt.plot(part.x)
    plt.show()

    plt.clf()
    for ii in range(0,200):
        part.update(dt=0.01)
        rho,x=part.get_density()
        plt.plot(x,rho)
        plt.draw()


#question bonus
import numpy
from matplotlib import pyplot as plt
class advect:
    def __init__(self,npart=300,u=1.0,dx=1.0):
        x=numpy.zeros(npart)
        x[npart/3:2*npart/3]=1.0;
        self.x=x
        self.u=u
        self.dx=dx
    def get_bc_periodic(self):
        self.x[0]=self.x[-2]
        self.x[-1]=self.x[1]
    def update(self,dt=1.0):
        self.get_bc_periodic()
        delt=self.x[1:]-self.x[0:-1]
        self.x[1:-1]+=self.u*dt/self.dx*delt[1:]

    def get_state_fft(self,t,dt=1.0):
        nstep=t/dt
        xft=numpy.fft.fft(self.x)
        kvec=numpy.arange(0,self.x.size)
        kvec=kvec/(self.x.size+0.0)*2*numpy.pi
        C=dt*self.u/self.dx
        i=numpy.complex(0,1.0)
        fac=1+C*numpy.exp(-i*kvec)-C
        newft=xft*(fac**nstep)
        return numpy.real(numpy.fft.ifft(newft))
    def get_half_k(self,t,dt=1.0):
        nstep=t/dt
        C=dt*self.u/self.dx
        alpha=2.0**(-0.5/nstep)
        cosk=(alpha-1+2*C-2*C*C)/(2*C*(1-C))
        kcrit=numpy.arccos(cosk)
        return kcrit

        

if __name__=='__main__':
    stuff=advect()
    tmax=stuff.x.size
    newx=stuff.get_state_fft(tmax)
    plt.clf()
    plt.plot(stuff.x)
    plt.plot(newx)
    plt.draw()
    

  
    tvec=numpy.array([10,20,50,100,200,300,600,900,1200])
    dt=0.5
    kcrit=stuff.get_half_k(tvec,dt)
    

    plt.clf()
    plt.plot(tvec,kcrit)    
    h=plt.gca()
    h.loglog()
    h.set_xlabel('Elapsed time')
    h.set_ylabel('Half-power K')
    h.set_title('Half power wavenumber as a function of time for dt=' + repr(dt))
    plt.draw()
    plt.savefig('half_power_k.png')
    
    tt=55
    kk=stuff.get_half_k(tt,dt)
    stuff.x[:]=0
    stuff.x[150]=1
    newx=stuff.get_state_fft(tt,dt)
    newft=numpy.abs(numpy.fft.fft(newx))
    oldft=numpy.abs(numpy.fft.fft(stuff.x))
    rat=newft/oldft
    i=0
    while rat[i]>0.5:
        i=i+1
    
    print 'k observed is ' + repr(i*numpy.pi/stuff.x.size) + ' and expected is ' + repr(kk)
