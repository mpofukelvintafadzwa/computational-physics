import numpy
from matplotlib import pyplot as plt
class Particles:
    def __init__(self,n=1000,G=1.0,soft=0.1):
        self.x=numpy.random.randn(n)
        self.y=numpy.random.randn(n)
        self.m=numpy.ones(n)*(1.0/n)
        self.vx=numpy.zeros(n)
        self.vy=numpy.zeros(n)
        self.opts={}
        self.opts['n']=n
        self.opts['G']=G
        self.opts['soft']=soft
    def get_forces(self):
        pot=0
        self.fx=numpy.zeros(self.opts['n'])
        self.fy=numpy.zeros(self.opts['n'])
        for i in range(0,self.opts['n']-1):
            dx=self.x[i]-self.x[i+1:]
            dy=self.y[i]-self.y[i+1:]
            rsqr=(dx*dx+dy*dy)
            rsqr[rsqr<self.opts['soft']]=self.opts['soft']
            r=numpy.sqrt(rsqr)
            r3inv=1.0/(r*rsqr)
            self.fx[i]-=numpy.sum(dx*r3inv*self.m[i+1:])
            self.fy[i]-=numpy.sum(dy*r3inv*self.m[i+1:])
            self.fx[i+1:]+=dx*r3inv*self.m[i]
            self.fy[i+1:]+=dy*r3inv*self.m[i]
            pot=pot+numpy.sum(self.opts['G']*self.m[i]*self.m[i+1:]*1.0/r)
        self.fx=self.fx*self.opts['G']
        self.fy=self.fy*self.opts['G']
        return pot
        
    def update_pos(self,dt=0.1):
        self.x+=self.vx*dt
        self.y+=self.vy*dt
        pot=self.get_forces()
        self.vx+=self.fx*dt
        self.vy+=self.fy*dt
        return pot
            
	#from previous assignment		
    def get_potential(self):
        pot=numpy.zeros(self.opts['n'])
        for i in range(0,self.opts['n']):
            dx=self.x[i]-self.x
            dy=self.y[i]-self.y
            r=numpy.sqrt(dx*dx+dy*dy)
            rinverse=1.0/r
            rinverse[i]=0  
            potential[i]=self.m[i]+numpy.sum(self.opts['G']*self.m[i]*self.m*rinverse)
        return potential

if __name__=='__main__':
    part=Particles()
    #pot=part.get_potential()
    pot_org=part.get_forces()
    print "original energy is " + repr(pot_org)
    plt.ion()
    nstep=200
    kk=numpy.zeros(nstep)
    pp=numpy.zeros(nstep)
    for xx in range(0,nstep):
        #pot=part.get_forces()
        pot=part.update_pos(0.05)        
        plt.clf()
        plt.plot(part.x,part.y,'*')
        plt.draw()
        kin=numpy.sum(part.m*(part.vx**2+part.vy**2))
        print 'total energy is ' + repr([pot,kin,kin-2.0*pot])
        pp[xx]=pot
        kk[xx]=kin

#question 2
import numpy
from matplotlib import pyplot as plt

x0=0;y0=0;vx0=0;vy0=0.5
x1=1;y1=0;vx1=0;vy1=-0.5;

dt=0.001
tmax=5
dstep=10

plt.ion()
plt.clf()
step=0
for t in numpy.arange(0,tmax,dt):    
    step+=1
    dx=x0-x1
    dy=y0-y1
    rsquare=dx*dx+dy*dy
    r=numpy.sqrt(rsquare)
    r3=r*r*r
   
    fx0=dx/r3;
    fy0=dy/r3
  
    fx1 =-fx0
    fy1 =-fy0
    
    x0 +=dt*vx0
    y0 +=dt*vy0
    vx0 +=-dt*fx0
    vy0 +=-dt*fy0
    
    x1+= dt*vx1
    y1+=dt*vy1
    vx1 -= dt*fx1
    vy1 -= dt*fy1
	
if step%dstep==0:
        plt.plot(x0,y0,'rx')
        plt.plot(x1,y1,'b*')
        plt.ylim(-1.5,1.5)
        plt.xlim(-1,2)
    
        plt.draw()
        pot=-1.0/r
        kin=0.5*(vx0*vx0+vy0*vy0+vx1*vx1+vy1*vy1)
        print 'kin and pot are '  + repr(kin) + '   ' + repr(pot) + '  ' + repr(pot+kin) 
