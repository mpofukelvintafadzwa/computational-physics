#question 1
class complexclass:
    def __init__(self,r=0,imag=0):
        self.r=r
        self.i=imag
    def copy(self):
        return complexclass(self.r,self.i)

    def __add__(self,val):
        ans=self.copy()
        if isinstance(val,complexclass):
            ans.r=ans.r+val.r
            ans.imag=ans.i+val.i
        else:
            ans.r=ans.r+val
        return ans
    def __mul__(self,val):
        ans=self.copy()
        if isinstance(val,complexclass):
            ans.r=self.r*val.r-self.i*val.i
            ans.i=self.r*val.i+self.i*val.r
        else:
            ans.r=ans.r*val
            ans.i=ans.i*val
        return ans
    def __sub__(self,val):
        ans=self.copy()
        if isinstance(val,complexclass):
            ans.r=ans.r+val.r
            ans.imag=ans.i-val.i
        else:
            ans.r=ans.r-val
        return ans

    def __div__(self,val):

        if isinstance(val,complexclass):
            val=val.copy()
            val.i=-1*val.i
            ans=self*val
            myabs=val.r**2+val.i**2
            ans=ans*(1.0/myabs)
        else:
            ans=self*(1.0/val)
        return ans
        
    def __repr__(self):
        if (self.i<0):
            return repr(self.r)+' - '+repr(-1*self.i) +'i'
        else:
            return repr(self.r)+' + '+repr(self.i) +'i'
    def __lshift__(self,crud):
        self.i=-1*self.i

#question2
import numpy
class Particles:
    def __init__(self,n=1000,G=1.0):
        self.x=numpy.random.randn(n)
        self.y=numpy.random.randn(n)
        self.m=numpy.ones(n)
        self.vx=numpy.zeros(n)
        self.vy=numpy.zeros(n)
        self.opts={}
        self.opts['n']=n
        self.opts['G']=G
    def calc_potential(self):
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
    pot=part.calc_potential()

#bonus
import math
class complexclass:
    def __init__(self,r=0,imag=0):
        self.r=r
        self.i=imag
    def copy(self):
        return complexclass(self.r,self.i)

    def __add__(self,val):
        ans=self.copy()
        if isinstance(val,complexclass):
            ans.r=ans.r+val.r
            ans.imag=ans.i+val.i
        else:
            ans.r=ans.r+val
        return ans
    def __mul__(self,val):
        ans=self.copy()
        if isinstance(val,complexclass):
            ans.r=self.r*val.r-self.i*val.i
            ans.i=self.r*val.i+self.i*val.r
        else:
            ans.r=ans.r*val
            ans.i=ans.i*val
        return ans
    def __sub__(self,val):
        ans=self.copy()
        if isinstance(val,complexclass):
            ans.r=ans.r+val.r
            ans.imag=ans.i-val.i
        else:
            ans.r=ans.r-val
        return ans
    def __div__(self,val):
        if isinstance(val,complexclass):
            val=val.copy()
            val.i=-1*val.i
            ans=self*val
            myabs=val.r**2+val.i**2
            ans=ans*(1.0/myabs)
        else:
            ans=self*(1.0/val)
        return ans
    def __pow_simple__(self,val):
        ang=math.atan2(self.i,self.r)
        abs=math.sqrt(self.i*self.i+self.r*self.r)
        ans=self.copy()
        newabs=abs**val
        newang=ang*val
        ans.r=newabs*math.cos(newang)
        ans.i=newabs*math.sin(newang)
        return ans
    def __pow__(self,val):
        if isinstance(val,complexclass):
            ang=math.atan2(self.i,self.r)
            myabs=math.sqrt(self.i*self.i+self.r*self.r)
            myexp=complexclass(math.log(myabs),ang)
            totexp=myexp*val
            newabs=math.exp(totexp.r)
            newang=complexclass(math.cos(totexp.i),math.sin(totexp.i))
            return newang*newabs
        else:
            return self.__pow_simples__(val)

    def __repr__(self):
        if (self.i<0):
            return repr(self.r)+' - '+repr(-1*self.i) +'i'
        else:
            return repr(self.r)+' + '+repr(self.i) +'i'
    def __lshift__(self,crud):
        self.i=-1*self.i
