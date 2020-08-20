# Dynamic Decouping of a Single Qubit

## DDC_Pulse tagging

```python
import math

def DD_Pulse(n,i,):
    
    A_i=(math.sin(math.pi*i/(2*(n+1))))**2
    
    return A_i

def DD_data(n,fidelity):
    
    DD_Pulse_List=[]

    i=1
    while i<n+1:
        Pulse=DD_Pulse(n,i)
        DD_Pulse_List.append(Pulse)
        
        i+=1
    return DD_Pulse_List
        
    
    
A=DD_data(4,0.9)
B=DD_data(5,0.5)
print(A)
print(B)

ci=2
n=0
ANS=0
List=[]
while n<len(A):
    List.append(B[n])
    List.append(A[n])
    n+=1
List.append(B[n])
print(List)

n=0
ANS+=(List[n])**ci
while n<len(List)-1:
    ANS+=(((-1)**(n+1))*((List[n+1])**ci-(List[n])**ci))
    print(ANS)
    n+=1

ANS+=(-1)*(1-(List[n]**ci))
print(ANS)
```

## Spin-Echo efficiency simulation

```python
from qutip import *
import math  
import numpy as np
import random
import matplotlib.pyplot as plt


class DD_Coupling:
    def __init__(self, view=[-60,30]):
        
        self.Phi=math.pi/2
        self.Phidefault=math.pi/2

        self.wX=-0.01
        self.wY=0.01
        self.wZ=0.01
        
        self.PhiList=[]
        self.ListX=[]
        self.ListY=[]
    
        
        
    def B_Random(self,Anv,Dev):
        ans=random.gauss(Anv,Dev)
        return ans
    
    # 产生随机磁场
    def RandomB(self,num):
        B=np.zeros(num)
        Tlist=np.arange(num)

        #选择各个频段
        for e_1 in Tlist:

            #选择每一项
            for e_2 in Tlist:

                if e_2%e_1==0:
                    B_A=random.gauss(0.000,0.002)

                #print(B_A)
                B[e_2]=B[e_2]+B_A

        return B     

        
# 厄米算符操作
    def SigX(self):
        self.Z=-self.Z
        self.Y=-self.Y
    
    def SigY(self):
        self.Z=-self.Z
        self.X=-self.X
        
    def SigZ(self):
        self.X=-self.X
        self.Y=-self.Y
        
        
        
        
# Dephasing自然演化------------------------------------------------------------------------------
    
    def U(self,times_all):
        #time_mark=50
        
        self.PhiDefault=self.Phi
        
        self.FAll_1=0
        self.FAll_2=0
        
        
        for e in range(100):
            
            self.Phi=math.pi/2
            
            B_list=self.RandomB(times_all)
            
        
            time=0

            
            
            while time<times_all:
                
                
                self.Phi+=B_list[time]
                

                time+=1
                
            self.PhiList.append(self.Phi)
                
            self.ListX.append(math.cos(self.Phi))
            self.ListY.append(math.sin(self.Phi))

            
            self.F2=math.sin(self.Phi)
            
            self.FAll_2+=self.F2
            
        self.FAll_2=self.FAll_2/100
        
    def USP(self,times_all,mark):
        time_mark=mark
        
        self.PhiDefault=self.Phi
        
        self.FAll_1=0
        self.FAll_2=0
        
        for e in range(100):
            
            self.Phi=math.pi/2
            
            B_list=self.RandomB(times_all)
            
            time=0
            
            
            while time<times_all:

                if time%time_mark==0:

                    self.Phi=-self.Phi
                    self.Phi+=B_list[time]

                else:
                
                    self.Phi+=B_list[time]
                

                time+=1
                
            self.PhiList.append(self.Phi)
                
            self.ListX.append(math.cos(self.Phi))
            self.ListY.append(math.sin(self.Phi))

            #print(self.AngZ)
            #self.F1=(1/4)*((1+math.cos(self.Phi-math.pi/2))**2+math.sin(self.Phi-math.pi/2)**2)
            #print(self.F)
            self.F2=math.sin(self.Phi)
            
            if self.F2<0:
                self.F2=-self.F2
                #self.F1=-self.F1

            #self.FAll_1+=self.F1
            self.FAll_2+=self.F2
            
        #self.FAll_1=self.FAll_1/1000
        self.FAll_2=self.FAll_2/100
        
    def F_TU(self):
        
        Tlist=np.arange(100)
        Flist_1=[]
        Flist_2=[]
        for e in Tlist:
            
            self.U(e)
            Flist_1.append(self.FAll_1)
            Flist_2.append(self.FAll_2)
        
        #plt.plot(Tlist,Flist_1)
        #plt.show()
        
        plt.plot(Tlist,Flist_2)
        #plt.show()
        
        #plt.scatter(self.ListX,self.ListY)
        #plt.show()

    def F_TS(self,mark):
        
        Tlist=np.arange(100)
        Flist_1=[]
        Flist_2=[]
        for e in Tlist:
            
            self.USP(e,mark)
            Flist_1.append(self.FAll_1)
            Flist_2.append(self.FAll_2)
            
        
        #plt.plot(Tlist,Flist_1)
        #plt.show()
        
        plt.plot(Tlist,Flist_2)
        #plt.show()
        
        #plt.scatter(self.ListX,self.ListY)
        #plt.show()

        

AA=DD_Coupling([-60,30])
AA.F_TU()
#AA.F_TS(10)
AA.F_TS(25)
AA.F_TS(50)
plt.show()
```
