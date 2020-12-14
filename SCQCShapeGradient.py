from qutip import *
import numpy as np
import matplotlib.pyplot as plt

def HX_coeff(t, args):
    
    global AxT
    global TmindxT
    global SigGaussxT
    return AxT * np.exp(-((t-TMindxT) / SigGaussxT) ** 2)

def HY_coeff(t, args):
    
    global AyT
    global TmindyT
    global SigGaussyT
    return AyT * np.exp(-((t-TMindyT) / SigGaussyT) ** 2)

def Run(XInit_1,XInit_2,XInit_3,YInit_1,YInit_2,YInit_3):
    
    global AxT
    AxT=XInit_1
    
    global TmindxT
    TMindxT=XInit_2
    
    global SigGaussxT
    SigGaussxT=XInit_3
    
    global AyT
    AyT=YInit_1
    
    global TmindyT
    TMindyT=YInit_2
    
    global SigGaussyT
    SigGaussyT=YInit_3

    UGoal = Qobj([[0,1,0],
              [1,0,0],
              [0,0,1]])

    H0 = Qobj([[0,0,0],
               [0,0,0],
               [0,0,1]])

    HX = Qobj([[0,1,0],
               [1,0,1.41421356],
               [0,1.41421356,0]])

    HY = Qobj([[0,0-1j,0],
               [0+1j,0,-1.41421356j],
               [0,1.41421356j,0]])
    
    H = [0.3*H0,[HX,HX_coeff],[HY,HY_coeff]]

    # initial state
    QInitList=[Qobj([[1],[0],[0]]),
               Qobj([[0],[1],[0]]),
               Qobj([[0],[0],[1]])]

    # list of times for which the solver should store the state vector
    tlist = np.linspace(0, 10, 1000)

    # Time Evolution Opreator
    UtList=[[],[],[]]

    E=[0,1,2]

    for QInit in QInitList:
        result = mesolve(H, QInit, tlist,  [])
        FinalState=result.states[-1]
        for e in E:

            UtList[e].append(FinalState[e][0][0])
        #print(result.states[-1])

    Ut=Qobj(UtList)
    
    global UtFinal
    UtFinal=Ut

    G=abs((UGoal.dag()*Ut).tr())
    return G


Ax=2
TMindx=2
SigGaussx=5
Ay=3
TMindy=8
SigGaussy=5

AxT=Ax
TMindxT=TMindx
SigGaussxT=SigGaussx
AyT=Ay
TMindyT=TMindy
SigGaussyT=SigGaussy

UtFinal=0

Delta=0.01

FuncPre=0
FuncLast=1
while FuncPre<FuncLast:
    
    Kx1=(Run(Ax+Delta,TMindx,SigGaussx,Ay,TMindy,SigGaussy)-Run(Ax,TMindx,SigGaussx,Ay,TMindy,SigGaussy))/Delta
    Kx2=(Run(Ax,TMindx+Delta,SigGaussx,Ay,TMindy,SigGaussy)-Run(Ax,TMindx,SigGaussx,Ay,TMindy,SigGaussy))/Delta
    Kx3=(Run(Ax,TMindx,SigGaussx+Delta,Ay,TMindy,SigGaussy)-Run(Ax,TMindx,SigGaussx,Ay,TMindy,SigGaussy))/Delta

    Ky1=(Run(Ax,TMindx,SigGaussx,Ay+Delta,TMindy,SigGaussy)-Run(Ax,TMindx,SigGaussx,Ay,TMindy,SigGaussy))/Delta
    Ky2=(Run(Ax,TMindx,SigGaussx,Ay,TMindy+Delta,SigGaussy)-Run(Ax,TMindx,SigGaussx,Ay,TMindy,SigGaussy))/Delta
    Ky3=(Run(Ax,TMindx,SigGaussx,Ay,TMindy,SigGaussy+Delta)-Run(Ax,TMindx,SigGaussx,Ay,TMindy,SigGaussy))/Delta
    
    print(Kx1,Kx2,Kx3,Ky1,Ky2,Ky3)

    Kr=(Kx1**2+Kx2**2+Kx3**2+Ky1**2+Ky2**2+Ky3**2)**0.5

    Run(Ax,TMindx,SigGaussx,Ay,TMindy,SigGaussy)

    FuncPre=Run(Ax,TMindx,SigGaussx,Ay,TMindy,SigGaussy)
    
    
    Ax+=(Kx1/Kr)*Delta
    TMindx+=(Kx2/Kr)*Delta
    SigGaussx+=(Kx3/Kr)*Delta

    Ay+=(Ky1/Kr)*Delta
    TMindy+=(Ky2/Kr)*Delta
    SigGaussy+=(Ky3/Kr)*Delta

    FuncLast=Run(Ax,TMindx,SigGaussx,Ay,TMindy,SigGaussy)

    print(Run(Ax,TMindx,SigGaussx,Ay,TMindy,SigGaussy))
    print(Ax,TMindx,SigGaussx,Ay,TMindy,SigGaussy)
    print('ok')

print('Done')
print(Run(Ax,TMindx,SigGaussx,Ay,TMindy,SigGaussy))
print(Ax,TMindx,SigGaussx,Ay,TMindy,SigGaussy)
print(UtFinal)

