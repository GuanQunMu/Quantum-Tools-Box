from qutip import *
import math  
import numpy as np
import random
import matplotlib.pyplot as plt


class BlockClass:
    def __init__(self, view=[-60,30], X=0,Y=0,Z=-1,thi=math.pi, phi=0):
        
        self.thi=thi
        self.phi=phi
        self.Block = Bloch()
        self.Block.view=view
        
        self.X = X
        self.Y = Y
        self.Z = Z

        self.XList=[]
        self.YList=[]
        self.ZList=[]
        
        self.wX=-0.01
        self.wY=0.01
        self.wZ=0.01
        
        
# |State> = cos(thi/2) |0> + exp(i * phi)sin(thi/2) |1>
# [x,y,z] = [sin(thi)cos(phi) , sin(thi)sin(phi) , cos(thi)]
    
# 基础工具箱----------------------------------------------------------------------
# 展示Block Sphere
    def Show(self):
        self.Block.show()
        
        
    def B_Random(self,Anv,Dev):
        ans=random.gauss(Anv,Dev)
        return ans
    
# 沿着Z为轴,右手系旋转
    def RotateZ(self,Ang):
        if self.X>0:
            self.AngZ=np.arctan(self.Y/self.X)

        elif self.X<0:
            self.AngZ=np.arctan(self.Y/self.X)+math.pi

        elif self.X==0:
            if self.Y>0:
                self.AngZ=math.pi/2

            elif self.Y<0:
                self.AngZ=-math.pi/2
                
            elif self.Y==0:
                self.AngZ=0
                
        self.R=(self.X**2+self.Y**2)**0.5
        self.X=math.cos(self.AngZ+Ang)*self.R
        self.Y=math.sin(self.AngZ+Ang)*self.R

# 沿着X为轴,右手系旋转
    def RotateX(self,Ang):
        if self.Y>0:
            self.AngX=np.arctan(self.Z/self.Y)

        elif self.Y<0:
            self.AngX=np.arctan(self.Z/self.Y)+math.pi

        elif self.Y==0:
            if self.Z>0:
                self.AngX=math.pi/2

            elif self.Z<0:
                self.AngX=-math.pi/2
            
            elif self.Z==0:
                self.AngX=0
                
        self.R=(self.Z**2+self.Y**2)**0.5
        self.Y=math.cos(self.AngX+Ang)*self.R
        self.Z=math.sin(self.AngX+Ang)*self.R
        
# 沿着Y为轴,右手系旋转
    def RotateY(self,Ang):
        if self.Z>0:
            self.AngY=np.arctan(self.X/self.Z)

        elif self.Z<0:
            self.AngY=np.arctan(self.X/self.Z)+math.pi

        elif self.Z==0:
            if self.X>0:
                self.AngY=math.pi/2

            elif self.X<0:
                self.AngY=-math.pi/2
            
            elif self.X==0:
                self.AngY=0
                
        self.R=(self.X**2+self.Z**2)**0.5
        self.Z=math.cos(self.AngY+Ang)*self.R
        self.X=math.sin(self.AngY+Ang)*self.R
        
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
        
        
# 工具箱结束-------------------------------------------------------------------------------------------
# 态的演化---------------------------------------------------------------------------------------------
    def RunPhaseLock(self,scale_1,scale_2):
        
        self.beta(500,scale_1,scale_2)
        
        #print(self.BetaList)
        #print(self.MtList)
        
        
        self.X=0
        self.Y=1
        self.Z=0
        
        for time in self.TimeList:
            
            # print([self.X,self.Y,self.Z])
            
            self.XList.append(self.X)
            self.YList.append(self.Y)
            self.ZList.append(self.Z)
            
            self.wZ = self.MtList[time]*math.cos(self.BetaList[time])
            self.wY = self.MtList[time]*math.sin(self.BetaList[time])
            
            # print(self.wZ)

            
            self.RotateY(self.wY)
            self.RotateZ(self.wZ)
            time+=1
        # print('X',self.X)
        # print('Y',self.Y)
            
        self.Block.add_points([self.XList,self.YList,self.ZList])
        
        if self.X>0:
            self.AngZ=np.arctan(self.Y/self.X)

        elif self.X<0:
            self.AngZ=np.arctan(self.Y/self.X)+math.pi

        elif self.X==0:
            if self.Y>0:
                self.AngZ=math.pi/2

            elif self.Y<0:
                self.AngZ=-math.pi/2

            elif self.Y==0:
                self.AngZ=0
                
        self.AngZ=self.AngZ-math.pi/2
        # self.AngZ<0:
        #    self.AngZ=-self.AngZ
        
        return self.AngZ
            
        # self.Show()
    
    def FinalPhaseRun(self):
        
        ScaleList=np.arange(8,12,0.01)
        SensativeList=[]
        
        for scale in ScaleList:
            SensativeList.append(self.RunPhaseLock(10,scale))
            
        plt.plot(ScaleList,SensativeList)
        plt.show()
        
            
        #self.Show()
            

        #self.RunPhaseLock(10,10)
        #self.RunPhaseLock(10,100)
    
        
# 新建Mt数组和Beta数组--------------------------------------------------------------------------------
    def beta(self,Times,scale_1,scale_2):
        self.TimeList=np.arange(Times)
        self.BetaList=[]
        self.MtList=[]
        
        for e in self.TimeList:
            self.BetaList.append(self.Omiga(e,scale_1))
            self.MtList.append(self.Mt(e,scale_2))
            
        
    def Omiga(self,Time,Scale):
        return ((Time//Scale)+1)*math.pi 
    
    def Mt(self,Time,Scale):
        
        #return 0.1
        
        return 0.004*math.cos((math.pi/Scale)*Time+math.pi/2)
        
'''   

# X逻辑门操作-----------------------------------------------------------------------------------------
    def RunX(self):
        times_all=100
        time=0
        
        while time<times_all:
            print([self.X,self.Y,self.Z])
            
            self.XList.append(self.X)
            self.YList.append(self.Y)
            self.ZList.append(self.Z)
            
            self.RotateX(self.wZ)
            time+=1
            
        self.Block.add_points([self.XList,self.YList,self.ZList])
            
        self.Show()

# Dephasing自然演化------------------------------------------------------------------------------
    
    def U(self,times_all):
        time_mark=50

        ListX=[]
        ListY=[]
        ListZ=[]
        self.Xdefault=self.X
        self.Ydefault=self.Y
        self.Zdefault=self.Z
        
        self.FAll=0
        for e in range(1):
            
            
            self.X=self.Xdefault
            self.Y=self.Ydefault
            self.Z=self.Zdefault
            
            time=1

            while time<=times_all:
                
                
                self.RotateZ(self.B_Random(0,0.05))
                

                time+=1
                 
            ListX.append(self.X)
            ListY.append(self.Y)
            ListZ.append(self.Z)
            
            # 找到Phi值
            if self.X>0:
                self.AngZ=np.arctan(self.Y/self.X)

            elif self.X<0:
                self.AngZ=np.arctan(self.Y/self.X)+math.pi

            elif self.X==0:
                if self.Y>0:
                    self.AngZ=math.pi/2

                elif self.Y<0:
                    self.AngZ=-math.pi/2

                elif self.Y==0:
                    self.AngZ=0
            
            
            #print(self.AngZ)
            self.F=(1/4)*((1+math.cos(self.AngZ-math.pi/2))**2+math.sin(self.AngZ-math.pi/2)**2)
            #print(self.F)
            self.FAll+=self.F
            
        self.FAll=self.FAll/1000
        
        self.Block.add_points([ListX,ListY,ListZ])
        
        
    def F_T(self):
        
        Tlist=np.arange(100)
        Flist=[]
        for e in Tlist:
            
            self.U(e)
            Flist.append(self.FAll)
            
        
        plt.plot(Tlist,Flist)
        self.Show()


#传统spin echo-----------------------------------------------------------------------------------

    def Spin_X(self):
        times_all=100
        time_mark=50
        time=1
        ListX=[]
        ListY=[]
        ListZ=[]
        self.Xdefault=self.X
        self.Ydefault=self.Y
        self.Zdefault=self.Z
        
        while time<=times_all:
            ListX.append(self.X)
            ListY.append(self.Y)
            ListZ.append(self.Z)
            
            if time==time_mark or time==time_mark*2:
                self.SigX()
                
            else:
                self.RotateY(self.wY)
                self.RotateZ(self.wZ)
                
            time+=1
            
        self.Block.add_points([ListX,ListY,ListZ])
        print("X",self.X-self.Xdefault)
        print("Y",self.Y-self.Ydefault)
        print("Z",self.Z-self.Zdefault)
        self.Show()
        
#新spin echo-------------------------------------------------------------------------------------

    def Spin_NZ(self):
        times_all=200
        time_mark=50
        time=1
        ListX=[]
        ListY=[]
        ListZ=[]
        self.Xdefault=self.X
        self.Ydefault=self.Y
        self.Zdefault=self.Z
        
        while time<=times_all:
            ListX.append(self.X)
            ListY.append(self.Y)
            ListZ.append(self.Z)
            
            if time==time_mark or time==time_mark*3:
                self.SigX()
                
            elif time==time_mark*2 or time==time_mark*4:
                self.SigY()
                
            else:
                self.RotateY(self.wY)
                self.RotateY(self.wY)
                self.RotateZ(self.wZ)
                
            time+=1
            
        self.Block.add_points([ListX,ListY,ListZ])
        print("X",self.X-self.Xdefault)
        print("Y",self.Y-self.Ydefault)
        print("Z",self.Z-self.Zdefault)
        self.Show()
        
'''
AA=BlockClass([-60,30],0,1,0)
AA.FinalPhaseRun()
