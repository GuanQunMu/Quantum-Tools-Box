import matplotlib.pyplot as plt


# For 171Y b+ ions
class UF_Gate:   
    def __init__(self):
        
        # Energy between |0> and |1> state
        self.OmigaHF = 12.6*10**9*3.1415*2

        # Coupling constants (U: up transimition; D: down transimition)
        self.C_P12_U = 1/(3**0.5)
        self.C_P12_D = 1/(3**0.5)
        
        self.C_P32_1U = (2/3)**0.5
        self.C_P32_1D = -(1/6)**0.5
        
        self.C_P32_2D = 1/(2**0.5)
        
        # Spontanrous Emission rate for every state
        self.PE_P12 = 1.238*10**8
        self.PE_P32 = 1.627*10**8
        
        # Satuation intensity
        self.Isat_P12 = 510.3
        self.Isat_P32 = 950.6
        
        # Initializing detuning
        self.Detu_1 = 0
        self.Detu_2 = -100*2*3.1416*10**12
        
        
    def run(self,DetuningList):
        
        # Building empty lists
        self.WavelenghList=[]
        self.STK_ShiftList=[]
        self.EM_RateList=[]
        self.OmigaList=[]
        
        # Start
        for e in DetuningList:
            
            # Scan the detuning 
            D=e[0]                    # Detuning from S1/2 state
            DEnd=e[1]                 # The end of the scanning
            DStep=e[2]                # The step of the scanning
            
            while D<DEnd:

                Wavelengh = 3*10**8/(D*10**12 + 3*10**8/(370*10**(-9)))*10**9

                self.Detu_1 = D*2*3.1416*10**12
                self.Detu_2 = (D-100)*2*3.1416*10**12

                # Rabi Frequency for one detuning
                Omiga = (1/4)*(self.PE_P12**2*self.C_P12_U*self.C_P12_D/(self.Isat_P12*self.Detu_1)+self.PE_P32**2*self.C_P32_1U*self.C_P32_1D/(self.Isat_P32*self.Detu_2))
                
                Delta_L0 = (1/8)*(self.PE_P12**2*self.C_P12_U**2/(self.Isat_P12*self.Detu_1)+self.PE_P32**2*self.C_P32_1U**2/(self.Isat_P32*self.Detu_2))
                Delta_L1 = (1/8)*(self.PE_P12**2*self.C_P12_D**2/(self.Isat_P12*(self.Detu_1-self.OmigaHF))+self.PE_P32**2*self.C_P32_1D**2/(self.Isat_P32*(self.Detu_2-self.OmigaHF))+self.PE_P32**2*self.C_P32_2D**2/(self.Isat_P32*(self.Detu_2-self.OmigaHF)))
                
                if Omiga<0:
                    Omiga=-Omiga
                
                # Spontanrous Emission Probability for one detuning
                SE_Rate = 3.1415*(1/12)*(self.PE_P12**3/(self.Isat_P12*self.Detu_1**2)+2*self.PE_P32**3/(self.Isat_P32*self.Detu_2**2))*3.1415/(2*Omiga)
                
                # Stark Shift for one detuning
                STK_Shift= (Delta_L1-Delta_L0)/Omiga
                
                if SE_Rate<0:
                    SE_Rate=-SE_Rate
                
                
                if STK_Shift<0:
                    STK_Shift=-STK_Shift


                self.WavelenghList.append(Wavelengh)
                self.STK_ShiftList.append(STK_Shift)
                self.EM_RateList.append(SE_Rate)
                self.OmigaList.append(Omiga)
                
                D+=DStep
        
        
        # drawing the figure of EM_probability & STK shift
        fig, ax1 = plt.subplots()

        ax2 = ax1.twinx()   
        #ax1.axis([320,380,0,0.00002])
        ax2.axis([320,380,0,0.002])
        
        
        #ax1.set_ylabel("Spontanrous Emission Probability", color='blue')
        ax2.set_ylabel("Stark Shift Frequency / Rabi Frequency",color='red')
        
        ax2.plot(self.WavelenghList,self.STK_ShiftList,color='red')
        #ax1.plot(self.WavelenghList,self.EM_RateList,color='blue')

        ax1.set_xlabel('nm')
        
        # drawing the figure of Rabi frequency
        plt.show()
        plt.axis([320,380,0,0.07])
        plt.ylabel("Rabi Frequency")
        plt.xlabel("nm")
        plt.plot(self.WavelenghList,self.OmigaList,color='green')
        plt.show()


A=UF_Gate()
A.run([[-500,-3,1],[0.1,90,1],[103,1000,1]])        # Set the scanning detuning list

