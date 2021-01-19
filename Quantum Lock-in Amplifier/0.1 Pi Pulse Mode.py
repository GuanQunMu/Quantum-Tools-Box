# Quantum Sensing

from qutip import *
import numpy as np
import matplotlib.pyplot as plt
import math

# Hamiltonian
H = sigmax()

# initial state
psi0 = Qobj([[0.707106],[0+0.707106j]])
print(psi0)

# list of times for which the solver should store the state vector

TimeFinal=50
Step=0.01

tlist = np.linspace(0, 50, int(TimeFinal/Step))
period=int(0.1/Step)

TAltered=[tlist[i:i + period] for i in range(0, len(tlist), period)]

delta=1

#def HX_coeff(t, args):
#    return 0 * math.sin(t*2*3.14159)

def HX_coeff(t, args):
    return 0 * math.cos(t*1*3.14159)

def HZ_coeff(t, args):
    return 0.1 * math.cos(t*1*3.14159)

H = [[sigmaz(),HZ_coeff],[sigmax(),HX_coeff]]

RotatedAngle=3.14159/10
XRotatePulse=Qobj([[complex(math.cos(RotatedAngle/2),0),complex(0,-math.sin(RotatedAngle/2))],
                 [complex(0,-math.sin(RotatedAngle/2)),complex(math.cos(RotatedAngle/2),0)]])

ResultFinal=[]
for TLip in TAltered:
    result = mesolve(H, psi0, TLip,[],[])       
    
    for e in result.states:
        ResultFinal.append(e)
        
    psi0=result.states[-1]
    psi0=XRotatePulse*psi0


# plot the expectation value of three Pauli operators
fig, axes = plt.subplots(1,1)

# result.expect[2] denotes the list of state expecation over Pauli-z operator.
axes.plot(tlist, expect(sigmaz(),ResultFinal), label=r'$\left<\sigma_z\right>$')
axes.plot(tlist, expect(sigmax(),ResultFinal), label=r'$\left<\sigma_x\right>$')
axes.plot(tlist, expect(sigmay(),ResultFinal), label=r'$\left<\sigma_y\right>$')
#axes.plot(tlist, result.expect[1], label=r'$\left<\sigma_y\right>$')
#axes.plot(tlist, result.expect[0], label=r'$\left<\sigma_x\right>$')

axes.set_xlabel(r'$t$', fontsize=20)
axes.legend(loc=2);
