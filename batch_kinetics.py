# This is a script for visualizing the solutions to the kinetic equations for the batch process.

import matplotlib.pyplot as pl # plotting
import numpy as np # numerics
from scipy.integrate import solve_ivp # differential equations

# Returns reaction rate as function of concentration
def dmdt(t,m,params):
    a=params[0] # specific surface area
    f=params[1] # void volume fraction
    m0=params[2] # initial concentration of metal
    Q=params[3] # total available binding surface
    ka=params[4] # forward reaction rate
    kd=params[5] # backwards reaction rate

    return -(a/f)*ka*Q*m+ka*m*(m0-m)+kd*(m0-m)

# timespan of solution
t=[0,86400/24] 

# parameters
a=3/(3e-3)
f=0.31
m0 = 0.015737*30
Q= (f/a)*0.015737*30
ka=0.015
KD=(10**(-13.7))
kd=ka*KD

fig,ax=pl.subplots()
for c in np.arange(-0.3,0.3, step=0.05):
# solving the IVP
    solution = solve_ivp(dmdt,t_span=t,y0=[m0],args=[[a,f,m0,Q+c*(f/a),ka,kd]],max_step=100)
    if abs(c)<=1e-10:
        ax.plot(solution.t,solution.y[0,:],color=(0,0,0),lw=2,label='[CupC] = '+'%4f'%(Q+c*(f/a))+' mol/m^2')
    else:
        ax.plot(solution.t,solution.y[0,:],lw=1,color=(0,0,0),ls='dotted',label='[CupC] = '+'%4f'%(Q+c*(f/a))+' mol/m^2')
ax.set_title('Concentration of copper versus time for various binding protein concentrations')
ax.grid()
ax.legend()
ax.set_xlabel(r'time $(s)$')
ax.set_ylabel(r'concentration $(mol/m^3)$')
pl.show()