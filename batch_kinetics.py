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

def q(t,m,params):
    f = params[0]
    a=params[1]
    Q=params[2]
    return ((f/a)*(m0-m))/Q

# timespan of solution
t=[0,5] 

# parameters
a=(33.4 * 1.5e6)/(1-0.31) 

# specific area = how much surface area per volume
# SA_density * mass_density / volume = m^2 / g * g/m^3 * m^3/m^3 = 1/m

f=0.31
m0 = 0.015*30
Q= 0.35e-6
ka=0.015
KD=(10**(-13.7))
kd=ka*KD

# fig,ax=pl.subplots()
# for c in np.arange(-0.3,0.3, step=0.05):
# # solving the IVP

solution = solve_ivp(dmdt,t_span=t,y0=[m0],args=[[a,f,m0,Q,ka,kd]],max_step=0.001)

# ax.plot(solution.t,solution.y[0,:],label="[Cu] remaining, [CopC]=0.35e-6")
# ax.set_title('Concentration of copper versus time for various binding protein concentrations')
# ax.grid()
# ax.legend()
# ax.set_xlabel(r'time $(s)$')
# ax.set_ylabel(r'concentration $(mol/m^3)$')
# pl.savefig('batch_reactor.png')
# pl.show()

fig,ax=pl.subplots()
ax.plot(solution.t,100*q(0,solution.y[0,:],[f,a,Q]),label="percentage [CopC] occupied")
ax.set_title('Proportion of occupied CopC versus time')
ax.grid()
ax.legend()
ax.set_xlabel(r'time $(s)$')
ax.set_ylabel(r'proportion (%)')
pl.savefig('batch_reactor.png')
pl.show()