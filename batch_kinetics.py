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

    return -(a/f)*ka*Q*m-ka*m*(m-m0)+kd*(m-m0)

# timespan of solution
t=[0,100] 

# parameters
a=0.5 
f=0.8
m0 = 1 
Q=0.8
ka=0.3
kd=0.01

# solving the IVP
solution = solve_ivp(dmdt,t_span=t,y0=[1],args=[[a,f,m0,Q,ka,kd]],max_step=0.1)
print(solution.y)
pl.plot(solution.t,solution.y[0,:])
pl.grid()
pl.xlabel('time (s)')
pl.ylabel('concentration (mol/m^3)')
pl.show()