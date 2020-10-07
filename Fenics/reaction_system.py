from fenics import *

T = 10            # final time
num_steps = 1000  # number of time steps
dt = T / num_steps  # time step size

# Read mesh from file
mesh = RectangleMesh(Point(0, 0), Point(1, 1), 20, 20)

# Define function space for system of concentrations
P1 = FiniteElement('P', triangle, 1)
element = MixedElement([P1, P1])
V = FunctionSpace(mesh, element)
W=FunctionSpace(mesh, P1)

# Define test functions
v_1, v_2 = TestFunctions(V)

# Define functions for velocity and concentrations
w = Constant((0, 0.1))
u = Function(V)  # metal concentration
u_n = Function(V)  # concentration at timestep n

# Split system functions to access components
u_1, u_2 = split(u)
u_n1, u_n2 = split(u_n)

# Define source terms
f_10 = Expression('pow(x[1],2)<1 ? 0.015737*30 : 0', degree=1)
f_1 = interpolate(f_10, W)

f_2 = Constant(0.35e-6)

# Define expressions used in variational forms
k = Constant(dt)
a = Constant((33.4 * 1.5e6)/(1-0.31))
f = Constant(0.31)
m0 = Constant(0.015737*30)
Q = Constant(0.35e-6)
ka = Constant(0.015)
KD = Constant(10**(-13.7))
kd = Constant(ka*KD)
D = Constant(1.297e-9)

# Define variational problem
F = ((u_1 - u_n1) / k)*v_1*dx - dot(w, grad(u_1))*v_1*dx \
  + D*dot(grad(u_1), grad(v_1))*dx  \
  + ((u_2 - u_n2) / k)*v_2*dx +\
  -(a/f)*(-ka*Q*u_1+ka*u_1*u_2+kd*u_2)*v_1*dx \
  + (-ka*Q*u_1+ka*u_1*u_2+kd*u_2)*v_2*dx \
  - f_1*v_1*dx-f_2*v_2*dx
# Create VTK files for visualization output
vtkfile_u_1 = File('reaction_system/u_1.pvd')
vtkfile_u_2 = File('reaction_system/u_2.pvd')

# Time-stepping
t = 0
for n in range(num_steps):

    # Update current time
    t += dt

    # Solve variational problem for time step
    solve(F == 0, u, solver_parameters={"newton_solver": {"relative_tolerance": 1e-8, "maximum_iterations": 500}})

    # Save solution to file (VTK)
    _u_1, _u_2 = u.split()
    vtkfile_u_1 << (_u_1, t)
    vtkfile_u_2 << (_u_2, t)

    # Update previous solution
    u_n.assign(u)
