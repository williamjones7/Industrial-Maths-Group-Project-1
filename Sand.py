import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation
%matplotlib inline

#0 dirichlet B.C's - no sand at either end of the belt
dir0 = 0.0  # u(0,t)= dir0
dir1 = 0.0  # u(1,t)= dir1

#  function for setting initial conditions in space
def I(x): # initial u(x,0) = 0, dump sand everywhere
    len_x = np.size(x)
    i_x = np.zeros(len_x)
    return i_x

#Function for exact solution, when we find it
# def U_exact(x,t):
#     M = np.size(x)
#     u_ex = np.zeros(M)  
#     u_ex[0] = dir0   # a1 is dirichlet boundary condition at u(0,t)
#     for i in range(1,M-1):
#         sum_u_ex = 0
#         for n in range(1,2000):
#             sum_u_ex += (((2*(-1)**n)/(np.pi*n))+(n%2)*(8/(np.pi*n)))*np.sin(n*np.pi*x[i])*np.exp((-1)*(n**2)*(np.pi**2)*t)
#         u_ex[i] = sum_u_ex
            
        
#     u_ex = u_ex + x  # add x for term in exact solution
#     u_ex[M-1] = dir1 # dirichlet boundary condition at u(L,t)
#     return u_ex

#plot function
def plotting(U,L):
    size_t = np.size(U,1)-1 # number of time steps
   
    fig, ax = plt.subplots(1, 1, figsize=(10, 8))
    markers =['X','.','+','o']

    colours = ['r','g','b','purple','yellow'] # make comparison easy
    colour_pos = 0;

    N_dots = 20
    x_dots = np.linspace(0, L, N_dots+1)    # spacial points to plot exact solution at

# some arbitarily chosen time points to compare at
    for plot_pos in [0, int(size_t/10),int(size_t/5), int(size_t/2), int(size_t)]:
    
    # exact solution
        # U_tplot = U_exact(x_dots,t[plot_pos]) 
        # label = "Exact, t=" + "%0.3f" % (t[plot_pos],)
        # ax.plot(x_dots,U_tplot,linestyle = ':',color = colours[colour_pos],marker = markers[0], label=label)

    # numerical solution
        label = "Numerical PDE, t=" + "%0.3f" % (t[plot_pos],)
        ax.plot(x_pde,U[:,plot_pos],color = colours[colour_pos],ls='--', label=label)
        colour_pos = colour_pos + 1;
    

    plt.xlim(0,L) # zoom in on area of interest
    ax.legend() # turn on legend 
    plt.show()



## ====== setting parameters =======

Nt_gaps = 30000    # number of timesteps
T = 1           # final time 
#h = t_max/(Nt_points)  # time step

Nt_points = Nt_gaps + 1

t = np.linspace(0.,T,Nt_points)  # times at each time step
Nx_spaces = 100; # number of spaces in x direction
Nx_points = Nx_spaces + 1 
L = 1; 
x_pde = np.linspace(0, L, Nx_points)    # mesh points in space
dx = x_pde[1] - x_pde[0] 
dt = t[1] - t[0]
C = dt/(dx**2)
print("Delta x =", dx, "Delta t = ", dt, "C =", C)


##----------intialise lists ------------------------

# Define the numerical solution 
# the first index is space and the second time
U = np.zeros((Nx_points,Nt_points))

# The initial condition

U[:,0]=I(x_pde)
# enforcing the boundary condition
U[0,0]  = dir0 
U[-1,0] = dir1

u_old = I(x_pde)

u = np.zeros(Nx_points)
# and to store the full solution
U = np.zeros((Nx_points,Nt_points))
U[:,0] = u_old

#compute numerical solution 
for n in range(1, Nt_points):
    
    # set Dirichlet boundary points here
    u[0] = dir0
    u[-1] = dir1
    
    # compute u at inner mesh points
    for i in range(1, Nx_points-1):
        u[i] = u_old[i] + (dt/(dx**2))*(u[i-1]-2*u[i]+u[i+1])
 
    # update u_old before next step
    u_old[:]= u

    # copy into full storage
    U[:,n] = u;
    
plotting(U, L)