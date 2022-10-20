import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation

## ====== setting parameters =======

Nt_gaps = 30000    # number of timesteps
T = 30           # final time - we're interested in time as it gets large  
Nt_points = Nt_gaps + 1

t = np.linspace(0.,T,Nt_points)  # times at each time step
Nx_spaces = 100; #number of spaces in x direction
Nx_points = Nx_spaces + 1 

#0 dirichlet B.C's - no sand at either end of the belt
dir0 = 0.0  # u(0,t)= dir0
dir1 = 0.0  # u(L,t)= dir1

length = 5 #default conveyor length
# #Set up spacial points for exact solution
N_dots = 50

## ===== defining our functions =====

#Allows us to change the length if needed 
def scale(L):
    x_pde = np.linspace(0, L, Nx_points) #mesh points in space
    dx = x_pde[1] - x_pde[0]
    dt = t[1] - t[0]
    return x_pde, dx, dt

#Function for setting initial conditions in space
def I(x): # initial u(x,0) = 0, dump sand at all x, assume we start with no sand
    len_x = np.size(x)
    i_x = np.zeros(len_x)
    return i_x

# #Function for exact solution
def exact(s, V, D, L):
    x = np.linspace(0, L, N_dots+1)
    M = np.size(x)
    u_ex = np.zeros(M) 
    for i in range(1,M-1):
        beta = (s*L)/(V*(1-np.exp((V*L)/D)))
        u_ex[i] = beta*(np.exp((V/D)*x[i])-1) + (s/V)*x[i]
    return u_ex

#Define the numerical solution for different belt speeds, and coefficients
#Define external constants V (speed of belt), D(Diffusion coefficient), and s (source term)
def numerical(s, V, D, L):
    x_pde, dx, dt = scale(L)

    #The first index is space and the second time
    U = np.zeros((Nx_points,Nt_points))

    #The initial condition
    U[:,0]=I(x_pde)

    #Boundary conditions
    U[0,0]  = dir0 
    U[-1,0] = dir1

    #Find coefficients for numerical solution 
    p = V*dt/dx
    r = D*dt/(dx**2)

    #For stability we require p <= 1 and r <= 1/2
    print("Delta x =", dx, "Delta t = ", dt, "p =", p, "r =", r)

    u_old = I(x_pde)

    u = np.zeros(Nx_points)
    # and to store the full solution
    U = np.zeros((Nx_points,Nt_points))
    U[:,0] = u_old

    #compute numerical solution 
    for n in range(1, Nt_points):
        
        #set Dirichlet boundary points here
        u[0] = dir0
        u[-1] = dir1
        
        #compute u at inner mesh points
        for i in range(1, Nx_points-1):
            u[i] = (p/2 + r)*u_old[i-1] + (1 - 2*r)*u_old[i] + (r - p/2)*u_old[i+1] + s*dt
    
        #update u_old before next step
        u_old[:]= u.copy()

        #copy into full storage
        U[:,n] = u.copy()

    return U


# ===== Plotting =====

# markers =['x', '+', 'o']
# linestyle = ['dotted', 'dashed', 'dashdot', 'solid']
# colours = ['r','g','b','purple','yellow','orange'] # make comparison easy

# x_pde, dx, dt = scale(length)
# x_dots = np.linspace(0, length, N_dots+1)

# plt.figure(1, figsize=(10,8), dpi =200)
# advection = [0.3, 0.4, 0.5]
# diffusion = [0.2, 0.6, 1.0]
# lengths = [5,6,7]

# for i in range(0, len(advection)):
#     label = "Exact, V=" + "%0.3f" % (advection[i])
#     plt.plot(x_dots, exact(1.0, advection[i], 0.5, length),color = colours[0], marker = markers[i], label = label)
#     label = "Numerical PDE, V=" + "%0.3f" % (advection[i])
#     plt.plot(x_pde, numerical(1.0, advection[i], 0.5, length)[:,-1], color = colours[0], linestyle = linestyle[i], label = label)

#     label = "Exact, D=" + "%0.3f" % (diffusion[i])
#     plt.plot(x_dots, exact(1.0, 0.25, diffusion[i], 5),color = colours[1], marker = markers[i], label = label)
#     label = "Numerical PDE, D=" + "%0.3f" % (diffusion[i])
#     plt.plot(x_pde, numerical(1.0, 0.25, diffusion[i], 5)[:,-1], color = colours[1], linestyle = linestyle[i], label = label)

#     x_pde_2, dx_2, dt_2 = scale(lengths[i])
#     x_dots_2 = np.linspace(0, lengths[i], N_dots+1)  
#     label = "Exact, L=" + "%0.3f" % (lengths[i])
#     plt.plot(x_dots_2, exact(1.0, 0.25, 0.5, lengths[i]),color = colours[2], marker = markers[i], label = label)
#     label = "Numerical PDE, L=" + "%0.3f" % (lengths[i])
#     plt.plot(x_pde_2, numerical(1.0, 0.25, 0.5, lengths[i])[:,-1], color = colours[2], linestyle = 'dashed', label = label)

# plt.legend()
# plt.xlim(0,max(lengths)) # zoom in on area of interest
# plt.xlabel('Position on conveyor belt (metres)', fontsize = 15)
# plt.ylabel('Height of sand (metres)', fontsize = 15)
# plt.savefig('Figure 1.png')


plt.figure(2, figsize=(10,8), dpi =200)
Nintervals = 50

max_height_adv = []
interval_size = (0.5-0.25)/Nintervals
for v in np.arange(0.25, 0.25 + interval_size*(Nintervals+1), interval_size):
    max_height_adv.append(max(exact(1.0, v, 0.5, length)))
plt.plot(max_height_adv, label = 'Advection')

max_height_diff = []
for d in np.arange(0.04, 0.04 + interval_size*(Nintervals+1), interval_size):
    max_height_diff.append(max(exact(1.0, 0.25, d, length)))
plt.plot(max_height_diff, label = 'Diffusion')

max_height_length = []
for l in np.arange(5, 5  + interval_size*(Nintervals+1), interval_size):
    max_height_length.append(max(exact(1.0, 0.25, 0.5, l)))
plt.plot(max_height_length, label = 'Length')

plt.legend()
plt.xlabel('Change in variable', fontsize = 15)
plt.ylabel('Max height', fontsize = 15)
plt.savefig('Figure 2.png')


# plt.figure(1, figsize=(10, 8), dpi=200)
# #define advection speeds
# advection = [0.25, 0.3, 0.35, 0.4, 0.45, 0.5]
# for i in range(0, len(advection)):
#     label = "Exact, V=" + "%0.3f" % (advection[i])
#     plt.plot(x_dots, exact(x_dots, 1.0, advection[i], 0.5, 5),color = colours[i], linestyle = ':', marker = markers[0], label = label)
#     label = "Numerical PDE, V=" + "%0.3f" % (advection[i])
#     plt.plot(x_pde, numerical(1.0, advection[i], 0.5, 5)[:,-1], color = colours[i], linestyle = 'dashed', label = label)

# plt.legend()
# plt.xlim(0,length) # zoom in on area of interest
# plt.xlabel('Position on conveyor belt (metres)', fontsize = 15)
# plt.ylabel('Height of sand (metres)', fontsize = 15)
# plt.savefig('Advection Figure 1.png')

# plt.figure(2, figsize=(10, 8), dpi=200)
# #define diffusion coefficients 
# diffusion = [0.2, 0.4, 0.6, 0.8, 1.0]
# for i in range(0, len(diffusion)):
#     label = "Exact, D=" + "%0.3f" % (diffusion[i])
#     plt.plot(x_dots, exact(x_dots, 1.0, 0.25, diffusion[i], 5),color = colours[i], linestyle = ':', marker = markers[1], label = label)
#     label = "Numerical PDE, D=" + "%0.3f" % (diffusion[i])
#     plt.plot(x_pde, numerical(1.0, 0.25, diffusion[i], 5)[:,-1], color = colours[i], linestyle = 'dashed', label = label)

# plt.legend()
# plt.xlim(0,5) # zoom in on area of interest
# plt.xlabel('Position on conveyor belt (metres)', fontsize = 15)
# plt.ylabel('Height of sand (metres)', fontsize = 15)
# plt.savefig('Diffusion Figure 1.png')

# plt.figure(3, figsize=(10, 8), dpi=200)
# #define lengths
# lengths = [5,6,7,8,9,10]
# for i in range(0, len(lengths)):
#     x_pde, dx, dt = scale(lengths[i])
#     x_dots = np.linspace(0, lengths[i], N_dots+1)  
#     label = "Exact, L=" + "%0.3f" % (lengths[i])
#     plt.plot(x_dots, exact(x_dots, 1.0, 0.25, 0.5, lengths[i]),color = colours[i], linestyle = ':', marker = markers[1], label = label)
#     label = "Numerical PDE, L=" + "%0.3f" % (lengths[i])
#     plt.plot(x_pde, numerical(1.0, 0.25, 0.5, lengths[i])[:,-1], color = colours[i], linestyle = 'dashed', label = label)

# plt.legend()
# plt.xlim(0,10) # zoom in on area of interest
# plt.xlabel('Position on conveyor belt (metres)', fontsize = 15)
# plt.ylabel('Height of sand (metres)', fontsize = 15)
# plt.savefig('Length Figure 1.png')

