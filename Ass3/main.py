import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pyplot as pl
# I made changes from the professors sample code, so it looks similar to it.
# Set up the grid and advection and diffusion parameters
Ngrid = 50
Nsteps = 5000
dt = 1 #(cs)
dx = 1 #(dm)
# our equation df/dt+u'(df/dx)=g+D1(d^2u/dx^2)
#difusion rate
D1 = 2.5# this  ratio is explained in the pdf
#effective gravity sin(a)*g
g = 0.017

beta1 = D1*dt/dx**2


x = np.arange(0, Ngrid*1., dx) / Ngrid # multiplying by 1. to make sure this is an array of floats not integers

# analytic solution: u(x) = -g/D1*(x^2/2-H*x) ,H=1dm
u =-g/D1*(x**2/2-1*x)*2350 # I scale it by 2350 to have the right height. I am not sure why the units didn't work out
f1 = 10*x#10 to make it have a closer starting value for it to run shorter, changing it should not have a different cruve in the end





# Set up plot
pl.ion()
fig, axes = pl.subplots(1,1)
axes.set_title('D=%3.1f'%D1)


# Plotting initial and final state in the background for reference
axes.plot(x, f1, 'k-')
axes.plot(x, u, 'k-')


# We will be updating these plotting objects
plt1, = axes.plot(x, f1, 'ro')


# Setting the axis limits for visualization

axes.set_xlim([0,1])
axes.set_ylim([0,10])

# this draws the objects on the plot
fig.canvas.draw()

for ct in range(Nsteps):

    ## Calculate diffusion first
    # Setting up matrices for diffusion operator
    A1 = np.eye(Ngrid) * (1.0 + 2.0 * beta1) + np.eye(Ngrid, k=1) * -beta1 + np.eye(Ngrid, k=-1) * -beta1

    ## Boundary conditions to keep the first element fixed
    # This ensures f in the first cell stays fixed at all times under diffusion
    A1[0][0] = 1.0

    A1[0][1] = 0

    # Stress-free boundary condition on the right
    A1[Ngrid - 1][Ngrid - 1] = 1.0 + beta1

    # Solving for the next timestep
    f1 = np.linalg.solve(A1, f1)

    ## Calculate accelaration due to g
    f1[1:Ngrid - 1]=f1[1:Ngrid - 1]+g*dt

    # update the plot
    plt1.set_ydata(f1)


    fig.canvas.draw()
    pl.pause(0.0001)

