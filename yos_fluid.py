#!/usr/bin/env python
# coding: utf-8

# In[2]:


# yangou Du

import numpy as np
get_ipython().run_line_magic('matplotlib', 'notebook')
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

A = np.array([-1,1])
B=np.array([-1,-1])
C=np.array([1,-1])
D=np.array([1,1])
dt=0.001
Set=[A,B,D,C]  # Set contains the four vortices

fig, ax = plt.subplots()

def u(p):
    u=np.array([0,0])
    j=0
    for i in Set:
        r=p-i
        if (j%2==0):  #top two spins one way, x is r rotates by 90Deg
            x1=-r[1]
            x2=r[0]
        else:          #bottom two spins the other
            x1=r[1]
            x2=-r[0] 
        j=j+1
        r2=x1**2+x2**2
        if r2==0 : #not dividing by zeros
            x=0
        else:
            x=np.array([x1/(r2),x2/(r2)])    #k/r k=1 here
        u=u+x
    return(u)


tx = np.arange(-1,11)
ty = np.arange(-2,3)

X, Y = np.meshgrid(tx, ty)

def step(): 
    j=0
    for i in Set:
        i=i+u(i)*dt
        Set[j]=i
        j=j+1
    

def animate(i):
    
    for _ in range(100): #This allows us to interate more than one step per frame
        step()
    
    x= [x for [x,y] in Set]  # x_val for vortixes
    y= [y for [x,y] in Set]  # y_val for vortixes
    
    ax.clear()
    ax.scatter(x,y)      #Plot vortixes

    ax.set_xlim([-1,10])
    ax.set_ylim([-2,2])
    
    ux = np.zeros(shape=X.shape)
    uy = np.zeros(shape=Y.shape)
    (lx, ly) = X.shape
    for i in range(lx): 
        for j in range(ly): 
            ux[i,j], uy[i,j] = u([X[i,j], Y[i,j]])
            
    
    ax.streamplot(X, Y, ux, uy)
    

anim = FuncAnimation(fig, animate, frames=100, interval=20, repeat=False)
def prog(i, n):
    if i % 10 == 0:
        print(f"Saving frame {i} of {n}. ({i/n * 100}%)")
anim.save("TLI.gif", dpi=300, writer=PillowWriter(fps=15), progress_callback=prog)
del(anim)


# In[ ]:




