# relaxationest.py
# Solve (dy/dt) = -y' + Dy'' + Q*exp(-x^2/L^2) and wait for (dy/dt) = 0
# coding: utf-8

# In[150]:

import numpy as np
import matplotlib.pyplot as plt

#   Define some global parameters
N = 100
D = 9.e-01
Q = 1.
L = 1.e-01
Norm = 1./(L*np.sqrt(np.pi))        #   Gaussian normalised to unity
xleft = -4.
xright = 2.

# In[151]:

# analytic solution on infinite domain the initial guess !!!

x,dx= np.linspace(xleft,xright, N,retstep=True)
yanalytic = np.zeros((x.size))
for i in range(N):
    if x[i] <= 0:
        yanalytic[i] = Q*np.exp(x[i]/D)
    else:
        yanalytic[i] = Q

y = yanalytic[0]*np.ones((x.size))          # Initial y(x,t) = uniform value of yanalytic(xleft)
yleft = y[0]
print('bcs ',yleft)

# Plot
plt.plot(x, yanalytic)
plt.show()

# In[152]:

T = 20  # total time
dt = .0005  # time step
n = int(T / dt)  # number of iterations

print(y)
print(n)

#y=np.ones((x.size))

for i in range(n):
    y[0]=yleft                              # b.c. at xleft, y = yanalytic(xleft)
    y[-1]=y[-2]                             # b.c. at xright, y' = 0

    adv=-(y[1:-1]-y[0:-2])/dx
    diff=D*(y[2:]-2.*y[1:-1]+y[0:-2])/(dx*dx)
    reac=Norm*Q*np.exp(-(x[1:-1]/L)**2)
    y[1:-1]=y[1:-1]+dt*(adv+diff+reac)

# In[153]:

# Plot solutions

plt.plot(x, y, label='y0')
plt.plot(x, yanalytic, label='exact')
plt.legend()
plt.show()
