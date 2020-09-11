import numpy as np
import matplotlib.pyplot as plt

#----------------GLOBAL PARAMETERS-----------------------------------------------------------------


N = 100

xleft = 0
xright = 6

x,dx= np.linspace(xleft,xright, N,retstep=True)

T = 20
dt = 0.0005
n = int(T / dt)

#----------SOLUTION FOR VELOCITY USING STEADY STATE METHOD---------------------------------------

v = np.linspace(1,1/7,N)


for i in range(n):
    v[0] = 1
    v[-1] = v[-2]

    A = -(v[1:-1]-v[0:-2])/dx
    B = -((7*v[1:-1]-1)*(1-v[1:-1])/v[1:-1])
    v[1:-1]=v[1:-1]+dt*(A+B)

#-------------------PLOTTING------------------------------------------------

plt.plot(x, v, label='$v(x)$')
plt.grid()
plt.legend()
plt.show()
