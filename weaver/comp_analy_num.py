import numpy as np
import matplotlib.pyplot as plt


N = 50
xleft = 0
xright = 6

x, dx = np.linspace(xleft, xright, N, retstep=True) #retstep yields the

T = 20
dt = 0.0005
n = int(T/dt)





v = np.linspace(1,1/7,N)     #initial guess for the velocity
itcount = 0
itcount2 = 0

for i in range(n):

    v[0] = 1        #the upstream value of the velocity is held constant
    v[-1] = v[-2]   #the last two values of the velocity are required to be identical

    dvdx = -(v[1:-1]-v[0:-2])/dx      #dv/dx is approximated
    nonderiv = -((7*v[1:-1]-1)*(1-v[1:-1])/v[1:-1])  #the non derivative term is implemented

    v[1:-1] = v[1:-1] + dt*(dvdx + nonderiv) #both terms are added and then applied to the velocity using the timestep length
    if i < 1000:
        if itcount2 == 100:
            print('plotted small')
            plt.plot(x,v,color = 'blue')
            itcount2 =0
    if itcount == 6000:
        plt.plot(x,v,color = 'blue')
        itcount = 0
    itcount2 += 1
    itcount += 1

plt.show()

v_left = 1 - 0.00000001
v_right = 1/7 + 0.00000001

v_a = np.linspace(v_left,v_right,10)

x_a  = 1/42 *(np.log((1-v_a)**7/(7*v_a-1))) - 1/42*(np.log(1/3*(3/7)**7))+2



plt.plot(x,v, label = 'Numerical Solution')
plt.scatter(x_a ,v_a, color = 'black', marker = 'x',label = 'Analytical Solution')
plt.xlabel('Distance in Shockframe $[10^6cm]$')
plt.ylabel('Dimensionless Velocity')
plt.xlim(left = 0)
plt.legend()
plt.grid()
plt.show()
