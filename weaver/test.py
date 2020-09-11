import numpy as np
import matplotlib.pyplot as plt

from matplotlib import animation, rc
from IPython.display import HTML, Image
from scipy.integrate import odeint


T = 20
dt = 0.005
n = int(T / dt)

def rhs(y,x):
    dydx = -(7*y - 1.)*(1. - y)/y
    return dydx

N = 100
x = np.linspace(0.,6., N)
y0 = 0.999999999            # Change this to y0 = 1. and soln will flatline
y = odeint(rhs, y0, x)

#
# Compare with analytic solution (5.10) - shifted to go through x[nhalf], y[nhalf]
# The factor 42 arises because of the scaling of numerical x above
#

nhalf = int(N/2)
ymin = 1./7. + 0.000001         # avoid log(0) below
ymax = 1. - 0.000001            # avoid log(0) below
yanalytic = np.linspace(ymin, ymax, 1000)
xanalytic = x[nhalf] + (1./42.)*np.log((1 - yanalytic)**7/(7.*yanalytic - 1.)) \
                     - (1./42.)*np.log((1 - y[nhalf])**7/(7.*y[nhalf] - 1.))



rc('animation', html='html5')

fig, ax = plt.subplots()

ax.set_xlim(( 0, 8))
ax.set_ylim(( 0, 1.05))


line, = ax.plot([], [], lw=2)

#Splt.plot(xanalytic,yanalytic)


def update(j):
    x,dx= np.linspace(0,8, 100,retstep=True)
    v = np.linspace(1,1/7,100)
    i = 0
    for i in range(j):
        v[0] = 1
        v[-1] = v[-2]

        deriv_velo = -(v[1:-1]-v[0:-2])/dx
        nonderiv_velo = -((7*v[1:-1]-1)*(1-v[1:-1])/v[1:-1])
        v[1:-1]=v[1:-1]+dt*(deriv_velo+nonderiv_velo)
    line.set_data(x,v)
    return (line,)

anim = animation.FuncAnimation(fig, update,frames=n, interval=0.000001, blit=True)

plt.ylabel('Velocity ')
plt.xlabel('X in Shockframe')
plt.grid()




plt.show()
