# radshock.py - programme to solve equation 5.6 of Weaver in the
# sanitised form (dy/dx) = -(7y - 1)*(1 - y)/y
# y = v/v0
# numerical x has absorbed the remaining constants in equation (5.6) - 2
# initial condition is y = 1 at numerical x = 0
#
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

def rhs(y,x):
    dydx = -(7*y - 1.)*(1. - y)/y
    return dydx

N = 1000
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
yanalytic = np.linspace(ymin, ymax, 10)
xanalytic = x[nhalf] + (1./42.)*np.log((1 - yanalytic)**7/(7.*yanalytic - 1.)) \
                     - (1./42.)*np.log((1 - y[nhalf])**7/(7.*y[nhalf] - 1.))


plt.plot(x,y,'-',xanalytic,yanalytic,'o')
plt.show()
