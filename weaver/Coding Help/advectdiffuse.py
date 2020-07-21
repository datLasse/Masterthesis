# advectdiffuse.py - programme to solve advection/diffusion eqn as a boundary value and initial value problem 
#   y' = Dy'' + Q*exp(-x^2/L^2)
#   Flow velocity constant from xleft to xright
#   D = diffusion coefficient, divided by constant flow velocity
#   Q = injection term, divided by constant velocity, and unit-normalised Gaussian
#   As L goes to zero, and infinite domain : (i) yanalytic = Q*exp(x/D) (x<0) (ii) yanalytic = Q (x>0)
#   Boundary values for bvp code are (normally) (i) y(xleft) = Q*exp(xleft/D) and (ii) y'(xright) = 0
#   Initial solution for bvp code is y = yanalytic
#   Initial conditions for ode code are y(left) and y'(left) given by analytic soln 
#
import numpy as np
import matplotlib.pyplot as plt 
from scipy.integrate import solve_bvp
from scipy.integrate import odeint

#   Define some global parameters

N = 100
D = 9.e-01
Q = 1.
L = 1.e-05
Norm = 2./(L*np.sqrt(np.pi))        #   Gaussian normalised to unity 
xleft = -4.
xright = 2.

#   Three functions : fun and bc for bvp solver; and derivs for odeint solver 

def fun(x,y):
    """ Returns RHS of 2 diff eqns in normalised units - for bvp 
        y[0] = n
        y[1] = n'
        x = position
    """
    rhs1 = y[1]
    rhs2 = D**(-1.)*(y[0] - Q*Norm*np.exp(-(x/L)**2))
    return np.array([rhs1, rhs2])

def bc(ya,yb):
    """ Boundary conditions (4) in terms of residual values
        a = left hand boundary
        b = right hand boundary 
    """
    nleft = Q*np.exp(xleft/D)
    return np.array([ya[0] - nleft, yb[1]])

def derivs(s,x):
    """ Returns RHS of 2 diff eqns in normalised units - for odeint
        s[0] = n
        s[1] = n'
        x = position
    """
    rhs1 = s[1]
    rhs2 = D**(-1.)*(s[0] - Q*Norm*np.exp(-(x/L)**2))
    deriv = np.array([rhs1, rhs2])
    return deriv

# Solve as a boundary value problem with analytic solution on infinite domain the initial guess !!!

x = np.linspace(xleft,xright, N)
y = np.zeros((2, x.size))
for i in range(N):
    if x[i] <= 0:
        y[0,i] = Q*np.exp(x[i]/D)
        y[1,i] = (Q/D)*np.exp(x[i]/D)
    else:
        y[0,i] = Q
        y[1,i] = 0.

yanalytic = y[0,]    
sol = solve_bvp(fun, bc, x, y)

# Solve as an initial value problem

dens_init = Q*np.exp(xleft/D)
grad_init = (Q/D)*np.exp(xleft/D) 
state = np.array([dens_init, grad_init])

solution = odeint(derivs, state, x) 

# Plot solutions 

x_plot = x
density = sol.sol(x_plot)[0]                                #   y from bvp
gradient = sol.sol(x_plot)[1]                               #   y' from bvp
density_ode = np.transpose(np.array([solution[:,0]]))       #   y from odeint
gradient_ode = np.array([solution[:,1]])                    #   y' from odeint 
plt.plot(x_plot, density, x_plot, yanalytic, x_plot, density_ode)
plt.show()
