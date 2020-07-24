import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

n = 100000
x = np.linspace(0,13, n)
y0 = 0.9

def velocity(y,x):
    dydx = -(7*y - 1)*(1 -y)/y
    return dydx
i = 0
for i in range(16):
    j = 0
    a= str(0.9)
    for j in range(i):
        a= a + str(9)
        y0 = float(a)

    sol = odeint(velocity,y0,x)
    plt.plot(x,sol)




plt.grid()

plt.show()
