import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(0,1,100)

rho = (1-x)**3

v = rho**(-0.19)

plt.plot(x,rho, label =  '$rho$')
plt.plot(x,v, label = 'v')
plt.legend()
plt.show()
