from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp

E_0 = 1 * 0.0000016021773 #erg: gcm^2/s^2
m_H = 1.6*10**(-24) #g
c = 3e11 #cm
sigma_c = 2*10**(-23)
n_0 = 1*10**(20)   #1/cm^3
v_0 = (2*E_0/m_H)**(0.5)  #cm/s

T = 10**7
b = 20.3
n_eq = b*T**3
n_s = 2.03*10**(19)
Q = 1#this is wrong; very clearly

n_num = 1000
v_ana = np.linspace(1/6*v_0,999/1000*v_0,10)
x_num = np.linspace(-1e6,3e6,n_num)

def place(v):
    x = (c/(21*sigma_c*n_0*v_0))*(np.log(np.power(v_0-v,7)/((7*v-v_0)*v_0**6))-np.log((1/3)*np.power(3/7,7)))
    return x


def velocity(v,x):
    dvdx = -sigma_c*n_0*v_0*((8*v_0*v-7*v**2-v_0**2)/(2*v*c))
    return dvdx



x_ana= np.zeros(len(v_ana))

i = 0
for i in range(len(v_ana)):
    x_ana[i] = place(v_ana[i])


sol = odeint(velocity,0.999999999999*v_0,x_num )

sol_new = np.reshape(sol,n_num)

def constants(v):
    D1 = (c*v/(3*n_0*v_0*sigma_c))
    D2 = ((v**2-8*v_0*v+v_0**2)/(6*v))
    D3 = sigma_c*n_0*v_0*((8*v_0*v-7*v**2-v_0**2)/(2*v*c))
    return D1,D2,D3

def numdens(y,x):
    v = np.interp(x,x_num,sol_new)
    D1,D2,D3 = constants(v)
    theta,omega = y
    return [omega,(-D2*omega-D3*theta+Q*((1-theta)/n_eq))/(D1)]

sol_num = odeint(numdens,(n_s,1000),x_num)


plt.scatter(x_ana,v_ana-v_0/7,label='$v(x) analytic$',color = 'red')
plt.plot(x_num,sol-v_0/7,label= '$v(x) numerical$')
#plt.plot(x_num, sol_num[:,0], label='$n(x)$')
#plt.yscale('log')
plt.grid(alpha=0.5)
plt.legend(framealpha=1)
plt.show()
