import numpy as np
from scipy.integrate import solve_bvp,odeint
import matplotlib.pyplot as plt


#------------CONSTANTS-----------------------------------------------------------------------

E_0 = 1 * 0.0000016021773 #erg: gcm^2/s^2
m_H = 1.6*10**(-24) #g
c = 3e11 #cm

n_0 = 1*10**(20)   #1/cm^3
v_0 = (2*E_0/m_H)**(0.5)  #cm/s
k = 1.38*10**(-16)
eps_c = 8.01*10**(-9)

T = 10**6
b = 20.3
n_eq = b*T**3/n_0
n_s = 2.03*10**(19)/n_0
sigma_c = 6.65*10**(-25)
#sigma_c = 3*T*k*sigma_t


n_num = 1000
x_num = np.linspace(0,6, n_num)

#------------SOLUTION FOR VELOCITY USING ODEINT----------------------------------------------


def velocity(v,x):
    dvdx = -((7*v-1)*(1-v)/(v))
    return dvdx

sol_velo = odeint(velocity,0.999,x_num)



#-------------SOLUTION FOR NUMBER DENSITY USING SOLVE_BVP-----------------------------------

sol_new = np.reshape(sol_velo,n_num)

def numdens(x,y):
    v = np.interp(x,x_num,sol_new)
    theta_e = (10/9)*v_0**2*m_H*(1-v)/y[0]
    #sigma_c = 3*theta_e*sigma_t
    n = (n_0)/(v)
    D1 = c*v/(3*n_0*sigma_c)
    D2 = v_0*n_0*((v**2-8*v+1)*v_0/(6*v))
    D3 = n_0**2*v_0**2**sigma_c*(((7*v-1)*(1-v)))/(v*2*c)
    T_e= theta_e/k
    lamb = eps_c/theta_e
    g = 1.226 - 0.475*np.log(lamb) + 0.0013*(np.log(lamb))**2
    E = -np.log(lamb) - 0.5772
    Q = 5.692*10**(-12)*T_e**(-0.5)*n**2*g*E
    n_eq = b*T_e**3
    return np.vstack((y[1],(-D2*y[1]-D3*y[0]-Q*(1-y[0]/n_eq))/(D1)))


def bc_num(ya, yb):
    return np.array([ya[0]-n_s,ya[1]])


y_num = np.array([np.linspace(n_s, 10000*n_eq, n_num),np.linspace(0, 0, n_num)])


sol_num = solve_bvp(numdens, bc_num, x_num, y_num)


#-------------SOLUTION FOR TEMPERATURE USING NUMBERDENSITY AND VELOCITY -------------------------

temp_e = np.zeros(len(sol_num.y[0]))

i = 0

for i in range(len(sol_num.y[0])):
    temp_e[i] = (10/9)*v_0**2*m_H*(1-sol_new[i])/(sol_num.y[0][i]*k)


#--------------PLOTTING SOLUTIONS ------------------------------------------------------------

plt.plot(sol_num.x, sol_num.y[0]*n_0, label='$n(x) $')
plt.plot(x_num, sol_velo, label='$v(x)$')
plt.plot(x_num,temp_e, label ='T(x)')

plt.yscale('log')
plt.grid()
plt.legend(framealpha=1)
plt.show()
