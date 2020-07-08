import numpy as np
from scipy.integrate import solve_bvp
import matplotlib.pyplot as plt

E_0 = 1 * 0.0000016021773 #erg: gcm^2/s^2
m_H = 1.6*10**(-24) #g
c = 3e11 #cm
sigma_c = 2*10**(-24)
n_0 = 2.03*10**(19)   #1/cm^3
v_0 = (2*E_0/m_H)**(0.5)  #cm/s

def velocity(x, v):
    return np.array([-sigma_c*n_0*v_0*((8*v_0*v[0]-7*v[0]**2-v_0**2)/(2*v[0]*c))])

def bc_velo(ya, yb):
    return np.array([ya[0]-v_0])


n_velo = 800
x_velo = np.linspace(-1e6,3e6, n_velo)
y_velo = np.array([np.linspace(v_0, v_0/7, n_velo)])


sol_velo = solve_bvp(velocity, bc_velo, x_velo, y_velo)

if sol_velo.status != 0:
    print("WARNING: sol.status is %d" % sol_velo.status)
print(sol_velo.message)

def place(v):
    x = (c/(21*sigma_c*n_0*v_0))*(np.log(np.power(v_0-v,7)/((7*v-v_0)*v_0**6))-np.log((1/3)*np.power(3/7,7)))
    return x

v = np.linspace(1/6*v_0,999/1000*v_0,100)#

x = np.zeros(len(v))
vv = np.zeros(len(v))
vn = np.zeros(len(v))
i = 0
for i in range(len(v)):
    x[i] = place(v[i])
    vv[i] = v[i] - v_0/7
print(x)



#-----------------------------------------------------------------------------------------------
T = 10**6
b = 20.3
n_eq = b*T**3
n_s = 2.03*10**(19)




def constants(x):
#    v (x)
    D1 = (c*v/(3*n_0*v_0*sigma_c))
    D2 = ((v**2-8*v_0*v+v_0**2)/(6*v))
    D3 = sigma_c*n_0*v_0*((8*v_0*v-7*v**2-v_0**2)/(2*v*c))
    return D1,D2,D3

def numberdens(x,n):
    D1,D2,D3 = constants(x)
    return np.vstack((n[1],(-D2*n[1]-D3*n[0]+Q*(1-n[0]/n_eq))/(D1)))

def bc_velo(ya, yb):
    return np.array([ya[0]-n_s,yb[0]-n_eq])


n_num = 700
x_num = np.linspace(-4.35*10**(-3),7.65*10**(-4), n_num)
y_num = np.array([np.linspace(n_s, n_eq, n_num)])


#sol_num = solve_bvp(numdens, bc_num, x_num, y_num)


#plt.plot(sol_num.x, sol_num.y[0], label='$n(x)$')
plt.plot(sol_velo.x, sol_velo.y[0]-v_0/7, label='$v(x) num$')
plt.plot(x,vv,label='$v(x)$')
plt.grid(alpha=0.5)
plt.legend(framealpha=1)
plt.show()
