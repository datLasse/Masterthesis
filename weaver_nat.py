import numpy as np
from scipy.integrate import solve_bvp,odeint
import matplotlib.pyplot as plt

E_0 = 1 * 0.0000016021773 #erg: gcm^2/s^2
m_H = 1.6*10**(-24) #g
c = 3e11 #cm
sigma_c = 2*10**(-23)
n_0 = 1*10**(20)   #1/cm^3
v_0 = (2*E_0/m_H)**(0.5)  #cm/s
k = 1.38*10**(-16)
eps_c = 8.01*10**(-9)

#def velocity(x, v):
#    return np.array([-sigma_c*n_0*v_0*((8*v_0*v[0]-7*v[0]**2-v_0**2)/(2*v[0]*c))])

#def bc_velo(ya, yb):
#    return np.array([ya[0]-v_0])


#n_velo = 800
#x_velo = np.linspace(-1e6,3e6, n_velo)
#y_velo = np.array([np.linspace(v_0, v_0/7, n_velo)])


#sol_velo = solve_bvp(velocity, bc_velo, x_velo, y_velo)

#if sol_velo.status != 0:
#    print("WARNING: sol.status is %d" % sol_velo.status)
#print(sol_velo.message)

#def place(v):
#    x = (c/(21*sigma_c*n_0*v_0))*(np.log(np.power(v_0-v,7)/((7*v-v_0)*v_0**6))-np.log((1/3)*np.power(3/7,7)))
#    return x

#v = np.linspace(1/6*v_0,999/1000*v_0,100)

#x = np.zeros(len(v))
#vv = np.zeros(len(v))
#vn = np.zeros(len(v))
#i = 0
#for i in range(len(v)):
#    x[i] = place(v[i])
#    vv[i] = v[i] - v_0/7



#-----------------------------------------------------------------------------------------------
#everything needs to be natural
T = 10**7
b = 20.3
n_eq = b*T**3/n_0
n_s = 2.03*10**(19)/n_0
#Q = 100#this is wrong; very clearly

print(n_s)
def velocity(v,x):
    dvdx = -((8*v-7*v**2-1**2)/(v))
    return dvdx


n_num = 400
x_num = np.linspace(0,10, n_num)

sol_velo = odeint(velocity,0.999999,x_num)

sol_new = np.reshape(sol_velo,n_num)



def constants(v):
    n = (n_0)/(v)
    D1 = c*v/(3*n_0*sigma_c)
    D2 = v_0*n_0*((v**2-8*v+1)*v_0/(6*v))
    D3 = n_0**2*v_0**2**sigma_c*(((7*v-1)*(1-v)))/(v*2*c)
    return D1,D2,D3,n

def numdens(x,y):
    v = np.interp(x,x_num,sol_new)
    D1,D2,D3,n = constants(v)
    theta_e = (10/9)*v_0**2*m_H*(1-v)/y[0]
    T_e= theta_e/k
    lamb = eps_c/theta_e
    g = 1.226 - 0.475*np.log(lamb) + 0.0013*(np.log(lamb))**2
    E = -np.log(lamb) - 0.5772
    Q = 5.692*10**(-12)*T_e**(-0.5)*n**2*g*E
    n_eq = b*T_e**3
    #print(y[0])
    return np.vstack((y[1],(-D2*y[1]-D3*y[0]-Q*(1-y[0]/n_eq))/(D1)))

def numdens1(y,x):
    v = np.interp(x,x_num,sol_new)
    D1,D2,D3,n = constants(v)
    theta,omega = y
    return [omega,(-D2*omega-D3*theta+Q*((1-theta)/n_eq))/(D1)]

def bc_num1(ya, yb):
    return np.array([ya[0]-n_s,ya[1]])

def bc_num2(ya, yb):
    return np.array([ya[0]-n_s,yb[0]-n_eq])

def bc_num3(ya, yb):
    return np.array([ya[0]-n_s,yb[1]])


y_num = np.array([np.linspace(n_s, n_eq, n_num),np.linspace(1000, 0, n_num)])

#y0= [1,-Q/v_0]


#sol_num1 = solve_bvp(numdens, bc_num1, x_num, y_num)
#sol_num2 = solve_bvp(numdens, bc_num2, x_num, y_num)
sol_num3 = solve_bvp(numdens, bc_num3, x_num, y_num)
#sol_ode = odeint(numdens1,y0,x_num)

#if sol_num.status != 0:
#    print("WARNING: sol.status is %d" % sol_velo.status)
#print(sol_num.message)
#temp_e = np.zeros(len(sol_num.y[0]))
#neq = np.zeros(len(sol_num.y[0]))
#i = 0
#for i in range(len(diff1)):
#    temp_e[i] = (10/9)*v_0**2*m_H*(1-sol_new[i])/(sol_num.y[0][i]*k)
#    neq[i] = n_eq = b*temp_e[i]**3


#plt.plot(sol_num1.x, sol_num1.y[0], label='n(x) boundary condition only upstream')
#plt.plot(sol_num2.x, sol_num2.y[0], label='n(x) boundary condition up- and downstream $n_{eq}$')
plt.plot(sol_num3.x, sol_num3.y[0], label='$n(x) $')
plt.plot(x_num, sol_velo, label='$v(x)$')
#plt.plot(x_num,temp_e, label ='T(x)')
#plt.plot(x_num,neq,label='$n_{eq}$')
#plt.plot(x_num,sol_ode[:,0], label='n(x) with odeint')
plt.yscale('log')
plt.grid(alpha=0.5)
plt.legend(framealpha=1)
plt.show()
