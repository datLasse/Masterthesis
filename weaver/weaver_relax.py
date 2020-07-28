import numpy as np
import matplotlib.pyplot as plt

#----------------GLOBAL PARAMETERS-----------------------------------------------------------------

N = 100

xleft = 0
xright = 6

x,dx= np.linspace(xleft,xright, N,retstep=True)

T = 20
dt = 0.0005
n = int(T / dt)

#----------SOLUTION FOR VELOCITY USING STEADY STATE METHOD---------------------------------------

v = np.linspace(1,1/7,N)


for i in range(n):
    v[0] = 1
    v[-1] = v[-2]

    deriv_velo = -(v[1:-1]-v[0:-2])/dx
    nonderiv_velo = -((7*v[1:-1]-1)*(1-v[1:-1])/v[1:-1])
    v[1:-1]=v[1:-1]+dt*(deriv_velo+nonderiv_velo)


#--------SOLUTION FOR NUMBER DENSITY USING STEADY STATE METHOD ----------------------------------

#neccesary parameters:
E_0 = 1 * 0.0000016021773 #erg: gcm^2/s^2
m_H = 1.6*10**(-24) #g
c = 3e11 #cm

n_0 = 1*10**(20)   #1/cm^3
v_0 = (2*E_0/m_H)**(0.5)  #cm/s
k = 1.38*10**(-16)
eps_c = 8.01*10**(-9)

T_e = 10**6
b = 20.3
n_eq = b*T_e**3/n_0
n_s = 2.03*10**(19)/n_0
sigma_c = 6.65*10**(-25)



y = np.linspace(n_s, 1000*n_eq, N)
v = np.linspace(1,1/7,N)

theta_0 = T_e*k

for j in range(n):
    y[0] = n_s
    y[-1] = y[-2]

    v[0] = 1
    v[-1] = v[-2]

    U = -(v[1:-1]-v[0:-2])/dx
    V = -((7*v[1:-1]-1)*(1-v[1:-1])/v[1:-1])


    theta_e = v_0**2*m_H*(1-v[1:-1])/y[1:-1] + theta_0

    D1 = c*v[1:-1]/(3*n_0*sigma_c)

    D2 = n_0*((v[1:-1]**2-8*v[1:-1]+1)*v_0/(6*v[1:-1]))
    D3 = n_0**2*v_0**2**sigma_c*(((7*v[1:-1]-1)*(1-v[1:-1])))/(v[1:-1]*2*c)

    T_e= theta_e/k

    lamb = eps_c/theta_e

    g = 1.226 - 0.475*np.log(lamb) + 0.0013*(np.log(lamb))**2
    E = -np.log(lamb) - 0.5772
    Q = 5.692*10**(-12)*((n_0)/(v[1:-1]))**2*g*E/(T_e**(0.5))

    n_eq = b*T_e**3
#    print(D1*((y[2:]-2.*y[1:-1]+y[0:-2])/(dx*dx)))
#    print(D2*(y[1:-1]-y[0:-2])/dx)
#    print(D3*(y[1:-1]-y[0:-2])/dx)
#    print(Q*(1-y[1:-1]/n_eq))
    A = ((y[2:]-2.*y[1:-1]+y[0:-2])/(dx*dx))
    B = (y[1:-1]-y[0:-2])/dx
    C = y[1:-1]
    R = (1-y[1:-1]/n_eq)

    y[1:-1] = y[1:-1] + dt*(A+B+C+R)
    v[1:-1]=v[1:-1]+dt*(U+V)

#-------------------PLOTTING------------------------------------------------

plt.plot(x, y, label='$n(x) $')
plt.plot(x, v, label='$v(x)$')
plt.plot(x[1:-1],T_e)

plt.yscale('log')
plt.grid()
plt.legend(framealpha=1)
plt.show()
