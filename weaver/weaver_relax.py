import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from matplotlib.animation import FuncAnimation
import sys
import matplotlib.animation as animation
#----------------GLOBAL PARAMETERS-----------------------------------------------------------------

N = 100

xleft = 0
xright = 8

x,dx= np.linspace(xleft,xright, N,retstep=True)

T = 20
dt = 0.0005
n = int(T / dt)


E_0 = 1 * 0.0000016021773 #erg: gcm^2/s^2
m_H = 1.6*10**(-24) #g
c = 2.99792458e10 #cm/s

n_0 = 1*10**(20)   #1/cm^3
v_0 = (2*E_0/m_H)**(0.5)  #cm/s
k = 1.38*10**(-16)
eps_c = 8.01*10**(-9)

Te = 10**6
b = 20.3
n_eq = b*Te**3/n_0
n_s = 2.03*10**(19)/n_0
sigma_c = 6.65*10**(-25)
m_e = 0.91*10**(-27) #g
r_0 = 2.8*10**(-13) #cm
alpha = 1/137
#----------SOLUTION FOR VELOCITY USING STEADY STATE METHOD---------------------------------------

v = np.linspace(1,1/7,N)
itcount = 0
itcount2 = 0

for i in range(n):
    v[0] = 1
    v[-1] = v[-2]

    deriv_velo = -(v[1:-1]-v[0:-2])/dx
    nonderiv_velo = -((7*v[1:-1]-1)*(1-v[1:-1])/v[1:-1])
    v[1:-1]=v[1:-1]+dt*(deriv_velo+nonderiv_velo)
    if i < 1000:
        if itcount2 == 50:
            print('plotted small')
            plt.plot(x,v,color = 'blue')
            itcount2 =0
    if itcount == 1000:
        plt.plot(x,v,color = 'blue')
        itcount = 0
    itcount2 += 1
    itcount += 1



plt.show()
#----------SOLUTION FOR VELOCITY ANALYTICALLY -------------------------------------------------

nhalf = int(N/2)
ymin = 1./7. + 0.000001         # avoid log(0) below
ymax = 1. - 0.000001            # avoid log(0) below
yanalytic = np.linspace(ymin, ymax, 1000)
xanalytic = x[nhalf] + (1./42.)*np.log((1 - yanalytic)**7/(7.*yanalytic - 1.)) \
                     - (1./42.)*np.log((1 - v[nhalf])**7/(7.*v[nhalf] - 1.))

#plt.plot(xanalytic*10**(-7),vanalytic-1/7*v_0)
#plt.show()



#-----------MAKE GIF-------------------------------------------------------------

#----------SOLUTION FOR VELOCITY USING ODEINT------------------------------------------------

def velocity(v,x):
    dvdx = -((8*v-7*v**2-1**2)/(v))
    return dvdx

sol = odeint(velocity,0.99999,x)
sol_velo = np.reshape(sol,N)

#--------SOLUTION FOR NUMBER DENSITY USING STEADY STATE METHOD ----------------------------------

#neccesary parameters:


print(sigma_c*n_0*v_0/(2*c))

y = np.linspace(n_s, 10000*n_eq, N)
v = np.linspace(1,1/7,N)

theta_0 = Te*k




for j in range(n):
    y[0] = n_s
    #y[-2] = n_s*10**4
    y[-1] = y[-2]

    v[0] = 1
    v[-1] = v[-2]

    U = -(v[1:-1]-v[0:-2])/dx
    V = -((7*v[1:-1]-1)*(1-v[1:-1])/v[1:-1])


    theta_e = (1-v[1:-1])/y[1:-1] + theta_0


    lamb = eps_c/theta_e


    tau = (theta_e*1.4*10**(-6))/(m_e*c**2)
    WNR = 16*((1)/(v[1:-1]))**2*m_e*c*(2*(theta_e*10**(-6))/(np.pi*m_e*c**2))**(0.5)*r_0**2*alpha
    WEI = np.zeros(len(WNR))
    WEE = np.zeros(len(WNR))
    QRCB = np.zeros(len(WNR))
    for u  in range(len(WNR)):
        if tau[u] <= 1.5:
            WEI[u] = WNR[u]*(19/24)*tau[u]*(1 + 1.022*tau[u]+0.221*tau[u]**2 - 0.239*tau[u]**3)
            WEE[u] = WNR[u]*3*tau[u]*(1-0.128*tau[u] + 0.898*tau[u]**2 - 0.439*tau[u]**3)
        else:
            WEI[u] = WNR[u] * (9/8)*(np.pi/2)**(0.5)*tau[u]**(0.5)*(np.log(2*tau[u]) + 0.923) - 1
            WEE[u] = 9/4*(np.pi/2)**(0.5)*tau[u]**(0.5)*(np.log(2*tau[u]) + 0.673)

        if lamb[u] <= 1/5:
            QRCB[u] = (WEI[u]+WEE[u])/(theta_e[u]*10**(-6))

        else:
            QRCB[u] = (WEI[u]+WEE[u])/(5*lamb[u]*theta_e[u]*10**(-6))

    #print(theta_e)
    D1 =v[1:-1]

    D2 = ((13*v[1:-1]**2-8*v[1:-1]+1)/(6*v[1:-1])) #changed after writeup
    D3 = ((7*v[1:-1]-1)*(1-v[1:-1]))/(v[1:-1]*2)

    T_e= theta_e/k
    n_eq = b*T_e**3






    g = 1.226 - 0.475*np.log(lamb) + 0.0013*(np.log(lamb))**2
    E = -np.log(lamb) - 0.5772
    QB = 5.692*10**(-12)*((1)/(v[1:-1]))**2*g*E/(T_e**(0.5))
    Q = 1*10**(-2)*QB + 10**5*QRCB

#    print(D1*((y[2:]-2.*y[1:-1]+y[0:-2])/(dx*dx)))
#    print(D2*(y[1:-1]-y[0:-2])/dx)
#    print(D3*(y[1:-1]-y[0:-2])/dx)
#    print(Q*(1-y[1:-1]/n_eq))
    #print(D2)
    A = D1*((y[2:]-2.*y[1:-1]+y[0:-2])/(dx*dx))#*c*sol_velo[1:-1]/(3*n_0*sigma_c)
    B = D2*(y[1:-1]-y[0:-2])/dx#*n_0*((sol_velo[1:-1]**2-8*sol_velo[1:-1]+1)*v_0/(6*sol_velo[1:-1]))
    C = D3*y[1:-1]#*n_0**2*v_0**2**sigma_c*(((7*sol_velo[1:-1]-1)*(1-sol_velo[1:-1])))/(sol_velo[1:-1]*2*c)
    R = Q#*5.692*10**(-12)*((n_0)/(sol_velo[1:-1]))**2*g*E/(T_e**(0.5))

    y[1:-1] = y[1:-1] + dt*(1.5*A+1*10**(-6)*B+sigma_c*1*10**(-12)*C+R)
    v[1:-1]=v[1:-1]+dt*(U+V)

print(y)
print('this is QRCB' + str(QRCB))
print('this is QB' + str(QB))
#print(v)

temp_e = np.zeros(len(y))
#neq = np.zeros(len(sol_num.y[0]))
i = 0
for i in range(len(y)):
    temp_e[i] = v_0**2*m_H*(1-v[i])/(y[i]*k) +Te
#-------------------PLOTTING------------------------------------------------



fig, ax1 = plt.subplots()

ax1.plot(x, (v)*v_0, label='$v(x)$', color = 'blue')
ax1.plot(x,temp_e, label ='T(x)', color = 'red')
ax1.set_yscale('log')
ax1.set_ylim(bottom = 1e5)
ax1.set_ylabel('Temperature [K]; Velocity [cm/s]')
ax1.tick_params(axis='y')
ax1.set_xlabel('Distance in Shockframe [$10^{6}$cm] ')


ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

#ax2.set_ylabel('sin', color=color)  # we already handled the x-label with ax1
ax2.plot(x, y*n_0, label='$n_{\gamma}(x) $',color = 'green')
ax2.set_yscale('log')
ax2.set_ylim(bottom = 1e19)
ax2.set_ylabel('Number Density [$cm^{-3}$]')
ax2.tick_params(axis='y')
ax2.grid()


lines, labels = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax2.legend(lines + lines2, labels + labels2, loc = 'lower right')
#fig.tight_layout()  # otherwise the right y-label is slightly clipped
ax1.grid()
plt.title('Shock Structure for $\epsilon_{0} = 1 MeV$')
plt.show()
