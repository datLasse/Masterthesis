from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp

E_0 = 1 * 0.0000016021773 #erg: gcm^2/s^2
m_H = 1.6*10**(-24) #g
c = 3e11 #cm
sigma_c = 2*10**(-24)
n_0 = 1*10**(20)   #1/cm^3
v_0 = (2*E_0/m_H)**(0.5)  #cm/s


v = np.linspace(1/6*v_0,999/1000*v_0,100)

def place(v):
    x = (c/(21*sigma_c*n_0*v_0))*(np.log(np.power(v_0-v,7)/((7*v-v_0)*v_0**6))-np.log((1/3)*np.power(3/7,7)))
    return x


def velocity(v,x):
    dvdx = -sigma_c*n_0*v_0*((8*v_0*v-7*v**2-v_0**2)/(2*v*c))
    return dvdx

def vat0(a):
    return odeint(velocity,a,x)

def parameter(v_0):
    itcount = 0
    terminator = 1
    b = v_0
    while terminator:
        sol = vat0(b)
        if int(sol[3]) in range(int(v_0*4/7-10),int(v_0*4/7+10)) :
            terminator = 0
            return b
        else:
            itcount += 1
            b = 1/itcount*v_0
            if itcount>10:
                b = itcount*v_0


def plotv(v,v_0):
    vv = v - (v_0/7)
    return vv

x = np.zeros(len(v))
vv = np.zeros(len(v))
vn = np.zeros(len(v))
i = 0
for i in range(len(v)):
    x[i] = place(v[i])
    vv[i] = v[i] - v_0/7
print(x)
x_num = np.linspace(x[99],x[0],1000)
sol = odeint(velocity,v_0,x_num,printmessag = )

print(x[99])
print(x[0])
print(sol)

plt.plot(x,vv,label='$v(x)$')
plt.plot(x_num,sol-v_0/7,label= '$v(x) numerical$')
plt.yscale('log')
plt.grid(alpha=0.5)
plt.legend(framealpha=1)
plt.show()
