import numpy as np
poly_index = 1.5
mu = 0.19
R = 800*6.9*10**(8)
Mass = 15*1.9*10**(30)
kappa = 0.034

def breakout_point(R,Mass,poly_index,kappa,mu):
    c = 3*10**8
    RHO = Mass/R**3
    r = R*(1-np.power((kappa*R*RHO)/(c*(poly_index+1)),(1)/(poly_index*(mu-1)-1))/1800000)
    return r


r = breakout_point(R,Mass,poly_index,kappa,mu)
print("We did it boyz d =",(R-r)/R,"R")
