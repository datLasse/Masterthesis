import numpy as np


#Solar mass, radius and the speed of light are in SI - Units. They are changed to
#gramms and centimetrs in the lines 16 and 17 so that the result is in (g/cm^3)


mass_solar = 1.989e30
radius_solar = 696340000
c = 299792458
n = 1.5                                   #RSG n = 1.5 BSG & WR n = 3
kappa = 0.034                            #RSG & BSG k = 0.034 WR n = 0.02  in SI- Units
mu = 0.19
X = 500                                   # RSG X = 500 BSG X = 50 WR X = 5

constant = (((3*15*mass_solar*10**3)
/(4*np.pi*X**3*(radius_solar*10**2)**3))
*np.power((kappa*1800e3*3*15*mass_solar)
/(c*(n+1)*4*np.pi*X**2*radius_solar**2),
(1/(mu-((n+1)/n)))))

print(constant)
