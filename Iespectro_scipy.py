import numpy as np
from scipy import integrate
from astropy import constants as const
from astropy import units as u
from math import pi



datos=np.loadtxt('sun_AM0.dat')
x=datos[:,0]*u.nm
y=datos[:,1]*u.W*((u.m)**-2)*((u.nm)**-1)
wavel=x.to('um')
flujo=y.to('erg / (s cm2 um)')
Itrap=integrate.trapz(flujo,wavel)
print Itrap  
