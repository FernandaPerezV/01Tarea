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

def f(x):
 return  x**(3)/(np.exp(x)-1)

I=integrate.quad(f,0,np.inf)
Ictes=    ( (2*pi*const.h/const.c**2) * ((const.k_B*(5778*u.K)/const.h)**4) )* I[0]
Iquad=Ictes.to('erg / (s cm2)')
print Iquad 

