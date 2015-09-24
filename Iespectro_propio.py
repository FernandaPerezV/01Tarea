#Este programa realiza la integral con el metodo de los trapecios la integral del 
#en longitud de onda.

import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from math import pi
from astropy import constants as const

#cargamos los datos
datos=np.loadtxt('sun_AM0.dat')
x=datos[:,0]*u.nm
y=datos[:,1]*u.W*((u.m)**-2)*((u.nm)**-1)


wavel=x.to('um')
flujo=y.to('erg / (s cm2 um)')
n=len(wavel)

integral=0
for i in range(n-1):
 dx=wavel[i+1]-wavel[i]
 sumy=flujo[i]+flujo[i+1]
 integral+=0.5*dx*sumy

print ("(cte solar) integral del espectro en longitud de onda=", integral)

