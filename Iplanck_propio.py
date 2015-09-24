#Este programa calcula la integral de planck para un cuerpo negro con T=5778K
#mediante metodo de simpson
import numpy as np
from math import pi
from astropy import constants as const
from astropy import units as u

T=5778*u.K

def f(y):
 return ((np.tan(y))**3)* (1+(np.tan(y))**2) / ( np.exp(np.tan(y))-1)


n=100 #aquí se escogio n=100
x0=0.0001
xf=pi/2

dx = (xf-x0)/n
x = x0
suma = 0

for i in range(n/2):
 suma += f(x) + 4.*f(x+dx) + f(x+2*dx)
 x += 2*dx

integral=(dx/3.)*suma
integralPlanck=( (2*pi*const.h/const.c**2) * ((const.k_B*(5778*u.K)/const.h)**4) ) * integral
print('energia por unidad de tiempo y area de un cuerpo negro con Teff_sol=5778K = ', integralPlanck.to('erg / (s cm2)'))




