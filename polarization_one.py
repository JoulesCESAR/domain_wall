# Script for calculation of ferroelectric wall domain profile

from math import *

alpha_0 = 7.5e5
T = 600

alpha = alpha_0 * (T - 752)
beta = -2.92e8
epsilon_0 = 8.8541878176e-12
t = 1e-9
D = 10e-9
d = 5e-9

epsilon_1 = 4.0/(epsilon_0 * (pi**3))

r_c = 0.5e-9

for x in [ 0.1*r_c, 0.2*r_c, 0.3*r_c, 0.4*r_c, 0.5*r_c]:

   P = sqrt(-alpha/beta - 0.5 * ((epsilon_1*t*d)/(beta*D**2))*(1-exp(-2 * pi*(D/d)))) * tanh(x/(2*r_c))
   print P

print alpha_0, beta, epsilon_0, epsilon_1
