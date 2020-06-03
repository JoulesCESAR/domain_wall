# Script for calculation of ferroelectric wall domain profile

from math import *
from numpy import *
import matplotlib.pyplot as plt

axes = plt.gca()
axes.set_xlim([-0.55,0.55])
axes.set_ylim([-0.9,0.9])


beta = -2.92e8
g = 0.54e-10
xi = 1.56e9
alpha0 = 7.6e5
q11 = 0.089
q12 = -0.026
s11 = -2.5
s12 = 9.0
e = 4.0/((8.8541878176e-12)*(pi**3))
n = 3.0
T = 50
t = 0.2e-9
D = 100e-9

sigma0 = 0.0
sigma1 = -6e9
sigma2 = -7e9
sigma3 = 0.7e9
sigma4 = 3e9

A0 = alpha0*(T-752.0) + 2.0*((e*t)/(n*D))*(1-exp(-2.0*pi*n))-2*(q11+2*q12)*sigma0
A1 = alpha0*(T-752.0) + 2.0*((e*t)/(n*D))*(1-exp(-2.0*pi*n))-2*(q11+2*q12)*sigma1
A2 = alpha0*(T-752.0) + 2.0*((e*t)/(n*D))*(1-exp(-2.0*pi*n))-2*(q11+2*q12)*sigma2
A3 = alpha0*(T-752.0) + 2.0*((e*t)/(n*D))*(1-exp(-2.0*pi*n))-2*(q11+2*q12)*sigma3
A4 = alpha0*(T-752.0) + 2.0*((e*t)/(n*D))*(1-exp(-2.0*pi*n))-2*(q11+2*q12)*sigma4

p0 = sqrt((-beta + sqrt(beta**2-4.0*A0*xi))/(2*xi))
p1 = sqrt((-beta + sqrt(beta**2-4.0*A1*xi))/(2*xi))
p2 = sqrt((-beta + sqrt(beta**2-4.0*A2*xi))/(2*xi))
p3 = sqrt((-beta + sqrt(beta**2-4.0*A3*xi))/(2*xi))
p4 = sqrt((-beta + sqrt(beta**2-4.0*A4*xi))/(2*xi))

w0 = sqrt(g) / (p0 * (xi * pi ** 2 + 0.5 * beta) ** 0.5)
C0 = (6.0 * xi * p0**2 + 3.0 * beta)/(4.0 * xi * p0**2 + 3.0 * beta)

w1 = sqrt(g) / (p1 * (xi * pi ** 2 + 0.5 * beta) ** 0.5)
C1 = (6.0 * xi * p1**2 + 3.0 * beta)/(4.0 * xi * p1**2 + 3.0 * beta)

w2 = sqrt(g) / (p2 * (xi * pi ** 2 + 0.5 * beta) ** 0.5)
C2 = (6.0 * xi * p2**2 + 3.0 * beta)/(4.0 * xi * p2**2 + 3.0 * beta)

w3 = sqrt(g) / (p3 * (xi * pi ** 2 + 0.5 * beta) ** 0.5)
C3 = (6.0 * xi * p3**2 + 3.0 * beta)/(4.0 * xi * p3**2 + 3.0 * beta)

w4 = sqrt(g) / (p4 * (xi * pi ** 2 + 0.5 * beta) ** 0.5)
C4 = (6.0 * xi * p4**2 + 3.0 * beta)/(4.0 * xi * p4**2 + 3.0 * beta)

x = arange(-7.0*w0*1e9, 7.0*w0*1e9, 0.01)
plt.xticks(fontsize = 10)
plt.yticks(arange(-0.9,0.9,0.1), fontsize = 9)

P0 = (p0 * sinh((x*1e-9)/w0))/sqrt(C0 + sinh((x*1e-9)/w0)**2)
P1 = (p1 * sinh((x*1e-9)/w1))/sqrt(C1 + sinh((x*1e-9)/w1)**2)
P2 = (p2 * sinh((x*1e-9)/w2))/sqrt(C2 + sinh((x*1e-9)/w2)**2)
P3 = (p3 * sinh((x*1e-9)/w3))/sqrt(C3 + sinh((x*1e-9)/w3)**2)
P4 = (p4 * sinh((x*1e-9)/w4))/sqrt(C4 + sinh((x*1e-9)/w4)**2)

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.plot(x, P0, '-',color = 'black',linewidth=1.2, label=r'\textit{$\sigma = 0$}')
plt.plot(x, P1, '--',color ='blue',linewidth=1.2, label=r'\textit{$\sigma = -6 GPa$}')
plt.plot(x, P2, '-.',color = 'blue', linewidth=1.2, label=r'\textit{$\sigma = -7 GPa$}')
plt.plot(x, P3, ':',color ='green',linewidth=1.2, label=r'\textit{$\sigma = 0.7 GPa$}')
plt.plot(x, P4, '+',color = 'green', linewidth=1.2, label=r'\textit{$\sigma = 3 GPa$}')

plt.xlabel(r'$x$ ($nm$)')
plt.ylabel(r'$P$ ($C/m^2$)',fontsize=16)
plt.legend(loc=2)

plt.savefig('/home/julio/Documents/Development/domain_wall/Figures/figure6.jpg', figsize=(3.30, 3.30), dpi=100)
plt.show()
