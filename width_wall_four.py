from math import *
from numpy import *
import matplotlib.pyplot as plt

axes = plt.gca()
axes.set_xlim([0.0,850.0])
axes.set_ylim([0.0,2.0])

beta = -2.92e8
g = 0.54e-10
xi = 1.56e9

alpha0 = 7.6e5
q11 = 0.089
q12 = -0.026
e = 4.0/((8.8541878176e-12)*(pi**3))
n = 3.0
D = 100e-9
t = 0.2e-9

sigma0 = 0.0
sigma1 = -2e9
sigma2 = -4e9
sigma3 = 0.5e9
sigma4 = 1e9

T = arange(0.0, 850.0, 1.5)

A0 = alpha0*(T-752.0) + 2.0*((e*t)/(n*D))*(1-exp(-2.0*pi*n))-2*(q11+2*q12)*sigma0
p0 = sqrt((-beta + sqrt(beta**2-4.0*A0*xi))/(2*xi))
w0 = sqrt(g) / (p0 * (xi * (p0 ** 2) + 0.5 * beta) ** 0.5)
w0 = w0 * 1e9

A1 = alpha0*(T-752.0) + 2.0*((e*t)/(n*D))*(1-exp(-2.0*pi*n))-2*(q11+2*q12)*sigma1
p1 = sqrt((-beta + sqrt(beta**2-4.0*A1*xi))/(2*xi))
w1 = sqrt(g) / (p1 * (xi * (p1 ** 2) + 0.5 * beta) ** 0.5)
w1 = w1 * 1e9

A2 = alpha0*(T-752.0) + 2.0*((e*t)/(n*D))*(1-exp(-2.0*pi*n))-2*(q11+2*q12)*sigma2
p2 = sqrt((-beta + sqrt(beta**2-4.0*A2*xi))/(2*xi))
w2 = sqrt(g) / (p2 * (xi * (p2 ** 2) + 0.5 * beta) ** 0.5)
w2 = w2 * 1e9

A3 = alpha0*(T-752.0) + 2.0*((e*t)/(n*D))*(1-exp(-2.0*pi*n))-2*(q11+2*q12)*sigma3
p3 = sqrt((-beta + sqrt(beta**2-4.0*A3*xi))/(2*xi))
w3 = sqrt(g) / (p3 * (xi * (p3 ** 2) + 0.5 * beta) ** 0.5)
w3 = w3 * 1e9

A4 = alpha0*(T-752.0) + 2.0*((e*t)/(n*D))*(1-exp(-2.0*pi*n))-2*(q11+2*q12)*sigma4
p4 = sqrt((-beta + sqrt(beta**2-4.0*A4*xi))/(2*xi))
w4 = sqrt(g) / (p4 * (xi * (p4 ** 2) + 0.5 * beta) ** 0.5)
w4 = w4 * 1e9

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.xticks(arange(0.0,850.0,50.0), fontsize = 10)
plt.yticks(arange(0.0,2.0,0.2), fontsize = 9)

plt.plot(T, w0, '-',color = 'black',linewidth=1.2, label=r'\textit{$\sigma = 0$}')
plt.plot(T, w1, '--',color = 'blue',linewidth=1.2, label=r'\textit{$\sigma = -2 GPa$}')
plt.plot(T, w2, '-.',color = 'blue', linewidth=1.2, label=r'\textit{$\sigma = -4 GPa$}')
plt.plot(T, w3, ':',color ='green',linewidth=1.2, label=r'\textit{$\sigma = 0.5 GPa$}')
plt.plot(T, w4, '+',color ='green',linewidth=1.2, label=r'\textit{$\sigma = 1 GPa$}')

#plt.legend(loc='upper left')
#plt.xlabel(r'\textbf{$P_0$} ($C/m^2$)')
plt.xlabel(r'$T$ ($K$)')
plt.ylabel(r'$\xi_{180^\circ}$ ($nm$)',fontsize=16)
plt.legend(loc=2)

plt.savefig('figure9.eps', figsize=(3.30, 3.30), dpi=100)
plt.show()
