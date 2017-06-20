# Script for calculation of ferroelectric wall domain profile

from math import *
from numpy import *
import matplotlib.pyplot as plt


p = 0.75
beta = -2.92e8
g = 0.54e-10
xi = 1.56e9


w = sqrt(g) / (p * (xi * pi ** 2 + 0.5 * beta) ** 0.5)
C = (6.0 * xi * p**2 + 3.0 * beta)/(4.0 * xi * p**2 + 3.0 * beta)

x = arange(-7.0*w*1e9, 7.0*w*1e9, 0.001)
P = (p * sinh((x*1e-9)/w))/sqrt(C + sinh((x*1e-9)/w)**2)
  
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.plot(x, P, 'b-', linewidth=2)

plt.xlabel(r'\textit{$x$} ($nm$)')
plt.ylabel(r'\textit{$P$} ($C/m^2$)',fontsize=16)


plt.savefig('figure2.eps', figsize=(3.30, 3.30), dpi=100)
plt.show()
