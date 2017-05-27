from math import *
from numpy import *
import matplotlib.pyplot as plt

beta = -2.92e8
g = 0.54e-10
xi = 1.56e9

p = arange(0.100, 0.750, 0.001)
w = sqrt(g) / (p * (xi * pi ** 2 + 0.5 * beta) ** 0.5)
w = w * 1e9

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.plot(p, w, 'b-', linewidth=2)
#plt.legend(loc='upper left')
#plt.xlabel(r'\textbf{$P_0$} ($C/m^2$)')
plt.xlabel(r'\textit{$P_0$} ($C/m^2$)')
plt.ylabel(r'\textit{$\xi_{180^\circ}$} ($nm$)',fontsize=16)

        
plt.savefig('figure1.eps')          
plt.show()
