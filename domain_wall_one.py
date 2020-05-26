from math import *
from numpy import *
import matplotlib.pyplot as plt

fig = plt.figure()
ax = fig.add_subplot(111)

axes = plt.gca()
axes.set_xlim([0.0,0.75])
axes.set_ylim([0.0,12.0])

beta = -2.92e8
g = 0.54e-10
xi = 1.56e9

p = arange(0.0, 0.750, 0.01)
w = sqrt(g) / (p * (xi * (p ** 2) + 0.5 * beta) ** 0.5)
w = w * 1e9

plt.rc('text', usetex=True)
plt.rc('font', family='serif')


plt.xticks(fontsize = 10)
plt.yticks(arange(0.0,12.0,1.0), fontsize = 9)


ax.text(0.27, 6,r'$P_{crit} = \sqrt{-\frac{\beta}{2\zeta}}=0.3052 C/m^2$', rotation='vertical', color='red',fontsize=12)
plt.axvline(x=0.3052,linewidth=1,linestyle = 'dashed', color='red')
plt.plot(p, w, 'b-', linewidth=1)
#plt.legend(loc='upper left')
#plt.xlabel(r'\textbf{$P_0$} ($C/m^2$)')
plt.xlabel(r'$P_0$ ($C/m^2$)')
plt.ylabel(r'$\xi_{180^\circ}$ ($nm$)',fontsize=16)


plt.savefig('figure1.eps', figsize=(3.30, 3.30), dpi=100)
plt.show()
