import numpy as np
from numpy import exp
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from uncertainties import ufloat
from scipy.stats import sem

xa, ya, za = np.genfromtxt('./Data/two_a.txt', unpack=True) #import data
xb, yb, zb = np.genfromtxt('./Data/two_b.txt', unpack=True) #import data
xc, yc, zc = np.genfromtxt('./Data/two_c.txt', unpack=True) #import data
delta = np.genfromtxt('./Data/two_delta.txt', unpack=True) #import data

# Definiere relative Abweichung
arel=abs(ya-za)/za*100
brel=abs(yb-zb)/zb*100
crel=abs(yc-zc)/abs(zc)*100

#Plot Funktionen aus Aufgabenteil a
plt.figure(1)
plt.plot (xa, ya, 'r-', label='$f_a(x)$') #Funktion
plt.plot (xa, za, 'k-', label='$\~f_a(x)$') #Funktion ohne Auslöschung
plt.xscale('log')
plt.yscale('log')
plt.grid()
plt.xlabel(r'x')
plt.ylim(10e-12, 2e-8)
plt.legend(loc="best")
plt.tight_layout()
plt.savefig('Plots/2a.pdf')

#Plot Funktionen aus Aufgabenteil b
plt.figure(2)
plt.plot (xb, yb, 'r-', label='$f_b(x)$')
plt.plot (xb, zb, 'k-', label='$\~f_b(x)$')
plt.grid()
plt.xscale('log')
plt.xlabel(r'x')
plt.legend(loc="best")
plt.tight_layout()
plt.savefig('Plots/2b.pdf')

#Plot Funktionen aus Aufgabenteil c
plt.figure(3)
plt.title('$\delta$ = '+str(delta))
plt.plot (xc, yc, 'r-', label='$f_c(x)$')
plt.plot (xc, zc, 'k-', label='$\~f_c(x)$')
plt.grid()
plt.xlabel(r'x')
plt.legend(loc="best")
plt.tight_layout()
plt.savefig('Plots/2c.pdf')

#Plot relative Abweichung aus Aufgabenteil a
plt.figure(4)
plt.plot (xa, arel, 'r-', label='Relative Abweichung')
plt.grid()
plt.xlabel(r'x')
plt.xscale('log')
plt.ylabel(r'Abweichung in %')
plt.legend(loc="best")
plt.tight_layout()
plt.savefig('Plots/2arel.pdf')

#Plot relative Abweichung aus Aufgabenteil b
plt.figure(5)
plt.plot (xb, brel, 'r-', label='Relative Abweichung')
plt.grid()
plt.xlabel(r'x')
plt.xscale('log')
plt.ylabel(r'Abweichung in %')
plt.legend(loc="best")
plt.tight_layout()
plt.savefig('Plots/2brel.pdf')

#Plot relative Abweichung aus Aufgabenteil c
plt.figure(6)
plt.plot (xc, crel, 'r-', label='Relative Abweichung')
plt.grid()
plt.xlabel(r'x')
plt.ylabel(r'Abweichung in %')
plt.legend(loc="best")
plt.tight_layout()
plt.savefig('Plots/2crel.pdf')
