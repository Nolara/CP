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

arel=abs(ya-za)/za*100
brel=abs(yb-zb)/zb*100
crel=abs(yc-zc)/zc*100

plt.figure(1)
plt.plot (xa, ya, 'r-', label='Funktion a mit Auslöschung')
plt.plot (xa, za, 'k-', label='Funktion a ohne Auslöschung')
plt.xscale('log')
plt.yscale('log')
plt.grid()
plt.xlabel(r'x')
plt.ylim(10e-12, 2e-8)
plt.legend(loc="best")
plt.tight_layout()
plt.savefig('Plots/2a.pdf')

plt.figure(2)
plt.plot (xb, yb, 'r-', label='Funktion b mit Auslöschung')
plt.plot (xb, zb, 'k-', label='Funktion b ohne Auslöschung')
plt.grid()
plt.xscale('log')
plt.xlabel(r'x')
plt.legend(loc="best")
plt.tight_layout()
plt.savefig('Plots/2b.pdf')

plt.figure(3)
plt.title('$\delta$ = '+str(delta))
plt.plot (xc, yc, 'r-', label='Funktion c mit Auslöschung')
plt.plot (xc, zc, 'k-', label='Funktion c ohne Auslöschung')
plt.grid()
plt.xlabel(r'x')
plt.legend(loc="best")
plt.tight_layout()
plt.savefig('Plots/2c.pdf')

plt.figure(4)
plt.plot (xa, arel, 'r-', label='Relativer Fehler')
plt.grid()
plt.xlabel(r'x')
plt.xscale('log')
plt.ylabel(r'Fehler in %')
plt.legend(loc="best")
plt.tight_layout()
plt.savefig('Plots/2arel.pdf')

plt.figure(5)
plt.plot (xb, brel, 'r-', label='Relativer Fehler')
plt.grid()
plt.xlabel(r'x')
plt.xscale('log')
plt.ylabel(r'Fehler in %')
plt.legend(loc="best")
plt.tight_layout()
plt.savefig('Plots/2brel.pdf')

plt.figure(6)
plt.plot (xc, crel, 'r-', label='Relativer Fehler')
plt.grid()
plt.xlabel(r'x')
plt.ylabel(r'Fehler in %')
plt.legend(loc="best")
plt.tight_layout()
plt.savefig('Plots/2crel.pdf')
