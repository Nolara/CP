import numpy as np
from numpy import exp
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from uncertainties import ufloat
from scipy.stats import sem

xa, ya = np.genfromtxt('./Data/two_a.txt', unpack=True) #import data
xb, yb = np.genfromtxt('./Data/two_b.txt', unpack=True) #import data
xc, yc = np.genfromtxt('./Data/two_c.txt', unpack=True) #import data

def f(x):
    return 1/(x)**(1/2)-1/(x+1)**(1/2)

xl=np.linspace(1000000,10000000, 50000)

plt.figure(1)
plt.plot (xa, ya, 'k-', label='Funktion a mit Auslöschung')
#plt.plot (xl, f(xl), 'r-', label='Exakter Wert')
plt.yscale('log')
plt.grid()
plt.xlabel(r't')
#plt.ylabel(r'')
plt.legend(loc="best")
plt.tight_layout()
plt.savefig('Plots/2a.pdf')
#plt.show()

plt.figure(2)
plt.plot (xb, yb, 'r-', label='Funktion b mit Auslöschung')
plt.plot (xb, ybo, 'k-', label='Funktion b ohne Auslöschung')
#plt.plot (xl, f(xl), 'r-', label='Exakter Wert')
plt.grid()
plt.yscale('log')
plt.xlabel(r'x')
#plt.ylabel(r'')
plt.legend(loc="best")
plt.tight_layout()
plt.savefig('Plots/2b.pdf')
plt.show()

plt.figure(3)
plt.plot (xc, yc, 'r-', label='Funktion c mit Auslöschung')
#plt.plot (xc, yco, 'k-', label='Funktion b ohne Auslöschung')
#plt.plot (xl, f(xl), 'r-', label='Exakter Wert')
plt.grid()
plt.yscale('log')
plt.xlabel(r'x')
#plt.ylabel(r'')
plt.legend(loc="best")
plt.tight_layout()
plt.savefig('Plots/2c.pdf')
plt.show()
