import numpy as np
from numpy import exp
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from uncertainties import ufloat
from scipy.stats import sem

x_inv = np.genfromtxt('./Data/three_t_inverse.txt', unpack=True) #import data
x_fullLU = np.genfromtxt('./Data/three_t_full_LU.txt', unpack=True) #import data
x_parLU = np.genfromtxt('./Data/three_t_partial_LU.txt', unpack=True) #import data

diff1, diff2, diff3 = np.genfromtxt('./Data/three_difference.txt', unpack=True) #import data

l=len(x_inv)
x0 = np.linspace(0, l, l)

#Plot Verfahren aus a)
plt.figure(1)
plt.plot (x0, x_inv, 'k.', label='Invertierung')
plt.plot (x0, x_fullLU, 'y.', label='Full LU')
plt.plot (x0, x_parLU, 'r.', label='Partial LU')
plt.grid()
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'Dimension N')
plt.ylabel(r'Laufzeit in ns')
plt.legend(loc="best")
plt.tight_layout()
plt.savefig('Plots/3.pdf')
plt.show()

plt.figure(2)
plt.plot (x0, diff1, 'k.', label='Differenz 1')
plt.plot (x0, diff2, 'y.', label='Differenz 2')
plt.plot (x0, diff3, 'r.', label='Differenz 3')
plt.grid()
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'Dimension N')
plt.ylabel(r'Laufzeit in ns')
plt.legend(loc="best")
plt.tight_layout()
plt.savefig('Plots/3_diff.pdf')
plt.show()
