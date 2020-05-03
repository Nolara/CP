import numpy as np
from numpy import exp
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from uncertainties import ufloat
from scipy.stats import sem

x1, x2 ,y,z, = np.genfromtxt('./Data/two_data.txt', unpack=True) #import data

l = np.linspace(0, 999, 999)

#Plot Verfahren aus a)
plt.figure(1)
plt.plot (l, x1, 'k.', label='Zufallsmatrix')
plt.plot (l, x2, 'y.', label='Zufallsvektor')
plt.plot (l, y, 'r.', label='LU Zerlegung')
plt.plot (l, z, 'b.', label='Lösen')
plt.grid()
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'Dimension N')
plt.ylabel(r'Laufzeit in ns')
plt.legend(loc="best")
plt.tight_layout()
plt.savefig('Plots/2.pdf')
plt.show()
