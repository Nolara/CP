import numpy as np
from numpy import exp
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from uncertainties import ufloat
from scipy.stats import sem

import numpy as np
from numpy import exp
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from uncertainties import ufloat
from scipy.stats import sem

I11, h11, a11 = np.genfromtxt('./Data/3_1_1.txt', unpack=True) #import data
I12, h12, a12 = np.genfromtxt('./Data/3_1_2.txt', unpack=True) #import data
I13, h13, a13 = np.genfromtxt('./Data/3_1_3.txt', unpack=True) #import data


#Plot für Integral 1
plt.figure(1)
plt.title(r"Integral 1")
plt.plot (h11, a11, 'x', label='Trapezregel', color='darkblue')
plt.plot (h12, a12, 'o', label='Mittelpunktsregel', color='darkorange')
plt.plot (h13, a13, 'x', label='Simpsonsregel', color='forestgreen')

plt.grid()
plt.xlabel(r'Schrittweite $h$')
plt.ylabel(r'$rel. Abweichung$')
plt.xscale('log')
plt.yscale('log')
plt.legend(loc="best")
plt.tight_layout()
plt.savefig('Plots/A3_Integral1.pdf')


I21, h21, a21 = np.genfromtxt('./Data/3_2_1.txt', unpack=True) #import data
I22, h22, a22 = np.genfromtxt('./Data/3_2_2.txt', unpack=True) #import data
I23, h23, a23 = np.genfromtxt('./Data/3_2_3.txt', unpack=True) #import data

#Plot für Integral 2
plt.figure(2)
plt.title(r"Integral 2")
plt.plot (h21, a21, 'x', label='Trapezregel', color='darkblue')
plt.plot (h22, a22, 'o', label='Mittelpunktsregel', color='darkorange')
plt.plot (h23, a23, 'x', label='Simpsonsregel', color='forestgreen')

plt.grid()
plt.xlabel(r'Schrittweite $h$')
plt.ylabel(r'$rel. Abweichung$')
plt.xscale('log')
plt.yscale('log')
plt.legend(loc="best")
plt.tight_layout()
plt.savefig('Plots/A3_Integral2.pdf')
