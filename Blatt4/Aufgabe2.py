import numpy as np
from numpy import exp
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from uncertainties import ufloat
from scipy.stats import sem

#import Data

x1,y1 = np.genfromtxt('./Data/two_1_in.txt', unpack=True)
x2,y2 = np.genfromtxt('./Data/two_1_out.txt', unpack=True)
x3,y3 = np.genfromtxt('./Data/two_2_in.txt', unpack=True)
x4,y4 = np.genfromtxt('./Data/two_2_out.txt', unpack=True)

l=np.linspace(2,8,1000)

#plot Aufgabeteil a)
plt.figure(1)
plt.plot (x1,y1, 'rx', label='|x|<a')
plt.plot (x2,y2, 'bx', label='|x|>a')
plt.plot (l,8/(l), '-', label='Monopolterm', color='orange')
plt.grid()
plt.xlabel('x/a')
plt.ylabel('$\Phi/(𝜌_0/4𝜋𝜖_0)$')
plt.legend(loc="best")
plt.tight_layout()
plt.savefig('Plots/2_a.pdf')

#plot Aufgabeteil b)
plt.figure(2)
plt.plot (x3,y3, 'rx', label='|x|<a')
plt.plot (x4,y4, 'bx', label='|x|>a')
plt.plot (l,8/(3*l**2), '-', label='Dipolterm', color='orange')
plt.grid()
plt.xlabel('x/a')
plt.ylabel('$\Phi/(𝜌_0/4𝜋𝜖_0)$')
plt.legend(loc="best")
plt.tight_layout()
plt.savefig('Plots/2_b.pdf')
