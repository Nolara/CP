import numpy as np
from numpy import exp
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from uncertainties import ufloat
from scipy.stats import sem

x,y = np.genfromtxt('./Data/two_data.txt', unpack=True) #import data
m,n = np.genfromtxt('./Data/two_x0.txt', unpack=True) #import data
def f(x, a,b):
    return a*x+b

l = np.linspace(-7.5, 14, 50000)

#Plot Verfahren aus a)
plt.figure(1)
plt.plot (x, y, 'kx', label='Datenpunkte')
plt.plot(l, f(l, m,n), 'b-', label='Lineare Ausgleichsgerade')
plt.grid()
plt.xlim(-7,14)
plt.xlabel(r'$x$')
plt.ylabel(r'$y$')
plt.legend(loc="best")
plt.tight_layout()
plt.savefig('Plots/2.pdf')
