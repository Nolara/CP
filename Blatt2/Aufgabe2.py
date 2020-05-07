import numpy as np
from numpy import exp
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from uncertainties import ufloat
from scipy.stats import sem

x1, x2 ,y,z, = np.genfromtxt('./Data/two_data.txt', unpack=True) #import data
l=len(x1)
x0 = np.linspace(0, l, l)

def f(x, a,b):
    return a*x+b

def g(x,a,b):
    return np.exp(b)*np.exp(a*x)

lx0=np.log(2**x0)
x1log=np.log(x1)
x2log=np.log(x2)
ylog=np.log(y)
zlog=np.log(z)

params1, cov1 = curve_fit(f, lx0[2:], x1log[2:])
covv1 = np.sqrt(np.diag(cov1))
print('a1 ist ',params1[0],'pm',covv1[0])
print('b1 ist ',params1[1],'pm',covv1[1])


#Plot Verfahren aus a)
plt.figure(1)
plt.plot (2**x0, x1, 'k.', label='Zufallsmatrix')
#plt.plot (2**x0, g(2**x0, *params1), 'r-', label='Lineare Ausgleichsgerade')
plt.plot (2**x0, x2, 'y.', label='Zufallsvektor')
plt.plot (2**x0, y, 'r.', label='LU Zerlegung')
plt.plot (2**x0, z, 'b.', label='Lösen')
plt.grid()
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'Dimension N')
plt.ylabel(r'Laufzeit in ns')
plt.legend(loc="best")
plt.tight_layout()
plt.savefig('Plots/2.pdf')
plt.show()
