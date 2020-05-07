import numpy as np
from numpy import exp
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from uncertainties import ufloat
from scipy.stats import sem

x_inv = np.genfromtxt('./Data/three_t_inverse.txt', unpack=True) #import data
x_fullLU = np.genfromtxt('./Data/three_t_full_LU.txt', unpack=True) #import data
x_parLU = np.genfromtxt('./Data/three_t_partial_LU.txt', unpack=True) #import data

diff1, diff2, diff3 = np.genfromtxt('./Data/three_relative_error.txt', unpack=True) #import data

l=len(x_inv)
x0 = np.linspace(0, l, l)

def f(x, a,b):
    return a*x+b

def f2(x, a,b,c):
    return a+b*x**c

def g(x, a,b):
    return np.exp(b)*x**a

lx0=np.log(2**x0)


lxinv=np.log(x_inv)
lxfull=np.log(x_fullLU)
lxpar=np.log(x_parLU)

x = np.linspace(2**3, 2**l, 2**l)

params1, cov1 = curve_fit(f2, lx0[2:], lxinv[2:])
covv1 = np.sqrt(np.diag(cov1))
print('a1 ist ',params1[0],'pm',covv1[0])
print('b1 ist ',params1[1],'pm',covv1[1])

params2, cov2 = curve_fit(f, lx0[2:], lxfull[2:])
covv2 = np.sqrt(np.diag(cov2))
print('a1 ist ',params2[0],'pm',covv2[0])
print('b1 ist ',params2[1],'pm',covv2[1])


params3, cov3 = curve_fit(f, lx0[2:], lxpar[2:])
covv3 = np.sqrt(np.diag(cov3))
print('a1 ist ',params3[0],'pm',covv3[0])
print('b1 ist ',params3[1],'pm',covv3[1])

#Plot Verfahren aus a)
plt.figure(1)
plt.plot (x, f2(x, *params1), 'k-')
plt.plot (x, g(x, *params2), 'y-')
plt.plot (x, g(x, *params3), 'r-')
plt.plot (2**x0, x_inv, 'kx', label='Invertierung')
plt.plot (2**x0, x_fullLU, 'yx', label='Full LU')
plt.plot (2**x0, x_parLU, 'rx', label='Partial LU')
plt.grid()
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'Dimension N')
plt.ylabel(r'Laufzeit in s')
plt.legend(loc="best")
plt.tight_layout()
plt.savefig('Plots/3.pdf')
plt.show()

ldiff=len(diff1)
xdiff = np.linspace(0, ldiff, ldiff)

plt.figure(2)
plt.plot (2**xdiff, diff1*100, 'kx', label='Invertierung')
plt.plot (2**xdiff, diff2*100, 'yx', label='LU (Vollst. Piv.)')
plt.plot (2**xdiff, diff3*100, 'rx', label='LU (Teilw. Piv.)')
plt.grid()
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'Dimension N')
plt.ylabel(r'Relative Abweichung in %')
plt.legend(loc="best")
plt.tight_layout()
plt.savefig('Plots/3_diff.pdf')
plt.show()
