import numpy as np
from numpy import exp
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from uncertainties import ufloat
from scipy.stats import sem


x1,y1 = np.genfromtxt('./Data/two_1_in.txt', unpack=True) #import Dichtematrix
x2,y2 = np.genfromtxt('./Data/two_1_out.txt', unpack=True)
x3,y3 = np.genfromtxt('./Data/two_2_in.txt', unpack=True)
x4,y4 = np.genfromtxt('./Data/two_2_out.txt', unpack=True)


plt.figure(1)
plt.plot (x1,y1, 'rx', label='|x|<a')
plt.plot (x2,y2, 'bx', label='|x|>a')
plt.plot (x2,8/(x2**2), 'bx', label='Theorie')
plt.grid()
plt.xlabel('x/a')
plt.ylabel('$\Phi$')
plt.legend(loc="best")
plt.tight_layout()
plt.savefig('Plots/2_a.pdf')

plt.figure(2)
plt.plot (x3,y3, 'rx', label='|x|<a')
plt.plot (x4,y4, 'bx', label='|x|>a')
plt.plot (x4,1/(3*x4**2), 'bx', label='Theorie')
plt.grid()
plt.xlabel('x/a')
plt.ylabel('$\Phi$')
plt.legend(loc="best")
plt.tight_layout()
plt.savefig('Plots/2_b.pdf')
