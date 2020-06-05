import numpy as np
from numpy import exp
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from uncertainties import ufloat
from scipy.stats import sem

t1, x1,v1= np.genfromtxt('./Data/2_1.txt', unpack=True)



plt.figure(1)
plt.plot (t1,x1, 'b-', label='Ort')
plt.plot (t1,v1, 'r-', label='Geschwindigkeit')
plt.grid()
plt.xlabel('t')
plt.ylabel('')
plt.legend(loc="best")
plt.tight_layout()
plt.savefig('Plots/2_1.pdf')
