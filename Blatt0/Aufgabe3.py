import numpy as np
from numpy import exp
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from uncertainties import ufloat
from scipy.stats import sem

dt, en, es, an = np.genfromtxt('./Data/three.txt', unpack=True) #import data

en_rel=abs(an-en)/an*100 #calculate the relative error
es_rel=abs(an-es)/an*100

plt.figure(1)
plt.plot (dt, en, 'k-', label='Normales Eulerverfahren')
plt.plot (dt, es, 'r-', label='Symmetrisches Eulerverfahren')
plt.plot (dt, an, 'g-', label='Analytisches Ergebnis')
plt.yscale('log')
plt.grid()
plt.xlabel(r't')
plt.legend(loc="best")
plt.tight_layout()
plt.savefig('Plots/3.pdf')

plt.figure(2)
plt.plot (dt, en_rel, 'k-', label='Relativer Fehler des normalen Eulerverfahren')
plt.plot (dt, es_rel, 'r-', label='Relativer Fehler des symmetrischen Eulerverfahren')
plt.yscale('log')
plt.xlim(0,10)
plt.grid()
plt.xlabel(r't')
plt.ylabel(r'Fehler in %')
plt.legend(loc="best")
plt.tight_layout()
plt.savefig('Plots/3rel.pdf')
