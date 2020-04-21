import numpy as np
from numpy import exp
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from uncertainties import ufloat
from scipy.stats import sem

dt, en, es, an = np.genfromtxt('Data/three.txt', unpack=True)


plt.plot (dt, en, 'k-', label='Normales Eulerverfahren')
plt.plot (dt, es, 'r-', label='Symmetrisches Eulerverfahren')
plt.plot (dt, an, 'g-', label='Analytisches Ergebnis')
plt.grid()
plt.xlabel(r't')
#plt.ylabel(r'')
plt.legend(loc="best")
plt.tight_layout()
plt.savefig('3.pdf')
plt.show()


en_rel=abs(an-en)/an*100
es_rel=abs(an-es)/an*100
plt.plot (dt, en_rel, 'k-', label='Relativer Fehler des normalen Eulerverfahren')
plt.plot (dt, es_rel, 'r-', label='Relativer Fehler des symmetrischen Eulerverfahren')
plt.grid()
plt.xlabel(r't')
plt.ylabel(r'Abweichung in %')
plt.legend(loc="best")
plt.tight_layout()
plt.savefig('3rel.pdf')
plt.show()
