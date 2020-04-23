import numpy as np
from numpy import exp
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from uncertainties import ufloat
from scipy.stats import sem

dt_a, en_a, es_a, an_a = np.genfromtxt('./Data/three_a.txt', unpack=True) #import data
dt_b, en_b, es_b, an_b = np.genfromtxt('./Data/three_b.txt', unpack=True) #import data

en_rel_a=abs(an_a-en_a)/an_a*100 #calculate the relative error
es_rel_a=abs(an_a-es_a)/an_a*100
en_rel_b=abs(an_b-en_b)/an_b*100 #calculate the relative error
es_rel_b=abs(an_b-es_b)/an_b*100

plt.figure(1)
plt.plot (dt_a, en_a, 'k-', label='Normales Eulerverfahren')
plt.plot (dt_a, es_a, 'r-', label='Symmetrisches Eulerverfahren')
plt.plot (dt_a, an_a, 'g-', label='Analytisches Ergebnis')
plt.yscale('log')
plt.grid()
plt.xlabel(r't')
plt.legend(loc="best")
plt.tight_layout()
plt.savefig('Plots/3_a.pdf')

plt.figure(2)
plt.plot (dt_a, en_rel_a, 'k-', label='Relativer Fehler des normalen Eulerverfahren')
plt.plot (dt_a, es_rel_a, 'r-', label='Relativer Fehler des symmetrischen Eulerverfahren')
plt.yscale('log')
plt.xlim(0,10)
plt.grid()
plt.xlabel(r't')
plt.ylabel(r'Fehler in %')
plt.legend(loc="best")
plt.tight_layout()
plt.savefig('Plots/3rel_a.pdf')

plt.figure(3)
plt.plot (dt_b, en_b, 'k-', label='Normales Eulerverfahren')
plt.plot (dt_b, es_b, 'r-', label='Symmetrisches Eulerverfahren')
plt.plot (dt_b, an_b, 'g-', label='Analytisches Ergebnis')
plt.yscale('log')
plt.grid()
plt.xlabel(r't')
plt.legend(loc="best")
plt.tight_layout()
plt.savefig('Plots/3_b.pdf')

plt.figure(4)
plt.plot (dt_b, en_rel_b, 'k-', label='Relativer Fehler des normalen Eulerverfahren')
plt.plot (dt_b, es_rel_b, 'r-', label='Relativer Fehler des symmetrischen Eulerverfahren')
plt.yscale('log')
plt.xlim(0,10)
plt.grid()
plt.xlabel(r't')
plt.ylabel(r'Fehler in %')
plt.legend(loc="best")
plt.tight_layout()
plt.savefig('Plots/3rel_b.pdf')


plt.figure(5)
plt.plot (dt_a, en_rel_a, 'k-', label='Relativer Fehler des normalen Eulerverfahren aus Aufgabenteil a)')
plt.plot (dt_a, es_rel_a, 'r-', label='Relativer Fehler des symmetrischen Eulerverfahren aus Aufgabenteil a)')
plt.plot (dt_b, en_rel_b, 'g-', label='Relativer Fehler des normalen Eulerverfahren aus Aufgabenteil b)')
plt.plot (dt_b, es_rel_b, 'c-', label='Relativer Fehler des symmetrischen Eulerverfahren aus Aufgabenteil b)')
plt.yscale('log')
plt.xlim(0,10)
plt.grid()
plt.xlabel(r't')
plt.ylabel(r'Fehler in %')
plt.legend(loc="best")
plt.tight_layout()
plt.savefig('Plots/3rel.pdf')
