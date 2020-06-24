import numpy as np
from numpy import exp
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from uncertainties import ufloat
from scipy.stats import sem
#Einlesen der Daten
t1, x1,v1, E1k, E1p= np.genfromtxt('./Data/2_1.txt', unpack=True)
t2, x2,v2, E2k, E2p= np.genfromtxt('./Data/2_2.txt', unpack=True)
t3, x3,v3, E3k, E3p= np.genfromtxt('./Data/2_3.txt', unpack=True)
t4, x4,v4, E4k, E4p= np.genfromtxt('./Data/2_4.txt', unpack=True)
t5, x5,v5, E5k, E5p= np.genfromtxt('./Data/2_5.txt', unpack=True)
tb, xb,vb, Ebk, Ebp= np.genfromtxt('./Data/2_b.txt', unpack=True)
tc, l_ab, f_ab, l_rk, f_rk= np.genfromtxt('./Data/2_Vergleich.txt', unpack=True)


l=np.linspace(0,10*np.pi, 1000)

#Plot für a=0.5
plt.figure(1)
plt.plot (t1,x1, 'b-', label='x')
plt.plot (t1,v1, 'r-', label='v')
plt.plot (l,np.exp(-0.25*l), 'k-', label='Einhüllende')
plt.plot (l,-np.exp(-0.25*l), 'k-')
plt.grid()
plt.xlabel('t')
plt.ylabel('')
plt.legend(loc="best")
plt.tight_layout()
plt.savefig('Plots/2_1.pdf')

#Plot für a=2
plt.figure(2)
plt.plot (t2,x2, 'b-', label='x')
plt.plot (t2,v2, 'r-', label='v')
plt.grid()
plt.xlabel('t')
plt.ylabel('')
plt.legend(loc="best")
plt.tight_layout()
plt.savefig('Plots/2_2.pdf')

#Plot für a=4
plt.figure(3)
plt.plot (l,np.exp(-(2-3**(0.5))*l), 'k-', label='Einhüllende')
plt.plot (l,-np.exp(-(2-3**(0.5))*l), 'k-')
plt.plot (t3,x3, 'b-', label='x')
plt.plot (t3,v3, 'r-', label='v')
plt.grid()
plt.xlabel('t')
plt.ylabel('')
plt.legend(loc="best")
plt.tight_layout()
plt.savefig('Plots/2_3.pdf')

#Plot für a=0
plt.figure(4)
plt.plot (t4,x4, 'b-', label='x')
plt.plot (t4,v4, 'r-', label='v')
plt.grid()
plt.xlabel('t')
plt.ylabel('')
plt.legend(loc="best")
plt.tight_layout()
plt.savefig('Plots/2_4.pdf')

#Plot für a=-0.1
plt.figure(5)
plt.plot (t5,x5, 'b-', label='x')
plt.plot (t5,v5, 'r-', label='v')
plt.plot (l,np.exp(0.05*l), 'k-', label='Einhüllende')
plt.plot (l,-np.exp(0.05*l), 'k-')
plt.grid()
plt.xlabel('t')
plt.ylabel('')
plt.legend(loc="best")
plt.tight_layout()
plt.savefig('Plots/2_5.pdf')

#Plot zur Energieerhaltung
plt.figure(6)
plt.plot (tb,xb, 'b-', label='x')
plt.plot (tb,vb, 'r-', label='v')
plt.plot (tb,Ebp, 'g-', label=' Potentielle Energie/m')
plt.plot (tb,Ebk, 'y-', label=' Kinetische Energie/m')
plt.plot (tb,Ebk+Ebp, 'k-', label=' Gesamte Energie/m')
plt.grid()
plt.xlabel('t')
plt.ylabel('')
plt.legend(loc="best")
plt.tight_layout()
plt.savefig('Plots/2_b.pdf')

#Plot zum Energieverhältniss
plt.figure(9)
plt.plot (tb,(Ebk+Ebp)/(Ebk[0]+Ebp[0]), 'g-', label='Energieverhältniss')
plt.grid()
plt.xlabel('t')
plt.ylabel('')
plt.legend(loc="best")
plt.tight_layout()
plt.savefig('Plots/2_b_rel.pdf')


#Lineare Ausgleichsgerade
def f(x, a,b):
    return a*x+b

params_l_ab, covv_l_ab = curve_fit(f, tc, l_ab)
covv_l_ab = np.sqrt(np.diag(covv_l_ab))
print('Steigung der Laufzeit von AB ist ',params_l_ab[0],'pm',covv_l_ab[0])
print('Achsenabschnitt der Laufzeit von AB ist ',params_l_ab[1],'pm',covv_l_ab[1])

params_l_rk, covv_l_rk = curve_fit(f, tc, l_rk)
covv_l_rk = np.sqrt(np.diag(covv_l_rk))
print('Steigung der Laufzeit von RK ist ',params_l_rk[0],'pm',covv_l_rk[0])
print('Achsenabschnitt der Laufzeit von RK ist ',params_l_rk[1],'pm',covv_l_rk[1])

params_f_ab, covv_f_ab = curve_fit(f, tc, f_ab)
covv_f_ab = np.sqrt(np.diag(covv_f_ab))
print('Steigung der Funktionsaufrufe von AB ist ',params_f_ab[0],'pm',covv_f_ab[0])
print('Achsenabschnitt der Funktionsaufrufe von AB ist ',params_f_ab[1],'pm',covv_f_ab[1])

params_f_rk, covv_f_rk = curve_fit(f, tc, f_rk)
covv_f_rk = np.sqrt(np.diag(covv_f_rk))
print('Steigung der Funktionsaufrufe von RK ist ',params_f_rk[0],'pm',covv_f_rk[0])
print('Achsenabschnitt der Funktionsaufrufe von RK ist ',params_f_rk[1],'pm',covv_f_rk[1])

lc=np.linspace(10,1000,10000)

#Plot der Laufzeiten
plt.figure(7)
plt.plot (tc,l_ab, 'rx', label='Adams Bashforth')
plt.plot (tc,l_rk, 'bx', label='Runge Kutta')
plt.plot (lc, f(lc, *params_l_ab), 'r-', label='Lineare Ausgleichsgerade')
plt.plot (lc, f(lc, *params_l_rk), 'b-', label='Lineare Ausgleichsgerade')
plt.grid()
plt.xlabel('t')
plt.xscale('log')
plt.yscale('log')
plt.ylabel('Laufzeit in s')
plt.legend(loc="best")
plt.tight_layout()
plt.savefig('Plots/2_c_laufzeit.pdf')

#Plot der Funktionsaufrufe
plt.figure(8)
plt.plot (tc,f_ab, 'rx', label='Adams Bashforth')
plt.plot (tc,f_rk, 'bx', label='Runge Kutta')
plt.plot (lc, f(lc, *params_f_ab), 'r-', label='Lineare Ausgleichsgerade')
plt.plot (lc, f(lc, *params_f_rk), 'b-', label='Lineare Ausgleichsgerade')
plt.grid()
plt.xlabel('t')
plt.ylabel('Funktionsaufrufe')
plt.xscale('log')
plt.yscale('log')
plt.legend(loc="best")
plt.tight_layout()
plt.savefig('Plots/2_c_funktionsaufrufe.pdf')
