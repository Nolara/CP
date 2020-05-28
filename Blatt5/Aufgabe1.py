import numpy as np
from numpy import exp
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from uncertainties import ufloat
from scipy.stats import sem

Re_direkt, Im_direkt, Re_FFT, Im_FFT= np.genfromtxt('./Data/1_1.txt', unpack=True) #import data
abs_exp, Re_exp,Im_exp,abs_recht,Re_recht,Im_recht= np.genfromtxt('./Data/1_2.txt', unpack=True)

k1=np.linspace(0,8,8)

plt.figure(1)
plt.plot (k1,Im_direkt[0:8], 'b+', label='Direkte Rechnung')
plt.plot (k1,Im_FFT[0:8], 'rx', label='FFT')
plt.grid()
plt.xlabel('x')
plt.ylabel('Im')
plt.legend(loc="best")
plt.tight_layout()
plt.savefig('Plots/1_1_im_3.pdf')

plt.figure(8)
plt.plot (k1,Re_direkt[0:8], 'b+', label='Direkte Rechnung')
plt.plot (k1,Re_FFT[0:8], 'rx', label='FFT')
plt.grid()
plt.xlabel('x')
plt.ylabel('Re')
plt.legend(loc="best")
plt.tight_layout()
plt.savefig('Plots/1_1_im_3.pdf')

k12=np.linspace(0,16,16)

plt.figure(9)
plt.plot (k12,Im_direkt[8:], 'b+', label='Direkte Rechnung')
plt.plot (k12,Im_FFT[8:], 'rx', label='FFT')
plt.grid()
plt.xlabel('x')
plt.ylabel('Im')
plt.legend(loc="best")
plt.tight_layout()
plt.savefig('Plots/1_1_im_4.pdf')

plt.figure(10)
plt.plot (k12,Re_direkt[8:], 'b+', label='Direkte Rechnung')
plt.plot (k12,Re_FFT[8:], 'rx', label='FFT')
plt.grid()
plt.xlabel('x')
plt.ylabel('Re')
plt.legend(loc="best")
plt.tight_layout()
plt.savefig('Plots/1_1_im_4.pdf')


ka=np.linspace(-10,10,1000)
k=np.linspace(-10,10,128)


plt.figure(2)
plt.plot (k,Re_exp, 'rx', label='FFT')
plt.plot (k,Re_exp, 'r-')
plt.plot(ka, 1/(2*np.pi)**(1/2)*np.exp(-ka**2*2), 'g-', label='Analytische Lösung')
plt.grid()
plt.xlabel('k')
plt.ylabel('Re')
plt.legend(loc="best")
plt.tight_layout()
plt.savefig('Plots/1_2_Re.pdf')

plt.figure(3)
plt.plot (k,Im_exp, 'rx', label='FFT')
plt.plot (k,Im_exp, 'r-')
plt.plot(ka,0*ka, 'g-', label='Analytische Lösung')
plt.grid()
plt.xlabel('k')
plt.ylabel('Im')
plt.legend(loc="best")
plt.tight_layout()
plt.savefig('Plots/1_2_Im.pdf')


plt.figure(6)
plt.plot (k,abs_exp, 'rx', label='FFT')
plt.plot (k,abs_exp, 'r-')
plt.plot(ka, 1/(2*np.pi)**(1/2)*np.exp(-ka**2/2), 'g-', label='Analytische Lösung')
plt.grid()
plt.xlabel('k')
plt.ylabel('F(k)')
plt.legend(loc="best")
plt.tight_layout()
plt.savefig('Plots/1_2_Abs.pdf')


ka2=np.linspace(-np.pi,np.pi,1000)
k2=np.linspace(-np.pi,np.pi,128)



plt.figure(4)
plt.plot (k2,Re_recht, 'rx', label='FFT')
plt.plot (k2,Re_recht, 'r-')
plt.plot(ka2,0*ka2 , 'g-', label='Analytische Lösung')
plt.grid()
plt.xlabel('k')
plt.ylabel('Re')
plt.legend(loc="best")
plt.tight_layout()
plt.savefig('Plots/1_3_Re.pdf')

plt.figure(5)
plt.plot (k2,Im_recht, 'rx', label='FFT')
plt.plot (k2,Im_recht, 'r-')
plt.plot(ka2,1/ka2*(2-2*np.cos(ka2*np.pi)), 'g-', label='Analytische Lösung')
plt.grid()
plt.xlabel('k')
plt.ylabel('Im')
plt.legend(loc="best")
plt.tight_layout()
plt.savefig('Plots/1_3_Im.pdf')

plt.figure(7)
plt.plot (k2,abs_recht, 'rx', label='FFT')
plt.plot (k2,abs_recht, 'r-')
plt.plot(ka2,1/ka2*(2-2*np.cos(ka2*np.pi)), 'g-', label='Analytische Lösung')
plt.grid()
plt.xlabel('k')
plt.ylabel('F(k)')
plt.legend(loc="best")
plt.tight_layout()
plt.savefig('Plots/1_3_Abs.pdf')
