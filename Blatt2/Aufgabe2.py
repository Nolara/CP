import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

x1, x2 ,y,z, = np.genfromtxt('./Data/two_data.txt', unpack=True) #import data

l=len(x1)
x0 = np.linspace(0, l, l)

#Plot Verfahren aus a)
plt.figure(1)
plt.plot (2**x0, x1, 'k.', label='Zufallsmatrix')
plt.plot (2**x0, x2, 'y.', label='Zufallsvektor')
plt.plot (2**x0, y, 'r.', label='LU Zerlegung')
plt.plot (2**x0, z, 'b.', label='Lösen')
plt.grid()
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'N')
plt.ylabel(r't [ns]')
plt.legend()
plt.savefig("./Data/plot2.jpg")
plt.close()


print("   ")
print("Abschätzung der Laufzeit für N=1.000.000")
N=1000000
LaufzeitM=f(N,*parM)
LaufzeitLU=f(N,*parLU)
Laufzeitx=f(N,*parx)
print("t_M",np.round(LaufzeitM*10**(-9),2),"s = ",np.round(LaufzeitM*10**(-9)/60/60,2),"h = ",np.round(LaufzeitM*10**(-9)/60/60/24/365,2),"a")
print("t_LU",np.round(LaufzeitLU*10**(-9),2),"s = ",np.round(LaufzeitLU*10**(-9)/60/60,2),"h = ",np.round(LaufzeitLU*10**(-9)/60/60/24/365,2),"a")
print("t_x",np.round(Laufzeitx*10**(-9),2),"s = ",np.round(Laufzeitx*10**(-9)/60/60,2),"h = ",np.round(Laufzeitx*10**(-9)/60/60/24/365,2),"a")

t=LaufzeitM+LaufzeitLU+Laufzeitx
print("Laufzeit insgesammt")
print("t",np.round(t*10**(-9),2),"s = ",np.round(t*10**(-9)/60/60,2),"h = ",np.round(t*10**(-9)/60/60/24,2),"d",np.round(t*10**(-9)/60/60/24/365,2),"a")
