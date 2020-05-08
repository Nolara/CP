import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

N, M, LU, x  = np.genfromtxt('./Data/.txt', unpack=True)



N_M  =np.zeros(len(N) -3)           ### Ignorieren der ersten paar Werte, da diese noch stark von
M_kl  =np.zeros(len(N) -3)          ### den eigentlichen Werten abweichen und daher den Fitt verfälschen

N_lu  =np.zeros(len(N) -3)
LU_kl =np.zeros(len(N) -3)

N_x  =np.zeros(len(N) -6)
x_kl  =np.zeros(len(N) -6)

i=0
while(i<len(N)-3):
    N_M[i]=N[i+3]
    N_lu[i]=N[i+3]
    M_kl[i]=M[i+3]
    LU_kl[i]=LU[i+3]
    i=i+1
i=0
while(i<len(N)-6):
    N_x[i]=N[i+6]
    x_kl[i]=x[i+6]
    i=i+1


def f(x,a,b):
    return a*x**b

n=np.linspace(8,5*10**4)

print("ax^b")
print("  ")
print('Erstellung von random NxN Matrix')
parM, cov = curve_fit(f,N_M,M_kl)
errM=np.sqrt(np.diag(cov))
print("a ",parM[0],"\pm", errM[0])
print("b ",parM[1],"\pm", errM[1])
plt.plot(n,f(n,*parM),label=r"Fit M")
plt.plot(N,M,'rx',label=r"Laufzeit M")


print("  ")
print('LU-Zerlegung von random NxN Matrix')
parLU, cov = curve_fit(f,N_lu,LU_kl)
errLU=np.sqrt(np.diag(cov))
print("a ",parLU[0],"\pm", errLU[0])
print("b ",parLU[1],"\pm", errLU[1])
plt.plot(n,f(n,*parLU),label=r"Fit LU")
plt.plot(N,LU,'yx',label=r"Laufzeit LU")


print("  ")
print('Lösung von random NxN Matrix')
parx, cov = curve_fit(f,N_x,x_kl)
errx=np.sqrt(np.diag(cov))
print("a ",parx[0],"\pm", errx[0])
print("b ",parx[1],"\pm", errx[1])
plt.plot(n,f(n,*parx),label=r"Fit x")
plt.plot(N,x,'gx',label=r"Laufzeit x")

plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'N')
plt.ylabel(r't [ns]')
plt.legend()
plt.savefig("./Plots/plot2.jpg")
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
