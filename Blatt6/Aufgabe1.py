import numpy as np
from numpy import exp
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from uncertainties import ufloat
from scipy.stats import sem

xrg,yrg,frg,abwrg= np.genfromtxt('./Data/1_gradient_rosenbrock.txt', unpack=True) #import data
xrc,yrc,frc,abwrc= np.genfromtxt('./Data/1_conjugate_rosenbrock.txt', unpack=True)

xlist = np.linspace(-1.1, 1.1, 100)
ylist = np.linspace(0.9, 1.1, 100)
X, Y = np.meshgrid(xlist, ylist)
Z = ((1-X)**2+100*(Y-X**2)**2)

plt.figure(1)
fig,ax=plt.subplots(1,1)
cp = ax.contour(X, Y, Z)
plt.plot (xrg,yrg, 'rx', label='Gradient')
plt.plot (xrc,yrc, 'bx', label='Conjugate Gradient')
fig.colorbar(cp) # Add a colorbar to a plot
#ax.set_title('Filled Contours Plot')
ax.set_xlabel('x1')
ax.set_ylabel('x2')
plt.tight_layout()
plt.savefig('Plots/1_Rosenbrock.pdf')

plt.figure(2)
fig,ax=plt.subplots(1,1)
cp = ax.contour(X, Y, Z)
plt.plot (xrg,yrg, 'rx', label='FFT')
fig.colorbar(cp) # Add a colorbar to a plot
ax.set_title('Filled Contours Plot')
#ax.set_xlabel('x (cm)')
ax.set_xlabel('x1')
ax.set_ylabel('x2')
plt.tight_layout()
plt.savefig('Plots/1_b.pdf')
