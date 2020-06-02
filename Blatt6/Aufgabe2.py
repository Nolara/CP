import numpy as np
from numpy import exp
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from uncertainties import ufloat
from scipy.stats import sem

x1,y1,f1,abw1= np.genfromtxt('./Data/2_1.txt', unpack=True) #import data
x2,y2,f2,abw2= np.genfromtxt('./Data/2_2.txt', unpack=True)
x3,y3,f3,abw3= np.genfromtxt('./Data/2_3.txt', unpack=True)

xlist = np.linspace(-1.1, 1.1, 100)
ylist = np.linspace(0.9, 1.1, 100)
X, Y = np.meshgrid(xlist, ylist)
Z = ((1-X)**2+100*(Y-X**2)**2)


plt.figure(1)
fig,ax=plt.subplots(1,1)
cp = ax.contour(X, Y, Z)
plt.plot (x1,y1, 'rx', label='FFT')
plt.plot (x2,y2, 'rx', label='FFT')
plt.plot (x3,y3, 'rx', label='FFT')
fig.colorbar(cp) # Add a colorbar to a plot
ax.set_title('Filled Contours Plot')
ax.set_ylabel('y (cm)')
plt.tight_layout()
plt.savefig('Plots/2.pdf')
