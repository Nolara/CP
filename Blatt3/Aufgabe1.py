import numpy as np
from numpy import exp
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from uncertainties import ufloat
from scipy.stats import sem


d = np.genfromtxt('./Data/one_Dichte.txt') #import data
#Logarithmieren zur besseren Sichtbarkeit
dlog=np.log(d)
plt.tight_layout()
extent=(1,50,-20,20)
plt.imshow(dlog, interpolation='nearest', extent=extent, origin='lower')
cbar = plt.colorbar(ticks=[np.log(10**(-33)),np.log(10**(-22)),np.log(10**(-11)),np.log(0.1)])
cbar.ax.set_yticklabels(['10^(-33)','10^(-22)','10^(-11)', '10^(-1)'])
plt.xticks(np.arange(0,50,5))
plt.yticks(np.arange(-20,20+2,2))
plt.xlabel('i')
plt.ylabel('epsilon')
plt.savefig('Plots/1.pdf')
