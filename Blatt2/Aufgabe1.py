import numpy as np
from numpy import exp
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from uncertainties import ufloat
from scipy.stats import sem
import seaborn as sns

a_Original = np.genfromtxt('./Data/one_original.txt') #import data
plt.imshow(a_Original, cmap='gray', interpolation='nearest')
plt.axis('off')
plt.tight_layout()
plt.savefig('Plots/1_Original.pdf')

a_10 = np.genfromtxt('./Data/Reduktion_10.txt') #import data
plt.imshow(a_10, cmap='gray', interpolation='nearest')
plt.axis('off')
plt.tight_layout()
plt.savefig('Plots/1_10.pdf')

a_20 = np.genfromtxt('./Data/Reduktion_20.txt') #import data
plt.imshow(a_20, cmap='gray', interpolation='nearest')
plt.axis('off')
plt.tight_layout()
plt.savefig('Plots/1_20.pdf')

a_50 = np.genfromtxt('./Data/Reduktion_50.txt') #import data
plt.imshow(a_50, cmap='gray', interpolation='nearest')
plt.axis('off')
plt.tight_layout()
plt.savefig('Plots/1_50.pdf')
