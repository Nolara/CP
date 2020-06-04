import numpy as np
from numpy import exp
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from uncertainties import ufloat
from scipy.stats import sem

k1,x1,y1,f1,abw1= np.genfromtxt('./Data/2_1.txt', unpack=True) #import data
k2,x2,y2,f2,abw2= np.genfromtxt('./Data/2_2.txt', unpack=True)
k3,x3,y3,f3,abw3= np.genfromtxt('./Data/2_3.txt', unpack=True)
