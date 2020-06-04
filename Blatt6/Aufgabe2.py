import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

k,x1,x2,f  = np.genfromtxt('./Data/two_a.txt', unpack=True)

plt.plot(k,f,'k.')
plt.xlabel('k')
plt.ylabel('f(xmin_k)')
#plt.ylim(-1,100)
plt.savefig("Plots/plot_A2_a.pdf")
#plt.show()
plt.close()

plt.plot(k,f,'k.')
plt.xlabel('k')
plt.ylabel('f(xmin_k)')
plt.ylim(-1,6)
plt.axhline(0)
plt.savefig("Plots/plot_A2_a1.pdf")
#plt.show()
plt.close()

k,x1,x2,f  = np.genfromtxt('./Data/two_b.txt', unpack=True)

plt.plot(k,f,'k.')
plt.xlabel('k')
plt.ylabel('f(xmin_k)')
#plt.ylim(-1,100)
plt.savefig("Plots/plot_A2_b.pdf")
#plt.show()
plt.close()

plt.plot(k,f,'k.')
plt.xlabel('k')
plt.ylabel('f(xmin_k)')
plt.ylim(-1,6)
plt.axhline(0)
plt.savefig("Plots/plot_A2_b1.pdf")
#plt.show()
plt.close()


k,x1,x2,f  = np.genfromtxt('./Data/two_c.txt', unpack=True)

plt.plot(k,f,'k.')
plt.xlabel('k')
plt.ylabel('f(xmin_k)')
#plt.ylim(-1,100)
plt.savefig("Plots/plot_A2_c.pdf")
#plt.show()
plt.close()

plt.plot(k,f,'k.')
plt.xlabel('k')
plt.ylabel('f(xmin_k)')
plt.ylim(-1,6)
plt.axhline(0)
plt.savefig("Plots/plot_A2_c1.pdf")
#plt.show()
plt.close()


