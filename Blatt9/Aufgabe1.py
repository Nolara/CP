import numpy as np
import matplotlib.pyplot as plt

H, m = np.genfromtxt('./Data/one.txt',unpack='True')

plt.plot(H,m,label=r'Numerische Rechnung', color='darkorange', linewidth=5.0)
plt.plot(H, np.tanh(H), label='Analytische Rechnung', color='darkblue')

plt.xlabel(r'Äußeres Magnetfeld $H$')
plt.ylabel(r'Magnetisierung $m$')
plt.grid()
plt.legend(loc='best')
plt.savefig("Plots/plot_1.pdf")
