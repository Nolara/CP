import numpy as np
from numpy import exp
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from uncertainties import ufloat
from scipy.stats import sem
from matplotlib.patches import ConnectionPatch

xrg,yrg,frg,abwrg= np.genfromtxt('./Data/1_gradient_rosenbrock.txt', unpack=True) #import data
xrc,yrc,frc,abwrc= np.genfromtxt('./Data/1_conjugate_rosenbrock.txt', unpack=True)

xb1g,yb1g,fb1g,abwb1g= np.genfromtxt('./Data/1_gradient_b_1.txt', unpack=True)
xb2g,yb2g,fb2g,abwb2g= np.genfromtxt('./Data/1_gradient_b_2.txt', unpack=True)
xb3g,yb3g,fb3g,abwb3g= np.genfromtxt('./Data/1_gradient_b_3.txt', unpack=True)
xb1c,yb1c,fb1c,abwb1c= np.genfromtxt('./Data/1_conjugate_b_1.txt', unpack=True)
xb2c,yb2c,fb2c,abwb2c= np.genfromtxt('./Data/1_conjugate_b_2.txt', unpack=True)
xb3c,yb3c,fb3c,abwb3c= np.genfromtxt('./Data/1_conjugate_b_3.txt', unpack=True)

len_ag=len(xrg)
len_ac=len(xrc)

len_bg1=len(xb1g)
len_bg2=len(xb2g)
len_bc1=len(xb1c)
len_bc2=len(xb2c)

xlista = np.linspace(-2, 2, 100)
ylista = np.linspace(-2, 2, 100)
Xa, Ya = np.meshgrid(xlista, ylista)
Za = (((1-Xa)**2+100*(Ya-Xa**2)**2))

xlistac = np.linspace(-3, 6, 100)
ylistac = np.linspace(-3, 27, 100)
Xac, Yac = np.meshgrid(xlistac, ylistac)
Zac = (((1-Xac)**2+100*(Yac-Xac**2)**2))


xlistb1 = np.linspace(1,2.0, 100)
ylistb1 = np.linspace(1.6, 2.5, 100)
Xb1, Yb1 = np.meshgrid(xlistb1, ylistb1)
Zb1 = ((1+exp(-10*(Xb1*Yb1-3)**2)/(Xb1**2+Yb1**2))**(-1))

xlistb2 = np.linspace(-2,-1, 100)
ylistb2 = np.linspace(-2.5, -1.6, 100)
Xb2, Yb2 = np.meshgrid(xlistb2, ylistb2)
Zb2 =((1+exp(-10*(Xb2*Yb2-3)**2)/(Xb2**2+Yb2**2))**(-1))

xlistb3 = np.linspace(-2,2.5, 100)
ylistb3 = np.linspace(-2, 2, 100)
Xb3, Yb3 = np.meshgrid(xlistb3, ylistb3)
Zb3 = ((1+exp(-10*(Xb3*Yb3-3)**2)/(Xb3**2+Yb3**2))**(-1))

xlistb4 = np.linspace(0.4,0.6, 100)
ylistb4 = np.linspace(0.5, 0.7, 100)
Xb4, Yb4 = np.meshgrid(xlistb4, ylistb4)
Zb4 = ((1+exp(-10*(Xb4*Yb4-3)**2)/(Xb4**2+Yb4**2))**(-1))


plt.figure(1)
fig,ax=plt.subplots(1,1)
cp = ax.contourf(Xa, Ya, Za, cmap='Spectral_r')
plt.plot (xrg,yrg, 'k.')
plt.plot (xrg[len_ag-1],yrg[len_ag-1], 'w.')
plt.plot (xrg[0],yrg[0],'.', color='limegreen')
for x in range(0, 3):
    xyA = (xrg[x], yrg[x])
    xyB = (xrg[x+1], yrg[x+1])
    coordsA = "data"
    coordsB = "data"
    con = ConnectionPatch(xyA, xyB,
                          coordsA, coordsB,
                          arrowstyle ="->",
                          shrinkA = 0, shrinkB = 0,
                          mutation_scale = 20,
                          fc ="w")
    ax.add_artist(con)
fig.colorbar(cp) # Add a colorbar to a plot
#ax.set_title('Filled Contours Plot')
ax.set_xlabel('x1')
ax.set_ylabel('x2')
plt.tight_layout()
plt.savefig('Plots/1_Rosenbrock_grad.pdf')

plt.figure(2)
fig,ax=plt.subplots(1,1)
cp = ax.contourf(Xac, Yac, Zac, cmap='Spectral_r')
plt.plot (xrc,yrc, 'k.')
plt.plot (xrc[185:len_ac-1],yrc[185:len_ac-1], '.', color='magenta')
plt.plot (xrc[len_ac-1],yrc[len_ac-1], 'w.')
plt.plot (xrc[0],yrc[0],'.', color='limegreen')
for x in range(0, 2):
    xyA = (xrc[x], yrc[x])
    xyB = (xrc[x+1], yrc[x+1])
    coordsA = "data"
    coordsB = "data"
    con = ConnectionPatch(xyA, xyB,
                          coordsA, coordsB,
                          arrowstyle ="->",
                          shrinkA = 0, shrinkB = 0,
                          mutation_scale = 20,
                          fc ="w")
    ax.add_artist(con)
for x in range(183, 185):
    xyA = (xrc[x], yrc[x])
    xyB = (xrc[x+1], yrc[x+1])
    coordsA = "data"
    coordsB = "data"
    con = ConnectionPatch(xyA, xyB,
                          coordsA, coordsB,
                          arrowstyle ="->",
                          shrinkA = 0, shrinkB = 0,
                          mutation_scale = 20,
                          fc ="w")
    ax.add_artist(con)
fig.colorbar(cp) # Add a colorbar to a plot
#ax.set_title('Filled Contours Plot')
ax.set_xlabel('x1')
ax.set_ylabel('x2')
plt.tight_layout()
plt.savefig('Plots/1_Rosenbrock_conj.pdf')

plt.figure(3)
fig,ax=plt.subplots(1,1)
cp = ax.contourf(Xb1, Yb1, Zb1, cmap='Spectral_r')
plt.plot (xb1g,yb1g, 'k.')
plt.plot (xb1g[len_bg1-1],yb1g[len_bg1-1], 'w.')
plt.plot (xb1g[0],yb1g[0], '.', color='limegreen')
for x in range(0, 2):
    xyA = (xb1g[x], yb1g[x])
    xyB = (xb1g[x+1], yb1g[x+1])
    coordsA = "data"
    coordsB = "data"
    con = ConnectionPatch(xyA, xyB,
                          coordsA, coordsB,
                          arrowstyle ="->",
                          shrinkA = 0, shrinkB = 0,
                          mutation_scale = 20,
                          fc ="w")
    ax.add_artist(con)
fig.colorbar(cp) # Add a colorbar to a plot
#ax.set_title('Filled Contours Plot')
ax.set_xlabel('x1')
ax.set_ylabel('x2')
plt.tight_layout()
plt.savefig('Plots/1_b1_grad.pdf')

plt.figure(4)
fig,ax=plt.subplots(1,1)
cp = ax.contourf(Xb1, Yb1, Zb1, cmap='Spectral_r')
plt.plot (xb1c,yb1c, 'k.')
plt.plot (xb1c[len_bc1-1],yb1c[len_bc1-1], 'w.')
plt.plot (xb1c[0],yb1c[0], '.', color='limegreen')
for x in range(0, 5):
    xyA = (xb1c[x], yb1c[x])
    xyB = (xb1c[x+1], yb1c[x+1])
    coordsA = "data"
    coordsB = "data"
    con = ConnectionPatch(xyA, xyB,
                          coordsA, coordsB,
                          arrowstyle ="->",
                          shrinkA = 0, shrinkB = 0,
                          mutation_scale = 20,
                          fc ="w")
    ax.add_artist(con)
fig.colorbar(cp) # Add a colorbar to a plot
#ax.set_title('Filled Contours Plot')
ax.set_xlabel('x1')
ax.set_ylabel('x2')
plt.tight_layout()
plt.savefig('Plots/1_b1_conj.pdf')

plt.figure(5)
fig,ax=plt.subplots(1,1)
cp = ax.contourf(Xb2, Yb2, Zb2, cmap='Spectral_r')
plt.plot (xb2g,yb2g, 'k.')
plt.plot (xb2g[len_bg2-1],yb2g[len_bg2-1], 'w.')
plt.plot (xb2g[0],yb2g[0], '.', color='limegreen')
for x in range(0, 2):
    xyA = (xb2g[x], yb2g[x])
    xyB = (xb2g[x+1], yb2g[x+1])
    coordsA = "data"
    coordsB = "data"
    con = ConnectionPatch(xyA, xyB,
                          coordsA, coordsB,
                          arrowstyle ="->",
                          shrinkA = 0, shrinkB = 0,
                          mutation_scale = 20,
                          fc ="w")
    ax.add_artist(con)
fig.colorbar(cp) # Add a colorbar to a plot
#ax.set_title('Filled Contours Plot')
ax.set_xlabel('x1')
ax.set_ylabel('x2')
plt.tight_layout()
plt.savefig('Plots/1_b2_grad.pdf')

plt.figure(6)
fig,ax=plt.subplots(1,1)
cp = ax.contourf(Xb3, Yb3, Zb3, cmap='Spectral_r')
plt.plot (xb2c,yb2c, 'k.')
plt.plot (xb2c[len_bc2-1],yb2c[len_bc2-1], 'w.')
plt.plot (xb2c[0],yb2c[0], '.', color='limegreen')
for x in range(0, 2):
    xyA = (xb2c[x], yb2c[x])
    xyB = (xb2c[x+1], yb2c[x+1])
    coordsA = "data"
    coordsB = "data"
    con = ConnectionPatch(xyA, xyB,
                          coordsA, coordsB,
                          arrowstyle ="->",
                          shrinkA = 0, shrinkB = 0,
                          mutation_scale = 20,
                          fc ="w")
    ax.add_artist(con)
fig.colorbar(cp) # Add a colorbar to a plot
#ax.set_title('Filled Contours Plot')
ax.set_xlabel('x1')
ax.set_ylabel('x2')
plt.tight_layout()
plt.savefig('Plots/1_b2_conj.pdf')

k1=np.linspace(0, len_ag, len_ag)
k2=np.linspace(0, len_ac, len_ac)
k3=np.linspace(0, len_bg1, len_bg1)
k4=np.linspace(0, len_bc1, len_bc1)
k5=np.linspace(0, len_bg2, len_bg2)
k6=np.linspace(0, len_bc2, len_bc2)

plt.figure(7)
plt.plot (k1,abwrg, 'b+', label='Gradienten-Verfahren')
plt.grid()
plt.xlabel('k')
plt.ylabel('Abw')
plt.legend(loc="best")
plt.tight_layout()
plt.savefig('Plots/1_abw_rosenbrock_g.pdf')

plt.figure(8)
plt.plot (k1,abwrg, 'b+', label='Gradienten-Verfahren')
plt.grid()
plt.xlabel('k')
plt.ylabel('Abw')
plt.legend(loc="best")
plt.tight_layout()
plt.savefig('Plots/1_abw_rosenbrock_g.pdf')

plt.figure(9)
plt.plot (k2,abwrc, 'rx', label='Konjugiertes-Gradienten-Verfahren')
plt.grid()
plt.xlabel('k')
plt.ylabel('Abw')
plt.legend(loc="best")
plt.tight_layout()
plt.savefig('Plots/1_abw_rosenbrock_c.pdf')

plt.figure(10)
plt.plot (k3,abwb1g, 'b+', label='Gradienten-Verfahren')
plt.grid()
plt.xlabel('k')
plt.ylabel('$\epsilon_k$')
plt.legend(loc="best")
plt.tight_layout()
plt.savefig('Plots/1_abw_b1_g.pdf')

plt.figure(11)
plt.plot (k4,abwb1c, 'rx', label='Konjugiertes-Gradienten-Verfahren')
plt.grid()
plt.xlabel('k')
plt.ylabel('$\epsilon_k$')
plt.legend(loc="best")
plt.tight_layout()
plt.savefig('Plots/1_abw_b1_c.pdf')

plt.figure(12)
plt.plot (k5,abwb2g, 'b+', label='Gradienten-Verfahren')
plt.grid()
plt.xlabel('k')
plt.ylabel('$\epsilon_k$')
plt.legend(loc="best")
plt.tight_layout()
plt.savefig('Plots/1_abw_b2_g.pdf')

plt.figure(13)
plt.plot (k6,abwb2c, 'rx', label='Konjugiertes-Gradienten-Verfahren')
plt.grid()
plt.xlabel('k')
plt.ylabel('$\epsilon_k$')
plt.legend(loc="best")
plt.tight_layout()
plt.savefig('Plots/1_abw_b2_c.pdf')
