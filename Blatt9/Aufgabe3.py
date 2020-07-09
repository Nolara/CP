import numpy as np
import matplotlib.pyplot as plt

x,y= np.genfromtxt('./Data/3_start.txt', unpack=True)
x_random, y_random= np.genfromtxt('./Data/3_start_random.txt', unpack=True)

x_0, y_0= np.genfromtxt('./Data/3_final_0.txt', unpack=True)
x_1, y_1= np.genfromtxt('./Data/3_final_1.txt', unpack=True)
x_2, y_2= np.genfromtxt('./Data/3_final_2.txt', unpack=True)
x_3, y_3= np.genfromtxt('./Data/3_final_3.txt', unpack=True)
x_4, y_4= np.genfromtxt('./Data/3_final_4.txt', unpack=True)
x_5, y_5= np.genfromtxt('./Data/3_final_5.txt', unpack=True)
x_6, y_6= np.genfromtxt('./Data/3_final_6.txt', unpack=True)
x_7, y_7= np.genfromtxt('./Data/3_final_7.txt', unpack=True)
x_8, y_8= np.genfromtxt('./Data/3_final_8.txt', unpack=True)
x_9, y_9= np.genfromtxt('./Data/3_final_9.txt', unpack=True)
x_10, y_10= np.genfromtxt('./Data/3_final_10.txt', unpack=True)
x_11, y_11= np.genfromtxt('./Data/3_final_11.txt', unpack=True)
x_12, y_12= np.genfromtxt('./Data/3_final_12.txt', unpack=True)
x_13, y_13= np.genfromtxt('./Data/3_final_13.txt', unpack=True)
x_14, y_14= np.genfromtxt('./Data/3_final_14.txt', unpack=True)

L_0= np.genfromtxt('./Data/3_0.txt', unpack=True)
L_1= np.genfromtxt('./Data/3_1.txt', unpack=True)
L_2= np.genfromtxt('./Data/3_2.txt', unpack=True)
L_3= np.genfromtxt('./Data/3_3.txt', unpack=True)
L_4= np.genfromtxt('./Data/3_4.txt', unpack=True)
L_5= np.genfromtxt('./Data/3_5.txt', unpack=True)
L_6= np.genfromtxt('./Data/3_6.txt', unpack=True)
L_7= np.genfromtxt('./Data/3_7.txt', unpack=True)
L_8= np.genfromtxt('./Data/3_8.txt', unpack=True)
L_9= np.genfromtxt('./Data/3_9.txt', unpack=True)
L_10= np.genfromtxt('./Data/3_10.txt', unpack=True)
L_11= np.genfromtxt('./Data/3_11.txt', unpack=True)
L_12= np.genfromtxt('./Data/3_12.txt', unpack=True)
L_13= np.genfromtxt('./Data/3_13.txt', unpack=True)
L_14= np.genfromtxt('./Data/3_14.txt', unpack=True)

x=np.append(x, x[0])
y=np.append(y, y[0])
x_random=np.append(x_random, x_random[0])
y_random=np.append(y_random, y_random[0])
x_0=np.append(x_0, x_0[0])
y_0=np.append(y_0, y_0[0])
x_1=np.append(x_1, x_1[0])
y_1=np.append(y_1, y_1[0])
x_2=np.append(x_2, x_2[0])
y_2=np.append(y_2, y_2[0])
x_3=np.append(x_3, x_3[0])
y_3=np.append(y_3, y_3[0])
x_4=np.append(x_4, x_4[0])
y_4=np.append(y_4, y_4[0])
x_5=np.append(x_5, x_5[0])
y_5=np.append(y_5, y_5[0])
x_6=np.append(x_6, x_6[0])
y_6=np.append(y_6, y_6[0])
x_7=np.append(x_7, x_7[0])
y_7=np.append(y_7, y_7[0])
x_8=np.append(x_8, x_8[0])
y_8=np.append(y_8, y_8[0])
x_9=np.append(x_9, x_9[0])
y_9=np.append(y_9, y_9[0])
x_10=np.append(x_10, x_10[0])
y_10=np.append(y_10, y_10[0])
x_11=np.append(x_11, x_11[0])
y_11=np.append(y_11, y_11[0])
x_12=np.append(x_12, x_12[0])
y_12=np.append(y_12, y_12[0])
x_13=np.append(x_13, x_13[0])
y_13=np.append(y_13, y_13[0])
x_14=np.append(x_14, x_14[0])
y_14=np.append(y_14, y_14[0])


plt.figure(15)
plt.plot (x,y, 'rx')
plt.plot (x,y, 'r-')
plt.grid()
plt.xlabel('x')
plt.ylabel('y')
plt.tight_layout()
plt.savefig('Plots/3_start.pdf')

plt.figure(16)
plt.plot (x_random,y_random, 'rx')
plt.plot (x_random,y_random, 'r-')
plt.grid()
plt.xlabel('x')
plt.ylabel('y')
plt.tight_layout()
plt.savefig('Plots/3_start_random.pdf')

plt.figure(0)
plt.plot (x_0,y_0, 'rx')
plt.plot (x_0,y_0, 'r-')
plt.grid()
plt.xlabel('x')
plt.ylabel('y')
plt.tight_layout()
plt.savefig('Plots/3_final_0.pdf')

plt.figure(1)
plt.plot (x_1,y_1, 'rx')
plt.plot (x_1,y_1, 'r-')
plt.grid()
plt.xlabel('x')
plt.ylabel('y')
plt.tight_layout()
plt.savefig('Plots/3_final_1.pdf')

plt.figure(2)
plt.plot (x_2,y_2, 'rx')
plt.plot (x_2,y_2, 'r-')
plt.grid()
plt.xlabel('x')
plt.ylabel('y')
plt.tight_layout()
plt.savefig('Plots/3_final_2.pdf')

plt.figure(3)
plt.plot (x_3,y_3, 'rx')
plt.plot (x_3,y_3, 'r-')
plt.grid()
plt.xlabel('x')
plt.ylabel('y')
plt.tight_layout()
plt.savefig('Plots/3_final_3.pdf')

plt.figure(4)
plt.plot (x_4,y_4, 'rx')
plt.plot (x_4,y_4, 'r-')
plt.grid()
plt.xlabel('x')
plt.ylabel('y')
plt.tight_layout()
plt.savefig('Plots/3_final_4.pdf')

plt.figure(5)
plt.plot (x_5,y_5, 'rx')
plt.plot (x_5,y_5, 'r-')
plt.grid()
plt.xlabel('x')
plt.ylabel('y')
plt.tight_layout()
plt.savefig('Plots/3_final_5.pdf')

plt.figure(6)
plt.plot (x_6,y_6, 'rx')
plt.plot (x_6,y_6, 'r-')
plt.grid()
plt.xlabel('x')
plt.ylabel('y')
plt.tight_layout()
plt.savefig('Plots/3_final_6.pdf')

plt.figure(7)
plt.plot (x_7,y_7, 'rx')
plt.plot (x_7,y_7, 'r-')
plt.grid()
plt.xlabel('x')
plt.ylabel('y')
plt.tight_layout()
plt.savefig('Plots/3_final_7.pdf')

plt.figure(8)
plt.plot (x_8,y_8, 'rx')
plt.plot (x_8,y_8, 'r-')
plt.grid()
plt.xlabel('x')
plt.ylabel('y')
plt.tight_layout()
plt.savefig('Plots/3_final_8.pdf')

plt.figure(9)
plt.plot (x_9,y_9, 'rx')
plt.plot (x_9,y_9, 'r-')
plt.grid()
plt.xlabel('x')
plt.ylabel('y')
plt.tight_layout()
plt.savefig('Plots/3_final_9.pdf')

plt.figure(10)
plt.plot (x_10,y_10, 'rx')
plt.plot (x_10,y_10, 'r-')
plt.grid()
plt.xlabel('x')
plt.ylabel('y')
plt.tight_layout()
plt.savefig('Plots/3_final_10.pdf')

plt.figure(11)
plt.plot (x_11,y_11, 'rx')
plt.plot (x_11,y_11, 'r-')
plt.grid()
plt.xlabel('x')
plt.ylabel('y')
plt.tight_layout()
plt.savefig('Plots/3_final_11.pdf')

plt.figure(12)
plt.plot (x_12,y_12, 'rx')
plt.plot (x_12,y_12, 'r-')
plt.grid()
plt.xlabel('x')
plt.ylabel('y')
plt.tight_layout()
plt.savefig('Plots/3_final_12.pdf')

plt.figure(13)
plt.plot (x_13,y_13, 'rx')
plt.plot (x_13,y_13, 'r-')
plt.grid()
plt.xlabel('x')
plt.ylabel('y')
plt.tight_layout()
plt.savefig('Plots/3_final_13.pdf')

plt.figure(14)
plt.plot (x_14,y_14, 'rx')
plt.plot (x_14,y_14, 'r-')
plt.grid()
plt.xlabel('x')
plt.ylabel('y')
plt.tight_layout()
plt.savefig('Plots/3_final_14.pdf')

n1=np.linspace(0,len(L_0),len(L_0))
plt.figure(17)
plt.plot (n1,L_0, '-', color='slateblue', label='s=10')
plt.plot (n1,L_1, '-', color='limegreen', label='s=100')
plt.plot (n1,L_2, '-', color='gold', label='s=1000')
plt.plot (n1,L_3, '-', color='darkorange', label='s=10000')
plt.plot (n1,L_4, '-', color='red', label='s=100000')
plt.legend(loc="best")
plt.grid()
plt.xlabel('N')
plt.ylabel('L')
plt.tight_layout()
plt.savefig('Plots/3_0_9.pdf')

n1=np.linspace(0,len(L_5),len(L_5))
plt.figure(18)
plt.plot (n1,L_5, '-', color='slateblue', label='s=10')
plt.plot (n1,L_6, '-', color='limegreen', label='s=100')
plt.plot (n1,L_7, '-', color='gold', label='s=1000')
plt.plot (n1,L_8, '-', color='darkorange', label='s=10000')
plt.plot (n1,L_9, '-', color='red', label='s=100000')
plt.legend(loc="best")
plt.grid()
plt.xlabel('N')
plt.ylabel('L')
plt.tight_layout()
plt.savefig('Plots/3_0_99.pdf')

n1=np.linspace(0,len(L_10),len(L_10))
plt.figure(19)
plt.plot (n1,L_10, '-', color='slateblue', label='s=10')
plt.plot (n1,L_11, '-', color='limegreen', label='s=100')
plt.plot (n1,L_12, '-', color='gold', label='s=1000')
plt.plot (n1,L_13, '-', color='darkorange', label='s=10000')
plt.plot (n1,L_14, '-', color='red', label='s=100000')
plt.legend(loc="best")
plt.grid()
plt.xlabel('N')
plt.ylabel('L')
plt.tight_layout()
plt.savefig('Plots/3_0_999.pdf')
