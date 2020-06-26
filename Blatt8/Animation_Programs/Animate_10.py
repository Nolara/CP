import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.animation import FFMpegWriter
writer = FFMpegWriter(fps=15, metadata=dict(artist='Me'), bitrate=1800)
x1,y1= np.genfromtxt('./Data/Werte_10.txt', unpack=True)


l=int(len(x1)/16)

x1_data=[]
y1_data=[]
x2_data=[]
y2_data=[]
x3_data=[]
y3_data=[]
x4_data=[]
y4_data=[]
x5_data=[]
y5_data=[]
x6_data=[]
y6_data=[]
x7_data=[]
y7_data=[]
x8_data=[]
y8_data=[]
x9_data=[]
y9_data=[]
x10_data=[]
y10_data=[]
x11_data=[]
y11_data=[]
x12_data=[]
y12_data=[]
x13_data=[]
y13_data=[]
x14_data=[]
y14_data=[]
x15_data=[]
y15_data=[]
x16_data=[]
y16_data=[]
fig, ax = plt.subplots()
ax.set_ylim(0,8)
ax.set_xlim(0,8)
line1, =ax.plot(x1[0],y1[0], 'o', color='maroon')
line2, =ax.plot(x1[1],y1[1], 'o', color='red')
line3, =ax.plot(x1[2],y1[2], 'o', color='darkorange')
line4, =ax.plot(x1[3],y1[3], 'o', color='gold')
line5, =ax.plot(x1[4],y1[4], 'o', color='yellow')
line6, =ax.plot(x1[5],y1[5], 'o', color='greenyellow')
line7, =ax.plot(x1[6],y1[6], 'o', color='limegreen')
line8, =ax.plot(x1[7],y1[7], 'o', color='aquamarine')
line9, =ax.plot(x1[8],y1[8], 'o', color='aqua')
line10, =ax.plot(x1[9],y1[9], 'o', color='deepskyblue')
line11, =ax.plot(x1[10],y1[10], 'o', color='royalblue')
line12, =ax.plot(x1[11],y1[11], 'o', color='slateblue')
line13, =ax.plot(x1[12],y1[12], 'o', color='blueviolet')
line14, =ax.plot(x1[13],y1[13], 'o', color='violet')
line15, =ax.plot(x1[14],y1[14], 'o', color='magenta')
line16, =ax.plot(x1[15],y1[15], 'o', color='deeppink')

def animate(i):
    x1_data=(x1[16*i])
    y1_data=(y1[16*i])
    line1.set_xdata(x1_data)
    line1.set_ydata(y1_data)
    x2_data=(x1[16*i+1])
    y2_data=(y1[16*i+1])
    line2.set_xdata(x2_data)
    line2.set_ydata(y2_data)
    x3_data=(x1[16*i+2])
    y3_data=(y1[16*i+2])
    line3.set_xdata(x3_data)
    line3.set_ydata(y3_data)
    x4_data=(x1[16*i+3])
    y4_data=(y1[16*i+3])
    line4.set_xdata(x4_data)
    line4.set_ydata(y4_data)
    x5_data=(x1[16*i+4])
    y5_data=(y1[16*i+4])
    line5.set_xdata(x5_data)
    line5.set_ydata(y5_data)
    x6_data=(x1[16*i+5])
    y6_data=(y1[16*i+5])
    line6.set_xdata(x6_data)
    line6.set_ydata(y6_data)
    x7_data=(x1[16*i+6])
    y7_data=(y1[16*i+6])
    line7.set_xdata(x7_data)
    line7.set_ydata(y7_data)
    x8_data=(x1[16*i+7])
    y8_data=(y1[16*i+7])
    line8.set_xdata(x8_data)
    line8.set_ydata(y8_data)
    x9_data=(x1[16*i+8])
    y9_data=(y1[16*i+8])
    line9.set_xdata(x9_data)
    line9.set_ydata(y9_data)
    x10_data=(x1[16*i+9])
    y10_data=(y1[16*i+9])
    line10.set_xdata(x10_data)
    line10.set_ydata(y10_data)
    x11_data=(x1[16*i+10])
    y11_data=(y1[16*i+10])
    line11.set_xdata(x11_data)
    line11.set_ydata(y11_data)
    x12_data=(x1[16*i+11])
    y12_data=(y1[16*i+11])
    line12.set_xdata(x12_data)
    line12.set_ydata(y12_data)
    x13_data=(x1[16*i+12])
    y13_data=(y1[16*i+12])
    line13.set_xdata(x13_data)
    line13.set_ydata(y13_data)
    x14_data=(x1[16*i+13])
    y14_data=(y1[16*i+13])
    line14.set_xdata(x14_data)
    line14.set_ydata(y14_data)
    x15_data=(x1[16*i+14])
    y15_data=(y1[16*i+14])
    line15.set_xdata(x15_data)
    line15.set_ydata(y15_data)
    x16_data=(x1[16*i+15])
    y16_data=(y1[16*i+15])
    line16.set_xdata(x16_data)
    line16.set_ydata(y16_data)

    return line1, line2, line3, line4, line5, line6, line7, line8, line9, line10, line11, line12, line13, line14, line15, line16,


anim = FuncAnimation(fig, animate, frames=np.arange(1,1000,1 ), save_count=100)

anim.save('Animations/Animation_10.mp4', writer=writer)
