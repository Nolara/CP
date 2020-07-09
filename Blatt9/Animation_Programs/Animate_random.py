import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import FFMpegWriter

Y1 = np.genfromtxt('./Data/2_1_random.txt')

nSeconds1 = 10
col1=len(Y1[0])
row1=len(Y1)//col1
fps1=row1//nSeconds1

fig1 = plt.figure( figsize=(8,8) )

im1 = plt.imshow(Y1[0:col1,0:col1])

def animate_func1(i):
    im1.set_array(Y1[i*col1:(i+1)*col1,0:col1].copy())
    return [im1]

anim1 = animation.FuncAnimation(fig1, animate_func1, frames = row1, interval = 1000 /fps1 )

anim1.save('Animations/Animation_1_random.mp4', fps=fps1, extra_args=['-vcodec', 'libx264'])

Y3 = np.genfromtxt('./Data/2_3_random.txt')

nSeconds3 = 10
col3=len(Y3[0])
row3=len(Y3)//col3
fps3=row3//nSeconds3

fig3 = plt.figure( figsize=(8,8) )

im3 = plt.imshow(Y3[0:col3,0:col3])

def animate_func3(i):
    im3.set_array(Y3[i*col3:(i+1)*col3,0:col3].copy())
    return [im3]

anim3 = animation.FuncAnimation(fig3, animate_func3, frames = row3, interval = 1000 /fps3 )

anim3.save('Animations/Animation_3_random.mp4', fps=fps3, extra_args=['-vcodec', 'libx264'])
