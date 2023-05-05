import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, fftfreq
import matplotlib.animation as animation

#constants and parameters
N=64
n=np.arange(-(N-1)/2, (N-1)/2 + 1)

def ind(i):
    return np.where(n==i)

def force(x, alpha, beta):
    return x + alpha*np.power(x, 2)

T=60000
dt=0.1
#FPU parameters
alpha = 0.25
beta = 0.0
one_period_wave = 2*np.pi/(N-1)
ddu=np.zeros(N)
du=np.zeros(N)
u=np.zeros(N)
omega = np.zeros(N)
energy = np.zeros(N)
a=np.zeros(N)
da=np.zeros(N)
#modes and time
time = []
first_mode = []
second_mode = []
third_mode = []
fourth_mode = []
fifth_mode = []
sixth_mode = []
#initial conditions
#set up initial conditions
for index in np.arange(0, N):
    u[index]= 0.0*np.sin(0.5*one_period_wave*index) + 0.5*np.sin(1*one_period_wave*index)
#integration loop
plot_counter = 0
#fig, (ax1, ax2) = plt.subplots(1, 2)
fig, ax = plt.subplots()
for t in np.arange(0, T, dt):
    #plot configuration and modes
    if t % 500== 0:
        #ax1.cla()
        #ax2.cla()
        ax.cla()
        #normal coordinates
        for k in np.arange(1, N+1):
            a[k-1] = np.dot(u, np.sin(k*np.pi/N*np.arange(1, N+1)))    
            da[k-1] = np.dot(du, np.sin(k*np.pi/N*np.arange(1, N+1)))    
            energy[k-1] = 0.5*np.power(da[k-1], 2) + 2*np.power(np.sin(np.pi*k/2/N)*a[k-1], 2)
        nf = fftfreq(N, 1)[:N//2]
        #ax1.plot(n, u)
        #ax2.plot(energy)
        #ax1.set_ylim([-0.1,0.1])
        #ax2.set_ylim([0,0.01])
        time.append(t)
        first_mode.append(energy[0])
        second_mode.append(energy[1])
        third_mode.append(energy[2])
        fourth_mode.append(energy[3])
        fifth_mode.append(energy[4])
        sixth_mode.append(energy[5])
        ax.plot(time, first_mode, label='1st')
        ax.plot(time, second_mode, label='2nd')
        ax.plot(time, third_mode, label='3rd')
        ax.plot(time, fourth_mode, label='4th')
        ax.plot(time, fifth_mode, label='5th')
        ax.plot(time, sixth_mode, label='6th')
        ax.legend()
    plt.pause(1e-5)
    #calculate forces
    #fixed boundary
    u[0]=0
    u[N-1]=0
    #periodic boundary
    #ddu[0]=force(u[1] - u[0], alpha, beta) + force(u[N-1] - u[0], alpha, beta)
    #ddu[N-1]=force(u[N-2] - u[N-1], alpha, beta) + force(u[0] - u[N-1], alpha, beta)
    #inner particles
    for i in np.arange(1, N-1):
        ddu[i]=force(u[i+1]-u[i], alpha, beta)-force(u[i]-u[i-1], alpha, beta)
    #integration
    for i in np.arange(0, N):
        du[i] = du[i]+ddu[i]*dt
        u[i] = u[i]+du[i]*dt
    #set to zero forces
    ddu=np.zeros(N)
    #increment step counter
    plot_counter+=1
plt.show()
