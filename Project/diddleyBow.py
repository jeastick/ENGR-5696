import numpy as np
import math as m
import matplotlib
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import scipy.io.wavfile

#### GLOBAL VARIABLES ####
T = 444           # Tension of material (Newtons)
mu = 8050          # Density of material (kg/m3)
yp = 0.01          # Height of "pluck" from neutral axis (m)
xp = 0.75           # Location of "pluck" along string (m)
v = (T*mu)**(0.5)
L = 1              # Period (string length) (m)
l = L/2            # Half Period (m)

pup = 0.75         # sound pickup location (m)

x_res = 501     # resolution of x

f_res = 44100        # time sampling frequency (Hz)
t_res = 1/f_res     # time step (s)

sampletime = 0.1 #total length of sample to take (s)

t_steps = int(sampletime*f_res)         # number of time samples

def y0(x):
    if(x<=xp):
        return yp/xp*x
    if(x>xp):  
        return (-yp/(L-xp))*x+yp+yp/(L-xp)*xp

def sin_term(n,x):
    return m.sin(n*m.pi*x/l)

# We will calculate Fourier sin series coefficients b_n by 


x_range = np.linspace(0,L,num=x_res)

y_init = np.zeros(len(x_range))

for i in range(len(x_range)):
    y_init[i] = y0(x_range[i])


integrand = np.zeros(len(x_range))

coefs = 5
bn = np.zeros(coefs) #i'll calculate this many bn terms

fig, ax = plt.subplots()
ax.plot(x_range,y_init)

for n in range(coefs):
    print("Looping with n = " + str(n))
    for i in range(len(x_range)): #generate the function i will take the integral of to find bn
        integrand[i] = y0(x_range[i])*m.sin(n*m.pi*x_range[i]/L)
    bn[n] = 1/l*np.trapz(integrand,x_range)
    print("b" + str(n) + " is " + str(bn[n]))
    # ax.plot(x_range,integrand)



y00 = np.zeros(len(x_range))

for i in range(len(x_range)):
    for n in range(len(bn)):
        y00[i] = y00[i] + bn[n]*m.sin((n)*m.pi*x_range[i]/L)
# ax.plot(x_range,y00)



t_range = np.arange(0,t_res*t_steps,t_res)
soundwave = np.zeros(t_steps,dtype = np.float32)

print("t_range is: " + str(t_range))

for t_step in range(len(t_range)):
    yxt = np.zeros(len(x_range))
    for i in range(len(x_range)):
        for n in range(len(bn)):
            yxt[i] = yxt[i] + bn[n]*m.sin((n)*m.pi*x_range[i]/L)*m.cos((n)*m.pi*v*t_range[t_step]/L)
        # print("x_range i is" + str(x_range[i]))
        if(x_range[i] == pup):
            soundwave[t_step] = yxt[i]
    ax.plot(x_range,yxt)
    
plt.show()

print("Datatype for soundfile number is " + str(type(soundwave[2])))

scipy.io.wavfile.write("wave1.wav", f_res, soundwave)

print("soundwave is " + str(soundwave))
fig2, ax = plt.subplots()
ax.plot(t_range,soundwave)
plt.show()









