import numpy as np
import math as m
import matplotlib
import matplotlib.pyplot as plt
import scipy.integrate as integrate


#### GLOBAL VARIABLES ####
T = 444           # Tension of material (Newtons)
mu = 8050          # Density of material (kg/m3)
yp = 0.005          # Height of "pluck" from neutral axis (m)
xp = 4           # Location of "pluck" along string (m)
v = (T*mu)**(0.5)
L = 10              # Period (string length) (m)
l = L/2            # Half Period (m)

x_res = 500 # resolution of x


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



t_range = np.linspace(0,1,num=10)



for t_step in range(len(t_range)):
    yxt = np.zeros(len(x_range))
    for i in range(len(x_range)):
        for n in range(len(bn)):
            yxt[i] = yxt[i] + bn[n]*m.sin((n)*m.pi*x_range[i]/L)*m.cos((n)*m.pi*v*t_step/L)
    ax.plot(x_range,yxt)
    




# print("bn array is " + str(bn))

print(y_init)
plt.show()









