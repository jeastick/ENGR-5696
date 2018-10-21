# ENGR 5696
# JEFF EASTICK
# 2018/9/25
# ASSIGNMENT 4 - QUESTION 1

import math as m
import matplotlib
import matplotlib.pyplot as plt
import numpy as np


l = 1
h = 0.1
M = 101
N_1 = 1
N_5 = 5
N_7 = 7
N_9 = 9
N_11 = 11

pi = m.pi

def function(x):
    if (x >= 0 and x < l/2):
        return 2*h*x/l
    if (x >= l/2 and x <= l):
        return -2*h*x/l+2*h

def bn(n):
    bn = 8*h/(n**2*pi**2)*m.sin(n*pi/2)
    print("When n = " + str(n) + "bn = " + str(bn))
    return bn

def function_fourier(x,N):
    temp = 0

    for n in range(N):
        n+=1
        print(str("N is " + str(N) + " and n is " + str(n)))
        temp += bn(n)*m.sin(n*pi*x/l)
    return temp
    



x_range = np.linspace(0,l,num=M)
fx_exact = np.zeros(M)

fx_fourier_1  = np.zeros(M)
fx_fourier_5  = np.zeros(M)
fx_fourier_7  = np.zeros(M)
fx_fourier_9  = np.zeros(M)
fx_fourier_11 = np.zeros(M)

error_numerator_1  = 0
error_numerator_5  = 0
error_numerator_7  = 0
error_numerator_9  = 0
error_numerator_11 = 0
error_denominator  = 0


for i in range(len(x_range)):
    fx_exact[i] = function(x_range[i])
    fx_fourier_1[i]  = function_fourier(x_range[i],N_1)
    fx_fourier_5[i]  = function_fourier(x_range[i],N_5)
    fx_fourier_7[i]  = function_fourier(x_range[i],N_7)
    fx_fourier_9[i]  = function_fourier(x_range[i],N_9)
    fx_fourier_11[i] = function_fourier(x_range[i],N_11)
    
    error_numerator_1  += (fx_fourier_1[i]  - fx_exact[i])**2
    error_numerator_5  += (fx_fourier_5[i]  - fx_exact[i])**2
    error_numerator_7  += (fx_fourier_7[i]  - fx_exact[i])**2
    error_numerator_9  += (fx_fourier_9[i]  - fx_exact[i])**2
    error_numerator_11 += (fx_fourier_11[i] - fx_exact[i])**2

    error_denominator += (fx_exact[i])**2


E_1 = error_numerator_1 /error_denominator
E_5 = error_numerator_5 /error_denominator
E_7 = error_numerator_7 /error_denominator
E_9 = error_numerator_9 /error_denominator
E_11 = error_numerator_11/error_denominator

error_trend = np.array([E_1,E_5,E_7,E_9,E_11])
error_N_vals = np.array([N_1,N_5,N_7,N_9,N_11])
print(error_trend)
print(error_N_vals)
print(E_1 )
print(E_5 )
print(E_7 )
print(E_9 )
print(E_11)


fig, ax = plt.subplots()
fig.set_figheight(8.5)
fig.set_figwidth(11)
ax.plot(x_range,fx_exact)
ax.plot(x_range,fx_fourier_1)
ax.plot(x_range,fx_fourier_5)
ax.plot(x_range,fx_fourier_7)
ax.plot(x_range,fx_fourier_9)
ax.plot(x_range,fx_fourier_11)
ax.legend(('f(x)_exact','N=1','N=5,','N=7', 'N=9', 'N=11'),loc = 'right')

title = 'Fourier Series Expansion of f(x)'
ax.set(xlabel='x', ylabel='f(x)',
       title=title)
ax.grid()

fig.savefig("A4Q1_Fourier_Plot.png")


fig,ax = plt.subplots()
fig.set_figheight(8.5)
fig.set_figwidth(11)
ax.plot(error_N_vals,error_trend)

title2 = 'Relative error of Fourier Expansion vs. N'
ax.set(xlabel='N', ylabel='Relative Error',
       title=title2)
ax.grid()
fig.savefig("A4Q1_Error_Plot.png")



plt.show()

