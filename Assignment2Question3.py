# ENGR 5696
# JEFF EASTICK
# 2018/9/25
# ASSIGNMENT 2 - QUESTION 3

import math as m
import matplotlib
import matplotlib.pyplot as plt
import numpy as np




def A_x_t(T_,k,x,w,t):
    Axt = T_*m.exp(1j*(k*x+w*t))
    return Axt

def Y_y_series(n, a0, lam, H, y):
    y_star = y-H/2
    y_sharp = 2*y_star/H
    Yn = (H**(2*n)*lam**(2*n)*a0*y_sharp**(2*n))/((4**n)*m.factorial(2*n))
    return Yn

def Y_y_numerical(iterations,a0,lam, H, y):
    Y_numerical = 0
    for i in range(iterations):
        Yn = Y_y_series(i,ao,Lambda,H,y)
        Y_numerical = Y_numerical + Yn
    return Y_numerical

def Y_y_exact(lam,H,y):
    a = m.exp(lam*H)
    b = m.exp(-lam*H)
    Yy = (1-(1-a)/(b-a))*m.exp(lam*y) + ((1-a)/(b-a))*m.exp(-lam*y)
    return Yy


L=10
k = 2*m.pi/L
w = 1
To = 1
H = 5
y_bc = 1
Lambda = 1


# First we solve for a_0 using the boundary condition Y(0)=Y(H)=1, using n iterations

ao = 1 # initial guess for ao

tolerance = 0.00000001
converged = False
N = 0
Nmax = 1000
n = 10

while converged == False and N < Nmax: 

    Y_numerical = Y_y_numerical(n,ao,Lambda,H,H) 

    error = (Y_numerical - y_bc)

    if abs(error) < tolerance:
        print ("Solution for a_o has converged!: Y_(H) = " +str(Y_numerical) + " when ao = " + str(ao))
        converged = True
    else:
        ao = ao - 0.01*error
        # print ("new ao is " + str(ao)) #tracer
        N+=1




# Then we will have both a numerical and analytical solution for Y(y) which we can plot and compare

resolution = 0.1
y = np.arange(0.0,H+resolution,resolution)

numerical_solution = np.zeros(len(y)) 
analytical_solution =np.zeros(len(y))

for i in range(len(y)):
    numerical_solution[i]=Y_y_numerical(n,ao,Lambda,H,y[i])
    analytical_solution[i]=Y_y_exact(Lambda,H,y[i]) 
    #print("i = " + str(i) + " y[i] = " + str(y[i]) + " Num = " + str(numerical_solution[i]) + "     Ana = " + str(analytical_solution[i])) #tracer


fig, ax = plt.subplots()
ax.plot(numerical_solution,y)

numtitle = 'Numerical Solution of Y(y) with H = ' + str(H)
title = 'Y(y) with H = ' + str(H)
ax.set(xlabel='Y(y)_numerical', ylabel='y',
       title=numtitle)
ax.grid()
fig.savefig("A2Q3_Numerical_Solution.png")



fig,ax = plt.subplots()
ax.plot(analytical_solution,y)

anatitle = 'Analytical Solution of Y(y) with H = ' + str(H)
ax.set(xlabel='Y(y)_analytical', ylabel='y',
       title=anatitle)
ax.grid()

fig.savefig("A2Q3_Analytical_Solution.png")
plt.show()

