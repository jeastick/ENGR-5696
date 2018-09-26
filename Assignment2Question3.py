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
    Yn = (H**(2*n)*lam**(2*n)*a0*y**(2*n))/((4**n)*m.factorial(2*n))
    print("made it through Y_y_series")
    return Yn

def Y_y_numerical(iterations,a0,lam, H, y):
    Y_numerical = 0
    for i in range(iterations):
        Yn = Y_y_series(i,ao,Lambda,H,y_bc)
        Y_numerical = Y_numerical + Yn
    return Y_numerical

def Y_y_exact(lam,H,y):
    a = m.exp(lam*H)
    b = m.exp(-lam*H)
    Yy = (1-(1-a)/(b-a))*a + ((1-a)/(b-a))*b
    print("made it through Y_y_exact")
    return Yy


L=1
k = 2*m.pi/L
w = 1
To = 1
H = 10
y_bc = 1

Lambda = 1

testyyexact = Y_y_exact(Lambda,H,y_bc)
print ("Y_y_exact: " + str(testyyexact))

# First we solve for a_0 using the boundary condition Y(1) = 1, using n iterations




ao = 1 # initial guess for ao
Y1 = 1 
tolerance = 0.00000001
converged = False
N = 0
Nmax = 100
n = 10

while converged == False and N < Nmax: 

    Y_numerical = Y_y_numerical(n,ao,Lambda,H,y_bc) 

    error = (Y_numerical - y_bc)

    if abs(error) < tolerance:
        print ("Solution for a_o has converged!: Y_(1) = " +str(Y_numerical) + " when ao = " + str(ao))
        converged = True
    else:
        ao = ao - 0.01*error
        print ("new ao is " + str(ao))
        N+=1




# Then we will have both a numerical and analytical solution for Y(y) which we can plot and compare



# y = np.arange(0.0,10.0,0.01)

# # numerical_solution = Y_y_numerical(n,ao,Lambda,H,y)
# # analytical_solution = Y_y_exact(Lambda,H,y)

# fig, ax = plt.subplots()
# ax.plot(y, Y_y_numerical(n,ao,Lambda,H,y))
# ax.plot(y, Y_y_exact(Lambda,H,y))

# ax.set(xlabel='y', ylabel='Y(y)',
#        title='Comparison')
# ax.grid()

# fig.savefig("test.png")
# plt.show()
