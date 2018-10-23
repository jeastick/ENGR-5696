# ENGR 5696
# JEFF EASTICK
# 2018/9/25
# ASSIGNMENT 4 - QUESTION 1

import math as m
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

r = 0.75
e = 1      #y
h = e*r      #z
res = 100
coefs = 100

pi = np.pi

print("e = " + str(e))
print("h = " + str(h))



def series_term(N,Y,Z):
    a = 0
    for n in range(N):
        n+=1
        a += (-1)**n*(32/(((2*n-1)**3)*pi**3)*(np.cosh((2*n-1)*(pi*Z/e)))/(np.cosh((2*n-1)*(pi/2)*(h/e)))*np.cos((2*n-1)*pi*(Y/e)))
    return a


y = np.arange(-e/2, e/2, 1/res)
z = np.arange(-h/2, h/2, 1/res)

print("y = " + str(y))
print("z = " + str(z))

uy, uz = np.meshgrid(y,z)


u = (1-(2*uy/e)**2 + series_term(coefs,uy,uz))
utest = series_term(coefs,uy,uz)


fig, ax = plt.subplots()
plot = plt.contour(y/e,z/h,u,10)
plt.clabel(plot, inline=1, fontsize=10)

plt.xlim(-0.5, 0.5)
plt.ylim(-0.5, 0.5)

plt.xlabel('y/e')
plt.ylabel('z/h')
plt.show()



