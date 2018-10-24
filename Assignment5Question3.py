# ENGR 5696
# JEFF EASTICK
# 2018/9/25
# ASSIGNMENT 4 - QUESTION 1

import math as m
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

r = 10
e = 1      #y
h = e*r      #z
res = 100
coefs = 10
mu=2
G=1

pi = np.pi

print("e = " + str(e))
print("h = " + str(h))



def series_term(N,Y,Z):
    a = 0
    for n in range(N):
        n+=1
        a += (-1)**n*(32/(((2*n-1)**3)*pi**3)*(np.cosh((2*n-1)*(pi*Z/e)))/(np.cosh((2*n-1)*(pi/2)*(h/e)))*np.cos((2*n-1)*pi*(Y/e)))
    return a


y = np.linspace(-e/2, e/2, res)
z = np.linspace(-h/2, h/2, res)

print("y = " + str(y))
print("z = " + str(z))

uy, uz = np.meshgrid(y,z)

u_mag = ((G*e**2)/(8*mu))

ustar = (1-(2*uy/e)**2)
ustar2 = series_term(coefs,uy,uz)
print("U star 2 = " + str(ustar2))
u = u_mag*(ustar + ustar2)

Q =  np.trapz(np.trapz(u,y),z)
Q_ = Q/((G*e**4)/(8*mu))
Rhe = (h/e)/(1+(h/e))
Q_rh = (Rhe)**4

print("Series term is " + str(series_term(19,0,0)))

print("Calculated Flow Rate Q  = " + str(Q))
print("Normalized Flow rate Q_ = " + str(Q_))
# print("(Rhe)^4 = " + str(Rhe**4))
print("Flow Rate Qrh = " + str(Q_rh))

# levels = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1
levels = 15
fig, ax = plt.subplots()
plot = plt.contour(y/e,z/h,u/u_mag,levels)
plt.clabel(plot, inline=1, fontsize=10)

plt.xlim(-0.5, 0.5)
plt.ylim(-0.5, 0.5)

plt.xlabel('y/e')
plt.ylabel('z/h')
plt.show()

plot2 = plt.contour(y/e,z/h,ustar,levels)
plt.clabel(plot2, inline=1, fontsize=10)
plt.xlabel('y/e')
plt.ylabel('z/h')
plt.show()

plot3 = plt.contour(y/e,z/h,ustar2,levels)
plt.clabel(plot3, inline=1, fontsize=10)
plt.xlabel('y/e')
plt.ylabel('z/h')
plt.show()
# plot = plt.contour(y,z,j,15)



