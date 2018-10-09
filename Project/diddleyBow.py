import numpy as np
import math as m
import matplotlib
import matplotlib.pyplot as plt
import scipy.integrate as integrate


T = 1000
mu = 8700
yp = 1
xp = 8
L = 10

x = np.linspace(0,L,num=1000)
print(x)

fx = np.piecewise(x, [ x <= xp, x > xp], [lambda x: yp/xp*x, lambda x: (-yp/(L-xp))*x+yp+yp/(L-xp)*xp])
print(fx)


terms = 10
bn_coefs = np.zeros(terms) 
l = L/2

#fucked still:
# for n in range(terms):
#     b = 1/l*integrate.quad(lambda n: fx*m.sin(n*m.pi*x/l), 0, 2*l)
#     print(b)




fig, ax = plt.subplots()
ax.plot(x,fx)
plt.show()


