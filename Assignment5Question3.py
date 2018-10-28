# ENGR 5696
# JEFF EASTICK
# 2018/9/25
# ASSIGNMENT 4 - QUESTION 1

import math as m
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate

class duct:

    def __init__(self, name, R, E, MU, g):
        self.name = name
        self.r = R
        self.e = E      #y
        self.h = self.e*self.r      #z
        self.res = 100
        self.coefs = 10
        self.mu = MU
        self.G = g
        self.u_mag = ((self.G*self.e**2)/(8*self.mu))
        pi = np.pi


        def u_y_z(Z,Y,N):
            a = 0
            for n in range(N):
                n+=1
                a += (-1)**n*(32/(((2*n-1)**3)*pi**3)*(np.cosh((2*n-1)*(pi*Z/self.e)))/(np.cosh((2*n-1)*(pi/2)*(self.h/self.e)))*np.cos((2*n-1)*pi*(Y/self.e)))

            return self.u_mag*((1-(2*Y/self.e)**2) + a)
    
    
        def series_term(N,Y,Z):
            a = 0
            for n in range(N):
                n+=1
                a += (-1)**n*(32/(((2*n-1)**3)*pi**3)*(np.cosh((2*n-1)*(pi*Z/self.e)))/(np.cosh((2*n-1)*(pi/2)*(self.h/self.e)))*np.cos((2*n-1)*pi*(Y/self.e)))
            return a
    
    
        self.y = np.linspace(-(self.e/2), self.e/2, self.res)
        self.z = np.linspace(-(self.h/2), self.h/2, self.res)
    
        # print("y = " + str(self.y))
        # print("z = " + str(self.z))
    
        self.uy, self.uz = np.meshgrid(self.y,self.z)
    

    
        # self.ustar = (1-(2*self.uy/self.e)**2)
        # self.ustar2 = series_term(self.coefs,self.uy,self.uz)
        # self.u = self.u_mag*(self.ustar + self.ustar2)

        self.u_test = u_y_z(self.uz,self.uy,self.coefs)
        self.Q_test = integrate.dblquad(u_y_z,-self.e/2, self.e/2, -self.h/2, self.h/2, args = (self.coefs,))
        self.Q_test = self.Q_test[0]
        print("Q_test is: " + str(self.Q_test))

        # self.Q =  np.trapz(integrate.cumtrapz(self.u,self.y,initial=0),self.z)
        # self.Q_ = self.Q_test/((self.G*self.e**4)/(8*self.mu))
        self.Q_ = self.Q_test/((self.G*self.e**3*self.h)/(8*self.mu))
        self.Rhe = (self.h/self.e)/(1+(self.h/self.e))
        self.Q_Rh = (self.Rhe)**4
        print("For " + self.name + " flow rates are as follows:")
        print("Calculated Flow Rate Q  = " + str(self.Q_test))
        print("Normalized Flow rate Q_ = " + str(self.Q_))
        print("Flow Rate Q_Rh = " + str(self.Q_Rh))
        print("Where R = " + str(self.r))
        print("Where e = " + str(self.e))
        print("Where h = " + str(self.h))
        print("")
    
        self.levels = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]
        fig, ax = plt.subplots()
        plot = plt.contour(self.y/self.e,self.z/self.h,self.u_test/self.u_mag,self.levels)
        plt.clabel(plot, inline=1, fontsize=10)
        ax.grid()
        fig.set_figheight(10)
        fig.set_figwidth(10)
        plt.xlim(-0.5, 0.5)
        plt.ylim(-0.5, 0.5)
        plt.title('Normalized flow contours for ' + str(self.name))
        plt.xlabel('y/e')
        plt.ylabel('z/h')
        plt.savefig("A5Q3_" + str(self.name + ".png"), bbox_inches='tight')
        plt.show()

    def getQ_exact(self):
        return self.Q_test

    def getQ_normal(self):
        return self.Q_

    def getQ_Rh(self):
        return self.Q_Rh

    def getR(self):
        return self.r

    def name(self):
        return self.name


# def __init(self, name, R, E, MU, g)
Duct1 = duct("r = 1",  1,  1, 1, 1)
Duct2 = duct("r = 3",  3,  1, 1, 1)
Duct3 = duct("r = 5",  5,  1, 1, 1)
Duct4 = duct("r = 9",  9,  1, 1, 1)
Duct5 = duct("r = 11", 11, 1, 1, 1)

Q_normals = np.array([Duct1.getQ_normal(), Duct2.getQ_normal(), Duct3.getQ_normal(), Duct4.getQ_normal(), Duct5.getQ_normal()])
Q_exacts = np.array([Duct1.getQ_exact(), Duct2.getQ_exact(), Duct3.getQ_exact(), Duct4.getQ_exact(), Duct5.getQ_exact()])

Q_Rhes = np.array([Duct1.getQ_Rh(), Duct2.getQ_Rh(), Duct3.getQ_Rh(), Duct4.getQ_Rh(), Duct5.getQ_Rh()])

ratios = np.array([Duct1.getR(), Duct2.getR(), Duct3.getR(), Duct4.getR(), Duct5.getR()])

print(Q_normals)
print(Q_Rhes)
print(ratios)


fig, ax = plt.subplots()
ax.plot(ratios,Q_normals)
ax.plot(ratios,Q_Rhes)
ax.grid()
fig.set_figheight(10)
fig.set_figwidth(10)
plt.xlabel('r = h/e')
plt.ylabel('Flow Rate')
plt.title('Normalized Flow Rate vs. Cross Sectional Aspect Ratio')
ax.legend(('Q_normal', 'Q_Rhe'),loc = 'right')
plt.show()







