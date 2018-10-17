# ======================================
# =        ENGR-5696 PROJECT           =                 
# =                                    = 
# =         JEFF EASTICK               =             
# =         OCTOBER 2018               =             
# ======================================

import numpy as np
import math as m
import matplotlib
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import scipy.io.wavfile

from mpl_toolkits.mplot3d import axes3d

class bow:

    def __init__(self, name, Tension, LinDens, PluckHeight, PluckLocation, StringLength, PickupLocation, SampleTime, Coefficients):
        
        self.name = name
        self.T = Tension                            # Tension of material (Newtons)
        self.mu = LinDens                           # Linear Density of material (kg/m)
        self.yp = PluckHeight                       # Height of "pluck" from neutral axis (m)
        self.L = StringLength                       # Period (string length) (m)
        self.xp = PluckLocation/100*self.L          # Location of "pluck" along string (m)
        self.coefs = Coefficients

        self.v = (self.T/self.mu)**(0.5)
        self.l = self.L/2                           # Half Period (m)

        self.x_res = 500                            # Resolution of x
        self.pup = PickupLocation/100*self.L        # Sound pickup location (m)


        self.f_res = 5000                           # Time sampling frequency (Hz)
        self.t_res = 1/self.f_res                   # Time step (s)
    
        self.sampletime = SampleTime                # Total length of sample to take (s)

        self.t_steps = int(self.sampletime*self.f_res)         # Number of time samples

        self.generateBow()

    def y0(self,x):
        if(x<=self.xp):
            return self.yp/self.xp*x
        if(x>self.xp):  
            return (-self.yp/(self.L-self.xp))*x+self.yp+self.yp/(self.L-self.xp)*self.xp

    def sin_term(self,n,x):
            return m.sin(n*m.pi*x/self.l)

    def generateBow(self):

        self.x_range = np.linspace(0,self.L,num=self.x_res+1)

        self.y_init = np.zeros(len(self.x_range))

        for i in range(len(self.x_range)):
            self.y_init[i] = self.y0(self.x_range[i])


        self.integrand = np.zeros(len(self.x_range))

        self.bn = np.zeros(self.coefs) 

        for n in range(self.coefs):
            print("Calculating " + str(self.coefs) + " fourier coefficients with n = " + str(n))
            for i in range(len(self.x_range)):
                self.integrand[i] = self.y0(self.x_range[i])*m.sin(n*m.pi*self.x_range[i]/self.L)
            self.bn[n] = 1/self.l*np.trapz(self.integrand,self.x_range)
            print("b" + str(n) + " is " + str(self.bn[n]))



        self.y00 = np.zeros(len(self.x_range))

        for i in range(len(self.x_range)):
            for n in range(len(self.bn)):
                self.y00[i] = self.y00[i] + self.bn[n]*m.sin((n)*m.pi*self.x_range[i]/self.L)

        self.t_range = np.arange(0,self.t_res*self.t_steps,self.t_res)
        self.soundwave = np.zeros(self.t_steps,dtype = np.float32)

        # print("self.t_range is: " + str(self.t_range))
   
        self.yxt = np.zeros((len(self.t_range),len(self.x_range)))

        for t_step in range(self.yxt.shape[0]):
            yx = np.zeros(self.yxt.shape[1])
            for i in range(self.yxt.shape[1]):
                for n in range(len(self.bn)):
                    yx[i] = yx[i] + self.bn[n]*m.sin((n)*m.pi*self.x_range[i]/self.L)*m.cos((n)*m.pi*self.v*self.t_range[t_step]/self.L)
                if(self.x_range[i] == self.pup):
                    self.soundwave[t_step] = yx[i]*100
                self.yxt[t_step,i] = yx[i]

        scipy.io.wavfile.write(str(self.name) + ".wav", self.f_res, self.soundwave)

## x_range and y_init
    def plot_y_init(self):
        fig, ax = plt.subplots()
        ax.plot(self.x_range,self.y_init)
        plt.title(str(self.name + ": exact value of y_init (boundary case)"))
        plt.show()

## x_range and self.integrand
    def plot_integrand(self):
        fig, ax = plt.subplots()
        ax.plot(self.x_range,self.integrand)
        plt.title(str(self.name + ": Fourier integrand function for " + self.name))
        plt.show()

## x_range and y00
    def plot_y00(self):
        fig, ax = plt.subplots()
        ax.plot(self.x_range,self.y00)
        plt.title(str(self.name + ": Fourier series approximation of y_init"))        
        plt.show()

    def get_bn(self):
        return self.bn

# ## x_range and yxt over i time steps, skipping j time steps
    def plot_yxt(self,numSteps,skipSteps):
        fig, ax = plt.subplots()
        counter = 0
        for t_step in range(self.yxt.shape[0]):
            if(skipSteps != 0 and counter < numSteps and t_step%skipSteps == 0):
                ax.plot(self.x_range,self.yxt[t_step,])
                counter += 1
            elif(counter < numSteps):
                ax.plot(self.x_range,self.yxt[t_step,])
                counter += 1   
        plt.title(str(self.name + ": First and every other " + str(skipSteps) + " time steps of y(x,t) up to " + str(numSteps) + " steps"))        
        plt.show()

## t_range and soundwave
    def plot_soundwave(self):
        fig, ax = plt.subplots()
        ax.plot(self.t_range,self.soundwave)
        plt.title(str(self.name + ": Soundwave of string vibration with sound pick-up at x = " + str(self.pup)))        
        plt.show()

#MAIN

# We will create a few digital guitar and bass strings:

t_global    = 0.1
c_global    = 10
pluck_height_global = 0.01 # m

L_guitar    = 0.640 # m
L_bass      = 0.860 # m
mu_guitar   = 0.001 # kg/m
mu_bass     = 0.016 # kg/m




#               Name, Tension (N), LinDens (kg/m), PluckHeight (m), PluckLocation (% of L), StringLength (m) , PickupLocation (% of L), SampleTime (s), Coefficients (number):

# Control Guitar
Guitar1 = bow("Guitar1", 45, mu_guitar, pluck_height_global, 15, L_guitar, 20, t_global, c_global)
Guitar1.plot_y_init()
Guitar1.plot_integrand()
Guitar1.plot_y00()
Guitar1.plot_soundwave()
Guitar1.plot_yxt(10,0)

# Change the pluck location for Guitar2:
Guitar2 = bow("Guitar2", 45, mu_guitar, pluck_height_global, 5, L_guitar, 20, t_global, c_global)
Guitar2.plot_y_init()
Guitar2.plot_integrand()
Guitar2.plot_y00()
Guitar2.plot_soundwave()
Guitar2.plot_yxt(10,0)

# Change the pick-up location for Guitar3:
Guitar3 = bow("Guitar3", 45, mu_guitar, pluck_height_global, 15, L_guitar, 40, t_global, c_global)
Guitar3.plot_y_init()
Guitar3.plot_integrand()
Guitar3.plot_y00()
Guitar3.plot_soundwave()
Guitar3.plot_yxt(10,0)

# Control Bass
Bass1 = bow("Bass1", 100, mu_bass, pluck_height_global, 15, L_bass, 20, t_global, c_global)
Bass1.plot_y_init()
Bass1.plot_integrand()
Bass1.plot_y00()
Bass1.plot_soundwave()
Bass1.plot_yxt(10,0)

# Change the tension for Bass2:
Bass2 = bow("Bass2", 130, mu_bass, pluck_height_global, 15, L_bass, 20, t_global, c_global)
Bass2.plot_y_init()
Bass2.plot_integrand()
Bass2.plot_y00()
Bass2.plot_soundwave()
Bass2.plot_yxt(10,0)

# Change the string density for Bass3 by + 15%:
Bass3 = bow("Bass3", 100, mu_bass*1.15, pluck_height_global, 15, L_bass, 20, t_global, c_global)
Bass3.plot_y_init()
Bass3.plot_integrand()
Bass3.plot_y00()
Bass3.plot_soundwave()
Bass3.plot_yxt(10,0)