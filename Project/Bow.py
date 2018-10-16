import numpy as np
import math as m
import matplotlib
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import scipy.io.wavfile

class bow:

    def __init__(self,Tension, LinDens, PluckHeight, PluckLocation, StringLength, PickupLocation, SampleTime, Coefficients):
        #### GLOBAL VARIABLES ####
        self.T = Tension           # Tension of material (Newtons)
        self.mu = LinDens          # Linear Density of material (kg/m)
        self.yp = PluckHeight          # Height of "pluck" from neutral axis (m)
        self.xp = PluckLocation           # Location of "pluck" along string (m)
        self.L = StringLength              # Period (string length) (m)
        self.coefs = Coefficients

        self.v = (self.T/self.mu)**(0.5)
        self.l = self.L/2            # Half Period (m)

        self.pup = PickupLocation         # sound pickup location (m)

        self.x_res = 501     # resolution of x

        self.f_res = 2000        # time sampling frequency (Hz)
        self.t_res = 1/self.f_res     # time step (s)

        self.sampletime = SampleTime #total length of sample to take (s)

        self.t_steps = int(self.sampletime*self.f_res)         # number of time samples

    def y0(self,x):
                if(x<=self.xp):
                    return self.yp/self.xp*x
                if(x>self.xp):  
                    return (-self.yp/(self.L-self.xp))*x+self.yp+self.yp/(self.L-self.xp)*self.xp

    def sin_term(self,n,x):
            return m.sin(n*m.pi*x/self.l)


    def generateBow(self):

        # We will calculate Fourier sin series coefficients b_n by 

        self.x_range = np.linspace(0,self.L,num=self.x_res)

        self.y_init = np.zeros(len(self.x_range))

        for i in range(len(self.x_range)):
            self.y_init[i] = self.y0(self.x_range[i])


        integrand = np.zeros(len(self.x_range))

        self.bn = np.zeros(self.coefs) #i'll calculate this many self.bn terms

        fig, ax = plt.subplots()
        ax.plot(self.x_range,self.y_init)

        for n in range(self.coefs):
            print("Looping with n = " + str(n))
            for i in range(len(self.x_range)): #generate the function i will take the integral of to find self.bn
                integrand[i] = self.y0(self.x_range[i])*m.sin(n*m.pi*self.x_range[i]/self.L)
            self.bn[n] = 1/self.l*np.trapz(integrand,self.x_range)
            print("b" + str(n) + " is " + str(self.bn[n]))
            # ax.plot(self.x_range,integrand)



        self.y00 = np.zeros(len(self.x_range))

        for i in range(len(self.x_range)):
            for n in range(len(self.bn)):
                self.y00[i] = self.y00[i] + self.bn[n]*m.sin((n)*m.pi*self.x_range[i]/self.L)
        ax.plot(self.x_range,self.y00)



        self.t_range = np.arange(0,self.t_res*self.t_steps,self.t_res)
        self.soundwave = np.zeros(self.t_steps,dtype = np.float32)

        print("self.t_range is: " + str(self.t_range))

        for t_step in range(len(self.t_range)):
            yxt = np.zeros(len(self.x_range))
            for i in range(len(self.x_range)):
                for n in range(len(self.bn)):
                    yxt[i] = yxt[i] + self.bn[n]*m.sin((n)*m.pi*self.x_range[i]/self.L)*m.cos((n)*m.pi*self.v*self.t_range[t_step]/self.L)
                # print("self.x_range i is" + str(self.x_range[i]))
                if(self.x_range[i] == self.pup):
                    self.soundwave[t_step] = yxt[i]*100
            if(t_step < 300 and t_step%2==0):
                ax.plot(self.x_range,yxt)
        
        plt.show()

        # print("Datatype for soundfile number is " + str(type(self.soundwave[2])))

        scipy.io.wavfile.write("wave1.wav", self.f_res, self.soundwave)

        print("self.soundwave is " + str(self.soundwave))
        fig2, ax = plt.subplots()
        ax.plot(self.t_range,self.soundwave)
        plt.show()


            # Tension, LinDens, PluckHeight, PluckLocation, StringLength, PickupLocation, SampleTime, Coefficients):
bow1 = bow(80, 0.020, 0.01, 0.20, 1, 0.5, 0.1, 10)
bow1.generateBow()




