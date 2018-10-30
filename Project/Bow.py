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
import gc



from mpl_toolkits.mplot3d import axes3d

class bow:

    def __init__(self, name, Tension, LinDens, PluckHeight, PluckLocation, StringLength, PickupLocation, Fret, SampleTime, Coefficients):
        
        self.name = name
        print(self.name + ": Initializing bow object.")
        self.T = Tension                            # Tension of material (Newtons)
        self.mu = LinDens                           # Linear Density of material (kg/m)
        self.yp = PluckHeight                       # Height of "pluck" from neutral axis (m)
        self.fret = Fret
        self.L = StringLength*(1-1/17.817)**self.fret                      # Period (string length) (m)
        self.xp = PluckLocation/100*StringLength          # Location of "pluck" along string (m)
        self.coefs = Coefficients+1
        self.amp = 1                              # Amplification factor used on the guitar pick-up signal to increase amplitude of sound wave

        self.v = (self.T/self.mu)**(0.5)
        self.l = self.L/2                           # Half Period (m)

        self.x_res = 100                            # Resolution of x
        self.pup = PickupLocation/100*StringLength        # Sound pickup location (m)
        self.pup_x_range_index = int(round(PickupLocation/100*self.x_res))

        self.f_res = 44100                          # Time sampling frequency (Hz)
        self.t_res = 1/self.f_res                   # Time step (s)
    
        self.sampletime = SampleTime                # Total length of sample to take (s)

        self.t_steps = int(self.sampletime*self.f_res)         # Number of time samples
        print("Generating bow object...")
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
        print("x_range has been defined.")


        self.integrand = np.zeros(len(self.x_range))

        self.bn = np.zeros(self.coefs) 

        for n in range(self.coefs):
            print("Calculating " + str(self.coefs) + " fourier coefficients with n = " + str(n))
            for i in range(len(self.x_range)):
                self.integrand[i] = self.y0(self.x_range[i])*m.sin(n*m.pi*self.x_range[i]/self.L)
            self.bn[n] = 1/self.l*np.trapz(self.integrand,self.x_range)
            print("b" + str(n) + " is " + str(self.bn[n]))

        print("Coefficients b_n have been calculated.")

        self.y00 = np.zeros(len(self.x_range))

        for i in range(len(self.x_range)):
            for n in range(len(self.bn)):
                self.y00[i] += self.bn[n]*m.sin((n)*m.pi*self.x_range[i]/self.L)

        print("Boundary condition y(x,t) at t = 0 has been calculated using the fourier coefficients bn")

        self.t_range = np.arange(0,self.t_res*self.t_steps,self.t_res)
        self.soundwave = np.zeros(self.t_steps,dtype = np.float32)
        # print("len(self.t_range) is: " + str(len(self.t_range)))
        # print("len(self.soundwave) is: " + str(len(self.soundwave)))
        # print("self.t_range.shape[0] is: " + str(self.t_range.shape[0]))
        # print("self.soundwave.shape[0] is: " + str(self.t_range.shape[0]))

        self.yxt = np.zeros((len(self.t_range),len(self.x_range)))

        for t_step in range(self.yxt.shape[0]):
            yx = np.zeros(self.yxt.shape[1])
            for i in range(self.yxt.shape[1]):
                for n in range(len(self.bn)):
                    yx[i] = yx[i] + self.bn[n]*m.sin((n)*m.pi*self.x_range[i]/self.L)*m.cos((n)*m.pi*self.v*self.t_range[t_step]/self.L)
                if(i == self.pup_x_range_index):
                    self.soundwave[t_step] = yx[i]*self.amp
                self.yxt[t_step,i] = yx[i]

        print("y(x,t) has been generated.")

        scipy.io.wavfile.write(str(self.name) + ".wav", self.f_res, self.soundwave)

        self.frequencyRange = np.fft.fftfreq(self.soundwave.shape[0], d = self.t_res)
        self.frequencyRange = self.frequencyRange[range(int((round(self.soundwave.shape[0]))/2))]
        self.transform = abs(np.fft.fft(self.soundwave)*2/self.soundwave.shape[0])
        self.transform = self.transform[range(int(round(self.soundwave.shape[0])/2))]

        print(".wav file and FFT has been generated")
        print(self.name + " - Bow generation complete.")

## x_range and y_init
    def plot_y_init(self):
        fig, ax = plt.subplots()
        ax.plot(self.x_range,self.y_init)
        plt.title(str(self.name + ": exact value of y_init (boundary case)"))
        plt.savefig(str(self.name + "_y_init.png"), bbox_inches='tight')
        plt.close()
        # plt.show()

## x_range and self.integrand
    def plot_integrand(self):
        fig, ax = plt.subplots()
        ax.plot(self.x_range,self.integrand)
        plt.title(str(self.name + ": Fourier integrand function for " + self.name))
        plt.savefig(str(self.name + "_integrand.png"), bbox_inches='tight')
        plt.close()

        # plt.show()

## x_range and y00
    def plot_y00(self):
        fig, ax = plt.subplots()
        ax.plot(self.x_range,self.y00)
        plt.title(str(self.name + ": Fourier series approximation of y_init")) 
        plt.savefig(str(self.name + "_y00.png"), bbox_inches='tight')
        plt.close()
        # plt.show()




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
        plt.savefig(str(self.name + "_yxt.png"), bbox_inches='tight')       
        plt.close()
        # plt.show()

    def plot_fft(self):
        fig, ax = plt.subplots()
        ax.plot(frequencyRange,self.transform)
        plt.xlim(1,20000)
        ax.set_xscale('log')
        plt.title(str(self.name + ": Discrete Fourier Transform"))
        plt.savefig(str(self.name + "_FFT.png"), bbox_inches='tight')
        plt.close()       


    def plot_soundwave(self):
        fig, ax = plt.subplots()
        ax.plot(self.t_range,self.soundwave)
        plt.title(str(self.name + ": Soundwave of string vibration with sound pick-up at x = " + str(self.pup))) 
        plt.savefig(str(self.name + "_wave.png"), bbox_inches='tight')
        plt.close()
        # plt.show()

    def plot_all(self):
        self.plot_y_init()
        self.plot_integrand()
        self.plot_y00()
        self.plot_soundwave()
        self.plot_yxt(20,4)

    def get_bn(self):
        return self.bn

    def soundwave(self):
        return self.soundwave

    def t_range(self):
        return self.t_range

    def transform(self):
        return self.transform

    def frequencyRange(self):
        return self.frequencyRange





#MAIN

# We will create a few digital guitar and bass strings:

t_global    = 0.1
c_global    = 14

pluck_height_global = 0.005 # m

T_guitar    = 80   # N
L_guitar    = 0.640 # m
L_bass      = 0.860 # m
mu_guitar   = 0.001 # kg/m
mu_bass     = 0.016 # kg/m

xlimlow  = 10
xlimhigh = 22100

# Relationships to show:

# 1. Tension and frequency
# 2. Density and frequency
# 3. Length and frequency
# 4. Harmonics and pick-up location



# TEST 1 - SHOW RELATIONSHIP BETWEEN TENSION AND FREQUENCY

Test_1_String1 = bow("Test_1_String1", T_guitar    , mu_guitar, pluck_height_global, 20, L_guitar, 15, 0, t_global, c_global)
Test_1_String2 = bow("Test_1_String2", T_guitar*1.2, mu_guitar, pluck_height_global, 20, L_guitar, 15, 0, t_global, c_global)
Test_1_String3 = bow("Test_1_String3", T_guitar*1.4, mu_guitar, pluck_height_global, 20, L_guitar, 15, 0, t_global, c_global)
Test_1_String4 = bow("Test_1_String4", T_guitar*1.6, mu_guitar, pluck_height_global, 20, L_guitar, 15, 0, t_global, c_global)

Test_1_String1.plot_all()
Test_1_String2.plot_all()
Test_1_String3.plot_all()
Test_1_String4.plot_all()


fig1, ax1 = plt.subplots()
fig1.set_figheight(10)
fig1.set_figwidth(10)
ax1.plot(Test_1_String1.frequencyRange,Test_1_String1.transform)
ax1.plot(Test_1_String2.frequencyRange,Test_1_String2.transform)
ax1.plot(Test_1_String3.frequencyRange,Test_1_String3.transform)
ax1.plot(Test_1_String4.frequencyRange,Test_1_String4.transform)
plt.xlim(xlimlow,xlimhigh)
ax1.set_xscale('log')
plt.title(str("Test 1 - Relationship between tension and frequency"))
ax1.legend(('String1', 'String2', 'String3', 'String4'),loc = 'right')
plt.savefig(str("Test1.png"), bbox_inches='tight')
plt.close()
gc.collect()

# # TEST 2 - SHOW RELATIONSHIP BETWEEN DENSITY AND FREQUENCY

# Test_2_String1 = bow("Test_2_String1", 45, mu_guitar    , pluck_height_global, 50, L_guitar, 50, 0, t_global, c_global)
# Test_2_String2 = bow("Test_2_String2", 45, mu_guitar*1.2, pluck_height_global, 50, L_guitar, 50, 0, t_global, c_global)
# Test_2_String3 = bow("Test_2_String3", 45, mu_guitar*1.4, pluck_height_global, 50, L_guitar, 50, 0, t_global, c_global)
# Test_2_String4 = bow("Test_2_String4", 45, mu_guitar*1.6, pluck_height_global, 50, L_guitar, 50, 0, t_global, c_global)

# Test_2_String1.plot_all()
# Test_2_String2.plot_all()
# Test_2_String3.plot_all()
# Test_2_String4.plot_all()

# fig2, ax2 = plt.subplots()
# fig2.set_figheight(10)
# fig2.set_figwidth(10)
# ax2.plot(Test_2_String1.frequencyRange,Test_2_String1.transform)
# ax2.plot(Test_2_String2.frequencyRange,Test_2_String2.transform)
# ax2.plot(Test_2_String3.frequencyRange,Test_2_String3.transform)
# ax2.plot(Test_2_String4.frequencyRange,Test_2_String4.transform)
# plt.xlim(xlimlow,xlimhigh)
# ax2.set_xscale('log')
# plt.title(str("Test 2 - Relationship between linear density and frequency"))
# ax2.legend(('String1', 'String2', 'String3', 'String4'),loc = 'right')
# plt.savefig(str("Test2.png"), bbox_inches='tight')
# plt.close()
# gc.collect()


# # TEST 3 - SHOW RELATIONSHIP BETWEEN LENGTH AND FREQUENCY

# Test_3_String1 = bow("Test_3_String1", 45, mu_guitar, pluck_height_global, 50, L_guitar    , 50, 0, t_global, c_global)
# Test_3_String2 = bow("Test_3_String2", 45, mu_guitar, pluck_height_global, 50, L_guitar*1.2, 50, 0, t_global, c_global)
# Test_3_String3 = bow("Test_3_String3", 45, mu_guitar, pluck_height_global, 50, L_guitar*1.4, 50, 0, t_global, c_global)
# Test_3_String4 = bow("Test_3_String4", 45, mu_guitar, pluck_height_global, 50, L_guitar*1.6, 50, 0, t_global, c_global)

# Test_3_String1.plot_all()
# Test_3_String2.plot_all()
# Test_3_String3.plot_all()
# Test_3_String4.plot_all()

# fig3, ax3 = plt.subplots()
# fig3.set_figheight(10)
# fig3.set_figwidth(10)
# ax3.plot(Test_3_String1.frequencyRange,Test_3_String1.transform)
# ax3.plot(Test_3_String2.frequencyRange,Test_3_String2.transform)
# ax3.plot(Test_3_String3.frequencyRange,Test_3_String3.transform)
# ax3.plot(Test_3_String4.frequencyRange,Test_3_String4.transform)
# plt.xlim(xlimlow,xlimhigh)
# ax3.set_xscale('log')
# plt.title(str("Test 3 - Relationship between length and frequency"))
# ax3.legend(('String1', 'String2', 'String3', 'String4'),loc = 'right')
# plt.savefig(str("Test3.png"), bbox_inches='tight')
# plt.close()
# gc.collect()


# # TEST 3 - SHOW RELATIONSHIP BETWEEN HARMONICS AND PICKUP LOCATION

# Test_4_String1 = bow("Test_4_String1", 45, mu_guitar, pluck_height_global, 50, L_guitar, 20, 0, t_global, c_global)
# Test_4_String2 = bow("Test_4_String2", 45, mu_guitar, pluck_height_global, 50, L_guitar, 30, 0, t_global, c_global)
# Test_4_String3 = bow("Test_4_String3", 45, mu_guitar, pluck_height_global, 50, L_guitar, 40, 0, t_global, c_global)
# Test_4_String4 = bow("Test_4_String4", 45, mu_guitar, pluck_height_global, 50, L_guitar, 50, 0, t_global, c_global)

# Test_4_String1.plot_all()
# Test_4_String2.plot_all()
# Test_4_String3.plot_all()
# Test_4_String4.plot_all()

# fig4, ax4 = plt.subplots()
# fig4.set_figheight(10)
# fig4.set_figwidth(10)
# ax4.plot(Test_4_String1.frequencyRange,Test_4_String1.transform)
# ax4.plot(Test_4_String2.frequencyRange,Test_4_String2.transform)
# ax4.plot(Test_4_String3.frequencyRange,Test_4_String3.transform)
# ax4.plot(Test_4_String4.frequencyRange,Test_4_String4.transform)
# plt.xlim(xlimlow,xlimhigh)
# ax4.set_xscale('log')
# plt.title(str("Test 4 - Relationship between harmonics and pickup location"))
# ax4.legend(('String1', 'String2', 'String3', 'String4'),loc = 'right')
# plt.savefig(str("Test4.png"), bbox_inches='tight')
# plt.close()
# gc.collect()













# def __init__(self, name, Tension, LinDens, PluckHeight, PluckLocation, StringLength, PickupLocation, Fret, SampleTime, Coefficients):


# # Control Guitar
# Guitar1 = bow("Guitar1", 45, mu_guitar, pluck_height_global, 15, L_guitar, 20, 0, t_global, c_global)
# Guitar1.plot_fft()
# Guitar1.plot_y_init()
# Guitar1.plot_integrand()
# Guitar1.plot_y00()
# Guitar1.plot_soundwave()
# Guitar1.plot_yxt(10,0)

# gc.collect()

# # Change the pluck location for Guitar2:
# Guitar2 = bow("Guitar2", 45, mu_guitar, pluck_height_global, 5, L_guitar, 20, 0, t_global, c_global)
# Guitar2.plot_y_init()
# Guitar2.plot_integrand()
# Guitar2.plot_y00()
# Guitar2.plot_soundwave()
# Guitar2.plot_yxt(10,0)
# gc.collect()

# # Change the pick-up location for Guitar3:
# Guitar3 = bow("Guitar3", 45, mu_guitar, pluck_height_global, 15, L_guitar, 40, 0, t_global, c_global)
# Guitar3.plot_y_init()
# Guitar3.plot_integrand()
# Guitar3.plot_y00()
# Guitar3.plot_soundwave()
# Guitar3.plot_yxt(10,0)
# gc.collect()

# # Control Bass
# Bass1 = bow("Bass1", 100, mu_bass, pluck_height_global, 15, L_bass, 20, 0, t_global, c_global)
# Bass1.plot_y_init()
# Bass1.plot_integrand()
# Bass1.plot_y00()
# Bass1.plot_soundwave()
# Bass1.plot_yxt(10,0)
# gc.collect()

# # Change the tension for Bass2:
# Bass2 = bow("Bass2", 130, mu_bass, pluck_height_global, 15, L_bass, 20, 0, t_global, c_global)
# Bass2.plot_y_init()
# Bass2.plot_integrand()
# Bass2.plot_y00()
# Bass2.plot_soundwave()
# Bass2.plot_yxt(10,0)
# gc.collect()

# # Change the string density for Bass3 by + 15%:
# Bass3 = bow("Bass3", 100, mu_bass*1.15, pluck_height_global, 15, L_bass, 20, 0, t_global, c_global)
# Bass3.plot_y_init()
# Bass3.plot_integrand()
# Bass3.plot_y00()
# Bass3.plot_soundwave()
# Bass3.plot_yxt(10,0)
# gc.collect()



# Guitar4_0  = bow("Guitar4_0",  45, mu_guitar, pluck_height_global, 5, L_guitar, 20, 0 , t_global, c_global)
# gc.collect()
# Guitar4_1  = bow("Guitar4_1",  45, mu_guitar, pluck_height_global, 5, L_guitar, 20, 1 , t_global, c_global)
# gc.collect()
# Guitar4_2  = bow("Guitar4_2",  45, mu_guitar, pluck_height_global, 5, L_guitar, 20, 2 , t_global, c_global)
# gc.collect()
# Guitar4_3  = bow("Guitar4_3",  45, mu_guitar, pluck_height_global, 5, L_guitar, 20, 3 , t_global, c_global)
# gc.collect()
# Guitar4_4  = bow("Guitar4_4",  45, mu_guitar, pluck_height_global, 5, L_guitar, 20, 4 , t_global, c_global)
# gc.collect()
# Guitar4_5  = bow("Guitar4_5",  45, mu_guitar, pluck_height_global, 5, L_guitar, 20, 5 , t_global, c_global)
# gc.collect()
# Guitar4_6  = bow("Guitar4_6",  45, mu_guitar, pluck_height_global, 5, L_guitar, 20, 6 , t_global, c_global)
# gc.collect()
# Guitar4_7  = bow("Guitar4_7",  45, mu_guitar, pluck_height_global, 5, L_guitar, 20, 7 , t_global, c_global)
# gc.collect()
# Guitar4_8  = bow("Guitar4_8",  45, mu_guitar, pluck_height_global, 5, L_guitar, 20, 8 , t_global, c_global)
# gc.collect()
# Guitar4_9  = bow("Guitar4_9",  45, mu_guitar, pluck_height_global, 5, L_guitar, 20, 9 , t_global, c_global)
# gc.collect()
# Guitar4_10 = bow("Guitar4_10", 45, mu_guitar, pluck_height_global, 5, L_guitar, 20, 10, t_global, c_global)
# gc.collect()
# Guitar4_11 = bow("Guitar4_11", 45, mu_guitar, pluck_height_global, 5, L_guitar, 20, 11, t_global, c_global)
# gc.collect()
# Guitar4_12 = bow("Guitar4_12", 45, mu_guitar, pluck_height_global, 5, L_guitar, 20, 12, t_global, c_global)
# gc.collect()

# gc.collect()
