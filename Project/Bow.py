#=================================#
#       ENGR-5696 PROJECT         #
#                                 #
#         JEFF EASTICK            #
#          FALL 2018              #
#=================================#

import numpy as np
import math as m
import matplotlib
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import scipy.io.wavfile
import gc
import time
import matplotlib.animation as animation


from mpl_toolkits.mplot3d import axes3d

class bow:

    def y0(self,x):

        if(np.where(x<=self.xp)):
            return self.yp/self.xp*x
        if(np.where(x>self.xp)):  
            return self.yp/(self.L-self.xp)*(self.L-x)

    def bn_exact(self,N):
        return 2*self.yp*self.L*(np.sin(m.pi*(N)*self.xp/self.L))*(1/self.xp + 1/(self.L-self.xp))/(np.pi**2)/(N**2)


    def yxt_exact(self,X,T,N):
        a=0
        for n in range(N):
            n += 1
            a += self.bn_exact(n)*np.sin((n)*np.pi*X/self.L)*np.cos((n)*np.pi*self.v*T/self.L)
        return a



    def __init__(self, name, Tension, LinDens, PluckHeight, PluckLocation, StringLength, PickupLocation, Fret, SampleTime, Coefficients):
        
        self.name = name
        print("-----------------START " + self.name + "---------------------")
        print("Initializing + " + self.name)
        self.T = Tension                            # Tension of material (Newtons)
        self.mu = LinDens                           # Linear Density of material (kg/m)
        self.yp = PluckHeight                       # Height of "pluck" from neutral axis (m)
        self.fret = Fret
        self.L = StringLength*(1-1/17.817)**self.fret                      # Period (string length) (m)
        self.xp = PluckLocation/100*StringLength          # Location of "pluck" along string (m)
        self.coefs = Coefficients
        self.amp = 5                                # Amplification factor used on the guitar pick-up signal to increase amplitude of sound wave

        self.v = (self.T/self.mu)**(0.5)
        self.l = self.L/2                           # Half Period (m)
        self.fund = self.v/(2*self.L)
        self.wavelength = self.v/self.fund
        self.x_res = 500                            # Resolution of x
        self.pup = PickupLocation/100*StringLength        # Sound pickup location (m)
        self.pup_x_range_index = int(round(PickupLocation/100*self.x_res))

        self.f_res = 5000                          # Time sampling frequency (Hz)
        self.t_res = 1/self.f_res                   # Time step (s)
    
        self.sampletime = SampleTime                # Total length of sample to take (s)

        self.t_steps = int(self.sampletime*self.f_res)         # Number of time samples


        self.bns = np.fromfunction(lambda n: self.bn_exact(n),(self.coefs,))
        # MAKING MESHGRID ARRAYS


        self.space   = np.linspace(0, self.L         , self.x_res)
        self.time    = np.linspace(0, self.sampletime, self.f_res)


        self.bc = np.piecewise(self.space,[self.space<=self.xp,self.space>self.xp],[lambda x: self.yp/self.xp*x, lambda x: self.yp/(self.L-self.xp)*(self.L-x)])

        self.ux, self.ut = np.meshgrid(self.space,self.time)

        self.realyxt = self.yxt_exact(self.ux,self.ut,self.coefs)

        print("Generating .wav file.")

        self.generateWav()


        print("Generating FFT.")

        self.frequencyRange = np.fft.fftfreq(self.wave.shape[0], d = 1/44100)
        self.frequencyRange = self.frequencyRange[range(int((round(self.wave.shape[0]))/2))]
        self.transform = abs(np.fft.fft(self.wave)*2/self.wave.shape[0])
        self.transform = self.transform[range(int(round(self.wave.shape[0])/2))]



        # print("Plotting graphs.")        
        # self.plot_all()

        print("-------------------END " + self.name + "---------------------")
        print("")

    def generateWav(self):
        length = 10

        self.wavetime    = np.linspace(0, length, length*44100)
        self.wave = self.yxt_exact(self.pup,self.wavetime,self.coefs)
        self.wave = self.wave.astype(np.float32)
        scipy.io.wavfile.write(str(self.name) + ".wav", 44100, self.wave*self.amp)


    def plot_all(self):
        self.plot_bn()
        self.plot_wav()
        self.plot_bc()
        # self.plot_wires_waves()
        # self.plot_wires_string()
        # self.plot_fft()


    def plot_bn(self):
        fig, ax = plt.subplots()
        ax.plot(self.bns)
        plt.title(str(self.name + ": Fourier coefficients")) 
        plt.savefig(str(self.name + "_bn.png"), bbox_inches='tight')
        plt.close()

        # TODO - fix formatting and axes

    def plot_wav(self):
        self.generateWav()
        fig, ax = plt.subplots()
        ax.plot(self.wavetime,self.wave)
        ax.set_xlim(0,0.01)
        plt.title(str(self.name + ": Sound wave")) 
        plt.savefig(str(self.name + "_wav.png"), bbox_inches='tight')
        plt.close()


    def plot_bc(self):
        fig, ax = plt.subplots()
        ax.plot(self.space,self.yxt_exact(self.space,0,self.coefs))
        ax.plot(self.space,self.bc)
        plt.title(str(self.name + ": Initial Boundary Condition")) 
        plt.savefig(str(self.name + "_bc.png"), bbox_inches='tight')
        plt.close()

    def plot_wires_waves(self):

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        plt.title(str(self.name + ": Sound waves at various pick-up points"))
        ax.plot_wireframe(self.ux, self.ut, self.realyxt,rstride = self.f_res, cstride = int(self.x_res/5))
        ax.set_xlabel('Position along string (x)')
        ax.set_ylabel('Time (t)')
        ax.set_zlabel('Amplitude y(x,t)')
        plt.show()


    def plot_wires_string(self):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        plt.title(str(self.name + ": String position y(x,t) at various times t"))
        ax.plot_wireframe(self.ux, self.ut, self.realyxt,rstride = 150, cstride=0)
        ax.set_xlabel('Position along string (x)')
        ax.set_ylabel('Time (t)')
        ax.set_zlabel('Amplitude y(x,t)')
        plt.show()

    def plot_fft(self):
        fig, ax = plt.subplots()
        ax.plot(self.frequencyRange,self.transform)
        plt.xlim(1,20000)
        ax.set_xlim(20,20000)
        ax.set_ylim(0.00001,1)
        ax.set_yscale('log')
        ax.set_xscale('log')
        plt.title(str(self.name + ": Discrete Fourier Transform"))
        plt.savefig(str(self.name + "_FFT.png"), bbox_inches='tight')
        plt.close()       





# legacy code that calculates Bn numerically and is also ugly and may not work anymore:
    def generateBow(self):

        self.x_range = np.linspace(0,self.L,num=self.x_res+1)

        self.y_init = np.zeros(len(self.x_range))

        for i in range(len(self.x_range)):
            self.y_init[i] = self.y0(self.x_range[i])
        print("y_0 has been defined.")


        self.integrand = np.zeros(len(self.x_range))

        self.bn = np.zeros(self.coefs+1) 

        for n in range(self.coefs):
            n+=1
            print("Calculating " + str(self.coefs) + " fourier coefficients with n = " + str(n))
            for i in range(len(self.x_range)):
                self.integrand[i] = self.y0(self.x_range[i])*m.sin(n*m.pi*self.x_range[i]/self.L)
            self.bn[n] = 2/self.L*np.trapz(self.integrand,self.x_range)
            
            print("b_approx_" + str(n) + " is " + str(self.bn[n]))
            self.bn[n] = self.bn_exact(n)
            print("b_exact_" + str(n) + " is " + str(self.bn[n]))

        print("Coefficients b_n have been calculated.")

        self.y00 = np.zeros(len(self.x_range))

        for i in range(len(self.x_range)):
            for n in range(len(self.bn)):
                self.y00[i] += self.bn[n]*m.sin((n)*m.pi*self.x_range[i]/self.L)

        print("Boundary condition y(x,t) at t = 0 has been calculated using the fourier coefficients bn")

        self.plot_y00()
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



        print(".wav file and FFT has been generated")
        print(self.name + " - Bow generation complete.")








#MAIN

# We will create a few digital guitar and bass strings:

t_global    = 0.01
c_global    = 10

pluck_height_global = 0.1 # m


T_guitar    = 120   # N
L_guitar    = 0.640 # m
L_bass      = 0.860 # m
mu_guitar   = 0.001 # kg/m
mu_bass     = 0.016 # kg/m

xlimlow  = 10
xlimhigh = 22050
ylimlow  = 0.00001
ylimhigh = pluck_height_global

fs = 'large'
# Relationships to show:

# 1. Tension and frequency
# 2. Density and frequency
# 3. Length and frequency
# 4. Harmonics and pick-up location
# 5. Pluck location and fourier series coefficients



# # TEST 1 - RELATIONSHIP BETWEEN TENSION AND FREQUENCY

Test_1_String1 = bow("Test_1_String1", T_guitar    , mu_guitar, pluck_height_global, 50, L_guitar, 50, 0, t_global, c_global)
Test_1_String2 = bow("Test_1_String2", T_guitar*1.2, mu_guitar, pluck_height_global, 50, L_guitar, 50, 0, t_global, c_global)
Test_1_String3 = bow("Test_1_String3", T_guitar*1.4, mu_guitar, pluck_height_global, 50, L_guitar, 50, 0, t_global, c_global)
Test_1_String4 = bow("Test_1_String4", T_guitar*1.6, mu_guitar, pluck_height_global, 50, L_guitar, 50, 0, t_global, c_global)

fig1, ax1 = plt.subplots()
fig1.set_figheight(10)
fig1.set_figwidth(10)
ax1.plot(Test_1_String1.frequencyRange,Test_1_String1.transform)
ax1.plot(Test_1_String2.frequencyRange,Test_1_String2.transform)
ax1.plot(Test_1_String3.frequencyRange,Test_1_String3.transform)
ax1.plot(Test_1_String4.frequencyRange,Test_1_String4.transform)
ax1.set_xlim(xlimlow,xlimhigh)
ax1.set_ylim(ylimlow,ylimhigh)
ax1.set_xscale('log')
plt.xlabel('Hz', fontsize = fs)
plt.ylabel('Amplitude', fontsize = fs)
plt.title(str("Test 1 - Relationship between tension and frequency"), fontsize = fs)
ax1.legend(('String1 - T = ' + str(Test_1_String1.T) + ' N, f = ' + str(round(Test_1_String1.fund,2)) + ' Hz', 
            'String2 - T = ' + str(Test_1_String2.T) + ' N, f = ' + str(round(Test_1_String2.fund,2)) + ' Hz',
            'String3 - T = ' + str(Test_1_String3.T) + ' N, f = ' + str(round(Test_1_String3.fund,2)) + ' Hz', 
            'String4 - T = ' + str(Test_1_String4.T) + ' N, f = ' + str(round(Test_1_String4.fund,2)) + ' Hz'),
loc = 'right', fontsize = fs)
plt.savefig(str("Test1_T_vs_freq.png"), bbox_inches='tight')
plt.close()
del Test_1_String1
del Test_1_String2
del Test_1_String3
del Test_1_String4
gc.collect()

# # # TEST 2 - SHOW RELATIONSHIP BETWEEN DENSITY AND FREQUENCY

Test_2_String1 = bow("Test_2_String1", T_guitar, mu_guitar    , pluck_height_global, 50, L_guitar, 50, 0, t_global, c_global)
Test_2_String2 = bow("Test_2_String2", T_guitar, mu_guitar*1.2, pluck_height_global, 50, L_guitar, 50, 0, t_global, c_global)
Test_2_String3 = bow("Test_2_String3", T_guitar, mu_guitar*1.4, pluck_height_global, 50, L_guitar, 50, 0, t_global, c_global)
Test_2_String4 = bow("Test_2_String4", T_guitar, mu_guitar*1.6, pluck_height_global, 50, L_guitar, 50, 0, t_global, c_global)

fig2, ax2 = plt.subplots()
fig2.set_figheight(10)
fig2.set_figwidth(10)
ax2.plot(Test_2_String1.frequencyRange,Test_2_String1.transform)
ax2.plot(Test_2_String2.frequencyRange,Test_2_String2.transform)
ax2.plot(Test_2_String3.frequencyRange,Test_2_String3.transform)
ax2.plot(Test_2_String4.frequencyRange,Test_2_String4.transform)
plt.xlim(xlimlow,xlimhigh)
ax2.set_xscale('log')
plt.xlabel('Hz', fontsize = fs)
plt.ylabel('Amplitude', fontsize = fs)
plt.title(str("Test 2 - Relationship between linear density and frequency"), fontsize = fs)
ax2.legend(('String1 - density = ' + str(Test_2_String1.mu) + ' kg/m, f = ' + str(round(Test_2_String1.fund,2)) + ' Hz', 
            'String2 - density = ' + str(Test_2_String2.mu) + ' kg/m, f = ' + str(round(Test_2_String2.fund,2)) + ' Hz', 
            'String3 - density = ' + str(Test_2_String3.mu) + ' kg/m, f = ' + str(round(Test_2_String3.fund,2)) + ' Hz', 
            'String4 - density = ' + str(Test_2_String4.mu) + ' kg/m, f = ' + str(round(Test_2_String4.fund,2)) + ' Hz'),loc = 'right', fontsize = fs)
plt.savefig(str("Test2_density_vs_freq.png"), bbox_inches='tight')
plt.close()

del Test_2_String1
del Test_2_String2
del Test_2_String3
del Test_2_String4
gc.collect()


# # # TEST 3 - SHOW RELATIONSHIP BETWEEN LENGTH AND FREQUENCY

Test_3_String1 = bow("Test_3_String1", T_guitar, mu_guitar, pluck_height_global, 50, L_guitar    , 50, 0, t_global, c_global)
Test_3_String2 = bow("Test_3_String2", T_guitar, mu_guitar, pluck_height_global, 50, L_guitar*1.2, 50, 0, t_global, c_global)
Test_3_String3 = bow("Test_3_String3", T_guitar, mu_guitar, pluck_height_global, 50, L_guitar*1.4, 50, 0, t_global, c_global)
Test_3_String4 = bow("Test_3_String4", T_guitar, mu_guitar, pluck_height_global, 50, L_guitar*1.6, 50, 0, t_global, c_global)

fig3, ax3 = plt.subplots()
fig3.set_figheight(10)
fig3.set_figwidth(10)
ax3.plot(Test_3_String1.frequencyRange,Test_3_String1.transform)
ax3.plot(Test_3_String2.frequencyRange,Test_3_String2.transform)
ax3.plot(Test_3_String3.frequencyRange,Test_3_String3.transform)
ax3.plot(Test_3_String4.frequencyRange,Test_3_String4.transform)
plt.xlim(xlimlow,xlimhigh)
ax3.set_xscale('log')
plt.xlabel('Hz', fontsize = fs)
plt.ylabel('Amplitude', fontsize = fs)
plt.title(str("Test 3 - Relationship between length and frequency"), fontsize = fs)
ax3.legend(('String1 - L = ' + str(round(Test_3_String1.L,3)) + ' m, f = ' + str(round(Test_3_String1.fund,2)) + ' Hz', 
            'String2 - L = ' + str(round(Test_3_String2.L,3)) + ' m, f = ' + str(round(Test_3_String2.fund,2)) + ' Hz', 
            'String3 - L = ' + str(round(Test_3_String3.L,3)) + ' m, f = ' + str(round(Test_3_String3.fund,2)) + ' Hz', 
            'String4 - L = ' + str(round(Test_3_String4.L,3)) + ' m, f = ' + str(round(Test_3_String4.fund,2)) + ' Hz'),loc = 'right', fontsize = fs)
plt.savefig(str("Test3_length_vs_freq.png"), bbox_inches='tight')
plt.close()

del Test_3_String1
del Test_3_String2
del Test_3_String3
del Test_3_String4
gc.collect()


# # TEST 4 - SHOW SOUNDWAVE AMPLITUDE AND HARMONIC CONTENT CHANGES WITH PICKUP LOCATION

Test_4_String1 = bow("Test_4_String1", T_guitar, mu_guitar, pluck_height_global, 20, L_guitar, 20, 0, t_global, c_global)
Test_4_String2 = bow("Test_4_String2", T_guitar, mu_guitar, pluck_height_global, 20, L_guitar, 30, 0, t_global, c_global)
Test_4_String3 = bow("Test_4_String3", T_guitar, mu_guitar, pluck_height_global, 20, L_guitar, 40, 0, t_global, c_global)
Test_4_String4 = bow("Test_4_String4", T_guitar, mu_guitar, pluck_height_global, 20, L_guitar, 50, 0, t_global, c_global)

fig4 = plt.figure()
plt.axis('off')

ax41 = fig4.add_subplot(4,1,1)
ax41.plot(Test_4_String1.frequencyRange,Test_4_String1.transform)
ax41.set_xlim(xlimlow,xlimhigh)
ax41.set_ylim(ylimlow,ylimhigh)
ax41.set_xscale('log')

ax42 = fig4.add_subplot(4,1,2)
ax42.plot(Test_4_String2.frequencyRange,Test_4_String2.transform)
ax42.set_xlim(xlimlow,xlimhigh)
ax42.set_ylim(ylimlow,ylimhigh)
ax42.set_xscale('log')


ax43 = fig4.add_subplot(4,1,3)
ax43.plot(Test_4_String3.frequencyRange,Test_4_String3.transform)
ax43.set_xlim(xlimlow,xlimhigh)
ax43.set_ylim(ylimlow,ylimhigh)
ax43.set_xscale('log')


ax44 = fig4.add_subplot(4,1,4)
ax44.plot(Test_4_String4.frequencyRange,Test_4_String4.transform)
ax44.set_xlim(xlimlow,xlimhigh)
ax44.set_ylim(ylimlow,ylimhigh)
ax44.set_xscale('log')

ax41.set_title(('Test 4, String1 - pup @ x = ' + str(round(Test_4_String1.pup,3)) + " m, L = " + str(round(Test_4_String1.L,3))), fontsize = fs)
ax42.set_title(('Test 4, String2 - pup @ x = ' + str(round(Test_4_String2.pup,3)) + " m, L = " + str(round(Test_4_String2.L,3))), fontsize = fs)
ax43.set_title(('Test 4, String3 - pup @ x = ' + str(round(Test_4_String3.pup,3)) + " m, L = " + str(round(Test_4_String3.L,3))), fontsize = fs)
ax44.set_title(('Test 4, String4 - pup @ x = ' + str(round(Test_4_String4.pup,3)) + " m, L = " + str(round(Test_4_String4.L,3))), fontsize = fs)

fig4.set_figheight(10)
fig4.set_figwidth(10)
fig4.subplots_adjust(hspace = 0.5)
plt.savefig(str("Test4_pickup.png"), bbox_inches='tight')
plt.close()

del Test_4_String1
del Test_4_String2
del Test_4_String3
del Test_4_String4
gc.collect()


# TEST 5 - SHOW bn CHANGES WITH PLUCK LOCATION
Test_5_String0 = bow("Test_5_String0", T_guitar, mu_guitar, pluck_height_global, 10, L_guitar, 50, 0, t_global, c_global)
Test_5_String1 = bow("Test_5_String1", T_guitar, mu_guitar, pluck_height_global, 20, L_guitar, 50, 0, t_global, c_global)
Test_5_String2 = bow("Test_5_String2", T_guitar, mu_guitar, pluck_height_global, 30, L_guitar, 50, 0, t_global, c_global)
Test_5_String3 = bow("Test_5_String3", T_guitar, mu_guitar, pluck_height_global, 40, L_guitar, 50, 0, t_global, c_global)
Test_5_String4 = bow("Test_5_String4", T_guitar, mu_guitar, pluck_height_global, 50, L_guitar, 50, 0, t_global, c_global)

fig5, ax5 = plt.subplots()
plt.title(str("Test 5 - Relationship between bn and pluck location"), fontsize = fs)
plt.xlabel('n', fontsize = fs)
plt.ylabel('Fourier coefficient \'Bn\'', fontsize = fs)
# plt.axis('off')
ax5.plot(Test_5_String0.bns)
ax5.plot(Test_5_String1.bns)
ax5.plot(Test_5_String2.bns)
ax5.plot(Test_5_String3.bns)
ax5.plot(Test_5_String4.bns)
ax5.grid(True,which='both')
plt.xticks(np.arange(1,c_global,1))
ax5.legend(('String0 - pluck @ x = ' + str(round(Test_5_String0.xp,3)) + ' m', 
            'String1 - pluck @ x = ' + str(round(Test_5_String1.xp,3)) + ' m', 
            'String2 - pluck @ x = ' + str(round(Test_5_String2.xp,3)) + ' m', 
            'String3 - pluck @ x = ' + str(round(Test_5_String3.xp,3)) + ' m', 
            'String4 - pluck @ x = ' + str(round(Test_5_String4.xp,3)) + ' m'),loc = 'right', fontsize = fs)

fig5.set_figheight(10)
fig5.set_figwidth(10)

plt.savefig(str("Test5_bn_vs_xp.png"), bbox_inches='tight')
plt.close()

Test_5_String0.plot_wires_waves()
Test_5_String0.plot_wires_string()
Test_5_String4.plot_wires_waves()
Test_5_String4.plot_wires_string()

del Test_5_String0
del Test_5_String1
del Test_5_String2
del Test_5_String3
del Test_5_String4
gc.collect()









# Piano sound generator demo - generates sample sets for three 88-key pianos. 
hammer = 15
pickup = 20

for string in range (88):
    string+=1
    f = ((2**(1/12))**(string-49)*440)
    print("Key number " + str(string) + " frequency is " + str(f))
    play = bow("3string_"+str(string), mu_guitar*(2*L_guitar*f)**2, mu_guitar, pluck_height_global, hammer, L_guitar, pickup, 0, t_global, 3)

for string in range (88):
    string+=1
    f = ((2**(1/12))**(string-49)*440)
    print("Key number " + str(string) + " frequency is " + str(f))
    play = bow("10string_"+str(string), mu_guitar*(2*L_guitar*f)**2, mu_guitar, pluck_height_global, hammer, L_guitar, pickup, 0, t_global, 10)


for string in range (88):
    string+=1
    f = ((2**(1/12))**(string-49)*440)
    print("Key number " + str(string) + " frequency is " + str(f))
    play = bow("30string_"+str(string), mu_guitar*(2*L_guitar*f)**2, mu_guitar, pluck_height_global, hammer, L_guitar, pickup, 0, t_global, 30)



# Diddley bow sound generator demo - generates sample set for single-string instrument with 12 frets. 
for fret in range (13):
    f = 440 #Hz
    play = bow("fret_"+str(fret), mu_guitar*(2*L_guitar*f)**2, mu_guitar, pluck_height_global, hammer, L_guitar, pickup, fret, t_global, c_global)

