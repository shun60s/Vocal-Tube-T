#coding:utf-8

#
# T Three Tube Model (3-way junction model) 
#


import numpy as np
from matplotlib import pyplot as plt


# Check version
#  Python 3.10.4, 64bit on Win32 (Windows 10)
#  numpy 1.21.6


class Class_T_ThreeTube(object):
    def __init__(self, L1, L2, L3, A1, A2, A3, rg0=0.95, rl0=0.9 ,sampling_rate=48000):
        # initalize Tube length and Tube area
        self.L1= L1 # set list of 1st tube's length by unit is [cm]
        self.A1= A1 # set list of 1st tube's area by unit is [cm^2]
        self.L2= L2 # set list of 2nd tube's length by unit is [cm]
        self.A2= A2 # set list of 2nd tube's area by unit is [cm^2]
        self.L3= L3 # set list of 3rd tube's length by unit is [cm]
        self.A3= A3 # set list of 3rd tube's area by unit is [cm^2]
        C0=35000.0  # speed of sound in air, round 35000 cm/second
        self.sr= sampling_rate
        self.tu1=self.L1 / C0   # delay time in 1st tube
        self.tu2=self.L2 / C0   # delay time in 2nd tube
        self.tu3=self.L3 / C0   # delay time in 3rd tube
        self.r12=((self.A3+self.A2) - self.A1)/(self.A3+self.A2+self.A1)  # reflection coefficient between 1st tube and others
        self.r21=(self.A2 - (self.A3+self.A1))/(self.A3+self.A2+self.A1)  # reflection coefficient between 2nd tube and others
        self.r31=(self.A3 - (self.A2+self.A1))/(self.A3+self.A2+self.A1)  # reflection coefficient between 3rd tube and others


        self.rg0=rg0 # rg is reflection coefficient between glottis and 1st tube
        self.rl0=rl0 # reflection coefficient between 2nd tube and mouth
        self.rl3=-0.97 #set certain a value, beside ideal value is -1
        
        self.att_norm=1.0 # attenuation constant per one step ahead
        
        self.num_of_tube=3
        
    def fone(self, xw):
        # calculate one point of frequecny response
        yi= 0.5 * ( 1.0 + self.rg0 ) * ( 1.0 + self.r21)  * ( 1.0 + self.rl0 ) * np.exp( -1.0j * ( self.tu1 + self.tu2 ) * xw) 
        yb1= 1.0 + self.r12 * self.rg0 *  np.exp( -2.0j * self.tu1 * xw ) 
        yb1= yb1 + self.r21  * self.rl0 *  np.exp( -2.0j * self.tu2 * xw ) 
        yb1= yb1 + self.rg0 * self.rl0 * (1.0 - self.r12 + self.r21) * np.exp( -2.0j * ( self.tu1 + self.tu2 ) * xw)
        yc1= self.rl3  * ( 1.0 + self.r31)
        yc1= yc1 * (1.0 + self.rg0 *  np.exp( -2.0j * self.tu1 * xw ))
        yc1= yc1 * (1.0 - self.rl0 * np.exp( -2.0j * self.tu2 * xw ))
        yc1= yc1 * np.exp( -2.0j *  self.tu3 * xw) #np.exp( -1.0j *  self.tu3 * xw)
        yd1= 1.0 - self.rl3 * np.exp( -2.0j * self.tu3 * xw ) 
        """
        if yd1 == 0.0:
            val=0.  # because yc/yd is infinity when yd = 0
        else:
            yc1= yc1 / yd1
            val= yi / ( yb1 + yc1)
        """
        yc1= yc1 / yd1
        val= yi / ( yb1 + yc1)
        return np.sqrt(val.real ** 2 + val.imag ** 2)
        
    def H0(self, freq_low=100, freq_high=5000, Band_num=256):
        # get Log scale frequecny response, from freq_low to freq_high, Band_num points
        amp=[]
        freq=[]
        bands= np.zeros(Band_num+1)
        fcl=freq_low * 1.0    # convert to float
        fch=freq_high * 1.0   # convert to float
        delta1=np.power(fch/fcl, 1.0 / (Band_num)) # Log Scale
        bands[0]=fcl
        #print ("i,band = 0", bands[0])
        for i in range(1, Band_num+1):
            bands[i]= bands[i-1] * delta1
            #print ("i,band =", i, bands[i]) 
        for f in bands:
            amp.append(self.fone(f * 2.0 * np.pi))
        return   np.log10(amp) * 20, bands # = amp value, freq list
        
    def process(self, yg ):
        # process reflection transmission of resonance tube: yg is input, y2tm is output
        # three T tube
        #
        #
        #             ------------------------------------------------
        # yg->  ya1-> |                                yb1->          |->  y2tm
        #          rg |     tube1        <-ya2                tube2   |rl0
        #          <-  -------------------             yc1            |<- yb2
        #                                |              |     -------
        #                                |              V     |
        #                                |      tube3         |
        #                                |                    |
        #                                |                    |
        #                                |                    |
        #                                ---------------------
        #                                 ^
        #                                 |       rl3
        #                                yc2
        #
        #
        M1= round( self.tu1 * self.sr ) + 1  # for precision, higher sampling_rate is better
        M2= round( self.tu2 * self.sr ) + 1  # for precision, higher sampling_rate is better
        M3= round( self.tu3 * self.sr ) + 1  # for precision, higher sampling_rate is better
        M1= int(M1)
        M2= int(M2)
        M3= int(M3)
        ya1=np.zeros(M1)
        ya2=np.zeros(M1)
        yb1=np.zeros(M2)
        yb2=np.zeros(M2)
        yc1=np.zeros(M3)
        yc2=np.zeros(M3)
        y2tm=np.zeros(len(yg))
        
        for tc0 in range(len(yg)):
            for i in range((M1-1),0,-1): # process one step
                ya1[i]=ya1[i-1] * self.att_norm
                ya2[i]=ya2[i-1] * self.att_norm
            for i in range((M2-1),0,-1): # process one step
                yb1[i]=yb1[i-1] * self.att_norm
                yb2[i]=yb2[i-1] * self.att_norm
            for i in range((M3-1),0,-1): # process one step
                yc1[i]=yc1[i-1] * self.att_norm
                yc2[i]=yc2[i-1] * self.att_norm
            # calculate reflection
            #  tube1 input
            ya1[0]= ((1. + self.rg0 ) / 2.) * yg[tc0] + self.rg0 * ya2[-1]
            ya2[0]= -1. * self.r12 * ya1[-1] + ( 1. - self.r12 ) * yb2[-1] + ( 1.  - self.r12) * yc2[-1]
            #  tube2  output
            yb1[0]= ( 1 + self.r21 ) * ya1[-1] + self.r21 * yb2[-1] + ( 1. + self.r21) * yc2[-1]
            yb2[0]=  -1. * self.rl0  * yb1[-1]
            y2tm[tc0]= (1 + self.rl0) * yb1[-1]
            
            #  tube 3 closed
            yc1[0]= ( 1. + self.r31 ) * ya1[-1] + self.r31 * yc2[-1] + ( 1. + self.r31 ) * yb2[-1]
            yc2[0]=  -1. * self.rl3  * yc1[-1]
            
            
        return y2tm


    def check1(self,):
        # check if sum of output, ya1, yb2, and yc2 is 1
        a1= -1. * self.r12  +  ( 1 + self.r21 ) + ( 1. + self.r31 )
        b2= ( 1. - self.r12 )   +  self.r21 + ( 1. + self.r31 )
        c2= ( 1. - self.r12 )   +  ( 1. + self.r21)  +  self.r31 
        print('check1',a1,b2,c2)




if __name__ == '__main__':
    
    # Length & Area value, from problems 3.8 in "Digital Processing of Speech Signals" by L.R.Rabiner and R.W.Schafer
    #
    # /a/
    L1_a=9.0    # set list of 1st tube's length by unit is [cm]
    A1_a=1.0    # set list of 1st tube's area by unit is [cm^2]
    L2_a=8.0    # set list of 2nd tube's length by unit is [cm]
    A2_a=7.0    # set list of 2nd tube's area by unit is [cm^2]
    
    # /u/
    L1_u=10.0   # set list of 1st tube's length by unit is [cm]
    A1_u=7.0    # set list of 1st tube's area by unit is [cm^2]
    L2_u=7.0    # set list of 2nd tube's length by unit is [cm]
    A2_u=3.0    # set list of 2nd tube's area by unit is [cm^2]
    
    # /o/: L3,A3 is  extend factor to /a/ connecting as /u/
    L3_o= L2_a * (L2_u / L1_u)     # set list of 3rd tube's length by unit is [cm]
    A3_o= A2_a * (A2_u / A1_u)     # set list of 3rd tube's area by unit is [cm^2]
    
    # insatnce
    tube  =  Class_T_ThreeTube(L1_a, L2_a, L3_o, A1_a, A2_a, A3_o, sampling_rate=48000)
    
    tube.check1()
    
    
    