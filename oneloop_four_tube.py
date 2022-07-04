#coding:utf-8

#
# one-loop Four Tube Model, combination of T three tube, with attenuation
#


import numpy as np
from matplotlib import pyplot as plt


# Check version
#  Python 3.10.4, 64bit on Win32 (Windows 10)
#  numpy 1.21.6


class Class_1loop_FourTube(object):
    def __init__(self, L1, L2, L3, L4, A1, A2, A3, A4, rg0=0.95, rl0=0.9 ,sampling_rate=48000):
        # initalize Tube length and Tube area
        self.L1= L1 # set list of 1st tube's length by unit is [cm]
        self.A1= A1 # set list of 1st tube's area by unit is [cm^2]
        self.L2= L2 # set list of 2nd tube's length by unit is [cm]
        self.A2= A2 # set list of 2nd tube's area by unit is [cm^2]
        self.L3= L3 # set list of 3rd tube's length by unit is [cm]
        self.A3= A3 # set list of 3rd tube's area by unit is [cm^2]
        self.L4= L4 # set list of 4th tube's length by unit is [cm]
        self.A4= A4 # set list of 4th tube's area by unit is [cm^2]
        C0=35000.0  # speed of sound in air, round 35000 cm/second
        self.sr= sampling_rate
        self.tu1=self.L1 / C0   # delay time in 1st tube
        self.tu2=self.L2 / C0   # delay time in 2nd tube
        self.tu3=self.L3 / C0   # delay time in 3rd tube
        self.tu4=self.L4 / C0   # delay time in 4th tube
        
        # loop-in portion
        self.r12=((self.A3+self.A2) - self.A1)/(self.A3+self.A2+self.A1)  # reflection coefficient between 1st tube and others
        self.r21=(self.A2 - (self.A3+self.A1))/(self.A3+self.A2+self.A1)  # reflection coefficient between 2nd tube and others
        self.r31=(self.A3 - (self.A2+self.A1))/(self.A3+self.A2+self.A1)  # reflection coefficient between 3rd tube and others
        # loop-out portion
        self.r23=((self.A4+self.A3) - self.A2)/(self.A4+self.A3+self.A2)  # reflection coefficient between 2nd tube and others
        self.r32=((self.A4+self.A2) - self.A3)/(self.A4+self.A3+self.A2)  # reflection coefficient between 3rd tube and others
        self.r42=(self.A4 - (self.A3+self.A2))/(self.A4+self.A3+self.A2)  # reflection coefficient between 4th tube and others
        
        
        self.rg0=rg0 # rg is reflection coefficient between glottis and 1st tube
        self.rl0=rl0 # reflection coefficient between 4th tube and mouth
        
        self.att_norm=1.0 # attenuation constant per one step ahead  out of loop
        self.att_loop=0.998 #1.0 # attenuation constant per one step ahead  in loop
        
        self.num_of_tube=4
        
    def fone(self, xw):
        pass
        
    def H0(self, freq_low=100, freq_high=5000, Band_num=256):
        pass
        
    def process(self, yg ):
        # process reflection transmission of resonance tube: yg is input, y2tm is output
        # one-loop consists of tube2 and tube3
        #
        #
        #             -----------------------------------------------------------------------------------------
        # yg->  ya1-> |                       yb1->                                     yd1->                  |->    y2tm
        #          rg |     tube1                          tube2                                  tube4        | rl0
        #             |                                                         <-yb2                          |<- yd2
        #             |                            -----------------------------             ------------------|
        #           <-|                   <-ya2                                              |
        #              -------------------                                                   |
        #                                |     yc1->       tube3                             |
        #                                |                                      <-yc2        |
        #                                -----------------------------------------------------
        #                               
        #
        #
        M1= round( self.tu1 * self.sr ) + 1  # for precision, higher sampling_rate is better
        M2= round( self.tu2 * self.sr ) + 1  # for precision, higher sampling_rate is better
        M3= round( self.tu3 * self.sr ) + 1  # for precision, higher sampling_rate is better
        M4= round( self.tu4 * self.sr ) + 1  # for precision, higher sampling_rate is better
        M1= int(M1)
        M2= int(M2)
        M3= int(M3)
        M4= int(M4)
        ya1=np.zeros(M1)
        ya2=np.zeros(M1)
        yb1=np.zeros(M2)
        yb2=np.zeros(M2)
        yc1=np.zeros(M3)
        yc2=np.zeros(M3)
        yd1=np.zeros(M4)
        yd2=np.zeros(M4)
        y2tm=np.zeros(len(yg))
        
        for tc0 in range(len(yg)):
            for i in range((M1-1),0,-1): # process one step
                ya1[i]=ya1[i-1] * self.att_norm
                ya2[i]=ya2[i-1] * self.att_norm
            for i in range((M2-1),0,-1): # process one step
                yb1[i]=yb1[i-1] * self.att_loop
                yb2[i]=yb2[i-1] * self.att_loop
            for i in range((M3-1),0,-1): # process one step
                yc1[i]=yc1[i-1] * self.att_loop
                yc2[i]=yc2[i-1] * self.att_loop
            for i in range((M4-1),0,-1): # process one step
                yd1[i]=yd1[i-1] * self.att_norm
                yd2[i]=yd2[i-1] * self.att_norm
            # calculate reflection
            #  tube1 input
            ya1[0]= ((1. + self.rg0 ) / 2.) * yg[tc0] + self.rg0 * ya2[-1]
            ya2[0]= -1. * self.r12 * ya1[-1] + ( 1. - self.r12 ) * yb2[-1] + ( 1.  - self.r12) * yc2[-1]
            #  tube2  loop 1/2
            yb1[0]= ( 1 + self.r21 ) * ya1[-1] + self.r21 * yb2[-1] + ( 1. + self.r21) * yc2[-1]
            yb2[0]=  -1. * self.r23  * yb1[-1] + ( 1. - self.r23 ) * yc1[-1] + ( 1.  - self.r23) * yd2[-1]
            #  tube3  loop 2/2
            yc1[0]= ( 1. + self.r31 ) * ya1[-1] + ( 1. + self.r31 ) * yb2[-1] + self.r31 * yc2[-1]
            yc2[0]= ( 1. - self.r32 ) * yb1[-1] +  -1. * self.r32  * yc1[-1] + ( 1.  - self.r32) * yd2[-1]
            #  tube 4 output
            yd1[0]=  ( 1. + self.r42 ) * yb1[-1] + ( 1. + self.r42 ) * yc1[-1] + self.r42 * yd2[-1]
            yd2[0]=  -1. * self.rl0  * yd1[-1]
            y2tm[tc0]= (1 + self.rl0) * yd1[-1]
            
        return y2tm


    def check1(self,):
        # check if sum of output, ya1, yb2, and yc2 is 1
        a1= -1. * self.r12  +  ( 1 + self.r21 ) + ( 1. + self.r31 )
        b2= ( 1. - self.r12 )   +  self.r21 + ( 1. + self.r31 )
        c2= ( 1. - self.r12 )   +  ( 1. + self.r21)  +  self.r31 
        print('check1',a1,b2,c2)
        # check if sum of output, yb1, yc1, and yd2 is 1
        b1= -1. * self.r23 + ( 1. - self.r32 ) +  ( 1. + self.r42 )
        c1= ( 1. - self.r23 ) +  -1. * self.r32  + ( 1. + self.r42 )
        d2= ( 1.  - self.r23) + ( 1.  - self.r32) + self.r42
        print('check2',b1,c1,d2)




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
    
    #
    L4= 1.0 # L1_a
    A4= A1_a
    
    # insatnce
    tube  =  Class_1loop_FourTube(L1_a, L2_a, L3_o, L4, A1_a, A2_a, A3_o, A4, sampling_rate=48000)
    
    tube.check1()
    
    
    