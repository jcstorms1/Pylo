# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 13:00:54 2015

@author: jordan
"""
from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
#from scipy.interpolate import interp1d

class Analysis(object):
    
    def __init__(self, time, voltage, current):    
        self.t = time
        self.v = voltage
        self.c = current
        self.indexes = []
        
        self.burst_start=[]
        self.burst_end = []
        self.bursts = []                
        self._single_burst = []
        self.isi = []
        
    def av_isi(self, time=None, isi=300):
        
        self._2nd_to_last_burst(time=time, isi=isi)
    
        return np.mean(self.isi[self._single_burst])
        
    
    def burstduration(self, time=None, voltage=None, isi=300):
                    
        try:
            if not self.isi:
                 self._get_isi(time=time)
                 self._bursts(isi=isi)
        except:
            self.isi, self.burst_end, self.burst_start
        
        return self.t[self.burst_end[-1]]-self.t[self.burst_start[-2]]
              

    def dv_v(self, time=None, voltage=None, plot=False, **kwargs):
        if time is None:
            time = self.t
        if voltage is None:
            voltage = self.v
        
        dvdt = np.diff(voltage)/np.diff(time)
        V = voltage[:-1]
        
        if plot:
            plt.figure()
            plt.plot(V, dvdt, 'og')
            plt.title('Phase-Plane Density')
            plt.xlabel('V (mV)')
            plt.ylabel('dV/dt (V/s)')
            plt.show()
        return dvdt, V
    
    def find_peaks(self, voltage = None, thresh=-45, dist=1):
        if voltage is None:
            voltage = self.v
        #take first order difference of voltage array 
        diff = np.diff(voltage)
        #Indexes that are < 0 when adding a 0 to the end of diff, 
        # > 0 when adding a 0 to the beginning of diff, and > thresh are peaks
        indexes = np.where((np.hstack([diff, 0.])< 0.) 
                            & (np.hstack([0., diff]) >0.) 
                            & (voltage > thresh))[0]    
        #sort for finding only max peaks given a distance (isi)              
        if dist > 1 and len(indexes) > 1:
            #resort the indexes so that their voltage value is from min to max
            order = indexes[np.argsort(voltage[indexes])]
            #create a boolean array the same length os voltage (will be all True)
            b = np.ones_like(voltage, dtype=bool)
            #change the indexes to False
            b[indexes] = False        
            #iterate over each ordered index
            for index in order:
                #create a slice that is the index + or - dist
                s = np.s_[max(0, index-dist):index+dist+1]
                #all indexes in the slice are now True
                # Note: indexes with smaller voltages will be overwritten with True
                b[s] = True
                #only the index is false
                b[index] = False
            #peaks are the indexes that correspond to False in the boolean array
            indexes = np.where(b==False)[0]
        self.indexes = indexes
    
    def fi_curve(self, amp, stim_start, stim_end, time=None, voltage=None):  
        if time is None:
            time = self.t
        if voltage is None:
            voltage = self.v
        
        t, v, s0, s1, = time, voltage, stim_start, stim_end
        #find array index of start and end time of current pulse
        t0 = int(np.where(np.isclose(t,s0,rtol=1e-8))[0])
        t1 = int(np.where(np.isclose(t,s0+s1,rtol=1e-8))[0])
        #slice out the portion between t0 and t1
        v_slice = v[t0:t1]
        
        index = self.find_peaks(v_slice)
        
        freq = (len(index)/(((s0+s1)-s0)/1000))
        
        #amp = [v.max() if v.max() > 0 else v.min()]
        
        plt.plot(amp, freq,'o-')
        plt.title('FI Curve')
        plt.xlabel('Amplitude (nA)')
        plt.ylabel('Frequency (Hz)')
        plt.show
        
    def impedance_cycle(self, time, voltage, current, window_len=11):
        if time is None:
            time = self.t
        if voltage is None:
            voltage = self.v
            
        zc = np.where(np.sign(current[1:]) != np.sign(current[:-1]))[0][::2]
        indexes = np.column_stack((zc[0:-1], zc[1:]))
        
        V = np.array([np.ptp(voltage[x[0]:x[1]]) for x in indexes])
        I = np.array([np.ptp(current[x[0]:x[1]]) for x in indexes])
        freq = np.array([1000./np.diff(time[x]) for x in indexes])
    
        Z = np.abs(V/I)
          
        fit = self.smooth(Z, window_len)
        return freq, Z, fit
        
    def impedance_fft(self, time, voltage, current):
        if time is None:
            time = self.t
        if voltage is None:
            voltage = self.v
            
        v = voltage - np.mean(voltage)
        c = current - np.mean(current)    
        dt = np.diff(time)[0]    
        L = len(time)
        T = int(L*dt)
        m = 2**int(np.ceil(np.log2(L)))
        sample = int((L/T)*1e3)
        f = (sample/2*(np.linspace(0,1,m/2)))[1:1000]
        
        V = np.fft.fft(v, n=m)/L
        I = np.fft.fft(c, n=m)/L
        Z = (V/I)
        Z = abs(Z[0:m])[1:1000]
    
        fit = self.smooth(Z, window_len=len(Z)/10)
            
        return f, Z, fit
    
    def mean_cur(self, time=None, voltage=None, current=None, plot=False):
        if time is None:
            time = self.t
        if voltage is None:
            voltage = self.v
        if current is None:
            current = self.c
                    
#        if not self.indexes:
#            self.indexes = self.find_peaks(voltage)
        
        a, b = self.indexes[0:-1], self.indexes[1:]
        
        t_delta = np.diff(time[np.column_stack((a,b))])
        hz = 1000*(1/t_delta)
    
        i_mean = []
        
        for j, k in zip(a, b):
            i_mean.append(np.mean(current[j:k]))
        
        if plot:
            plt.plot(i_mean, hz, 'bo')
            plt.xlabel('Mean Current (nA)')
            plt.ylabel('Instantaneous Frequency (Hz)')
            plt.title('Current Analysis')
        
        return hz, i_mean
    
    def period(self, time=None):
        "returns time difference between last two spikes"
               
        if not self.bursts:
            self._get_isi()
            self._bursts()
        
        return float(self.burst_end[-2]-self.burst_end[-3])
        
    def resample(self, time, original_array, *args):         
        L = len(time)
        num = 2**int(np.ceil(np.log2(L)))
        newtime = np.linspace(time[0], time[-1], num)
        resampled = np.zeros_like(newtime)  
        
        f = interp1d(time, original_array)
        resampled = f(newtime)
           
        if args:
            for arg in args:
                f = interp1d(time, arg)
                resampled = np.column_stack((resampled, f(newtime)))
                
        return newtime, resampled.T
    
    def smooth(self, v, window_len=21):
        if v.ndim != 1:
            raise ValueError, "smooth only accepts 1 dimension arrays."
        elif v.size < window_len:
            raise ValueError, "Input vector needs to be bigger than window size."
    
        if window_len<3:
            return v
    
        slices = np.r_[v[window_len-1:0:-1],v,v[-1:-window_len:-1]]
    
        window = np.hamming(window_len)
    
        y = np.convolve(window/np.sum(window),slices,mode='valid')
    
        index = window_len/2
    
        return y[index:-index]
    
    def spikes_per_burst(self, time=None):
        
        if not self.bursts:
            self._get_isi()
            self._bursts()
        
        self._2nd_to_last_burst()
        
        return self._single_burst.size
            
    def _get_isi(self, time=None):
        #
        if time is None:
            time = self.t
        
        try:
            if not self.indexes:
                self.find_peaks()
        except:
            self.indexes
        
        pad = np.hstack([0, self.indexes])
        self.isi = np.diff(time[pad])
        
    def _bursts(self, time=None, isi=300):
        if time is None:
            time = self.t
        
        if not self.indexes:
            self.find_peaks()
               
        if not self.isi:
            self._get_isi()
                        
        s = np.where(self.isi >= isi)[0]
        e = s-1
        self.burst_start = self.indexes[s]
        self.burst_end = self.indexes[e]
        self.bursts = np.column_stack((self.burst_start, self.burst_end))          
    
    def _2nd_to_last_burst(self, time=None, isi=300):
        if time is None:
            time = self.t
        
        if not self.isi:
           self._get_isi()
           self._bursts()
                
            
        self._single_burst = np.where((self.burst_start[-2] <= self.indexes) &
                            (self.indexes <= self.burst_end[-1]))[0]
        
        

        
        
    
        
    #from scipy.stats import binning_statistic_2d as bstat
    #d, v = Dv_v(t, soma)
    #a, b, c, num = bstat(v, d, bins=100, statistic='mean')
    #x = [~np.isnan(a)]
    #rmse = np.sqrt(((x-x1)**2).mean())
    
    #def impedance_cycle(time, voltage, current, delta_freq, 
    #                                    freq_step=0.1, window_len=11):
    #    '''Returns: freq, Z, fit'''
    #    #Takes the frequency recorded from the zap mechanism
    #    t, v, c, f = time, voltage, current, delta_freq
    #    #Create array that uses a freq step to find the end of a period by each
    #    #step size change in frequency (e.g. 0.1 to 4.0 with a 0.1 step size)
    #    ints = np.arange(f[0], np.ceil(np.max(f)), freq_step)
    #    sort = np.searchsorted(f,ints)
    #    indexes = np.column_stack((sort[0:-1], sort[1:]))
    #        
    #    V = np.array([np.ptp(v[x[0]:x[1]]) for x in indexes])
    #    I = np.array([np.ptp(c[x[0]:x[1]]) for x in indexes])
    #    freq = np.array([1000./np.diff(t[x]) for x in indexes])
    #
    #    Z = np.abs(V/I)
    #      
    #    fit = smooth(Z, window_len)
    #    return freq, Z, fit