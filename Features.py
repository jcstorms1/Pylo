# -*- coding: utf-8 -*-
"""
Created on Tue Sep 29 14:27:07 2015

@author: jordan
"""

from __future__ import division
import numpy as np

#from scipy.interpolate import interp1d

class Analyze_Features(object):
    """
    Analysis class designed for integration with the GA module.

    Only uses the last 3/4 of the time, soma, and axon arrays
    (can be changed in __init__). All results are stored in the self.results
    dict as the average of a specific feature for the entire array.

    Currently the following methods are available for feature extraction:
            1. **ahp**                      *(class method)*
            2. **av_isi**                   *(class method)*
            3. **burst_duration**           *(class method)*
            4. **period**                   *(class method)*
            5. **slow_wave_amp**            *(class method)*
            6. **spike_amp**                *(class method)*
            7. **time_to_first_spike**      *(class method)*
            8. **vrange**                   *(class variable)*

    Other semi-private variables store indexing information, and in some cases,
    contain the value of each occurence (`e.g.` self._ahp_list (*indexes*) and
    self._ahp_voltage (*actual ahp voltage*)). However, most of these variables
    will only appear after a method is called.

    Each method is decorated with the _check_std method. The first argument is
    the feature name, which is used as the key for the self.results dictionary.
    The second argument is an estimate of the maximum STD possible for that
    specific feature in a good model. As such, each method being used for the
    GA should have the threshold argument set independently in the @_check_std
    decorator, above that specific feature's method.

    The @_check_std decorator allows for additional control over the analysis,
    by discarding features that do not lie within one standerd deviation from
    the mean. The _check_std method can be fine tuned via the threshold to
    allow for cases where a small STD may be calculated. Without the threshold
    set properly, good models may be rejected, because a good feature will be
    is considered an outlier. Too large of a threshold will accept bad models
    as good models, and too low of a threshold may reject good models.

    Notes:
    ------
    The `self._flag` (*Boolean*) stops analysis when certain methods return
    empty arrays or values. The default is `self._flag = False`.

    Attributes:
    ----------
    time : array
        Time array of a simulation
    v_soma : array
        Voltage array of soma
    v_axon : array
        Voltage array of axon
    max_isi : int
        Maximum isi allowed for spike detection.
        **Caution:** Too large of an isi will interfere with burst and
        intraburst detection.
    """

    def __init__(self, time, v_soma, v_axon, max_isi=300,
                 max_thresh=-45, max_dist=1, min_dist=1):
        #use only last 3/4 of arrays
        qtr = np.abs(time-np.max(time*0.25)).argmin()
        self.t = time[qtr:]
        self.v = v_soma[qtr:]
        self.a = v_axon[qtr:]
        self.max_isi = max_isi
        #get maximum (soma & axon), minimum points, and vrange
        self.s_max = self.find_maxima(self.v, thresh=max_thresh, dist=max_dist)
        self.a_max= self.find_maxima(self.a, thresh=0.0, dist=100)
        self.minima = self.find_minima(dist=min_dist)
        self.vrange = np.ptp(self.v)
        #get somatic spike indexes based off axon peaks
        try:
            self.spikes = self.s_max[np.searchsorted(
                               self.t[self.s_max], self.t[self.a_max],'right')]
        except:
            self.spikes = np.array([])

        self.results = {}
        #if we have a bad model skip analysis and set the result to 0.0
        self._flag = True #Good models keep this True bad

    def find_maxima(self, arr, thresh=-80, dist=1):
        #take first order difference of voltage array
        diff = np.diff(arr)
        #Indexes that are < 0 when adding a 0 to the end of diff,
        # > 0 when adding a 0 to the beginning of diff, and > thresh are peaks
        indexes = np.where((np.hstack([diff, 0.])< 0.)
                            & (np.hstack([0., diff]) >0.)
                            & (arr > thresh))[0]
        #sort for finding only max peaks given a distance (isi)
        if dist > 1 and len(indexes) > 1:
            #resort the indexes so that their self.v value is from min to max
            order = indexes[np.argsort(arr[indexes])]
            #create a boolean array the same length os voltage (will be all True)
            b = np.ones_like(arr, dtype=bool)
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

        return indexes

    def find_minima(self, dist=1):
        """
        Same as find_maxima, but we flip the signs of the `indexes=` statement
        and we don't use a threshold. Creates a list of AHP's as well.
        """
        #take first order difference of voltage array
        diff = np.diff(self.v)
        #Indexes that are < 0 when adding a 0 to the end of diff,
        # > 0 when adding a 0 to the beginning of diff, and > thresh are peaks
        indexes = np.where((np.hstack([diff, 0.])> 0.)
                            & (np.hstack([0., diff]) <0.))[0]
        #sort for finding only max peaks given a distance (isi)
        if dist > 1 and len(indexes) > 1:
            #resort the indexes so that their voltage value is from min to max
            order = indexes[np.argsort(self.v[indexes])]
            #create a boolean array the same length os voltage (will be all True)
            b = np.ones_like(self.v, dtype=bool)
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
            self.minima = np.where(b==False)[0]
        return indexes

    def _check_std(name,std_thresh):
        #wrapper to check all elements of a feature are within 1 standard dev.
        def wrap(f):
            def inner(self):
                x = f(self)
                p_std = x.mean() + x.std()
                m_std = x.mean() - x.std()
                if x.std() < std_thresh:
                    self.results[name] = x.mean()
                elif m_std <= x.all() <= p_std:
                    self.results[name] = x.mean()
                else:
                    self.results[name] = 0.0
                return self.results[name]
            return inner
        return wrap

    @_check_std('ahp', 5)
    def ahp(self):

        try:
            self._ahp_voltage
        except:
            self._get_ahps()


        if self._flag:
            return self._ahp_voltage

        else:
            return np.array([])

    @_check_std('avisis', 15)
    def av_isi(self):
        #check we have the intraburst list
        try:
            self._intrabursts
        #if not get it
        except:
            self._get_intraburst()

        if self._flag:
            try:
                #get the index of the intrabursts in the isi list to delete them
                intra_idx = [np.where(np.isclose(self._isi_list, x))[0]
                                   for x in self._intrabursts.ravel()]

                #delete the first isi and the intrabursts in the isi_list
                isi = np.delete(self._isi_list, intra_idx)[1:]

                return  isi

            except:
                return np.array([])

        else:
            return np.array([])

    @_check_std('burstduration', 25)
    def burst_duration(self):
        #check we have the burst times
        try:
            self.bursts
        #if not call the _get_bursts method
        except:
            self._get_bursts()

        if self._flag:
            #get all the burst lengths and flatten into a single array
            self.burst_lengths = np.diff(self.t[self.bursts]).ravel()
            return self.burst_lengths

        else:
            return np.array([])

    @_check_std('period', 50)
    def period(self):
        "returns time difference between last two spikes"
        #check that we have an array of each cycle time
        try:
            self.cycles
        #if not call the appropriate methods
        except:
            self._get_bursts()

        if self._flag:
            return self.cycles

        else:
            return np.array([0])

    @_check_std('slowwaveamp',1)
    def slow_wave_amp(self):
        try:
            self._spike_mins
        except:
            self._get_spike_mins()

        try:
            self._ahp_voltage
        except:
            self._get_ahps()

        if self._flag:
            try:
                #take the min of all intraspike minima
                mean_spike_min = np.mean(self.v[self._spike_mins],axis=1)
                #make sure we have the same amount of mean spike mins and ahps
                #take the difference of the voltage
                if self._last_spike[-1] > self._ahp_list[-1]:
                    swa = np.abs(mean_spike_min[:-1] - self.v[self._ahp_list])
                elif self._last_spike[0] < self._ahp_list[0]:
                    swa = np.abs(mean_spike_min - self.v[self._ahp_list[1:]])
                else:
                    swa = np.abs(mean_spike_min - self.v[self._ahp_list])

                self._slow_wave_amps = swa

                return self._slow_wave_amps
            except:
                return np.array([])

        else:
            return np.array([])

    @_check_std('spikeamp', 1)
    def spike_amp(self):

        try:
            self._spike_mins

        except:
            self._get_spike_mins()

        if self._flag:
            try:
                #get all spike indexes minus the last in the burst
                rmv_last_spk = np.delete(self.spikes,
                                np.searchsorted(self.spikes, self._last_spike))

                #stack the indexes for the next part
                peak_to_min = np.column_stack(
                      (self.v[rmv_last_spk], self.v[self._spike_mins.ravel()]))
                #take 1st order diff to get the spike amplitude for every burst
                self._spike_amps = np.abs(np.diff(peak_to_min)).ravel()

                return self._spike_amps
            except:
                return np.array([])
        else:
            return np.array([])
            
    @_check_std('timetofirstspike', 30)
    def time_to_first_spike(self):
        try:
            self._spike_mins
        except:
            self._get_spike_mins()

        try:
            self._ahp_voltage
        except:
            self._get_ahps()

        if self._flag:
            if self._last_spike[-1] > self._ahp_list[-1]:
                tts = np.column_stack((self._last_spike[:-1], self._ahp_list))
            elif self._last_spike[0] < self._ahp_list[0]:
                tts = np.column_stack((self._last_spike, self._ahp_list[1:]))
            else:
                tts = np.column_stack((self._last_spike, self._ahp_list))

            self._times_to_first_spikes = np.diff(self.t[tts]).ravel()

            return self._times_to_first_spikes

        else:
            return np.array([])

    def _get_ahps(self):
        #check we have our intraburst index list
        try:
            self._intra_indexes
        #check
        except:
            self._get_intraburst()

        if self._flag:

            try:
                #get the indexes of all minimas that fall within each intraburst
                ahp_list = [self.minima[np.where(
                            (self.t[self.minima]>self.t[x[0]]) &
                            (self.t[self.minima]<self.t[x[1]]))[0]].tolist()
                                                 for x in self._intra_indexes]

                #get indexes of the absolute minimums in the first ahp_list
                self._ahp_list = np.array(
                                    [x[np.argmin(self.v[x])] for x in ahp_list])

                #take the minimum voltage if mutliple indexes are in the intraburst
                self._ahp_voltage = self.v[self._ahp_list]
            except:
                self._ahp_list = np.array([])
                self._ahp_voltage = np.array([])

    def _get_isi(self):
        #check that we get an array returned
        try:
            #pad the a_max indexes with zero...
            pad = np.hstack([0, self.a_max])
            #take the first order diff of the time array with those indexes
            self._isi_list = np.diff(self.t[pad]).ravel()
        #if not return an empty array
        except:
            self._isi_list = np.array([])

    #check
    def _get_bursts(self):
        #make sure we have an isi_list
        try:
            self._isi_list
        #if we don't, call the _get_isi method
        except:
            self._get_isi()
        #find the vaues in isi_list that are greater than allowed maximum isi
        #define s as the start of the burst
        s = np.where(self._isi_list >= self.max_isi)[0]
        if s.size:
            #define e as the end of the burst
            e = s-1
            #create start and end arrays of those indexes
            self._burst_start = self.a_max[s]
            #sort since 'e = s-1' gives us the last index as the first element
            self._burst_end = np.sort(self.a_max[e])
            #create stacked array: [burst_start[0], burst_end[0]]...[bs[n],be[n]]
            self.bursts = np.column_stack((self._burst_start, self._burst_end))
            #get the cycle times
            self.cycles = np.diff(self.t[self._burst_start])
        else:
            self._flag=False

    def _get_intraburst(self):
        #makr sure we have our burst list
        try:
            self.bursts
        #if not run _get_bursts method
        except:
            self._get_bursts()

        if self._flag:
            #roll burst array forward one element and cut off first row
            self._intra_indexes = np.roll(self.bursts, 1)[1:,:]
            #get itraburst time difference for using the reorded indexes
            self._intrabursts = np.diff(self.t[self._intra_indexes])

    def _get_spike_mins(self):

        try:
            self._bursts
        except:
            self._get_bursts()

        if self._flag:
            try:
                #get first spike indexes of each burst
                self._first_spike = self.s_max[np.searchsorted(
                                       self.s_max, self._burst_start)]
                #get last spike indexex of each burst
                self._last_spike = self.s_max[np.searchsorted(
                                       self.s_max, self._burst_end)]

                #stack the start and end indexes for the next part...
                start_end = np.column_stack((self._first_spike, self._last_spike))
                #get all minima between the start and end indexes
                self._spike_mins = np.array([self.minima[np.where(
                                    (self.t[self.minima]>self.t[x[0]]) &
                                    (self.t[self.minima]<self.t[x[1]]))[0]]
                                                for x in start_end])
            except:
                self._first_spike = np.array([])
                self._last_spike = np.array([])
                self._spike_mins = np.array([])
