 # -*- coding: utf-8 -*-
"""
Created on Tue Sep 15 11:57:26 2015

@author: jordan
"""
import os
import Features as f
import pandas as pd
from time import time, strftime
from inspyred import ec


class Ga(object):

    def __init__(self, cell_name, outfile_dir):
        '''
        what_to_optimize = list of conductances or optimized parameters
        bounds = constraints on each parameter for solution generation
            can check bounds by typing self.bounder.upper_bounds (or lower)
        '''
        self.cell = cell_name
        
        self.myfile = outfile_dir + strftime("GA_%Y_%b_%d") + ".h5"
        
        self.gen_number = 0
        self.n = 0
        self.optimize = []

        #Create or rename our hdf file
        while os.path.exists(self.myfile):
            self.n += 1
            self.myfile = outfile_dir+strftime("GA_%Y_%b_%d") + "_" + str(self.n) + ".h5"
        
#        init = self.myfile.rstrip('.h5')
#        index = init.lstrip(outfile_dir).rstrip('.h5')
#        df = pd.DataFrame(self.cell._initials)
#        print df
#        df.to_csv(outfile_dir+"GA_Initials.csv")

        self.hd5 = pd.HDFStore(self.myfile, 'a')

        self._features = {}

        #Get boundaries for inpsyred's bounder & generator functions
        bounds=[]
        for key, val in self.cell._optimize.items():
            for key2, val2 in val.items():
                self.optimize.append(key+'_'+key2)
                bounds.append(val2)
        self._lb, self._ub = [list(x) for x in zip(*bounds)]
        #Set inspyred's bounder function with lower and upper bounds
        self.bounder = ec.Bounder(self._lb, self._ub)

        print "\nOptimizing the following:\n",  '\n'.join(self.optimize)
    def __call__(self, *args, **kwargs):
        candidate = [a for a in args]
        fit = self.evaluator([candidate], kwargs)
        return fit[0]

    def generator(self, random, args):
        #Function to generate solutions

        return [random.uniform(self._lb[x],self._ub[x])
                                    for x in range(len(self._lb))]

    def membership(self, lower, upper, value):

        _range = upper-lower
        midpoint = lower + (_range / 2.0)
        pt_a = lower - (((midpoint - lower) / 4.0) * 3.0)
        pt_b = lower + (((midpoint - lower) / 4.0) * 1.0)
        pt_c = upper - (((upper - midpoint) / 4.0) * 1.0)
        pt_d = upper + (((upper - midpoint) / 4.0) * 3.0)

        if value >= pt_b and value <= pt_c:
            nMembership = 1.0

        elif value > pt_a and value < pt_b:
            nMembership = (value - pt_a) / (pt_b - pt_a)

        elif value > pt_c and value < pt_d:

            nMembership = (pt_d - value) / (pt_d - pt_c)
        else:
            nMembership = 0.0

        return nMembership

    def analyze(self, time, v_soma, v_axon):

        mem = []
        an = f.Analyze_Features(time, v_soma, v_axon)

        #Objective 1: Period 800.0, 1200.0
        mem.append(self.membership(800., 1200., an.period()))

        #objective 2: BurstDuration 100.0, 225.0
        mem.append(self.membership(100., 225., an.burst_duration()))

        #objective 3: ISI 20.0, 45.0
        mem.append(self.membership(20., 45., an.av_isi()))

        #objective 4: Spikeamp 3.0, 8.0
        mem.append(self.membership(3., 8., an.spike_amp()))

        #objective 5: AHP -65., -50
        mem.append(self.membership(-65., -50., an.ahp()))

        #objective 6: Slow wave amp
        mem.append(self.membership(10., 25., an.slow_wave_amp()))

        #objective 6: Time to first spike
        #mem.append(self.membership(#,#, an.time_to_first_spike()))

        #objective 8: Voltage Range
        #mem.append(self.membership(#,#, an.vrange))

        self._features = an.results

        return mem

    def sim_protocol(self):
#
#        self.cell.zero_synapses(group='ne')
#        self.cell.iclamp(sec='soma')
#        self.cell.add_vector(
#                'i', self.cell.clamp[str(self.cell._clampcount-1)]._ref_i)
        results = self.cell.run()
        return results

    def evaluator(self, candidates, args):
        #Evaluate fitness
        fitness = []

        params = pd.DataFrame(columns=self.optimize)
        fit_val = pd.DataFrame()
        analysis = pd.DataFrame()
        s = time()
        #print 'Evaluating generation %d (values : fitness)' %self.gen_number

        #iterate over all candidate solutions
        for cs in candidates:
            #join 'what_to_optimize' with a single candidate solution
            opt = dict(zip(self.optimize, cs))
            #append fitness for internal GA functions#append current solution dictionary to the pandas dataframe
            params = params.append(opt, ignore_index=True)
            #update the parameters to the new candidate values/percentage
            for x in self.optimize:
                p = x.partition('_')
                self.cell._optimize[p[0]][p[2]] = opt[x]

            self.cell.ga_change_gbar(self.cell._optimize)

            #run the sim
            d = self.sim_protocol()

            #analyze their fitness
            fit = self.analyze(d['time'], d['vs'], d['va'])

            #append fitness and feature output to pandas dataframes
            fzip = dict(zip(self._features.keys(),fit))
            fit_val = fit_val.append(fzip, ignore_index=True)
            fit_val['total_fit'] = fit_val.sum(axis=1)
            analysis = analysis.append(self._features, ignore_index=True)

            #append fitness for GA
            fitness.append(ec.emo.Pareto(fit))
            print ', '.join(map(str, cs)) + ' :\nFit ' + ', '.join(map(str,fit))
        print '\nCompleted generation %d in %d seconds\n' %(
                self.gen_number,(time()-s))

        self._hd5_write('/gbars', params)
        self._hd5_write('/fitness', fit_val)
        self._hd5_write('/feature', analysis)

        self.gen_number += 1
        self.cell._clampcount = 0
        return fitness

    #def observer(self, candidates, args):

    def _hd5_write(self, folder, data):
        fname = 'gen'+str(self.gen_number)
        return data.to_hdf(self.hd5,
                           fname+folder,
                           format='table',
                           data_columns=True)
