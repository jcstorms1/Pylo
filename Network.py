# -*- coding: utf-8 -*-
"""
Created on Wed Nov  4 14:18:42 2015

@author: jordan
"""
from __future__ import division
from neuron import h
from collections import OrderedDict
import matplotlib.pyplot as plt
import numpy as np
import inspect

class Network(object):
    def __init__(self, *args):
        
        self.vecs = {}
        self._ref = {}
        
        self.cells = []
        
        for arg in args:
            setattr(self, arg.name.lower(), arg)
            
            for key, val in arg._ref.items():
                if key not in self._ref:
                    self.vecs[key] = h.Vector()                    
                    self._ref[key] = val
            
            arg.shunt()
            
    def run(self, v_init=-60, tstop=20000., dt=0.1,
                                        cvode=True, ga_use_half=False):
        '''
        Simulates this cell and all desired vectors are recorded. Uses fixed
        or variable timestep depending on the `cvode=` boolean.

        Parameters:
        ----------
        v_init : int, float
            The starting voltage of the simulation.
        tstop : int, float
            The maximum time of integration
        dt : float
            The desired integration step.
        cvode : bool
            Selects variable time step integration. Default is False.
        ga_use_half : bool
            Will only use the 2nd have of recordings for GA
        '''
        
        h.load_file('stdrun.hoc')
        h.v_init = v_init
        h.tstop = tstop
        h.dt = dt
        #set the recording into the vecs dictionary
        #the _ref dictionary contain the hoc object attribute references
        for key in self._ref.keys():
            #This makes sure we overwrite any vectors if we already ran a sim
            if isinstance(self.vecs['time'], np.ndarray):
                self.vecs[key] = h.Vector()
                self.vecs[key].record(self._ref[key])
            else:
                self.vecs[key].record(self._ref[key])

        
        
        if cvode:
            solver = h.CVode()
            solver.active(1)
            h.finitialize(h.v_init)
            solver.solve(h.tstop)
        else:
            h.CVode().active(0)
            h.finitialize()
            for t in range(0, int(h.tstop/h.dt)):
                h.fadvance()

        for key, val in self.vecs.iteritems():
            self.vecs[key] = np.array(val)

        if ga_use_half:
            for key, val in self.vecs.iteritems():
                self.vecs[key] = val[(val.size)/2:]

        return self.vecs

    def plot(self, *args, **kwargs):
        """
        Function for simple 'time' vs anything plots.

        Takes any user defined vectors except `Time`. To check what was
        recorded type `cell.vecs` in the console. An example of how to use
        this function once a simulation has finished is as follows::

                cell.plot('vsoma')
                cell.plot('vsoma','vaxon','izap')

        """

        n = len(args)

        self.fig, ax = plt.subplots(n,1)
        if 'title' in kwargs:
            self.fig.canvas.set_window_title(kwargs['title'])
            self.fig.suptitle(kwargs['title'], fontsize=11, fontweight='bold')
        if n == 1:
            ax.plot(self.vecs['time'], self.vecs[args[0]])
            ax.set_title('Time vs. ' + args[0].title())

            ax.set_ylabel(args[0].title())
            ax.set_xlim([self.vecs['time'][0], self.vecs['time'][-1]])

        else:
            for i in range(n):
                ax[i].plot(self.vecs['time'], self.vecs[args[i]])
                ax[i].set_title('Time vs. ' + args[i].title())
                ax[i].set_ylabel(args[i].title())
                ax[i].set_xlim([self.vecs['time'][0], self.vecs['time'][-1]])
                if i != (n-1):
                    plt.setp(ax[i].get_xticklabels(), visible=False)
                else:
                    ax[i].set_xlabel('Time')

            plt.tight_layout(h_pad=0.2)
        plt.subplots_adjust(top=0.85)
        plt.show()    