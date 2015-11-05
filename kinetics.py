# -*- coding: utf-8 -*-
"""
Created on Sat Oct 31 19:56:47 2015

@author: jordan
"""
import re
import numpy as np
import sympy.parsing.sympy_parser as sp
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
from sympy import lambdify, symbols

class Kinetics(object):

    def __init__(self,filename):
        self.eqs = {}
        self.vars = ['v']
        self.initials = []
        self._slider_vals = []
        self._sliders = []
        self._started = 0
        self._args = {}
        self._lams = {}
        self._keys = []
        self._vrange = np.arange(-100,100,1)
        self._button = None
        
        with open(filename,'r') as f:
            for line in f:

                if "}" in line:
                    self._started=0

                elif line.startswith("PARAMETER"):
                    self._started = 1

                elif line.startswith("PROCEDURE"):
                    self._started = 2

                elif (self._started==1 and not line.startswith(":")
                        and not "{" in line):

                    if "scale" in line or "shift" in line:
                        x = line.strip(' \n')
                        x = x.partition('=')
                        self.vars.append(x[0].strip())
                        y = x[2].rstrip(' (ms)\n')
                        y = x[2].rstrip(' (mV)\n')
                        self.initials.append(float(y))
                elif (self._started==2 and not "{" in line
                        and not "TABLE" in line and not "FROM" in line):
                      x = line.strip(' \n')
                      x = x.partition('=')
                      y = re.sub('\(ms\)','',x[2])
                      y = re.sub('\(mv\)','',y)
                      args = [i for i in self.vars if i in y]
                      sympy_exr = sp.parse_expr(y,
                                  transformations=sp.standard_transformations)
                      self._args[x[0]] = args[1:]
                      self._keys.append(x[0])
                      self._lams[x[0]] = [symbols(i) for i in args]
                      self.eqs[x[0]] =  lambdify(self._lams[x[0]], sympy_exr, 'numpy')

        self._create_canvas()
        self._get_sliders()
        self._dynamic_plot()
        for slider in self._slider_vals:
            slider.on_changed(self._dynamic_plot)
        plt.show()
        self._button.on_clicked(self._reset)
    def _create_canvas(self):

        if len(self._keys)>3 and len(self._keys)%2==0:
            self.fig, self.axes = plt.subplots(
                                    len(self._keys)/2,len(self._keys)/2)

            for i, ax in enumerate(self.axes.flatten()):
                ax.set_title(self._keys[i])

        elif len(self._keys)>3 and len(self._keys)%2==1:
            fig, axes = plt.subplots(len(self._keys),1)

            for i, ax in enumerate(self.axes.flatten()):
                ax.set_title(self._keys[i])
        else:
            fig, axes = plt.subplots(len(self._keys))
            for i, ax in enumerate(self.axes):
                ax.set_title(self._keys[i])
        
        self.fig.canvas.set_window_title( "Kinetics")
    
        
    def _dynamic_plot(self, *args):

        n, idx = 0, 0

        step = (len(self._sliders)/len(self._keys))

        for ax in self.axes.flatten():
            ax.clear()
            ax.set_title(self._keys[n])
            ax.plot(self._vrange, self.eqs[self._keys[n]](self._vrange,
                    self._slider_vals[idx].val, self._slider_vals[idx+1].val)) 
            
            self.fig.canvas.draw_idle()

            n+=1
            idx+=step


    def _get_sliders(self):
        #get slider properties
        n = 0
        for eq in self._keys:
            for var in self._args[eq]:
                if 'v' != var:
                    self._sliders.append({'label': var,
                                         'valmin': -80.,
                                         'valmax': 80.,
                                         'valinit': self.initials[n]})
                    n+=1
        #Use properties to build the sliders                
        nVars = len(self._sliders)
        
        #Create canvas
        fig = plt.figure()
        fig.canvas.set_window_title( "Kinetic Inputs")

        #CREATE THE SLIDERS
        self._slider_vals = []
        for i, properties in enumerate(self._sliders):
            ax = plt.subplot2grid((nVars+1, 3), (i, 0), colspan=3)
            self._slider_vals.append(Slider(ax=ax, **properties))
        button_ax = plt.subplot2grid((nVars+1, 3),(nVars, 1))
        self._button = Button(button_ax, 'Reset', hovercolor='0.975')
    
    def _reset(self,*args):
        for s in self._slider_vals:
            s.reset()