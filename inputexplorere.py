# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 14:02:38 2015

@author: jordan
"""

import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
 
def inputExplorer(f, sliders_properties, wait_for_validation = False):
    """ A light GUI to manually explore and tune the outputs of 
        a function.
        slider_properties is a list of dicts (arguments for Slider )
        whose keys are in ( label, valmin, valmax, valinit=0.5, 
        valfmt='%1.2f', closedmin=True, closedmax=True, slidermin=None,
        slidermax=None, dragging=True)
         
        def volume(x,y,z):
            return x*y*z
     
        intervals = [ { 'label' :  'width',  'valmin': 1 , 'valmax': 5 },
                  { 'label' :  'height',  'valmin': 1 , 'valmax': 5 },
                  { 'label' :  'depth',  'valmin': 1 , 'valmax': 5 } ]
        inputExplorer(volume,intervals)
    """
         
    nVars = len(sliders_properties)
       
    # CREATE THE CANVAS
     
    figure = plt.figure()
    figure.canvas.set_window_title( "Inputs for '%s'"%(f.func_name) )
     
    
    # CREATE THE SLIDERS
    sliders = []
     
    for i, properties in enumerate(sliders_properties):
        ax = plt.subplot2grid((nVars,1),(i, 0))
        sliders.append( Slider(ax=ax, **properties) )
     
     
    # CREATE THE CALLBACK FUNCTIONS
     
    def on_changed(event) : 
         
        res = f(*(s.val for s in sliders))
        print [s.val for s in sliders]
        if res is not None:
             
            print res
     
    def on_key_press(event):
         
        if event.key is 'enter':
             
            on_changed(event)   
     
    figure.canvas.mpl_connect('key_press_event', on_key_press)
     
    # AUTOMATIC UPDATE ?
     
    if not wait_for_validation:
         
        for s in sliders :
             
            s.on_changed(on_changed)
     
     
    # DISPLAY THE SLIDERS
     
    plt.show()