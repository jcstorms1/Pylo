# -*- coding: utf-8 -*-
"""
Created on Tue Jul 14 13:41:16 2015

@author: jordan
"""
from __future__ import division
from neuron import h
from collections import OrderedDict
import matplotlib.pyplot as plt
import numpy as np
import inspect


class Create_Cell(object):
    """
    This class creates a model of a cell using NEURON [H]_.

    Passing a text file following the specified format is parsed using a
    finite-state machine when instantiating the class object. Additional
    user-defined options can be added to the text file if necessary. However,
    the additional state must be added to the __init__ method to correctly
    create them in NEURON.

    Note:
    ----
    Do not include the `self` parameter in the ``Parameters`` section

    Attributes:
    ----------
    filename : str
        Textfile following the appropriate format (see example file)

    Additional Attributes:
    ----------
    tstop : float
        This is the total simulation time and is passed to NEURON's h.tstop.
    compartments : dict
        Used for referencing each compartment. Contains all compartments as
        keys and the associated channels (neuron mechanisms) as values.
    groups : dict
        Used for referencing groups of compartments as defined in `filename`.
    vecs : dict
        Contains all vectors recorded at simulation time. Begins as hoc
        vectors, but are converted to numpy arrays following simulation.

    Examples:
    --------
        >>> #Create a new cell
        >>> myfile = '/mnt/01D09A2178FAB940/Dropbox/Taylor/lpstuff'
        >>> import CreateCells as cc
        >>> lp = cc.Create_Cell(myfile)
        Building a LP cell...
        Creating compartments...
        Inserting channels...
        Connecting compartments...
        Setting paramaters...
        Creating vectors for simulation...
        Cell is completed <('-')> (>'o')> ^(^o^)^ !

        >>> #Check compartments
        >>> lp.compartments
        {'axon': ['leak', 'na', 'kda', 'aa'],
         'cap_near': ['synab', 'synpd', 'synpy'],
         'neur_far': ['synab', 'synpd', 'synpy'],
         'res_near': [],
         'soma': ['leak', 'kd', 'af', 'aas', 'h', 'caint', 'ca', 'kca', 'pr']}

        >>> #Check groups
        >>> lp.groups
        {'all': ['soma', 'cap_near', 'neur_far', 'cap_near', 'neur_far'],
         'ne': ['cap_near', 'neur_far'],
         'sn': ['soma', 'cap_near', 'neur_far']}

        >>> #Check vecs
        >>> lp.vecs
        {'time': <hoc.HocObject at 0x7f334690ac48>,
         'va': <hoc.HocObject at 0x7f334691d420>,
         'vs': <hoc.HocObject at 0x7f334690aae0>}

        >>> #Check conductance and reversal of a channel
        >>> lp.soma.gbar_leak
        0.01564
        >>> lp.soma.e_leak
        -17.96435
        >>> lp.soma.e_leak = -15.
        >>> lp.soma.e_leak
        -15.0

    References:
    ----------
    .. [H] https://www.neuron.yale.edu/neuron/
    """

    def __init__(self, filename):

        self.myfile = filename
        self.compartments = {}
        self.groups = {}
        self.vecs = {}
        self.clamp = {}
        
        #Semi-private attributes
        self._h = h
        self._ref = {}
        self._started = 0
        self._sh = None
        self._clamp = None
        self._clampcount = 0
        self._optimize = OrderedDict()
        self._initials = OrderedDict()
        self._li = lambda l: [x.strip() for x in l.split(',')]

        with open(self.myfile, 'r') as f:
            for line in f:
                if '##' in line:
                    self._started = 0

                elif '#Name' in line:
                    self._started = 1

                elif '#Groups' in line:
                    self._started = 2

                elif '#Morphology' in line:
                    print 'Creating compartments...'
                    self._started = 3

                elif '#Mechanisms' in line:
                    print 'Inserting channels...'
                    self._started = 4

                elif '#Topology' in line:
                    print 'Connecting compartments...'
                    self._started = 5

                elif '#Universals' in line:
                    self._started = 6
		
                elif '#Globals' in line:
                    self._started = 7

                elif '#Parameters' in line:
                    print 'Setting paramaters...'
                    self._started = 8

                elif '#Record' in line:
                    print 'Creating vectors for simulation...'
                    self._started = 9

                elif (self._started == 1 and line.strip()
                                         and not line.startswith('#')):
                    print 'Building a %s cell...' %line.strip()
                    self.name = line.strip()

                elif (self._started == 2 and line.strip()
                                         and not line.startswith('#')):
                    li = self._li(line)
                    self.groups[li[0]] = li[1:]

                #This creates the compartments with the desired dimensions
                elif (self._started == 3 and line.strip()
                                         and not line.startswith('#')):
                    li = self._li(line)
                    setattr(self, li[0], h.Section())
                    self.__dict__[li[0]].L = float(li[1])
                    self.__dict__[li[0]].diam = float(li[2])
                    self.compartments[li[0]] = []

                #This will insert the channels (density mechanisms/mod files)
                elif (self._started == 4 and line.strip()
                                         and not line.startswith('#')):
                    li = self._li(line)
                    if li[0] in self.groups.keys():
                        for sec in self.groups[li[0]]:
                            map(lambda x: self.__dict__[sec].insert(x), li[1:])
                            self.compartments[sec] = li[1:]
                    else:
                        map(lambda x: self.__dict__[li[0]].insert(x), li[1:])
                        self.compartments[li[0]] = li[1:]

                #This will connect the compartments (Topology)
                elif (self._started == 5 and line.strip()
                                         and not line.startswith('#')):
                    li = self._li(line)
                    self.__dict__[li[0]].connect(self.__dict__[li[2]],
                                                    float(li[3]), float(li[1]))
                elif (self._started == 6 and line.strip()
                                         and not line.startswith('#')):
                    li = self._li(line)
                    if 'all' in li:
                        for key in self.compartments.keys():
                            setattr(self.__dict__[key], li[1], float(li[2]))
                    else:
                        setattr(self.__dict[li[0]], li[1], float(li[2]))
                
                #set globals
                elif (self._started == 7 and line.strip()
                                         and not line.startswith('#')):
                    li = self._li(line)
                    setattr(self._h, li[0], float(li[1]))

                #set params
                elif (self._started == 8 and line.strip()
                                         and not line.startswith('#')):
                    li = self._li(line)
                    opt = True if len(li) > 3 else False
                    if li[0] in self.groups.keys():
                        for comp in self.groups[li[0]]:
                            setattr(self.__dict__[comp], li[1], float(li[2]))

                        if opt:
                            try:
                                self._optimize[li[0]][li[1]] = (
                                                    float(li[3]), float(li[4]))
                                self._initials[li[0]][li[1]] = float(li[2])
                            except:
                                self._initials[li[0]] = OrderedDict()
                                self._initials[li[0]][li[1]] = float(li[2])
                                self._optimize[li[0]] = OrderedDict()
                                self._optimize[li[0]][li[1]] = (
                                                    float(li[3]), float(li[4]))
                    elif li[0] == 'h':
                        setattr(self._h, li[1], float(li[2]))
                        if opt:
                            raise AttributeError('h object not set for ga')
                    else:
                        setattr(self.__dict__[li[0]], li[1], float(li[2]))

                        if opt:
                            try:
                                self._optimize[li[0]][li[1]] = (
                                                    float(li[3]), float(li[4]))
                                self._initials[li[0]][li[1]] = float(li[2])
                            except:
                                self._initials[li[0]] = OrderedDict()
                                self._initials[li[0]][li[1]] = float(li[2])
                                self._optimize[li[0]] = OrderedDict()
                                self._optimize[li[0]][li[1]] = (
                                                    float(li[3]), float(li[4]))
                #set vectors for recording
                elif (self._started == 9 and line.strip()
                                         and not line.startswith('#')):
                    li = self._li(line)
                    self.vecs[li[0]] = h.Vector()
                    #check if neuron has that attribute e.g. t
                    if hasattr(self._h, li[-1]):
                        #if True, set _ref dict with the attribute
                        self._ref[li[0]] = getattr(self._h, '_ref_%s' %li[-1])
                    #if False, section has that attribute e.g. v
                    elif hasattr(self.__dict__[li[1]], li[-1]):
                        #if True, set _ref dict with neuron syntax...
                        #e.g. 'soma(0.5)._ref_v'
                        self._ref[li[0]] = getattr(
                         self.__dict__[li[1]](float(li[2])), '_ref_%s' %li[-1])
                    #if both are false raise an exception
                    else:
                        raise AttributeError(
                                        "Can't find the attribute %s!" %li[3])

            print "Cell is completed <('-')> (>'o')> ^(^o^)^ !"

    def add_vector(self, *args):
        """
        This method quickly adds vectors for recording that were not specified
        in the initial file.

        The order of arguments should be::

                cell.add_vector('name', cell.clamp['#']._ref_i)

        After using the `iclamp` or `zap` methods, this function is required
        to add the current to the `cell.vecs` dictionary for recording.

        Examples:
        --------
        A single addition:

        >>> cell.iclamp(sec='soma')
        Created iclamp as clamp 0 in LP.clamp...
        Don't forget to add the vector(s) for recording!!!

        >>> cell.add_vector('i1', cell.clamp['0']._ref_i)
        >>> lp.vecs
        {'i1': <hoc.HocObject at 0x7f57bbb3f9c0>,
         'time': <hoc.HocObject at 0x7f57bbb553d8>,
         'va': <hoc.HocObject at 0x7f57bbb555d0>,
         'vs': <hoc.HocObject at 0x7f57bbb55468>}

        Multiple additions:

        >>> cell.add_vector('i1', cell.clamp['0']._ref_i, 'i2',\
        >>> cell.clamp['1']._ref_i)
        >>> lp.vecs
        {'i1': <hoc.HocObject at 0x7f57bbb3f9c0>,
         'i2': <hoc.HocObject at 0x7f57bbb563d0>,
         'time': <hoc.HocObject at 0x7f57bbb553d8>,
         'va': <hoc.HocObject at 0x7f57bbb555d0>,
         'vs': <hoc.HocObject at 0x7f57bbb55468>}
        """
        args = list(args)
        zargs = zip(args, args[1:])[::2]
        for z in zargs:
            self.vecs[z[0]] = h.Vector()
            self.vecs[z[0]].record(z[1])

    def iclamp(self, amp=-2.0, delay=1000, dur=3000, loc=0.5, sec=None):
        """
        Creates a current clamp NEURON object.

        Note:
        ----
        Loop to simulate multiple step injections

        Parameters:
        ----------
        amp : float
            Amplitude of the current step. Default is -2.0
        delay : int, float
            Time to delay the onset of the current injection. Defalt is 1000.
        dur : int, float
            Duration of the injection. Default is 3000.
        loc : float
            Location of injection in compartment. Default is 0.5
        sec : str
            Compartment where the injection should occur. Default is None

        """
        self._clamp = h.IClamp(loc, sec=self.__dict__[sec])
        self._clamp.amp = amp
        self._clamp.delay = delay
        self._clamp.dur = dur
        self._clamp_setter()

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

        self.shunt()
        
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

    def shunt(self, g=1.4999e-2, e=-35, loc=0.5, sec='soma'):
        """
        An additional leak current in Taylor's model. Possibly to simulate the
        increase in leak caused by the insertation of an electrode?
        """
        self._sh = h.Shunt(loc, sec=self.__dict__[sec])
        self._sh.G = g
        self._sh.e = e

    def zap(self, amp=2, delay=0, dur=None,
                    fmin=0.1, fmax=4.0, loc=0.5, sec='soma'):
        """
        Creates a `Zap` current clamp NEURON object.

        Usage is the same as the `iclamp` method.

        Parameters:
        ----------
        amp : float
            Amplitude of the current step. Default is -2.0
        delay : int, float
            Time to delay the onset of the current injection. Defalt is 1000.
        dur : int, float
            Duration of the injection. Default is 3000.
        loc : float
            Location of injection in compartment. Default is 0.5
        fmin : float
            Minumum zap frequency
        fmax : float
            Maximum zap frequency
        sec : str
            Compartment where the injection should occur. Default is None

        """

        self._clamp = h.zap(loc, self.__dict__[sec])
        self._clamp.amp = amp

        if dur is None:
            self._clamp.dur = self.tstop
        else:
            self._clamp.dur = dur

        self._clamp.delay = delay
        self._clamp.f0 = fmin
        self._clamp.f1 = fmax

        #self.clamp.fmin = fmin/1000
        #self.clamp.fmax = fmax/1000
        #print ("Don't forget to add the 'i' and 'f' vectors "
        #        "if you want to reord them!")
        self._clamp_setter()

    def zero_synapses(self, group=None, compartment=None):
        '''
        Function for quickly zeroing synapses prior to any clamp
        '''
        if group:
            for i in self.groups[group]:
                for x in self.compartments[i]:
                    setattr(self.__dict__[i], 'gbar_' + x, 0.)
        else:
            for x in self.compartments[compartment]:
                setattr(self.__dict__[compartment], 'gbar_' + x, 0.)


    def ga_change_gbar(self, runvals, percent=False, per_val = 100.):
        '''
        Function for GA.
        '''
        for key, val in runvals.items():
            #check if a group of compartments
            if key in self.groups:
                #if it is, iterate over each compartment in the group...
                for comp in self.groups[key]:
                    #for every parameter and value...
                    for k, v in val.items():
                        #reset the components value to the new GA value
                        setattr(self.__dict__[comp], k, v)
            #if not in a group, do individually
            else:
                for k, v in val.items():
                    setattr(self.__dict__[key], k, v)

    def _clamp_setter(self):
        """
        Semi-private function that handles the addition of single or multiple
        clamps or zaps.
        """
        self.clamp[str(self._clampcount)] = self._clamp
        cframe = inspect.currentframe()
        func = inspect.getframeinfo(cframe.f_back).function
        print ('Created ' + func + ' as clamp ' +str(self._clampcount)
                               + ' in ' + self.name + '.clamp...')
        print "Don't forget to add the vector(s) for recording!!!"
        self._clampcount +=1
