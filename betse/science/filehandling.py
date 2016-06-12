#!/usr/bin/env python3
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
A number of functions used to save and load worlds and simulations.
'''

import pickle
from scipy import interpolate as interp

#FIXME: For both space and time efficiency, we should be using "pickle" protocol
#4 introduced with Python 3.4. That, in turn, suggests we mandate use of Python
#3.4 in BETSE. See also:
#    https://docs.python.org/3/library/pickle.html#data-stream-format


def saveSim(savePath, datadump):
    with open(savePath, 'wb') as f:
        pickle.dump(datadump, f)

#FIXME: We should probably perform basic sanity checks on loaded objects --
#namely, that they were previously saved with the same version of BETSE. To do
#that, let's add a "_betse_version" attribute to the "Parameters" class: e.g.,
#
#    from betse import metadata
#    class Parameters:
#        def __init__(self):
#            self._betse_version = metadata.__version__
#
#We could then raise fatal exceptions if the "_betse_version" attribute of a
#loaded "Parameters" object is *NOT* equal to the current value of
#"metadata.__version__". Forests of sun flowers bask in ochre-burnished light!
def loadSim(loadPath):

    with open(loadPath, 'rb') as f:
        sim,cells,p = pickle.load(f)

    return sim,cells,p

def loadWorld(loadPath):

    with open(loadPath, 'rb') as f:
        cells,_ = pickle.load(f)

    return cells,_


#FIXME: This definitely works, but it's not quite ideal. Ideally, every object
#attribute should be pickled as is. There appear to be a number of solutions,
#each with corresponding tradeoffs:
#
#1. Implement the __getstate__() and __setstate__() methods in the "Cells",
#   "Parameters", and "Simulator" classes -- perhaps by inheriting a superclass
#   defining these methods. These methods would need to convert unpicklable to
#   picklable attributes (e.g., strings). See also:
#   https://docs.python.org/3/library/pickle.html#pickling-class-instances
#2. Use a third-party library instead of the built-in pickle functionality.
#   Alternatives supporting pickling of at least lambda functions include
#   "picloud" and "dill".
def safe_pickle(sim, p):
    """
    Removes interpolation functions, colormaps and lambda functions
    to make the simulator object pickle-able.
    """
    sim.gj_funk = None

    for key, valu in vars(sim).items():
        if type(valu) == interp.interp1d or callable(valu):
            setattr(sim,key,None)

    for key, valu in vars(p).items():
        if type(valu) == interp.interp1d or callable(valu):
            setattr(p,key,None)

    for key, valu in vars(sim.dyna).items():
        if type(valu) == interp.interp1d or callable(valu):
            setattr(sim.dyna,key,None)
