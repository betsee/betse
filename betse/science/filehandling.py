#!/usr/bin/env python3
# Copyright 2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

"""

A number of functions used to save and load worlds and simulations.

"""

import pickle


#FIXME: For orthogonality, this function should probably be refactored to accept
#the same objects in the same order returned by the save functions below: e.g.,
#
#    def saveSim(savePath, sim, cells, p):
#       with open(savePath, 'wb') as f:
#           pickle.dump((sim, cells, p), f)
def saveSim(savePath,datadump):
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
