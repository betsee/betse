#!/usr/bin/env python3
# Copyright 2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Abstract base classes for timed event classes.
'''

# ....................{ IMPORTS                            }....................
from abc import ABCMeta #, abstractmethod, abstractstaticmethod
# from betse.util.type import types

# ....................{ BASE                               }....................
#FIXME: Actually use me.
class Animation(object, metaclass=ABCMeta):
    '''
    Abstract base class of all animation classes.

    Instances of this class animate the spatial distribution of modelled
    variables (e.g., Vmem) over all time steps of the simulation.

    Attributes
    ----------
    '''
    pass
