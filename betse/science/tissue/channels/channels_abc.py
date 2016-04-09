#!/usr/bin/env python3
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Abstract base classes of all channel classes.
'''

# ....................{ IMPORTS                            }....................
from abc import ABCMeta, abstractmethod

# ....................{ BASE                               }....................
class ChannelsABC(object, metaclass=ABCMeta):
    '''
    Abstract base class of all channel classes.

    Attributes
    ----------
    '''

    @abstractmethod
    def init(self, dyna,sim,cells,p):
        '''
        Do something.
        '''
        pass

    @abstractmethod
    def run(self, dyna,sim,cells,p):
        '''
        Do something.
        '''
        pass