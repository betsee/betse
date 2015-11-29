#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

#FIXME: Shift this module to the top-level "betse" package.

'''
`betse`-specific exception hierarchy.
'''

# ....................{ IMPORTS                            }....................
from abc import ABCMeta

# ....................{ EXCEPTIONS                         }....................
class BetseException(Exception, metaclass = ABCMeta):
    '''
    Abstract base class of all `betse`-specific exceptions.
    '''
    pass

class BetseExceptionLog(BetseException):
    '''
    Log-specific exception.
    '''
    pass

# ....................{ EXCEPTIONS ~ util                 }....................
class BetseExceptionModule(BetseException):
    '''
    Module-specific exception.
    '''
    pass

class BetseExceptionMethod(BetseException):
    '''
    Method-specific exception.
    '''
    pass

class BetseExceptionRegex(BetseException):
    '''
    Regular exception-specific exception.
    '''
    pass

# ....................{ EXCEPTIONS ~ science               }....................
class BetseExceptionParameters(BetseException):
    '''
    Parameters-specific exception.
    '''
    pass

class BetseExceptionSimulation(BetseException):
    '''
    Simulation-specific exception.
    '''
    pass

# ....................{ EXCEPTIONS ~ path                  }....................
class BetseExceptionPath(BetseException):
    '''
    Path-specific exception.
    '''
    pass

class BetseExceptionDir(BetseExceptionPath):
    '''
    Directory-specific exception.
    '''
    pass

class BetseExceptionFile(BetseExceptionPath):
    '''
    File-specific exception.
    '''
    pass

# --------------------( WASTELANDS                         )--------------------
