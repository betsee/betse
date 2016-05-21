#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
BETSE-specific exception hierarchy.
'''

# ....................{ IMPORTS                            }....................
from abc import ABCMeta

# ....................{ EXCEPTIONS                         }....................
class BetseException(Exception, metaclass=ABCMeta):
    '''
    Abstract base class of all BETSE-specific exceptions.
    '''
    pass


class BetseExceptionLog(BetseException):
    '''
    Log-specific exception.
    '''
    pass

# ....................{ EXCEPTIONS ~ python                }....................
class BetseExceptionInterpreter(BetseException):
    '''
    General-purpose low-level Python interpreter exception.

    This exception is appropriate for use in relation to low-level issues
    concerning the active Python interpreter (e.g., inability to retrieve this
    interpreter's absolute path).
    '''
    pass


class BetseExceptionModule(BetseException):
    '''
    Module-specific exception.
    '''
    pass


class BetseExceptionFunction(BetseException):
    '''
    Function-specific exception.
    '''
    pass


class BetseExceptionLambda(BetseException):
    '''
    Lambda-specific exception.
    '''
    pass


class BetseExceptionMethod(BetseException):
    '''
    Method-specific exception.
    '''
    pass


class BetseExceptionMethodUnimplemented(BetseException, NotImplementedError):
    '''
    Unimplemented method-specific exception.

    This exception is typically raised from **unimplemented optional methods**
    (i.e., non-mandatory methods _not_ intended to be called) of concrete
    subclasses of abstract base classes. While the optimal solution for
    defining **unimplemented mandatory methods** (i.e., non-optional methods
    also _not_ intended to be called) is via the canonical `@abc.abstractmethod`
    decorator, there currently exists no canonical alternative for defining
    optional methods. Hence, this exception.
    '''

    def __init__():
        # Avoid circular import dependencies.
        from betse.util.py import callers
        super().__init__(
            'Optional method {}() unimplemented.'.format(
                callers.get_caller_basename()))

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

# ....................{ EXCEPTIONS ~ type                  }....................
class BetseExceptionRegex(BetseException):
    '''
    Regular exception-specific exception.
    '''
    pass

# ....................{ EXCEPTIONS ~ lib                   }....................
class BetseExceptionMatplotlib(BetseException):
    '''
    Matplotlib-specific exception.
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

# --------------------( WASTELANDS                         )--------------------
