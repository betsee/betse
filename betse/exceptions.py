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
#FIXME: Define an __init__() method asserting that the passed exception message
#is non-None, which Python permits by default but which is functionally useless.
class BetseException(Exception, metaclass=ABCMeta):
    '''
    Abstract base class of all BETSE-specific exceptions.
    '''
    pass


class BetseLogException(BetseException):
    '''
    Log-specific exception.
    '''
    pass

# ....................{ EXCEPTIONS ~ python                }....................
class BetseInterpreterException(BetseException):
    '''
    General-purpose low-level Python interpreter exception.

    This exception is appropriate for use in relation to low-level issues
    concerning the active Python interpreter (e.g., inability to retrieve this
    interpreter's absolute path).
    '''
    pass


class BetseModuleException(BetseException):
    '''
    Module-specific exception.
    '''
    pass

# ....................{ EXCEPTIONS ~ callable              }....................
class BetseDecoratorException(BetseException):
    '''
    Decorator-specific exception.
    '''
    pass


class BetseFunctionException(BetseException):
    '''
    Function-specific exception.
    '''
    pass


class BetseLambdaException(BetseException):
    '''
    Lambda-specific exception.
    '''
    pass


class BetseMethodException(BetseException):
    '''
    Method-specific exception.
    '''
    pass


class BetseMethodUnimplementedException(BetseException, NotImplementedError):
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

    def __init__(self):
        # Avoid circular import dependencies.
        from betse.util.py import callers
        super().__init__(
            'Optional method {}() unimplemented.'.format(
                callers.get_caller_basename()))

# ....................{ EXCEPTIONS ~ path                  }....................
class BetsePathException(BetseException):
    '''
    Path-specific exception.
    '''
    pass


class BetseDirException(BetsePathException):
    '''
    Directory-specific exception.
    '''
    pass


class BetseFileException(BetsePathException):
    '''
    File-specific exception.
    '''
    pass


class BetseCommandException(BetseFileException):
    '''
    Command-specific exception.
    '''
    pass

# ....................{ EXCEPTIONS ~ type                  }....................
class BetseEnumException(BetseException):
    '''
    Enumeration-specific exception.
    '''
    pass


class BetseIntegerException(BetseException):
    '''
    Integer-specific exception.
    '''
    pass


class BetseStringException(BetseException):
    '''
    String-specific exception.
    '''
    pass


class BetseRegexException(BetseException):
    '''
    Regular exception-specific exception.
    '''
    pass

# ....................{ EXCEPTIONS ~ lib                   }....................
class BetseMatplotlibException(BetseException):
    '''
    Matplotlib-specific exception.
    '''
    pass

# ....................{ EXCEPTIONS ~ arg                  }....................
#FIXME: Raise this exception throughout the "betse.cli" subpackage.
class BetseArgumentException(BetseException):
    '''
    Argument-specific exception.
    '''

    pass


class BetseArgumentParserException(SystemExit):
    '''
    :class:`betse.script.argparse.ArgumentParser`-specific exception connoting
    the :meth:`betse.script.argparse.parse_args` to have unsuccessfully parsed
    the argument list passed to the current BETSE script.
    '''

    pass

# ....................{ EXCEPTIONS ~ science               }....................
class BetseParametersException(BetseException):
    '''
    Parameters-specific exception.
    '''
    pass


class BetseSimulationException(BetseException):
    '''
    Simulation-specific exception.
    '''
    pass


class BetseSimulationInstabilityException(BetseSimulationException):
    '''
    Simulation-specific exception connoting the current simulation to have
    gone computationally unstable.
    '''
    pass
