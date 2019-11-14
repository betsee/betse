#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Application-specific exception hierarchy.
'''

# ....................{ IMPORTS                           }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: To avoid race conditions during application startup, this module may
# import *ONLY* from modules guaranteed to exist at startup. This includes all
# standard Python and application modules but *NOT* third-party dependencies,
# which if unimportable will only be validated at some later time in startup.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

from abc import ABCMeta

# ....................{ EXCEPTIONS                        }....................
#FIXME: Define an __init__() method asserting that the passed exception message
#is non-None, which Python permits by default but which is functionally useless.
class BetseException(Exception, metaclass=ABCMeta):
    '''
    Abstract base class of all application-specific exceptions.
    '''

    pass


class BetseAttrException(BetseException):
    '''
    **Attribute** (i.e., variable or method bound to an object)-specific
    exception.
    '''

    pass


class BetseMetaAppException(BetseException):
    '''
    **Application metadata singleton** (i.e., application-wide object
    synopsizing application metadata via read-only properties)-specific
    exception.
    '''

    pass

# ....................{ EXCEPTIONS ~ call                 }....................
class BetseCallableException(BetseException):
    '''
    General-purpose exception applicable to all **callables** (e.g., functions,
    lambdas, methods, properties).
    '''

    pass


class BetseCallbackException(BetseCallableException):
    '''
    Callback-specific exception.
    '''

    pass


class BetseDecoratorException(BetseCallableException):
    '''
    Decorator-specific exception.
    '''

    pass


class BetseLambdaException(BetseCallableException):
    '''
    Lambda-specific exception.
    '''

    pass


class BetseParamException(BetseCallableException):
    '''
    **Parameter** (i.e., positional or keyword argument passed to a
    callable)-specific exception.
    '''

    pass

# ....................{ EXCEPTIONS ~ call : func          }....................
class BetseFunctionException(BetseCallableException):
    '''
    Function-specific exception.
    '''

    pass


class BetseFunctionUnimplementedException(
    BetseFunctionException, NotImplementedError):
    '''
    Unimplemented function-specific exception.

    This exception is typically raised from **unimplemented dependency- or
    platform-specific functions** (i.e., functions assuming a specific
    dependency or platform that have yet to be implemented).
    '''

# ....................{ EXCEPTIONS ~ call : method        }....................
class BetseMethodException(BetseCallableException):
    '''
    Method-specific exception.
    '''

    pass


class BetseMethodUnimplementedException(
    BetseMethodException, NotImplementedError):
    '''
    Unimplemented method-specific exception.

    This exception is typically raised from **unimplemented optional methods**
    (i.e., non-mandatory methods *not* intended to be called) of concrete
    subclasses of abstract base classes. While the optimal solution for
    defining **unimplemented mandatory methods** (i.e., non-optional methods
    also *not* intended to be called) is via the standard
    :call:`collections.abc.abstractmethod` decorator, there exists no standard
    alternative for defining optional methods. Hence, this exception.
    '''

    # ..................{ INITIALIZERS                      }..................
    def __init__(self) -> None:

        # Avoid circular import dependencies.
        from betse.util.type.call import callers

        # Initialize our superclass with a human-readable exception message.
        super().__init__('Method {}() unimplemented.'.format(
            callers.get_caller_basename()))

# ....................{ EXCEPTIONS ~ call : descriptor    }....................
class BetseDescriptorException(BetseCallableException):
    '''
    Descriptor-specific exception.
    '''
    pass


class BetseExprAliasException(BetseDescriptorException):
    '''
    Expression alias-specific exception.
    '''
    pass


# ....................{ EXCEPTIONS ~ cli                  }....................
class BetseCLIException(BetseException):
    '''
    General-purpose command-line interface (CLI) exception.
    '''

    pass


class BetseCLIArgException(BetseCLIException):
    '''
    Command-line interface (CLI) argument-specific exception.
    '''

    pass


class BetseCLIArgParserException(BetseCLIArgException):
    '''
    **Argument parser** (i.e., :class:`argparse.ArgumentParser`)-specific
    exception.
    '''

    pass

# ....................{ EXCEPTIONS ~ lib                  }....................
class BetseLibException(BetseException):
    '''
    General-purpose exception applicable to third-party dependencies.
    '''

    pass


class BetseMatplotlibException(BetseLibException):
    '''
    :mod:`matplotlib`-specific exception.
    '''

    pass


class BetsePyDotException(BetseLibException):
    '''
    :mod:`pydot`-specific exception.
    '''

    pass


class BetseYamlException(BetseException):
    '''
    Yet Another Markup Language (YAML)-specific exception.
    '''

    pass

# ....................{ EXCEPTIONS ~ log                  }....................
class BetseLogException(BetseException):
    '''
    Logging-specific exception.
    '''

    pass


class BetseLogRaceException(BetseLogException):
    '''
    Logging race condition-specific exception.

    Motivation
    ----------
    Concurrent attempts to log to the same physical file from multiple
    processes can and typically will produce fatal race conditions producing
    raised exceptions from one or more of these processes. On logfile rotation,
    each process will aggressively contend with each other process for write
    access to the same physical file to be rotated.

    The thread- *and* process-safe
    :class:`betse.util.io.log.conf.logconfhandle.RotatingFileHandler` handler
    obviates most of these concerns. Unfortunately, due to inadequacies in the
    Python ecosystem with respect to process-safe logging, that handler
    *cannot* constrain logfile access during rotation. Instead, on detecting
    exceptions produced by race conditions between multiple processes competing
    for access when attempting to emit log records, that handler temporarily
    halts the existing process for a negligible amount of the timeslice (e.g.,
    100ms) and repeats the attempt a negligible number of times (e.g., 8)
    *before* giving up and raising a fatal instance of this exception class.
    '''

    pass

# ....................{ EXCEPTIONS ~ math                 }....................
class BetseMathException(BetseException):
    '''
    General-purpose exception applicable to all low-level math algorithms.
    '''

    pass


class BetseMathLineException(BetseMathException):
    '''
    Line- and line segment-specific math exception.
    '''

    pass


class BetseMathPointException(BetseMathException):
    '''
    Point-specific math exception.
    '''

    pass


class BetseMathPolygonException(BetseMathException):
    '''
    Polygon-specific math exception.
    '''

    pass

# ....................{ EXCEPTIONS ~ module               }....................
class BetseModuleException(BetseException):
    '''
    Module-specific exception.
    '''

    pass


class BetsePackageException(BetseModuleException):
    '''
    Module-specific exception.
    '''

    pass

# ....................{ EXCEPTIONS ~ os                   }....................
class BetseOSException(BetseException):
    '''
    General-purpose low-level operating system (OS) exception.
    '''

    pass


class BetseShellEnvException(BetseOSException):
    '''
    Shell environment-specific exception.
    '''

    pass

# ....................{ EXCEPTIONS ~ os : proc            }....................
class BetseProcessException(BetseOSException):
    '''
    General-purpose low-level process exception.
    '''

    pass


class BetseProcessNotFoundException(BetseProcessException):
    '''
    Low-level exception raised on failing to find a given process.
    '''

    pass

# ....................{ EXCEPTIONS ~ path                 }....................
class BetsePathnameException(BetseException):
    '''
    Pathname-specific exception.
    '''

    pass


class BetsePathException(BetseException):
    '''
    Path-specific exception.
    '''

    pass


class BetseArchiveException(BetsePathException):
    '''
    Archive-specific exception.
    '''

    pass


class BetseDirException(BetsePathException):
    '''
    Directory-specific exception.
    '''

    pass


class BetseGitException(BetsePathException):
    '''
    Git-specific exception.
    '''

    pass

# ....................{ EXCEPTIONS ~ path : file          }....................
class BetseFileException(BetsePathException):
    '''
    File-specific exception.
    '''

    pass


class BetseCommandException(BetseFileException):
    '''
    **Command** (i.e., executable binary file)-specific exception.
    '''

    pass


class BetseImageException(BetseFileException):
    '''
    Image file-specific exception.
    '''

    pass

# ....................{ EXCEPTIONS ~ py                   }....................
class BetsePyException(BetseException):
    '''
    General-purpose low-level Python interpreter exception.

    This exception is appropriate for use in relation to low-level issues
    concerning the active Python interpreter (e.g., inability to retrieve this
    interpreter's absolute path).
    '''

    pass


class BetsePyFrozenException(BetsePyException):
    '''
    **Frozen executable** (i.e., Python codebase converted into a
    platform-specific executable)-specific exception.
    '''

    pass


class BetsePyIdentifierException(BetsePyException):
    '''
    **Python identifier** (i.e., class, module, or attribute name)-specific
    exception.
    '''

    pass

# ....................{ EXCEPTIONS ~ test                 }....................
class BetseTestException(BetseException):
    '''
    General-purpose exception pertaining to this application's test suite.
    '''

    pass


class BetseTestFixtureException(BetseTestException):
    '''
    Fixture-specific test exception.
    '''

    pass


class BetseTestHookException(BetseTestException):
    '''
    Hook-specific test exception.
    '''

    pass


class BetseTestParamException(BetseTestException):
    '''
    Parameter-specific test exception.
    '''

    pass

# ....................{ EXCEPTIONS ~ type                 }....................
class BetseTypeException(BetseException):
    '''
    General-purpose exception applicable to all types (i.e., classes).
    '''

    pass


class BetseEnumException(BetseTypeException):
    '''
    Enumeration-specific exception.
    '''

    pass


class BetseNumericException(BetseTypeException):
    '''
    **Numeric** (i.e., integer *and* float)-specific exception.
    '''

    pass


class BetseIntException(BetseTypeException):
    '''
    Integer-specific exception.
    '''

    pass

# ....................{ EXCEPTIONS ~ type : iter          }....................
class BetseIterableException(BetseTypeException):
    '''
    Iterable-specific exception.
    '''

    pass


class BetseSequenceException(BetseIterableException):
    '''
    Sequence-specific exception.
    '''

    pass

# ....................{ EXCEPTIONS ~ type : iter : map    }....................
class BetseMappingException(BetseIterableException):
    '''
    Dictionary-specific exception.
    '''

    pass


class BetseMappingKeyException(BetseMappingException):
    '''
    Dictionary key-specific exception.
    '''

    pass


class BetseMappingValueException(BetseMappingException):
    '''
    Dictionary value-specific exception.
    '''

    pass

# ....................{ EXCEPTIONS ~ type : str           }....................
class BetseStrException(BetseTypeException):
    '''
    String-specific  exception.
    '''

    pass


class BetseCharException(BetseStrException):
    '''
    **Character** (i.e., string of length 1)-specific  exception.
    '''

    pass



class BetseRegexException(BetseStrException):
    '''
    Regular exception-specific  exception.
    '''

    pass

# ....................{ EXCEPTIONS ~ sci : sim            }....................
class BetseSimException(BetseException):
    '''
    General-purpose simulation exception.
    '''

    pass


class BetseSimConfException(BetseSimException):
    '''
    Simulation configuration-specific exception.
    '''

    pass


class BetseSimEventException(BetseSimException):
    '''
    Simulation event-specific exception.
    '''

    pass


class BetseSimPhaseException(BetseSimException):
    '''
    Simulation phase-specific exception.
    '''

    pass


class BetseSimTissueException(BetseSimException):
    '''
    Simulation tissue-specific exception.
    '''

    pass


class BetseSimVectorException(BetseSimException):
    '''
    Vector-specific exception.
    '''

    pass

# ....................{ EXCEPTIONS ~ sci : sim : unstable }....................
class BetseSimUnstableException(BetseSimException):
    '''
    Simulation-specific exception indicating the current simulation to have
    unexpectedly failed due to computational instability.
    '''

    pass


class BetseSimUnstableNaNException(BetseSimUnstableException):
    '''
    Simulation-specific exception indicating the current simulation to have
    unexpectedly failed due to a **NaN-based computational instability** (i.e.,
    an operation producing a Numpy array containing at least one NaN value,
    typically due to division by zero or infinity).
    '''

    pass

# ....................{ EXCEPTIONS ~ sci : visual         }....................
class BetseSimVisualException(BetseSimException):
    '''
    Simulation visualization-specific exception, applicable to both plots and
    animations.
    '''

    pass


class BetseSimVisualLayerException(BetseSimVisualException):
    '''
    Simulation visualization layer-specific exception.
    '''

    pass

# ....................{ EXCEPTIONS ~ sci : pipe           }....................
class BetseSimPipeException(BetseSimException):
    '''
    Simulation pipeline-specific exception.
    '''

    pass


class BetseSimPipeRunnerUnsatisfiedException(BetseSimPipeException):
    '''
    Simulation pipeline-specific exception raised on attempting to run a runner
    with unsatisfied requirements (e.g., a post-simulation animation requiring
    extracellular spaces to be enabled by the current simulation
    configuration).

    Attributes
    ----------
    result : str
        Human-readable string justifying this failure. For generality, this
        string is neither capitalized *nor* punctuated.
    reason : str
        Human-readable string justifying this failure. For generality, this
        string is neither capitalized *nor* punctuated.
    '''

    # ..................{ INITIALIZERS                      }..................
    def __init__(self, result: str, reason: str) -> None:
        '''
        Initialize this exception.

        Parameters
        ----------
        result : str
            Human-readable string describing this failure. For generality, this
            string is expected to be neither capitalized *nor* punctuated.
        reason : str
            Human-readable string justifying this failure. For generality, this
            string is expected to be neither capitalized *nor* punctuated.
        '''

        # Capitalize all passed parameters.
        self.result = result
        self.reason = reason

        # Initialize our superclass, concatenating these strings into a single
        # human-readable exception message.
        super().__init__('{}: {}.'.format(self.result, self.reason))
