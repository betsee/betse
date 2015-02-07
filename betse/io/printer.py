#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Low-level facilities for transparently printing and logging objects.
'''

# ....................{ IMPORTS                            }....................
from betse.io.file.log import logger
import traceback, sys

# ....................{ OUTPUTTERS                         }....................
def output_error(*objects) -> None:
    '''
    Print all passed objects to standard error *without* logging such message.
    '''
    print(*objects, file = sys.stderr)

# ....................{ CLASSES                            }....................
class LoggingPrinter(object):
    '''
    Low-level printing and logging of arbitrary objects.

    Attributes
    ----------
    _logger : Logger
        Low-level object implementing such logging.
    '''

    #FIXME: Cache a logging implementation here as a new private field
    #"_logger".

    def __init__(self):
        #FIXME: Actually, this should just be a global variable of the new
        #"betse.io.file.log" module.
        #FIXME: Right. And since we no longer need this as a field, there's
        #really no justification for having this as a class. Revert back to
        #simple functions, please.

        self._logger = None

    def print_exception(self, exception: Exception) -> None:
        '''
        Print the passed exception to standard error *and* log such exception.
        '''
        # If such printing itself raises an exception...
        try:
            assert isinstance(exception, Exception),\
                '"{}" not an exception.'.format(exception)

            #FIXME: Print such exception.

            # Log such exception *AFTER* printing such exception to standard error,
            # as logging is substantially more fragile and hence likely to itself
            # raise further exceptions.
            logger.exception(exception)
        # ...catch and print such exception using standard Python facilities
        # guaranteed not to raise additional exceptions.
        except Exception:
            output_error('print_exception() recursively raised exception:\n')
            traceback.print_exc()

# ....................{ SINGLETONS                         }....................
printer = LoggingPrinter()
'''
Singleton instance of LoggingPrinter(), simplifying printing and logging of
arbitrary objects in both the CLI and GUI interfaces.
'''

# --------------------( WASTELANDS                         )--------------------
# Utility functions manipulating **standard error** (i.e., the canonical file
# handle to which human-readable errors are written).
