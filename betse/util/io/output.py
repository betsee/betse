#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Low-level output facilities.
'''

# ....................{ IMPORTS                            }....................
import sys

# ....................{ OUTPUTTERS                         }....................
def error(*objects) -> None:
    '''
    Print all passed objects to standard error *without* logging such objects.
    '''
    print(*objects, file = sys.stderr)

# --------------------( WASTELANDS                         )--------------------
# from betse.io.file.log import logger
# ....................{ CLASSES                            }....................
# class LoggingPrinter(object):
#     '''
#     Low-level printing and logging of arbitrary objects.
#
#     Attributes
#     ----------
#     '''
#
#     #FUXME: Cache a logging implementation here as a new private field
#     #"_logger".
#
#     def __init__(self):
#         #FUXME: Actually, this should just be a global variable of the new
#         #"betse.io.file.log" module.
#         #FUXME: Right. And since we no longer need this as a field, there's
#         #really no justification for having this as a class. Revert back to
#         #simple functions, please.

# # ....................{ SINGLETONS                         }....................
# printer = LoggingPrinter()
# '''
# Singleton instance of LoggingPrinter(), simplifying printing and logging of
# arbitrary objects in both the CLI and GUI interfaces.
# '''

# Utility functions manipulating **standard error** (i.e., the canonical file
# handle to which human-readable errors are written).
