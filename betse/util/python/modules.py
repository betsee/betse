#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Low-level module facilities.
'''

# ....................{ IMPORTS                            }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: To raise human-readable exceptions on missing mandatory dependencies,
# the top-level of this module may import *ONLY* from packages guaranteed to
# exist at installation time (e.g., stock Python packages).
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
from betse.util.exceptions import BetseExceptionModule
import importlib

# ....................{ EXCEPTIONS                         }....................
def die_unless(
    module_name: str, exception_message: str = None) -> None:
    '''
    Raise an exception with the passed message if the module with the passed
    name is *not* importable.

    If no message is passed, a default message is synthesized from the passed
    module name.
    '''
    # If such module is missing, raise an exception.
    if not is_module(module_name):
        # If no exception message was passed, synthesize one from such name.
        if not exception_message:
            exception_message = 'Module "{}" not found.'.format(module_name)
        assert isinstance(exception_message, str),\
            '"{}" not a string.'.format(exception_message)

        # Raise such exception.
        raise BetseExceptionModule(exception_message)

# ....................{ TESTERS                            }....................
def is_module(module_name: str) -> bool:
    '''
    True if the module with the passed name is importable under the active
    Python 3 interpreter.
    '''
    assert isinstance(module_name, str),\
        '"{}" not a string.'.format(module_name)
    return importlib.find_loader(module_name)

# --------------------( WASTELANDS                         )--------------------
    # Such module is assumed to signify a mandatory `betse` dependency. Likewise,
    # such requirements are assumed to be in `setuptools` format (e.g.,
    # `numpy >= 1.9.0`).

#  Ideally, such
    # exception would be an instance of the setuptools-specific
    # "DistributionNotFound" class. Yet, as setuptools and hence such class is
    # *NOT* guaranteed to be importable, the conventional Python exception for
    # import errors is raised instead.
        # if not importlib.find_loader(requirement.project_name):
        #     raise DistributionNotFound(
        #         'Mandatory dependency {} not found.'.format(requirement))
    # for module_name in {
    #     'scipy'
    # }:
    # For each such dependency, this function also attempts to validate such
    # dependency's version. Specifically, if such dependency both exists *and*,
    # an exception is raised if

    # Caveats
    # ----------
    # Dependency versions are *not* validated. This is subject to change
