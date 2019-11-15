#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level module and package importation facilities.

This submodule *only* defines functions implementing non-standard and
occasionally risky "black magic" fundamentally modifying Python's standard
importation semantics and mechanics. This submodule does *not* define
commonplace functions for dynamically importing modules or testing or
validating that importation.

See Also
----------
:mod:`betse.util.py.module.pymodname`
    Related submodule defining functions importing modules by name as well as
    testing and validating that importation.
'''

# ....................{ IMPORTS                           }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: To raise human-readable exceptions on missing mandatory dependencies,
# the top-level of this module may import *ONLY* from packages guaranteed to
# exist at installation time -- which typically means *ONLY* BETSE packages and
# stock Python packages.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

import sys
from betse.util.io.log import logs
from betse.util.type.types import type_check

# ....................{ CONVERTERS                        }....................
def to_str_modules_imported_name() -> str:
    '''
    Human-readable string of the fully-qualified names of all previously
    imported modules for the active Python interpreter.

    For readability, these names are sorted in lexicographic order and
    delimited by a newline.
    '''
    # Avoid circular import dependencies.
    from betse.util.type.iterable.mapping import mappings

    # Return the human-readable string produced by joining on newline a
    # generator comprehension yielding the colon-delimited name and value of
    # all environment variables sorted in lexicographic order.
    #
    # Note that the "environ" object is of non-standard type "os._Environ",
    # which the pprint.pformat() fails to recognize and hence format as a
    # "dict"-compatible mapping. Ergo, passing "environ" directly to the
    # mappings.to_str_flat() function would yield a non-human-readable string.
    # While this could be ameliorated by converting "environ" to a "dict" first
    # (e.g., "mappings.to_str_flat(dict(environ))"), doing so still produces
    # less human-readable output than the current approach.
    return mappings.to_str_flat(sys.modules)

# ....................{ REGISTRARS                        }....................
@type_check
def register_dir(dirname: str) -> None:
    '''
    Register all files and subdirectories of the directory with the passed
    dirname to be importable modules and packages (respectively) for the
    remainder of the active Python interpreter if this directory has yet to be
    registered *or* reduce to a noop otherwise (i.e., if this directory is
    registered already).

    Specifically, this function appends this dirname to the current
    :data:`sys.path` listing (in order) the dirnames of all directories to be
    iteratively searched for any module or package on first importing that
    module or package. To comply with Python standards in which the first item
    of this list is either the dirname of the directory containing the script
    from which this process was invoked *or* the empty string (signifying the
    current directory), this list is appended to rather than prepended to.

    Parameters
    ----------
    dirname : str
        Absolute or relative path of the directory to be registered.
    '''

    # Avoid circular import dependencies.
    from betse.util.path import dirs

    # Log this addition.
    logs.log_debug('Registering import directory: %s', dirname)

    # If this directory does *NOT* exist or is unreadable, raise an exception.
    dirs.die_unless_dir(dirname)

    # If the current PYTHONPATH already contains this directory...
    if dirname in sys.path:
        # Log this edge case.
        logs.log_debug('Ignoring already registered import directory.')

        # Reduce to a noop.
        return
    # Else, the current PYTHONPATH does *NOT* already contain this directory.

    # Append this directory to the current PYTHONPATH.
    sys.path.append(dirname)
