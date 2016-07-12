#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level dependency validation facilities.

Functions defined by this module typically validate the existence of both
mandatory and optional dependencies.
'''

# ....................{ IMPORTS                            }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: To raise human-readable exceptions on missing mandatory dependencies,
# the top-level of this module may import *ONLY* from packages guaranteed to
# exist at installation time (i.e., stock Python and BETSE packages).
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# ....................{ INITIALIZERS                       }....................
def init() -> None:
    '''
    Raise an exception if one or more mandatory runtime dependencies of BETSE
    are _not_ satisfiable.

    See Also
    ----------
    die_unless_satisfiable_all()
        Function to which this function defers.
    '''

    die_unless_satisfiable_all()

# ....................{ EXCEPTIONS                         }....................
def die_unless_satisfiable_all() -> None:
    '''
    Raise an exception unless all mandatory runtime dependencies of BETSE are
    **satisfiable** (i.e., importable and of a satisfactory version).

    Equivalently, this function raises an exception if at least one such
    dependency is unsatisfied. For importable unsatisfied dependencies with
    `setuptools`-specific metadata (e.g., `.egg-info/`-suffixed subdirectories
    of the `site-packages/` directory for the active Python 3 interpreter,
    typically created by `setuptools` at install time), this function
    additionally validates the versions of such dependencies to satisfy `betse`
    requirements.
    '''

    # Avoid circular import dependencies.
    from betse.util.py import modules

    # If the "pkg_resources" setuptools dependency is missing, raise an
    # exception *BEFORE* importing such dependency below.
    modules.die_unless_module(
        module_name='pkg_resources',
        exception_message='Mandatory dependency "pkg_resources" not found.',
    )

    # Validate these dependencies via "pkg_resources".
    from betse.lib import setuptool
    setuptool.die_unless_requirements_satisfiable_all()
