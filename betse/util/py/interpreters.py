#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Low-level **Python implementation** (i.e., Python interpreter conforming to the
established Python syntax and semantics of the CPython reference implementation)
facilities.

Caveats
----------
Implementation-specific functions are generally considered poor form. Call these
functions _only_ where necessary.
'''

# ....................{ IMPORTS                            }....................
import platform
from betse.util.type.mappings import OrderedArgsDict

# ....................{ TESTERS                            }....................
def is_cpython() -> bool:
    '''
    `True` only if the active Python interpreter is an instance of the official
    CPython implementation.
    '''

    # Depressingly, this actually appears to be the most efficacious test.
    return get_name() == 'CPython'


def is_pypy() -> bool:
    '''
    `True` only if the active Python interpreter is an instance of the
    third-party PyPy implementation.
    '''

    # Depressingly, this actually appears to be the most efficacious test.
    return get_name() == 'PyPy'

# ....................{ GETTERS                            }....................
def get_name() -> str:
    '''
    Human-readable name of the active Python interpreter's implementation (e.g.,
    `CPython`, `IronPython`, `Jython`, `PyPy`).
    '''

    return platform.python_implementation()

# ....................{ GETTERS ~ metadata                 }....................
def get_metadata() -> OrderedArgsDict:
    '''
    Ordered dictionary synopsizing the active Python interpreter's
    implementation.

    This function aggregates the metadata reported by the reasonably
    cross-platform module `platform` into a simple dictionary.
    '''

    # This dictionary.
    metadata = OrderedArgsDict(
        'name', get_name(),
        'vcs revision', platform.python_revision() or 'none',
        'vcs branch', platform.python_branch() or 'none',
        'compiler', platform.python_compiler(),
    )

    # 2-tuple providing this interpreter's build number and date as strings.
    python_build = platform.python_build()

    # Append this metadata.
    metadata['build number'] = python_build[0]
    metadata['build data'] = python_build[1]

    # Return this dictionary.
    return metadata
