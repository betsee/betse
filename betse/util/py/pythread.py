#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level **Python-based multithreading** (i.e., multithreading implemented at
the Python API level rather than by third-party C extensions, libraries, or
packages) facilities.
'''

# ....................{ IMPORTS                            }....................
from betse.exceptions import BetsePyException

# ....................{ EXCEPTIONS                         }....................
def die_unless_threadable() -> None:
    '''
    Raise an exception unless the active Python interpreter supports
    multithreading.

    Raises
    ----------
    BetsePyException
        If this interpreter does *not* support multithreading.

    See Also
    ----------
    :func:`is_threadable`
        Further details.
    '''

    # If this Python interpreter is *NOT* multithreadable, raise an exception.
    if not is_threadable():
        raise BetsePyException(
            'Python interpreter not multithreadable '
            '(i.e., "_thread" module not importable).')

# ....................{ TESTERS                            }....................
def is_threadable() -> bool:
    '''
    ``True`` only if the active Python interpreter supports multithreading.

    Equivalently, this function returns ``True`` only if the low-level standard
    :mod:`_thread` module is importable under this interpreter. This module is
    strictly optional and hence available only for systems both supporting
    multithreading *and* compiling Python with such support.

    Systems supporting multithreading includes:

    * Windows.
    * macOS.
    * Linux.
    * SGI IRIX.
    * Solaris 2.x.
    * Systems implementing the POSIX-compliant ``pthread`` API.
    '''

    # Avoid circular import dependencies.
    from betse.util.type import modules

    # Return true only if a low-level optional module installed *ONLY* when
    # Python supports multithreading is importable.
    return modules.is_module('_thread')
