#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level **Python-based multithreading** (i.e., multithreading implemented at
the Python API level rather than by third-party C extensions, libraries, or
packages) facilities.
'''

# ....................{ IMPORTS                           }....................
from betse.exceptions import BetsePyException
from betse.util.type.decorator.decmemo import func_cached

# ....................{ EXCEPTIONS                        }....................
def die_unless_gil() -> None:
    '''
    Raise an exception unless the active Python interpreter prohibits genuine
    multithreading via a Global Interpreter Lock (GIL).

    Equivalently, this function raises an exception if this interpreter is
    **GIL-less** (i.e., permits genuine multithreading and hence is *not*
    encumbered by a GIL).

    Raises
    ----------
    BetsePyException
        If this interpreter is GIL-less.

    See Also
    ----------
    :func:`is_gil`
        Further details.
    '''

    # If this Python interpreter is GIL-less, raise an exception.
    if not is_gil():
        raise BetsePyException(
            'Python interpreter GIL-less '
            '(i.e., permits genuine multithreading).')


def die_unless_threadable() -> None:
    '''
    Raise an exception unless the active Python interpreter supports
    multithreading (in principal) at the API level.

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

# ....................{ TESTERS                           }....................
@func_cached
def is_gil() -> bool:
    '''
    ``True`` only if the active Python interpreter prohibits genuine
    multithreading via a Global Interpreter Lock (GIL).

    Caveats
    ----------
    This function returning ``True`` does *not* imply this interpreter to
    prohibit multithreading (in principal) at the API level. Interpreters that
    prohibit genuine multithreading via a GIL typically still at least support
    standard multithreading APIs, albeit in a manner silently reducing multi-
    to single-core processing. Consider calling the :func:`is_threadable`
    function to clarify this edge case.

    This function returns ``True`` only if this interpreter is either CPython
    (i.e., the reference Python implementation) *or* PyPy; all other
    interpreters (e.g., IronPython, Jython) are well-known to *not* prohibit
    genuine multithreading via a GIL. While an admittedly non-future-proofed
    hack, the Python API currently provides no portable means of querying for
    the "GIL-ness" of an interpreter.
    '''

    # Avoid circular import dependencies.
    from betse.util.py import pyimpl

    # Return true only if this is either CPython or PyPy. See above.
    return pyimpl.is_cpython() or pyimpl.is_pypy()


@func_cached
def is_threadable() -> bool:
    '''
    ``True`` only if the active Python interpreter supports multithreading (in
    principal) at the API level.

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

    Caveats
    ----------
    This function returning ``True`` does *not* imply this interpreter to be
    **GIL-less** (i.e., permit genuine multithreading and hence is *not*
    encumbered by a GIL). Interpreters that prohibit genuine multithreading via
    a GIL typically still at least support standard multithreading APIs, albeit
    in a manner silently reducing multi- to single-core processing. Consider
    calling the :func:`is_gil` function to clarify this edge case.
    '''

    # Avoid circular import dependencies.
    from betse.util.py.module import pymodname

    # Return true only if a low-level optional module installed *ONLY* when
    # Python supports multithreading is importable.
    return pymodname.is_module('_thread')
