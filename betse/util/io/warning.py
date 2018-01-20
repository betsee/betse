#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level **warning** (i.e., non-fatal errors with associated types emitted by
the standard :mod:`warnings` module) facilities.
'''

# ....................{ IMPORTS                            }....................
import warnings
from betse.util.type.types import type_check, ClassType, GeneratorType
from contextlib import contextmanager

# ....................{ MANAGERS                           }....................
#FIXME: For orthogonality, rename this function to ignoring_warnings().
@contextmanager
@type_check
def warnings_ignored(warning_class: ClassType = Warning) -> GeneratorType:
    '''
    Single-shot context manager temporarily ignoring all warnings of the passed
    type emitted by the :mod:`warnings` module for the duration of this context.

    Caveats
    -----------
    This context manager is single-shot and hence *not* safely reusable across
    multiple ``with`` blocks. Why? Because the low-level
    :class:`warnings.catch_warnings` class to which this context manager defers
    is itself single-shot, explicitly raising exceptions on reuse. While
    lamentable, this overhead appears to be unavoidable.

    Parameters
    -----------
    warning_class : optional[ClassType]
        Base class of all warnings to be ignored by this context manager.
        Defaults to :class:`Warning`, the superclass of all warning
        subclasses (builtin and user-defined), in which case this context
        manager will ignore *all* warnings.
    '''

    # Isolate all side effects produced by the "warnings" module to this block.
    with warnings.catch_warnings():
        # Ignore all warnings emitted by this module that are instances of the
        # previously passed base class.
        warnings.simplefilter(action='ignore', category=warning_class)

        # Yield control to the body of the caller's "with" block.
        yield
