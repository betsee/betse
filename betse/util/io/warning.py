#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Low-level **warning** (i.e., non-fatal errors with associated types emitted by
the standard :mod:`warnings` module) facilities.
'''

# ....................{ IMPORTS                            }....................
import warnings
from betse.util.type.types import GeneratorType
from contextlib import contextmanager

# ....................{ MANAGERS                           }....................
@contextmanager
def warnings_ignored() -> GeneratorType:
    '''
    Context manager temporarily ignoring all warnings emitted for the duration
    of this context.
    '''

    # Defer to another lower-level context manager permitting specific types of
    # warnings to be ignored.
    with warnings.catch_warnings():
        # Unconditionally ignore all warnings.
        warnings.simplefilter("ignore")

        # Yield control to the body of the caller's "with" block.
        yield
