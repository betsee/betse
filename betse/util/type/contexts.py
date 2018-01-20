#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level **context manager** (i.e., classes implementing the context manager
protocol and hence the ``__enter__`` and ``__exit__`` special methods)
facilities.
'''

# ....................{ IMPORTS                            }....................
from betse.util.type.types import GeneratorType
from contextlib import contextmanager

# ....................{ MANAGERS                           }....................
@contextmanager
def noop_context() -> GeneratorType:
    '''
    Empty context manager exiting immediately after being entered.
    '''

    yield
