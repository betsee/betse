#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2025 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level **context manager** (i.e., classes implementing the context manager
protocol and hence the ``__enter__`` and ``__exit__`` special methods)
facilities.
'''

# ....................{ IMPORTS                            }....................
from beartype.typing import Iterator
from contextlib import contextmanager

# ....................{ MANAGERS                           }....................
@contextmanager
def noop_context() -> Iterator[None]:
    '''
    Empty context manager exiting immediately after being entered.
    '''

    yield
