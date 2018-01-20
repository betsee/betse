#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level **generator** (i.e., objects satisfying the standard generator API)
facilities.
'''

# ....................{ IMPORTS                            }....................
from betse.util.type.types import GeneratorType

# ....................{ GENERATORS                         }....................
def empty_generator() -> GeneratorType:
    '''
    Empty generator yielding... absolutely nothing.

    See Also
    ----------
    https://stackoverflow.com/a/36658865/2809027
        StackOverflow answer strongly inspiring this implementation.
    '''

    yield from ()
