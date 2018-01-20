#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level **iterator** (i.e., objects satisfying the standard iterator API)
facilities.
'''

# ....................{ IMPORTS                            }....................
from betse.util.type.types import IterableTypes

# ....................{ GENERATORS                         }....................
def empty_iterator() -> IterableTypes:
    '''
    Empty iterator iterating over... absolutely nothing.

    See Also
    ----------
    https://stackoverflow.com/a/10621647/2809027
        StackOverflow answer strongly inspiring this implementation.
    '''

    # It is clever. It is Python.
    return iter(())
