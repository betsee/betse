#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Low-level **decorator** (i.e., callables dynamically wrapping arbitrary other
callables) facilities.
'''

# ....................{ IMPORTS                            }....................

# ....................{ DECORATORS                         }....................
def noop(obj: object) -> object:
    '''
    Identity decorator returning the decorated callable unmodified.
    '''

    return obj
