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
    Noop decorator, also referred to as the **identity decorator.**

    This decorator returns the decorated callable (e.g., function) unmodified.
    '''

    return obj
