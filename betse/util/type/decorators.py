#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level **decorator** (i.e., callables dynamically wrapping arbitrary other
callables) facilities.
'''

# ....................{ IMPORTS                            }....................
from betse.util.type.types import CallableTypes  # type_check

# ....................{ DECORATORS                         }....................
# While type-checking these types would probably be advisable, this decorator is
# currently passed non-callables by "py.test". While resolving that issue would
# itself probably be advisable, we simply cannot be bothered. Hence, these types
# remain blithely unchecked.
def identity_decorator(callable: CallableTypes) -> CallableTypes:
    '''
    Identity decorator returning the decorated callable unmodified.
    '''

    return callable
