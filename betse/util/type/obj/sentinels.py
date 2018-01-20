#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
**Sentinel** (i.e., objects of arbitrary unique value, typically used to
distinguish edge-case results from the standard ``None`` singleton) facilities.
'''

# ....................{ CLASSES                            }....................
class Sentinel(object):
    '''
    Class encapsulating sentinel objects of arbitrary (albeit unique) value.

    Instances of this class are intended to be used as placeholder objects in
    iterables, typically to identify erroneous or edge-case algorithm input.
    '''

    def __repr__(self) -> str:
        '''
        Human- and machine-readable representation of this sentinel.

        This method has been overridden purely to improve the debuggability of
        algorithms requiring instances of this class.
        '''

        return 'Sentinel()'

# ....................{ CONSTANTS                          }....................
SENTINEL = Sentinel()
'''
Sentinel object of arbitrary value.

This object is internally leveraged by various utility functions (e.g.,
:func:`betse.util.type.iterables.zip_isometric`) to identify erroneous and
edge-case input (e.g., iterables of insufficient length).
'''
