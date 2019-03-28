#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Unit tests exercising the :mod:`betse.util.type.iterable.iterables` submodule.
'''

# ....................{ IMPORTS                           }....................

# ....................{ TESTS ~ default                   }....................
#FIXME: Exercise all remaining edge cases... of which there are many. See the
#iterables.to_iterable() docstring for further details.
def test_to_iterable() -> None:
    '''
    Unit test the
    :func:`betse.util.type.iterable.iterables.to_iterable` function.
    '''

    # Defer heavyweight imports.
    from betse.util.type.iterable import iterables

    # Input iterable to be converted from.
    howl = (
        'incomparable blind streets of shuddering cloud',
        'and lightning in the mind leaping toward poles',
        'of Canada & Paterson,',
        'illuminating all the motionless world',
        'of Time between,',
    )

    # Exercise the edge case when this input iterable remains unconverted.
    assert iterables.to_iterable(iterable=howl) is howl

    # Exercise the edge case when this input iterable is to be converted into
    # an output iterable of differing type.
    assert iterables.to_iterable(iterable=howl, cls=set) == set(howl)
