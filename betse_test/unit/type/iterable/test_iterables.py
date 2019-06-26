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

    # Input tuple to be converted from.
    howl_tuple = (
        'incomparable blind streets of shuddering cloud',
        'and lightning in the mind leaping toward poles',
        'of Canada & Paterson,',
        'illuminating all the motionless world',
        'of Time between,',
    )

    # Input generator to be converted from.
    howl_generator = (stanza for stanza in howl_tuple)

    # Test whether this conversion preserves this input tuple as is.
    assert iterables.to_iterable(iterable=howl_tuple) is howl_tuple

    # Test whether this conversion converts this input tuple into an output
    # iterable of differing type.
    assert iterables.to_iterable(iterable=howl_tuple, cls=set) == set(
        howl_tuple)

    # Test whether this conversion converts this input generator into an output
    # tuple identical to this input tuple.
    assert iterables.to_iterable(iterable=howl_generator) == howl_tuple
