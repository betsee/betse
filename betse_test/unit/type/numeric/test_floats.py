#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Unit tests exercising floating point functionality.
'''

# ....................{ IMPORTS                            }....................
# import pytest

# ....................{ TESTS                              }....................
def test_get_precision() -> None:
    '''
    Unit test the :func:`betse.util.type.numeric.floats.get_precision` getter.
    '''

    # Defer heavyweight imports.
    from betse.util.type.numeric import floats

    # Assert this getter to behave as expected on common edge cases, including:
    #
    # * A large float formatted in scientific notation by the str() builtin.
    # * A small float formatted in scientific notation by the str() builtin.
    # * A small float formatted in decimal notation by the str() builtin.
    assert floats.get_precision(7.77e+66) == 66
    assert floats.get_precision(6.66e-77) == 77
    assert floats.get_precision(0.000123456789) == 12
