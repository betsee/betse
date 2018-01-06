#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Unit tests exercising floating point functionality.
'''

# ....................{ IMPORTS                            }....................
# import pytest

# ....................{ TESTS ~ base 10                    }....................
def test_get_base_10_exponent() -> None:
    '''
    Unit test the :func:`betse.util.type.numeric.floats.get_base_10_exponent`
    getter.
    '''

    # Defer heavyweight imports.
    from betse.util.type.numeric import floats

    # Assert this getter to behave as expected on common edge cases.
    assert floats.get_base_10_exponent(7.77e+66) == 66
    assert floats.get_base_10_exponent(6.66e-77) == -77
    assert floats.get_base_10_exponent(0.000123456789) == -4
    assert floats.get_base_10_exponent(100001234567.9) == 11


def test_get_base_10_precision() -> None:
    '''
    Unit test the :func:`betse.util.type.numeric.floats.get_base_10_precision`
    getter.
    '''

    # Defer heavyweight imports.
    from betse.util.type.numeric import floats

    # Assert this getter to behave as expected on common edge cases, including:
    #
    # * A large float formatted in scientific notation by the str() builtin.
    # * A small float formatted in scientific notation by the str() builtin.
    # * A small float with a leading digit formatted in decimal notation by the
    #   str() builtin.
    # * A small float with no leading digit formatted in decimal notation by the
    #   str() builtin.
    #
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # CAUTION: Python implicitly reformats all floating point numbers smaller
    # than 0.0001 when converted to strings in scientific rather than decimal
    # notation. Hence, small floats should be no smaller than 0.0001.
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    assert floats.get_base_10_precision(7.77e+66) == 66
    assert floats.get_base_10_precision(6.66e-77) == 77
    assert floats.get_base_10_precision(0.000123456789) == 12
    assert floats.get_base_10_precision( .000123456789) == 12
