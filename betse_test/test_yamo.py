#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Temporary parametrization-based tests.
'''

# ....................{ IMPORTS                            }....................
from betse_test.util.mark.param import parametrize_test

# ....................{ TESTS                              }....................
@parametrize_test(
    aba=('ABA', 'BAB',),
    cac=('CAC', 'ACA',),
)
def test_params(aba, cac):
    assert aba in ('ABA', 'BAB',)
    assert cac in ('CAC', 'ACA',)
