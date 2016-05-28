#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Functional tests for BETSE's CLI testing all simulation-specific subcommands
(e.g., `betse try`).
'''

# ....................{ IMPORTS                            }....................
import pytest
# from betse_test.mark.skip import skip

# ....................{ TESTS                              }....................
#FIXME: If this actually works, generalize into a parametrized fixture and try
#again.

@pytest.mark.serial_parametrized
# @pytest.mark.parametrize(('func', 'cls'), [(sum, list), (len, int)])
@pytest.mark.parametrize(('func', 'cls'), [(sum, int), (len, int)])
def test_serial(func, cls) -> None:
    assert isinstance(func([1, 2]), cls)
