#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Functional tests for BETSE's CLI testing all simulation-specific subcommands
(e.g., `betse try`).
'''

#FIXME: Excise this module entirely, which is no longer relevant.

# ....................{ IMPORTS                            }....................
from pytest import fixture
from betse_test.mark.param import parametrize_fixture_serial

# ....................{ TESTS                              }....................
@parametrize_fixture_serial
@fixture(params=([(sum, int), (len, int)]), ids=['sum', 'len'])
def serial(request):
    return {
        'func': request.param[0],
        'cls': request.param[1],
    }

def test_serial(serial) -> None:
    assert isinstance(serial['func']([1, 2]), serial['cls'])
