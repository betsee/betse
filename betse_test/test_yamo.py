#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Temporary parametrization-based tests.
'''

# ....................{ IMPORTS                            }....................
from betse_test.util.mark.param import parametrize_test
from pytest import fixture

# ....................{ FIXTURES                           }....................
@fixture()  # scope='session'
def betse_zyz(request) -> str:
    # print('\n!!!!!!!!!!betse_zyz() returning {}'.format(pqp))
    return request.param.lower()

# ....................{ TESTS                              }....................
@parametrize_test(
    params={
        'western_dragon': ('Celedyr', 'Hestaby',),
        'eastern_dragon': ('Masaru', 'Ryumyo',),
    },
    ids=('bad-dragons', 'good-dragons',),
)
def test_yamo_unfix(western_dragon, eastern_dragon):
    assert western_dragon in ('Celedyr', 'Hestaby',)
    assert eastern_dragon in ('Masaru', 'Ryumyo',)


@parametrize_test(
    params={
        'aba': ('ABA', 'BAB',),
        'cac': ('CAC', 'ACA',),
    },
    fixtures={
        'betse_zyz': ('ZYZ', 'PQP',)
    },
)
def test_yamo_fix(aba, cac, betse_zyz):
    assert aba in ('ABA', 'BAB',)
    assert cac in ('CAC', 'ACA',)
    assert betse_zyz in ('zyz', 'pqp',)
