#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

#FIXME: Excise this temporary file, please.

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
def test_yamo_fix(request, aba, cac, betse_zyz):
    assert aba in ('ABA', 'BAB',)
    assert cac in ('CAC', 'ACA',)
    assert betse_zyz in ('zyz', 'pqp',)

# ....................{ TESTS ~ more                       }....................
# @parametrize_test(fixtures={
#     'betse_sim_config': ({},)
# })
def test_yamo_conf(
    # betse_cli_sim: 'CLITester',
    betse_sim_config: 'SimTestState',
) -> None:
    '''
    Test simulating with the default simulation configuration.

    Parameters
    ----------
    betse_cli_sim : CLITester
        Object encapsulating CLI-driven simulation testing.
    betse_sim_config_default : SimTestState
        Object encapsulating the default simulation configuration.
    '''

    # Test the currently parametrized simulation-specific BETSE CLI subcommand.
    # betse_cli_sim()
    pass
