#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Functional tests for BETSE's CLI testing all subcommands _except_ for
simulation-specific subcommands (e.g., `betse info`).
'''

# ....................{ TESTS                              }....................
def test_cli_sans_args(betse_cli) -> None:
    '''
    Test the `betse` command given no arguments.
    '''

    betse_cli()


def test_cli_info(betse_cli) -> None:
    '''
    Test the `betse info` subcommand.
    '''

    betse_cli('info')
