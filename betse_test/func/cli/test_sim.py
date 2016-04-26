#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Functional tests for BETSE's CLI testing all simulation-specific subcommands
(e.g., `betse try`).
'''

# ....................{ TESTS                              }....................
#FIXME: As the first order of business, refactor:
#
#* This to accept *ONLY* a "betse_sim_*"-style fixture.
#* The "_betse_sim" fixture to:
#  * Accept a "betse_cli" fixture.
#  * Classify the "CLITestRunner" instance returned by the "betse_cli" fixture
#    to the "SimTestContext" instance returned by the "_betse_sim" fixture.
#  * Disable all interactive plotting configuration options.
#FIXME: While useful, this should be split apart into its constituent steps,
#most of which should probably become fixtures for efficient reuse.

# def test_cli_try(betse_cli) -> None:
#     '''
#     Test the `betse try` subcommand.
#     '''
#
#     betse_cli('try')
