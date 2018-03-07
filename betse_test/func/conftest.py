#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

# ....................{ IMPORTS ~ fixture                  }....................
'''
Global functional test configuration for BETSE's command-line interface (CLI).

``py.test`` implicitly imports all functionality defined by this module into all
CLI-specific functional test modules. As this functionality includes all
publicly declared functional fixtures in this ``fixture`` subpackage, these
tests may reference these fixtures without explicit imports.
'''

# ....................{ IMPORTS ~ fixture                  }....................
from betse_test.func.fixture.clier import betse_cli
from betse_test.func.fixture.sim.clisimer import (
    betse_cli_sim,
    betse_cli_sim_default,
    betse_cli_sim_compat,
)
