#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Global functional test configuration common to all interfaces (e.g., CLI, GUI).

`py.test` implicitly imports all functionality defined by this module into all
functional test modules. As this functionality includes all publicly declared
functional fixtures in this `fixture` subpackage, these tests may reference
these fixtures without explicit imports.
'''

# ....................{ IMPORTS ~ fixture : private        }....................
# Register private fixtures *BEFORE* public fixtures requiring these private
# fixtures. For unknown reasons, py.test requires this import in this "conftest"
# module rather than in the modules defining these public fixtures.

# ....................{ IMPORTS ~ fixture : public         }....................
# Register all public fixtures requiring that private fixture. (Order is
# insignificant here.)
from betse_test.func.fixture.sim.config import betse_sim_config_default
