#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Public simulation configuration fixtures.

These fixtures provide temporary simulation configurations and metadata on these
configurations, including:

* The absolute path of this configuration's on-disk YAML file.
* This configuration's in-memory dictionary deserialized from this file.
'''

# ....................{ IMPORTS                            }....................
from pytest import fixture

# ....................{ IMPORTS ~ fixtures                 }....................
import betse_test.func.fixture.sim.configbase

# ....................{ FIXTURES                           }....................
@fixture(scope='session')
def betse_sim_config_default(_betse_sim_config) -> 'SimTestConfig':
    '''
    Context manager-driven fixture creating a temporary default simulation
    configuration _and_ returning a test-specific object encapsulating this
    configuration.

    This fixture copies BETSE's default simulation configuration file,
    complete with all external assets (e.g., geometry masks) referenced and
    required by this file, into a temporary directory with basename `default`.

    Parameters
    ----------
    _betse_sim_config : SimTestConfig
        Test-specific object encapsulating this simulation configuration file.

    Returns
    ----------
    SimTestConfig
        Passed `_betse_sim_config` parameter unmodified.

    See Also
    ----------
    _betse_sim_config()
        Further details on return type and method of construction.
    '''

    # Write the default configuration to disk with only the requisite
    # modifications performed by this fixture (e.g., disabling interactivity).
    _betse_sim_config.config.overwrite()

    # Return the object encapsulating this default configuration.
    return _betse_sim_config
