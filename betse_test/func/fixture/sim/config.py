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
from betse_test.func.fixture.sim import configapi
from pytest import fixture

# ....................{ FIXTURES                           }....................
#FIXME: Non-ideal repetition. Ideally, we would define a new BETSE-specific
#utility fixture decorator @sim_config_fixture that:
#
#* Forced use of session scope, as repeated below.
#* Forced use of the "request" and "tmpdir_factory" fixtures, as repeated below.
#
#Python decoration is non-trivial. This is left as an exercise for the reader.
#Since py.test only supports callable- rather than class-based fixtures,
#fixture inheritance is right out -- which is what we'd *REALLY* prefer. One
#option could be to implement this as an xunit-style class-based fixture, which
#(in theory) would permit inheritance. Anyway, a discourse for another time.

@fixture(scope='session')
def betse_sim_config_default(
    request: '_pytest.python.FixtureRequest',
    tmpdir_factory: '_pytest.tmpdir.tmpdir_factory',
) -> 'SimTestState':
    '''
    Context manager-driven fixture creating a temporary default simulation
    configuration _and_ returning a test-specific object encapsulating this
    configuration.

    This fixture copies BETSE's default simulation configuration file,
    complete with all external assets (e.g., geometry masks) referenced and
    required by this file, into a temporary directory with basename `default`.

    Parameters
    ----------
    request : _pytest.python.FixtureRequest
        Builtin fixture required by the `configapi.make()` factory function.
    tmpdir_factory : _pytest.tmpdir.tmpdir_factory
        Builtin fixture required by the `configapi.make()` factory function.

    Returns
    ----------
    SimTestState
        Test-specific object encapsulating this simulation configuration file.

    See Also
    ----------
    _betse_sim_config
        For further details on return type and method of construction.
    '''

    # # Test-specific object encapsulating this simulation configuration file.
    sim_state = configapi.make(request, tmpdir_factory)

    # Write the default configuration to disk with only the requisite
    # modifications performed by this fixture (e.g., disabling interactivity).
    sim_state.config.overwrite()

    # Return this encapsulation object.
    # print('Here!!!!!!!!!!!!!')
    return sim_state
