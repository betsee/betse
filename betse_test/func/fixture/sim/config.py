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
# from betse.util.type.types import CallableTypes, NoneType
from betse_test.func.fixture.sim import configapi
from betse_test.func.fixture.sim.configapi import SimTestState
# from betse_test.util import requests
from pytest import fixture

# ....................{ FIXTURES                           }....................
@fixture(scope='session')
def betse_sim_config_default(
    request: '_pytest.python.FixtureRequest',
    tmpdir_factory: '_pytest.tmpdir.tmpdir_factory',
) -> SimTestState:
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
        Further details on this fixture's return value.
    '''

    return configapi.make(
        request, tmpdir_factory,
        sim_config_dir_basename='default',
    )


# #FIXME: Refactor to use indirect fixture parametrization both here and above,
# #ultimately permitting these fixtures to be entirely eliminated. See also:
# #    http://docs.pytest.org/en/latest/example/parametrize.html#deferring-the-setup-of-parametrized-resources
@fixture(scope='session')
def betse_sim_config_visuals(
    request: '_pytest.python.FixtureRequest',
    tmpdir_factory: '_pytest.tmpdir.tmpdir_factory',
) -> SimTestState:
    '''
    Context manager-driven fixture creating a temporary simulation configuration
    enabling all supported animations and all features required by these
    animations _and_ returning a test-specific object encapsulating this
    configuration.

    See Also
    ----------
    :func:`betse_sim_config_anims`
        Further details on fixture parameters and return value.
    '''

    return configapi.make(
        request, tmpdir_factory,
        sim_config_dir_basename='visuals',
        config_modifier=lambda config: config.enable_visuals(),
    )


#FIXME: Parametrize this in the manner expected by config.enable_video().
# @fixture(scope='session')
# def betse_sim_config_video(
#     request: '_pytest.python.FixtureRequest',
#     tmpdir_factory: '_pytest.tmpdir.tmpdir_factory',
# ) -> SimTestState:
#     '''
#     Context manager-driven fixture creating a temporary simulation configuration
#     enabling encoding of at least one animation and all features required by
#     this animation as compressed video of the passed filetype with the preferred
#     matplotlib animation writer of the passed name _and_ returning a
#     test-specific object encapsulating this configuration.
#
#     See Also
#     ----------
#     :func:`betse_sim_config_anims`
#         Further details on fixture parameters and return value.
#     '''
#
#     return configapi.make(
#         request, tmpdir_factory,
#         config_modifier=lambda config: config.enable_video())
