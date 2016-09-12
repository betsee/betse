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
from betse.util.type.types import MappingType
from betse_test.func.fixture.sim import configapi
from betse_test.func.fixture.sim.configapi import SimTestState
from betse_test.util import requests
from pytest import fixture

# ....................{ FIXTURES                           }....................
#FIXME: This reusable fixture would seem to obviate the need for a separate
#configapi.make() factory function. That being the case, shift the body of the
#configapi.make() factory function into the body of this fixture function and
#remove the former entirely.
#FIXME: Document this fixture to optionally require indirect parametrization,
#where the parametrized value if any *MUST* be a dictionary. If no parameter is
#supplied, this parameter defaults to the empty dictionary. Else:
#
#* The value of the "config_modifier" key if any *MUST* be either a callable or
#  "None".

# @fixture(scope='session')
@fixture()
def betse_sim_config(
    request: '_pytest.python.FixtureRequest',
    tmpdir_factory: '_pytest.tmpdir.tmpdir_factory',
) -> SimTestState:

    # Keyword arguments with which this fixture is indirectly parametrized if
    # any or the empty dictionary otherwise.
    kwargs = requests.get_fixture_param(
        request,
        default_value={},
        type_expected=MappingType,
    )

    return configapi.make(
        request, tmpdir_factory,
        config_modifier=kwargs.get('config_modifier', None),
    )


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


#FIXME: Refactor to use indirect fixture parametrization both here and above,
#ultimately permitting these fixtures to be entirely eliminated. See also:
#    http://docs.pytest.org/en/latest/example/parametrize.html#deferring-the-setup-of-parametrized-resources
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
