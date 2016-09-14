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
from betse.util.type import strs
# from betse.util.type.types import CallableTypes, NoneType
from betse_test.func.fixture.sim import configapi
from betse_test.func.fixture.sim.configapi import SimTestState
from betse_test.util import requests
from pytest import fixture

# ....................{ FIXTURES                           }....................
#FIXME: Refactor as follows:
#
#* Shift the entirety of the "configapi" submodule into this submodule.
#* Rename this submodule to "betse_test.fixture.sim_config". Note the shift
#  upwards into a fixture subdirectory accessible by both functional and
#  unit tests.
#* Remove the "configapi" submodule.
#* Remove the now-empty "betse_test.func.fixture.sim" subpackage.

@fixture()
def betse_sim_config(
    request: '_pytest.python.FixtureRequest',
    tmpdir_factory: '_pytest.tmpdir.tmpdir_factory',
) -> SimTestState:
    '''
    Per-test fixture creating a temporary default simulation configuration file
    and returning an object encapsulating the contents of this file.

    Configuration Modifications (On-disk)
    ----------
    This fixture copies BETSE's default simulation configuration file,
    complete with all external assets (e.g., geometry masks) referenced and
    required by this file, into a temporary directory whose basename is the name
    of the test requesting this fixture excluding the prefixing substring
    `test_`. As example, when requested by the `test_cli_sim_default` test,
    this fixture creates a temporary simulation configuration file
    `{tmpdir}/cli_sim_default/sim_config.yaml` (e.g.,
    `/tmp/pytest-0/cli_sim_default/sim_config.yaml`), where `{tmpdir}` is the
    absolute path of this test session's root temporary directory (e.g.,
    `/tmp/pytest-0/`).

    This temporary directory and hence simulation configuration will be
    accessible _only_ for the duration of the current test. Subsequently run
    tests and fixtures may _not_ safely reuse this configuration.

    Configuration Modifications (In-memory)
    ----------
    This fixture also transforms the in-memory instance of the
    :class:`betse.science.config.wrapper.SimConfigWrapper` class encapsulating
    this configuration as follows:

    * All configuration options either requiring interactive input _or_
      displaying interactive output are disabled (e.g., plots, animations).
    * The space and time costs associated with simulating this configuration
      are safely minimized in a manner preserving all features.

    Since this fixture does _not_ write these changes back to this file, the
    parent fixture or test is expected to do so manually (e.g., by calling the
    :meth:`betse.science.config.wrapper.SimConfigWrapper.overwrite` method on
    the `config` attribute of the object returned by this fixture).

    Parameters
    ----------
    request : _pytest.python.FixtureRequest
        Builtin fixture parameter describing the parent fixture or test of this
        fixture (and similar contextual metadata).
    tmpdir_factory : _pytest.tmpdir.tmpdir_factory
        Builtin session-scoped fixture whose `mktemp()` method returns a
        :class:`py.path.local` instance encapsulating a new temporary directory.

    Returns
    ----------
    SimTestState
        Test-specific object encapsulating a temporary simulation configuration
        file specific to the current test.

    See Also
    ----------
    :class:`betse_test.fixture.sim_config.SimTestState`
        Class of the object returned by this fixture.
    '''

    # Name of the current test.
    test_name = requests.get_tested_name(request)
    # print('    request.node: {}'.format(request.node))
    # print('    test_name: {}'.format(test_name))

    # Basename of the temporary directory containing this configuration file,
    # set to the name of the current test excluding the prefixing "test_".
    sim_config_dir_basename = strs.remove_prefix(
        text=test_name,
        prefix='test_',
        exception_message=(
            'Test name "{}" not prefixed by "test_".'.format(test_name)),
    )

    # Create this temporary directory and wrap this directory's absolute path
    # with a high-level "py.path.local" object. See also:
    #     http://pytest.org/latest/tmpdir.html#the-tmpdir-factory-fixture
    sim_config_dirpath = tmpdir_factory.mktemp(sim_config_dir_basename)

    # Absolute path of this configuration file in this temporary directory.
    sim_config_filepath = sim_config_dirpath.join('sim_config.yaml')

    # Test-specific object encapsulating this simulation configuration file.
    sim_state = SimTestState(config_filepath=sim_config_filepath)

    # Return this object.
    return sim_state

# ....................{ FIXTURES ~ obsolete                }....................
#FIXME: Excise.
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

    sim_state = configapi.make(
        request, tmpdir_factory,
        sim_config_dir_basename='default',
    )
    sim_state.config.overwrite()
    return sim_state



#FIXME: Excise.
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

    sim_state = configapi.make(
        request, tmpdir_factory,
        sim_config_dir_basename='visuals',
        config_modifier=lambda config: config.enable_visuals(),
    )
    sim_state.config.overwrite()
    return sim_state


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
