#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
**Global test configuration** (i.e., early-time configuration guaranteed to be
run by :mod:`pytest` *after* passed command-line arguments are parsed) for
this test suite.

:mod:`pytest` implicitly imports *all* functionality defined by this module
into *all* submodules of this subpackage.

See Also
----------
:mod:`conftest`
    Root test configuration applied before this configuration.
'''

# ....................{ IMPORTS                           }....................
import pytest

# ....................{ IMPORTS ~ fixture : autouse       }....................
# Import fixtures automatically run at the start of the current test session
# and hence typically *NOT* manually required by specific tests.

from betse_test.fixture.metaapper import betse_app_meta

# ....................{ IMPORTS ~ fixture : manual        }....................
# Import fixtures required to be manually required by specific tests.

from betse_test.fixture.tempdirer import betse_temp_dir
from betse_test.fixture.simconf.simconfer import (
    betse_sim_conf,
    betse_sim_conf_default,
    betse_sim_conf_compat,
)

# ....................{ GLOBALS                           }....................
EXPORT_SIM_CONF_DIRNAME = None
'''
Absolute or relative dirname of the target directory to export (i.e.,
recursively copy) each source simulation configuration directory into if any
*or* ``None`` otherwise.

See Also
----------
:func:`pytest_configure`
    Further details.
'''

# ....................{ HOOKS ~ plugin                    }....................
def pytest_configure(config) -> None:
    '''
    Hook run immediately *after* both parsing all :mod:`pytest` command-line
    options and loading all third-party :mod:`pytest` plugins (including
    application-specific ``conftest`` scripts).

    Specifically:

    * The global :attr:`EXPORT_SIM_CONF_DIRNAME` variable is defined as follows
      for subsequent lookup from module scope (e.g., pytest markers):

      * If the custom ``--export-sim-conf-dir`` command-line option was passed
        by the caller (as parsed by the root :mod:`conftest` module), this
        variable's value is that of this option's (i.e., the absolute or
        relative dirname of the target directory to export and hence
        recursively copy each source simulation configuration directory into).
      * Else, this variable's value is ``None``.

    * If the external ``${DISPLAY}`` environment variable is currently set
      (e.g., to the X11-specific socket to be connected to display GUI
      components), unset this variable. Permitting this variable to remain set
      would permit tests erroneously attempting to connect to an X11 server to
      locally succeed but remotely fail, as headless continuous integration
      (CI) typically has no access to an X11 server. Unsetting this variable
      ensures orthogonality between these cases by coercing the former to fail
      as well.
    '''

    # Defer heavyweight imports.
    from betse.util.os.shell import shellenv

    # Global variables to be set below.
    global EXPORT_SIM_CONF_DIRNAME

    # Globalize the value of the application-specific "--export-sim-conf-dir"
    # command-line option if any for subsequent lookup from module scope.
    #
    # This is required as the "pytest.config" object is no longer safely
    # accessible from module scope. Attempting to do so now results in a
    # deprecation warning resembling:
    #     PytestDeprecationWarning: the `pytest.config` global is deprecated.  Please use `request.config` or `pytest_configure` (if you're a pytest plugin) instead.
    #
    # This ad-hoc circmvention is shamelessly inspired by the following
    # exhaustive StackOverflow treatise on this subject:
    #     https://stackoverflow.com/a/51884507/2809027
    EXPORT_SIM_CONF_DIRNAME = config.getoption('export_sim_conf_dirname')

    #FIXME: This operation should be converted into autouse fixtures defined in
    #this plugin above. Such fixtures should require the builtin fixture
    #permitting us to temporarily change environment variables, which should
    #then be used to temporarily undefine the ${DISPLAY} variable. What was
    #that called again... "monkeypatch"? Contemplate eternity and see:
    #
    #    http://pytest.org/latest/fixture.html#autouse-fixtures-xunit-setup-on-steroids

    # Unset the external `${DISPLAY}` environment variable if currently set.
    # Technically, this operation needs to be performed:
    #
    # * Only once for the entire test suite when py.test is *NOT* parallelized
    #   with "xdist", in which case all tests run in the same process and hence
    #   share the same environment variables.
    # * Once for each test when py.test is parallelized with "xdist", in which
    #   case each test is run in a distinct subprocess and hence does *NOT*
    #   share the same environment variables.
    #
    # Since unsetting environment variables is fast, doing so here
    # transparently supports both use cases detailed above with no discernable
    # downside. See the docstring for additional commentary.
    shellenv.unset_var_if_set('DISPLAY')


def pytest_unconfigure(config) -> None:
    '''
    Hook run immediately *before* exiting the current :mod:`pytest` test
    session.
    '''

    pass

# ....................{ HOOKS ~ test                      }....................
def pytest_runtest_setup(item: 'pytest.main.Item') -> None:
    '''
    Hook run immediately *before* running the passed test.

    Specifically:

    * If this is a **serial test** (i.e., test method bound to an instance of
      the :class:`betse.util.test.pytest.pytabc.SerialTestABC` superclass) for
      which a prior serial test in the same test class was recorded as failing
      by the :func:`pytest_runtest_makereport` hook, this hook marks this test
      as xfailing.

    Parameters
    ----------
    item : pytest.main.Item
        Metadata encapsulating this test callable (e.g., function, method).

    See Also
    ----------
    https://pytest.org/latest/example/simple.html#incremental-testing-test-steps
        Official py.test code snippet inspiring this implementation.
    '''

    # Defer heavyweight imports.
    from betse.exceptions import BetseTestHookException
    from betse.util.test.pytest.pytabc import SerialTestABC

    # For each list of fixtures requested by this test...
    for fixture_defs in item._fixtureinfo.name2fixturedefs.values():
        # For each such fixture...
        for fixture_def in fixture_defs:
            # True only if this fixture is parametrized serially (i.e., is
            # decorated by the @serialize_parametrized_fixture decorator).
            is_fixture_parametrized_serially = getattr(
                fixture_def.func,
                '_betse_is_fixture_parametrized_serially',
                False)
            # print('funcinfo: {}'.format(fixture_def))
            # print('funcinfo.func: {}'.format(fixture_def.func))
            # print('funcinfo.func.is: {}'.format(is_fixture_parametrized_serially))

            # If this fixture is parametrized serially, propagate this setting
            # to this test by dynamically marking this test accordingly.
            if is_fixture_parametrized_serially:
                item.keywords['serialize_parametrized_test'] = True

    # If this test is parametrized serially, do so *AFTER* propagating this
    # setting from fixtures requested by this test above.
    if 'serialize_parametrized_test' in item.keywords:
        # If this test is *NOT* parametrized, raise an exception.
        if not hasattr(item, '_genid'):
            raise BetseTestHookException(
                'Unparametrized test "{}" marked as serial.'.format(item.name))

        # Test callable (i.e., underlying function or method object).
        test_callable = item.obj

        # Unique identifier of the first failing parameters passed to this test
        # if any or "None" otherwise.
        first_failing_param_id = getattr(
            test_callable, '_betse_first_failing_param_id', None)

        # If any parameters passed to this test failed, mark these subsequently
        # passed parameters as also xfailing.
        if first_failing_param_id is not None:
            pytest.xfail(
                'Prior serial test parametrization "{}" failed.'.format(
                    first_failing_param_id))
    #FIXME: This should *REALLY* simply be implemented as a new
    #@serialize_class decorator, much like @serialize_parametrized_test above.

    # If this is an unparametrized test intended to be run serially, do so.
    elif SerialTestABC.is_test_serial(item):
        # Object to which this test method is bound.
        test_instance = item.parent

        # Name of the first previously declared test method to fail in this
        # test method's class if any or "None" otherwise.
        first_failure_method_name = getattr(
            test_instance, '_first_failure_method_name', None)

        # If a previously declared test failed, mark this subsequently declared
        # test as also xfailing.
        if first_failure_method_name is not None:
            pytest.xfail('Prior serial test "{}" failed.'.format(
                test_instance._first_failure_method_name))


def pytest_runtest_makereport(
    item: 'pytest.main.Item', call: 'pytest.runner.CallInfo'):
    '''
    Hook run immediately *after* the passed test returned the passed result.

    Specifically:

    * If this is a **serial test** (i.e., test method bound to an instance of
      the the :class:`betse.util.test.pytest.pytabc.SerialTestABC` superclass)
      that failed, this hook records this failure for subsequent analysis by
      the :func:`pytest_runtest_setup` hook.

    Parameters
    ----------
    item : pytest.main.Item
        Metadata encapsulating this test callable (e.g., function, method).
    call : pytest.runner.CallInfo
        Metadata encapsulating the value returned or exception raised by
        calling this test callable.

    See Also
    ----------
    https://pytest.org/latest/example/simple.html#incremental-testing-test-steps
        Official py.test code snippet inspiring this implementation.
    '''

    # Defer heavyweight imports.
    from betse.util.test.pytest.pytabc import SerialTestABC

    # If this test failed...
    if call.excinfo is not None:
        # If this is a parametrized test intended to be run serially...
        if 'serialize_parametrized_test' in item.keywords:
            # Test callable (i.e., underlying function or method object).
            test_callable = item.obj

            # Record the unique identifiers of these test parameters for
            # subsequent use by the pytest_runtest_setup() hook.
            test_callable._betse_first_failing_param_id = item._genid
        #FIXME: This should *REALLY* simply be implemented as a new
        #@class_serial decorator, much like @parametrize_serial above.
        # If this is an unparametrized test intended to be run serially...
        elif SerialTestABC.is_test_serial(item):
            # Object to which this test method is bound.
            test_parent = item.parent

            # Record the name of this test for subsequent use by the
            # pytest_runtest_setup() hook.
            test_parent._first_failure_method_name = item.name
