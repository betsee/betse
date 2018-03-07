#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
CLI-specific functional tests exporting simulation configurations to external
directories, typically as a subprocess of a larger functional test (e.g., to
validate application backward compatibility).
'''

# ....................{ IMPORTS                            }....................
import pytest
from betse_test.util import requests
from betse_test.util.mark.skip import skip_if

# ....................{ OPTIONS                            }....................
EXPORT_SIM_CONF_DIRNAME = pytest.config.option.export_sim_conf_dirname
'''
Absolute or relative path of the target directory to export (i.e., recursively
copy) each source simulation configuration directory into if any *or* ``None``
otherwise.

See Also
----------
:func:`betse_test.conftest.pytest_addoption`
    Further details.
'''

# ....................{ DECORATORS                         }....................
skip_unless_export_sim_conf_dir = skip_if(
    EXPORT_SIM_CONF_DIRNAME is None,
    reason='CLI option "--export-sim-conf-dir" unpassed.')
'''
Decorator skipping the decorated test if the ``--export-sim-conf-dir``
command-line option specific to this test suite was *not* passed to the
``py.test`` command.
'''

# ....................{ TESTS                              }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# CAUTION: The name of this test is assumed to *NOT* change across Git commits.
# Changing this name would violate this assumption and hence forwards
# compatibility with future versions of this test suite.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

@skip_unless_export_sim_conf_dir
def test_cli_sim_export(
    request: '_pytest.python.FixtureRequest',
    betse_cli_sim: 'CLISimTester',
) -> None:
    '''
    Functional test exporting the simulation configuration directory produced by
    running all BETSE subcommands exporting pickled objects (i.e., ``seed``,
    ``init``, and ``sim``) to the ``test_sim_export`` subdirectory of the
    ``--export-sim-conf-dir`` command-line option passed to the ``py.test``
    command.

    Usage
    ----------
    This pseudo-test is typically run explicitly from the CLI as follows:

    .. code:: bash

        $ ./test -k test_cli_sim_export --export-sim-conf-dir=~/some/where

    Parameters
    ----------
    request : _pytest.python.FixtureRequest
        Builtin fixture describing this test.
    betse_cli_sim : CLISimTester
        Object running BETSE CLI simulation subcommands.
    '''

    # Defer heavyweight imports.
    from betse.util.path import dirs, pathnames

    # Name of the current test.
    test_name = requests.get_tested_name(request)

    # Absolute or relative pathname of the target directory to export into,
    # canonicalized for caller convenience (e.g., to ensure tilde expansion) and
    # created if this directory does *NOT* already exist.
    trg_sim_conf_parent_dirname = dirs.canonicalize_and_make_unless_dir(
        EXPORT_SIM_CONF_DIRNAME)

    # Absolute or relative pathname of the test-specific target directory to
    # copy the test-specific source simulation configuration directory to.
    trg_sim_conf_dirname = pathnames.join(
        trg_sim_conf_parent_dirname, test_name)

    # If this test-specific target directory already exists, raise an exception.
    dirs.die_if_dir(trg_sim_conf_dirname)

    # Absolute pathname of the test-specific source simulation configuration
    # directory to copy from.
    src_sim_conf_dirname = betse_cli_sim.sim_state.conf_dirname

    # Test all subcommands exporting pickled objects. Note that plotting
    # subcommands export no such objects and hence are intentionally omitted.
    betse_cli_sim.run_subcommands_sim()

    # Export (i.e., recursively copy) this source to target directory.
    dirs.copy(
        src_dirname=src_sim_conf_dirname,
        trg_dirname=trg_sim_conf_dirname)
