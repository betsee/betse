#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
CLI-specific functional tests exporting simulation configurations to external
directories, typically as a subprocess of a larger functional test (e.g., to
validate backwards compatibility).
'''

# ....................{ IMPORTS                            }....................
import pytest
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
    reason='CLI option "--export-sim-conf-dir" not passed')
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

#FIXME: Add support to "betse_setup/test.py" for a new "--pytest-options" option
#simply passing its argument as additional options to "py.test".
@skip_unless_export_sim_conf_dir
def test_sim_export(betse_cli_sim: 'CLISimTester') -> None:
    '''
    Export the simulation configuration directory produced by running all BETSE
    subcommands exporting pickled objects (i.e., ``seed``, ``init``, and
    ``sim``) to the ``test_sim_export`` subdirectory of the
    ``--export-sim-conf-dir`` command-line option passed to the ``py.test``
    command.

    Usage
    ----------
    This pseudo-test is typically run explicitly from the CLI as follows:

    .. code:: bash

        $ ./test -k test_sim_export --pytest-options="--export-sim-conf-dir=~/some/where"

    Parameters
    ----------
    betse_cli_sim : CLISimTester
        Object running BETSE CLI simulation subcommands.
    '''

    # Defer heavyweight imports.
    from betse.util.path import dirs, pathnames

    # Test all subcommands exporting pickled objects. Note that plotting
    # subcommands export no such objects and hence are intentionally omitted.
    betse_cli_sim.run_subcommands(('seed',))
    # betse_cli_sim.run_subcommands(('seed',), ('init',), ('sim',))

    # Absolute or relative pathname of the test-specific target directory to
    # copy the test-specific source simulation configuration directory to.
    trg_sim_conf_dirname = pathnames.join(
        EXPORT_SIM_CONF_DIRNAME, 'test_sim_export')

    # Export (i.e., recursively copy) this source to target directory.
    dirs.copy(
        src_dirname=betse_cli_sim.sim_state.conf_filename,
        trg_dirname=trg_sim_conf_dirname)
