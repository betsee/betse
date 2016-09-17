#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Functional tests for BETSE's CLI testing all subcommands _except_ for
simulation-specific subcommands (e.g., `betse info`).
'''

# ....................{ TESTS                              }....................
def test_cli_no_arg(betse_cli: 'CLITester') -> None:
    '''
    Test the `betse` command passed no arguments.

    Parameters
    ----------
    betse_cli : CLITester
        Test-specific object encapsulating the BETSE CLI.
    '''

    betse_cli.run()


def test_cli_info(betse_cli: 'CLITester') -> None:
    '''
    Test the `betse info` subcommand.

    Parameters
    ----------
    betse_cli : CLITester
        Test-specific object encapsulating the BETSE CLI.
    '''

    betse_cli.run('info')


def test_cli_config(
    betse_cli: 'CLITester',
    tmpdir_factory: '_pytest.tmpdir.tmpdir_factory',
) -> None:
    '''
    Test the `betse config` subcommand, creating a configuration in a temporary
    directory specific to this test.

    While other tests requesting any fixture whose name is prefixed by
    `betse_sim_config_` (e.g., the `betse_sim_config_default` fixture) already
    create such configurations, they do so by directly calling a comparatively
    low-level function of the BETSE API rather than by running a comparatively
    high-level subcommand of the BETSE CLI. Since these two logic paths are
    different (albeit related), this test ensures the latter to be exercised.

    Parameters
    ----------
    betse_cli : CLITester
        Test-specific object encapsulating the BETSE CLI.
    tmpdir_factory : _pytest.tmpdir.tmpdir_factory
        Builtin session-scoped fixture whose `mktemp()` method returns a
        `py.path.local` instance encapsulating a new temporary directory.
    '''

    # Create this temporary directory and wrap this directory's absolute path
    # with a high-level "py.path.local" object.
    sim_config_dirpath = tmpdir_factory.mktemp('config')

    # Absolute path of this configuration file in this temporary directory.
    sim_config_filepath = sim_config_dirpath.join('sim_config.yaml')

    # Create this file. (Avoid shell-quoting this path. Doing so unnecessarily
    # adds an additional level of quoting... which is bad.)
    betse_cli.run('config', str(sim_config_filepath))

    # Assert this file to have been created.
    assert sim_config_filepath.check(file=1)
