#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
CLI-specific functional tests exercising all miscellaneous subcommands (e.g.,
`betse info`), thus excluding those specific to networks and simulations.
'''

# ....................{ IMPORTS                           }....................
import pytest
from betse.util.test.pytest.mark.pytskip import (
    skip_unless_lib_runtime_optional)

# ....................{ TESTS                             }....................
def test_cli_no_arg(betse_cli: 'CLITester') -> None:
    '''
    Test the ``betse`` command passed no arguments.

    Parameters
    ----------
    betse_cli : CLITester
        Object encapsulating the BETSE CLI.
    '''

    betse_cli.run()


#FIXME: Something has gone horrifyingly wrong. For unknown reasons, "py.test"
#appears to have removed the standard "${DISPLAY}" environment variable. We
#suspect that the "pytest-qt" plugin is the culprit, implying that our external
#test driver (e.g., "betse.lib.setuptool.commands.supcmdtest") will need to
#explicitly disable this plugin whenever present.
#
#Further testing is warranted. As an initial simplistic test, consider simply
#uninstalling "pytest-qt" and retest us up.
def test_cli_info(betse_cli: 'CLITester') -> None:
    '''
    Test the ``betse info`` subcommand.

    Parameters
    ----------
    betse_cli : CLITester
        Object encapsulating the BETSE CLI.
    '''

    from betse.util.os.shell import shellenv
    print('Environment variables:\n{}'.format(shellenv.to_str()))
    betse_cli.run('info')


def test_cli_config(
    betse_cli: 'CLITester',
    betse_temp_dir: 'LocalPath',
) -> None:
    '''
    Test the ``betse config`` subcommand, creating a configuration in a
    temporary directory specific to this test.

    While other tests requesting any fixture whose name is prefixed by
    ``betse_sim_config_`` (e.g., the ``betse_sim_config_default`` fixture)
    already create such configurations, they do so by directly calling a
    comparatively low-level function of the BETSE API rather than by running a
    comparatively high-level subcommand of the BETSE CLI. Since these two logic
    paths are different (albeit related), this test ensures the latter to be
    exercised.

    Parameters
    ----------
    betse_cli : CLITester
        Object encapsulating the BETSE CLI.
    betse_temp_dir : LocalPath
        Object encapsulating a temporary directory isolated to the current
        test.
    '''

    # Absolute path of this configuration file in this temporary directory.
    sim_config_filepath = betse_temp_dir.join('sim_config.yaml')

    # Create this file. (Avoid shell-quoting this path. Doing so unnecessarily
    # adds an additional level of quoting... which is bad.)
    betse_cli.run('config', str(sim_config_filepath))

    # Assert this file to have been created.
    assert sim_config_filepath.check(file=1)

# ....................{ TESTS ~ parametrized              }....................
@pytest.mark.parametrize(
    ('profile_type',), (
        # Profile types leveraging stdlib functionality guaranteed to exist.
        ('none',),
        ('call',),

        # Profile types leveraging third-party packages necessitating checking.
        pytest.param(
            'size', marks=skip_unless_lib_runtime_optional('pympler')),
    ),
)
def test_cli_profile(
    betse_cli: 'CLITester',
    betse_temp_dir: 'LocalPath',
    profile_type: str,
) -> None:
    '''
    Test profiling of an arbitrary ``betse`` command with the parametrized
    profiling type, serialized to a temporary file specific to this test.

    Parameters
    ----------
    betse_cli : CLITester
        Object encapsulating the BETSE CLI.
    betse_temp_dir : LocalPath
        Object encapsulating a temporary directory isolated to the current
        test.
    profile_type : str
        Type of profiling to perform -- either:

        * ``none``, performing no profiling.
        * ``call``, performing call-granularity profiling.
        * ``size``, performing memory profiling.
    '''
        # * `line`, performing line-granularity profiling.

    # CLI option enabling this type of profiling.
    profile_type_option = '--profile-type={}'.format(profile_type)

    # If this is a profiler outputing no profile file...
    if profile_type in ('none', 'size',):
        # Profile the basic "betse" command to only the logfile. While the
        # absolute path of a profile file could also be passed, doing so
        # requires inefficiently creating a temporary directory. Since these
        # profilers output no profile file, this directory would remain empty.
        betse_cli.run(profile_type_option)
    else:
        # Absolute path of the profile file isolated to a test-specific
        # temporary directory output by the following command.
        profile_filepath = betse_temp_dir.join('profile.prof')

        # Profile the basic "betse" command to the logfile *AND* this file.
        betse_cli.run(
            profile_type_option, '--profile-file={}'.format(profile_filepath))

        # Assert this file to have been created.
        assert profile_filepath.check(file=1)
