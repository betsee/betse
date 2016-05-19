#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
External command fixtures.

These fixtures automate testing of the current versions of BETSE's externally
runnable CLI and GUI commands (e.g., `betse`, `betse-qt4`) in subprocesses of
the current test process, regardless of whether these commands have been
editably installed (i.e., as synchronized symlinks rather than desynchronized
copies) into the current Python environment or not.
'''

#FIXME: We've added a new "--dist-dir" option to our "freeze_*" family of
#setuptools subcommands. Now, let's add support for fixtures running these
#subcommands. Is there any means of programmatically running setuptools
#subcommands in the current Python process? We suspect not, sadly.

#FIXME: Add seemless support for exercising all CLI-specific functional tests
#leveraging the "betse_cli" fixture against a frozen rather than unfrozen
#version of BETSE's CLI. Do *NOT* parametrize this fixture to unconditionally
#run each such test against both frozen and unfrozen versions of BETSE's CLI, as
#doing so would substantially increase space and time complexity with little to
#no tangible gain. Instead:
#
#* Define a new "--is-frozen" command-line option specific to our setuptools
#  "test" subcommand, defaulting to disabled.
#* When this option is *NOT* passed:
#  * One and only one PyInstaller test should be run -- say, test_pyi_try().
#    This test is intended only to provide a coarse-grained sanity check of
#    whether or not:
#    * BETSE is actually freezable under PyInstaller.
#    * The resulting executable is successfully simulatable with the default
#      simulation configuration via "betse try".
#  * Consequently, this test should (in order):
#    1. Freeze BETSE in one-dir mode (which requires no additional compression
#       and decompression steps and hence is slightly more time efficient for
#       our simple purposes) into a temporary subdirectory.
#    2. Run "betse try" against this frozen executable.
#  * All other CLI tests should be run as is against the unfrozen version of
#    BETSE importable by the active Python interpreter.
#* When this option is passed:
#  * The aforementioned test_pyi_try() test should *NOT* be run, as doing so
#    would be wholly redundant.
#  * BETSE should be frozen in one-dir mode, as discussed above.
#  * All other CLI tests should be run against this frozen executable rather
#    than the unfrozen, importable version of BETSE.
#
#Hence, the "--is-frozen" CLI option serves as a high-level switch fundamentally
#changing testing behaviour. The intention, of course, is that this option would
#only be passed immediately before producing an official new frozen version of
#BETSE. This provides a sanity check on frozen executable behaviour otherwise
#difficult (if not infeasible) to do manually.

# ....................{ IMPORTS                            }....................
from pytest import fixture

# ....................{ CLASSES                            }....................
class CLITestRunner(object):
    '''
    CLI-specific test runner, efficiently executing the external command for the
    BETSE CLI (i.e., `betse`) in the active Python interpreter.

    For both efficiency and reliably, this runner does _not_ actually execute
    this command. Doing so would introduce installation complications and
    portability concerns (e.g., conflicting versions of the `betse` command in
    the current `${PATH}`). Instead, this runner:

    * Imports the `betse.cli.__main__` module implementing the BETSE CLI.
    * Passes this module's `run()` method the passed arguments.

    CLI-specific fixtures typically return instances of this class as a
    means of communicating this runner to other fixtures and tests.
    '''

    def __call__(self, *args) -> None:
        '''
        Call the entry point for BETSE's CLI with the passed positional
        arguments.

        This special method is a convenience permitting this fixture to be
        called as is rather than via the `run()` method.

        See Also
        ----------
        `run()`
            For further details, which this special method internally defers to.
        '''

        return self.run(*args)


    @staticmethod
    def run(*args) -> None:
        '''
        Call the entry point for BETSE's CLI with the passed arguments.

        To improve debuggability for failing tests, this function
        unconditionally passes these options to this call as well:

        * `--verbose`, logging low-level debugging messages to stdout, which
          `py.test` captures for all tests and displays for all failing tests.
        * `--log-type=none`, redirecting all log messages to either stdout or
          stderr but _not_ a logfile. While `py.test` can be configured to
          capture logfile messages, doing so sanely is complicated by the fact
          that `py.test` already captures both stdout and stderr by default.  As
          there is no benefit in recapturing logfile messages already logged to
          either stdout or stderr, we avoid doing so entirely.

        Parameters
        ----------
        args : list
            List of zero or more arguments to be passed to this entry point,
            corresponding exactly to the set of command-line arguments accepted
            by the external command for the BETSE CLI (i.e., `betse`).

        See Also
        ----------
        `betse --help`
            Further details on these arguments.
        '''

        # Defer heavyweight imports to their point of use.
        from betse.cli.__main__ import main
        from betse.util.command import exits

        # List of arguments:
        #
        # * Converted from this tuple of arguments.
        # * Prefixed by failure-friendly options.
        arg_list = ['--verbose', '--log-type=none'] + list(args)

        # Exit status of the entry point for BETSE's CLI passed these arguments.
        exit_status = main(arg_list)

        # If this exit status connotes failure, fail the caller test or fixture.
        assert exits.is_success(exit_status), (
            'BETSE CLI failed with exit status {} '
            'given argument list {}.'.format(exit_status, arg_list))

# ....................{ FIXTURES ~ low-level               }....................
@fixture(scope='session')
def betse_cli() -> CLITestRunner:
    '''
    Fixture returning a singleton instance of the `CLITestRunner` class.

    For efficiency, this instance is shared by all invocations of this fixture
    for the current test session.
    '''

    return CLITestRunner()
