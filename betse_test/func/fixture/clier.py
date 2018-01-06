#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Fixtures and fixture classes efficiently exercising a single subcommand of the
BETSE CLI in the active Python interpreter.
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
#FIXME: This is great. We'd only like to make one improvement to the
#aforementioned PyInstaller testing support: marks. Use them. Explicit is
#typically better than implicit. In particular:
#
#* Define a new "betse_test.mark.able" submodule.
#* Define a new "@freezable" mark decorator in this submodule.
#* Decorate all tests suitable for running under a frozen copy of BETSE with
#  this decorator. This should probably *ONLY* include:
#  * All functional CLI tests, excluding "test_pyi_try". Everything else (e.g.,
#    unit tests) are *NOT* suitable for running frozen.
#FIXME: Urgh. That's great, but it's still not quite enough. Why? Simulation
#configurations. We currently create them by importing BETSE packages, which
#clearly won't work for frozen commands. Instead, if running frozen, we'll need
#to just call "betse config". Simple, if tedious.

# ....................{ IMPORTS                            }....................
from betse.util.path.command import cmdexit
from betse.util.type.types import type_check
from pytest import fixture

# ....................{ CONSTANTS                          }....................
_CLI_OPTIONS_MANDATORY = (
    '--verbose',
    '--log-level=none',
    '--matplotlib-backend=agg',
)
'''
Tuple of all failure-friendly command-line options unconditionally passed to all
invocations of the BETSE CLI by functional tests.

To improve debuggability for failing tests, these options are unconditionally
passed by :class:`CLITester` instances created by the :func:`betse_cli` fixture:

* ``--verbose``, logging low-level debugging messages to stdout, which `py.test`
  captures for all tests and displays for all failing tests.
* ``--log-level=none``, redirecting all log messages to either stdout or stderr
  but *not* a logfile. While ``py.test`` can be configured to capture logfile
  messages, doing so sanely is complicated by the fact that ``py.test`` already
  captures both stdout and stderr by default. As there is no benefit in
  recapturing logfile messages already logged to either stdout or stderr, tests
  avoid doing so entirely.
* ``--matplotlib-backend=agg``, enabling the default non-interactive matplotlib
  backend *guaranteed* to be usable on all platforms. By default, matplotlib
  enables an interactive backend (e.g., ``Qt5Agg``) inhibiting sane test
  automation.
'''

# ....................{ CLASSES                            }....................
class CLITester(object):
    '''
    BETSE CLI test runner, efficiently testing a single subcommand of the
    official BETSE CLI (i.e., `betse`) in the active Python interpreter.

    Simple functional fixtures (e.g., `betse_cli`) typically return instances of
    this class to other fixtures and tests exercising a single facet of the
    BETSE CLI.

    Command Execution
    ----------
    For both efficiency and reliably, this runner does _not_ actually fork this
    subcommand as a separate process. Doing so would introduce installation
    complications, portability concerns, and non-reproduceable edge-cases (e.g.,
    conflicting versions of the same `betse` command in the current `${PATH}`).
    Instead, this runner:

    * Imports the :mod:`betse.cli.__main__` module implementing the BETSE CLI.
    * Passes this module's :func:`betse.cli.__main__.run()` function the passed
      arguments extended by the mandatory arguments defined by the
      :data:`_CLI_OPTIONS_MANDATORY` tuple global.
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(self) -> None:
        '''
        Initialize this test runner.
        '''

        pass  # Well, that was easy.

    # ..................{ RUNNERS                            }..................
    @type_check
    def run(self, *args: str) -> None:
        '''
        Run the BETSE CLI with the passed positional string arguments, extended
        by the mandatory positional string arguments defined by the
        :data:`_CLI_OPTIONS_MANDATORY` tuple global.

        Parameters
        ----------
        args : Tuple[str]
            Tuple of zero or more arguments to be passed to this entry point,
            corresponding exactly to the set of command-line arguments accepted
            by the external command for the BETSE CLI (i.e., `betse`).

        See Also
        ----------
        `betse --help`
            Further details on arguments accepted by the BETSE CLI.
        '''

        # Defer heavyweight imports to their point of use.
        from betse.__main__ import main

        # Prefixed this argument list by failure-friendly options.
        args_evolved = _CLI_OPTIONS_MANDATORY + args
        # print('BETSE arg list: {}'.format(arg_list))

        # Run the BETSE CLI subcommand corresponding to these arguments,
        # capturing the exit status of that subcommand for testing.
        exit_status = main(args_evolved)

        # If this exit status signifies failure, fail the current test.
        assert cmdexit.is_success(exit_status), (
            'BETSE CLI failed with exit status {} '
            'given arguments: {}'.format(exit_status, args_evolved))

# ....................{ FIXTURES                           }....................
# Test-scope fixture creating and returning a new object for each discrete test.
@fixture
def betse_cli() -> CLITester:
    '''
    Fixture returning a test-specific object suitable for running arbitrary
    BETSE CLI subcommands required by the current fixture or test.

    Returns
    ----------
    CLITester
        Object running arbitrary BETSE CLI subcommands.
    '''

    # Print a newline, ensuring that the first line of stderr or stdout printed
    # by the BETSE CLI is visually offset from current py.test output.
    print()

    # Create and return a new BETSE CLI test runner.
    return CLITester()
