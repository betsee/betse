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
from contextlib import ExitStack
from betse.util.command import exits
from betse.util.type import types, objects
from betse_test.util import requests
from pytest import fixture

# ....................{ CLASSES                            }....................
class CLITestRunner(object):
    '''
    BETSE interface test runner, efficiently testing the external command for
    either the official BETSE CLI (e.g., `betse`) in the active Python
    interpreter.

    Functional test fixtures typically return instances of this class to other
    functional test fixtures and tests testing BETSE's CLI.

    Command Execution
    ----------
    For both efficiency and reliably, this runner does _not_ actually execute
    this command. Doing so would introduce installation complications and
    portability concerns (e.g., conflicting versions of the `betse` command in
    the current `${PATH}`). Instead, this runner:

    * Imports the `betse.cli.__main__` module implementing the BETSE CLI.
    * Passes this module's `run()` method the passed arguments.

    Attributes
    ----------
    contexts : list
        List of all context managers with which to call the
        `betse.cli.__main__.run()` method when this object's `run()` method is
        called.
    '''


    # corresponding to the context managers returned by the optional
    # get_command_context() methods defined by the instances of these fixtures.
    def __init__(self, contexts: list) -> None:
        '''
        Initialize this test runner with the passed `request` fixture object.

        Parameters
        ----------
        contexts : list
            List of all context managers under which to subsequently call the
            `betse.cli.__main__.run()` method when this object's `run()` method
            is called.
        '''
        assert types.is_sequence_nonstr(contexts), (
            types.assert_not_sequence_nonstr(contexts))

        self._contexts = contexts


    def __call__(self, *args) -> None:
        '''
        Call the entry point for this BETSE interface with the passed positional
        arguments.

        This special method is a convenience permitting this fixture to be
        called as is rather than via the `run()` method.

        See Also
        ----------
        `run()`
            For further details, which this special method internally defers to.
        '''

        return self.run(*args)


    def run(self, *args) -> None:
        '''
        Call the entry point for this BETSE interface with the passed positional
        arguments.

        To improve debuggability for failing tests, this function
        unconditionally passes these command-line options to this interface:

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
            Further details on such arguments.
        '''

        # Defer heavyweight imports to their point of use.
        from betse.cli.__main__ import main

        # List of arguments:
        #
        # * Converted from this tuple of arguments.
        # * Prefixed by failure-friendly options.
        arg_list = ['--verbose', '--log-type=none'] + list(args)
        # print('BETSE arg list: {}'.format(arg_list))

        #FIXME: If the current test also requires a fixture whose name is
        #prefixed by "betse_sim_config_", the following call to main() should be
        #embedded in a "with"-style block calling paths.change_current():
        #
        #1. Refactor the "betse_cli" fixture to do the following *BEFORE*
        #   instantiating and returning this class:
        #   1. Accept a "request" fixture.
        #   2. Search the "request.fixturenames" list for all fixtures whose
        #      names are prefixed by "betse_sim_config_". There should only be
        #      exactly one. This logic already resides in the
        #      configbase._betse_sim_config() fixture, suggesting we generalize
        #      that logic to a new test utility getter returning the name of
        #      this single fixture.
        #   3. Call request.getfuncargvalue(fixture_name) to obtain the
        #      "SimTestConfig" instance returned by that fixture.
        #2. Pass this instance to CLITestRunner.__init__().
        #3. Refactor CLITestRunner.__init__() to classify the passed
        #   "SimTestConfig" instance if any as a private attribute. While None
        #   is an acceptable value, this parameter should *NOT* be optional.
        #4. If that attribute is non-None, embed this call in a "with"
        #   statement changing to the directory containing this attribute's
        #   configuration file.
        #
        #The alternative, of course, would be to permanently change the CWD for
        #the entirety of this test session in the _betse_sim_config fixture.
        #Doing so strikes us as a bad idea.
        #FIXME: See also the builtin "monkeypatch" fixture, which provides
        #chdir() and undo() methods -- the latter of which *SHOULD* be
        #implicitly called on teardown. However, is there any point? Our
        #context-manager is almost certainly safer and already exists.
        #FIXME: Ah. Actually, the simplest (and almost certainly most correct)
        #way comes to mind:
        #
        #* Refactor the "betse_cli" fixture into a new "betse_command" fixture
        #  calling *EITHER* the "betse" CLI or appropriate GUI command
        #  conditionally based on the name of the calling test. Yup! By far the
        #  simplest way. If the calling test's name is prefixed by:
        #  * "test_cli_", then call the CLI.
        #  * "test_gui_", then call the GUI.
        #  * Else, raise an exception.
        #* Refactor this class so support such conditionality.
        #* Move this module to "betse_test/func/fixture/command.py".
        #* Refactor the _betse_sim_config fixture to *ALWAYS* require the new
        #  "betse_command" fixture, ensuring the former has access to the
        #  current instance of this class.
        #* Define a new "_current_dirname = None" attribute in this class'
        #  __init__() method.
        #* Define a new set_current_dirname() method on this class, setting this
        #  attribute to the passed non-empty string.
        #* If that attribute is non-None here, use a "with" statement as
        #  outlined above.
        #
        #Done! Pretty awesome, actually.
        #FIXME: O.K.; the above conception would work, but has a few critical
        #design issues. The most significant is that:
        #
        #* It doesn't generalize to multiple transformations. How would we add
        #  support for additional conditional requirements on behalf of other
        #  optional fixtures required by tests?
        #* The "with" statement logic should, logically speaking, reside in the
        #  class that is actually concerned with that logic: namely,
        #  "SimTestConfig".
        #
        #Happily, we can hit both stones with one bird as follows:
        #
        #* Search "request.fixturenames" for all fixtures with names prefixed by
        #  "betse_". Since this repeats logic in "sim/configbase.py", generalize
        #  this into a new function in "betse_test.util.requests".
        #* Within a "contextlib.ExitStack" context manager:
        #  * For each matched fixture name (ignoring the name of the "betse_cli"
        #    fixture, obviously):
        #    * Obtain each such fixture's object by calling
        #      request.getfuncargvalue() with that fixture's name. Again,
        #      generalize this logic into a new function.
        #    * If that fixture object defines a get_command_context() method:
        #      * Call that method and pass that method's return value to
        #        exit_stack.enter_context().
        #  * Call main() as below.
        #
        #Surprisingly, this is actually *EASIER* than the approach outlined
        #above. No new properties are required or coordination between multiple
        #fixtures. Indeed, the new approach outlined here is far more general.

        # Create and enter a context manager of context managers: that is, a
        # context manager permitting a variable number of other context managers
        # to be dynamically entered.
        with ExitStack() as stack:
            # Enter each context manager passed to this object's __init__()
            # method, equivalent to specifying each such manager in a "with"
            # statement in a nested manner.
            for context in self._contexts:
                stack.enter_context(context)

            # Exit status of BETSE's CLI passed these arguments.
            exit_status = main(arg_list)

        # If this exit status connotes failure, fail the current test.
        assert exits.is_success(exit_status), (
            'BETSE CLI failed with exit status {} '
            'given argument list {}.'.format(exit_status, arg_list))

# ....................{ FIXTURES                           }....................
# To force this fixture to return a new object for all parent fixtures and
# tests, this fixture is declared with default scope (i.e., test).
@fixture
def betse_cli(request: '_pytest.python.FixtureRequest') -> CLITestRunner:
    '''
    Fixture returning an instance of the `CLITestRunner` class, suitable for
    running the BETSE CLI command required by the current fixture or test.

    Parameters
    ----------
    request : _pytest.python.FixtureRequest
        Builtin fixture parameter describing the parent fixture or test of this
        fixture (and similar contextual metadata).

    Returns
    ----------
    CLITestRunner
        Object running the BETSE CLI command.
    '''

    # Names of all BETSE-specific fixtures required by the current test,
    # excluding the current fixture.
    betse_fixture_names = requests.get_fixture_names_prefixed_by(
        request=request, fixture_name_prefix='betse_')

    # List of the command-specific context managers provided by these fixtures,
    # corresponding to the context managers returned by the optional
    # get_command_context() methods defined by the instances of these fixtures.
    betse_fixture_command_contexts = []

    # For the name of each such fixture...
    for betse_fixture_name in betse_fixture_names:
        # Fixture object returned by this fixture.
        betse_fixture = requests.get_fixture(
            request=request, fixture_name=betse_fixture_name)

        # get_command_context() method defined by this fixture object if any.
        betse_fixture_get_command_context = objects.get_method_or_none(
            obj=betse_fixture, method_name='get_command_context')

        # If this fixture object defines this method, append the context manager
        # returned by call this method to this list.
        if betse_fixture_get_command_context is not None:
            betse_fixture_command_contexts.append(
                betse_fixture_get_command_context())

    # Return a new CLI runner specific to the current test.
    return CLITestRunner(contexts=betse_fixture_command_contexts)
