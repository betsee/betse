#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
External command fixture classes.
'''

# ....................{ IMPORTS                            }....................
from contextlib import ExitStack
from betse.util.command import exits
from betse.util.type import types

# ....................{ CLASSES                            }....................
#FIXME: Implement me!
class CLIMultiTester(object):
    pass

#FIXME: Rename to "CLITester" for conciseness.
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
