#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
External command fixture classes.
'''

# ....................{ IMPORTS                            }....................
from contextlib import ExitStack
from betse.util.path.command import exits
from betse.util.type.types import type_check, SequenceTypes

# ....................{ CLASSES ~ single                   }....................
class CLITester(object):
    '''
    BETSE CLI test runner, efficiently testing a single subcommand of the
    official BETSE CLI (i.e., `betse`) in the active Python interpreter.

    Simple functional fixtures (e.g., `betse_cli`) typically return instances of
    this class to other fixtures and tests exercising a single facet of the
    BETSE CLI.

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
        List of all context managers with which to subsequently call the
        `betse.cli.__main__.run()` method when this object's `run()` method is
        called.
    '''

    @type_check
    def __init__(self, contexts: SequenceTypes) -> None:
        '''
        Initialize this test runner with the passed `request` fixture object.

        Parameters
        ----------
        contexts : SequenceType
            List of all context managers with which to subsequently call the
            `betse.cli.__main__.main()` method when this object's `run()` method
            is called.
        '''

        # Classify the passed parameters.
        self._contexts = contexts


    def __call__(self, *args) -> None:
        '''
        Run the BETSE CLI with the passed positional arguments.

        This special method is a convenience permitting this fixture to be
        called as is rather than via the `run()` method.

        See Also
        ----------
        run
            Method to which this special method internally defers.
        '''

        return self.run(*args)


    def run(self, *args) -> None:
        '''
        Run the BETSE CLI with the passed positional arguments.

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

# ....................{ CLASSES ~ multi                    }....................
class CLITesterPreArged(object):
    '''
    BETSE CLI test runner, efficiently testing a single subcommand of the
    official BETSE CLI (i.e., `betse`) in the active Python interpreter with an
    argument list passed to this object's `__init__()` rather than `run()`
    method.

    Complex functional fixtures (e.g., `betse_cli_sim`) typically return
    instances of this class to other fixtures and tests exercising a predefined
    facet of the BETSE CLI.

    Attributes
    ----------
    _cli : CLITester
        BETSE CLI test runner, testing a single subcommand of the official
        BETSE CLI (i.e., `betse`) in the active Python interpreter.
    _subcommand_args : SequenceType
        Argument list comprising the BETSE CLI subcommand to be tested.

    See Also
    ----------
    CLITestRunner
        Further details on BETSE CLI execution.
    '''


    @type_check
    def __init__(
        self,
        cli: CLITester,
        subcommand_args: SequenceTypes,
    ) -> None:
        '''
        Initialize this test runner with the passed `request` fixture object.

        Parameters
        ----------
        cli : CLITester
            BETSE CLI test runner, testing a single subcommand of the official
            BETSE CLI (i.e., `betse`) in the active Python interpreter.
        subcommand_args : collections.Sequence
            Argument list to be subsequently passed to the `cli.run()` method
            when this object's `run()` method is called. These arguments should
            comprise the BETSE CLI subcommand to be tested.
        '''

        # Classify the passed parameters.
        self._cli = cli
        self._subcommand_args = subcommand_args


    def __call__(self, *args) -> None:
        '''
        Run the BETSE CLI with the argument list previously passed to the
        `__init__()` method.

        This special method is a convenience permitting this fixture to be
        called as is rather than via the `run()` method.

        See Also
        ----------
        run
            Method to which this special method internally defers.
        '''

        return self.run(*args)


    def run(self) -> None:
        '''
        Run the BETSE CLI with the argument list previously passed to the
        `__init__()` method.
        '''

        self._cli(*self._subcommand_args)
