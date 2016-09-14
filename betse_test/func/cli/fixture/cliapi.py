#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
External command fixture classes.
'''

# ....................{ IMPORTS                            }....................
from betse.util.path.command import exits
from betse.util.type.types import type_check, SequenceTypes
from betse_test.func.fixture.sim.configapi import SimTestState
from contextlib import ExitStack

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


    #FIXME: Excise the "contexts" argument and "_contexts" attribute.
    @type_check
    def __init__(self, contexts: SequenceTypes) -> None:
        '''
        Initialize this test runner.

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
        args : tuple
            Tuple of zero or more arguments to be passed to this entry point,
            corresponding exactly to the set of command-line arguments accepted
            by the external command for the BETSE CLI (i.e., `betse`).

        See Also
        ----------
        `betse --help`
            Further details on arguments accepted by the BETSE CLI.
        '''

        # Defer heavyweight imports to their point of use.
        from betse.cli.__main__ import main

        # List of arguments:
        #
        # * Converted from this tuple of arguments.
        # * Prefixed by failure-friendly options.
        arg_list = ['--verbose', '--log-type=none'] + list(args)
        # print('BETSE arg list: {}'.format(arg_list))

        #FIXME: Reduce this entire block to simply:
        #
        #    exit_status = main(arg_list)
        #
        #Context management is now handled elsewhere.

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

# ....................{ CLASSES ~ sim                      }....................
class CLISimTester(object):
    '''
    BETSE CLI simulation test runner, exercising multiple subcommands of the
    BETSE CLI (i.e., `betse`) in the active Python interpreter with the same
    temporary simulation configuration.

    Complex functional fixtures (e.g., `betse_cli_sim`) typically return
    instances of this class to other fixtures and tests exercising multiple
    facets of the BETSE CLI.

    Attributes
    ----------
    _cli : CLITester
        BETSE CLI test runner, testing a single subcommand of the official
        BETSE CLI (i.e., `betse`) in the active Python interpreter.
    _sim_state: SimTestState
        Test-specific object encapsulating a temporary simulation
        configuration file specific to the current test.

    See Also
    ----------
    :class:`CLITester`
        Further details on BETSE CLI execution.
    '''


    @type_check
    def __init__(
        self,
        cli: CLITester,
        sim_state: SimTestState,
    ) -> None:
        '''
        Initialize this test runner.

        Parameters
        ----------
        cli : CLITester
            BETSE CLI test runner, testing a single subcommand of the official
            BETSE CLI (i.e., `betse`) in the active Python interpreter.
        sim_state: SimTestState
            Test-specific object encapsulating a temporary simulation
            configuration file specific to the current test.
        '''

        # Classify the passed parameters.
        self._cli = cli
        self._sim_state = sim_state


    def run_subcommands(self, *subcommands_args: SequenceTypes) -> None:
        '''
        Run all BETSE CLI subcommands signified by the passed argument lists in
        the active Python process (_in the passed order_).

        To ensure that each such subcommand runs the the same simulation, this
        method implicitly:

        * Appends each such argument list by the absolute path of the simulation
          configuration file with which this test runner was initialized.
        * Temporarily changes the current working directory (CWD) to the
          directory containing this file.

        Parameters
        ----------
        subcommands_args : tuple
            Tuple of sequences of **subcommand arguments** (i.e., one or more
            shell words comprising the BETSE CLI subcommand to be tested).
        '''

        for subcommand_args in subcommands_args:
            self.run_subcommand(*subcommand_args)


    def run_subcommand(self, *subcommand_args: str) -> None:
        '''
        Run the BETSE CLI subcommand signified by the passed argument list in
        the active Python process.

        This method implicitly:

        * Appends this argument list by the absolute path of the simulation
          configuration file with which this test runner was initialized.
        * Temporarily changes the current working directory (CWD) to the
          directory containing this file.

        Parameters
        ----------
        subcommand_args : tuple
            Tuple of **subcommand arguments** (i.e., one or more shell words
            comprising the BETSE CLI subcommand to be tested).
        '''

        # Temporarily change the CWD to this simulation file's directory.
        with self._sim_state.get_command_context():
            # String preceding and following the argument list printed below.
            banner_str = '=' * 30

            # Notify users of this argument list.
            print('\n{} {} {}'.format(
                banner_str, ' '.join(subcommand_args), banner_str,))

            # Append the absolute path of this runner's configuration file to
            # the passed tuple of arguments. While inefficient, converting this
            # tuple into a list would be even more inefficient.
            subcommand_args += (self._sim_state.config.filename,)

            # Run this subcommand.
            self._cli(*subcommand_args)

# ....................{ CLASSES ~ obsolete                 }....................
#FIXME: Excise.
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
