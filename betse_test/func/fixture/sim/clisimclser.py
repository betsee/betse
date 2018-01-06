#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Fixture classes running BETSE CLI subcommands in the active Python interpreter.
'''

# ....................{ IMPORTS                            }....................
from betse.util.type.types import type_check, SequenceTypes
from betse_test.fixture.simconf.simconfclser import SimConfTestABC
from betse_test.func.fixture.clier import CLITester

# ....................{ CLASSES                            }....................
class CLISimTester(object):
    '''
    BETSE CLI simulation test runner, exercising multiple subcommands of the
    BETSE CLI (i.e., ``betse``) in the active Python interpreter with the same
    temporary simulation configuration.

    Complex functional fixtures (e.g., :func:`betse_cli_sim`) typically return
    instances of this class to other fixtures and tests exercising multiple
    facets of the BETSE CLI.

    Attributes
    ----------
    cli_tester : CLITester
        BETSE CLI test runner, testing a single subcommand of the official
        BETSE CLI (i.e., ``betse``) in the active Python interpreter.
    sim_state: SimConfTestABC
        Test-specific object encapsulating a temporary simulation
        configuration file specific to the current test.

    See Also
    ----------
    :class:`CLITester`
        Further details on BETSE CLI execution.
    '''

    # ..................{ CONSTANTS                          }..................
    _SUBCOMMANDS_ARGS_DEFAULT = (
        ('seed',),
        ('init',),
        ('sim',),
        ('plot', 'seed'),
        ('plot', 'init'),
        ('plot', 'sim'),
    )
    '''
    Tuple of the argument lists comprising all default simulation-specific BETSE
    CLI subcommands.

    This tuple is a convenience simplifying common-case tests exercising *only*
    these subcommands.
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(
        self, cli_tester: CLITester, sim_state: SimConfTestABC) -> None:
        '''
        Initialize this test runner.

        Parameters
        ----------
        cli_tester: CLITester
            BETSE CLI test runner, testing a single subcommand of the official
            BETSE CLI (i.e., ``betse``) in the active Python interpreter.
        sim_state: SimConfTestABC
            Test-specific object encapsulating a temporary simulation
            configuration file specific to the current test.
        '''

        # Classify these parameters.
        self.cli_tester = cli_tester
        self.sim_state = sim_state

    # ..................{ RUNNERS                            }..................
    def run_subcommands_default(self) -> None:
        '''
        Run all default simulation-specific BETSE CLI subcommands in the active
        Python process.

        See Also
        ----------
        :attr:`_SUBCOMMANDS_ARGS_DEFAULT`
            Constant listing the argument lists comprising these subcommands.
        '''

        self.run_subcommands(*self._SUBCOMMANDS_ARGS_DEFAULT)


    @type_check
    def run_subcommands(
        self,
        *subcommands_args: SequenceTypes,
        is_config_overwrite: bool = True
    ) -> None:
        '''
        Run all simulation-specific BETSE CLI subcommands signified by the
        passed argument lists in the active Python process (in the passed
        order).

        **Order is significant.** subcommands producing output required by
        subsequent subcommands as input should be passed first.

        To guarantee that each such subcommand efficiently reuses the same
        underlying simulation, this method implicitly:

        * Appends each such argument list with the absolute path of the
          simulation configuration file with which this test runner was
          initialized.
        * Temporarily changes the current working directory (CWD) to the
          directory containing this file.

        Parameters
        ----------
        subcommands_args : tuple
            Tuple of sequences of **subcommand arguments** (i.e., one or more
            shell words comprising the BETSE CLI subcommand to be tested).
        is_config_overwrite : optional[bool]
            If ``True``, all in-memory changes to the current simulation
            configuration are persisted back to disk. Defaults to ``True``.
        '''

        # If persisting all in-memory configuration changes back to disk, do so
        # *BEFORE* running the passed.subcommands requiring these changes.
        if is_config_overwrite:
            self._overwrite_config()

        # For each such subcommand...
        for subcommand_args in subcommands_args:
            # Run this subcommand, preventing the configuration file from being
            # re-overwritten. While technically safe, doing so would impose
            # unnecessary I/O inefficiencies.
            self.run_subcommand(
                *subcommand_args, is_config_overwrite=False)


    @type_check
    def run_subcommand(
        self, *subcommand_args: str, is_config_overwrite: bool = True
    ) -> None:
        '''
        Run the simulation-specific BETSE CLI subcommand signified by the passed
        argument list in the active Python process.

        This method implicitly:

        * Temporarily changes the current working directory (CWD) to the
          directory containing the desired simulation configuration file.
        * Appends this argument list by the basename of this file.

        The basename rather than absolute path of this file is passed to ensure
        that a fatal error is raised if the test requesting this fixture object
        failed to change the current working directory (CWD) to the directory
        containing this file. Why? Because operating outside of this directory
        encourages accidental permanent modification of the filesystem by tests
        and hence _must_ be discouraged.

        Parameters
        ----------
        subcommand_args : tuple
            Tuple of **subcommand arguments** (i.e., one or more shell words
            comprising the BETSE CLI subcommand to be tested).
        is_config_overwrite : optional[bool]
            If ``True``, all in-memory changes to the current simulation
            configuration are persisted back to disk. Defaults to ``True``.
        '''

        # Defer heavyweight imports.
        from betse.util.io.log import logs
        from betse.util.type.text import strs

        # If persisting all in-memory configuration changes back to disk, do so
        # *BEFORE* running the passed.subcommands requiring these changes.
        if is_config_overwrite:
            self._overwrite_config()

        # Human-readable title to be embedded in a single-line terminal banner.
        subcommand_banner_title = strs.join_on(subcommand_args, delimiter=' ')

        # Log a single-line terminal banner embedding this title.
        logs.log_banner(title=subcommand_banner_title, padding='=')

        # Temporarily change the CWD to this simulation file's directory.
        with self.sim_state.context():
            # Append the absolute path of this runner's configuration file to
            # the passed tuple of arguments. While inefficient, converting this
            # tuple into a list would be even more inefficient.
            #
            # Avoid shell-quoting this path. Doing so unnecessarily adds an
            # additional level of quoting... which is bad.
            subcommand_args += (self.sim_state.conf_filename,)

            # Run this subcommand.
            self.cli_tester.run(*subcommand_args)

    # ..................{ PRIVATE                            }..................
    def _overwrite_config(self) -> None:
        '''
        Persist all in-memory changes to the current simulation configuration
        wrapper back to disk.
        '''

        self.sim_state.p.save_inplace()
