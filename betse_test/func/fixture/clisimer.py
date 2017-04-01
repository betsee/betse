#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Fixtures and fixture classes efficiently exercising multiple subcommands of the
BETSE CLI in the active Python interpreter.
'''

# ....................{ IMPORTS                            }....................
from betse.util.type.types import type_check, SequenceTypes
from betse_test.fixture.simconfig.simconfer import SimTestState
from betse_test.func.fixture.clier import CLITester
from pytest import fixture

# ....................{ CLASSES ~ sim                      }....................
class CLISimTester(object):
    '''
    BETSE CLI simulation test runner, exercising multiple subcommands of the
    BETSE CLI (i.e., `betse`) in the active Python interpreter with the same
    temporary simulation configuration.

    Complex functional fixtures (e.g., `betse_cli_sim_old`) typically return
    instances of this class to other fixtures and tests exercising multiple
    facets of the BETSE CLI.

    Attributes
    ----------
    cli_tester : CLITester
        BETSE CLI test runner, testing a single subcommand of the official
        BETSE CLI (i.e., `betse`) in the active Python interpreter.
    sim_state: SimTestState
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

    This tuple is a convenience simplifying common-case tests exercising _only_
    these subcommands.
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(self, cli_tester: CLITester, sim_state: SimTestState) -> None:
        '''
        Initialize this test runner.

        Parameters
        ----------
        cli_tester: CLITester
            BETSE CLI test runner, testing a single subcommand of the official
            BETSE CLI (i.e., `betse`) in the active Python interpreter.
        sim_state: SimTestState
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
        is_overwriting_config: bool = True
    ) -> None:
        '''
        Run all simulation-specific BETSE CLI subcommands signified by the
        passed argument lists in the active Python process (_in the passed
        order_).

        **Order is significant.** subcommands producing output required by
        subsequent subcommands as input should be passed first.

        To guarantee that each such subcommand efficiently reuses the same
        underlying simulation, this method implicitly:

        * Appends each such argument list by the absolute path of the simulation
          configuration file with which this test runner was initialized.
        * Temporarily changes the current working directory (CWD) to the
          directory containing this file.

        Parameters
        ----------
        subcommands_args : tuple
            Tuple of sequences of **subcommand arguments** (i.e., one or more
            shell words comprising the BETSE CLI subcommand to be tested).
        is_overwriting_config : optional[bool]
            If `True`, this method persists all in-memory changes to the current
            simulation configuration wrapper back to disk. Defaults to `True`.
        '''

        # If persisting all in-memory configuration changes back to disk, do so
        # *BEFORE* running the passed.subcommands requiring these changes.
        if is_overwriting_config:
            self._overwrite_config()

        # For each such subcommand...
        for subcommand_args in subcommands_args:
            # Run this subcommand, preventing the configuration file from being
            # re-overwritten. While technically safe, doing so would impose
            # unnecessary I/O inefficiencies.
            self.run_subcommand(
                *subcommand_args, is_overwriting_config=False)


    @type_check
    def run_subcommand(
        self, *subcommand_args: str, is_overwriting_config: bool = True
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
        is_overwriting_config : optional[bool]
            If `True`, this method persists all in-memory changes to the current
            simulation configuration wrapper back to disk. Defaults to `True`.
        '''

        # If persisting all in-memory configuration changes back to disk, do so
        # *BEFORE* running the passed.subcommands requiring these changes.
        if is_overwriting_config:
            self._overwrite_config()

        # Temporarily change the CWD to this simulation file's directory.
        with self.sim_state.context():
            # Print these arguments as a single-line banner.
            self._print_subcommand_args(subcommand_args)

            # Append the absolute path of this runner's configuration file to
            # the passed tuple of arguments. While inefficient, converting this
            # tuple into a list would be even more inefficient.
            #
            # Avoid shell-quoting this path. Doing so unnecessarily adds an
            # additional level of quoting... which is bad.
            subcommand_args += (self.sim_state.config.filename,)

            # Run this subcommand.
            self.cli_tester.run(*subcommand_args)

    # ..................{ PRIVATE                            }..................
    def _overwrite_config(self) -> None:
        '''
        Persist all in-memory changes to the current simulation configuration
        wrapper back to disk.
        '''

        self.sim_state.config.overwrite()


    #FIXME: Consider integrating this method the main codebase, ideally by
    #calling this method before running each simulation phase. This banner is
    #sufficiently aesthetic that I'd certainly appreciate seeing it everywhere.
    @type_check
    def _print_subcommand_args(self, subcommand_args: tuple) -> None:
        '''
        Print the passed tuple of subcommand arguments as a single-line
        human-readable banner.
        '''

        # Defer heavyweight imports.
        from betse.util.type import strs

        # Maximum length of the single-line banner to be printed.
        BANNER_LEN_MAX = 80

        # Text declaring the passed argument list.
        BANNER_TEXT = strs.join_on(subcommand_args, delimiter=' ')
        # print('\n    subcommand_args: {}'.format(subcommand_args))
        # BANNER_TEXT = ' '.join(subcommand_args)

        # Character both preceding and following this text.
        BANNER_CHAR = '='

        # Minimal-length banner containing this text and the smallest
        # possible prefixing and suffixing banner characters.
        BANNER_MIN = '{char} {text} {char}'.format(
            text=BANNER_TEXT, char=BANNER_CHAR,)

        # If the length of this minimal-length banner exceeds the maximum
        # length, print this text with no such banner.
        if len(BANNER_MIN) > BANNER_LEN_MAX:
            print(BANNER_TEXT)
        # Else, print this text within a banner.
        else:
            # Length of the string preceding this text. To produce a uniform
            # rather than ragged left margin for aesthetic uniformity, this
            # length is a constant.
            BANNER_PREFIX_LEN = 30

            # String preceding this text.
            BANNER_PREFIX = BANNER_CHAR * BANNER_PREFIX_LEN

            # Banner to be printed excluding all suffixing banner characters.
            BANNER_SANS_SUFFIX = '{} {} '.format(BANNER_PREFIX, BANNER_TEXT)

            # Length of the string following this text. To dynamically fill all
            # remaining line space, this length contextually depends on the
            # lengths of all other strings. By the above conditional, this
            # length is guaranteed to be non-zero and need *NOT* be tested.
            BANNER_SUFFIX_LEN = BANNER_LEN_MAX - len(BANNER_SANS_SUFFIX)

            # Nonetheless, test this length for sanity.
            assert BANNER_SUFFIX_LEN > 0, (
                'Banner suffix length {} not positive.'.format(
                    BANNER_SUFFIX_LEN))

            # String following this text.
            BANNER_SUFFIX = BANNER_CHAR * BANNER_SUFFIX_LEN

            # Print this text as a single-line banner.
            print('\n{}{}'.format(BANNER_SANS_SUFFIX, BANNER_SUFFIX,))

# ....................{ FIXTURES                           }....................
# Test-scope fixture creating and returning a new object for each discrete test.
@fixture
def betse_cli_sim(
    request: '_pytest.python.FixtureRequest',
    betse_cli: CLITester,
    betse_sim_config: SimTestState,
) -> CLISimTester:
    '''
    Fixture returning a test-specific object suitable for running one or more
    simulation-specific BETSE CLI subcommands (e.g., `seed`, `sim`) with the
    temporary simulation configuration required by the current fixture or test.

    Parameters
    ----------
    request : _pytest.python.FixtureRequest
        Builtin fixture describing this fixture's parent fixture or test.
    betse_cli : CLITester
        Object running a single simulation-specific BETSE CLI subcommand.
    betse_sim_config : SimTestState
        Object encapsulating a temporary simulation configuration file.

    Returns
    ----------
    CLISimTester
        Object running multiple simulation-specific BETSE CLI subcommands.
    '''

    return CLISimTester(
        cli_tester=betse_cli,
        sim_state=betse_sim_config,
    )
