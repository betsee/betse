#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Top-level abstract base class of all **subcommandable command line interface
(CLI)** (i.e., CLI accepting one or more subcommands) subclasses.
'''

# ....................{ IMPORTS                            }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: To raise human-readable exceptions on application startup, the
# top-level of this module may import *ONLY* from submodules guaranteed to:
# * Exist, including standard Python and application modules.
# * Never raise exceptions on importation (e.g., due to module-level logic).
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

from betse.cli import cliinfo, cliutil
from betse.cli.api.cliabc import CLIABC
from betse.util.io.log import logs
from betse.util.path import files, pathnames
from betse.util.py import pyident, pys
from betse.util.type.call.memoizers import property_cached
from betse.util.type.obj import objects
from betse.util.type.types import SequenceTypes

# ....................{ SUBCLASS                           }....................
class CLISubcommandableABC(CLIABC):
    '''
    Top-level abstract base class of all **subcommandable command line interface
    (CLI)** (i.e., CLI accepting one or more subcommands) subclasses, suitable
    for use by both CLI and GUI front-ends for BETSE.

    Unlike the parent :class:`CLIABC` superclass, this superclass provides
    explicit support for subcommands. Concrete subclasses implementing
    subcommands should directly subclass this rather than that superclass.

    Attributes
    ----------
    _arg_subparsers_top : ArgParserType
        Container of all argument subparsers parsing all top-level subcommands
        accepted by this application (e.g., ``init``, ``sim``).
    '''

    # ..................{ INITIALIZERS                       }..................
    def __init__(self) -> None:

        # Initialize our superclass.
        super().__init__()

        # Nullify attributes for safety.
        self._arg_subparsers_top = None

    # ..................{ SUBCLASS ~ property                }..................
    @property
    def _help_subcommands_prefix(self) -> str:
        '''
        Help string template preceding the list of all subcommands.
        '''

        return '''
Exactly one of the following subcommands must be passed:
'''

    # ..................{ SUPERCLASS ~ args                  }..................
    def _config_arg_parsing(self) -> None:

        # Container of all top-level argument subparsers for this application.
        self._arg_subparsers_top = self._arg_parser.add_subparsers(
            # Name of the attribute storing the passed subcommand name.
            dest='subcommand_name_top',

            # Title of the subcommand section in help output.
            title='subcommands',

            # Description to be printed *BEFORE* subcommand help.
            description=self.expand_help(self._help_subcommands_prefix),
        )

        #FIXME: Consider refactoring this logic as follows:
        #
        #* The local "subcommand_name_to_subparser" dictionary should be
        #  replaced by performing the following in the body of the
        #  cliutil.add_arg_subparsers_subcommands() function:
        #  * If a global "help.SUBCOMMANDS_{subcommand.name}" dictionary
        #    exists, non-recursively loop over that dictionary in the same
        #    manner as well. Hence, this generalizes support of subcommand
        #    subcommands. Nice!
        #* Likewise, if a "self._arg_parser_{subcommand.name}" instance
        #  variable exists, set that variable to this parser.
        #FIXME: The above is an intrepid starting point. Ultimately, however,
        #what would be required to genuinely refactor this for the better would
        #be to abstract the tuples returned by the _make_subcommands_top() and
        #_make_subcommands_plot() methods into instances of a new
        #"betse.cli.clicmd.CLISubcommands" class whose __init__() method for
        #this class should accept a sequence of all subcommands as well as *ALL*
        #parameters required to fully create argument subparsers from the
        #subcommands listed by this class. This includes:
        #
        #     @type_check
        #     def __init__(
        #         self,
        #
        #         # Sequence of all subcommands provided by this command.
        #         subcommands: SequenceTypes,
        #
        #         # Name of the attribute storing the passed subcommand name.
        #         dest: str = 'subcommand_name_top',
        #
        #         # Title of the subcommand section in help output.
        #         title: str = 'subcommands',
        #
        #         # Description to be printed *BEFORE* subcommand help.
        #         description: str = cliutil.expand_help(SUBCOMMANDS_PREFIX),
        #
        #         # Tuple of all "CLISubcommandABC" instances.
        #         subcommands: SequenceTypes,
        #     )
        #
        #Given such a class, it should then be feasible for instances of this
        #class to nest instances of this class. How? By generalizing the
        #existing "CLISubcommandParentObsolete" class (which is used to create
        #subcommands themselves containing subcommands, like the "plot"
        #subcommand) to provide a "subcommands" instance variable of type
        #"CLISubcommands", presumably passed to the
        #CLISubcommandParentObsolete.__init__() method.
        #
        #While performing this generalization currently exceeds our capacity for
        #sanity, it would be fairly sweet (if a little overkill, considering we
        #only nest one level deep at the moment). *shrug*

        # Dictionary mapping from the name of each top-level subcommand to the
        # argument subparser parsing that subcommand.
        subcommand_name_to_subparser = cliutil.add_arg_subparsers_subcommands(
            subcommands=self._make_subcommands_top(),
            arg_subparsers=self._arg_subparsers_top,
            arg_subparser_kwargs=self._arg_parser_kwargs,
        )

        # Subparser parsing arguments passed to the "plot" subcommand.
        self._arg_parser_plot = subcommand_name_to_subparser['plot']

        # Container of all sub-level argument subparsers for the "plot"
        # subcommand.
        self._arg_subparsers_plot = self._arg_parser_plot.add_subparsers(
            # Name of the attribute storing the passed subcommand name.
            dest='subcommand_name_plot',

            # Title of the subcommand section in help output.
            title='plot subcommands',

            # Description to be printed *BEFORE* subcommand help.
            description=self.expand_help(self._help_subcommands_prefix),
        )

        # Dictionary mapping from the name of each "plot" subcommand to the
        # argument subparser parsing that subcommand. Note that this dictionary
        # is non-essential and hence garbage-collected immediately.
        cliutil.add_arg_subparsers_subcommands(
            subcommands=self._make_subcommands_plot(),
            arg_subparsers=self._arg_subparsers_plot,
            arg_subparser_kwargs=self._arg_parser_kwargs,
        )

    # ....................{ SUBCLASS ~ mandatory               }....................
    # The following methods *MUST* be implemented by subclasses.

    @abstractmethod
    def _make_subcommands_top(self) -> SequenceTypes:
        '''
        Sequence of all :class:`CLISubcommandABC` instances defining the
        top-level subcommands accepted by this application.

        For each such subcommand, a corresponding argument subparser is
        subsequently created and added to the container of all top-level
        argument subparsers (i.e., :attr:`_arg_subparsers_top`).

        **Order is significant,** defining the order that the ``--help`` option
        synopsizes these subcommands in. Subcommands omitted here are *not*
        parsed by argument subparsers and are thus effectively ignored.

        Returns
        ----------
        SequenceTypes
            Sequence of all such :class:`CLISubcommandABC` instances.
        '''

        pass


    def _make_subcommands_plot(self) -> SequenceTypes:
        '''
        Sequence of all :class:`CLISubcommandABC` instances defining the
        sub-level subcommands accepted by the top-level ``plot`` subcommand
        accepted by this application.

        See Also
        ----------
        :meth:`_make_subcommands_plot`
            Further details
        '''

        return (
            CLISubcommandYAMLOnly(
                name='seed',
                synopsis='plot a seeded cell cluster for a config file',
                description='''
Plot the previously seeded cell cluster defined by the passed configuration
file. Plot results will be saved to output files defined by this configuration,
while the previously seeded cell cluster will be loaded from input files
defined by this configuration.
''',
            ),


            CLISubcommandYAMLOnly(
                name='init',
                synopsis='plot an initialized cell cluster for a config file',
                description='''
Plot the previously initialized cell cluster defined by the passed configuration
file. Plot results will be saved to output files defined by this configuration,
while the previously initialized cell cluster will be loaded from input files
defined by this configuration.
''',
            ),


            CLISubcommandYAMLOnly(
                name='sim',
                synopsis='plot a simulated cell cluster for a config file',
                description='''
Plot the previously simulated cell cluster defined by the passed configuration
file. Plot results will be saved to output files defined by this configuration,
while the previously simulated cell cluster will be loaded from input files
defined by this configuration.
''',
            ),


            CLISubcommandYAMLOnly(
                name='sim-grn',
                synopsis=(
                    'plot a simulated gene regulatory network '
                    'for a config file'
                ),
                description='''
Plot the previously simulated gene regulatory network (GRN) defined by the
passed configuration file. Plot results will be saved to output files defined by
this configuration, while the previously simulated cell cluster will be loaded
from input files defined by this configuration.
''',
            ),
        )

    # ..................{ SUPERCLASS ~ cli                   }..................
    def _do(self) -> object:
        '''
        Implement this command-line interface (CLI).

        If a subcommand was passed, this method runs this subcommand and returns
        the result of doing so; else, this method prints help output and returns
        the current instance of this object.
        '''

        # If no subcommand was passed...
        if not self._args.subcommand_name_top:
            #,Print help output. Note that this common case constitutes neither
            # a fatal error nor a non-fatal warning.
            print()
            self._arg_parser.print_help()

            # Return the current instance of this object. While trivial, this
            # behaviour simplifies memory profiling of this object.
            return self

        # Else, a subcommand was passed.
        #
        # Sanitized name of this subcommand.
        subcommand_name_top = pyident.sanitize_snakecase(
            self._args.subcommand_name_top)

        # Name of the method running this subcommand.
        subcommand_method_name = '_do_' + subcommand_name_top

        # Method running this subcommand. If this method does *NOT* exist,
        # get_method() will raise a non-human-readable exception. Usually, that
        # would be bad. In this case, however, argument parsing coupled with a
        # reliable class implementation guarantees this method to exist.
        subcommand_method = objects.get_method(
            obj=self, method_name=subcommand_method_name)

        # Run this subcommand and return the result of doing so (if any).
        return subcommand_method()

    # ..................{ SUBCOMMANDS ~ info                 }..................
    def _show_header(self) -> None:

        cliinfo.log_header()


    def _do_info(self) -> None:
        '''
        Run the ``info`` subcommand.
        '''

        cliinfo.log_info()

    # ..................{ SUBCOMMANDS ~ sim                  }..................
    def _do_try(self) -> object:
        '''
        Run the ``try`` subcommand and return the result of doing so.
        '''

        # Basename of the sample configuration file to be created.
        config_basename = 'sample_sim.yaml'

        # Relative path of this file, relative to the current directory.
        self._args.conf_filename = pathnames.join(
            'sample_sim', config_basename)

        #FIXME: Insufficient. We only want to reuse this file if this file's
        #version is identical to that of the default YAML configuration file's
        #version. Hence, this logic should (arguably) be shifted elsewhere --
        #probably into "betse.science.sim_config".

        # If this file already exists, reuse this file.
        if files.is_file(self._args.conf_filename):
            logs.log_info(
                'Reusing simulation configuration "%s".', config_basename)
        # Else, create this file.
        else:
            self._do_config()

        # Run all general-purposes phases, thus excluding network-isolated
        # phases (e.g., "_do_sim_grn"), in the expected order.
        self._do_seed()
        self._do_init()
        self._do_sim()
        self._do_plot_seed()
        self._do_plot_init()

        # Return the value returned by the last such phase, permitting this
        # subcommand to be memory profiled. While any value would technically
        # suffice, the value returned by the last such phase corresponds to a
        # complete simulation run and hence is likely to consume maximal memory.
        return self._do_plot_sim()


    def _do_config(self) -> None:
        '''
        Run the ``config`` subcommand.
        '''

        # Avoid importing modules importing dependencies at the top level.
        from betse.science.config import confio
        confio.write_default(self._args.conf_filename)


    def _do_seed(self) -> object:
        '''
        Run the ``seed`` subcommand and return the result of doing so.
        '''

        return self._sim_runner.seed()


    def _do_init(self) -> object:
        '''
        Run the ``init`` subcommand and return the result of doing so.
        '''

        return self._sim_runner.init()


    def _do_sim(self) -> object:
        '''
        Run the ``sim`` subcommand and return the result of doing so.
        '''

        return self._sim_runner.sim()


    def _do_sim_grn(self) -> object:
        '''
        Run the ``sim-grn`` subcommand and return the result of doing so.
        '''

        return self._sim_runner.sim_grn()


    def _do_plot(self) -> object:
        '''
        Run the ``plot`` subcommand and return the result of doing so.
        '''

        # If no subcommand was passed, print help output and return. See the
        # _do() method for similar logic and commentary.
        if not self._args.subcommand_name_plot:
            print()
            self._arg_parser_plot.print_help()
            return

        # Run this subcommand's subcommand and return the result of doing so.
        # See the _run() method for similar logic and commentary.
        subcommand_name_plot = pyident.sanitize_snakecase(
            self._args.subcommand_name_plot)
        subcommand_method_name = '_do_plot_' + subcommand_name_plot
        subcommand_method = getattr(self, subcommand_method_name)
        return subcommand_method()


    def _do_plot_seed(self) -> object:
        '''
        Run the ``plot`` subcommand's ``seed`` subcommand and return the result
        of doing so.
        '''

        return self._sim_runner.plot_seed()


    def _do_plot_init(self) -> object:
        '''
        Run the ``plot`` subcommand's ``init`` subcommand and return the result
        of doing so.
        '''

        return self._sim_runner.plot_init()


    def _do_plot_sim(self) -> object:
        '''
        Run the ``plot`` subcommand's ``sim`` subcommand and return the result
        of doing so.
        '''

        return self._sim_runner.plot_sim()


    def _do_plot_sim_grn(self) -> object:
        '''
        Run the ``plot`` subcommand's ``sim-grn`` subcommand and return the
        result of doing so.
        '''

        return self._sim_runner.plot_grn()


    def _do_repl(self) -> None:
        '''
        Run the ``repl`` subcommand.
        '''

        # In the unlikely edge-case of the "repl" subcommand being erroneously
        # run by a functional test, prohibit this by raising an exception.
        # Permitting this would probably cause tests to indefinitely hang.
        if pys.is_testing():
            from betse.exceptions import BetseTestException
            raise BetseTestException(
                'REPL unavailable during testing for safety.')

        # Defer heavyweight imports until *AFTER* possibly failing above.
        from betse.cli.repl import repls

        # Start the desired REPL.
        repls.start_repl()

    # ..................{ GETTERS                            }..................
    @property_cached
    def _sim_runner(self):
        '''
        Simulation runner preconfigured with sane defaults.
        '''

        # Defer heavyweight imports.
        from betse.science.simrunner import SimRunner

        # Return this runner.
        return SimRunner(conf_filename=self._args.conf_filename)
