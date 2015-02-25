#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
`betse`'s command line interface (CLI).
'''

# ....................{ IMPORTS                            }....................
from argparse import ArgumentParser
from betse import metadata, pathtree
from betse.cli import help
from betse.cli.cli import CLI
from betse.util.io import loggers
from betse.util.path import dirs, files
from betse.util.type import strs
from betse.util.system import processes
from collections import OrderedDict

# ....................{ MAIN                               }....................
def main() -> int:
    '''Run `betse`'s command line interface (CLI).

    This function is provided as a convenience to callers requiring procedural
    functions rather than classical methods (e.g., `setuptools`).

    Returns
    ----------
    int
        Exit status of such interface. This is a non-negative integer in
        `[0, 255]` where 0 signifies success and all other values failure.
    '''
    return CLICLI().run()

# ....................{ CLASS                              }....................
class CLICLI(CLI):
    '''`betse`'s command line interface (CLI).

    Attributes
    ----------
    _arg_subparsers : ArgumentParser
        `argparse`-specific list of subparsers parsing top-level subcommands
        (e.g., `sim`).
    _arg_parser_sim : ArgumentParser
        `argparse`-specific subparser of command-line arguments passed to the
        `sim` subcommand.
    '''
    # ..................{ SUPERCLASS                         }..................
    def _configure_arg_parsing(self):
        # List of argument subparsers parsing arguments for root subcommands.
        self._arg_subparsers = self._arg_parser.add_subparsers(
            # Name of the attribute storing the passed command name.
            dest = 'command_name',

            # Title of the subcommand section in help output.
            title = 'subcommands',

            # Description of the subcommand section in help output.
            description = self._format_help_template(
                help.TEMPLATE_SUBCOMMANDS),
        )

        # ................{ SUBPARSER ~ sim                    }................
        #FIXME: Implement me.

        self._add_subparser(
            name = 'try',
            help = 'run a sample tissue simulation',
            description = self._format_help_template(
                help.TEMPLATE_SUBCOMMAND_TRY),
        )
        self._configure_arg_parsing_sim()

        # ................{ SUBPARSER ~ info                   }................
        self._add_subparser(
            name = 'info',
            help = 'print informational metadata about the program',
            description = self._format_help_template(
                help.TEMPLATE_SUBCOMMAND_INFO),
        )

    def _run(self) -> None:
        '''
        Run `betse`'s command line interface (CLI).
        '''
        # If no subcommand was passed, print help output and return.
        if not self._args.command_name:
            self._arg_parser.print_help()
            return

        # Else, a subcommand was passed.
        #
        # Name of the method running such subcommand.
        subcommand_method_name = '_run_' + self._args.command_name

        # Method running such subcommand. If such method does *NOT* exist,
        # getattr() will raise a non-layman-readable exception. Typically, this
        # would be bad. In this case, however, argument parsing coupled with a
        # reliable class implementation guarantees such method to exist.
        subcommand_method = getattr(self, subcommand_method_name)

        # Run such subcommand.
        subcommand_method()

    # ..................{ SUBPARSERS                         }..................
    def _add_subparser(self, *args, **kwargs) -> ArgumentParser:
        '''
        Add and return a top-level command-line argument subparser, initialized
        with the passed positional and keyword arguments.
        '''
        # Extend the passed dictionary of keyword arguments with the dictionary
        # of preinitialized keyword arguments.
        kwargs.update(self._arg_parser_kwargs)

        # Add such subparser.
        return self._arg_subparsers.add_parser(*args, **kwargs)

    def _configure_arg_parsing_sim(self):
        '''
        Configure argument parsing for the `sim` subcommand.
        '''
        # Create and return a new simulation-specific argument subparser
        # initialized with the set of passed keyword parameters.
        #
        # Such subparser will be preconfigured to parse options `-c` and
        # `--config-file`, specifying such simulation's configuration file.
        self._arg_parser_sim = self._add_subparser(
            name = 'sim',
            help = 'run tissue simulation subcommand(s)',
            description = self._format_help_template(
                help.TEMPLATE_SUBCOMMAND_SIM),
        )
        self._arg_parser_sim.add_argument(
            'sim_subcommand_names',
            nargs = '+',
            metavar = 'SUBCOMMAND_NAME',
            help = 'simulation subcommand(s) to be run',
        )
        self._arg_parser_sim.add_argument(
            'sim_config_filename',
            metavar = 'CONFIG_FILE',
            help = 'simulation configuration file',
        )

    # ..................{ SUBCOMMANDS ~ info                 }..................
    def _run_info(self) -> None:
        '''
        Run the `info` subcommand.
        '''
        #FIXME: Also print the versions of installed mandatory dependencies.
        #FIXME; For aesthetics, convert to yppy-style "cli.memory_table" output.

        # Dictionary of string keys and string values to be output below,
        # ordered so as to preserve the specified order.
        info_key_to_value = OrderedDict((
            ('script basename', processes.get_current_basename()),
            ('program version', metadata.__version__),
            ('program authors', metadata.AUTHORS),
            ('home directory', pathtree.HOME_DIRNAME),
            ('dot directory',  pathtree.DOT_DIRNAME),
            ('data directory', pathtree.DATA_DIRNAME),
            ('log file', pathtree.LOG_DEFAULT_FILENAME),
            ('default simulation config file',
             pathtree.SIMULATION_CONFIG_DEFAULT_FILENAME),
        ))

        # String to be output by this subcommand.
        info_output = '\n' + strs.join_on_newline(
            '{}: {}'.format(info_key, info_value)
            for info_key, info_value in info_key_to_value.items()
        )

        # Log rather than merely output such string, as such logging simplifies
        # cliest-side bug reporting.
        loggers.log_info(info_output)

    # ..................{ SUBCOMMANDS ~ sim                  }..................
    def _run_sim(self) -> None:
        '''
        Run the `sim` subcommand.
        '''
        # loggers.log_debug(
        #     'sim subcommand names: %s', str(self._args.sim_subcommand_names))
        # loggers.log_debug(
        #     'sim_config_filename: %s', self._args.sim_config_filename)

        # Run each simulation-specific subcommand.
        for sim_subcommand_name in self._args.sim_subcommand_names:
            # Name of the method running such subcommand.
            sim_subcommand_method_name = '_run_sim_' + sim_subcommand_name

            # Method running such subcommand, passing "None" to prevent
            # getattr() from raising a non-layman-readable exception on
            # unrecognized subcommands.
            sim_subcommand_method = getattr(
                self, sim_subcommand_method_name, None)

            # If such subcommand is unrecognized, print help output and exit the
            # current process with failure.
            if not sim_subcommand_method:
                # Error message to be printed.
                exit_message = (
                    '"sim" subcommand "{}" unrecognized. '
                    'See the following usage documentation for '
                    'the list of subcommands recognized by "sim".\n\n{}'
                ).format(
                    sim_subcommand_name,
                    self._arg_parser_sim.format_help(),
                )

                # Exit with such message.
                processes.exit_with_failure(exit_message)

            # Else, run such subcommand.
            sim_subcommand_method()

    def _run_sim_cfg(self) -> None:
        '''
        Run the `sim` subcommand's `cfg` subcommand.
        '''
        # Create the directory to which such file will be written if needed.
        dirs.make_parent_unless_found(self._args.sim_config_filename)

        # Copy the default configuration file to such file.
        files.copy(
            pathtree.SIMULATION_CONFIG_DEFAULT_FILENAME,
            self._args.sim_config_filename,
        )

# --------------------( WASTELANDS                         )--------------------
    #FUXME: Get me working correctly *BEFORE* implementing any others.

            # description = (
            #     'Run the passed tissue simulation subcommand(s) '
            #     'configured by the passed configuration file. For example, '
            #     'to initialize, run, and plot a tissue simulation '
            #     'configured by a file "my_sim.yaml" '
            #     'in the current directory:\n\n'
            #     ':    {script_basename} sim init run plot my_sim.yaml\n\n'
            #     'Valid subcommands include:\n\n'
            #     'cfg\n'
            #     ' ----------\n'
            #     'Write a default tissue simulation configuration to '
            #     'the passed output file, '
            #     'whose filename should (ideally) be suffixed by ".yaml". '
            #     'For portability, this configuration will save '
            #     'simulation results and plots '
            #     'to the directory in which this file resides. '
            #     'You are welcome to modify this file at any time.\n\n'
            #     'init\n'
            #     ' ----------\n'
            #     'Initialize the tissue simulation '
            #     'configured by the passed input configuration file. '
            #     'Initialization results will be saved to the output file '
            #     'configured in such configuration.\n\n'
            #     'run\n'
            #     ' ----------\n'
            #     'Run the previously initialized tissue simulation '
            #     'configured by the passed input configuration file. '
            #     'Simulation results will be saved to the output file '
            #     'configured in such configuration. '
            #     'Likewise, the previously run initialization '
            #     'will be loaded from the input file '
            #     'configured in such configuration. '
            #     'If such file does not exist, an error is raised.\n\n'
            #     'plot\n '
            #     ' ----------\n'
            #     'Plot the previously run tissue simulation '
            #     'configured by the passed input configuration file. '
            #     'Plot results will be saved to the output files '
            #     'configured in such configuration. '
            #     'Likewise, the previously run simulation '
            #     'will be loaded from the input file '
            #     'configured in such configuration. '
            #     'If such file does not exist, an error is raised.'
            # ).format(
            #     script_basename = self._script_basename,
            # ),

            # Pass preinitialized keyword arguments.
            # **self._arg_parser_kwargs
    # def __init__(self):
    #     super().__init__()

            # description = 'Subcommand to be performed.',
            # description = (
            #     'Run a sample tissue simulation by '
            #     'automatically creating a default configuration file and '
            #     'initializing, running, and plotting the tissue simulation '
            #     'configured by such file. '
            #     'This convenience command is equivalent to the following:\n\n'
            #     '    {} sim cfg init run plot sim_config.yaml\n\n'
            # ).format(self._script_basename),

            #FUXME: Ugh. "\n" escapes are ignored in descriptions, so we'll need
            #to manually embed newlines. Hmm. Would even that work? It might be
            #that argparse squelches *ALL* newlines in descriptions. Google up.

                # 'will instruct '
                # 'other tissue simulation subcommands (e.g., "run", "plot") '
                # 'to the current directory.\n\n'
# from argparse import ArgumentParser, _SubParsersAction
                    # 'Run "{} sim --help" for a list of supported'
    # def _configure_arg_parsing_simulation(
    #     self, arg_subparsers_main: ArgumentParser):
    #     '''
    #     Configure subclass-specific argument parsing for the `sim` subcommand
    #     given the passed list of top-level argument parsers.
    #     '''
    #     assert isinstance(arg_subparsers_main, ArgumentParser),\
    #         '"{}" not an argument parser.'.format(arg_subparsers_main)

    # def _configure_arg_parsing(self, arg_parser: ArgumentParser):
    #     assert isinstance(arg_parser, ArgumentParser),\
    #         '"{}" not an argument parser.'.format(arg_parser)
                # program_name = metadata.NAME,
        # ................{ SUBPARSER ~ sim                    }................
        # self._add_subparser_simulation(
        #     name = 'sim.init,run,plot',
        #     help = 'initialize, run, and plot a tissue simulation',
        #     description = (
        #         'Initialize, run, and plot the tissue simulation '
        #         'specified by the passed configuration file. '
        #         'This subcommand aggregates the behaviour of '
        #         'the "sim.init", "sim.run", and "sim.plot" subcommands. '
        #         'Simulation results and plots will be saved to '
        #         'the output files '
        #         'specified by such configuration. '
        #         'The simulation will be initialized before being run by '
        #         'loading the input initialization file '
        #         'specified by such configuration. '
        #         'If such file does not exist, '
        #         'it will be automatically created on your behalf by '
        #         'performing an initialization. '
        #         'caching simulation results to the passed output file. '
        #         '(In short, this subcommand always tries to '
        #         '"do the right thing.")'
        #     )
        # )
        #
        # # ................{ SUBPARSER ~ sim.init               }................
        # self._add_subparser_simulation(
        #     name = 'sim.cfg',
        #     help = 'make a default tissue simulation configuration',
        #     description = (
        #         'Write a default tissue simulation configuration to '
        #         'the passed output file. Such configuration will instruct {} '
        #         'to save simulation results and plots '
        #         'to sensibly named files in the current directory.'
        #     ).format(metadata.NAME),
        # )
        #
        # # ................{ SUBPARSER ~ sim.init               }................
        # self._add_subparser_simulation(
        #     name = 'sim.init',
        #     help = 'initialize a tissue simulation',
        #     description = (
        #         'Initialize the tissue simulation '
        #         'specified by the passed configuration file. '
        #         'Initialization results will be saved to '
        #         'the output file '
        #         'specified by such configuration. '
        #     ),
        # )
        #
        # # ................{ SUBPARSER ~ sim.run                }................
        # self._add_subparser_simulation(
        #     name = 'sim.run',
        #     help = 'run a previously initialized tissue simulation',
        #     description = (
        #         'Run the previously initialized tissue simulation '
        #         'specified by the passed configuration file. '
        #         'Simulation results will be saved to '
        #         'the output file '
        #         'specified by such configuration. '
        #         'The simulation will be initialized before being run by '
        #         'loading the input initialization file '
        #         'specified by such configuration. '
        #         'If such file does not exist, '
        #         'this subcommand will fail with an error.'
        #     )
        # )
        #
        # # ................{ SUBPARSER ~ sim.plot               }................
        # self._add_subparser_simulation(
        #     name = 'sim.plot',
        #     help = 'plot a previously run tissue simulation',
        #     description = (
        #         'Plot the previously run tissue simulation '
        #         'specified by the passed configuration file. '
        #         'Plot results will be saved to '
        #         'the output files '
        #         'specified by such configuration. '
        #         'The simulation will be loaded before being plotted by '
        #         'loading the input simulation file '
        #         'specified by such configuration. '
        #         'If such file does not exist, '
        #         'this subcommand will fail with an error.'
        #     ),
        # )
#
    #FUXME: This method should arguably simply be passed "self._arg_parser",
    #which would then permit us to localize such field in the "CLI" class.

# ..................{ SUBPARSERS                         }..................
# def add_arg_subparsers_parser(
#     arg_subparsers: ArgumentParser, **kwargs) -> _SubParsersAction:
#     '''
#     Create and add a new argument subparser to the passed list of argument
#     subparsers, initializing such subparser with the set of passed keyword
#     parameters and returning such subparser.
#
#     Additionally, if the `help` keyword parameter is passed *and* the
#     `description` keyword parameter is not, the latter will be implicitly
#     synthesized from the former.
#     '''
#     assert isinstance(arg_subparsers, ArgumentParser),\
#         '"{}" not an argument parser.'.format(arg_subparsers)
#
#     # If the "help" parameter was passed *AND* the "description" parameter
#     # was not, synthesize the latter from the former.
#     if 'help' in kwargs and 'description' not in kwargs:
#         kwargs['description'] = kwargs['help'].capitalize() + '.'
#
#     # Make such subparser.
#     arg_subparser = arg_subparsers.add_parser(**kwargs)
#
#     # Get such subparser.
#     return arg_subparser
#
    # _arg_subparsers : argparse._SubParsersAction
    #     `argparse`-specific object storing all *argument subparsers* (i.e.,
    #     parsers parsing subcommand-specific arguments).
        # Create and return a new argument subparser initialized with the set of
        # passed keyword parameters.
    # def _add_subparser_simulation(self, **kwargs) -> argparse._SubParsersAction:
    #     '''
    #     Create and return a new simulation-specific argument subparser
    #     initialized with the set of passed keyword parameters.
    #
    #     Such subparser will be preconfigured to parse options `-c` and
    #     `--config-file`, specifying such simulation's configuration file.
    #     '''
    #     subparser_sim = self._add_subparser(**kwargs)
    #     subparser_sim.add_argument(
    #         'sim_config_filename',
    #         metavar = 'CONFIG_FILE',
    #         help = 'simulation configuration file'
    #     )
    #     return subparser_sim

        # ................{ SUBPARSER ~ sim                    }................
        # self._add_subparser_simulation(
        #     name = 'sim',
        #     help = 'initialize, run, and plot a tissue simulation',
        #     description = (
        #         'Initialize, run, and plot the tissue simulation '
        #         'specified by the passed configuration file. '
        #         'This subcommand aggregates the behaviour of '
        #         'the "sim.init", "sim.run", and "sim.plot" subcommands. '
        #         'Simulation results and plots will be saved to '
        #         'the output files '
        #         'specified by such configuration. '
        #         'The simulation will be initialized before being run by '
        #         'loading the input initialization file '
        #         'specified by such configuration. '
        #         'If such file does not exist, '
        #         'it will be automatically created on your behalf by '
        #         'performing an initialization. '
        #         'caching simulation results to the passed output file. '
        #         '(In short, this subcommand always tries to '
        #         '"do the right thing.")'
        #     )
        # )
#
                # 'If such file does not exist, '
                # 'This subcommand does *NOT* plot such simulation but is '
                # 'otherwise identical to the "run" subcommand.'
                # 'See help for the "run" subcommand for further details.'
                # 'Unlike the "run" subcommand, '
                # 'this subcommand does *NOT* plot such simulation. '
    # def _add_subparser_simulation(self, **kwargs) -> argparse._SubParsersAction:
    #     '''
    #     Create and return a new simulation-specific argument subparser
    #     initialized with the set of passed keyword parameters.
    #
    #     Such subparser will be preconfigured to parse options `-c` and
    #     `--config-file`, specifying such simulation's configuration file.
    #     '''
    #     subparser_sim = self._add_subparser_simulation(**kwargs)
    #     subparser_sim.add_argument(
    #         '-c', '--config-file',
    #         dest = 'sim_config_filename',
    #         help = 'simulation configuration file'
    #     )
    #     return subparser_sim

        # ................{ SUBPARSER ~ sim                    }................
        # subparser_sim = self._add_subparser_simulation(
        #     name = 'sim',
        #     help = 'run and plot a tissue simulation',
        #     description = (
        #         'Run and plot a tissue simulation, '
        #         'cache the simulation results to the passed output file. '
        #         'If an optional input file is passed, '
        #         'the simulation will be initialized from '
        #         'the previously cached contents of such file; '
        #         'otherwise, an initialization will be '
        #         'automatically performed on your behalf and '
        #         'the simulation then initialized from '
        #         'the contents of the file produced by such initialization. '
        #         '(In short, this subcommand always tries to '
        #         '"do the right thing.")'
        #     )
        # )
        # subparser_sim.add_argument(
        #     '-c', '--config-file',
        #     dest = 'sim_config_filename',
        #     help = 'simulation configuration file'
        # )
        #
        # # ................{ SUBPARSER ~ sim.run                }................
        # subparser_sim = self._add_subparser_simulation(
        #     name = 'sim.run',
        #     help = 'run a tissue simulation',
        #     description = (
        #         'Run a tissue simulation and '
        #         'cache the results to the passed output file. '
        #         'If an optional input file is passed, '
        #         'the simulation will be initialized from '
        #         'the previously cached contents of such file; '
        #         'otherwise, an initialization will be '
        #         'automatically performed on your behalf and '
        #         'the simulation then initialized from '
        #         'the contents of the file produced by such initialization. '
        #         '(In short, this subcommand always tries to '
        #         '"do the right thing.")'
        #     )
        # )
        # subparser_sim.add_argument(
        #     '-c', '--config-file',
        #     dest = 'sim_config_filename',
        #     help = 'simulation configuration file'
        # )
        #
        # # ................{ SUBPARSER ~ sim.init               }................
        # subparser_sim_init = self._add_subparser_simulation(
        #     name = 'sim.init',
        #     help = 'initialize a tissue simulation (to be subsequently run)',
        #     description = (
        #         'Initialize a tissue simulation and '
        #         'cache the results to the passed output file.'
        #     ),
        # )

                # 'initialized on your behalf before being run. '
                # 'If this simulation has *NOT* yet been '
                # 'initialized, such simulation will be automatically '
                # 'initialized on your behalf before being run. '
                # '(This is a good thing.)'
        # Add argument subparsers accepting options.
        # Add argument subparsers accepting no options.
# from betse.util.io import stdout
        # stdout.output_lines(
        #FUXME: Contemplate localizing.

# import inspect
        # Dictionary from subcommand name to _run_*() method running such
        # subcommand.
        # subcommand_name_to_method = dict(
        #     (method_name[len('_run_'):], method)
        #     for (method_name, method) in
        #         inspect.getmembers(self, inspect.ismethod)
        #     if method_name.startswith('_run_')
        # )
        # assert self._args.command_name in subcommand_name_to_method,\
        #     '"{}" not a recognized subcommand'.format(self._args.command_name)
        #
        # # Run such subcommand.
        # subcommand_name_to_method[self._args.command_name]()
        # subcommand_method = getattr(self, subcommand_method_name, None)
        # assert callable(subcommand_method),\
        #     '"{}" not callable'.format(subcommand_method_name)
        #FUXME: Implement me.
        # Run the command specified by such arguments.

        #FUXME: Display a default help message, when the user passes no
        #arguments. Does such parser already do so? This is trivial for us
        #to do as follows:
        #
        #    if len(sys.argv) == 1:
        #         self._arg_parser.print_help()
        #         return

        #FUXME: Not entirely clear as to why we require or want this dictionary.
        #Smacks of overkill, frankly. Also unclear why we require
        #_add_subparser_sim()-style methods. Just define all such subcommands
        #here, for now.

        # # Dictionary from command name to subparser object.
        # command_parsers = {}
        # command_parsers['sim'] = self._add_subparser_sim()
        #
        # # Add an identifying name and description for each command parser.
        # for command_parser_name, command_parser in command_parsers.items():
        #     command_parser.set_defaults(command_name = command_parser_name)

        # Parse command-line arguments into object attributes.
        # self._parse_args()
        # self._parse_common_args()

#FUXME; The following snippet courtesy Matthew Leingan affords an elegant means
#of integrating built-in Python argument parsing and logging;
#
#    import argparse
#    import logging
#
#    parser = argparse.ArgumentParser()
#    parser.add_argument('-d','--debug',
#        help='Print lots of debugging statements',
#        action="store_const",dest="loglevel",const=logging.DEBUG,
#        default=logging.WARNING
#    )
#    parser.add_argument('-v','--verbose',
#        help='Be verbose',
#        action="store_const",dest="loglevel",const=logging.INFO
#    )
#    args = parser.parse_args()
#    logging.basicConfig(level=args.loglevel)
#
#In his own words: "So if --debug is set, the logging level is set to DEBUG. If
#--verbose, logging is set to INFO. If neither, the lack of --debug sets the
#logging level to the default of WARNING."
#
#That said, devoting a separate command-line option to "--debug" strikes us as
#overkill when we can simply repeat "-vv" to achive the same, ala this snippet;
#
#    # This isn't quite right, as it uses "--verbose=1"-style integer
#    # assignment. We just want to repeat "-v" and/or "--verbose". Much simpler.
#    # Nonetheless, this should be of some use.
#    parser = argparse.ArgumentParser()
#    parser.add_argument("-v", "--verbose", const=1, default=0, type=int, nargs="?",
#                        help="increase verbosity: 0 = only warnings, 1 = info, 2 = debug. No number means info. Default is no verbosity.")
#    args = parser.parse_args()
#
#    logger = logging.getLogger()
#    if args.verbose == 0:
#        logger.setLevel(logging.WARN)
#    elif args.verbose == 1:
#        logger.setLevel(logging.INFO)
#    elif args.verbose == 2:
#        logger.setLevel(logging.DEBUG)

# from betse.util.io import loggers
            # assert action_parser, '"{}" is None'.format(action_parser_name)
        #FUXME: This doesn't seem quite right. Don't we want to parse arguments
        #only after defining all argument syntax for actions as well?

    # @property
    # def _script_basename(self) -> str:
    #     return metadata.SCRIPT_NAME_CLI

        # return exception.errno if hasattr(exception, 'errno') else 1
#FUXME: Define a new module "betse/dependency.py" performing validation of
#external dependencies, both Python and non-Python. Although we believe "yppy"
#implemented such functionality, google about for the optimum Python 3 solution
#to this presumably commonplace problem.

    # _args : list
    #     List of zero or more arguments passed to such interface (e.g., from the
    #     command line).
        # Initialize such arguments to the current argument list, excluding
        # such list's first item. By cross-platform consent, such item is
        # *ALWAYS* the command name for the current process (e.g., "betse") and
        # hence ignorable.
        # self._args = sys.argv[1:]

#List of zero or more external arguments passed from the command line.
# def main(args = None):
#     CLI().run()
#     '''Run betse`'s command line interface (CLI).
#
#     Parameters
#     ----------
#     args : list, optional
#         List of zero or more arguments passed to such interface (e.g., from the
#         command line) or `None` if called as the entry point in an external
#         script installed by `setuptools`.
#     '''
#     # If called from a setuptools-installed script, copy such arguments from the
#     # argument list excluding the first item of such list. By cross-platform
#     # agreement, such item is *ALWAYS* the command name of the current process
#     # (e.g., "betse") and hence ignorable.
#     if args is None:
#         args = sys.argv[1:]

# if __name__ == '__main__':
#     main()

#FUXME; Configure me for CLI usage. Note that I'm no longer convinced that the
#way we launched "yppy" (e.g., "bin/yppy.bash") was ideal. We really want to do
#the "Pythonic" thing here. ruamel.yaml, for example, installs a Python wrapper
#"/usr/lib/yaml" which (in order):
#
#* Finds an appropriate Python interpreter.
#* Replaces the current process with the result of interpreting
#  "/usr/lib/python-exec/python${PYTHON_VERSION}/yaml". Such file appears to be
#  autogenerated by setuptools at installation time.
#FUXME; Hmm: it looks like we want a new file "betse/__main__.py" resembling:
#    from betse.main import main
#    main()
#This then permits betse to be run as follows:
#    # Yes, you either have to be in the parent directory of the directory
#    # containing such "__main__.py" file *OR* you have to fiddle with
#    # ${PYTHONPATH}.
#    >>> cd ~/py/betse
#    >>> python -m betse
#Naturally, this lends itself well to shell scripting. (Yay!)
#FUXME; Wo! Even nicer. setuptools has implicit support for "__main__.py"-style
#entry points. We just need a "setup.py" resembling:
#    setup(
#        # [...]
#        entry_points={
#            'betse': ['betse = betse.main:main'],
#        },
#    )
#What's sweet about this is that we can define additional separate scripts with
#deeper entry points if we need and or want to.
