#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
`betse`'s command line interface (CLI).
'''

#FIXME: The "~/.betse/cache" directory grows fairly large fairly quickly. It'd
#be great to emit non-fatal warnings if its size exceeds some reasonable
#threshold (e.g., 1MB).

# ....................{ IMPORTS                            }....................
from argparse import ArgumentParser

from betse import metadata
from betse.cli import help, info
from betse.cli.cli import CLI
from betse.util.io.log import logs
from betse.util.path import files, paths


# ....................{ CLASS                              }....................
class CLICLI(CLI):
    '''
    `betse`'s command line interface (CLI).

    Attributes
    ----------
    _arg_subparsers_top : ArgumentParser
        Subparsers parsing top-level subcommands (e.g., `plot`).
    _arg_subparsers_plot : ArgumentParser
        Subparsers parsing `plot` subcommands (e.g., `plot seed`).
    _arg_parser_plot : ArgumentParser
        Subparser parsing arguments passed to the `plot` subcommand.
    '''
    def __init__(self):
        super().__init__()

        # Nullify attributes for safety.
        self._arg_subparsers_top = None
        self._arg_subparsers_plot = None
        self._arg_parser_plot = None

    # ..................{ SUPERCLASS ~ args                  }..................
    def _get_arg_parser_top_kwargs(self):
        # Keyword arguments passed to the top-level argument parser.
        return {
            'epilog': help.TEMPLATE_SUBCOMMANDS_SUFFIX,
        }


    def _config_arg_parsing(self):
        # Top-level subcommands.
        self._arg_subparsers_top = self._arg_parser.add_subparsers(
            # Name of the attribute storing the passed subcommand name.
            dest='subcommand_name_top',

            # Title of the subcommand section in help output.
            title='subcommands',

            # Description to be printed *BEFORE* subcommand help.
            description=self._format_help_template(
                help.TEMPLATE_SUBCOMMANDS_PREFIX),
        )

        self._add_arg_subparser_top_configured(
            name='config',
            help='create a new {} simulation configuration'.format(
                metadata.NAME),
            description=self._format_help_template(
                help.TEMPLATE_SUBCOMMAND_CONFIG),
        )
        self._add_arg_subparser_top_configured(
            name='seed',
            help='create the cell cluster defined by a config file',
            description=self._format_help_template(
                help.TEMPLATE_SUBCOMMAND_SEED),
        )
        self._add_arg_subparser_top_configured(
            name='init',
            help='init the created cell cluster defined by a config file',
            description=self._format_help_template(
                help.TEMPLATE_SUBCOMMAND_INIT),
        )
        self._add_arg_subparser_top_configured(
            name='sim',
            help='simulate the initted cell cluster defined by a config file',
            description=self._format_help_template(
                help.TEMPLATE_SUBCOMMAND_SIM),
        )
        self._config_arg_parsing_plot()
        self._add_arg_subparser_top(
            name='info',
            help='show information about {} and the current system'.format(
                metadata.NAME),
            description=self._format_help_template(
                help.TEMPLATE_SUBCOMMAND_INFO),
        )
        self._add_arg_subparser_top(
            name='try',
            help='create, init, simulate, and plot a sample simulation',
            description=self._format_help_template(
                help.TEMPLATE_SUBCOMMAND_TRY),
        )

    # ..................{ SUPERCLASS ~ cli                   }..................
    def _do(self) -> None:
        '''
        Command-line interface (CLI) for `betse`.
        '''

        # If no subcommand was passed, print help output and return. Note that
        # this does *NOT* constitute a fatal error.
        if not self._args.subcommand_name_top:
            print()
            self._arg_parser.print_help()
            return

        # Else, a subcommand was passed.
        #
        # Name of the method running this subcommand.
        subcommand_method_name = '_do_' + self._args.subcommand_name_top

        # Method running this subcommand. If this method does *NOT* exist,
        # getattr() will raise a non-human-readable exception. Usually, that
        # would be bad. In this case, however, argument parsing coupled with a
        # reliable class implementation guarantees this method to exist.
        subcommand_method = getattr(self, subcommand_method_name)

        # Run this subcommand.
        subcommand_method()

    # ..................{ SUBCOMMAND ~ plot                  }..................
    def _config_arg_parsing_plot(self):
        '''
        Configure argument parsing for the `plot` subcommand.
        '''

        # This subcommand.
        self._arg_parser_plot = self._add_arg_subparser_top(
            name='plot',
            help='plot previously created, initted, or simulated simulations',
            description=self._format_help_template(
                help.TEMPLATE_SUBCOMMAND_PLOT),
        )

        # Collection of all subcommands of this subcommand.
        self._arg_subparsers_plot = self._arg_parser_plot.add_subparsers(
            # Name of the attribute storing the passed subcommand name.
            dest='subcommand_name_plot',

            # Title of the subcommand section in help output.
            title='plot subcommands',

            # Description to be printed *BEFORE* subcommand help.
            description=self._format_help_template(
                help.TEMPLATE_SUBCOMMANDS_PREFIX),
        )
        self._add_arg_subparser_plot_configured(
            name='seed',
            help='plot the created cell cluster defined by a config file',
            description=self._format_help_template(
                help.TEMPLATE_SUBCOMMAND_PLOT_SEED),
        )
        self._add_arg_subparser_plot_configured(
            name='init',
            help='plot the initted cell cluster defined by a config file',
            description=self._format_help_template(
                help.TEMPLATE_SUBCOMMAND_PLOT_INIT),
        )
        self._add_arg_subparser_plot_configured(
            name='sim',
            help='plot the simulated cell cluster defined by a config file',
            description=self._format_help_template(
                help.TEMPLATE_SUBCOMMAND_PLOT_SIM),
        )

    # ..................{ SUBPARSER ~ top                    }..................
    def _add_arg_subparser_top_configured(
        self, *args, **kwargs) -> ArgumentParser:
        '''
        Create a new argument subparser requiring a configuration filename, add
        this subparser to the collection of top-level argument subparsers, and
        return this subparser.
        '''
        return self._add_arg_subparser_configured(
            self._arg_subparsers_top, *args, **kwargs)


    def _add_arg_subparser_top(self, *args, **kwargs) -> ArgumentParser:
        '''
        Create a new argument subparser, add such subparser to the collection
        of top-level argument subparsers, and return this subparser.
        '''
        return self._add_arg_subparser(
            self._arg_subparsers_top, *args, **kwargs)

    # ..................{ SUBPARSER ~ plot                   }..................
    def _add_arg_subparser_plot_configured(self, *args, **kwargs) -> ArgumentParser:
        '''
        Create a new argument subparser requiring a configuration filename, add
        such subparser to the subparser corresponding to the `plot` subcommand,
        and return such subparser.
        '''
        return self._add_arg_subparser_configured(
            self._arg_subparsers_plot, *args, **kwargs)

    # ..................{ SUBPARSER                          }..................
    def _add_arg_subparser_configured(
        self, arg_subparsers: ArgumentParser, *args, **kwargs) -> ArgumentParser:
        '''
        Create a new argument subparser requiring a configuration filename, add
        such subparser to the passed set of argument subparsers, and return such
        subparser.

        Parameters
        ----------
        arg_subparsers : ArgumentParser
            Set of argument subparsers, typically corresponding to a top-level
            subcommand (e.g., `plot`).
        '''
        # Create such subparser.
        arg_subparser = self._add_arg_subparser(
            arg_subparsers, *args, **kwargs)

        # Configure such subparser to require a passed configuration file.
        arg_subparser.add_argument(
            'config_filename',
            metavar = 'CONFIG_FILE',
            help = 'simulation configuration file',
        )

        # Get such subparser.
        return arg_subparser


    def _add_arg_subparser(
        self, arg_subparsers, *args, **kwargs) -> ArgumentParser:
        '''
        Create a new **argument subparser** (i.e., an `argparse`-specific object
        responsible for parsing a single command-line argument), add such
        subparser to the passed **argument subparsers** (i.e., another
        `argparse`-specific object cotaining multiple subparsers), and return
        such subparser.

        All remaining positional and keyword arguments are passed as is to such
        subparser's `__init__()` method.
        '''
        # Extend the passed dictionary of keyword arguments with the dictionary
        # of preinitialized keyword arguments.
        kwargs.update(self._arg_parser_kwargs)

        # Add such subparser.
        return arg_subparsers.add_parser(*args, **kwargs)

    # ..................{ SUBCOMMANDS ~ info                 }..................
    def _do_info(self) -> None:
        '''
        Run the `info` subcommand.
        '''

        info.output_info()

    # ..................{ SUBCOMMANDS ~ sim                  }..................
    def _do_try(self) -> None:
        '''
        Run the `try` subcommand.
        '''

        # Basename of the sample configuration file to be created.
        config_basename = 'sample_sim.yaml'

        # Relative path of this file, relative to the current directory.
        self._args.config_filename = paths.join(
            'sample_sim', config_basename)

        #FIXME: Insufficient. We only want to reuse this file if this file's
        #version is identical to that of the default YAML configuration file's
        #version. Hence, this logic should (arguably) be shifted elsewhere --
        #probably into "betse.science.sim_config".

        # If this file already exists, reuse this file.
        if files.is_file(self._args.config_filename):
            logs.log_info(
                'Reusing simulation configuration "{}".'.format(
                    config_basename))
        # Else, create this file.
        else:
            self._do_config()

        # Initialize and run such simulation.
        self._do_seed()
        self._do_init()
        self._do_sim()


    def _do_config(self) -> None:
        '''
        Run the `config` subcommand.
        '''

        # Delay importing this module, which imports heavy-weight dependencies.
        from betse.science.config import sim_config
        sim_config.write_default(self._args.config_filename)


    def _do_seed(self) -> None:
        '''
        Run the `seed` subcommand.
        '''
        self._get_sim_runner().makeWorld()


    def _do_init(self) -> None:
        '''
        Run the `init` subcommand.
        '''
        self._get_sim_runner().initialize()


    def _do_sim(self) -> None:
        '''
        Run the `sim` subcommand.
        '''
        self._get_sim_runner().simulate()


    def _do_plot(self) -> None:
        '''
        Run the `plot` subcommand.
        '''
        # Run such subcommand's passed subcommand. See _run() for details.
        subcommand_method_name = '_do_plot_' + self._args.subcommand_name_plot
        subcommand_method = getattr(self, subcommand_method_name)
        subcommand_method()


    def _do_plot_seed(self) -> None:
        '''
        Run the `plot` subcommand's `seed` subcommand.
        '''
        self._get_sim_runner().plotWorld()


    def _do_plot_init(self) -> None:
        '''
        Run the `plot` subcommand's `init` subcommand.
        '''
        self._get_sim_runner().plotInit()


    def _do_plot_sim(self) -> None:
        '''
        Run the `plot` subcommand's `sim` subcommand.
        '''
        self._get_sim_runner().plotSim()

    # ..................{ GETTERS                            }..................
    def _get_sim_runner(self):
        '''
        Get a new simulation runner configured with sane defaults.
        '''
        # Import from "betse.science" in a just-in-time manner. Why? This
        # importation imports heavy-weight dependencies and hence is slow.
        from betse.science.simrunner import SimRunner

        # Get such runner.
        return SimRunner(config_filename = self._args.config_filename)
