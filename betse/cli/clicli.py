#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
BETSE's command line interface (CLI).
'''

# ....................{ IMPORTS                            }....................
from argparse import ArgumentParser, _SubParsersAction
from betse.cli import clisubcommands, info
from betse.cli.cliabc import expand_help, CLIABC
from betse.cli.clisubcommands import CLISubcommand
from betse.exceptions import BetseTestException
from betse.util.io.log import logs
from betse.util.path import files, paths
from betse.util.py import pys
from betse.util.type.types import type_check

# ....................{ CLASS                              }....................
class CLICLI(CLIABC):
    '''
    BETSE's command line interface (CLI).

    Attributes
    ----------
    _arg_parser_plot : ArgumentParser
        Subparser parsing arguments passed to the `plot` subcommand.
    _arg_subparsers_top : ArgumentParser
        Subparsers parsing top-level subcommands (e.g., `plot`).
    _arg_subparsers_plot : ArgumentParser
        Subparsers parsing `plot` subcommands (e.g., `plot seed`).
    '''

    # ..................{ INITIALIZERS                       }..................
    def __init__(self):
        super().__init__()

        # Nullify attributes for safety.
        self._arg_parser_plot = None
        self._arg_subparsers_top = None
        self._arg_subparsers_plot = None

    # ..................{ SUPERCLASS ~ args                  }..................
    def _get_arg_parser_top_kwargs(self):
        # Keyword arguments passed to the top-level argument parser constructor.
        return {
            'epilog': clisubcommands.SUBCOMMANDS_SUFFIX,
        }


    def _config_arg_parsing(self):
        # Collection of top-level argument subparsers.
        self._arg_subparsers_top = self._arg_parser.add_subparsers(
            # Name of the attribute storing the passed subcommand name.
            dest='subcommand_name_top',

            # Title of the subcommand section in help output.
            title='subcommands',

            # Description to be printed *BEFORE* subcommand help.
            description=expand_help(clisubcommands.SUBCOMMANDS_PREFIX),
        )

        # Dictionary mapping from top-level subcommand name to the argument
        # subparser parsing this subcommand.
        subcommand_name_to_subparser = {}

        # For each top-level subcommand...
        for subcommand in clisubcommands.SUBCOMMANDS:
            #FIXME: Refactor this logic as follows:
            #
            #* The local "subcommand_name_to_subparser" dicitonary should be
            #  replaced by performing the following in the body of this loop:
            #  * If a global "help.SUBCOMMANDS_{subcommand.name}" dictionary
            #    exists, non-recursively loop over that dictionary in the same
            #    manner as well. Hence, this generalizes support of subcommand
            #    subcommands. Nice!
            #* Likewise, if a "self._arg_parser_{subcommand.name}" instance
            #  variable exists, set that variable to this parser. Such logic
            #  should, arguably, be performed by _add_subcommand().

            # Add an argument subparser parsing this subcommand.
            subcommand_name_to_subparser[subcommand.name] = self._add_subcommand(
                subcommand=subcommand, arg_subparsers=self._arg_subparsers_top)

        # Configure arg parsing for subcommands of the "plot" subcommand.
        self._arg_parser_plot = subcommand_name_to_subparser['plot']
        self._config_arg_parsing_plot()

    # ..................{ SUBCOMMAND ~ plot                  }..................
    def _config_arg_parsing_plot(self) -> None:
        '''
        Configure argument parsing for subcommands of the `plot` subcommand.
        '''

        #FIXME: Refactor into iteration resembling that performed by the
        #_config_arg_parsing() method.

        # Collection of all subcommands of the "plot" subcommand.
        self._arg_subparsers_plot = self._arg_parser_plot.add_subparsers(
            # Name of the attribute storing the passed subcommand name.
            dest='subcommand_name_plot',

            # Title of the subcommand section in help output.
            title='plot subcommands',

            # Description to be printed *BEFORE* subcommand help.
            description=expand_help(clisubcommands.SUBCOMMANDS_PREFIX),
        )

        # For each subcommand of the "plot" subcommand...
        for subcommand in clisubcommands.SUBCOMMANDS_PLOT:
            # Add an argument subparser parsing this subcommand.
            self._add_subcommand(
                subcommand=subcommand, arg_subparsers=self._arg_subparsers_plot)

    # ..................{ SUBCOMMANDS                        }..................
    @type_check
    def _add_subcommand(
        self,
        subcommand: CLISubcommand,
        arg_subparsers: _SubParsersAction,
        *args, **kwargs
    ) -> ArgumentParser:
        '''
        Create a new **argument subparser** (i.e., an `argparse`-specific object
        parsing command-line arguments) for the passed subcommand, add this
        subparser to the passed collection of **argument subparsers** (i.e.,
        another `argparse`-specific object cotaining multiple subparsers), and
        return this subparser.

        This subparser will be configured to:

        * If this subcommand accepts a configuration filename, require such an
          argument be passed.
        * Else, require no arguments be passed.

        Parameters
        ----------
        subcommand: CLISubcommand
            Subcommand to be added.
        arg_subparsers : _SubParsersAction
            Collection of sibling subcommand argument parsers to which the
            subcommand argument parser created by this method is added. This
            collection is owned either by:
            * A top-level subcommand (e.g., `plot`), in which case the
              subcommand created by this method is a child of that subcommand.
            * No subcommand, in which case the subcommand created by this method
              is a top-level subcommand.

        All remaining positional and keyword arguments are passed as is to this
        subparser's `__init__()` method.

        Returns
        ----------
        ArgumentParser
            Subcommand argument parser created by this method.
        '''

        # Extend the passed dictionary of keyword arguments with (in any order):
        #
        # * The dictionary of globally applicable keyword arguments.
        # * The dictionary of subcommand-specific keyword arguments.
        kwargs.update(self._arg_parser_kwargs)
        kwargs.update(subcommand.get_arg_parser_kwargs())

        # Subcommand argument subparser added to this collection.
        arg_subparser = arg_subparsers.add_parser(*args, **kwargs)

        # If this subcommand requires a configuration file, configure this
        # subparser accordingly.
        if subcommand.is_passed_yaml:
            arg_subparser.add_argument(
                'config_filename',
                metavar='CONFIG_FILE',
                help='simulation configuration file',
            )

        # Return this subparser.
        return arg_subparser

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
        # Sanitized name of this subcommand.
        subcommand_name_top = clisubcommands.sanitize_name(
            self._args.subcommand_name_top)

        # Name of the method running this subcommand.
        subcommand_method_name = '_do_' + subcommand_name_top

        # Method running this subcommand. If this method does *NOT* exist,
        # getattr() will raise a non-human-readable exception. Usually, that
        # would be bad. In this case, however, argument parsing coupled with a
        # reliable class implementation guarantees this method to exist.
        subcommand_method = getattr(self, subcommand_method_name)

        # Run this subcommand.
        subcommand_method()

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

        # Run all general-purposes phases, thus excluding network-specific
        # phases (e.g., "_do_sim_brn", "_do_sim_grn"), in the expected order.
        self._do_seed()
        self._do_init()
        self._do_sim()
        self._do_plot_seed()
        self._do_plot_init()
        self._do_plot_sim()


    def _do_config(self) -> None:
        '''
        Run the `config` subcommand.
        '''

        # Avoid importing modules importing dependencies at the top level.
        from betse.science.config import confdefault
        confdefault.write(self._args.config_filename)


    def _do_seed(self) -> object:
        '''
        Run the `seed` subcommand and return the result of doing so.
        '''

        return self._get_sim_runner().seed()


    def _do_init(self) -> object:
        '''
        Run the `init` subcommand and return the result of doing so.
        '''

        return self._get_sim_runner().init()


    def _do_sim(self) -> object:
        '''
        Run the `sim` subcommand and return the result of doing so.
        '''

        return self._get_sim_runner().sim()


    def _do_sim_brn(self) -> object:
        '''
        Run the `sim-brn` subcommand and return the result of doing so.
        '''

        return self._get_sim_runner().sim_brn()


    def _do_sim_grn(self) -> object:
        '''
        Run the `sim-grn` subcommand and return the result of doing so.
        '''

        return self._get_sim_runner().sim_grn()


    def _do_plot(self) -> object:
        '''
        Run the `plot` subcommand and return the result of doing so.
        '''

        # If no subcommand was passed, print help output and return. Note that
        # this does *NOT* constitute a fatal error.
        if not self._args.subcommand_name_plot:
            print()
            self._arg_parser_plot.print_help()
            return

        # Run this subcommand's passed subcommand and return the result
        # of doing so. See _run() for details.
        subcommand_name_plot = clisubcommands.sanitize_name(
            self._args.subcommand_name_plot)
        subcommand_method_name = '_do_plot_' + subcommand_name_plot
        subcommand_method = getattr(self, subcommand_method_name)
        return subcommand_method()


    def _do_plot_seed(self) -> object:
        '''
        Run the `plot` subcommand's `seed` subcommand and return the result
        of doing so.
        '''

        return self._get_sim_runner().plot_seed()


    def _do_plot_init(self) -> object:
        '''
        Run the `plot` subcommand's `init` subcommand and return the result
        of doing so.
        '''

        return self._get_sim_runner().plot_init()


    def _do_plot_sim(self) -> object:
        '''
        Run the `plot` subcommand's `sim` subcommand and return the result
        of doing so.
        '''

        return self._get_sim_runner().plot_sim()


    def _do_plot_sim_brn(self) -> object:
        '''
        Run the `plot` subcommand's `sim-brn` subcommand and return the result
        of doing so.
        '''

        return self._get_sim_runner().plot_brn()


    def _do_plot_sim_grn(self) -> object:
        '''
        Run the `plot` subcommand's `sim-grn` subcommand and return the result
        of doing so.
        '''

        return self._get_sim_runner().plot_grn()


    def _do_repl(self) -> None:
        '''
        Run the `repl` subcommand.
        '''

        # In the unlikely edge-case of the "repl" subcommand being erroneously
        # run by a functional test, prohibit this by raising an exception.
        # Permitting this would probably cause tests to indefinitely hang.
        if pys.is_testing():
            raise BetseTestException('REPL unavailable for testing.')

        # Defer heavyweight imports until *AFTER* possibly failing above.
        from betse.cli.repl import repls

        # Start the desired REPL.
        repls.start_repl()

    # ..................{ GETTERS                            }..................
    #FIXME: Inefficient, particularly when running _do_try(). Instead, refactor
    #this method into a new "_sim_runner" property cached by @property_cached,
    #ensuring that one and only one SimRunner is instantiated for each CLI
    #instance.
    def _get_sim_runner(self):
        '''
        BETSE simulation runner preconfigured with sane defaults.
        '''

        # Defer heavyweight imports.
        from betse.science.simrunner import SimRunner

        # Return this runner.
        return SimRunner(config_filename=self._args.config_filename)
