#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Concrete subclasses defining this application's command line interface (CLI).
'''

#FIXME: Refactor the class defined below to subclass "CLISubcommandableABC"
#instead, which should dramatically reduce the size and scope of this submodule.

# ....................{ IMPORTS                            }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: To raise human-readable exceptions on application startup, the
# top-level of this module may import *ONLY* from submodules guaranteed to:
# * Exist, including standard Python and application modules.
# * Never raise exceptions on importation (e.g., due to module-level logic).
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

from betse.cli import cliinfo
from betse.cli.api.cliabc import CLIABC
from betse.cli.api.clicmd import CLISubcommander
from betse.util.io.log import logs
from betse.util.path import files, pathnames
from betse.util.py import pyident, pys
from betse.util.type.call.memoizers import property_cached
from betse.util.type.obj import objects
from betse.util.type.types import ModuleType

# ....................{ SUBCLASS                           }....................
class BetseCLI(CLIABC):
    '''
    Command line interface (CLI) for this application.

    Attributes
    ----------
    _arg_parser_plot : ArgParserType
        Argument parser parsing arguments passed to the ``plot`` subcommand.
    '''

    # ..................{ INITIALIZERS                       }..................
    def __init__(self) -> None:

        # Initialize our superclass.
        super().__init__()

        # Nullify attributes for safety.
        self._arg_parser_plot = None

    # ..................{ SUPERCLASS ~ property              }..................
    @property
    def _help_epilog(self) -> str:

        return '''
subcommand help:
For help with specific subcommands, pass either the "-h" or "--help" argument to
the desired subcommand. For example, for help with both the "plot" subcommand
and that subcommand's "seed" subsubcommand, run:

;    {script_basename} plot --help
;    {script_basename} plot seed --help
'''


    @property
    def _module_ignition(self) -> ModuleType:

        from betse import ignition
        return ignition


    @property
    def _module_metadata(self) -> ModuleType:

        from betse import metadata
        return metadata

    # ..................{ SUPERCLASS ~ args                  }..................
    def _config_arg_parsing(self) -> None:

        # Container of all top-level argument subparsers for this application.
        subcommander_top = self._make_subcommander_top()
        subcommander_top.add(cli=self, arg_parser=self._arg_parser_top)

        # Argument parser parsing arguments passed to the "plot" subcommand.
        self._arg_parser_plot = (
            subcommander_top.subcommand_name_to_arg_parser['plot'])

    # ....................{ SUBCOMMANDS                        }....................
    #FIXME: Document us up.
    def _make_subcommander_top(self) -> CLISubcommander:

        # Avoid circular import dependencies.
        from betse.cli.api.clicmd import (
            CLISubcommandNoArg, CLISubcommandParent, CLISubcommandYAMLOnly,)

        # Return a container of all top-level subcommands accepted by this
        # application's CLI command.
        return CLISubcommander(
            subcommand_var_name='subcommand_name_top',
            help_title='subcommands',
            subcommands=(

                CLISubcommandYAMLOnly(
                    name='config',
                    help_synopsis=(
                        'create a default config file for '
                        '{program_name} simulations'
                    ),
                    help_description='''
Write a default tissue simulation configuration to the passed output file. While
not strictly necessary, this file should have filetype ".yaml" . If this file
already exists, an error will be printed.

You may edit this file at any time. By default, this file instructs
{program_name} to save simulation results (e.g., plots) to the directory
containing this file.
''',),


                CLISubcommandYAMLOnly(
                    name='seed',
                    help_synopsis='seed a new cell cluster for a config file',
                    help_description='''
Create the cell cluster defined by the passed configuration file. The results
will be saved to output files defined by this configuration.
''',),


                CLISubcommandYAMLOnly(
                    name='init',
                    help_synopsis=(
                        'initialize a seeded cell cluster for a config file'),
                    help_description='''
Initialize (i.e., calculate steady-state concentrations for) the previously
created cell cluster defined by the passed configuration file. Initialization
results will be saved to output files defined by this configuration, while the
previously created cell cluster will be loaded from input files defined by this
configuration.
''',),


                CLISubcommandYAMLOnly(
                    name='sim',
                    help_synopsis=(
                        'simulate an initialized cell cluster '
                        'for a config file'
                    ),
                    help_description='''
Simulate the previously initialized cell cluster defined by the passed
configuration file. Simulation results will be saved to output files defined by
this configuration, while the previously initialized cell cluster will be loaded
from input files defined by this configuration.
''',),


                CLISubcommandYAMLOnly(
                    name='sim-grn',
                    help_synopsis=(
                        'simulate a gene regulatory network for a config file'),
                    help_description='''
Simulate a gene regulatory network (GRN) for the previously initialized cell
cluster defined by the passed configuration file, whose "gene regulatory network
config" option specifies the path of the configuration file defining this
network. All other simulation features and options will be ignored.

Simulation results will be saved to output files defined by this configuration,
while the previously initialized cell cluster will be loaded from input files
defined by this configuration.
''',),


                CLISubcommandParent(
                    name='plot',
                    help_synopsis=(
                        'plot a seeded, initialized, or simulated simulation'),
                    help_description='''
Run the passed plotting subcommand. For example, to plot the previous
simulation defined by a configuration file "my_sim.yaml" in the current
directory:

;    betse plot sim my_sim.yaml
''',
                    subcommander=self._make_subcommander_plot(),
                ),


                CLISubcommandNoArg(
                    name='repl',
                    help_synopsis=(
                        'enter an interactive {program_name}-aware REPL'),
                    help_description='''
Initialize the {program_name} environment and immediately open a
Read-Evaluate-Print Loop (REPL). This allows interactive manipulation of the
simulations and analyses.
''',),


                CLISubcommandNoArg(
                    name='info',
                    help_synopsis=(
                        'print metadata synopsizing '
                        '{program_name} and current system'
                    ),
                    help_description='''
Print informational metadata in ":"-delimited key-value format, including:

* Program name, version, and principal authors.

* Absolute paths of critical files and directories used by {program_name},
including:

* {program_name}'s data directory (i.e., the program-specific directory to
    which non-Python files intended for use by external users are stored).

* {program_name}'s dot directory (i.e., the user-specific directory to which
    files and directories intended for internal program use are stored).

* {program_name}'s log file (i.e., the user-specific file to which all runtime
    messages are appended, including low-level debug statements, non-fatal
    warnings, and fatal errors).
''',),


                CLISubcommandNoArg(
                    name='try',
                    help_synopsis=(
                        'create, init, simulate, and plot a sample simulation'),
                    help_description='''
Run a sample tissue simulation. This subcommand (A) creates a default YAML
configuration file, (B) creates the cell cluster defined by that file, and
(C) initializes, (D) simulates, and (E) plots the tissue simulation defined by
that file given that cluster. All files and directories created by these
operations will be preserved (rather than deleted on subcommand completion).

Equivalently, this subcommand is shorthand for the following:

;    betse config   sample_sim/sample_sim.yaml
;    betse seed     sample_sim/sample_sim.yaml
;    betse init     sample_sim/sample_sim.yaml
;    betse sim      sample_sim/sample_sim.yaml
;    betse plot sim sample_sim/sample_sim.yaml
''',),

            ))


    def _make_subcommander_plot(self) -> CLISubcommander:
        '''
        Sequence of all :class:`CLISubcommandABC` instances defining the
        low-level subcommands accepted by the top-level ``plot`` subcommand
        accepted by this application.
        '''

        # Avoid circular import dependencies.
        from betse.cli.api.clicmd import CLISubcommandYAMLOnly

        # Return a container of all child subcommands accepted by this parent
        # "plot" subcommand.
        return CLISubcommander(
            subcommand_var_name='subcommand_name_plot',
            help_title='plot subcommands',
            subcommands=(
                CLISubcommandYAMLOnly(
                    name='seed',
                    help_synopsis=(
                        'plot a seeded cell cluster for a config file'),
                    help_description='''
Plot the previously seeded cell cluster defined by the passed configuration
file. Plot results will be saved to output files defined by this configuration,
while the previously seeded cell cluster will be loaded from input files
defined by this configuration.
''',),


                CLISubcommandYAMLOnly(
                    name='init',
                    help_synopsis=(
                        'plot an initialized cell cluster for a config file'),
                    help_description='''
Plot the previously initialized cell cluster defined by the passed configuration
file. Plot results will be saved to output files defined by this configuration,
while the previously initialized cell cluster will be loaded from input files
defined by this configuration.
''',),


                CLISubcommandYAMLOnly(
                    name='sim',
                    help_synopsis=(
                        'plot a simulated cell cluster for a config file'),
                    help_description='''
Plot the previously simulated cell cluster defined by the passed configuration
file. Plot results will be saved to output files defined by this configuration,
while the previously simulated cell cluster will be loaded from input files
defined by this configuration.
''',
            ),


                CLISubcommandYAMLOnly(
                    name='sim-grn',
                    help_synopsis=(
                        'plot a simulated gene regulatory network '
                        'for a config file'
                    ),
                    help_description='''
Plot the previously simulated gene regulatory network (GRN) defined by the
passed configuration file. Plot results will be saved to output files defined by
this configuration, while the previously simulated cell cluster will be loaded
from input files defined by this configuration.
''',),

            ))

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
            self._arg_parser_top.print_help()

            # Return the current instance of this object. While trivial, this
            # behaviour simplifies memory profiling of this object.
            return self
        # Else, a subcommand was passed.

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
