#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Concrete subclasses defining this application's command line interface (CLI).
'''

# ....................{ IMPORTS                           }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: To raise human-readable exceptions on application startup, the
# top-level of this module may import *ONLY* from submodules guaranteed to:
# * Exist, including standard Python and application modules.
# * Never raise exceptions on importation (e.g., due to module-level logic).
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

from betse.cli import cliinfo
from betse.util.cli.clicmd import (
    CLISubcommander,
    CLISubcommandNoArg,
    CLISubcommandParent,
    CLISubcommandYAMLOnly,
)
from betse.util.cli.clicmdabc import CLISubcommandableABC
from betse.util.cli.cliopt import CLIOptionBoolTrue
from betse.util.io.log import logs
from betse.util.os import displays
from betse.util.path import files, pathnames
from betse.util.type.decorator.decmemo import property_cached
from betse.util.type.types import SequenceTypes

# ....................{ SUBCLASS                          }....................
class BetseCLI(CLISubcommandableABC):
    '''
    Command line interface (CLI) for this application.

    Attributes (of :attr:`_args`)
    ----------
    is_headless : bool
        ``True`` only if the active Python interpreter is to be coerced into
        running **headless** (i.e., with *no* access to a GUI display).
        Defaults to ``False``.
    '''

    # ..................{ SUPERCLASS ~ property : commands  }..................
    @property_cached
    def _subcommander_top(self) -> CLISubcommander:

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
                        'plot an initialized or simulated simulation'),
                        # 'plot a seeded, initialized, or simulated simulation'),
                    help_description='''
Run the passed plotting subcommand. For example, to plot the previous
simulation defined by a configuration file "my_sim.yaml" in the current
directory:

;    betse plot sim my_sim.yaml
''',
                    subcommander=self._subcommander_plot,
                ),


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
configuration file, (B) creates the cell cluster defined by that file, (C)
initializes that cell cluster, (D) plots that initialization, (E) simulates that
initialization, and (F) plots that simulation. All files and directories created
by these operations will be preserved (rather than deleted on subcommand
completion).

Equivalently, this subcommand is shorthand for the following:

;    betse config    sample_sim/sample_sim.yaml
;    betse seed      sample_sim/sample_sim.yaml
;    betse init      sample_sim/sample_sim.yaml
;    betse sim       sample_sim/sample_sim.yaml
;    betse plot init sample_sim/sample_sim.yaml
;    betse plot sim  sample_sim/sample_sim.yaml
''',),

            ))


    @property
    def _subcommander_plot(self) -> CLISubcommander:
        '''
        Container of all child subcommands accepted by the ``plot`` subcommand.
        '''

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

    # ..................{ SUPERCLASS ~ property : options   }..................
    @property
    def _options_top(self) -> SequenceTypes:

        # Sequence of all default top-level options.
        options_top = super()._options_top

        # Return a list extending this sequence with subclass-specific options.
        return options_top + [
            CLIOptionBoolTrue(
                long_name='--headless',
                synopsis=(
                    'enable headless mode '
                    '(display no plots or animations) '
                    '[forces "--matplotlib-backend=agg"]'
                ),
            ),
        ]


    def _parse_options_top(self) -> None:

        # Parse all default top-level options.
        super()._parse_options_top()

        # If coercing this Python interpreter into running headless, do so.
        if self._args.is_headless:
            displays.set_headless(True)
        # Else, allow the displays_is_headless() function to implicitly detect
        # whether this interpreter is actually running headless.

    # ..................{ SUPERCLASS ~ help                 }..................
    def _log_header(self) -> None:

        cliinfo.log_header()

    # ..................{ PROPERTIES                        }..................
    @property_cached
    def _sim_runner(self) -> 'betse.science.simrunner.SimRunner':
        '''
        Simulation runner running simulation subcommands on the YAML-formatted
        simulation configuration file passed by the user at the CLI.
        '''

        # Defer heavyweight imports.
        from betse.science.parameters import Parameters
        from betse.science.simrunner import SimRunner

        # Simulation configuration loaded from this YAML-formatted file.
        p = Parameters.make(conf_filename=self._args.conf_filename)

        # Create and return a simulation runner for this configuration.
        return SimRunner(p=p)

    # ..................{ SUBCOMMANDS                       }..................
    def _do_info(self) -> None:
        '''
        Run the ``info`` subcommand.
        '''

        cliinfo.log_info()

    # ..................{ SUBCOMMANDS ~ sim                 }..................
    def _do_try(self) -> object:
        '''
        Run the ``try`` subcommand and return the result of doing so.
        '''

        # Basename of the sample configuration file to be created.
        config_basename = 'sample_sim.yaml'

        # Relative filename of this file, relative to the current directory.
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

        # Defer heavyweight imports.
        from betse.science.parameters import Parameters

        # Copy the default simulation configuration to this target file.
        Parameters().copy_default(trg_conf_filename=self._args.conf_filename)


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
