#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Help strings to be printed by `betse`'s command line interface (CLI).
'''

# ....................{ TEMPLATES ~ subcommands            }....................
TEMPLATE_SUBCOMMANDS_PREFIX = '''
Exactly one of the following subcommands must be passed:
'''
'''
Help string template for the set of subcommands.
'''

TEMPLATE_SUBCOMMANDS_SUFFIX = '''
subcommand help:

For help with specific subcommands, pass the "-h" or "--help" arguments to each.
For example, for help with both the "plot" subcommand and that subcommand's
"world" subsubcommand:

;    betse plot --help
;    betse plot world --help
'''
'''
Help string template for the **program epilog** (i.e., string printed after
*all* other text in top-level help output).
'''

# ....................{ TEMPLATES ~ subcommand             }....................
TEMPLATE_SUBCOMMAND_INFO = '''
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
'''
'''
Help string template for the `info` subcommand.
'''

TEMPLATE_SUBCOMMAND_TRY = '''
Run a sample tissue simulation. This subcommand (A) creates a default YAML
configuration file, (B) creates the cellular world defined by that file, and
(C) initializes, (D) runs, and (E) plots the tissue simulation defined by that
file within that world. For safety, this subcommand preserves (i.e., does not
delete on completion) files and directories created by these operations.

Equivalently, this is shorthand for the following:

;    betse config   sample_sim/sample_sim.yaml
;    betse world    sample_sim/sample_sim.yaml
;    betse prep     sample_sim/sample_sim.yaml
;    betse run      sample_sim/sample_sim.yaml
;    betse plot run sample_sim/sample_sim.yaml
'''
'''
Help string template for the `try` subcommand.
'''

TEMPLATE_SUBCOMMAND_CONFIG = '''
Write a default tissue simulation configuration to the passed output file. While
not strictly necessary, such file should have filetype ".yaml" . If such file
already exists, an error will be printed.

You may freely edit this file at any time. Note that the default configuration
instructs {program_name} to save simulation output (e.g., plots) into the
directory in which this file resides.
'''
'''
Help string template for the `config` subcommand.
'''

TEMPLATE_SUBCOMMAND_WORLD = '''
Create the cellular world defined by the passed configuration file. The results
will be saved to output files defined by such configuration.

If this configuration file does not exist, an error will be printed.
'''
'''
Help string template for the `world` subcommand.
'''

TEMPLATE_SUBCOMMAND_PREP = '''
Prep (i.e., calculate steady-state concentrations for) the previously created
cellular world defined by the passed configuration file. Prep results will be
saved to output files defined by such configuration. Likewise, the previously
created world will be loaded from input files defined by this configuration.

If this configuration file does not exist, an error will be printed.
'''
'''
Help string template for the `prep` subcommand.
'''

TEMPLATE_SUBCOMMAND_RUN = '''
Run the previously prepped tissue simulation defined by the passed configuration
file. Simulation results will be saved to output files defined by this
configuration. Likewise, the previously prepped simulation will be loaded from
input files defined by this configuration.

If this configuration file does not exist, an error will be printed.
'''
'''
Help string template for the `run` subcommand.
'''

# ....................{ TEMPLATES ~ subcommand : plot      }....................
TEMPLATE_SUBCOMMAND_PLOT = '''
Run the passed plotting subcommand. For example, to plot the previously run
simulation defined by a configuration file "my_sim.yaml" in the current
directory:

;    betse plot run my_sim.yaml
'''
'''
Help string template for the `plot` subcommand.
'''

TEMPLATE_SUBCOMMAND_PLOT_WORLD = '''
Plot the previously created cellular world defined by the passed configuration
file. Plot results will be saved to output files defined by this configuration.
Likewise, the previously created world will be loaded from input files defined
by this configuration.

If this configuration file does not exist, an error will be printed.
'''
'''
Help string template for the `plot` subcommand's `world` subcommand.
'''

TEMPLATE_SUBCOMMAND_PLOT_PREP = '''
Plot the previously prepped tissue simulation defined by the passed
configuration file. Plot results will be saved to output files defined by this
configuration. Likewise, the previously prepped simulation will be loaded from
input files defined by this configuration.

If this configuration file does not exist, an error will be printed.
'''
'''
Help string template for the `plot` subcommand's `prep` subcommand.
'''

TEMPLATE_SUBCOMMAND_PLOT_RUN = '''
Plot the previously run tissue simulation defined by the passed configuration
file. Plot results will be saved to output files defined by this configuration.
Likewise, the previously run simulation will be loaded from input files defined
by this configuration.

If this configuration file does not exist, an error will be printed.
'''
'''
Help string template for the `plot` subcommand's `run` subcommand.
'''

# --------------------( WASTELANDS                         )--------------------
# TEMPLATE_SUBCOMMAND_SIM = '''
# Run the passed tissue simulation subcommand. For example, to initialize the
# simulation defined by an existing configuration file "my_sim.yaml" in the
# current directory:
#
# ;    betse setup run my_sim.yaml
# '''
# '''
# Help string template for the `sim` subcommand.
# '''

# For further help with subcommand behavior and expected arguments, pass the "-h"
# Template with which to create **epilog templates** (i.e., strings printed after
# *all* other text at some level of help output).
# TEMPLATE_SUBCOMMANDS_TOP_SUFFIX = TEMPLATE_EPILOG_TEMPLATE.format(
#     example_subcommand_name = 'sim',
#     example_subcommand = 'sim --help')
# '''
# Help string template for the **program epilog** (i.e., string printed after
# *all* other text in top-level help output).
# '''
#
# TEMPLATE_SUBCOMMANDS_SIM_SUFFIX = TEMPLATE_EPILOG_TEMPLATE.format(
#     example_subcommand_name = 'cfg',
#     example_subcommand = 'sim cfg --help')
# '''
# Help string template for the `sim` subcommand's epilog.
# '''

#'''
# For further help on subcommand behaviour and expected arguments, pass the "-h"
# or "--help" argument to any subcommand. For example, for further help with the
# "sim" subcommand, run:
#
# ;    {script_basename} sim --help
# '''
# from betse.util.path import paths
#.format(
    # simulation_config_basename = paths.get_basename(
    #     pathtree.SIMULATION_CONFIG_DEFAULT_FILENAME))
# from betse import pathtree
#Subcommand to be run.

# TEMPLATE_SUBCOMMAND_SIM = '''
# Run the passed tissue simulation subcommand(s) configured by the passed
# configuration file. For example, to initialize, run, and plot a simulation
# configured by an existing file "my_sim.yaml" in the current directory:
#
# ;    {script_basename} sim init run plot_run my_sim.yaml
#
# The last passed argument must be the name of a simulation configuration file. If
# the "cfg" subcommand is run, such file need not already exist; else, such file
# must already exist (e.g., by a previous run of the "cfg" subcommand).
#
# All other passed arguments are the subcommands to be run. Since subcommands will
# be run in the passed order, order is significant. At least one of the following
# subcommands must be passed:
#
# ;cfg
# ;----------
# Write a default tissue simulation configuration to the passed output file. Such
# filename should (ideally) be suffixed by ".yaml". By default, this file will
# instruct {program_name} to save simulation output (e.g., plots) to the same
# directory containing this file. This file may be freely edited at any time.
#
# ;init
# ;----------
# Initialize the tissue simulation defined by the passed configuration file.
# Initialization results will be saved to the output file defined by such
# configuration.
#
# ;run
# ;----------
# Run the previously initialized tissue simulation defined by the passed
# configuration file. Simulation results will be saved to the output file
# defined by this configuration. Likewise, the previously run initialization
# will be loaded from the input file defined by this configuration. (If such
# file does not exist, an error will be printed.)
#
# ;plot_init
# ;----------
# Plot the previously initialized tissue simulation defined by the passed
# configuration file. Plot results will be saved to the output file(s) defined
# in this configuration. Likewise, the previously run initialization will be
# loaded from the input file defined by this configuration. (If such file does
# not exist, an error will be printed.)
#
# ;plot_run
# ;----------
# Plot the previously run tissue simulation defined by the passed configuration
# file. Plot results will be saved to the output file(s) defined by such
# configuration. Likewise, the previously run simulation will be loaded from the
# input file defined by this configuration. (If such file does not exist, an
# error will be printed.)
# '''
# '''
# Help string template for the `sim` command.
# '''
