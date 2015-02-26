#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Help strings to be printed by `betse`'s command line interface (CLI).
'''

# ....................{ IMPORTS                            }....................
from betse.util.path import paths

# ....................{ TEMPLATES                          }....................
TEMPLATE_SUBCOMMANDS = '''
Exactly one of the following subcommands must be passed:
'''
'''
Help string template for the set of subcommands.
'''

TEMPLATE_EPILOG = '''
For further help on subcommand behaviour and expected arguments, pass the "-h"
or "--help" argument to any subcommand. For example, for further help with the
"sim" subcommand, run:

;    {script_basename} sim --help
'''
'''
Help string template for *epilog* (i.e., string printed after *all* other text).
'''

# ....................{ TEMPLATES ~ command                }....................
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
Help string template for the `info` command.
'''

TEMPLATE_SUBCOMMAND_TRY = '''
Run a sample tissue simulation by automatically creating a default configuration
file and initializing, running, and plotting the simulation configured by such
file. This command is shorthand for the following:

;    {script_basename} sim cfg init run plot my_sim.yaml
'''
'''
Help string template for the `try` command.
'''

TEMPLATE_SUBCOMMAND_SIM = '''
Run the passed tissue simulation subcommand(s) configured by the passed
configuration file. For example, to initialize, run, and plot a simulation
configured by an existing file "my_sim.yaml" in the current directory:

;    {script_basename} sim init run plot my_sim.yaml

The last passed argument must be the name of a simulation configuration file. If
the "cfg" subcommand is run, such file need not already exist; else, such file
must already exist (e.g., by a previous run of the "cfg" subcommand).

All other passed arguments are the subcommands to be run. Since subcommands will
be run in the passed order, order is significant. At least one of the following
subcommands must be passed:

;cfg
;----------
Write a default tissue simulation configuration to the passed output file. Such
filename should (ideally) be suffixed by ".yaml". By default, this file will
instruct {program_name} to save simulation output (e.g., plots) to the same
directory containing this file. This file may be freely edited at any time.

;init
;----------
Initialize the tissue simulation specified by the passed configuration file.
Initialization results will be saved to the output file specified in such
configuration.

;run
;----------
Run the previously initialized tissue simulation specified by the passed
configuration file. Simulation results will be saved to the output file
specified in such configuration. Likewise, the previously run initialization
will be loaded from the input file specified in such configuration. (If such
file does not exist, an error is raised.)

;plot
;----------
Plot the previously run tissue simulation specified by the passed configuration
file. Plot results will be saved to the output file(s) specified in such
configuration. Likewise, the previously run simulation will be loaded from the
input file specified in such configuration. (If such file does not exist, an
error is raised.)
''' #.format()
'''
Help string template for the `sim` command.
'''

# --------------------( WASTELANDS                         )--------------------
#.format(
    # simulation_config_basename = paths.get_basename(
    #     pathtree.SIMULATION_CONFIG_DEFAULT_FILENAME))
# from betse import pathtree
#Subcommand to be run.
