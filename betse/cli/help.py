#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Help strings printed by `betse`'s command line interface (CLI).
'''

# ....................{ TEMPLATES ~ subcommands            }....................
OPTION_VERSION = '''
print program version and exit
'''
'''
Help string template synopsizing the `--version` option.
'''


OPTION_VERBOSE = '''
print low-level debugging messages
'''
'''
Help string template synopsizing the `--verbose` option.
'''


OPTION_LOG_TYPE = '''
type of logging to perform (defaults to "{default}"):
;* "none", logging to stdout and stderr
;* "file", logging to the "--log-file" file
'''
'''
Help string template synopsizing the `--log-type` option.
'''


OPTION_LOG_FILE = '''
file to log to if "--log-type" is "file" (defaults to "{default}")
'''
'''
Help string template synopsizing the `--log-file` option.
'''


OPTION_PROFILE_TYPE = '''
type of profiling to perform (defaults to "{default}"):
;* "none", disabling profiling
;* "call", profiling method and function calls
'''
'''
Help string template synopsizing the `--profile-type` option.
'''


OPTION_PROFILE_FILE = '''
file to profile to if "--profile-type" is "call" (defaults to "{default}")
'''
'''
Help string template synopsizing the `--profile-file` option.
'''

# ....................{ TEMPLATES ~ subcommands            }....................
SUBCOMMANDS_PREFIX = '''
Exactly one of the following subcommands must be passed:
'''
'''
Help string template preceding the list of all subcommands.
'''

SUBCOMMANDS_SUFFIX = '''
subcommand help:

For help with specific subcommands, either pass the "-h" or "--help" argument to
the desired subcommand. For example, for help with both the "plot" subcommand
and that subcommand's "seed" subsubcommand:

;    betse plot --help
;    betse plot seed --help
'''
'''
Help string template for the **program epilog** (i.e., string printed after
*all* other text in top-level help output).
'''

# ....................{ TEMPLATES ~ subcommand             }....................
SUBCOMMAND_INFO = '''
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


SUBCOMMAND_TRY = '''
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
'''
'''
Help string template for the `try` subcommand.
'''


SUBCOMMAND_CONFIG = '''
Write a default tissue simulation configuration to the passed output file. While
not strictly necessary, this file should have filetype ".yaml" . If this file
already exists, an error will be printed.

You may edit this file at any time. By default, this file instructs
{program_name} to save simulation results (e.g., plots) to the directory
containing this file.
'''
'''
Help string template for the `config` subcommand.
'''


SUBCOMMAND_SEED = '''
Create the cell cluster defined by the passed configuration file. The results
will be saved to output files defined by this configuration.
'''
'''
Help string template for the `seed` subcommand.
'''


SUBCOMMAND_INIT = '''
Initialize (i.e., calculate steady-state concentrations for) the previously
created cell cluster defined by the passed configuration file. Initialization
results will be saved to output files defined by this configuration, while the
previously created cell cluster will be loaded from input files defined by this
configuration.
'''
'''
Help string template for the `init` subcommand.
'''


SUBCOMMAND_SIM = '''
Simulate the previously initialized cell cluster defined by the passed
configuration file. Simulation results will be saved to output files defined by
this configuration, while the previously initialized cell cluster will be loaded
from input files defined by this configuration.
'''
'''
Help string template for the `sim` subcommand.
'''

# ....................{ TEMPLATES ~ subcommand : plot      }....................
SUBCOMMAND_PLOT = '''
Run the passed plotting subcommand. For example, to plot the previous
simulation defined by a configuration file "my_sim.yaml" in the current
directory:

;    betse plot run my_sim.yaml
'''
'''
Help string template for the `plot` subcommand.
'''


SUBCOMMAND_PLOT_SEED = '''
Plot the previously created cell cluster defined by the passed configuration
file. Plot results will be saved to output files defined by this configuration,
while the previously created cell cluster will be loaded from input files
defined by this configuration.
'''
'''
Help string template for the `plot` subcommand's `seed` subcommand.
'''


SUBCOMMAND_PLOT_INIT = '''
Plot the previously initialized cell cluster defined by the passed configuration
file. Plot results will be saved to output files defined by this configuration,
while the previously initialized cell cluster will be loaded from input files
defined by this configuration.
'''
'''
Help string template for the `plot` subcommand's `init` subcommand.
'''


SUBCOMMAND_PLOT_SIM = '''
Plot the previously simulated cell cluster defined by the passed configuration
file. Plot results will be saved to output files defined by this configuration,
while the previously simulated cell cluster will be loaded from input files
defined by this configuration.
'''
'''
Help string template for the `plot` subcommand's `sim` subcommand.
'''
