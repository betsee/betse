#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Metadata describing subcommands accepted by BETSE's command line interface
(CLI).
'''

# ....................{ IMPORTS                            }....................
from betse.cli.cliabc import expand_help
from betse.util.type.types import type_check

# ....................{ FUNCTIONS                          }....................
@type_check
def sanitize_name(subcommand_name: str) -> str:
    '''
    Sanitize the passed subcommand name (e.g., from `sim-gnr` to `sim_gnr`).

    Specifically, this utility function:

    * Replaces all hyphens in this name with underscores, as method names
        are generated from subcommand names but cannot contain hyphens.
    '''

    return subcommand_name.replace('-', '_')

# ....................{ CLASSES                            }....................
class CLISubcommand(object):
    '''
    Metadata encapsulating a **CLI subcommand** (i.e., a name passed to the
    external `betse` command identifying an action to be performed).

    This metadata encapsulates all human-readable help strings for this
    subcommand as well as additional options and arguments accepted by this
    subcommand.

    Attributes
    ----------
    name : str
        Machine-readable name of this subcommand (e.g., `plot`).
    synopsis : str
        Human-readable synopsis of this subcommand.
    description : str
        Human-readable description of this subcommand.
    is_passed_yaml : bool
        `True` if this subcommand only accepts a configuration filename _or_
        `False` if this subcommand accepts no passed arguments.
    '''

    @type_check
    def __init__(
        self,
        name: str,
        synopsis: str,
        description: str,
        is_passed_yaml: bool,
    ) -> None:
        '''
        Describe this CLI subcommand.

        Params
        ----------
        name : str
            Machine-readable name of this CLI subcommand (e.g., `plot`),
            typically only a single word.
        synopsis : str
            Human-readable synopsis of this CLI subcommand, typically only one
            to three lines of lowercase, unpunctuated text. All `{`- and `}`-
            delimited format substrings (e.g., `{program_name}`) supported by
            the `_format_help()` method will be globally replaced.
        description : str
            Human-readable description of this CLI subcommand, typically one
            to several paragraphs of grammatical sentences. All `{`- and `}`-
            delimited format substrings (e.g., `{program_name}`) supported by
            the `_format_help()` method will be globally replaced.
        is_passed_yaml : bool
            `True` if this subcommand only accepts a configuration filename _or_
            `False` if this subcommand accepts no passed arguments.
        '''

        # Classify these parameters, expanding all default keywords in the
        # human-readable parameters.
        self.name = name
        self.synopsis = expand_help(synopsis)
        self.description = expand_help(description)
        self.is_passed_yaml = is_passed_yaml


    def get_arg_parser_kwargs(self) -> dict:
        '''
        Keyword arguments suitable for initializing the argument subparser
        parsing this subcommand.
        '''

        return {
            'name': self.name,
            'help': self.synopsis,
            'description': self.description,
        }

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

# ....................{ SUBCOMMANDS                        }....................
SUBCOMMANDS = (
    CLISubcommand(
        name='config',
        synopsis='create a default config file for {program_name} simulations',
        description='''
Write a default tissue simulation configuration to the passed output file. While
not strictly necessary, this file should have filetype ".yaml" . If this file
already exists, an error will be printed.

You may edit this file at any time. By default, this file instructs
{program_name} to save simulation results (e.g., plots) to the directory
containing this file.
''',
        is_passed_yaml=True,
    ),


    CLISubcommand(
        name='seed',
        synopsis='seed a new cell cluster for a config file',
        description='''
Create the cell cluster defined by the passed configuration file. The results
will be saved to output files defined by this configuration.
''',
        is_passed_yaml=True,
    ),


    CLISubcommand(
        name='init',
        synopsis='initialize a seeded cell cluster for a config file',
        description='''
Initialize (i.e., calculate steady-state concentrations for) the previously
created cell cluster defined by the passed configuration file. Initialization
results will be saved to output files defined by this configuration, while the
previously created cell cluster will be loaded from input files defined by this
configuration.
''',
        is_passed_yaml=True,
    ),


    CLISubcommand(
        name='sim',
        synopsis='simulate an initialized cell cluster for a config file',
        description='''
Simulate the previously initialized cell cluster defined by the passed
configuration file. Simulation results will be saved to output files defined by
this configuration, while the previously initialized cell cluster will be loaded
from input files defined by this configuration.
''',
        is_passed_yaml=True,
    ),


    CLISubcommand(
        name='sim-brn',
        synopsis='simulate a biochemical reaction network for a config file',
        description='''
Simulate a biochemical reaction network (BRN) for the previously initialized
cell cluster defined by the passed configuration file, whose "metabolism config"
option specifies the path of the configuration file defining this network.  All
other simulation features and options will be ignored.

Simulation results will be saved to output files defined by this configuration,
while the previously initialized cell cluster will be loaded from input files
defined by this configuration.
''',
        is_passed_yaml=True,
    ),


    CLISubcommand(
        name='sim-grn',
        synopsis='simulate a gene regulatory network for a config file',
        description='''
Simulate a gene regulatory network (GRN) for the previously initialized cell
cluster defined by the passed configuration file, whose "gene regulatory network
config" option specifies the path of the configuration file defining this
network. All other simulation features and options will be ignored.

Simulation results will be saved to output files defined by this configuration,
while the previously initialized cell cluster will be loaded from input files
defined by this configuration.
''',
        is_passed_yaml=True,
    ),


    CLISubcommand(
        name='plot',
        synopsis='plot a seeded, initialized, or simulated simulation',
        description='''
Run the passed plotting subcommand. For example, to plot the previous
simulation defined by a configuration file "my_sim.yaml" in the current
directory:

;    betse plot sim my_sim.yaml
''',
        is_passed_yaml=False,
    ),


    CLISubcommand(
        name='info',
        synopsis='show information about {program_name} and the current system',
        description='''
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
''',
        is_passed_yaml=False,
    ),


    CLISubcommand(
        name='try',
        synopsis='create, init, simulate, and plot a sample simulation',
        description='''
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
''',
        is_passed_yaml=False,
    ),
)
'''
Tuple of `CLISubcommand` instances describing top-level subcommands.

**Order is significant.** The order in which these instances are listed defines
the order in which the `betse --help` command synopsizes these subcommands.
Moreover, top-level subcommands _not_ listed here will _not_ be parsed by
argument subparsers and hence will be effectively ignored.
'''

# ....................{ SUBCOMMANDS                        }....................
SUBCOMMANDS_PLOT = (
    CLISubcommand(
        name='seed',
        synopsis='plot a seeded cell cluster for a config file',
        description='''
Plot the previously seeded cell cluster defined by the passed configuration
file. Plot results will be saved to output files defined by this configuration,
while the previously seeded cell cluster will be loaded from input files
defined by this configuration.
''',
        is_passed_yaml=True,
    ),


    CLISubcommand(
        name='init',
        synopsis='plot an initialized cell cluster for a config file',
        description='''
Plot the previously initialized cell cluster defined by the passed configuration
file. Plot results will be saved to output files defined by this configuration,
while the previously initialized cell cluster will be loaded from input files
defined by this configuration.
''',
        is_passed_yaml=True,
    ),


    CLISubcommand(
        name='sim',
        synopsis='plot a simulated cell cluster for a config file',
        description='''
Plot the previously simulated cell cluster defined by the passed configuration
file. Plot results will be saved to output files defined by this configuration,
while the previously simulated cell cluster will be loaded from input files
defined by this configuration.
''',
        is_passed_yaml=True,
    ),


    CLISubcommand(
        name='sim-brn',
        synopsis='plot a simulated biochemical reaction network for a config file',
        description='''
Plot the previously simulated biochemical reaction network (BRN) defined by the
passed configuration file. Plot results will be saved to output files defined by
this configuration, while the previously simulated cell cluster will be loaded
from input files defined by this configuration.
''',
        is_passed_yaml=True,
    ),


    CLISubcommand(
        name='sim-grn',
        synopsis='plot a simulated gene regulatory network for a config file',
        description='''
Plot the previously simulated gene regulatory network (BRN) defined by the
passed configuration file. Plot results will be saved to output files defined by
this configuration, while the previously simulated cell cluster will be loaded
from input files defined by this configuration.
''',
        is_passed_yaml=True,
    ),
)
'''
Tuple of `CLISubcommand` instances describing subcommands of the `plot`
subcommand.

See Also
----------
SUBCOMMANDS
    Further details on tuple structure and significance.
'''
