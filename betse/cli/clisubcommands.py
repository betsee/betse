#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Metadata describing subcommands accepted by BETSE's command line interface
(CLI).
'''

# ....................{ IMPORTS                            }....................
from abc import ABCMeta
from betse.cli import cliutil
from betse.util.type.types import (
    type_check, ArgParserType, ArgSubparsersType, IterableTypes, MappingType)

# ....................{ SUPERCLASSES                       }....................
class CLISubcommandABC(object, metaclass=ABCMeta):
    '''
    Abstract base class of all **CLI subcommand** (i.e., shell word passed to
    the external `betse` command signifying a high-level action to be
    performed) subclasses.

    This class encapsulates all metadata pertaining to this subcommand,
    including:

    * Human-readable help strings describing this subcommand.
    * All options and arguments accepted by this subcommand.

    Attributes
    ----------
    name : str
        Machine-readable name of this CLI subcommand (e.g., `plot`), typically
        only a single word.
    synopsis : str
        Human-readable synopsis of this CLI subcommand, typically only one to
        three lines of lowercase, unpunctuated text. All `{`- and `}`- delimited
        format substrings (e.g., `{program_name}`) supported by the
        :meth:`cliutil.expand_help` function will be globally replaced.
    description : str
        Human-readable description of this CLI subcommand, typically one to
        several paragraphs of grammatical sentences. All `{`- and `}`- delimited
        format substrings (e.g., `{program_name}`) supported by the
        :meth:`cliutil.expand_help` function will be globally replaced.
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(
        self,
        name: str,
        synopsis: str,
        description: str,
    ) -> None:
        '''
        Define this CLI subcommand.

        Parameters
        ----------
        See the class docstring for parameter documentation. All parameters
        accepted by this method are instance variables of the same name.
        '''

        # Classify these parameters, expanding all default keywords in the
        # human-readable parameters.
        self.name = name
        self.synopsis = cliutil.expand_help(synopsis)
        self.description = cliutil.expand_help(description)

    # ..................{ ADDERS                             }..................
    @type_check
    def add(
        self,
        arg_subparsers: ArgSubparsersType,
        arg_subparser_kwargs: MappingType,
    ) -> (
        ArgParserType):
        '''
        Create a new **argument subparser** (i.e., :mod:`argparse`-specific
        object parsing command-line arguments) parsing this subcommand, add this
        subparser to the passed collection of **argument subparsers** (i.e.,
        another :mod:`argparse`-specific object cotaining multiple subparsers),
        and return this subparser.

        This subparser is configured to:

        * If this subcommand accepts a configuration filename, require such an
          argument be passed.
        * Else, require no arguments be passed.

        Parameters
        ----------
        arg_subparsers : ArgSubparsersType
            Collection of sibling subcommand argument parsers to which the
            subcommand argument parser created by this method is added. This
            collection is owned either by:
            * A top-level subcommand (e.g., `plot`), in which case the
              subcommand created by this method is a child of that subcommand.
            * No subcommand, in which case the subcommand created by this method
              is a top-level subcommand.
        arg_subparser_kwargs : MappingType
            Dictionary of all keyword arguments to be passed to theh
            :meth:`ArgumentParser.init` method creating this argument subparser.

        Returns
        ----------
        ArgParserType
            Subcommand argument parser created by this method.
        '''

        # Initialize this parser with subcommand-specific keyword arguments.
        kwargs = {
            'name':        self.name,
            'help':        self.synopsis,
            'description': self.description,
        }
        kwargs.update(arg_subparser_kwargs)

        # Create and return this parser, added to this container of subparsers.
        return arg_subparsers.add_parser(**kwargs)

# ....................{ SUBCLASSES                         }....................
class CLISubcommandNoArg(CLISubcommandABC):
    '''
    CLI subcommand accepting _no_ passed arguments.
    '''

    pass


class CLISubcommandParent(CLISubcommandABC):
    '''
    CLI subcommand that is itself the parent of one or more CLI subcommands,
    accepting _only_ the name of a child subcommand as a passed argument.
    '''

    # We almost don't believe it either.
    pass


class CLISubcommandYAMLOnly(CLISubcommandABC):
    '''
    CLI subcommand accepting _only_ a configuration filename as a passed
    argument.
    '''

    # ..................{ ADDERS                             }..................
    def add(self, *args, **kwargs) -> ArgParserType:

        # Subcommand argument subparser added by our superclass.
        arg_subparser = super().add(*args, **kwargs)

        # Configure this subparser to require a configuration file argument.
        arg_subparser.add_argument(
            'config_filename',
            metavar='CONFIG_FILE',
            help='simulation configuration file',
        )

        # Return this subparser.
        return arg_subparser

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

# ....................{ ADDERS ~ top                       }....................
@type_check
def add_top(
    arg_subparsers: ArgSubparsersType,
    arg_subparser_kwargs: MappingType,
) -> dict:
    '''
    Add one new argument subparser parsing each top-level subcommand to the
    passed container of argument subparsers.

    Parameters
    ----------
    arg_subparsers : ArgSubparsersType
        Container of argument subparsers to add these subparsers to.
    arg_subparser_kwargs : MappingType
        Dictionary of all keyword arguments to be passed to the
        :meth:`CLISubcommandABC.add` method of each such subcommand.

    Returns
    ----------
    dict
        Dictionary mapping from the name of each such subcommand to the new
        argument subparser parsing this subcommand.
    '''

    # Tuple of "CLISubcommandABC" instances describing top-level subcommands.
    #
    # Order is significant, defining that the `betse --help` command synopsizes
    # these subcommands in. Subcommands *NOT* listed here are *NOT* be parsed by
    # argument subparsers and hence effectively ignored.
    SUBCOMMANDS_TOP = (
        CLISubcommandYAMLOnly(
            name='config',
            synopsis=(
                'create a default config file for {program_name} simulations'),
            description='''
Write a default tissue simulation configuration to the passed output file. While
not strictly necessary, this file should have filetype ".yaml" . If this file
already exists, an error will be printed.

You may edit this file at any time. By default, this file instructs
{program_name} to save simulation results (e.g., plots) to the directory
containing this file.
''',
        ),

        CLISubcommandYAMLOnly(
            name='seed',
            synopsis='seed a new cell cluster for a config file',
            description='''
Create the cell cluster defined by the passed configuration file. The results
will be saved to output files defined by this configuration.
''',
        ),

        CLISubcommandYAMLOnly(
            name='init',
            synopsis='initialize a seeded cell cluster for a config file',
            description='''
Initialize (i.e., calculate steady-state concentrations for) the previously
created cell cluster defined by the passed configuration file. Initialization
results will be saved to output files defined by this configuration, while the
previously created cell cluster will be loaded from input files defined by this
configuration.
''',
        ),

        CLISubcommandYAMLOnly(
            name='sim',
            synopsis='simulate an initialized cell cluster for a config file',
            description='''
Simulate the previously initialized cell cluster defined by the passed
configuration file. Simulation results will be saved to output files defined by
this configuration, while the previously initialized cell cluster will be loaded
from input files defined by this configuration.
''',
        ),

        CLISubcommandYAMLOnly(
            name='sim-brn',
            synopsis=(
                'simulate a biochemical reaction network for a config file'),
            description='''
Simulate a biochemical reaction network (BRN) for the previously initialized
cell cluster defined by the passed configuration file, whose "metabolism config"
option specifies the path of the configuration file defining this network.  All
other simulation features and options will be ignored.

Simulation results will be saved to output files defined by this configuration,
while the previously initialized cell cluster will be loaded from input files
defined by this configuration.
''',
        ),

        CLISubcommandYAMLOnly(
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
        ),

        CLISubcommandParent(
            name='plot',
            synopsis='plot a seeded, initialized, or simulated simulation',
            description='''
Run the passed plotting subcommand. For example, to plot the previous
simulation defined by a configuration file "my_sim.yaml" in the current
directory:

;    betse plot sim my_sim.yaml
''',
        ),

        CLISubcommandNoArg(
            name='repl',
            synopsis='drop into a REPL within an initialized BETSE environment',
            description='''
Initialize the BETSE environment and immediately open a REPL. This allows
interactive manipulation of the simulations and analyses.
''',
        ),

        CLISubcommandNoArg(
            name='info',
            synopsis=(
                'show information about {program_name} and the current system'),
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
        ),

        CLISubcommandNoArg(
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
        ),
    )

    # Defer to this low-level function to add subparsers for these subcommands.
    return _add_subcommands(
        subcommands=SUBCOMMANDS_TOP,
        arg_subparsers=arg_subparsers,
        arg_subparser_kwargs=arg_subparser_kwargs,
    )

# ....................{ ADDERS ~ plot                      }....................
@type_check
def add_plot(
    arg_subparsers: ArgSubparsersType,
    arg_subparser_kwargs: MappingType,
) -> dict:
    '''
    Add one new argument subparser parsing each subcommand of the top-level
    `plot` subcommand to the passed container of argument subparsers.

    Parameters
    ----------
    arg_subparsers : ArgSubparsersType
        Container of argument subparsers to add these subparsers to.
    arg_subparser_kwargs : MappingType
        Dictionary of all keyword arguments to be passed to the
        :meth:`CLISubcommandABC.add` method of each such subcommand.

    Returns
    ----------
    dict
        Dictionary mapping from the name of each such subcommand to the new
        argument subparser parsing that subcommand.
    '''

    # Tuple of "CLISubcommandABC" instances describing plot subcommands.
    #
    # Order is significant, defining that the `betse --help` command synopsizes
    # these subcommands in. Subcommands *NOT* listed here are *NOT* be parsed by
    # argument subparsers and hence effectively ignored.
    SUBCOMMANDS_PLOT = (
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
            name='sim-brn',
            synopsis=(
                'plot a simulated biochemical reaction network '
                'for a config file'),
            description='''
Plot the previously simulated biochemical reaction network (BRN) defined by the
passed configuration file. Plot results will be saved to output files defined by
this configuration, while the previously simulated cell cluster will be loaded
from input files defined by this configuration.
''',
        ),

        CLISubcommandYAMLOnly(
            name='sim-grn',
            synopsis=(
                'plot a simulated gene regulatory network for a config file'),
            description='''
Plot the previously simulated gene regulatory network (BRN) defined by the
passed configuration file. Plot results will be saved to output files defined by
this configuration, while the previously simulated cell cluster will be loaded
from input files defined by this configuration.
''',
        ),
    )

    # Defer to this low-level function to add subparsers for these subcommands.
    return _add_subcommands(
        subcommands=SUBCOMMANDS_PLOT,
        arg_subparsers=arg_subparsers,
        arg_subparser_kwargs=arg_subparser_kwargs,
    )

# ....................{ ADDERS ~ private                   }....................
@type_check
def _add_subcommands(
    subcommands: IterableTypes,
    arg_subparsers: ArgSubparsersType,
    arg_subparser_kwargs: MappingType,
) -> dict:
    '''
    Add one new argument subparser parsing each subcommand in the passed
    iterable of subcommands to the passed container of argument subparsers and
    return a dictionary mapping from the name of each such subcommand to the
    corresponding argument subparser.

    Parameters
    ----------
    subcommands : IterableTypes
        Iterable of all subcommands to add subparsers for.
    arg_subparsers : ArgSubparsersType
        Container of argument subparsers to add these subparsers to.
    arg_subparser_kwargs : MappingType
        Dictionary of all keyword arguments to be passed to the
        :meth:`CLISubcommandABC.add` method of each such subcommand.

    Returns
    ----------
    dict
        Dictionary mapping from the name of each such subcommand to the new
        argument subparser parsing that subcommand.
    '''

    # Dictionary to be returned.
    subcommand_name_to_subparser = {}

    # For each top-level subcommand, add an argument parser parsing this
    # subcommand to this container of argument subparsers.
    for subcommand in subcommands:
        subcommand_name_to_subparser[subcommand.name] = subcommand.add(
            arg_subparsers=arg_subparsers,
            arg_subparser_kwargs=arg_subparser_kwargs)

    # Return this dictionary.
    return subcommand_name_to_subparser
