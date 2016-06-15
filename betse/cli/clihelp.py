#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Help strings printed by `betse`'s command line interface (CLI).
'''

#FIXME: Rename this submodule to "clihelp", as this module's name conflicts with
#that of the help() builtin.
#FIXME: Shift all "subcommand"-specific classes and globals into a new
#"betse.cli.subcommand" module. On doing so:
#
#* Rename the "CLISubcommand" class to merely "CLISubcommand".

# ....................{ IMPORTS                            }....................
from betse import metadata
from betse.util.command import commands
from betse.util.type import strs, types

# ....................{ FUNCTIONS                          }....................
#FIXME: Rename to merely sanitize_name() after shifting into a new
#"clisubcommand" module.
def sanitize_subcommand_name(subcommand_name: str) -> str:
    '''
    Sanitize the passed subcommand name (e.g., from `sim-gnr` to `sim_gnr`).

    Specifically, this utility function:

    * Replaces all hyphens in this name with underscores, as method names
        are generated from subcommand names but cannot contain hyphens.
    '''
    assert types.is_str_nonempty(subcommand_name), (
        types.assert_not_str_nonempty(subcommand_name, 'Subcommand name'))

    return subcommand_name.replace('-', '_')


def expand(text: str, **kwargs) -> str:
    '''
    Interpolate the passed keyword arguments into the passed help string
    template, stripping all prefixing and suffixing whitespace from this
    template.

    For convenience, the following default keyword arguments are
    unconditionally interpolated into this template:

    * `{script_basename}`, expanding to the basename of the current script
        (e.g., `betse`).
    * `{program_name}`, expanding to this script's human-readable name
        (e.g., `BETSE`).
    '''
    assert types.is_str(text), types.assert_not_str(text)

    return strs.remove_presuffix_whitespace(text.format(
        program_name=metadata.NAME,
        script_basename=commands.get_current_basename(),
        **kwargs
    ))

# ....................{ CLASSES                            }....................
class CLISubcommand(object):
    '''
    Collection of all human-readable help strings for one CLI subcommand.

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
        assert types.is_str_nonempty(name), (
            types.assert_not_str_nonempty(name, 'Subcommand name'))
        assert types.is_str_nonempty(synopsis), (
            types.assert_not_str_nonempty(synopsis, 'Subcommand synopsis'))
        assert types.is_str_nonempty(description), (
            types.assert_not_str_nonempty(description, 'Subcommand description'))
        assert types.is_bool(is_passed_yaml), (
            types.assert_not_bool(is_passed_yaml))

        # Classify these parameters, expanding all default keywords in the
        # human-readable parameters.
        self.name = name
        self.synopsis = expand(synopsis)
        self.description = expand(description)
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

# ....................{ TEMPLATES ~ OPTIONS                }....................
#FIXME: Refactor these string globals into a single dictionary mapping from
#option name to help string. The current approach is *MUCH* too heavyweight.

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


    #FIXME: Define help strings, please.
    CLISubcommand(
        name='sim-brn',
        synopsis='simulate a biochemical reaction network for a config file',
        description='''
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                  UNDEFINED                                 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
''',
        is_passed_yaml=True,
    ),


    #FIXME: Define help strings, please.
    CLISubcommand(
        name='sim-grn',
        synopsis='simulate a gene regulatory network for a config file',
        description='''
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                  UNDEFINED                                 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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


    #FIXME: Define help strings, please.
    CLISubcommand(
        name='sim-brn',
        synopsis='plot a simulated biochemical reaction network for a config file',
        description='''
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                  UNDEFINED                                 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
''',
        is_passed_yaml=True,
    ),


    #FIXME: Define help strings, please.
    CLISubcommand(
        name='sim-grn',
        synopsis='plot a simulated gene regulatory network for a config file',
        description='''
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                  UNDEFINED                                 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
