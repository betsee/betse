#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Metadata describing options accepted by BETSE's command line interface (CLI).
'''

# ....................{ IMPORTS                            }....................
from argparse import ArgumentParser
from betse.util.type.types import type_check, StrOrNoneTypes

# ....................{ SUPERCLASSES                       }....................
class CLIOption(object):
    '''
    Metadata encapsulating a **CLI option** (i.e., `-` and/or `--`-prefixed
    option passed to the external `betse` command identifying an optional
    configuration setting).

    This metadata encapsulates all human-readable help strings for this option
    as well as the conversion of this option from a low-level command-line
    string externally passed by users into a high-level Python object internally
    accessed by the codebase.

    Attributes
    ----------
    name : str
        `--`-prefixed and `-`-delimited machine-readable name of the long
        variant of this CLI option (e.g., `--matplotlib-backend`), typically
        either a single word or two words.
    char : str
        `-`-prefixed machine-readable name of the short variant of this CLI
        option (e.g., `-v`), typically only a single character.
    synopsis : str
        Human-readable synopsis of this CLI option, typically only one to three
        lines of lowercase, unpunctuated text. All `{`- and `}`- delimited
        format substrings (e.g., `{program_name}`) supported by the
        :meth:`cliutil.expand_help` method will be globally replaced.
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(
        self,
        name: StrOrNoneTypes,
        char: StrOrNoneTypes,
        synopsis: str,
    ) -> None:
        '''
        Define this CLI option.

        Parameters
        ----------
        See the class docstring for parameter documentation. All parameters
        accepted by this method are instance variables of the same name.
        '''

        #FIXME: To avoid circularity issues, shift expand_help() into a new
        #"cliutil" submodule of this subpackage.

        # Avoid circular import dependencies.
        from betse.cli.cliabc import expand_help

        # Classify these parameters, expanding all default keywords in the
        # human-readable parameters.
        self.name = name
        self.char = char
        self.synopsis = expand_help(synopsis)

        #FIXME: Raise a human-readable exception if "name" and "char" are both
        #None. (Urgh.)

    # ..................{ ADDERS                             }..................
    @type_check
    def add(self, arg_parser: ArgumentParser, **kwargs) -> None:
        '''
        Add a new argument parsing this option to the passed argument parser.

        Parameters
        ----------
        arg_parser : _SubParsersAction
            Argument parsers to add an argument parsing this option to.

        All remaining keyword arguments are passed as is to the
        :meth:`ArgumentParser.add_argument` method of this parser.
        '''

        # Positional arguments to be passed to the add_argument() method below.
        args = []

        # Add the short and long variants of this option to these arguments.
        if self.char is not None:
            args.append(self.char)
        if self.name is not None:
            args.append(self.name)

        # Add an argument parsing this option to this parser.
        arg_parser.add_argument(*args, **kwargs)

# ....................{ SUBCLASSES                         }....................

# ....................{ OPTIONS                            }....................
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
;* "none", logging to stdout and stderr only
;* "file", logging to stdout, stderr, and "--log-file"
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
;* "call", profiling callables (functions, methods)
;* "line", profiling code lines (requires "pprofile")
;* "size", profiling object sizes (requires "pympler")
'''
'''
Help string template synopsizing the `--profile-type` option.
'''

OPTION_PROFILE_FILE = '''
file to profile to if "--profile-type" is not "none" (defaults to "{default}")
'''
'''
Help string template synopsizing the `--profile-file` option.
'''
