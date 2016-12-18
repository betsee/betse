#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Metadata describing options accepted by BETSE's command line interface (CLI).
'''

# ....................{ IMPORTS                            }....................
from abc import ABCMeta
from argparse import ArgumentParser
from betse.exceptions import BetseCLIArgException
from betse.util.py import identifiers
from betse.util.type import strs
from betse.util.type.types import (
    type_check, EnumType, EnumMemberType, MappingOrNoneTypes, StrOrNoneTypes)

# ....................{ SUPERCLASSES                       }....................
class CLIOptionABC(object, metaclass=ABCMeta):
    '''
    **CLI option** (i.e., `-` and/or `--`-prefixed option passed to the external
    `betse` command identifying an optional configuration setting).

    This class encapsulates all metadata pertaining to this option, including:

    * Human-readable help strings describing this option.
    * Logic required to convert this option from a low-level command-line string
      externally passed by end users into a high-level Python object internally
      accessed throughout this codebase.

    Attributes
    ----------
    _add_argument_args : list
        List of all positional arguments to be passed to the
        :meth:`ArgumentParser.add_argument` method, specifying the short and/or
        long variants of this option.
    _add_argument_kwargs : dict
        Dictionary of all keyword arguments to be passed to the
        :meth:`ArgumentParser.add_argument` method, converting this option from
        a command-line string into a Python object.
    _identifier : str
        Machine-readable string suitable for use as a **Python identifier**
        (i.e., class, function, module, or variable name), derived from this
        option's long and short variants.
    _name : StrOrNoneTypes
        Long variant of this CLI option if any _or_ `None` otherwise.
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(
        self,
        char: StrOrNoneTypes,
        name: StrOrNoneTypes,
        synopsis: str,
        synopsis_kwargs: MappingOrNoneTypes,
        add_argument_kwargs: MappingOrNoneTypes,
    ) -> None:
        '''
        Define this CLI option.

        Parameters
        ----------
        char : optional[str]
            `-`-prefixed machine-readable name of the short variant of this CLI
            option (e.g., `-v`) if defined _or_ `None` otherwise. If non-`None`,
            this is typically only a single character. Defaults to `None`.
        name : optional[str]
            `--`-prefixed and `-`-delimited machine-readable name of the long
            variant of this CLI option (e.g., `--matplotlib-backend`) if defined
            _or_ `None` otherwise. If non-`None`, this is typically either one
            word _or_ two words delimited by `-`. Defaults to `None`.
        synopsis : str
            Human-readable synopsis of this CLI option, typically only one to
            three lines of lowercase, unpunctuated text. All `{`- and `}`-
            delimited format substrings (e.g., `{program_name}`) supported by
            the :meth:`cliutil.expand_help` method will be globally replaced.
        synopsis_kwargs : optional[MappingType]
            Dictionary of all keyword arguments to be interpolated into this
            synopsis by the :meth:`betse.cli.cliabc.expand_Help` function if any
            _or_ `None` otherwise. Defaults to `None`.
        add_argument_kwargs : optional[MappingType]
            Dictionary of all keyword arguments to be passed to the
            :meth:`ArgumentParser.add_argument` method if any _or_ `None`
            otherwise. If non-`None`, these arguments typically define the
            onversion of this option from a command-line string into a Python
            object. Defaults to `None`.
        '''

        #FIXME: To avoid circularity issues, shift expand_help() into a new
        #"cliutil" submodule of this subpackage.

        # Avoid circular import dependencies.
        from betse.cli.cliabc import expand_help

        # If neither a short nor long variant of this option is passed, fail.
        if name is None and char is None:
            raise BetseCLIArgException(
                'Short and long option variants both undefined.')

        # Default all unpassed parameters *BEFORE* classifying these parameters.
        if add_argument_kwargs is None:
            add_argument_kwargs = {}
        if synopsis_kwargs is None:
            synopsis_kwargs = {}

        # Positional and keyword arguments to be passed to add_argument().
        self._add_argument_args = []
        self._add_argument_kwargs = add_argument_kwargs

        # Synopsize this option, interpolating these keyword arguments if any
        # into this synopsis.
        self._add_argument_kwargs['help'] = expand_help(
            synopsis, **synopsis_kwargs)

        # Python identifier derived from this option's long and short variants.
        self._identifier = None

        # If a short variant is passed, validate this variant. To preferentially
        # assign this option's identifier to the long variant validated below,
        # do so *BEFORE* validating the long variant below.
        if char is not None:
            # Remove the "-" prefixing this variant if present.
            char = strs.remove_prefix_if_found(text=name, prefix='-')

            # If the remainder is *NOT* a single character, fail.
            if len(char) != 1:
                raise BetseCLIArgException(
                    'Short option variant "{}" not a character.'.format(char))

            # Default this option's identifier to this short variant. Since this
            # variant is guaranteed to be only a single character, no string
            # munging is needed.
            self._identifier = char

            # (Re)prefix this variant by "-".
            char = '-' + char

            # Add this variant to these positional arguments.
            self._add_argument_args(char)

        # If a long variant is passed...
        if name is not None:
            # Remove the "--" prefixing this variant if present.
            name = strs.remove_prefix_if_found(text=name, prefix='--')

            # If the remainder is less than two characters long, fail.
            if len(name) < 2:
                raise BetseCLIArgException(
                    'Long option variant "{}" empty or '
                    'only one character.'.format(name))

            # Replace this option's identifier by this long variant, sanitizing.
            # all hyphens with underscores to generate a conformant identifier.
            self._identifier = identifiers.sanitize(name)

            # (Re)prefix this variant by "--".
            name = '--' + name

            # Add this variant to these positional arguments.
            self._add_argument_args(name)

    # ..................{ ADDERS                             }..................
    @type_check
    def add(self, arg_parser: ArgumentParser) -> None:
        '''
        Add a new argument parsing this option to the passed argument parser.

        Parameters
        ----------
        arg_parser : _SubParsersAction
            Argument parsers to add an argument parsing this option to.
        '''

        # Add an argument parsing this option to this parser.
        arg_parser.add_argument(
            *self._add_argument_args, **self._add_argument_kwargs)

# ....................{ SUBCLASSES                         }....................
class CLIOptionVersion(CLIOptionABC):
    '''
    CLI option reporting the version of this application.
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(
        self,
        *args,
        version: str,
        **kwargs
    ) -> None:
        '''
        Define this CLI option.

        Parameters
        ----------
        version : str
            Version specifier to be printed when this option is passed.

        All remaining parameters are passed as is to the superclass method.
        '''

        # Initialize our superclass with the passed arguments.
        super().__init__(*args, **kwargs)

        # Print this version specifier when this option is passed.
        self._add_argument_kwargs.update({
            'action': 'version',
            'version': version,
        })


class CLIOptionBool(CLIOptionABC):
    '''
    CLI option persisting a mandatory string argument in boolean format (i.e.,
    either `true` or `false`) into an instance variable of type `bool`.
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(self, *args, **kwargs) -> None:

        # Initialize our superclass with the passed arguments.
        super().__init__(*args, **kwargs)

        # Validate and convert this option's string argument.
        self._add_argument_kwargs.update({
            # Persist this option's string argument to a boolean variable.
            'action': 'store_true',

            # Prefix this variable by "is_", conforming to codebase nomenclature
            # (e.g., from "verbose" to "is_verbose").
            'dest': 'is_' + self._identifier,
        })


class CLIOptionEnum(CLIOptionABC):
    '''
    CLI option persisting a mandatory string argument matching one of several
    alternatives into an instance variable whose value is the name of the
    corresponding member of an enumeration type.

    Caveats
    ----------
    Unfortunately, the :class:`ArgumentParser` API significantly complicates
    what should otherwise be a simplistic conversion from string arguments to
    enumeration members. Ideally, this class would instead persist into an
    instance variable of high-level type :class:`Enum` rather than of low-level
    type :class:`str`. This is currently infeasible.

    Why? Entirely because of the overly simplistic `choices` parameter accepted
    by the :meth:`ArgumentParser.add_argument` method, which validates this
    option's argument _after_ converting that argument from a string into the
    desired with the type or callable supplied by the `type` parameter. Whereas
    the names of enumeration members are uppercase throughout this codebase
    (e.g., `ProfileType['SIZE']`), their corresponding arguments are lowercase
    (e.g., `--profile-type=size`). Since `choices` provides no means of mapping
    between the two, callers should do so manually instead.
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(
        self,
        *args,
        enum_type: EnumType,
        default: EnumMemberType,
        synopsis_kwargs: MappingOrNoneTypes,
        **kwargs
    ) -> None:
        '''
        Define this CLI option.

        Parameters
        ----------
        enum_type : EnumType
            Enumeration constraining this option's string argument.
        default : EnumMemberType
            Member of this enumeration to default to if this option is unpassed.

        All remaining parameters are passed as is to the superclass method.
        '''

        # Avoid circular import dependencies.
        from betse.util.type import enums

        # Default all unpassed parameters *BEFORE* classifying these parameters.
        if synopsis_kwargs is None:
            synopsis_kwargs = {}

        # Lowercased name of the default member of this enumeration.
        default_name = default.name.lower()

        # Interpolate this name into this option's synopsis.
        synopsis_kwargs['default'] = default_name

        # Initialize our superclass with the passed arguments.
        super().__init__(*args, synopsis_kwargs=synopsis_kwargs, **kwargs)

        # Validate and convert this option's string argument.
        self._add_argument_kwargs.update({
            # Persist this option's argument to a variable with this name.
            'action': 'store',
            'dest': self._identifier,

            # Set of the lowercased names of all members of this enumeration.
            'choices': enums.get_names_lowercase(enum_type),

            # Lowercased name of the default member of this enumeration.
            'default': default_name,
        })

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
