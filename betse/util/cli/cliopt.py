#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Classes (both abstract and concrete) for defining ``-``- and ``--``-prefixed
command line interface (CLI) options.
'''

# ....................{ IMPORTS                            }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: To raise human-readable exceptions on application startup, the
# top-level of this module may import *ONLY* from submodules guaranteed to:
# * Exist, including standard Python and application modules.
# * Never raise exceptions on importation (e.g., due to module-level logic).
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

from abc import ABCMeta
from betse.exceptions import BetseCLIArgException
from betse.util.py import pyident
from betse.util.type.text import strs
from betse.util.type.types import (
    type_check,
    ArgParserType,
    EnumType,
    EnumMemberType,
    MappingOrNoneTypes,
    StrOrNoneTypes,
)

# ....................{ EXPANDERS                          }....................
#FIXME: Replace with usage of the CLIABC.expand_help() method; then excise this.
#To do so, we'll probably need to:
#
#* Define a new "CLIOptioner" container class analogous to the
#  "clicmd.CLISubcommander" container class.
#* Define a CLIOptioner.add() method.
#* Refactor the CLIOptionABC.__init__() method to cease calling expand_help().
#* Refactor all CLIOptionABC.add() methods to:
#  * Accept an initial "cli" parameter.
#  * Call the cli.expand_help() method to expand help strings as needed.
#
#Trivial, of course, given that we've already done so for subcommands. *sigh*

@type_check
def expand_help(text: str, **kwargs) -> str:
    '''
    Interpolate the passed keyword arguments into the passed help string
    template, stripping all prefixing and suffixing whitespace from this
    template.

    For convenience, the following default keyword arguments are unconditionally
    interpolated into this template:

    * ``{script_basename}``, expanding to the basename of the current script
        (e.g., ``betse``).
    * ``{program_name}``, expanding to this script's human-readable name
        (e.g., ``BETSE``).
    '''

    from betse import metadata
    from betse.util.path.command import cmds

    return strs.remove_whitespace_presuffix(text.format(
        program_name=metadata.NAME,
        script_basename=cmds.get_current_basename(),
        **kwargs
    ))

# ....................{ SUPERCLASSES                       }....................
class CLIOptionABC(object, metaclass=ABCMeta):
    '''
    **CLI option** (i.e., ``-`` and/or ``--``-prefixed option passed to the
    external ``betse`` command identifying an optional configuration setting).

    This class encapsulates all metadata pertaining to this option, including:

    * Human-readable help strings describing this option.
    * Logic required to convert this option from a low-level command-line string
      externally passed by end users into a high-level Python object internally
      accessed throughout this codebase.

    Attributes
    ----------
    _add_argument_args : list
        List of all positional arguments to be passed to the
        :meth:`ArgParserType.add_argument` method, specifying the short and/or
        long variants of this option.
    _add_argument_kwargs : dict
        Dictionary of all keyword arguments to be passed to the
        :meth:`ArgParserType.add_argument` method, converting this option from
        a command-line string into a Python object.
    _var_name : str
        Machine-readable string suitable for use as a **Python identifier**
        (i.e., class, function, module, or variable name), derived from this
        option's long and short variants.
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(
        self,

        # Mandatory parameter.
        synopsis: str,

        # Optional parameters.
        short_name: StrOrNoneTypes = None,
        long_name: StrOrNoneTypes = None,
        var_name: StrOrNoneTypes = None,
        synopsis_kwargs: MappingOrNoneTypes = None,
        add_argument_kwargs: MappingOrNoneTypes = None,
    ) -> None:
        '''
        Define this CLI option.

        Parameters
        ----------
        synopsis : str
            Human-readable synopsis of this CLI option, typically only one to
            three lines of lowercase, unpunctuated text. All ``{``- and ``}``-
            delimited format substrings (e.g., ``{program_name}``) supported by
            the :meth:`cliutil.expand_help` function will be globally replaced.
        short_name : optional[str]
            ``-``-prefixed machine-readable name of the short variant of this
            CLI option (e.g., ``-v``) if defined *or* ``None`` otherwise. If
            non-``None``, this is ideally only a single character. Defaults to
            ``None``.
        long_name : optional[str]
            ``--``-prefixed and ``-``-delimited machine-readable name of the
            long variant of this CLI option (e.g., ``--matplotlib-backend``) if
            defined *or* ``None`` otherwise. If non-``None``, this is typically
            either one word *or* two words delimited by ``-``. Defaults to
            ``None``.
        var_name : optional[str]
            Name of the instance variable to which the converted value of this
            option's optional argument is to be persisted. Defaults to ``None``,
            in which case this name defaults to a Python identifier derived from
            this option's long and short variants.
        synopsis_kwargs : optional[MappingType]
            Dictionary of all keyword arguments to be interpolated into this
            synopsis by the :meth:`cliutil.expand_help` function if any
            *or* ``None`` otherwise. Defaults to ``None``, in which case this
            dictionary is empty.
        add_argument_kwargs : optional[MappingType]
            Dictionary of all keyword arguments to be passed to the
            :meth:`ArgParserType.add_argument` method if any *or* ``None``
            otherwise. If non-``None``, these arguments typically define the
            onversion of this option from a command-line string into a Python
            object. Defaults to ``None``.
        '''

        # If neither a short nor long variant of this option is passed, fail.
        if short_name is None and long_name is None:
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
        var_name_default = None

        # If a short variant is passed, validate this variant. To preferentially
        # assign this option's identifier to the long variant validated below,
        # do so *BEFORE* validating the long variant below.
        if short_name is not None:
            # Remove the "-" prefixing this variant if present.
            short_name = strs.remove_prefix_if_found(
                text=short_name, prefix='-')

            # If the remainder is *NOT* a single character, fail.
            if len(short_name) != 1:
                raise BetseCLIArgException(
                    'Short option variant "{}" not a character.'.format(
                        short_name))

            # Default this option's identifier to this short variant. Since this
            # variant is guaranteed to be only a single character, no string
            # munging is needed.
            var_name_default = short_name

            # (Re)prefix this variant by "-".
            short_name = '-' + short_name

            # Add this variant to these positional arguments.
            self._add_argument_args.append(short_name)

        # If a long variant is passed...
        if long_name is not None:
            # Remove the "--" prefixing this variant if present.
            long_name = strs.remove_prefix_if_found(text=long_name, prefix='--')

            # If the remainder is less than two characters long, fail.
            if len(long_name) < 2:
                raise BetseCLIArgException(
                    'Long option variant "{}" empty or '
                    'only one character.'.format(long_name))

            # Replace this option's identifier by this long variant, sanitizing.
            # all hyphens with underscores to generate a conformant identifier.
            var_name_default = pyident.sanitize_snakecase(long_name)

            # (Re)prefix this variant by "--".
            long_name = '--' + long_name

            # Add this variant to these positional arguments.
            self._add_argument_args.append(long_name)

        # Name of the variable to which this option is persisted, set to the
        # passed name if any or defaulting to the name set above otherwise.
        self._var_name = var_name if var_name is not None else var_name_default

        # If this vraiable name is invalid, raise an exception.
        pyident.die_unless_var_name(self._var_name)

    # ..................{ ADDERS                             }..................
    @type_check
    def add(self, arg_parser: ArgParserType) -> None:
        '''
        Add an argument parsing this option to the passed argument parser.

        Parameters
        ----------
        arg_parser : ArgParserType
            Argument parser to add this argument to.
        '''
        # print('\nAdding option... args: {}; kwargs: {}'.format(self._add_argument_args, self._add_argument_kwargs))

        # Add an argument parsing this option to this parser.
        arg_parser.add_argument(
            *self._add_argument_args, **self._add_argument_kwargs)

# ....................{ CLASSES ~ argless                  }....................
class CLIOptionBoolTrue(CLIOptionABC):
    '''
    CLI option accepting *no* arguments setting an instance variable to ``True``
    when passed and ``False`` when unpassed.
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
            'dest': 'is_' + self._var_name,
        })


class CLIOptionVersion(CLIOptionABC):
    '''
    CLI option accepting *no* arguments printing a version specifier and halting
    the current process when passed.
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(self, *args, version: str, **kwargs) -> None:
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

# ....................{ CLASSES ~ arg                      }....................
class CLIOptionArgABC(CLIOptionABC):
    '''
    Abstract base class of all CLI option subclasses accepting an optional
    argument.persisted to an instance variable whose name derives from this
    option's long and short variants.
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(
        self,
        *args,

        # Mandatory parameters.
        default_value: object,

        # Optional parameters.
        synopsis_kwargs: MappingOrNoneTypes = None,
        **kwargs
    ) -> None:
        '''
        Define this CLI option.

        Parameters
        ----------
        default_value : object
            Object to default this option's optional argument to if unpassed.
        synopsis_kwargs : optional[MappingType]
            Dictionary of all keyword arguments to be interpolated into this
            synopsis by the :meth:`betse.cli.cliutil.expand_help` function if
            any *or* ``None`` otherwise. Defaults to ``None``, in which case
            this dictionary is empty. This method adds the following keys to
            this dictionary:
            * ``default``, whose value is the passed ``default_value``
              parameter. Each substring ``{default}`` in the passed
              ``synopsis`` parameter is globally replaced with this value.

        All remaining parameters are passed as is to the superclass method.
        '''

        # Default all unpassed parameters *BEFORE* classifying these parameters.
        if synopsis_kwargs is None:
            synopsis_kwargs = {}

        # Interpolate this name into this option's synopsis.
        synopsis_kwargs['default'] = default_value

        # Initialize our superclass with the passed arguments.
        super().__init__(*args, synopsis_kwargs=synopsis_kwargs, **kwargs)

        # Validate and convert this option's string argument.
        self._add_argument_kwargs.update({
            # Persist this option's argument to a variable with this name.
            'action': 'store',
            'dest': self._var_name,

            # Lowercased name of the default member of this enumeration.
            'default': default_value,
        })


class CLIOptionArgEnum(CLIOptionArgABC):
    '''
    CLI option persisting an optional string argument matching one of several
    alternatives into an instance variable whose value is the name of the
    corresponding member of an enumeration type.

    Caveats
    ----------
    Unfortunately, the :class:`ArgParserType` API significantly complicates
    what should otherwise be a simplistic conversion from string arguments to
    enumeration members. Ideally, this class would instead persist into an
    instance variable of high-level type :class:`Enum` rather than of low-level
    type :class:`str`. This is currently infeasible.

    Why? Entirely thanks to the overly simplistic ``choices`` parameter accepted
    by the :meth:`ArgParserType.add_argument` method, which validates this
    option's argument *after* converting that argument from a string into the
    desired with the type or callable supplied by the `type` parameter. Whereas
    the names of enumeration members are uppercase throughout this codebase
    (e.g., ``ProfileType['SIZE']``), their corresponding arguments are lowercase
    (e.g., ``--profile-type=size``). Since `choices` provides no means of
    mapping between the two, callers should do so manually instead.
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(
        self,
        *args,
        enum_type: EnumType,
        enum_default: EnumMemberType,
        **kwargs
    ) -> None:
        '''
        Define this CLI option.

        Parameters
        ----------
        enum_type : EnumType
            Enumeration constraining this option's string argument.
        enum_default : EnumMemberType
            Member of this enumeration to default to if this option is unpassed.

        All remaining parameters are passed as is to the superclass method.
        '''

        # Avoid circular import dependencies.
        from betse.util.type import enums

        # If this default is *NOT* a member of this enumeration, fail.
        enums.die_unless_member(enum_type, enum_default)

        # Initialize our superclass with the passed arguments.
        super().__init__(
            *args,
            # Default this option's argument if unpassed to the lowercased name
            # of this default member of this enumeration.
            default_value=enum_default.name.lower(),
            **kwargs)

        # Validate and convert this option's string argument.
        self._add_argument_kwargs.update({
            # Set of the lowercased names of all members of this enumeration.
            'choices': enums.get_member_names_lowercase(enum_type),
        })


class CLIOptionArgStr(CLIOptionArgABC):
    '''
    CLI option persisting an optional string argument of arbitrary format into
    an instance variable of type :class:`str`.
    '''

    # The best things in life are always free.
    pass
