#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Subcommands accepted by this application's command line interface (CLI).
'''

# ....................{ IMPORTS                            }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: To raise human-readable exceptions on application startup, the
# top-level of this module may import *ONLY* from submodules guaranteed to:
# * Exist, including standard Python and application modules.
# * Never raise exceptions on importation (e.g., due to module-level logic).
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

from abc import ABCMeta
from betse.util.type.types import (
    type_check,
    ArgParserType,
    ArgSubparsersType,
    CallableTypes,
    MappingType,
    SequenceTypes,
)

# ....................{ CLASSES ~ container                }....................
class CLISubcommands(object):
    '''
    Container of all CLI subcommands accepted by either a top-level CLI command
    (e.g., ``betse``) or a subcommand of such command (e.g., ``betse plot``).

    This class encapsulates all metadata pertaining to this container,
    including:

    * Human-readable help strings describing this container.
    * All options and arguments accepted by this container.

    Attributes
    ----------
    _subcommand_var_name : str
        Name of the instance variable to which the :class:`ArgumentParser`
        parsing these subcommands from command-line options and arguments
        passed by the external user stores either:
        * The name of the currently passed subcommand for this container if
            the user passed such a subcommand.
        * ``None`` otherwise.
    _subcommands : SequenceTypes
        Sequence of all **subcommands** (i.e., :class:`CLISubcommandABC`
        instances) to add to this container.
    _help_description : str
        Description to be printed *before* subcommand help.
    _help_title : optional[str]
        Human-readable title of the subcommand section in help output,
        typically only one to three words of lowercase, unpunctuated text.
        Defaults to simply ``subcommands``.
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(
        self,

        # Mandatory arguments.
        subcommand_var_name: str,
        subcommands: SequenceTypes,
        help_description: str,

        # Optional arguments.
        help_title: str = 'subcommands',
    ) -> None:
        '''
        Initialize this container of CLI subcommands.

        Parameters
        ----------
        subcommand_var_name : str
            Name of the instance variable to which the :class:`ArgumentParser`
            parsing these subcommands from command-line options and arguments
            passed by the external user stores either:
            * The name of the currently passed subcommand for this container if
              the user passed such a subcommand.
            * ``None`` otherwise.
        subcommands : SequenceTypes
            Sequence of all **subcommands** (i.e., :class:`CLISubcommandABC`
            instances) to add to this container.
        help_description : str
            Description to be printed *before* subcommand help. For convenience,
            all format substrings supported by the :meth:`CLIABC.expand_help`
            method (e.g., ``{program_name}``) are globally replaced as expected.
        help_title : optional[str]
            Human-readable title of the subcommand section in help output,
            typically only one to three words of lowercase, unpunctuated text.
            Defaults to simply ``subcommands``. As in the ``description``
            parameter, all format substrings are globally replaced as expected.
        '''

        # Classify all passed parameters.
        self._subcommand_var_name = subcommand_var_name
        self._subcommands = subcommands
        self._help_description = help_description
        self._help_title = help_title

    # ..................{ ADDERS                             }..................
    @type_check
    def add(
        self,

        # Avoid circular import dependencies.
        cli: 'betse.cli.api.cliabc.CLIABC',
        arg_parser: ArgParserType,
        arg_subparsers_kwargs: MappingType,
    ) -> ArgSubparsersType:
        '''
        Create a new **argument subparsers container** (i.e., child
        :mod:`argparse`-specific object containing each argument subparser)
        parsing all subcommands in this container, add this container to the
        passed **argument parser** (i.e., parent :mod:`argparse`-specific
        object parsing command-line arguments), and return this container.

        Parameters
        ----------
        cli : betse.cli.api.cliabc.CLIABC
            High-level command-line interface (CLI) owning this container. To
            avoid circular references, neither this method nor any method
            transitively called by this method retains this reference.
        arg_parser : ArgParserType
            Parent argument parser to which the container of subcommand argument
            subparsers created by this method is added. If this parser is:
            * The root argument parser for this CLI, all subcommands in this
              container will be created as top-level subcommands.
            * Any other argument parser, all subcommands in this container will
              be created as child subcommands of another parent subcommand.

        Returns
        ----------
        ArgSubparsersType
            Container of subcommand argument subparsers created by this method.
        '''

        # Keyword arguments with which to initialize this container,
        # interpolating all format substrings in all human-readable arguments.
        kwargs = {
            'dest':        self._subcommand_var_name,
            'title':       cli.expand_help(self._help_title),
            'description': cli.expand_help(self._help_description),
        }
        kwargs.update(arg_subparsers_kwargs)

        # Container of subcommand argument subparsers to be returned.
        arg_subparsers = arg_parser.add_subparsers(**kwargs)

        # For each contained subcommand, add an argument parser parsing this
        # subcommand to this container of argument subparsers.
        for subcommand in self._subcommands:
            subcommand.add(cli=cli, arg_subparsers=arg_subparsers)

        # Return this container.
        return arg_subparsers

# ....................{ SUPERCLASSES                       }....................
class CLISubcommandABC(object, metaclass=ABCMeta):
    '''
    Abstract base class of all **CLI subcommand** (i.e., argument passed to
    the external ``betse`` command signifying a high-level action to perform)
    subclasses.

    This class encapsulates all metadata pertaining to this subcommand,
    including:

    * Human-readable help strings describing this subcommand.
    * All options and arguments accepted by this subcommand.

    Attributes
    ----------
    _name : str
        Machine-readable name of this CLI subcommand (e.g., ``plot``), typically
        only a single word.
    _help_synopsis : str
        Human-readable synopsis of this CLI subcommand, typically only one to
        three lines of lowercase, unpunctuated text.
    _help_description : str
        Human-readable description of this CLI subcommand, typically one to
        several paragraphs of grammatical sentences.
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
        Initialize this CLI subcommand.

        Parameters
        ----------
        name : str
            Machine-readable name of this CLI subcommand (e.g., ``plot``), typically
            only a single word.
        synopsis : str
            Human-readable synopsis of this CLI subcommand, typically only one
            to three lines of lowercase, unpunctuated text. For convenience, all
            format substrings supported by the :meth:`CLIABC.expand_help` method
            (e.g., ``{program_name}``) are globally replaced as expected.
        description : str
            Human-readable description of this CLI subcommand, typically one to
            several paragraphs of grammatical sentences. As in the ``synopsis``
            parameter, all format substrings are globally replaced as expected.
        '''

        # Classify all passed parameters.
        self._name = name
        self._help_synopsis = synopsis
        self._help_description = description

    # ..................{ ADDERS                             }..................
    @type_check
    def add(
        self,

        # Avoid circular import dependencies.
        cli: 'betse.cli.api.cliabc.CLIABC',
        arg_subparsers: ArgSubparsersType,
    ) -> ArgParserType:
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
        cli : betse.cli.api.cliabc.CLIABC
            High-level command-line interface (CLI) owning this subcommand. To
            avoid circular references, neither this method nor any method
            transitively called by this method retains this reference.
        arg_subparsers : ArgSubparsersType
            Collection of sibling subcommand argument parsers to which the
            subcommand argument parser created by this method is added. This
            collection is owned either by:
            * A top-level subcommand (e.g., ``plot``), in which case the
              subcommand created by this method is a child of that subcommand.
            * No subcommand, in which case the subcommand created by this method
              is a top-level subcommand.

        Returns
        ----------
        ArgParserType
            Subcommand argument parser created by this method.
        '''

        # Keyword arguments with which to initialize this container,
        # interpolating all format substrings in all human-readable arguments.
        kwargs = {
            'name':        self._name,
            'help':        cli.expand_help(self._help_synopsis),
            'description': cli.expand_help(self._help_description),
        }

        # Merge these subcommand-specific arguments with all default arguments.
        kwargs.update(cli.arg_parser_kwargs)

        # Create and return this parser, added to this container of subparsers.
        return arg_subparsers.add_parser(**kwargs)

# ....................{ SUBCLASSES                         }....................
class CLISubcommandNoArg(CLISubcommandABC):
    '''
    CLI subcommand accepting *no* passed arguments.
    '''

    pass


#FIXME: Replace all usage of this class by "CLISubcommandParent".
class CLISubcommandParentObsolete(CLISubcommandABC):
    '''
    CLI subcommand that is itself the parent of one or more CLI subcommands,
    accepting *only* the name of a child subcommand as a passed argument.
    '''

    pass


class CLISubcommandParent(CLISubcommandABC):
    '''
    CLI subcommand which itself is the parent of one or more CLI subcommands,
    accepting *only* the name of a child subcommand as a passed argument.

    Parameters
    ----------
    _subcommands : CLISubcommands
        Container of child subcommands accepted by this parent subcommand.
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(
        self,
        subcommands: CLISubcommands,
        *args, **kwargs
    ) -> None:
        '''
        Initialize this CLI subcommand.

        Parameters
        ----------
        subcommands : CLISubcommands
            Container of child subcommands accepted by this parent subcommand.

        All remaining parameters are passed as is to the
        :meth:`CLISubcommandABC.__init__` method.
        '''

        # Initialize our superclass with all passed parameters.
        super().__init__(*args, **kwargs)

        # Classify all remaining parameters.
        self._subcommands = subcommands

    # ..................{ ADDERS                             }..................
    @type_check
    def add(
        self,

        # Avoid circular import dependencies.
        cli: 'betse.cli.api.cliabc.CLIABC',
        *args, **kwargs
    ) -> ArgParserType:

        # Argument parser parsing this subcommand.
        arg_parser = self.add(*args, **kwargs)

        # Add this container of argument subparsers to this parser.
        self._subcommands.add(cli=cli, arg_parser=arg_parser)

        # Return this argument parser.
        return arg_parser


class CLISubcommandYAMLOnly(CLISubcommandABC):
    '''
    CLI subcommand accepting *only* a configuration filename as a passed
    argument.
    '''

    # ..................{ ADDERS                             }..................
    def add(self, *args, **kwargs) -> ArgParserType:

        # Subcommand argument subparser added by our superclass.
        arg_subparser = super().add(*args, **kwargs)

        # Configure this subparser to require a configuration file argument.
        arg_subparser.add_argument(
            'conf_filename',
            metavar='CONFIG_FILE',
            help='simulation configuration file',
        )

        # Return this subparser.
        return arg_subparser
