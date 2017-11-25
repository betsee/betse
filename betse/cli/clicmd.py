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
from betse.cli import cliutil
from betse.util.type.types import (
    type_check, ArgParserType, ArgSubparsersType, IterableTypes, MappingType)

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
    name : str
        Machine-readable name of this CLI subcommand (e.g., ``plot``), typically
        only a single word.
    synopsis : str
        Human-readable synopsis of this CLI subcommand, typically only one to
        three lines of lowercase, unpunctuated text. All ``{``- and ``}``-
        delimited format substrings (e.g., ``{program_name}``) supported by the
        :meth:`cliutil.expand_help` function will be globally replaced.
    description : str
        Human-readable description of this CLI subcommand, typically one to
        several paragraphs of grammatical sentences. All ``{``- and ``}``-
        delimited format substrings (e.g., ``{program_name}``) supported by the
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
            * A top-level subcommand (e.g., ``plot``), in which case the
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
    CLI subcommand accepting *no* passed arguments.
    '''

    pass


class CLISubcommandParent(CLISubcommandABC):
    '''
    CLI subcommand that is itself the parent of one or more CLI subcommands,
    accepting *only* the name of a child subcommand as a passed argument.
    '''

    # We almost don't believe it either.
    pass


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
