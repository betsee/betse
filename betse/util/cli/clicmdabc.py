#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Top-level abstract base class of all **subcommandable command line interface
(CLI)** (i.e., CLI accepting one or more subcommands) subclasses.
'''

# ....................{ IMPORTS                            }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: To raise human-readable exceptions on application startup, the
# top-level of this module may import *ONLY* from submodules guaranteed to:
# * Exist, including standard Python and application modules.
# * Never raise exceptions on importation (e.g., due to module-level logic).
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

from abc import abstractmethod
from betse.util.cli.cliabc import CLIABC
from betse.util.cli.clicmd import CLISubcommander, CLISubcommandParent
from betse.util.py import pyident
from betse.util.type.obj import objects
from betse.util.type.types import type_check, ArgParserType

# ....................{ SUBCLASS                           }....................
class CLISubcommandableABC(CLIABC):
    '''
    Top-level abstract base class of all **subcommandable command line interface
    (CLI)** (i.e., CLI accepting one or more subcommands) subclasses, suitable
    for use by both CLI and GUI front-ends for BETSE.

    Unlike the parent :class:`CLIABC` superclass, this superclass provides
    explicit support for subcommands. Concrete subclasses implementing
    subcommands should directly subclass this rather than that superclass.

    Attributes
    ----------
    _subcommander_top : CLISubcommander
        Container of all argument subparsers parsing all top-level subcommands
        accepted by this application (e.g., ``init``, ``sim``).
    '''

    # ..................{ INITIALIZERS                       }..................
    def __init__(self) -> None:

        # Initialize our superclass.
        super().__init__()

        # Nullify all instance variables for safety.
        self._subcommander_top = None

    # ..................{ SUPERCLASS ~ args                  }..................
    def _config_arg_parsing(self) -> None:

        # Container of all top-level argument subparsers for this application.
        self._subcommander_top = self._make_subcommander_top()
        self._subcommander_top.add(cli=self, arg_parser=self._arg_parser_top)

    # ....................{ SUBCLASS ~ mandatory               }....................
    # The following methods *MUST* be implemented by subclasses.

    @abstractmethod
    def _make_subcommander_top(self) -> CLISubcommander:
        '''
        Container of all top-level subcommands accepted by this CLI command.

        **Order is significant,** defining the order that the ``--help`` option
        synopsizes these subcommands in. Subcommands omitted from this container
        will *not* be parsed by argument subparsers and thus ignored.

        For each such subcommand, the :meth:`_config_arg_parsing` method creates
        and adds a corresponding argument subparser to the lower-level container
        of all top-level argument subparsers.
        '''

        pass

    # ..................{ SUPERCLASS ~ cli                   }..................
    def _do(self) -> object:
        '''
        Implement this command-line interface (CLI).

        If the caller passed a top-level subcommand to this CLI command, this
        method runs this subcommand and returns the result of doing so; else,
        this method prints help and returns the current instance of this CLI.
        '''

        self._run_subcommand(
            arg_parser=self._arg_parser_top,
            subcommander=self._subcommander_top,
            subcommand_method_name_prefix='_do_',
        )


    @type_check
    def _run_subcommand(
        self,
        arg_parser: ArgParserType,
        subcommander: CLISubcommander,
        subcommand_method_name_prefix: str,
    ) -> object:
        '''
        Recursively run the child subcommand parsed by the passed argument
        parser, contained by the passed parent subcommands container, and
        selected by the external caller on the command line and return the
        result of doing so if the caller passed a child subcommand for this
        parent subcommand *or* print help for this parent subcommand and return
        the current instance of this CLI.

        To avoid conflict with the ``_do_``-prefixed names of subcommand
        methods, this method is intentionally *not* prefixed by ``_do_``.

        Parameters
        ----------
        arg_parser : ArgParserType
            Argument parsing parsing the parent subcommand of these child
            subcommands.
        subcommander : CLISubcommander
            Container of these child subcommands.
        subcommand_method_name_prefix : str
            Substring prefixing the names of all subclass methods implementing
            these child subcommands.
        '''

        # Name of the child subcommand for this parent subcommand passed by the
        # external user on the command line if a child subcommand was passed
        # *OR* "None" otherwise.
        subcommand_name = getattr(self._args, subcommander.subcommand_var_name)

        # If no child subcommand was passed...
        if not subcommand_name:
            #,Print help for this parent subcommand. Note that this common case
            #constitutes neither a fatal error nor a non-fatal warning.
            print()
            arg_parser.print_help()

            # Return the current instance of this CLI. While trivial, doing so
            # simplifies memory profiling of this CLI.
            return self
        # Else, a child subcommand was passed.

        # Child subcommand and argument parser corresponding to this name.
        subcommand = self._subcommander_top.subcommand_name_to_subcommand[
            subcommand_name]
        subcommand_arg_parser = (
            self._subcommander_top.subcommand_name_to_arg_parser[
                subcommand_name])

        # Snake_case-style name of this subcommand, suitable for use in
        # constructing the syntactically valid Python method name.
        subcommand_name_snakecase = pyident.sanitize_snakecase(subcommand_name)

        # If this subcommand is itself the parent of child subcommands...
        if isinstance(subcommand, CLISubcommandParent):
            # Append the prefix of the names of all subclass methods
            # implementing these child subcommands with this subcommand's name.
            subcommand_method_name_prefix += subcommand_name_snakecase + '_'

            # Call this function recursively and return the result of doing so.
            return self._run_subcommand(
                arg_parser=subcommand_arg_parser,
                subcommander=subcommand.subcommander,
                subcommand_method_name_prefix=subcommand_method_name_prefix,
            )
        # Else, this subcommand is a leaf rather than parent. In this case,
        # terminate recursion by running this subcommand's subclass method.
        else:
            # Name of the method running this subcommand.
            subcommand_method_name = (
                subcommand_method_name_prefix + subcommand_name_snakecase)

            # Method running this subcommand. If this method does *NOT* exist,
            # get_method() will raise a non-human-readable exception. Usually,
            # that would be bad. In this case, argument parsing coupled with a
            # reliable class implementation ensures this method to exist.
            subcommand_method = objects.get_method(
                obj=self, method_name=subcommand_method_name)

            # Run this subcommand and return the result of doing so (if any).
            return subcommand_method()
