#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Low-level command line interface (CLI) utilities.
'''

# ....................{ IMPORTS                            }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: To raise human-readable exceptions on application startup, the
# top-level of this module may import *ONLY* from submodules guaranteed to:
# * Exist, including standard Python and application modules.
# * Never raise exceptions on importation (e.g., due to module-level logic).
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

from betse import metadata
from betse.util.path.command import cmds
from betse.util.type.text import strs
from betse.util.type.types import (
    type_check, ArgSubparsersType, IterableTypes, MappingType,)

# ....................{ GETTERS                            }....................
def get_version() -> str:
    '''
    Human-readable version specifier suitable for printing to end users.
    '''

    return '{} {}'.format(cmds.get_current_basename(), metadata.VERSION)

# ....................{ ADDERS                             }....................
@type_check
def add_arg_subparsers_subcommands(
    # Avoid circular import dependencies.
    cli: 'betse.cli.api.cliabc.CLIABC',
    subcommands: IterableTypes,
    arg_subparsers: ArgSubparsersType,
) -> MappingType:
    '''
    Add one new argument subparser parsing each subcommand in the passed
    iterable of subcommands to the passed container of argument subparsers and
    return a dictionary mapping from the name of each such subcommand to the
    corresponding argument subparser.

    Parameters
    ----------
    cli : betse.cli.api.cliabc.CLIABC
        High-level command-line interface (CLI) owning this subcommand. To avoid
        circular references, neither this method nor any method transitively
        called by this method retains this reference.
    subcommands : IterableTypes
        Iterable of all subcommands to add subparsers for.
    arg_subparsers : ArgSubparsersType
        Container of argument subparsers to add these subparsers to.

    Returns
    ----------
    MappingType
        Dictionary mapping from the name of each such subcommand to the new
        argument subparser parsing that subcommand.
    '''

    # Dictionary to be returned.
    subcommand_name_to_subparser = {}

    # For each passed subcommand, create, add, and map an argument parser
    # parsing this subcommand to this container of argument subparsers.
    for subcommand in subcommands:
        subcommand_name_to_subparser[subcommand._name] = subcommand.add(
            cli=cli, arg_subparsers=arg_subparsers)

    # Return this dictionary.
    return subcommand_name_to_subparser

# ....................{ EXPANDERS                          }....................
#FIXME: Replace with usage of the CLIABC.expand_help() method; then excise this.
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

    return strs.remove_whitespace_presuffix(text.format(
        program_name=metadata.NAME,
        script_basename=cmds.get_current_basename(),
        **kwargs
    ))
