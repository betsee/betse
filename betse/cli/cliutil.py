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
from betse.util.type.types import type_check

# ....................{ EXPANDERS                          }....................
def get_version() -> str:
    '''
    Human-readable version specifier suitable for printing to end users.
    '''

    return '{} {}'.format(
        cmds.get_current_basename(), metadata.VERSION)

# ....................{ EXPANDERS                          }....................
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
