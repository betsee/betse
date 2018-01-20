#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Implementation of the ``info`` subcommand for this command line interface (CLI).
'''

#FIXME; For aesthetics, convert to yppy-style "cli.memory_table" output.

# ..................{ IMPORTS                                }..................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: To avoid non-trivial delays on importing this module, this module
# should import *ONLY* from modules and packages whose importation is unlikely
# to incur such delays. This includes all standard Python packages and all BETSE
# packages required by the log_info_header() function, which is called
# sufficiently early in application startup as to render these imports
# effectively mandatory.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

from betse import metadata
from betse.util.io.log import logconfig, logs
from betse.util.os import oses
from betse.util.py import pyimpl, pys
from betse.util.type.mapping.mapcls import OrderedArgsDict
from collections import OrderedDict

# ..................{ GETTERS                                }..................
def get_metadata() -> OrderedArgsDict:
    '''
    Ordered dictionary synopsizing the current installation of this application.
    '''

    # Defer heavyweight imports.
    from betse import pathtree
    from betse.util.path.command import cmds

    # Return this dictionary.
    return OrderedArgsDict(
        'basename', cmds.get_current_basename(),
        'version', metadata.VERSION,
        'codename', metadata.CODENAME,
        'authors', metadata.AUTHORS,
        'license', metadata.LICENSE,
        'home directory',  pathtree.get_home_dirname(),
        'dot directory',   pathtree.get_dot_dirname(),
        'data directory',  pathtree.get_data_dirname(),
    )

# ..................{ LOGGERS                                }..................
def log_header() -> None:
    '''
    Log a single-line human-readable sentence synopsizing the state of the
    current application (e.g., name, codename, version).
    '''

    logs.log_info(
        'Welcome to <<'
        '{program_name} {program_version} | '
        '{py_name} {py_version} | '
        '{os_name} {os_version}'
        '>>.'.format(
            program_name=metadata.NAME,
            program_version=metadata.VERSION,
            py_name=pyimpl.get_name(),
            py_version=pys.get_version(),
            os_name=oses.get_name(),
            os_version=oses.get_version(),
        ))


def log_info() -> None:
    '''
    Log all metadata required by the ``info`` subcommand.
    '''

    # Notify the current user of a possible wait *BEFORE* importing modules
    # whose importation contributes to this wait.
    logs.log_info('Harvesting system metadata... (This may take a moment.)')

    # Defer heavyweight imports.
    from betse.lib import libs
    from betse.util.os import displays, kernels

    #FIXME: Shift into a more appropriate general-purpose submodule.
    # Tuple of BETSE-specific metadata.
    BETSE_METADATAS = (
        # Application metadata.
        (metadata.NAME.lower(), get_metadata()),

        # Logging metadata.
        ('logging', logconfig.get_metadata()),
    )

    #FIXME: Shift into a more appropriate general-purpose submodule.
    # Tuple of system-specific metadata.
    SYSTEM_METADATAS = (
        # Python metadata.
        ('python', pys.get_metadata()),
        ('python interpreter', pyimpl.get_metadata()),

        # Operating system (OS) metadata.
        ('os', oses.get_metadata()),
        ('os kernel', kernels.get_metadata()),
        ('os display', displays.get_metadata()),
    )

    # Dictionary of human-readable labels to dictionaries of all
    # human-readable keys and values categorized by such labels. All such
    # dictionaries are ordered so as to preserve order in output.
    info_type_to_dict = OrderedDict(
        BETSE_METADATAS +
        libs.get_metadatas() +
        SYSTEM_METADATAS
    )

    # String formatting this information.
    info_buffer = 'Harvested system information:\n'

    # Format each such dictionary under its categorizing label.
    for info_type, info_dict in info_type_to_dict.items():
        # Format this label.
        info_buffer += '\n{}:'.format(info_type)

        # Format this label's dictionary.
        for info_key, info_value in info_dict.items():
            info_buffer += '\n  {}: {}'.format(info_key, info_value)

    # Log rather than merely output this string, as logging simplifies
    # cliest-side bug reporting.
    logs.log_info(info_buffer)
