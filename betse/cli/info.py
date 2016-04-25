#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
`info` subcommand for `betse`'s command line interface (CLI).
'''

#FIXME; For aesthetics, convert to yppy-style "cli.memory_table" output.

# ....................{ IMPORTS                            }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: To avoid non-trivial delays on importing this module, this module
# should import *ONLY* from modules and packages whose importation is unlikely
# to incur such delays. This includes all standard Python packages and all BETSE
# packages required by the log_info_header() function, which is called
# sufficiently early in application startup as to render these imports
# effectively mandatory.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

from collections import OrderedDict
from io import StringIO

from betse import metadata
from betse.util.io.log import log_config, logs
from betse.util.os import oses
from betse.util.py import pys


# ..................{ LOGGERS                                }..................
def log_header() -> None:
    '''
    Log a one-line synopsis of metadata logged by the `info` subcommand.
    '''

    logs.log_info(
        'Welcome to <<'
        '{program_name} {program_version} | '
        '{py_name} {py_version} | '
        '{os_name} {os_version}'
        '>>.'.format(
            program_name=metadata.NAME,
            program_version=metadata.__version__,
            py_name=pys.get_name(),
            py_version=pys.get_version(),
            os_name=oses.get_name(),
            os_version=oses.get_version(),
        ))

# ..................{ OUTPUTTERS                             }..................
def output_info() -> None:
    '''
    Print all output for the `info` subcommand.
    '''

    # Notify the current user of a possible wait *BEFORE* importing modules
    # whose importation contributes to this wait.
    logs.log_info(
        'Harvesting system information. (This may take a moment.)')

    # Avoid inefficient import delays.
    from betse import pathtree
    from betse.lib import libs
    from betse.lib.matplotlib.matplotlibs import mpl_config
    from betse.util.os import kernels
    from betse.util.process import processes

    # Dictionary of human-readable labels to dictionaries of all
    # human-readable keys and values categorized by such labels. All such
    # dictionaries are ordered so as to preserve order in output.
    info_type_to_dict = OrderedDict((
        # Application metadata.
        (metadata.NAME.lower(), OrderedDict((
            ('basename', processes.get_current_basename()),
            ('version', metadata.__version__),
            ('authors', metadata.AUTHORS),
            ('home directory', pathtree.HOME_DIRNAME),
            ('dot directory',  pathtree.DOT_DIRNAME),
            ('data directory', pathtree.DATA_DIRNAME),
            ('default config file', pathtree.CONFIG_DEFAULT_FILENAME),
            ('default config geometry directory',
             pathtree.DATA_GEOMETRY_DIRNAME),
        ))),

        # Logging metadata.
        ('logging', log_config.get_metadata()),

        # Dependencies metadata.
        ('dependencies', libs.get_metadata()),

        # "matplotlib" metadata.
        ('matplotlib', mpl_config.get_metadata()),

        # Python metadata.
        ('python', pys.get_metadata()),

        # Operating system (OS) metadata.
        ('kernel', kernels.get_metadata()),
        ('os', oses.get_metadata()),
    ))

    # String buffer formatting such information.
    info_buffer = StringIO()

    # True if this is the first label to be output.
    is_info_type_first = True

    # Format each such dictionary under its categorizing label.
    for info_type, info_dict in info_type_to_dict.items():
        # If this is *NOT* the first label, delimit this label from the
        # prior label.
        if is_info_type_first:
            is_info_type_first = False
        else:
            info_buffer.write('\n')

        # Format such label.
        info_buffer.write('{}:\n'.format(info_type))

        # Format such dictionary.
        for info_key, info_value in info_dict.items():
            info_buffer.write('  {}: {}\n'.format(info_key, info_value))

    # Log rather than merely output such string, as logging simplifies
    # cliest-side bug reporting.
    logs.log_info(info_buffer.getvalue())
