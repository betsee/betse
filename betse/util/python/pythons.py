#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
High-level Python facilities.

This module provides functionality pertaining to the active Python interpreter
as a whole.
'''

# ....................{ IMPORTS                            }....................
from collections import OrderedDict
import platform

# ....................{ GETTERS                            }....................
def get_metadata() -> str:
    '''
    Get an ordered dictionary synopsizing the active Python interpreter.

    This function aggregates the metadata reported by the reasonably
    cross-platform module `platform` into a simple dictionary.
    '''
    # Such dictionary.
    metadata = OrderedDict((
        ('implementation', platform.python_implementation()),
        ('version', platform.python_version()),
        ('vcs revision', platform.python_revision()),
        ('vcs branch', platform.python_branch()),
        ('compiler', platform.python_compiler()),
    ))

    # 2-tuple providing such interpreter's build number and date as strings.
    python_build = platform.python_build()

    # Append such metadata.
    metadata['build number'] = python_build[0]
    metadata['build data'] = python_build[1]

    # Get such dictionary.
    return metadata

# --------------------( WASTELANDS                         )--------------------
    #FUXME: Not terribly human-readable and hence hardly ideal.
    # ordered_dict[' platform.platform(aliased = True)

# def format() -> str:
#     '''
#     Get a human-readable string synopsizing the current system.
#     '''
