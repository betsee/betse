#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
High-level system facilities.

This module provides functionality pertaining to the system as a whole and hence
comprising both hardware and software. Hardware- and software-specific
functionality resides elsewhere.
'''

# ....................{ IMPORTS                            }....................
from collections import OrderedDict
import platform

# ....................{ GETTERS                            }....................
def get_metadata() -> OrderedDict:
    '''
    Get an ordered dictionary synopsizing the current system.

    This function expands the metadata reported by the reasonably cross-platform
    function `platform.uname()` with additional diagnostics.
    '''
    # Such dictionary.
    metadata = vars(platform.uname())

    #FIXME: Expand such dictionary here with additioral metadata. See:
    #    https://docs.python.org/3/library/platform.html

    # Get such dictionary.
    return metadata

# --------------------( WASTELANDS                         )--------------------
    #FUXME: Not terribly human-readable and hence hardly ideal.
    # metadata[' platform.platform(aliased = True)

# def format() -> str:
#     '''
#     Get a human-readable string synopsizing the current system.
#     '''
