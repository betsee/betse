#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Low-level dictionary facilities.
'''

# ....................{ IMPORTS                            }....................
import pprint

# ....................{ JOINERS                            }....................
def format(map: dict) -> str:
    '''
    Convert the passed dictionary to a human-readable string.
    '''
    assert isinstance(map, dict), '"{}" not a dictionary.'.format(map)
    return pprint.pformat(map)

# --------------------( WASTELANDS                         )--------------------
