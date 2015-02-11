#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Low-level string facilities.
'''

# ....................{ IMPORTS                            }....................
# import sys

# ....................{ JOINERS                            }....................
def join(*strs) -> str:
    '''
    Concatenate the passed strings with no separating delimiter.

    This is a convenience function wrapping the standard `"".join(())` method,
    whose syntax is arguably overly verbose.
    '''
    return ''.join((*strs))

# --------------------( WASTELANDS                         )--------------------
