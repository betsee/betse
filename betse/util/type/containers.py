#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Low-level container facilities.
'''

# ....................{ SORTERS                            }....................
def sort_as_lexicographic_ascending(container):
    '''
    Sort the passed container in **ascending lexicographic order** (i.e., in the
    order of conventional dead-tree dictionaries).
    '''
    return sorted(container)

# --------------------( WASTELANDS                         )--------------------
    # return hasattr(obj, '__iter__')
