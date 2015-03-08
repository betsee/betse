#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Low-level container facilities.
'''

# ....................{ IMPORTS                            }....................
from collections import Iterable

# ....................{ TESTERS                            }....................
def is_iterable(obj) -> bool:
    '''
    True if the passed object is an *interable* (i.e., implements the abstract
    interface defined by `collections.Iterable`).
    '''
    return isinstance(obj, Iterable)

def is_iterable_nonstring(obj) -> bool:
    '''
    True if the passed object is a *non-string interable* (i.e., implements the
    abstract interface defined by `collections.Iterable` *and* is not a string).
    '''
    return not isinstance(obj, str) and is_iterable(obj)

# ....................{ SORTERS                            }....................
def sort_as_lexicographic_ascending(container):
    '''
    Sort the passed container in **ascending lexicographic order** (i.e., in the
    order of conventional dead-tree dictionaries).
    '''
    return sorted(container)

# --------------------( WASTELANDS                         )--------------------
    # return hasattr(obj, '__iter__')
