#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Low-level dictionary facilities.
'''

# ....................{ IMPORTS                            }....................
import pprint
from betse.util.type import types

# ....................{ TESTERS                            }....................
def is_keys(mapping: dict, *key_names) -> str:
    '''
    `True` only if the passed dictionary contains _all_ passed keys.

    Parameters
    ----------
    mapping : dict
        Dictionary to be tested.
    key_names : list
        List of all keys to be tested for.
    '''
    assert types.is_mapping(mapping), types.assert_not_mapping(mapping)

    # Yes, this is ridiculously awesome.
    return set(key_names).issubset(mapping)

# ....................{ FORMATTERS                         }....................
def format(mapping: dict) -> str:
    '''
    Convert the passed dictionary to a human-readable string.
    '''
    assert types.is_mapping(mapping), types.assert_not_mapping(mapping)
    return pprint.pformat(mapping)
