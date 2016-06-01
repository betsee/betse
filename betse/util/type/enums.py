#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Low-level **enumeration** (i.e., enumerable types created by the `Enum()` class)
facilities.
'''

# ....................{ IMPORTS                            }....................
from betse.util.type import types

# ....................{ GETTERS                            }....................
def get_names_lowercase(enum: "EnumMeta") -> 'collections.Sequence':
    '''
    Non-string sequence of the lowercased names of all members of the passed
    enumeration type in **declaration order** (i.e., the order in which these
    members were originally declared).

    Parameters
    ----------
    enum : EnumMeta
        Enumeration type to be inspected.

    Returns
    ----------
    collections.Sequence
        Lowercased names of all enumeration members.
    '''
    assert types.is_enum(enum), types.assert_not_enum(enum)

    return tuple(enum_member.name.lower() for enum_member in enum)
