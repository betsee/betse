#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
**Pickling** (i.e., serialization and deserialization of arbitrarily complex
objects to and from on-disk files) facilities.
'''

# ....................{ IMPORTS                            }....................
import pickle
from betse.util.type.types import type_check

# ....................{ CONSTANTS                          }....................
# The improved pickle-ability of protocol 4 appears to be required to pickle
# C-based data structures (e.g., "scipy.spatial.KDTree").
PROTOCOL = 4
'''
Pickle protocol used by the `save()` and `load()` functions.

This protocol is the most recent pickle protocol supported by the minimum
version of Python supported by this application, satisfying the following two
competing tradeoffs:

* Compatibility with all versions of Python supported by this application.
* Maximal **pickle-ability** (i.e., the capacity to pickle objects), improving
  support for such edge cases as very large objects and edge-case object types.
'''

# ....................{ SAVERS                             }....................
@type_check
def save(obj: object, filename: str) -> None:
    '''
    Save (i.e., write, pickle, serialize) the passed object to the file with the
    passed path.

    Parameters
    ----------
    obj : object
        Arbitrarily complex object to be serialized. This object and all objects
        transitively referenced by this object will be serialized to this file.
    filename : str
        Absolute or relative path of this file.
    '''

    with open(filename, 'wb') as pickle_file:
        pickle.dump(obj, file=pickle_file, protocol=PROTOCOL)

# ....................{ LOADERS                            }....................
@type_check
def load(filename: str) -> object:
    '''
    Load (i.e., read, unpickle, deserialize) the object previously saved to the
    file with the passed path.

    Parameters
    ----------
    filename : str
        Absolute or relative path of this file.

    Returns
    ----------
    object
        Arbitrarily complex object and all objects transitively referenced by
        this object loaded from this file.
    '''

    with open(filename, 'rb') as unpickle_file:
        return pickle.load(file=unpickle_file)
