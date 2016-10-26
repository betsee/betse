#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Low-level **pickling** (i.e., serialization and deserialization of arbitrarily
complex objects to and from on-disk files) facilities.

Caveats
----------
**This submodule leverages the third-party `dill` package rather than the
standard `pickle` package.** The former conforms to the API of the latter with
additional support for so-called "exotic" types required by this application,
including:

* `lambda` expressions.
* Generators.
* Ranges.
* Slices.
* Numpy `ndarray` subclass objectes.
* Numpy `unfunc` objects.
'''

#FIXME: Add support to:
#
#* The save() function for transparently compressing a pickled file immediately
#  after pickling that file *IF* the passed filename has a filetype signifying a
#  supported compression format (e.g., "bz2", "gz", "zip").
#* The load() function for transparently decompressing a pickled file
#  immediately before unpickling that file *IF* the passed filename again has a
#  filetype signifying a supported compression format.
#
#For relevant examples, see:
#    https://stackoverflow.com/questions/18474791/decreasing-the-size-of-cpickle-objects

# ....................{ IMPORTS                            }....................
import dill as pickle
from betse.util.io.log import logs
from betse.util.path import files
from betse.util.type.types import type_check

# ....................{ CONSTANTS                          }....................
# The improved pickle-ability of protocol 4 appears to be required to pickle
# C-based data structures (e.g., "scipy.spatial.KDTree").
PROTOCOL = 4
'''
Pickle protocol used by the :func:`save` and :func:`load` functions.

This protocol is the most recent pickle protocol supported by the minimum
version of Python supported by this application, satisfying the following two
competing tradeoffs:

* Compatibility with all versions of Python supported by this application.
* Maximal **pickle-ability** (i.e., the capacity to pickle objects), improving
  support for such edge cases as very large objects and edge-case object types.
'''

# ....................{ SAVERS                             }....................
@type_check
def save(*objects: object, filename: str) -> None:
    '''
    Save (i.e., write, pickle, serialize) the tuple of all passed objects to the
    file with the passed path if two or more objects are passed _or_ the single
    passed object if only one object is passed.

    Parameters
    ----------
    objects : tuple
        One or more arbitrarily complex object to be serialized. These objects
        and all objects transitively referenced by this object will be
        serialized to this file. If:
        * Only one object is passed, only that object will be saved.
        * Two or more objects are passed, the tuple of all such objects will be
          saved.
    filename : str
        Absolute or relative path of this file.
    '''

    # If only one object is passed, save only that object rather than the
    # 1-tuple consisting only of that object.
    if len(objects) == 1:
        objects = objects[0]

    # Log this pickling.
    logs.log_debug('Saving "%s"...', filename)

    # Save these objects to this file.
    with files.write_bytes(filename) as pickle_file:
        pickle.dump(
            objects,
            file=pickle_file,
            protocol=PROTOCOL,

            #FIXME: We may need to increase the maximum recursion depth *BEFORE*
            #calling pickle.dump() above. Let's ignore this until issues arise.

            # Physically pickle the contents of all objects transitively
            # referring to globals via stack-based recursion. By default, dill
            # non-physically pickles these objects by reference. The latter
            # approach has the benefit of avoiding recursion and hence fatal
            # errors on exceeding Python's maximum recursion depth, but the
            # distinct disadvantage of preventing pickled objects referring to
            # modules whose fully-qualified names have since changed (e.g., due
            # to those modules having since been moved, renamed, or removed)
            # from being unpickled.
            #
            # In our case, the latter is of considerably more concern than the
            # former. Whereas the former is trivially solved by increasing the
            # maximum recursion depth, the latter has no sane solution. An
            # unpicklable object is simply unpicklable. Since this codebase is
            # under constant scrutiny and construction, the fully-qualified
            # names of modules frequently change.
            #
            # To safeguard backward compatibility with previously pickled
            # objects, recursive object discovery is *STRONGLY* preferred.
            recurse=True,
        )

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

    # Log this unpickling.
    logs.log_debug('Loading "%s"...', filename)

    # Load and return all objects saved to this file.
    with files.read_bytes(filename) as unpickle_file:
        return pickle.load(file=unpickle_file)
