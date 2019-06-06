#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level **pickling** (i.e., serialization and deserialization of arbitrarily
complex objects to and from on-disk files) facilities.

Caveats
----------
**This submodule leverages the third-party :mod:`dill` package rather than the
standard :mod:`pickle` package.** The former conforms to the API of the latter
with additional support for so-called "exotic" types required by this
application, including:

* Generators.
* Lambda expressions.
* Ranges.
* Slices.
* Numpy :class:`ndarray` subclass instances.
* Numpy :class:`ufunc` objects.
'''

# ....................{ IMPORTS                           }....................
import dill
from betse.util.io import iofiles
from betse.util.io.log import logs
from betse.util.type.decorator.decmemo import CALLABLE_CACHED_VAR_NAME_PREFIX
from betse.util.type.obj import objtest
from betse.util.type.types import type_check

# ....................{ CONSTANTS                         }....................
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

# ....................{ CLASSES                           }....................
#FIXME: Given the unreliability of the private "dill._dill" API, it would be
#both substantially safer *AND* simpler to directly monkey-patch the existing
#dill.Pickler.save() method with that defined below rather than attempting to
#entirely replace the "dill._dill.Pickler" class as we currently do. Food for
#contemplative thought -- especially if dill breaks backward compatibility in
#this respect again.
class BetsePickler(dill.Pickler):
    '''
    Application-specific :mod:`dill`-based custom pickler.

    This pickler augments :mod:`dill` with additional support for
    application-specific constructs, including:

    * Preventing temporary in-memory cached data from being pickled to disk,
      including all private instance variables cached by decorators defined by
      the :mod:`betse.util.type.decorator.decmemo` module (e.g.,
      :func:`property_cached`). To do so efficiently, this pickler uncaches
      *all* previously cached data from *all* objects pickled to disk. This
      data is guaranteed to be transparently re-cached on the next in-memory
      access of this data and is thus safely uncachable. While technically
      avoidable (e.g., by saving and restoring uncached instance variables into
      a local dictionary internal to the :meth:`save` method), doing so would
      incur additional space, time, and maintenance penalties. In short, the
      lazy way still remains the best way.

    See Also
    ----------
    https://github.com/uqfoundation/dill/issues/225#issuecomment-294286518
        Feature request response by GitHub user matsjoyce_ on the :mod:`dill`
        issue tracker from which this implementation was strongly inspired.
        Thanks a metric ton, matsjoyce_!

    .. _matsjoyce:
        https://github.com/matsjoyce
    '''

    # ..................{ SAVERS                            }..................
    def save(self, obj, *args, **kwargs):
        '''
        Prepare the passed object to be pickled.

        This method is recursively called ala the Visitor pattern for each
        object to be pickled reachable in the current object graph.
        '''
        # logs.log_debug('Pickling object "%s"...', obj.__class__.__name__)

        # If this object defines a dictionary mapping from the name to value of
        # each attribute defined on this object and hence is unslotted, prevent
        # cached attributes from being pickled. Slotted objects cannot have
        # attributes dynamically added or removed at runtime and hence *MUST*
        # be ignored here.
        #
        # Since this guarantees the "obj.__dict__" attribute to exist, this
        # attribute is accessed directly below rather than indirectly via the
        # vars() builtin. While feasible, the latter is mildly less efficient.
        if hasattr(obj, '__dict__'):
            # For the name of each such attribute...
            for obj_attr_name in obj.__dict__.keys():
                # If this attribute is prefixed by a substring implying this
                # attribute to be a private instance variable to which some
                # caching decorators (e.g., @property_cached) has cached the
                # value returned by a decorated callable...
                if obj_attr_name.startswith(CALLABLE_CACHED_VAR_NAME_PREFIX):
                    # Log this deletion attempt.
                    logs.log_debug(
                        'Uncaching transient data "%s.%s"...',
                        obj.__class__.__name__, obj_attr_name)

                    # Delete this variable from this object. Again, do so
                    # directly rather than via the del() built. While feasible,
                    # the latter is slightly less efficient.
                    obj.__dict__.pop(obj_attr_name)

        # Pickle this object.
        super().save(obj, *args, **kwargs)

# ....................{ LOADERS                           }....................
@type_check
def load(filename: str) -> object:
    '''
    Load (i.e., read, unpickle, deserialize) the object previously saved to the
    file with the passed path.

    This function transparently decompresses these objects into this file when
    this filename's filetype is that of a supported archive format.

    Parameters
    ----------
    filename : str
        Absolute or relative path of this file. If this filename is suffixed by
        a supported archive filetype (i.e., if the
        :func:`betse.util.path.archives.is_filetype` function returns
        ``True`` when passed this filename), this file is automatically
        decompressed as an archive of that filetype.

    Returns
    ----------
    object
        Arbitrarily complex object and all objects transitively referenced by
        this object loaded from this file.
    '''

    # Load and return all objects saved to this file, silently decompressing
    # this file if compressed.
    with iofiles.reading_bytes(filename) as unpickle_file:
        return dill.load(file=unpickle_file)

# ....................{ SAVERS                            }....................
@type_check
def save(
    *objs,
    filename: str,
    is_overwritable: bool = False
) -> None:
    '''
    Save (i.e., write, pickle, serialize) the tuple of all passed objects to
    the file with the passed path if two or more objects are passed *or* the
    single passed object if only one object is passed.

    This function transparently compresses these objects into this file when
    this filename's filetype is that of a supported archive format.

    Parameters
    ----------
    objs : tuple
        One or more arbitrarily complex object to be serialized. These objects
        and all objects transitively referenced by this object will be
        serialized to this file. If:

        * Only one object is passed, only that object will be saved.
        * Two or more objects are passed, the tuple of all such objects will be
          saved.
    filename : str
        Absolute or relative path of this file. If this filename is suffixed by
        a supported archive filetype (i.e., if the
        :func:`betse.util.path.archives.is_filetype` function returns
        ``True`` when passed this filename), this file is automatically
        compressed into an archive of that filetype.
    is_overwritable : optional[bool]
        ``True`` if overwriting this file when this file already exists *or*
        ``False`` if raising an exception when this file already exists.
        Defaults to `False` for safety.
    '''

    # If only one object is passed, save only that object rather than the
    # 1-tuple consisting only of that object.
    if len(objs) == 1:
        objs = objs[0]

    # Save these objects to this file, silently compressing this file if this
    # filename is suffixed by an archive filetype.
    with iofiles.writing_bytes(
        filename=filename, is_overwritable=is_overwritable) as pickle_file:
        dill.dump(
            objs,
            file=pickle_file,
            protocol=PROTOCOL,

            #FIXME: We may need to increase the maximum recursion depth *BEFORE*
            #calling dill.dump() above. Let's ignore this until issues arise.

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

# ....................{ INITIALIZERS                      }....................
def init() -> None:
    '''
    Initialize both this submodule *and* the :mod:`dill` package.

    Specifically, this function instructs :mod:`dill` to pickle with our
    application-specific :class:`BetsePickler` subclass.
    '''

    # Log this initialization.
    logs.log_debug('Initializing dill...')

    # Core "dill" submodule. Sadly, dill 0.2.8.0 and 0.2.8.1 but (confusingly)
    # *NOT* dill >= 0.2.8.2 fundamentally broke backward compatibility by
    # renaming this public API -- necessitating conditional logic here to
    # properly find this submodule. For further details, see:
    #     https://github.com/uqfoundation/dill/issues/268
    dill_core_submodule = None

    # For dill >= 0.2.8.2, attempt to find the new "dill._dill" submodule.
    try:
        dill_core_submodule = dill._dill
    # For dill < 0.2.8.2, fallback to find the older "dill.dill" submodule.
    except:
        dill_core_submodule = dill.dill

    # If this submodule does *NOT* contain the expected "Pickler" base class,
    # raise an exception.
    objtest.die_unless_has_class(dill_core_submodule, 'Pickler')

    # Instruct "dill" to pickle with our application-specific pickler subclass.
    dill_core_submodule.Pickler = BetsePickler
