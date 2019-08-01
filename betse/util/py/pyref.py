#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level **reference** (i.e., pointer to the address of an object) facilities.
'''

# ....................{ IMPORTS                           }....................
import weakref
from betse.util.type.types import WeakRefProxyTypes, WeakRefTypes

# ....................{ PROXIERS                          }....................
def proxy_weak(obj: object) -> WeakRefProxyTypes:
    '''
    Weak reference to the passed object as a proxy object, such that the latter
    transparently masquerades as the former.

    The reference returned by this function raises a :data:`ReferenceError`
    exception on each attempt to access this reference *after* Python collects
    and hence deletes the passed object. If this object is:

    * Already a proxy object, this object is returned as is.
    * Callable, a :class:`CallableProxyType` object proxying this object is
      returned.
    * Uncallable, a :class:`ProxyType` object proxying this object is returned.

    Caveats
    ----------
    This function should *always* be called in lieu of the lower-level
    :func:`weakref.proxy` function, which unsafely raises the following
    exception when passed an existing proxy object:

        TypeError: cannot create weak reference to 'weakproxy' object

    This function should typically (but *not* necessarily always) be called in
    lieu of the lower-level :func:`get_weak` function, which silently and hence
    unsafely returns ``None`` rather than safely raising an exception after
    Python collects and hence deletes the passed object.

    Parameters
    ----------
    obj : object
        Object to be weakly proxied.

    Returns
    ----------
    WeakRefProxyTypes
        Weak reference to the passed object as a proxy object. When the passed
        object is garbage-collected, this proxy object raises a
        :data:`ReferenceError` exception rather than yielding the passed
        object. This proxy object is unhashable regardless of the hashability
        of the passed object, avoiding issues arising from the mutability of
        weak references and preventing their use as dictionary keys.

    See Also
    ---------
    :func:`get_weak`
        Lower-level getter returning ``None`` rather than raising an exception.
    '''

    # Return...
    return (
        # If this object is already a weak proxy, this object as is.
        obj if isinstance(obj, WeakRefProxyTypes) else
        # Else, a weak reference proxying this object.
        weakref.proxy(obj)
    )

# ....................{ REFERERS                          }....................
def refer_weak(obj: object) -> WeakRefTypes:
    '''
    Weak reference to the passed object as an unproxied callable object, such
    that calling the latter returns a strong (i.e., non-weak and hence
    standard) reference to the former if the former has yet to be collected
    *or* ``None`` otherwise (i.e., if the former has already been collected).

    This function explicitly supports common edge-cases unsupported by the
    underlying :func:`weakref.ref` function. Notably, if the passed object is:

    * Already a **weak reference** (i.e., object whose type is in the
      :class:`WeakRefTypes` tuple), that object is returned as is.
    * A **bound method** (i.e., function whose first parameter is *always*
      implicitly passed the same parent object), that method is wrapped with a
      weak reference specific to bound methods.
    * Any other object, that object is wrapped with a general-purpose unproxied
      weak reference.

    Caveats
    ----------
    This function should typically (but *not* necessarily always) be avoided in
    favour of the higher-level :func:`proxy_weak` function, which safely raises
    an exception rather than silently and hence unsafely returning ``None``
    after Python collects and hence deletes the passed object.

    This function should *always* be called in lieu of the lower-level
    :func:`weakref.ref` function, which unsafely raises the following
    exception when passed an existing weak reference:

        TypeError: cannot create weak reference to 'weakref' object

    Parameters
    ----------
    obj : object
        Object to be weakly refered to.

    Returns
    ----------
    WeakRefTypes
        Weak reference to the passed object as an unproxied callable object.
        When the passed object is garbage-collected, calling this weak
        reference returns ``None`` rather than the passed object. This weak
        reference is unhashable regardless of the hashability of the passed
        object, avoiding issues arising from the mutability of weak references
        and preventing their use as dictionary keys.

    Examples
    ----------
    Calling the weak reference created and returned by this function yields a
    strong reference to the object initially passed to this function. To avoid
    subtle race conditions between pure-Python code leveraging this weak
    reference and the garbage collector collecting the underlying object, this
    reference should *only* be called in the following usage pattern:

        >>> from betse.util.io.log import logs
        >>> from betse.util.io.log.conf import logconf
        >>> from betse.util.py import pyref
        >>> log_config_weak = pyref.refer_weak(logconf.get_log_conf())
        >>> log_config_strong = log_config_weak()
        >>> if log_config_strong is not None:
        ...     logs.log_debug('Verbose logging enabled: %r',
        ...                    log_config_strong.is_verbose)
        ... else:
        ...     logs.log_debug('Logging configuration no longer exists.')

    That is to say, the following recommendations should be observed:

    * Whenever this weak reference is called, the resulting strong reference
      should be persisted to a variable for subsequent reference. Conversely,
      this weak reference itself should *never* be directly referenced, as
      doing so invites race conditions with the garbage collector.
    * The first access of this strong reference should be embedded in an if
      conditional testing whether the underlying object:

      * Is still alive, in which case this reference is non-``None``.
      * Has already been collected, in which case this reference is ``None``.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.types import BoundMethodTypes

    # Return...
    return (
        # If this object is already a weak reference, this object as is.
        obj if isinstance(obj, WeakRefTypes) else
        # Else if this object is a bound method, a weak reference specific
        # to bound methods (whose obscure nature warrants special handling).
        weakref.WeakMethod(obj) if isinstance(obj, BoundMethodTypes) else
        # Else, a general-purpose weak reference referring to this object.
        weakref.ref(obj)
    )
