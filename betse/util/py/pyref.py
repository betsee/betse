#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level **reference** (i.e., pointer to the address of an object) facilities.
'''

# ....................{ IMPORTS                           }....................
import weakref
from betse.util.type.types import WeakRefProxyTypes, WeakRefType

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
    unsafely returns ``None`` rather than safely raising an exception after the
    passed object is collected and hence deleted.

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
def refer_weak(obj: object) -> WeakRefType:
    '''
    Weak reference to the passed object as an unproxied callable object, such
    that calling the latter returns a strong (i.e., non-weak and hence
    standard) reference to the former if the former has yet to be collected
    *or* ``None`` otherwise (i.e., if the former has already been collected).

    If this object is already a weak reference, this object is returned as is.

    Caveats
    ----------
    This function should *always* be called in lieu of the lower-level
    :func:`weakref.ref` function, which unsafely raises the following
    exception when passed an existing weak reference:

        TypeError: cannot create weak reference to 'weakref' object

    This function should typically (but *not* necessarily always) be avoided in
    favour of the higher-level :func:`proxy_weak` function, which safely raises
    an exception rather than silently and hence unsafely returning ``None``
    after the passed object is collected and hence deleted.

    Parameters
    ----------
    obj : object
        Object to be weakly refered to.

    Returns
    ----------
    WeakRefType
        Weak reference to the passed object as an unproxied callable object. When the passed
        object is garbage-collected, calling this weak reference returns
        ``None`` rather than the passed object. This weak reference is
        unhashable regardless of the hashability of the passed object, avoiding
        issues arising from the mutability of weak references and preventing
        their use as dictionary keys.
    '''

    # Return...
    return (
        # If this object is already a weak reference, this object as is.
        obj if isinstance(obj, WeakRefType) else
        # Else, a weak reference referring to this object.
        weakref.ref(obj)
    )
