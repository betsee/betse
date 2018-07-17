#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level **reference** (i.e., pointer to the address of an object) facilities.
'''

# ....................{ IMPORTS                           }....................
import weakref
from betse.util.type.types import WeakRefProxyTypes

# ....................{ GETTERS                           }....................
def proxy_weak(obj: object) -> WeakRefProxyTypes:
    '''
    Weak reference to the passed object as a proxy object, such that the latter
    transparently masquerades as the former.

    The reference returned by this function raises a :data:`ReferenceError`
    exception when Python collects and hence deletes the passed object. If this
    object is:

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

    Parameters
    ----------
    obj : object
        Object to be proxied.

    Returns
    ----------
    WeakRefProxyTypes
        Weak reference to the passed object as a proxy object. When the passed
        object is garbage-collected, this proxy object raises a
        :data:`ReferenceError` exception rather than yielding the passed
        object. This proxy object is unhashable regardless of the hashability
        of the passed object, avoiding issues arising from the mutability of
        weak references and preventing their use as dictionary keys.
    '''

    # Return...
    return (
        # If this object is already a weak reference, this object as is.
        obj if isinstance(obj, WeakRefProxyTypes) else
        # Else, a weak reference proxying this object.
        weakref.proxy(obj)
    )
