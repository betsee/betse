#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level object facilities.
'''

# ....................{ IMPORTS                            }....................
import inspect
from betse.util.type.types import type_check, CallableTypes, GeneratorType

# ....................{ TESTERS                            }....................
@type_check
def is_method(obj: object, method_name: str) -> bool:
    '''
    `True` only if a method with the passed name is bound to the passed object.

    Parameters
    ----------
    obj : object
        Object to test for this method.
    method_name : str
        Name of the method to test this object for.

    Returns
    ----------
    bool
        `True` only if a method with this name is bound to this object.
    '''

    # Attribute with this name in this object if any or None otherwise.
    method = getattr(obj, method_name, None)

    # Return whether this attribute is a method.
    return callable(method)

# ....................{ GETTERS                            }....................
@type_check
def get_method_or_none(obj: object, method_name: str) -> CallableTypes:
    '''
    Method with the passed name bound to the passed object if any _or_ `None`
    otherwise.

    Parameters
    ----------
    obj : object
        Object to obtain this method from.
    method_name : str
        Name of the method to be obtained.

    Returns
    ----------
    callable, None
        Method with this name in this object if any _or_ `None` otherwise.
    '''

    # Attribute with this name in this object if any or None otherwise.
    method = getattr(obj, method_name, None)

    # If this attribute is a method, return this attribute; else, return None.
    return method if method is not None and callable(method) else None

# ....................{ ITERATORS                          }....................
def iter_fields_simple_custom(obj: object) -> GeneratorType:
    '''
    Generator yielding a 2-tuple of the name and value of each **non-builtin
    non-property field** (i.e., variable whose name is _not_ both prefixed and
    suffixed by `__` and whose value _not_ dynamically defined by the
    `@property` decorator  to be the implicit result of a method call) bound to
    the passed object, in lexicographically sorted field name order.

    Only fields registered in this object's internal dictionary (e.g.,
    `__dict__` in standard unslotted objects) will be yielded. Fields defined
    by this object's `__getattr__()` method or related runtime magic will _not_
    be yielded.

    Parameters
    ----------
    obj : object
        Object to yield the non-builtin fields of.

    Yields
    ----------
    (field_name, field_value)
        2-tuple of the name and value of each non-builtin field bound to this
        object, in lexicographically sorted field name order.
    '''

    # Ideally, this function would be internally reimplemented in terms of
    # the canonical inspect.getmembers() function. Dynamic inspection is
    # surprisingly non-trivial in the general case, particularly when virtual
    # base classes rear their diamond-studded faces. Moreover, doing so would
    # support edge-case attributes when passed class objects, including:
    #
    # * Metaclass attributes of the passed class.
    # * Attributes decorated by "@DynamicClassAttribute" of the passed class.
    #
    # Sadly, inspect.getmembers() internally accesses attributes via the
    # dangerous getattr() builtin rather than the safe inspect.getattr_static()
    # function. Since the codebase explicitly requires the latter, this function
    # reimplements rather than defers to inspect.getmembers(). (Sadness reigns.)
    #
    # For the same reason, the unsafe vars() builtin cannot be called either.
    # Since that builtin fails for builtin containers (e.g., "dict", "list"),
    # this is not altogether a bad thing.
    for attr_name in dir(obj):
        # If this attribute is *NOT* a builtin...
        if not (attr_name.startswith('__') and attr_name.endswith('__')):
            # Value of this attribute guaranteed to be statically rather than
            # dynamically retrieved. The getattr() builtin performs the latter,
            # dynamically calling this attribute's getter if this attribute is
            # a property. Since such call could conceivably raise unwanted
            # exceptions *AND* since this function explicitly ignores
            # properties, static attribute retrievable is strongly preferable.
            attr_value = inspect.getattr_static(obj, attr_name)

            # If this value is neither callable nor a property, this attribute
            # is a field. In this case, yield this field.
            if not (callable(attr_value) or isinstance(attr_value, property)):
                yield attr_name, attr_value

# ....................{ DECORATORS                         }....................
