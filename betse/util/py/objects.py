#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level object facilities.
'''

# ....................{ IMPORTS                            }....................
# from betse.util.type import types

# ....................{ ITERATORS                          }....................
def iter_fields_nonbuiltin(obj: object):
    '''
    Generator yielding a 2-tuple of the name and value of each **non-builtin
    field** (i.e., variable with name _not_ both prefixed and suffixed by `__`)
    bound to the passed object, in lexicographically sorted field name order.

    Only fields registered in this object's internal dictionary (e.g.,
    `__dict__` in standard unslotted objects) will be yielded. Fields defined
    by this object's `__getattr__()` method or related runtime magic will _not_
    be yielded.
    '''
    # Note that:
    #
    # * Calling the dir() wrapper inspect.getmembers() here would allow to us to
    #   support edge-case attributes when passed class objects, including:
    #   * Metaclass attributes of the passed class.
    #   * Attributes decorated by "@DynamicClassAttribute" of the passed class.
    #   Since BETSE currently requires neither, we prefer calling the slightly
    #   lower-level but substantially faster dir() builtin.
    # * Calling vars() rather than dir() here would allow us to avoid calling
    #   getattr() and hence slightly increase time efficiency at a cost of
    #   failing for builtin containers (e.g., "dict", "list") *AND* copying
    #   object attribute values into a new dictionary. You do the ugly math.
    for attr_name in dir(obj):
        # If this attribute is *NOT* a builtin...
        if not (attr_name.startswith('__') and attr_name.endswith('__')):
            # ...and is a field, yield this field.
            attr_value = getattr(obj, attr_name)
            if not callable(attr_value):
                yield attr_name, attr_value
