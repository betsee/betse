#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level **decorator** (i.e., classes and callables dynamically wrapping other
classes and callables at runtime) facilities.
'''

# ....................{ IMPORTS                            }....................
from abc import ABCMeta  #, abstractmethod
from betse.util.type.types import (
    type_check, CallableTypes, ClassType, MethodType)

# ....................{ DECORATORS                         }....................
class MethodDecorator(object, metaclass=ABCMeta):
    '''
    Abstract base class of all **method decorators** (i.e., decorators
    decorating *only* methods bound to class instances).

    This superclass efficiently caches bound methods on behalf of subclasses,
    guaranteeing all subclasses to be efficiently callable as proper methods.

    Attributes
    ----------
    _method: CallableTypes
        Unbound method (i.e., function) to be decorated.
    _obj_id_to_method_bound : dict
        Dictionary mapping from the unique identifier associated with each
        class instance containing a method decorated by this subclass to the
        same method bound to that class instance. While technically optional,
        the cache implemented by this dictionary avoids the need to recreate
        bound methods on each call to the :meth:`__get__` method.
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(self, method: CallableTypes) -> None:
        '''
        Initialize this method decorator.

        Parameters
        ----------
        method: CallableTypes
            Unbound method (i.e., function) to be decorated.
        '''

        # Classify all passed parameters.
        self._method = method

        # Initialize all remaining instance variables.
        self._obj_id_to_method_bound = {}

    # ..................{ DESCRIPTORS                        }..................
    def __get__(self, obj: object, cls: ClassType) -> MethodType:
        '''
        Create, cache, and return a decorated method bound to the passed object.
        '''

        # If this descriptor is accessed as a class rather than instance
        # variable, return this low-level descriptor rather than the high-level
        # value to which this expression evaluates.
        if obj is None:
            return self

        # Unique identifier associated with this object.
        obj_id = id(obj)

        # Attempt to...
        try:
            # Return the previously bound decorated method.
            return self._obj_id_to_method_bound[obj_id]
        except KeyError:
            #FIXME: If the line below fails, consider this common alternative:
            #    ... = functools.partial(self.__call__, obj)
            method_bound = self._obj_id_to_method_bound[obj_id] = MethodType(
                self.__call__, obj)
            return method_bound

    # ..................{ CALLERS                            }..................
    def __call__(self, obj, *args, **kwargs) -> object:
        '''
        Call the decorated method previously passed to the :meth:`__init__`
        method bound to the passed object with the passed positional and keyword
        arguments, returning the value returned by this call.

        This special method is typically overriden by subclass implementations
        wrapping the decorated method with additional functionality.
        '''

        return self._method(obj, *args, **kwargs)
