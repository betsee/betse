#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level **data descriptor** (i.e., objects satisfying the data descriptor
protocol, typically defined at class scope) facilities.
'''

# ....................{ IMPORTS                            }....................
# from betse.util.type.types import type_check

# ....................{ CLASSES                            }....................
class DataDescriptorBound(object):
    '''
    High-level object encapsulating a low-level user-defined **data descriptor**
    (i.e., object satisfying the standard data descriptor protocol), bound at
    instantiation time to an instance of a user-defined class containing this
    descriptor and hence servicing exactly one instance.

    Motivation
    ----------
    This wrapper is principally intended to circumvent the conflict between
    immutable types in Python (e.g., :class:`int`, :class:`str`) and Python's
    pass-by-reference semantics.

    Consider two similar scenarios:

    * Passing an integer variable to a function as an integer parameter. In this
      case, modifying the latter does *not* modify the former.
    * Passing an instance of this class encapsulating an integer variable to a
      function as a parameter expecting only an instance of this class. In this
      case, modifying this instance actually modifies the integer variable
      encapsulated by this instance -- thereby circumventing this conflict.

    Design
    ----------
    Instances of this class access descriptors outside of class scope. Doing so
    is technically feasible but circumvents Python's default
    :meth:`object.__getattribute__` implementation transforming every attribute
    access ``b.x`` where ``x`` is a data descriptor declared at class scope into
    ``type(b).__dict__['x'].__get__(b, type(b))``. Circumventing this machinery
    typically defeats the purpose of instantiating a descriptor, but is
    warranted in this case.

    Why? Because this class neither requires nor desires this descriptor to ever
    be accessed via special machinery. Since the ``__get__()`` and ``__set__()``
    methods defined by this descriptor already implement all functionality
    required by this class, explicitly calling these methods as is suffices.

    Attributes
    ----------
    data_desc : object
        Low-level data descriptor encapsulated by this higher-level object.
    _obj : object
        Object to be temporarily "bound" to this data descriptor for the
        duration of each call to the :meth:`get` and :meth:`set` methods.

    See Also
    ----------
    :class:`DataDescriptorUnbound`
        Data descriptor *not* bound to its parent object at instantiation time.
    :func:`betse.util.descriptor.expralias.expr_alias`
        Lower-level function creating data descriptors encapsulated by this
        higher-level class, intended to be used where declaring a descriptor at
        class scope is appropriate (e.g., where the data descriptor to be
        aliased is known at class definition time).
    '''

    # ..................{ INITIALIZERS                       }..................
    def __init__(self, obj: object, data_desc: object) -> None:
        '''
        Encapsulate the passed data descriptor against the passed object.

        Parameters
        ----------
        obj : object
            Object to temporarily "bind" to this data descriptor for the
            duration of each call to the :meth:`get` and :meth:`set` methods.
        data_desc : None
            **Data descriptor** (i.e., object satisfying the standard data
            descriptor protocol) to be encapsulated.
        '''

        # Raise an exception unless this object is a data descriptor.
        die_unless_data_desc(data_desc)

        # Classify all passed parameters.
        self._obj = obj
        self.data_desc = data_desc

    # ..................{ GETTERS                            }..................
    def get(self) -> object:
        '''
        Value this data descriptor currently evaluates to.
        '''

        # Yes, the descriptor protocol explicitly supports direct method calls.
        # In particular, the "Invoking Descriptors" subsection of the official
        # "Descriptor HowTo Guide" notes:
        #
        #     A descriptor can be called directly by its method name. For
        #     example, d.__get__(obj).
        #
        # Since this data descriptor is accessed with respect to an instance
        # variable rather than a class, the second parameter passed to this
        # method call is "None".
        #
        # See: https://docs.python.org/3/howto/descriptor.html
        return self.data_desc.__get__(self._obj, None)

    # ..................{ SETTERS                            }..................
    def set(self, value: object) -> None:
        '''
        Set this data descriptor's current value to the passed value.
        '''

        self.data_desc.__set__(self._obj, value)


class DataDescriptorUnbound(object):
    '''
    High-level object encapsulating a low-level user-defined **data descriptor**
    (i.e., object satisfying the standard data descriptor protocol), *not* bound
    at instantiation time to an instance of a user-defined class containing this
    descriptor and hence servicing multiple different instances.

    Attributes
    ----------
    data_desc : object
        Low-level data descriptor encapsulated by this higher-level object.

    See Also
    ----------
    :class:`DataDescriptorBound`
        Data descriptor bound to its parent object at instantiation time.
    '''

    # ..................{ INITIALIZERS                       }..................
    def __init__(self, data_desc: object) -> None:
        '''
        Encapsulate the passed data descriptor.

        Parameters
        ----------
        data_desc : None
            **Data descriptor** (i.e., object satisfying the standard data
            descriptor protocol) to be encapsulated.
        '''

        # Raise an exception unless this object is a data descriptor.
        die_unless_data_desc(data_desc)

        # Classify all passed parameters.
        self.data_desc = data_desc

    # ..................{ GETTERS                            }..................
    def get(self, obj: object) -> object:
        '''
        Value this data descriptor currently evaluates to.

        Parameters
        ----------
        obj : object
            Object to be temporarily bound to this data descriptor for the
            duration of this call.

        Returns
        ----------
        object
            This data descriptor's current value.
        '''

        return self.data_desc.__get__(obj, None)

    # ..................{ SETTERS                            }..................
    def set(self, obj: object, value: object) -> None:
        '''
        Set this data descriptor's current value to the passed value.

        Parameters
        ----------
        obj : object
            Object to be temporarily bound to this data descriptor for the
            duration of this call.
        value : object
            Value to set this data descriptor to.
        '''

        self.data_desc.__set__(obj, value)

# ....................{ EXCEPTIONS                         }....................
def die_unless_data_desc(data_desc: object) -> None:
    '''
    Raise an exception unless the passed object is a **data descriptor** (i.e.,
    object defining both the special ``__get__`` and ``__set__`` methods and
    hence satisfying the standard data descriptor protocol).

    Parameters
    ----------
    data_desc : object
        Object to be tested.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.obj import objects

    # Raise an exception unless all special methods required by the standard
    # data descriptor protocol are bound to this object.
    objects.die_unless_method(data_desc, '__get__', '__set__')
