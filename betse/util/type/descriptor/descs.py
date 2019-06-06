#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level **non-data descriptor** (i.e., objects satisfying the non-data
descriptor protocol, typically defined at class scope) facilities.
'''

# ....................{ IMPORTS                           }....................
from betse.util.type.types import (
    type_check,
    CallableTypes,
    ClassOrNoneTypes,
    ClassBoundMethodTypes,
)

# ....................{ CLASSES                           }....................
class ClassPropertyReadOnly(object):
    '''
    **Read-only class property** (i.e., property bound to a class rather than
    an instance of a class with no corresponding setter and hence read-only)
    descriptor.

    Decorators defined by this submodule (e.g., :func:`classproperty_readonly`)
    wrap their decorated callables with instances of this class.

    Caveats
    ----------
    While feasible, implementing writable class properties requires
    metaclasses, which invites subtle conflicts between uncooperative
    metaclasses in non-trivial class hierarchies. For maintainability, consider
    implementing writable class properties as customary class setters instead
    (e.g., ``def set_something(cls, something: object) -> None``).

    See also this `authoritative StackOverflow commentary`_ on the necessity of
    metaclasses in implementing writable class properties.

    .. _authoritative StackOverflow commentary:
       https://stackoverflow.com/questions/5189699/how-to-make-a-class-property

    Attributes
    ----------
    _getter : ClassBoundMethodTypes
        Class-bound getter method underlying this property.
    '''

    # ..................{ INITIALIZERS                      }..................
    @type_check
    def __init__(self, getter: ClassBoundMethodTypes) -> None:
        '''
        Initialize this read-only class property.

        Attributes
        ----------
        getter : ClassBoundMethodTypes
            Class-bound getter method underlying this property.
        '''

        # Classify this getter.
        self._getter = getter

    # ..................{ GETTERS                           }..................
    def __get__(self, obj: object, cls: ClassOrNoneTypes = None) -> object:
        '''
        Value of this read-only class property.
        '''

        # If this descriptor is accessed as an instance rather than class
        # variable, default to the class of this instance.
        if cls is None:
            cls = type(obj)

        # Return the result of calling the unbound function underlying the
        # underlying class-bound method descriptor with this object and class.
        return self._getter.__get__(obj, cls)()

# ....................{ DECORATORS                        }....................
@type_check
def classproperty_readonly(getter: CallableTypes) -> ClassPropertyReadOnly:
    '''
    **Read-only class property** (i.e., property bound to a class rather than
    an instance of a class with no corresponding setter and hence read-only),
    implemented by the passed unbound class method.

    Parameters
    ----------
    getter : CallableTypes
        **Unbound class method** (i.e., function whose first parameter is the
        subclass binding the property returned by this decorator).

    Returns
    ----------
    ClassPropertyReadOnly
        Read-only class property implemented by this unbound class method.

    See Also
    ----------
    :class:`ClassPropertyReadOnly`
        Further details.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.obj import objtest

    # If this getter is already a class-bound method wrapping an underlying
    # unbound function, raise an exception.
    objtest.die_if_instance(obj=getter, cls=ClassBoundMethodTypes)

    # Class method wrapping this getter.
    abstract_getter = classmethod(getter)

    # Return this method wrapped by a read-only class property descriptor.
    return ClassPropertyReadOnly(abstract_getter)


def abstractclassproperty_readonly(getter: CallableTypes) -> None:
    '''
    **Abstract read-only class property** (i.e., property bound to a class
    rather than an instance of a class with no corresponding setter and hence
    read-only, required to be reimplemented by subclasses), implemented by the
    passed unbound class method.

    See Also
    ----------
    :func:`classproperty_readonly`
        Further details.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.decorator.deccls import abstractclassmethod
    from betse.util.type.obj import objtest

    # If this getter is already a class-bound method wrapping an underlying
    # unbound function, raise an exception.
    objtest.die_if_instance(obj=getter, cls=ClassBoundMethodTypes)

    # Abstract class method wrapping this getter.
    abstract_getter = abstractclassmethod(getter)

    # Return this method wrapped by a read-only class property descriptor.
    return ClassPropertyReadOnly(abstract_getter)
