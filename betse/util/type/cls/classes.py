#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level class facilities.
'''

# ....................{ IMPORTS                            }....................
from betse.exceptions import BetseTypeException
from betse.util.type.types import (
    type_check,
    CallableTypes,
    ClassType,
    GeneratorType,
    MappingType,
    SequenceTypes,
)

# ....................{ EXCEPTIONS                         }....................
@type_check
def die_unless_subclass(subclass: ClassType, superclass: ClassType) -> None:
    '''
    Raise an exception unless the first passed class inherits and is thus a
    subclass of the second passed class.

    Parameters
    ----------
    subclass : ClassType
        Subclass to be validated.
    superclass : ClassType
        Superclass to be validated.
    '''

    if not issubclass(subclass, superclass):
        raise BetseTypeException(
            'Class {!r} not a subclass of class {!r}'.format(
                subclass, superclass))

# ....................{ GETTERS                            }....................
@type_check
def get_module_name(cls: ClassType) -> str:
    '''
    Fully-qualified name of the module defining the passed class.

    Parameters
    ----------
    cls : ClassType
        Class to retrieve this module name for.

    Returns
    ----------
    str
        Fully-qualified name of this module.

    Raises
    ----------
    BetseTypeException
        If this class has no ``__module__`` attribute, which should ideally
        *never* happen.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.obj import objects

    # Defer to this existing function, which suffices.
    return objects.get_class_module_name(cls)

# ....................{ GETTERS ~ name                     }....................
@type_check
def get_name(cls: ClassType) -> str:
    '''
    Unqualified name of the passed class.

    Parameters
    ----------
    cls : ClassType
        Class to retrieve this unqualified name for.

    Returns
    ----------
    str
        Unqualified name of this class.
    '''

    # Elegant simplicity diminishes aggressive tendencies.
    return cls.__name__


@type_check
def get_name_snakecase(cls: ClassType) -> str:
    '''
    Unqualified name of the passed class converted from CamelCase to snake_case
    (e.g., from ``MyClassName`` to ``my_class_name``).

    Parameters
    ----------
    cls : ClassType
        Class to retrieve this unqualified snake_case name for.

    Returns
    ----------
    str
        Unqualified snake_case name of this class.
    '''

    # Avoid circular import dependencies.
    from betse.util.py import pyident

    # Unqualified CamelCase name of this class.
    name_camelcase = get_name(cls)

    # Can it be? But it can.
    return pyident.convert_camelcase_to_snakecase(name_camelcase)

# ....................{ GETTERS ~ method                   }....................
@type_check
def get_method(cls: ClassType, method_name: str) -> CallableTypes:
    '''
    Possibly unbound method with the passed name defined by the passed class if
    any *or* raise an exception otherwise.

    If this method is:

    * A classmethod, this method is bound to this class.
    * A normal method, this method is unbound.
    * A static method, this method is unbound.

    Parameters
    ----------
    cls : ClassType
        Class to retrieve this method for.
    method_name : str
        Name of the method to be obtained.

    Returns
    ----------
    CallableTypes
        Possibly unbound method with this name.

    Raises
    ----------
    BetseMethodException
        If no such method is defined by this class.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.obj import objects

    # It's so simple, my eyes bleed.
    return objects.get_method(obj=cls, method_name=method_name)

# ....................{ ITERATORS                          }....................
@type_check
def iter_methods(cls: ClassType) -> GeneratorType:
    '''
    Generator yielding 2-tuples of the name and value of each method defined by
    the passed class (in ascending lexicographic order of method name).

    This generator *only* yields methods statically registered in the internal
    dictionary for this class (e.g., ``__dict__`` in unslotted classes),
    including:

    * Builtin methods, whose names are both prefixed and suffixed by ``__``.
    * Custom methods, whose names are *not* prefixed and suffixed by ``__``,
      including:
      * Custom standard methods.
      * Custom property methods (i.e., methods decorated by the builtin
        :func:`property` decorator).

    Parameters
    ----------
    cls : ClassType
        Class to iterate all methods of.

    Yields
    ----------
    (method_name, method_value)
        2-tuple of the name and value of each method bound to this object (in
        ascending lexicographic order of method name).
    '''

    # Avoid circular import dependencies.
    from betse.util.type.obj import objects

    # Well, isn't that special?
    yield from objects.iter_methods(obj=cls)


@type_check
def iter_methods_matching(
    cls: ClassType, predicate: CallableTypes) -> GeneratorType:
    '''
    Generator yielding 2-tuples of the name and value of each method defined by
    the passed class whose method name matches the passed predicate (in
    ascending lexicographic order of method name).

    Parameters
    ----------
    cls: ClassType
        Class to yield all matching methods of.
    predicate : CallableTypes
        Callable iteratively passed the name of each method bound to this
        object, returning ``True`` only if that name matches this predicate.

    Yields
    ----------
    (method_name, method_value)
        2-tuple of the name and value of each matching method bound to this
        object (in ascending lexicographic order of method name).

    See Also
    ----------
    :func:`iter_methods`
        Further details.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.obj import objects

    # Well, isn't that special?
    yield from objects.iter_methods_matching(obj=cls, predicate=predicate)

# ....................{ DEFINERS                           }....................
@type_check
def define_class(
    # Mandatory parameters.
    class_name: str,

    # Optional parameters.
    class_attr_name_to_value: MappingType = {},
    base_classes: SequenceTypes = (),
) -> ClassType:
    '''
    Dynamically define a new class with the passed name subclassing all passed
    base classes and providing all passed class attributes (e.g., class
    variables, methods).

    Parameters
    ----------
    class_name : str
        Name of the class to be created.
    class_attr_name_to_value : MappingType
        Mapping from the name to the initial value of each class attribute
        (e.g., class variable, method) to declare this class to contain.
        Defaults to the empty dictionary, equivalent to declaring a class with
        the trivial body ``pass``.
    base_classes : optional[SequenceTypes]
        Sequence of all base classes to subclass this class from. Defaults to
        the empty tuple, equivalent to the 1-tuple ``(object,)`` containing only
        the root base class of all classes.

    Returns
    ----------
    ClassType
        Class dynamically defined with this name from these base classes and
        class attributes.
    '''

    # Thank you, bizarre 3-parameter variant of the type.__init__() method.
    return type(class_name, base_classes, class_attr_name_to_value)
