#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level class facilities.
'''

# ....................{ IMPORTS                           }....................
import inspect
from betse.exceptions import BetseTypeException
from betse.util.type.types import (
    type_check,
    CallableTypes,
    ClassType,
    GeneratorType,
    MappingOrNoneTypes,
    SequenceOrNoneTypes,
    StrOrNoneTypes,
)

# ....................{ EXCEPTIONS                        }....................
@type_check
def die_unless_subclass(subclass: ClassType, superclass: ClassType) -> None:
    '''
    Raise an exception unless the first passed class inherits and is thus a
    subclass of the second passed class.

    For simplicity, note that any class is considered to be a subclass of
    itself within the context of this function. (See the example below.)

    Parameters
    ----------
    subclass : ClassType
        Subclass to be validated.
    superclass : ClassType
        Superclass to be validated.

    Raises
    ----------
    BetseTypeException
        If the first such class is *not* a subclass of the second such class.

    Examples
    ----------
        >>> from betse.util.type.cls import classes
        >>> class antiihnii(object): pass
        >>> class yee_naaldlooshii(antiihnii): pass
        >>> classes.die_unless_subclass(
        ...     subclass=antiihnii, superclass=antiihnii)
        >>> classes.die_unless_subclass(
        ...     subclass=yee_naaldlooshii, superclass=yee_naaldlooshii)
        >>> classes.die_unless_subclass(
        ...     subclass=yee_naaldlooshii, superclass=antiihnii)
        >>> classes.die_unless_subclass(
        ...     subclass=antiihnii, superclass=yee_naaldlooshii)
        Traceback (most recent call last):
          File "/home/leycec/tmp/yilo.py", line 12, in <module>
            subclass=antiihnii, superclass=yee_naaldlooshii)
          File "<string>", line 24, in __die_unless_subclass_type_checked__
          File "/home/leycec/py/betse/betse/util/type/cls/classes.py", line 63, in die_unless_subclass
            subclass, superclass))
        betse.exceptions.BetseTypeException: Class <class '__main__.antiihnii'> not a subclass of class <class '__main__.yee_naaldlooshii'>
    '''

    if not issubclass(subclass, superclass):
        raise BetseTypeException(
            'Class {!r} not a subclass of class {!r}'.format(
                subclass, superclass))

# ....................{ TESTERS                           }....................
@type_check
def is_abstract(cls: ClassType) -> bool:
    '''
    ``True`` only if the passed class is an **abstract base class** (i.e., a
    class whose metaclass is :class:`abc.ABCMeta` such that at least one
    abstract method decorated by the :func:`abc.abstractmethod` decorator
    remains unimplemented).

    Parameters
    ----------
    cls : ClassType
        Class to be tested.

    Returns
    ----------
    bool
        ``True`` only if this class is an abstract base class (ABC).
    '''

    # One-liners for ponderous victory.
    return inspect.isabstract(cls)

# ....................{ GETTERS                           }....................
@type_check
def get_module_name_qualified(cls: ClassType) -> str:
    '''
    Fully-qualified name of the module defining the passed class if this class
    is defined by a module *or* raise an exception otherwise.

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
        If this class is *not* defined by a module. Specifically, if either:

        * This class defines no special ``__module__`` attribute.
        * This class defines a special ``__module__`` attribute whose value is
          either:

          * ``None``.
          * A special identifier. Notably, builtin classes defined by no
            modules (e.g., :class:`str`) define a special ``__module__``
            attribute whose value is ``__builtin__``.

    See Also
    ----------
    https://stackoverflow.com/a/13653312/2809027
        StackOverflow answer strongly inspiring this implementation.
    '''

    # Avoid circular import dependencies.
    from betse.util.py import pyident

    # If this class does *NOT* define the special "__module__" attribute,
    # raise an exception. In theory, all classes including builtin classes
    # should define this. In practice, some do not. Thanks alot, everybody!
    if not hasattr(cls, '__module__'):
        raise BetseTypeException(
            'Class "{class_name}" module undefined '
            '(i.e., "{class_name}.__module__" attribute not found).'.format(
                class_name=get_name_unqualified(cls)))
    # Else, this class defines this attribute.

    # Fully-qualified name of the module defining this class if any *OR* "None"
    # otherwise.
    module_name = cls.__module__

    # If this name is either undefined *OR* prefixed and suffixed by both "__"
    # (e.g., as is the case with "str.__module__ == '__builtin__'"), no module
    # defines this class. In this case, raise an exception.
    if module_name is None or pyident.is_special(module_name):
        raise BetseTypeException(
            'Class "{class_name}" module undefined '
            '(i.e., "{class_name}.__module__" attribute is '
            '"{module_name}").'.format(
                class_name=get_name_unqualified(cls),
                module_name=str(module_name)))

    # Else, return this name.
    return module_name


@type_check
def get_module_name_qualified_or_none(cls: ClassType) -> StrOrNoneTypes:
    '''
    Fully-qualified name of the module defining the passed class if this class
    is defined by a module *or* ``None`` otherwise.

    Parameters
    ----------
    cls : ClassType
        Class to retrieve this module name for.

    Returns
    ----------
    StrOrNoneTypes
        Fully-qualified name of this module if any *or* ``None`` otherwise.
        Specifically, ``None`` is returned if either:

        * This class defines no special ``__module__`` attribute.
        * This class defines a special ``__module__`` attribute whose value is
          either:

          * ``None``.
          * A special identifier. Notably, builtin classes defined by no
            modules (e.g., :class:`str`) define a special ``__module__``
            attribute whose value is ``__builtin__``.

    See Also
    ----------
    https://stackoverflow.com/a/13653312/2809027
        StackOverflow answer strongly inspiring this implementation.
    '''

    # Avoid circular import dependencies.
    from betse.util.py import pyident

    # If this class does *NOT* define the special "__module__" attribute,
    # return "None". In theory, all classes including builtin classes should
    # define this. In practice, some do not. Thanks for nuthin', somebody!
    if not hasattr(cls, '__module__'):
        return None
    # Else, this class defines this attribute.

    # Fully-qualified name of the module defining this class if any *OR* "None"
    # otherwise.
    module_name = cls.__module__

    # Return "None" if this name is either undefined *OR* prefixed and suffixed
    # by both "__" (e.g., as is the case with "str.__module__ ==
    # '__builtin__'") or this name otherwise
    return (
        None
        if module_name is None or pyident.is_special(module_name) else
        module_name
    )

# ....................{ GETTERS ~ name                    }....................
@type_check
def get_name_qualified(cls: ClassType) -> str:
    '''
    Fully-qualified name of the passed class.

    Parameters
    ----------
    cls : ClassType
        Class to retrieve this fully-qualified name for.

    Returns
    ----------
    str
        Fully-qualified name of this class.
    '''

    # Unqualified name of this class.
    class_name = get_name_unqualified(cls)

    # Fully-qualified name of the module defining this class if this class is
    # defined by a module *OR* "None" otherwise.
    module_name = get_module_name_qualified_or_none(cls)

    # Return either...
    return (
        # The concatenation of this class and module name if this module name
        # is well-defined.
        '{}.{}'.format(module_name, class_name)
        if module_name is not None else
        # This class name as is otherwise.
        class_name
    )


@type_check
def get_name_unqualified(cls: ClassType) -> str:
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

# ....................{ GETTERS ~ method                  }....................
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
    return objects.get_callable(obj=cls, callable_name=method_name)

# ....................{ ITERATORS                         }....................
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
    (method_name : str, method : MethodType)
        2-tuple of the name and value of each method bound to this object (in
        ascending lexicographic order of method name).
    '''

    # Avoid circular import dependencies.
    from betse.util.type.obj import objiter

    # Well, isn't that special?
    yield from objiter.iter_methods(obj=cls)


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
        Callable iteratively passed the name of each method and that method
        bound to this object, returning ``True`` only if that name matches this
        predicate.

    Yields
    ----------
    (method_name : str, method : MethodType)
        2-tuple of the name and value of each matching method bound to this
        object (in ascending lexicographic order of method name).

    See Also
    ----------
    :func:`iter_methods`
        Further details.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.obj import objiter

    # Well, isn't that special?
    yield from objiter.iter_methods_matching(obj=cls, predicate=predicate)

# ....................{ DEFINERS                          }....................
@type_check
def define_class(
    # Mandatory parameters.
    class_name: str,

    # Optional parameters.
    class_attr_name_to_value: MappingOrNoneTypes = None,
    base_classes: SequenceOrNoneTypes = None,
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
        the empty tuple, equivalent to the 1-tuple ``(object,)`` containing
        only the root base class of all classes.

    Returns
    ----------
    ClassType
        Class dynamically defined with this name from these base classes and
        class attributes.
    '''

    if class_attr_name_to_value is None:
        class_attr_name_to_value = {}
    if base_classes is None:
        base_classes = ()

    # Thank you, bizarre 3-parameter variant of the type.__init__() method.
    return type(class_name, base_classes, class_attr_name_to_value)
