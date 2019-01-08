#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level object facilities.
'''

#FIXME: For maintainability, extract all "die_"- and "is_"-prefixed functions
#into a newly defined "objtest" submodule of this subpackage.

# ....................{ IMPORTS                           }....................
import platform
from betse.exceptions import (
    BetseAttributeException, BetseMethodException, BetseTypeException)
from betse.util.type import types
from betse.util.type.types import (
    type_check,
    CallableTypes,
    CallableOrNoneTypes,
    ClassType,
    ModuleType,
    TestableTypes,
)

# ....................{ EXCEPTIONS                        }....................
@type_check
def die_if_instance(obj: object, cls: TestableTypes) -> None:
    '''
    Raise an exception if the passed object is an instance of either the passed
    class or one or more classes of the passed tuple of classes.

    Parameters
    ----------
    obj : object
        Object to be validated as instance of this class or tuple of classes.
    cls : TestableTypes
        Class or tuple of classes to validate this object to be an instance of.

    Raises
    ----------
    BetseTypeException
        If this object is an instance of either this class or one or more
        classes of this tuple of classes.
    '''

    if isinstance(obj, cls):
        raise BetseTypeException(
            'Object {!r} instances {!r}.'.format(obj, cls))


@type_check
def die_unless_instance(obj: object, cls: TestableTypes) -> None:
    '''
    Raise an exception unless the passed object is an instance of either the
    passed class or one or more classes of the passed tuple of classes.

    Parameters
    ----------
    obj : object
        Object to be validated as instance of this class or tuple of classes.
    cls : TestableTypes
        Class or tuple of classes to validate this object to be an instance of.

    Raises
    ----------
    BetseTypeException
        If this object is an instance of neither this class nor one or more
        classes of this tuple of classes.
    '''

    if not isinstance(obj, cls):
        raise BetseTypeException(
            'Object {!r} not an instance of {!r}.'.format(obj, cls))

# ....................{ EXCEPTIONS ~ attr                 }....................
@type_check
def die_unless_class(obj: object, *class_names: str) -> None:
    '''
    Raise an exception unless the passed object provides all classes with the
    passed names.

    Parameters
    ----------
    obj : object
        Object to test for these classes.
    class_names : tuple[str]
        Tuple of the names of all classes to test this object for.

    Raises
    ----------
    BetseMethodException
        If one or more such classes are *not* bound to this object.
    '''

    for class_name in class_names:
        if not is_class(obj=obj, class_name=class_name):
            raise BetseMethodException(
                'Object "{}" class "{}" undefined.'.format(obj, class_name))


@type_check
def die_unless_method(obj: object, *method_names: str) -> None:
    '''
    Raise an exception unless the passed object provides all methods with the
    passed names.

    Parameters
    ----------
    obj : object
        Object to test for these methods.
    method_names : tuple[str]
        Tuple of the names of all methods to test this object for.

    Raises
    ----------
    BetseMethodException
        If one or more such methods are *not* bound to this object.
    '''

    for method_name in method_names:
        if not is_method(obj=obj, method_name=method_name):
            raise BetseMethodException(
                'Object "{}" method {}() undefined.'.format(obj, method_name))

# ....................{ TESTERS                           }....................
@type_check
def is_attr(obj: object, attr_name: str) -> bool:
    '''
    ``True`` only if the passed object provides an attribute of arbitrary type
    with the passed name.

    This function is a convenient synonym of the :func:`hasattr` builtin,
    provided entirely as a caller convenience.

    Parameters
    ----------
    obj : object
        Object to be tested for this attribute.
    attr_name : str
        Name of the attribute to be tested for.

    Returns
    ----------
    bool
        ``True`` only if this attribute is bound to this object.
    '''

    return hasattr(obj, attr_name)


@type_check
def is_class(obj: object, class_name: str) -> bool:
    '''
    ``True`` only if the passed object defines a class with the passed name.

    Parameters
    ----------
    obj : object
        Object to test for this class.
    class_name : str
        Name of the class to test this object for.

    Returns
    ----------
    bool
        ``True`` only if a class with this name is bound to this object.
    '''

    # Attribute with this name in this object if any or "None" otherwise.
    cls = getattr(obj, class_name, None)

    # Return whether this attribute is a class. Since "None" is guaranteed to
    # *NOT* be a class, this simple test produces the expected result in all
    # possible cases.
    return types.is_class(cls)


@type_check
def is_method(obj: object, method_name: str) -> bool:
    '''
    ``True`` only if the passed object provides a method with the passed name.

    Parameters
    ----------
    obj : object
        Object to test for this method.
    method_name : str
        Name of the method to test this object for.

    Returns
    ----------
    bool
        ``True`` only if a method with this name is bound to this object.
    '''

    # Method with this name bound to this object if any or "None" otherwise.
    method = get_method_or_none(obj=obj, method_name=method_name)

    # Return true only if this method exists.
    return method is not None

# ....................{ TESTERS ~ pure                    }....................
# Tester distinguishing pure-Python from C-based class instances. Doing so is
# surprisingly non-trivial and, indeed, technically feasible with 100% accuracy
# *ONLY* under the official CPython interpreter. Why? Because the most accurate
# means of portably implementing this test in a cross-interpreter manner is to
# test whether the passed class instance defines the reserved "__dict__" or
# "__slots__" attributes. Whereas all pure-Python class instances are
# guaranteed to declare one or the other, most C-based class instances define
# neither. Of course, "most" is hardly "all." A C-based class may optionally
# request a "__dict__" attribute via "tp_dictoffset", resulting in false
# negatives in uncommon edge cases.
#
# If *NOT* running under the official CPython interpreter, this is the
# strongest possible implementation of this test; if running under the official
# CPython interpreter, however, this implementation may be strengthened by
# detecting whether or not the class of the passed class instance is a "heap
# type" (as defined by the "_Py_TPFLAGS_HEAPTYPE" docstring below). Only
# pure-Python classes are heap types, preserving a one-to-one correspondence
# between tester results and underlying interpreter reality.
#
# If the active Python interpreter is the official CPython implementation,
# prefer a more reliable CPython-specific solution guaranteed to succeed. Note
# that, to avoid circular import dependencies, this conditional unavoidably
# duplicates the existing betse.util.py.interpreters.is_cpython() tester.
if platform.python_implementation() == 'CPython':
    _Py_TPFLAGS_HEAPTYPE = (1<<9)
    '''
    Magic number defined by the `Include/object.h` C header in the official
    CPython codebase.

    CPython ORs the `__flags__` bit field of the class object for each **heap
    type** (i.e., pure-Python class whose type structure is dynamically
    allocated at interpreter runtime rather than a C-based class whose type
    structure is statically preallocated at either interpreter or C extension
    compile time).
    '''

    # To reduce duplication, this tester's docstring is dynamically set below.
    @type_check
    def is_pure_python(obj: object) -> bool:

        # If the passed object is *NOT* a class, this object *MUST* be a class
        # instance. (By Python 3.x design, all objects -- including primitive
        # builtins such as "int" and "bool" -- are one or the other). In this
        # case, obtain the class this instance is an instance of.
        cls = obj if isinstance(obj, type) else type(obj)

        # Return True only if this class is a heap type, as indicated by its
        # heap type bit flag being non-zero.
        return bool(cls.__flags__ & _Py_TPFLAGS_HEAPTYPE)
# Else, fallback to a CPython-agnostic solution typically but *NOT* necessarily
# succeeding. For all real-world objects of interest, this is effectively
# successful. Edge cases exist but are suitably rare. See the commentary above.
else:
    # To reduce duplication, this tester's docstring is dynamically set below.
    @type_check
    def is_pure_python(obj: object) -> bool:

        # If the passed object is a class, return True only if this class
        # defines either the "__dict__" or "__slots__" attributes.
        # Unfortunately, the trivial test "hasattr(cls, '__dict__')" does *NOT*
        # suffice to decide whether this class defines the "__dict__"
        # attribute.  Why? Because this attribute unconditionally exists for
        # *ALL* classes, pure-Python and C-based alike. Hence:
        #
        #     >>> hasattr(int, '__dict__')
        #     True
        #
        # For unknown and presumably banal reasons, the dir() builtin strips
        # the "__dict__" attribute name from its returned list only for C-based
        # classes. Hence, detecting whether a class is pure-Python or C-based
        # in a cross-interpreter manner reduces to iteratively searching the
        # list returned by dir() for this name.
        #
        # This constraint does *NOT* extend to the "__slots__" attribute, which
        # exists if and only if this class explicitly defines this attribute.
        # Are we having fun yet?
        if isinstance(obj, type):
            return '__dict__' in dir(obj) or hasattr(obj, '__slots__')
        # Else, this object *MUST* be a class instance. In this case, return
        # True only if similar conditions hold. Since a class instance defines
        # the "__dict__" attribute only if this instance is a pure-Python
        # instance whose class does *NOT* define the "__slots__" attribute,
        # this conditional may efficiently lookup the "__dict__" attribute
        # directly rather than inefficiently defer to the dir() builtin. Fun!
        else:
            return hasattr(obj, '__dict__') or hasattr(obj, '__slots__')

# Docstring dynamically set for the tester defined above.
is_pure_python.__doc__ = '''
``True`` if the passed object is either a pure-Python class or instance of such
a class *or* ``False`` otherwise (i.e., if this object is either a C-based
class or instance of such a class, either builtin or defined by a C
extension).

Parameters
----------
obj : object
    Object to be tested.

Returns
----------
bool
    ``True`` only if this object is a pure-Python class *or* instance of such a
    class.

See Also
----------
https://stackoverflow.com/a/41012823/2809027
    Stackoverflow answer strongly inspiring this implementation.
'''


def is_c_based(obj: object) -> bool:
    '''
    ``True`` if the passed object is either a C-based class or instance of such
    a class (either builtin or defined by a C extension) *or* ``False``
    otherwise (i.e., if this object is either a pure-Python class or instance
    of such a class).

    Parameters
    ----------
    obj : object
        Object to be tested.

    Returns
    ----------
    bool
        ``True`` only if this object is a C-based class *or* instance of such a
        class.
    '''

    return not is_pure_python(obj)

# ....................{ GETTERS : attr                    }....................
@type_check
def get_attr(obj: object, attr_name: str) -> object:
    '''
    Attribute with the passed name bound to the passed object if any *or* raise
    an exception otherwise.

    Parameters
    ----------
    obj : object
        Object to obtain this attribute from.
    attr_name : str
        Name of the attribute to be obtained.

    Returns
    ----------
    object
        Attribute with this name bound to this object.

    Raises
    ----------
    BetseAttributeException
        If no such attribute is bound to this object.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.obj.sentinels import SENTINEL

    # Attribute with this name in this object if any or the sentinel otherwise.
    attr = get_attr_or_sentinel(obj, attr_name)

    # If no such attribute exists, raise an exception.
    if attr is SENTINEL:
        raise BetseAttributeException(
            'Object attribute "{}.{}" undefined.'.format(
                obj.__class__.__name__, attr_name))

    # Else, this attribute exists. Return this attribute as is.
    return attr


@type_check
def get_attr_or_none(obj: object, attr_name: str) -> object:
    '''
    Attribute with the passed name bound to the passed object if any *or*
    ``None`` otherwise.

    Caveats
    ----------
    For disambiguity, consider calling the :func:`get_attr_or_sentinel`
    function instead. Whereas this function fails to distinguish between
    existing attributes whose values are ``None`` and non-existing attributes
    for which this function returns ``None``, the :func:`get_attr_or_sentinel`
    function trivially disambiguates between these two common edge cases.

    Parameters
    ----------
    obj : object
        Object to obtain this attribute from.
    attr_name : str
        Name of the attribute to be obtained.

    Returns
    ----------
    object
        Attribute with this name bound to this object if any *or* ``None``
        otherwise.
    '''

    # All that which glitters is not gold.
    return getattr(obj, attr_name, None)


@type_check
def get_attr_or_sentinel(obj: object, attr_name: str) -> object:
    '''
    Value of the attribute with the passed name bound to the passed object if
    any *or* :attr:`betse.util.type.obj.sentinels.SENTINEL` otherwise,
    permitting callers to distinguish between attributes that do *not* exist
    and attributes that do exist but whose values are ``None``.

    Parameters
    ----------
    obj : object
        Object to obtain this attr from.
    attr_name : str
        Name of the attribute to be obtained.

    Returns
    ----------
    object
        Value of this attribute if any *or* the sentinel singleton otherwise.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.obj.sentinels import SENTINEL

    # Return this attribute if any or the sentinel otherwise.
    return getattr(obj, attr_name, SENTINEL)

# ....................{ GETTERS ~ class                   }....................
def get_class(obj: object) -> ClassType:
    '''
    Passed object if this object is itself a class *or* the class of this
    object otherwise (i.e., if this object is *not* a class).

    Parameters
    ----------
    obj : object
        Object to retrieve this class for.

    Returns
    ----------
    ClassType
        This object if this object is a class *or* this object's class.
    '''

    # Simplicity is not a place in Simple City.
    return obj if isinstance(obj, ClassType) else type(obj)

# ....................{ GETTERS ~ class : name            }....................
def get_class_name_unqualified(obj: object) -> str:
    '''
    Unqualified name of either the passed object if this object is itself a
    class *or* the class of this object otherwise (i.e., if this object is
    *not* a class).

    Parameters
    ----------
    obj : object
        Object to retrieve this class name for.

    Returns
    ----------
    str
        Unqualified name of this class.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.cls import classes

    # This object if this object is a class *OR* this object's class otherwise.
    cls = get_class(obj)

    # Return the unqualified name of this class.
    return classes.get_name_unqualified(cls)


def get_class_module_name_qualified(obj: object) -> str:
    '''
    Fully-qualified name of the module defining either the passed object if
    this object is itself a class *or* the class of this object otherwise
    (i.e., if this object is *not* a class).

    Parameters
    ----------
    obj : object
        Object to retrieve this module name for.

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
    from betse.util.type.cls import classes

    # This object if this object is a class *OR* this object's class otherwise.
    cls = get_class(obj)

    # Return the fully-qualified name of this class.
    return classes.get_module_name_qualified(cls)

# ....................{ GETTERS : method                  }....................
@type_check
def get_method(obj: object, method_name: str) -> CallableTypes:
    '''
    Method with the passed name bound to the passed object if any *or* raise
    an exception otherwise.

    Parameters
    ----------
    obj : object
        Object to obtain this method from.
    method_name : str
        Name of the method to be obtained.

    Returns
    ----------
    CallableTypes
        Method with this name bound to this object.

    Raises
    ----------
    BetseMethodException
        If no such method is bound to this object.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.obj.sentinels import SENTINEL

    # Attribute with this name in this object if any or the sentinel otherwise.
    # To raise human-readable exceptions, this attribute is *NOT* retrieved via
    # the higher-level get_method_or_none() method; doing so would obscure
    # whether this attribute does not exist or does but is not a method.
    method = get_attr_or_sentinel(obj, method_name)

    # If no such attribute exists, raise an exception.
    if method is SENTINEL:
        raise BetseMethodException('Method {}.{}() undefined.'.format(
            obj.__class__.__name__, method_name))

    # If this attribute exists but is *NOT* a method, raise an exception.
    if not callable(method):
        raise BetseMethodException(
            'Object attribute "{}.{}" not a method: {!r}'.format(
            obj.__class__.__name__, method_name, method))

    # Else, this attribute is a method. Return this method as is.
    return method


@type_check
def get_method_or_none(obj: object, method_name: str) -> CallableOrNoneTypes:
    '''
    Method with the passed name bound to the passed object if any *or* ``None``
    otherwise.

    Parameters
    ----------
    obj : object
        Object to obtain this method from.
    method_name : str
        Name of the method to be obtained.

    Returns
    ----------
    CallableOrNoneTypes
        This method if any *or* ``None`` otherwise.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.obj.sentinels import SENTINEL

    # Attribute with this name in this object if any or the sentinel otherwise.
    method = get_attr_or_sentinel(obj, method_name)

    # If this attribute is a method, return this attribute; else, return
    # "None".
    return method if method is not SENTINEL and callable(method) else None
