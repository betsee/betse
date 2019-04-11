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
    BetseAttrException,
    BetseCallableException,
    BetseMethodException,
    BetseTypeException,
)
from betse.util.type import types
from betse.util.type.types import (
    type_check,
    CallableTypes,
    CallableOrNoneTypes,
    ClassType,
    TestableTypes,
    TestableOrNoneTypes,
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
def die_unless_has_attr(obj: object, *attr_names: str) -> None:
    '''
    Raise an exception unless the passed object defines all attributes with the
    passed names.

    Parameters
    ----------
    obj : object
        Object to be validated.
    attr_names : tuple[str]
        Tuple of the names of all attributes to validate this object to define.

    Raises
    ----------
    BetseAttrException
        If one or more such attributes are *not* bound to this object.
    '''

    # If this object fails to define one or more of these attributes...
    if not has_attr(obj, *attr_names):
        # For the name of each such attribute...
        for attr_name in attr_names:
            # If this object fails to define a attribute with this name, raise
            # an exception describing this failure. For simplicity, reuse the
            # existing exception raising implementation of get_attr().
            if not has_attr(obj, attr_name):
                get_attr(obj, attr_name)


@type_check
def die_unless_has_class(obj: object, *class_names: str) -> None:
    '''
    Raise an exception unless the passed object defines all classes with the
    passed names.

    Parameters
    ----------
    obj : object
        Object to be validated.
    class_names : tuple[str]
        Tuple of the names of all classes to validate this object to define.

    Raises
    ----------
    BetseTypeException
        If one or more such classes are *not* bound to this object.
    '''

    # If this object fails to define one or more of these classes...
    if not has_class(obj, *class_names):
        # For the name of each such class...
        for class_name in class_names:
            # If this object fails to define a class with this name, raise an
            # exception describing this failure.
            if not has_class(obj, class_name):
                raise BetseTypeException(
                    'Class "{}.{}" undefined.'.format(
                        get_class_name_unqualified(obj), class_name))


@type_check
def die_unless_has_method(obj: object, *method_names: str) -> None:
    '''
    Raise an exception unless the passed object defines all methods with the
    passed names.

    Parameters
    ----------
    obj : object
        Object to be validated.
    method_names : tuple[str]
        Tuple of the names of all methods to validate this object to define.

    Raises
    ----------
    BetseMethodException
        If one or more such methods are *not* bound to this object.
    '''

    # If this object fails to define one or more of these methods...
    if not has_method(obj, *method_names):
        # For the name of each such method...
        for method_name in method_names:
            # If this object fails to define a method with this name, raise an
            # exception describing this failure.
            if not has_method(obj, method_name):
                raise BetseMethodException(
                    'Method "{}.{}" undefined.'.format(
                        get_class_name_unqualified(obj), method_name))

# ....................{ TESTERS ~ attr                    }....................
@type_check
def has_attr(obj: object, *attr_names: str) -> bool:
    '''
    ``True`` only if the passed object defines all attributes with the passed
    names.

    Parameters
    ----------
    obj : object
        Object to be tested.
    attr_names : tuple[str]
        Tuple of the names of all attributes to test this object for.

    Returns
    ----------
    bool
        ``True`` only if this object defines all such attributes.
    '''

    # Return true only if, for each such attribute name, this object defines an
    # attribute with this name.
    return all(hasattr(obj, attr_name) for attr_name in attr_names)


@type_check
def has_class(obj: object, *class_names: str) -> bool:
    '''
    ``True`` only if the passed object defines all classes with the passed
    names.

    Parameters
    ----------
    obj : object
        Object to be tested.
    class_names : tuple[str]
        Tuple of the names of all classes to test this object for.

    Returns
    ----------
    bool
        ``True`` only if this object defines all such classes.
    '''

    # Return true only if...
    return all(
        # This object defines a class with this name... Dismantled, this is:
        #
        # * getattr(), returning an attribute with this name bound to this
        #   object if any *OR* "None" otherwise.
        # * is_class(), returning True only if this attribute is a class. Since
        #   "None" is guaranteed to *NOT* be a class, this test suffices.
        types.is_class(getattr(obj, class_name, None))
        # For each passed class name.
        for class_name in class_names)


@type_check
def has_method(obj: object, *method_names: str) -> bool:
    '''
    ``True`` only if the passed object defines all methods with the passed
    names.

    Parameters
    ----------
    obj : object
        Object to be tested.
    method_names : tuple[str]
        Tuple of the names of all methods to test this object for.

    Returns
    ----------
    bool
        ``True`` only if this object defines all such methods.
    '''

    # Return true only if...
    return all(
        # This object defines a method with this name...
        get_callable_or_none(obj=obj, callable_name=method_name) is not None
        # For each passed method name.
        for method_name in method_names)

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
    Magic number defined by the ``Include/object.h`` C header in the official
    CPython codebase.

    CPython ORs the ``__flags__`` bit field of the class object for each **heap
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
        # For unknown but probably banal reasons, the dir() builtin strips the
        # "__dict__" attribute name from its returned list only for C-based
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
def get_attr(obj: object, attr_name: str, **kwargs) -> object:
    '''
    Value of the attribute with the passed name bound to the passed object if
    this object defines this attribute *or* raise an exception otherwise (i.e.,
    if this object defines no such attribute), optionally validated to be of
    the passed type.

    Parameters
    ----------
    obj : object
        Object to be inspected.
    attr_name : str
        Name of the attribute to be retrieved.

    All remaining keyword arguments are passed as is to the
    :func:`get_attr_or_default` function.

    Returns
    ----------
    object
        Value of the attribute with this name bound to this object.

    Raises
    ----------
    BetseAttrException
        If no such attribute is bound to this object.
    BetseTypeException
        If the ``attr_type`` parameter is non-``None`` and the type of the
        current value of this attribute is *not* an instance of ``attr_type``.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.obj.sentinels import SENTINEL

    # Value of the attribute with this name defined by this object if any *OR*
    # the sentinel otherwise.
    attr_value = get_attr_or_sentinel(obj, attr_name, **kwargs)

    # If no such attribute exists, raise an exception.
    if attr_value is SENTINEL:
        raise BetseAttrException(
            'Attribute "{}.{}" undefined.'.format(
                get_class_name_unqualified(obj), attr_name))
    # Else, this attribute exists.

    # Return this value.
    return attr_value


@type_check
def get_attr_or_none(obj: object, attr_name: str, **kwargs) -> object:
    '''
    Attribute with the passed name bound to the passed object if this object
    defines this attribute *or* ``None`` otherwise (i.e., if this object
    defines no such attribute), optionally validated to be of
    the passed type.

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
        Object to be inspected.
    attr_name : str
        Name of the attribute to be obtained.

    All remaining keyword arguments are passed as is to the
    :func:`get_attr_or_default` function.

    Returns
    ----------
    object
        Attribute with this name bound to this object if any *or* ``None``
        otherwise.

    Raises
    ----------
    BetseTypeException
        If the ``attr_type`` parameter is non-``None`` and the type of the
        current value of this attribute is *not* an instance of ``attr_type``.
    '''

    # Return the current value of the attribute with this name defined by this
    # object if any *OR* "None" otherwise.
    return get_attr_or_default(
        obj=obj, attr_name=attr_name, attr_default=None, **kwargs)


@type_check
def get_attr_or_sentinel(obj: object, attr_name: str, **kwargs) -> object:
    '''
    Value of the attribute with the passed name bound to the passed object if
    this object defines this attribute *or* the sentinel singleton otherwise
    (i.e., if this object defines no such attribute), optionally validated to
    be of the passed type.

    This function enables callers to safely distinguish between non-existing
    attributes and existing attributes whose values are ``None``.

    Parameters
    ----------
    obj : object
        Object to be inspected.
    attr_name : str
        Name of the attribute to be obtained.

    All remaining keyword arguments are passed as is to the
    :func:`get_attr_or_default` function.

    Returns
    ----------
    object
        Either:

        * If this object declares this attribute, this attribute's value.
        * Else, the **sentinel singleton** (i.e.,
          :attr:`betse.util.type.obj.sentinels.SENTINEL`).

    Raises
    ----------
    BetseTypeException
        If the ``attr_type`` parameter is non-``None`` and the type of the
        current value of this attribute is *not* an instance of ``attr_type``.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.obj.sentinels import SENTINEL

    # Return the current value of the attribute with this name defined by this
    # object if any *OR* the sentinel otherwise.
    return get_attr_or_default(
        obj=obj, attr_name=attr_name, attr_default=SENTINEL, **kwargs)


@type_check
def get_attr_or_default(
    # Mandatory parameters.
    obj: object,
    attr_name: str,
    attr_default: object,

    # Optional parameters.
    attr_type: TestableOrNoneTypes = None,
) -> object:
    '''
    Value of the attribute with the passed name bound to the passed object if
    this object defines this attribute *or* the passed default value otherwise
    (i.e., if this object defines no such attribute), optionally validated to
    be of the passed type.

    Parameters
    ----------
    obj : object
        Object to be inspected.
    attr_name : str
        Name of the attribute to return the current value of.
    attr_default : object
        Default value to be returned if this object defines no such attribute.
    attr_type : TestableOrNoneTypes
        Expected type of the current value of this attribute. This function
        effectively performs the equivalent of the :meth:`type_check` decorator
        at runtime by raising an exception if all of the following apply:

        * This type is *not* ``None``.
        * This value is *not* this default value, implying this attribute to be
          defined by this object.
        * This value is *not* an instance of this type.

        Defaults to ``None``, in which case no such type checking is performed.

    Returns
    ----------
    object
        Either:

        * If this object defines this attribute, this attribute's value.
        * Else, this default value.

    Raises
    ----------
    BetseTypeException
        If the ``attr_type`` parameter is non-``None`` and the type of the
        current value of this attribute is *not* an instance of ``attr_type``.
    '''

    # Value of the attribute with this name defined by this object if any *OR*
    # this default value otherwise.
    attr_value = getattr(obj, attr_name, attr_default)

    # If this value is to be type-checked *AND* is *NOT* this default value
    # (which by definition already satisfies caller requirements regardless of
    # type), type-check this value.
    if attr_type is not None and attr_value is not attr_default:
        die_unless_instance(obj=attr_value, cls=attr_type)

    # Return this value.
    return attr_value

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

# ....................{ GETTERS : callable                }....................
@type_check
def get_callable(obj: object, callable_name: str) -> CallableTypes:
    '''
    Callable with the passed name defined by the passed object if this object
    defines such a callable *or* raise an exception otherwise (i.e., if this
    object defines no such callable).

    Parameters
    ----------
    obj : object
        Object to be queried.
    callable_name : str
        Name of the callable to be returned.

    Returns
    ----------
    CallableTypes
        Callable with this name defined by this object.

    Raises
    ----------
    BetseCallableException
        If this object defines no callable with this name.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.obj.sentinels import SENTINEL

    # Attribute with this name defined by this object if any *OR* the sentinel.
    # To raise human-readable exceptions, this attribute is *NOT* retrieved via
    # the higher-level get_callable_or_none() method; doing so would obscure
    # whether this attribute does not exist or does but is not callable.
    func = get_attr_or_sentinel(obj, callable_name)

    # If this attribute does *NOT* exist, raise an exception.
    if func is SENTINEL:
        raise BetseCallableException(
            'Method {}.{}() undefined.'.format(
                obj.__class__.__name__, callable_name))

    # If this attribute exists but is *NOT* callable, raise an exception.
    if not callable(func):
        raise BetseCallableException(
            'Object attribute "{}.{}" not callable: {!r}'.format(
                obj.__class__.__name__, callable_name, func))
    # Else, this attribute is callable.

    # Return this callable.
    return func


@type_check
def get_callable_or_none(obj: object, callable_name: str) -> (
    CallableOrNoneTypes):
    '''
    Callable with the passed name defined by the passed object if this object
    defines such a callable *or* ``None`` otherwise (i.e., if this object
    defines no such callable).

    Parameters
    ----------
    obj : object
        Object to be queried.
    callable_name : str
        Name of the callable to be returned.

    Returns
    ----------
    CallableOrNoneTypes
        Either:

        * If this object defines this callable, this callable.
        * Else, ``None``.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.obj.sentinels import SENTINEL

    # Attribute with this name defined by this object if any *OR* the sentinel.
    func = get_attr_or_sentinel(obj, callable_name)

    # Return this attribute if this attribute is callable *OR* "None".
    return func if func is not SENTINEL and callable(func) else None

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
