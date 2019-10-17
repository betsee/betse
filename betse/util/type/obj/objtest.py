#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level **object validation** (i.e., functions enabling callers to test and
validate arbitrary objects to expose various attributes) facilities.
'''

# ....................{ IMPORTS                           }....................
import platform
from betse.exceptions import BetseMethodException, BetseTypeException
from betse.util.type import types
from betse.util.type.types import type_check, TestableTypes

# ....................{ EXCEPTIONS                        }....................
def die_unless_is(*objs: object) -> None:
    '''
    Raise an exception unless all passed objects are identical.

    Equivalently, this function raises an exception if any passed object is
    *not* identical to every other passed object.

    Parameters
    ----------
    objs : tuple[object]
        Tuple of all objects to be validated as identical.

    Raises
    ----------
    BetseTypeException
        If any passed object is *not* identical to every other passed object.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.iterable import iteriter

    # Iterable of all possible pairs of the passed objects.
    obj_pairs = iteriter.iter_combinations(iterable=objs, k=2)

    # For each possible pair of the passed objects...
    for obj_a, obj_b in obj_pairs:
        # If these two objects differ, raise an exception.
        if obj_a is not obj_b:
            raise BetseTypeException('{!r} not {!r}.'.format(obj_a, obj_b))

# ....................{ EXCEPTIONS ~ instance             }....................
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

    # If this object is an instance of this class, raise an exception.
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

    # If this object is *NOT* an instance of this class, raise an exception.
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

    # Avoid circular import dependencies.
    from betse.util.type.obj import objects

    # If this object fails to define one or more of these attributes...
    if not has_attr(obj, *attr_names):
        # For the name of each such attribute...
        for attr_name in attr_names:
            # If this object fails to define a attribute with this name, raise
            # an exception describing this failure. For simplicity, reuse the
            # existing exception raising implementation of get_attr().
            if not has_attr(obj, attr_name):
                objects.get_attr(obj, attr_name)


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

    # Avoid circular import dependencies.
    from betse.util.type.obj import objects

    # If this object fails to define one or more of these classes...
    if not has_class(obj, *class_names):
        # For the name of each such class...
        for class_name in class_names:
            # If this object fails to define a class with this name, raise an
            # exception describing this failure.
            if not has_class(obj, class_name):
                raise BetseTypeException(
                    'Class "{}.{}" undefined.'.format(
                        objects.get_class_name_unqualified(obj), class_name))


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

    # Avoid circular import dependencies.
    from betse.util.type.obj import objects

    # If this object fails to define one or more of these methods...
    if not has_method(obj, *method_names):
        # For the name of each such method...
        for method_name in method_names:
            # If this object fails to define a method with this name, raise an
            # exception describing this failure.
            if not has_method(obj, method_name):
                raise BetseMethodException(
                    'Method "{}.{}" undefined.'.format(
                        objects.get_class_name_unqualified(obj), method_name))

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

    # Avoid circular import dependencies.
    from betse.util.type.obj import objects

    # Return true only if...
    return all(
        # This object defines a method with this name...
        objects.get_callable_or_none(
            obj=obj, callable_name=method_name) is not None
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
