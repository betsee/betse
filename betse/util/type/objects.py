#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level object facilities.
'''

# ....................{ IMPORTS                            }....................
import inspect, sys
from betse.util.type.types import (
    type_check,
    CallableTypes,
    MethodType,
    NumericTypes,
    GeneratorType,
    SequenceTypes,
)
from functools import wraps
if False: wraps  # silence IDE warnings

# ....................{ DECORATORS                         }....................
def property_cached(property_method: MethodType) -> CallableTypes:
    '''
    Decorate the passed property method to cache the value returned by the first
    implicit call of this method.

    On the first access of a property decorated with this decorator, the passed
    method implementing this property is called, the value returned by this
    property internally cached into a private attribute of this method, and this
    value returned. On all subsequent accesses of this property, the cached
    value is returned as is _without_ calling this method. Hence, this method is
    called at most once for each instance of the class containing this property.
    '''

    # Raw string of Python statements comprising the body of this wrapper,
    # including (in order):
    #
    # * A "@wraps" decorator propagating the name, docstring, and other
    #   identifying metadata of the original function to this wrapper.
    # * A private "__beartype_func" parameter initialized to this function.
    #   In theory, the "func" parameter passed to this decorator should be
    #   accessible as a closure-style local in this wrapper. For unknown
    #   reasons (presumably, a subtle bug in the exec() builtin), this is
    #   not the case. Instead, a closure-style local must be simulated by
    #   passing the "func" parameter to this function at function
    #   definition time as the default value of an arbitrary parameter. To
    #   ensure this default is *NOT* overwritten by a function accepting a
    #   parameter of the same name, this edge case is tested for below.
    #
    # While there exist numerous alternative implementations for caching
    # properties, the approach implemented below has been profiled to be the
    # most efficient. Alternatives include (in order of decreasing efficiency):
    #
    # * Dynamically getting and setting a property-specific key-value pair of
    #   the internal dictionary for the current object, timed to be
    #   approximately 1.5 times as slow as exception handling: e.g.,
    #
    #     if not {property_name!r} in self.__dict__:
    #         self.__dict__[{property_name!r}] = __property_method(self)
    #     return self.__dict__[{property_name!r}]
    #
    # * Dynamically getting and setting a property-specific attribute of
    #   the current object: e.g.,
    #   the internal dictionary for the current object, timed to be
    #   approximately 1.5 times as slow as exception handling: e.g.,
    #
    #     if not hasattr(self, {property_name!r}):
    #         setattr(self, {property_name!r}, __property_method(self))
    #     return getattr(self, {property_name!r})
    func_body = '''
@wraps(__property_method)
def property_method_cached(self, __property_method=__property_method):
    try:
        return self.{property_name}
    except AttributeError:
        self.{property_name} = __property_method(self)
        return self.{property_name}
'''.format(property_name='__property_cached_' + property_method.__name__)

    # Dictionary mapping from local attribute name to value. For efficiency,
    # only attributes required by the body of this wrapper are copied from the
    # current namespace. (See below.)
    local_attrs = {'__property_method': property_method}

    # Dynamically define this wrapper as a closure of this decorator. For
    # obscure and presumably uninteresting reasons, Python fails to locally
    # declare this closure when the locals() dictionary is passed; to capture
    # this closure, a local dictionary must be passed instead.
    exec(func_body, globals(), local_attrs)

    # Return this wrapper method wrapped by a property descriptor.
    return property(local_attrs['property_method_cached'])

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

# ....................{ PRINTERS                           }....................
#FIXME: Document us up.
def print_vars_size(obj: object) -> None:

    # Avoid circular import dependencies.
    from betse.util.type.ints import MiB

    for attr_name, attr_size in iter_vars_size(obj, size_divisor=MiB):
        print('{}: {:.02f} MiB'.format(attr_name, attr_size))

# ....................{ ITERATORS                          }....................
def iter_vars(obj: object) -> SequenceTypes:
    '''
    Sequence of 2-tuples of the name and value of each variable bound to the
    passed object (_in ascending lexicographic order of variable name_).

    This includes:

    * All variables statically registered in this object's internal dictionary
      (e.g., `__dict__` in unslotted objects), including both:
      * Builtin variables, whose names are both prefixed and suffixed by `__`.
      * Custom variables, whose names are _not_ prefixed and suffixed by `__`.
    * All variables dynamically defined by the :func:`property` decorator to be
      the implicit result of a method call.

    This excludes _only_ variables dynamically defined by this object's
    `__getattr__()` method or related runtime magic, which cannot (_by
    definition_) be inspected by this generator.

    Parameters
    ----------
    obj : object
        Object to iterate the variables of.

    Returns
    ----------
    SequenceTypes
        Each element of this sequence is a 2-tuple `(var_name, var_value)` of
        the name and value of each variable bound to this object (_in ascending
        lexicographic order of variable name_).

    Raises
    ----------
    Exception
        If any variable of this object is a property whose
        :func:`property`-decorated method raises an exception. To avoid handling
        such exceptions, consider instead calling the
        :func:`iter_vars_simple_custom` generator excluding all properties.
    '''

    # The buck stops here.
    return inspect.getmembers(obj)


def iter_vars_simple(obj: object) -> GeneratorType:
    '''
    Generator yielding 2-tuples of the name and value of each **non-builtin
    variable** (i.e., variable whose name is _not_ both prefixed and suffixed by
    `__`) bound to the passed object (_in ascending lexicographic order of
    variable name_).

    Parameters
    ----------
    obj : object
        Object to yield all non-builtin variables of.

    Yields
    ----------
    (var_name, var_value)
        2-tuple of the name and value of each non-builtin variable bound to this
        object (_in ascending lexicographic order of variable name_).

    See Also
    ----------
    :func:`iter_vars`
        Further details.
    '''

    # For the name and value of each variable of this object...
    for var_name, var_value in iter_vars(obj):
        # If this variable is *NOT* builtin, yield this name and value.
        if not (var_name.startswith('__') and var_name.endswith('__')):
            yield var_name, var_value


#FIXME: Reimplement in terms of .
def iter_vars_simple_custom(obj: object) -> GeneratorType:
    '''
    Generator yielding 2-tuples of the name and value of each **non-builtin
    non-property variable** (i.e., variable whose name is _not_ both prefixed
    and suffixed by `__` and whose value _not_ dynamically defined by the
    `@property` decorator to be the implicit result of a method call) bound to
    the passed object (_in ascending lexicographic order of variable name_).

    Only variables statically registered in this object's internal dictionary
    (e.g., `__dict__` in unslotted objects) are yielded. Variables dynamically
    defined by this object's `__getattr__()` method or related runtime magic
    are simply ignored.

    Parameters
    ----------
    obj : object
        Object to yield all non-builtin non-property variables of.

    Yields
    ----------
    (var_name, var_value)
        2-tuple of the name and value of each non-builtin non-property variable
        bound to this object (_in ascending lexicographic order of variable
        name_).
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
    for var_name in dir(obj):
        # If this attribute is *NOT* a builtin...
        if not (var_name.startswith('__') and var_name.endswith('__')):
            # Value of this attribute guaranteed to be statically rather than
            # dynamically retrieved. The getattr() builtin performs the latter,
            # dynamically calling this attribute's getter if this attribute is
            # a property. Since such call could conceivably raise unwanted
            # exceptions *AND* since this function explicitly ignores
            # properties, static attribute retrievable is strongly preferable.
            var_value = inspect.getattr_static(obj, var_name)

            # If this value is neither callable nor a property, this attribute
            # is a non-property variable.
            if not (callable(var_value) or isinstance(var_value, property)):
                yield var_name, var_value


#FIXME: *UGH*. The sys.getsizeof() function only returns the shallow size of the
#passed object (i.e., the size of this object excluding the size of other
#objects referred to by this object). To retrieve the deep size of the passed
#object, we sadly need to implement a stack-based depth-first search (DFS) of
#some sort. Yes, we are not kidding.
#FIXME: Two general-purpose third-party solutions also exist: pympler and guppy.
#Pympler appears to be the more appropriate for our purposes. Given that,
#perhaps we could implement a solution resembling:
#
#* If pympler is available, defer to that.
#* Else, fallback to a DFS implemented by hand. Since this is *ONLY* a fallback,
#  this DFS should be implemented as trivially as possible.

#FIXME: Call this method during profiling with our three principal objects:
#"Simulator", "Cells", and "Parameters".
#FIXME: Call this method as a new functional test detecting erroneously large
#attributes in these three principal objects. To do so, determine the largest
#attribute currently produced in these objects for the minified world
#environment used for all functional tests. Then raise an exception if any
#attribute of these objects is greater than or equal to 10 times this
#expected maximum.
def iter_vars_size(
    obj: object, size_divisor: NumericTypes = 1) -> SequenceTypes:
    '''
    Sequence of 2-tuples of the name and in-memory size of each variable bound
    to the passed object (_in descending order of size_).

    Parameters
    ----------
    obj : object
        Object to iterate the variable sizes of.
    size_divisor: optional[NumericTypes]
        Number (i.e., integer, float) by which to divide the size in bytes of
        each variable of this object, hence converting from sizes in bytes to
        sizes in a larger denomination (e.g., :data:`betse.util.type.ints.KiB`,
        yielding sizes in kibibytes). Defaults to 1, yielding sizes in bytes.

    Returns
    ----------
    SequenceTypes
        Each element of this sequence is a 2-tuple `(var_name, var_size)` of
        the name and in-memory size of each variable bound to this object (_in
        descending order of size_).

    Raises
    ----------
    TypeError
        If any variable of this object does _not_ define the `__sizeof__()`
        special method. While all pure-Python and builtin objects define this
        special method, objects defined by third-party extensions are _not_
        required to do so.
    Exception
        If any variable of this object is a property whose
        :func:`property`-decorated method raises an exception.

    See Also
    ----------
    :func:`iter_vars`
        Further details.
    '''

    # Avoid circular import dependencies.
    from betse.util.type import iterables

    # List of all 2-tuples to be returned.
    vars_size = []

    # For each name and value of this object (in arbitrary order)...
    for var_name, var_value in iter_vars(obj):
        #FIXME: Define a new get_size() function of this submodule internally
        #deferring to either Pympler's asizeof() or sys.getsizeof() as needed.

        # In-memory size of this object, converted from bytes to another
        # denomination by the passed divisor.
        var_size = sys.getsizeof(var_value) / size_divisor

        # Append a 2-tuple of this variable's name and size.
        vars_size.append((var_name, var_size))

    # Return this list sorted on the second element of each such 2-tuple (i.e.,
    # each variable's size in bytes) in descending order.
    return iterables.sort_by_index_descending(
        iterable=vars_size, subiterable_index=1)
