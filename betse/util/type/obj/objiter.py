#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level **object iteration** (i.e., functions enabling callers to iterate
over attributes on arbitrary objects, typically with generators) facilities.
'''

# ....................{ IMPORTS                           }....................
import inspect
from betse.util.type.types import (
    type_check, CallableTypes, GeneratorType, TestableTypes)

# ....................{ ITERATORS ~ attrs                 }....................
@type_check
def iter_attrs(obj: object) -> GeneratorType:
    '''
    Generator yielding 2-tuples of the name and value of each attribute bound
    to the passed object (in ascending lexicographic order of attribute name).

    Parameters
    ----------
    obj : object
        Object to yield all matching attributes of.

    Yields
    ----------
    (attr_name : str, attr_value : object)
        2-tuple of the name and value of each matching attribute bound to this
        object (in ascending lexicographic order of attribute name).
    '''

    # Defer to this generator.
    yield from iter_attrs_implicit_matching(
        obj=obj, predicate=noop_attr_predicate)


@type_check
def iter_attrs_instance_of(obj: object, cls: TestableTypes) -> GeneratorType:
    '''
    Generator yielding 2-tuples of the name and value of each attribute bound
    to the passed object whose value is an instance of the passed class or
    tuple of classes (in ascending lexicographic order of attribute name).

    Parameters
    ----------
    obj : object
        Object to yield all matching attributes of.
    cls : TestableTypes
        Class or tuple of classes to yield all instances of.

    Yields
    ----------
    (attr_name : str, attr_value : object)
        2-tuple of the name and value of each matching attribute bound to this
        object (in ascending lexicographic order of attribute name).
    '''

    # Defer to this generator.
    yield from iter_attrs_implicit_matching(
        obj=obj, predicate=lambda attr_name, attr_value: (
            isinstance(attr_value, cls)))

# ....................{ ITERATORS ~ attrs : matching      }....................
@type_check
def iter_attrs_implicit_matching(
    obj: object, predicate: CallableTypes) -> GeneratorType:
    '''
    Generator yielding 2-tuples of the name and **implicit value** (i.e., value
    retrieved by implicitly calling the ``@property``-decorated method
    implementing this attribute if this attribute is a property *or* the value
    of this attribute as is otherwise) of each attribute bound to the passed
    object whose name and/or value matches the passed predicate (in ascending
    lexicographic order of attribute name).

    Specifically, for each such attribute, the value of this attribute yielded
    by this generator is:

    * If this attribute is a **property** (i.e., method decorated by the
      standard ``@property`` decorator),  the high-level value returned by
      implicitly querying the low-level data descriptor underlying this
      property rather than that data descriptor.
    * Else, the value of this attribute as is.

    Caveats
    ----------
    This generator is substantially more dangerous than the comparable
    :func:`iter_attrs_simple_matching` generator. Whereas this generator
    implicitly calls the low-level method implementing each high-level property
    of the passed object and hence raises exceptions if any such method raises
    exceptions, that generator *never* raises unexpected exceptions. Unless
    properties are of interest, callers are strongly encouraged to call that
    rather than this generator.

    Parameters
    ----------
    obj : object
        Object to yield all matching attributes of.
    predicate : CallableTypes
        Callable iteratively passed both the name and value of each attribute
        bound to this object, returning ``True`` only if that name and/or value
        matches this predicate.

    Yields
    ----------
    (attr_name : str, attr_value : object)
        2-tuple of the name and value of each matching attribute bound to this
        object (in ascending lexicographic order of attribute name).
    '''

    # Return a generator comprehension...
    return (
        # Yielding this attribute's name and implicit value...
        (attr_name, attr_value)
        # For the name and definition of each attribute bound to this object...
        for attr_name, attr_value in inspect.getmembers(obj)
        # If this attribute matches this predicate.
        if predicate(attr_name, attr_value)
    )


@type_check
def iter_attrs_explicit_matching(
    obj: object, predicate: CallableTypes) -> GeneratorType:
    '''
    Generator yielding 2-tuples of the name and **explicit value** (i.e., value
    retrieved *without* implicitly calling the ``@property``-decorated method
    implementing this attribute if this attribute is a property) of each
    attribute bound to the passed object whose name and/or value also matches
    the passed predicate (in ascending lexicographic order of variable name).

    Specifically, for each such attribute, the value of this attribute yielded
    by this generator is:

    * If this attribute is a **property** (i.e., method decorated by the
      standard ``@property`` decorator), the low-level data descriptor
      underlying this property rather than the high-level value returned by
      implicitly querying that data descriptor.
    * Else, the value of this attribute as is.

    Caveats
    ----------
    This generator is substantially safer than the comparable
    :func:`iter_attrs_implicit_matching` generator, which implicitly calls the
    low-level method implementing each high-level property of the passed object
    and hence raises exceptions if any such method raises exceptions. By
    compare, this generator *never* raises unexpected exceptions. Unless
    properties are of interest, callers are strongly encouraged to call this
    rather than that generator.

    Only variables statically registered in this object's internal dictionary
    (e.g., ``__dict__`` in unslotted objects) are yielded. Variables
    dynamically defined by this object's ``__getattr__()`` method or related
    runtime magic are simply ignored.

    Parameters
    ----------
    obj : object
        Object to yield all matching non-property attributes of.
    predicate : CallableTypes
        Callable iteratively passed both the name and value of each attribute
        bound to this object, returning ``True`` only if that name and/or value
        matches this predicate.

    Yields
    ----------
    (var_name : str, var_value : object)
        2-tuple of the name and value of each matching non-property attribute
        bound to this object (in ascending lexicographic order of attribute
        name).
    '''

    # Ideally, this function would be reimplemented in terms of the
    # iter_attrs_implicit_matching() function calling the canonical
    # inspect.getmembers() function. Dynamic inspection is surprisingly
    # non-trivial in the general case, particularly when virtual base classes
    # rear their diamond-studded faces. Moreover, doing so would support
    # edge-case attributes when passed class objects, including:
    #
    # * Metaclass attributes of the passed class.
    # * Attributes decorated by "@DynamicClassAttribute" of the passed class.
    #
    # Sadly, inspect.getmembers() internally accesses attributes via the
    # dangerous getattr() builtin rather than the safe inspect.getattr_static()
    # function. This function explicitly requires the latter and hence *MUST*
    # reimplement rather than defer to inspect.getmembers(). (Sadness reigns.)
    #
    # For the same reason, the unsafe vars() builtin cannot be called either.
    # Since that builtin fails for builtin containers (e.g., "dict", "list"),
    # this is not altogether a bad thing.
    for attr_name in dir(obj):
        # Value of this attribute guaranteed to be statically rather than
        # dynamically retrieved. The getattr() builtin performs the latter,
        # dynamically calling this attribute's getter if this attribute is
        # a property. Since that call could conceivably raise unwanted
        # exceptions *AND* since this function explicitly ignores
        # properties, static attribute retrievable is unavoidable.
        attr_value = inspect.getattr_static(obj, attr_name)

        # If this attribute matches this predicate...
        if predicate(attr_name, attr_value):
            # Yield this attribute's name and explicit value. Note that, due to
            # the above assignment, this iteration *CANNOT* reasonably be
            # optimized into a generator comprehension.
            yield attr_name, attr_value

# ....................{ GETTERS : method                  }....................
@type_check
def iter_attr_names(obj: object) -> GeneratorType:
    '''
    Generator yielding the name of each attribute bound to the passed object
    (in ascending lexicographic order of attribute name).

    Parameters
    ----------
    obj : object
        Object to iterate the attribute names of.

    Yields
    ----------
    attr_name : str
        Name of each attribute bound to this object (in ascending lexicographic
        order of attribute name).
    '''

    # Defer to this generator.
    return (
        attr_name
        for attr_name, attr_value in iter_attrs(obj)
    )

# ....................{ ITERATORS ~ methods               }....................
def iter_methods(obj: object) -> GeneratorType:
    '''
    Generator yielding a 2-tuple of the name and value of each method bound to
    the passed object (in ascending lexicographic order of method name).

    This includes:

    * All methods statically registered in this object's internal dictionary
      (e.g., ``__dict__`` in unslotted objects), including both:

      * Builtin methods, whose names are both prefixed and suffixed by ``__``.
      * Custom methods, whose names are *not* prefixed and suffixed by ``__``.

    * All methods dynamically returned by **property methods** (i.e., methods
      decorated by the :func:`property` decorator).

    This excludes:

    * All property methods themselves.
    * All methods dynamically defined by this object's ``__getattr__()`` method
      or related runtime magic, which cannot (by definition) be inspected by
      this generator.

    Parameters
    ----------
    obj : object
        Object to iterate all methods of.

    Yields
    ----------
    (method_name : str, method : MethodType)
        2-tuple of the name and value of each method bound to this object (in
        ascending lexicographic order of method name).

    Raises
    ----------
    Exception
        If any method of this object is a property whose
        :func:`property`-decorated method raises an exception. To avoid
        handling such exceptions, consider instead calling the
        :func:`iter_methods_custom_simple` generator excluding all properties.
    '''

    # Defer to this generator.
    yield from iter_methods_matching(obj=obj, predicate=noop_attr_predicate)


@type_check
def iter_methods_matching(
    obj: object, predicate: CallableTypes) -> GeneratorType:
    '''
    Generator yielding 2-tuples of the name and value of each method bound to
    the passed object whose method name and value matches the passed predicate
    (in ascending lexicographic order of method name).

    Parameters
    ----------
    obj : object
        Object to yield all matching methods of.
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

    # Defer to this generator.
    yield from iter_attrs_explicit_matching(
        # For the name and definition of each attribute bound to this object,
        # exclude this attribute from consideration unless...
        obj=obj, predicate=lambda attr_name, attr_value: (
            # This attribute is a method (rather than variable) *AND*...
            callable(attr_value) and
            # This variable satisfies the passed predicate.
            predicate(attr_name, attr_value)))


def iter_methods_custom(obj: object) -> GeneratorType:
    '''
    Generator yielding 2-tuples of the name and value of each **non-builtin
    method** (i.e., method whose name is *not* both prefixed and suffixed by
    `__`) bound to the passed object (in ascending lexicographic order of
    method name).

    Parameters
    ----------
    obj : object
        Object to yield all non-builtin methods of.

    Yields
    ----------
    (method_name : str, method : MethodType)
        2-tuple of the name and value of each non-builtin method bound to this
        object (in ascending lexicographic order of method name).

    See Also
    ----------
    :func:`iter_methods`
        Further details.
    '''
    # Avoid circular import dependencies.
    from betse.util.py import pyident


    # Defer to this generator.
    yield from iter_methods_matching(
        obj=obj, predicate=lambda method_name, method: not (
            pyident.is_special(method_name)))


def iter_methods_prefixed(obj: object, prefix: str) -> GeneratorType:
    '''
    Generator yielding 2-tuples of the name and value of each method bound to
    the passed object whose method name is prefixed by the passed substring (in
    ascending lexicographic order of method name).

    Parameters
    ----------
    obj : object
        Object to yield all prefixed methods of.
    prefix : str
        Substring prefixing the names of all such methods.

    Yields
    ----------
    (method_name : str, method_value : object)
        2-tuple of the name and value of each prefixed method bound to this
        object (in ascending lexicographic order of method name).
    '''

    # Defer to this generator.
    yield from iter_methods_matching(
        obj=obj, predicate=lambda method_name, method_value: (
            method_name.startswith(prefix)))

# ....................{ ITERATORS ~ vars                  }....................
def iter_vars(obj: object) -> GeneratorType:
    '''
    Generator yielding a 2-tuple of the name and value of each variable bound
    to the passed object (in ascending lexicographic order of variable name).

    This includes:

    * All variables statically registered in this object's internal dictionary
      (e.g., ``__dict__`` in unslotted objects), including both:

      * Builtin variables, whose names are both prefixed and suffixed by
        ``__``.
      * Custom variables, whose names are *not* prefixed and suffixed by
        ``__``.

    * All variables dynamically returned by **property methods** (i.e., methods
      decorated by the :func:`property` decorator).

    This excludes *only* variables dynamically defined by this object's
    ``__getattr__()`` method or related runtime magic, which cannot (by
    definition) be inspected by this generator.

    Parameters
    ----------
    obj : object
        Object to iterate the variables of.

    Yields
    ----------
    (var_name : str, var_value : object)
        2-tuple of the name and value of each variable bound to this object (in
        ascending lexicographic order of variable name).

    Raises
    ----------
    Exception
        If any variable of this object is a property whose
        :func:`property`-decorated method raises an exception. To avoid
        handling such exceptions, consider instead calling the
        :func:`iter_vars_custom_simple` generator excluding all properties.
    '''

    # Defer to this generator.
    yield from iter_vars_matching(obj=obj, predicate=noop_attr_predicate)


@type_check
def iter_vars_matching(
    obj: object, predicate: CallableTypes) -> GeneratorType:
    '''
    Generator yielding 2-tuples of the name and value of each variable bound
    to the passed object whose name and/or value matches the passed predicate
    (in ascending lexicographic order of attribute name).

    Parameters
    ----------
    obj : object
        Object to yield all matching variables of.
    predicate : CallableTypes
        Callable iteratively passed both the name and value of each variables
        bound to this object, returning ``True`` only if that name and/or value
        matches this predicate.

    Yields
    ----------
    (var_name : str, var_value : object)
        2-tuple of the name and value of each matching variable bound to this
        object (in ascending lexicographic order of variable name).
    '''

    #FIXME: While satisfactory, the current approach fails to return the full
    #set of all instance variables in the edge case that the values of one or
    #more instance variables of the passed object are themselves methods.
    #
    #A more general-purpose (albeit currently untested) approach might be to
    #define the set of all instance variables as the set difference between the
    #set of all object attributes minus the set of all object class attributes.
    #In optimistic theory, the resulting set *SHOULD* be that desired. That
    #said, this is sufficiently fragile that we need to guarantee this with one
    #or ideally more unit tests. A preliminary implementation might resemble:
    #
    #    obj_type_attr_names = set(iter_attr_names(type(obj)))
    #
    #    yield from iter_attrs_implicit_matching(
    #        obj=obj, predicate=lambda attr_name, attr_value: (
    #            (not callable(attr_value) or
    #             attr_name not in obj_type_attr_names)) and
    #           predicate(attr_name, attr_value))

    # Defer to this generator.
    yield from iter_attrs_implicit_matching(
        # For the name and definition of each attribute bound to this object,
        # exclude this attribute from consideration unless...
        obj=obj, predicate=lambda attr_name, attr_value: (
            # This attribute is *NOT* a method and hence is a variable *AND*...
            not callable(attr_value) and
            # This variable satisfies the passed predicate.
            predicate(attr_name, attr_value)))


def iter_vars_custom(obj: object) -> GeneratorType:
    '''
    Generator yielding 2-tuples of the name and value of each **non-builtin
    variable** (i.e., variable whose name is *not* both prefixed and suffixed
    by `__`) bound to the passed object (in ascending lexicographic order of
    variable name).

    Parameters
    ----------
    obj : object
        Object to yield all non-builtin variables of.

    Yields
    ----------
    (var_name : str, var_value : object)
        2-tuple of the name and value of each non-builtin variable bound to
        this object (in ascending lexicographic order of variable name).

    See Also
    ----------
    :func:`iter_vars`
        Further details.
    '''

    # Avoid circular import dependencies.
    from betse.util.py import pyident

    # Defer to this generator.
    yield from iter_vars_matching(
        obj=obj, predicate=lambda var_name, var_value: not (
            pyident.is_special(var_name)))

# ....................{ ITERATORS ~ vars : custom simple  }....................
def iter_vars_custom_simple(obj: object) -> GeneratorType:
    '''
    Generator yielding 2-tuples of the name and value of each **non-builtin
    non-property variable** (i.e., variable whose name is *not* both prefixed
    and suffixed by ``__`` and whose value is *not* dynamically defined by the
    ``@property`` decorator to be the implicit result of a method call) bound
    to the passed object (in ascending lexicographic order of variable name).

    Only variables statically registered in this object's internal dictionary
    (e.g., ``__dict__`` in unslotted objects) are yielded. Variables
    dynamically defined by this object's ``__getattr__()`` method or related
    runtime magic are simply ignored.

    Parameters
    ----------
    obj : object
        Object to yield all non-builtin non-property variables of.

    Yields
    ----------
    (var_name : str, var_value : object)
        2-tuple of the name and value of each non-builtin non-property variable
        bound to this object (in ascending lexicographic order of variable
        name).
    '''

    # Defer to this generator with the null predicate unconditionally accepting
    # all such variables regardless of name or value.
    yield from iter_vars_custom_simple_matching(
        obj=obj, predicate=noop_attr_predicate)


def iter_vars_custom_simple_prefixed(
    obj: object, prefix: str) -> GeneratorType:
    '''
    Generator yielding 2-tuples of the name and value of each **non-builtin
    non-property prefixed variable** (i.e., variable whose name is prefixed by
    the passed substring and *not* both prefixed and suffixed by ``__`` and
    whose value is *not* dynamically defined by the ``@property`` decorator to
    be the implicit result of a method call) bound to the passed object (in
    ascending lexicographic order of variable name).

    Only variables statically registered in this object's internal dictionary
    (e.g., ``__dict__`` in unslotted objects) are yielded. Variables
    dynamically defined by this object's ``__getattr__()`` method or related
    runtime magic are simply ignored.

    Parameters
    ----------
    obj : object
        Object to yield all non-builtin non-property prefixed variables of.
    prefix : str
        Substring prefixing the names of all such variables.

    Yields
    ----------
    (var_name : str, var_value : object)
        2-tuple of the name and value of each non-builtin non-property prefixed
        variable bound to this object (in ascending lexicographic order of
        variable name).
    '''

    # Defer to this generator.
    yield from iter_vars_custom_simple_matching(
        obj=obj, predicate=lambda var_name, var_value: (
            var_name.startswith(prefix)))


#FIXME: Reimplement in terms of the newly defined
#iter_attrs_explicit_matching() generator.
@type_check
def iter_vars_custom_simple_matching(
    obj: object, predicate: CallableTypes) -> GeneratorType:
    '''
    Generator yielding 2-tuples of the name and value of each **non-builtin
    non-property variable** (i.e., variable whose name is *not* both prefixed
    and suffixed by ``__`` and whose value is *not* dynamically defined by the
    ``@property`` decorator to be the implicit result of a method call) bound
    to the passed object whose name and/or value also matches the passed
    predicate (in ascending lexicographic order of variable name).

    Only variables statically registered in this object's internal dictionary
    (e.g., ``__dict__`` in unslotted objects) are yielded. Variables
    dynamically defined by this object's ``__getattr__()`` method or related
    runtime magic are simply ignored.

    Parameters
    ----------
    obj : object
        Object to yield all non-builtin non-property variables of.
    predicate : CallableTypes
        Callable iteratively passed both the name and value of each attribute
        bound to this object, returning ``True`` only if that name and/or value
        matches this predicate.

    Yields
    ----------
    (var_name : str, var_value : object)
        2-tuple of the name and value of each matching non-builtin non-property
        variable bound to this object (in ascending lexicographic order of
        variable name).
    '''

    # Avoid circular import dependencies.
    from betse.util.py import pyident

    # Ideally, this function would be reimplemented in terms of the
    # iter_attrs_implicit_matching() function calling the canonical
    # inspect.getmembers() function. Dynamic inspection is surprisingly
    # non-trivial in the general case, particularly when virtual base classes
    # rear their diamond-studded faces. Moreover, doing so would support
    # edge-case attributes when passed class objects, including:
    #
    # * Metaclass attributes of the passed class.
    # * Attributes decorated by "@DynamicClassAttribute" of the passed class.
    #
    # Sadly, inspect.getmembers() internally accesses attributes via the
    # dangerous getattr() builtin rather than the safe inspect.getattr_static()
    # function. This function explicitly requires the latter and hence *MUST*
    # reimplement rather than defer to inspect.getmembers(). (Sadness reigns.)
    #
    # For the same reason, the unsafe vars() builtin cannot be called either.
    # Since that builtin fails for builtin containers (e.g., "dict", "list"),
    # this is not altogether a bad thing.
    for attr_name in dir(obj):
        # If this attribute is *NOT* a builtin...
        if not pyident.is_special(attr_name):
            # Value of this attribute guaranteed to be statically rather than
            # dynamically retrieved. The getattr() builtin performs the latter,
            # dynamically calling this attribute's getter if this attribute is
            # a property. Since that call could conceivably raise unwanted
            # exceptions *AND* since this function explicitly ignores
            # properties, static attribute retrievable is unavoidable.
            attr_value = inspect.getattr_static(obj, attr_name)

            # If all of the following conditions apply:
            if (
                # This value is neither callable nor a property, this attribute
                # must necessarily be a non-property variable.
                not (
                    callable(attr_value) or
                    isinstance(attr_value, property)
                ) and
                # This non-property variable matches this predicate. For
                # safety, this predicate is passed this attribute only *AFTER*
                # this attribute has been guaranteed to be such a variable;
                # this predicate may thus safely assume this to be the case.
                predicate(attr_name, attr_value)
            # ...yield this matching variable.
            ):
                yield attr_name, attr_value

# ....................{ PREDICATES                        }....................
def noop_attr_predicate(attr_name: object, attr_value: object) -> bool:
    '''
    Null attribute predicate unconditionally returning ``True`` for all passed
    attribute name and value combinations, suitable for passing as the value of
    the ``predicate`` parameter accepted by attribute iterators in this
    submodule (e.g., :func:`iter_attrs_implicit_matching`).
    '''

    # ...our work is done here, folks.
    return True
