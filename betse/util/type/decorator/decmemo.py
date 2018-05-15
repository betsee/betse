#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level **decorator-based memoization** (i.e., efficient caching and reuse of
the values returned by decorated callables) facilities.
'''

# ....................{ IMPORTS                            }....................
from betse.util.type.types import type_check, CallableTypes, PropertyType
from functools import wraps

if False: wraps  # silence contemptible IDE warnings

# ....................{ CONSTANTS                          }....................
CALLABLE_CACHED_VAR_NAME_PREFIX = '_betse_cached__'
'''
Substring prefixing the names of all private instance variables to which all
caching decorators (e.g., :func:`property_cached`) cache values returned by
decorated callables.

This prefix guarantees uniqueness across *all* instances -- including those
instantiated from official Python and unofficial third-party classes and those
internally defined by this application. Doing so permits logic elsewhere (e.g.,
pickling filtering) to uniquely match and act upon these variables.
'''


FUNCTION_CACHED_VAR_NAME = CALLABLE_CACHED_VAR_NAME_PREFIX + 'function_value'
'''
Name of the private instance variable to which the :func:`func_cached`
decorator statically caches the value returned by the decorated function.
'''


PROPERTY_CACHED_VAR_NAME_PREFIX = CALLABLE_CACHED_VAR_NAME_PREFIX + 'property_'
'''
Substring prefixing the names of all private instance variables to which the
:func:`property_cached` decorator dynamically caches the value returned by the
decorated property method.
'''

# ....................{ DECORATORS                         }....................
@type_check
def func_cached(func: CallableTypes) -> CallableTypes:
    '''
    Decorate the passed **non-property callable** (e.g., function, lambda,
    method) to cache the value returned by the first call of this callable.

    On the first call of a callable decorated with this decorator, the passed
    callable is called with all passed parameters, the value returned by this
    callable is cached into a private attribute of this callable, and this value
    is returned. On each subsequent call of this callable, the cached value is
    returned as is *without* calling this callable. Hence, this callable is
    called at most once for each instance of the class containing this property.

    Caveats (Memoization)
    ----------
    **This decorator does not memoize callables.** Memoization would map each
    permutation of parameters passed to the decorated callable to a unique
    return value conditionally cached (and hence returned) for that permutation.

    This decorator instead unconditionally caches (and hence returns) a single
    return value for *all* permutations of passed parameters. Hence, this
    decorator is principally intended to decorate callables accepting *no*
    parameters (e.g., simple testers and getters).

    Caveats (Methods)
    ----------
    **This decorator is only intended to cache return values of functions.**
    This decorator is *not* intended to cache return values of non-function
    callables (e.g., bound or unbound methods). Ergo, this decorator is named
    ``func_cached`` rather than the more general-purpose name
    ``callable_cached``.

    While this decorator *could* technically be used to decorate non-function
    callables, it never should be. Why? Because:

    - Decorating lambdas is highly non-trivial and hence impractical.
    - Attempting to cache return values of unbound methods with this decorator
      would do so globally rather than on an instance-specific basis and hence
      be functionally useless. To cache return values of bound method on an
      instance-specific basis, an entirely new decorator would need to be
      created -- presumably returning a descriptor whose ``__get__()`` special
      method creates and returns bound methods whose return values are cached
      onto instance variables of those bound methods themselves. For similar
      functionality, see the
      ``betse.util.type.descriptor.expralias.expr_alias`` data descriptor.

    See Also
    ----------
    :func:`property_cached`
        Decorator caching object properties.
    '''

    # Raw string of Python statements comprising the body of this wrapper.
    # See @property_cached for further details, including time profilings.
    #
    # Technically, note that this decorator *SHOULD* raise a
    # "BetseCallableException" if any parameter accepted by this callable has
    # the reserved name `__callable_cached`. However, validating this constraint
    # would require inspecting this callable's signature and iteratively
    # searching the parameter list of that signature for parameter names
    # violating this constraint. Since the likelihood of a callable accepting
    # such parameters is "None", doing so would increase time complexity
    # *WITHOUT* yielding tangible benefits. (Also, we're lazy.)
    func_body = '''
@wraps(__callable_cached)
def callable_cached_arged(
    *args, __callable_cached=__callable_cached, **kwargs):
    try:
        return __callable_cached.{callable_var_name}
    except AttributeError:
        __callable_cached.{callable_var_name} = __callable_cached(
            *args, **kwargs)
        return __callable_cached.{callable_var_name}
'''.format(callable_var_name=FUNCTION_CACHED_VAR_NAME)

    # Dictionary mapping from local attribute names to values. For efficiency,
    # only attributes required by the body of this wrapper are copied from the
    # current namespace. (See below.)
    local_attrs = {'__callable_cached': func}

    # Dynamically define this wrapper as a closure of this decorator. For
    # obscure and presumably uninteresting reasons, Python fails to locally
    # declare this closure when the locals() dictionary is passed; to capture
    # this closure, a local dictionary must be passed instead.
    exec(func_body, globals(), local_attrs)

    # Return this wrapper method.
    return local_attrs['callable_cached_arged']


@type_check
def property_cached(property_method: CallableTypes) -> PropertyType:
    '''
    Decorate the passed property method to cache the value returned by the first
    implicit call of this method.

    On the first access of a property decorated with this decorator, the passed
    method implementing this property is called, the value returned by this
    property is internally cached into a private attribute of the object to
    which this method is bound, and this value is returned. On each subsequent
    access of this property, this cached value is returned as is _without_
    calling this method. Hence, this method is called at most once for each
    object exposing this property.

    Caveats
    ----------
    **This decorator does not destroy bound property methods.** Technically, the
    most efficient means of caching a property value into an instance is to
    replace the property method currently bound to that instance with an
    instance variable initialized to that value (e.g., as documented by this
    `StackOverflow answer`_).

    Since a property should only ever be treated as an instance variable,
    there superficially exists little harm in dynamically changing the type of
    the former to the latter. Sadly, doing so introduces numerous subtle issues
    with no plausible workaround.

    In particular, replacing property methods by instance variables prevents
    pickling logic elsewhere from automatically excluding cached property
    values, forcing these values to *always* be pickled to disk. This is bad.
    Cached property values are *always* safely recreatable in memory (and hence
    need *not* be pickled) and typically space-consumptive in memory (and hence
    best *not* pickled). The slight efficiency gain from replacing property
    methods by instance variables is hardly worth the significant space loss
    from pickling these variables.

    .. _StackOverflow answer:
        https://stackoverflow.com/a/36684652/2809027
    '''

    # Name of the private instance variable to which this decorator caches the
    # value returned by the decorated property method.
    property_var_name = (
        PROPERTY_CACHED_VAR_NAME_PREFIX + property_method.__name__)

    # Raw string of Python statements comprising the body of this wrapper,
    # including (in order):
    #
    # * A "@wraps" decorator propagating the name, docstring, and other
    #   identifying metadata of the original function to this wrapper.
    # * A private "__property_method" parameter initialized to this function.
    #   In theory, the "property_method" parameter passed to this decorator
    #   should be accessible as a closure-style local in this wrapper. For
    #   unknown reasons (presumably, a subtle bug in the exec() builtin), this
    #   is not the case. Instead, a closure-style local must be simulated by
    #   passing the "property_method" parameter to this function at function
    #   definition time as the default value of an arbitrary parameter.
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
        return self.{property_var_name}
    except AttributeError:
        self.{property_var_name} = __property_method(self)
        return self.{property_var_name}
'''.format(property_var_name=property_var_name)

    # Dictionary mapping from local attribute names to values. For efficiency,
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
