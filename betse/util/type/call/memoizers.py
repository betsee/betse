#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Low-level **memoization** (i.e., efficient caching and reuse of the values
returned by callables, typically via decorators) facilities.
'''

# ....................{ IMPORTS                            }....................
from betse.util.type.types import type_check, CallableTypes, PropertyType
from functools import wraps

if False: wraps  # silence contemptible IDE warnings

# ....................{ DECORATORS                         }....................
@type_check
def callable_cached(func: CallableTypes) -> CallableTypes:
    '''
    Decorate the passed **non-property callable** (e.g., function, lambda,
    method) to cache the value returned by the first call of this callable.

    On the first call of a callable decorated with this decorator, the passed
    callable is called with all passed parameters, the value returned by this
    callable is cached into a private attribute of this callable, and this value
    is returned. On each subsequent call of this callable, the cached value is
    returned as is _without_ calling this callable. Hence, this callable is
    called at most once for each instance of the class containing this property.

    Caveats
    ----------
    **This decorator does not memoize callables.** Memoization would map each
    permutation of parameters passed to the decorated callable to a unique
    return value conditionally cached (and hence returned) for that permutation.

    This decorator instead unconditionally caches (and hence returns) a single
    return value for _all_ permutations of passed parameters. Hence, this
    decorator is principally intended to decorate callables accepting _no_
    parameters (e.g., simple testers and getters).

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
        return __callable_cached.__callable_cached_value
    except AttributeError:
        __callable_cached.__callable_cached_value = __callable_cached(
            *args, **kwargs)
        return __callable_cached.__callable_cached_value
'''

    # Dictionary mapping from local attribute name to value. For efficiency,
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


#FIXME: Reimplement this to leverage this far more efficient alternative
#overwriting itself with the value returned by the decorated function:
#
#    https://stackoverflow.com/a/36684652/2809027
#
#Quite clever and *MUCH* more succinct than the approach pursued below. Note
#that this approach does change the type of this property from a callable to an
#instance variable. For this reason, this approach *CANNOT* be used to similarly
#improve the callable_cached() decorator defined above. In this case, however,
#doing so is perfectly safe: a property should only ever be treated as an
#instance variable, so there exists no harm in dynamically changing the type of
#the former to the latter. (Awesome.)
#
#In fact, the above approach is essentially what Django does (albeit with a more
#cumbersome class decorator than a function decorator). Which means, in turn,
#that this is indeed the canonical solution to caching properties.
#FIXME: Actually, the solution below may indeed be superior -- but for one
#particularly sutble reasons not commonly discussed: pickling. The third-party
#pickling package "dill" explicitly supports pickling of @property-decorated
#methods by pickling the data descriptor produced by those methods rather than
#the underlying data -- which is good. Reducing this approach to the attribute
#inlining scheme outlined above, however, would result in "dill" pickling the
#underlying data -- which would be very bad.
#FIXME: That said, the current approach still does *NOT* suffice. Why? Because
#"dill" will still pickle the underlying data cached by all
#@property_cached-decorated methods that have been called at least once, as this
#data is cached to simple instance variables of their parent object.
#Fortunately, these variables are *ALL* prefixed by the same substring:
#"__property_cached_". This implies that "dill" may be instructed to ignore all
#instance variables prefixed by this substring when pickling. Let it be, please.
#FIXME: While "dill" doesn't appear to support such functionality at present, a
#feature request for doing so has been opened at:
#    https://github.com/uqfoundation/dill/issues/225

# Note that, for unknown reasons, the "property_method" parameter cannot be
# assumed to be a method. Property methods appear to be of type function rather
# than method, presumably due to being decorated *BEFORE* being bound to a class
# as a method.
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
    '''

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
