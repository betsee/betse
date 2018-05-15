#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level **decorator** (i.e., callables dynamically wrapping arbitrary other
callables) facilities.
'''

# ....................{ IMPORTS                            }....................
from betse.util.io.log import logs
from betse.util.type.types import (
    type_check, CallableTypes, DecoratorTypes)
from functools import wraps

# ....................{ DECORATORS                         }....................
# While type-checking these types would probably be advisable, this decorator
# is currently passed non-callables by "py.test". While resolving that issue
# would itself probably be advisable, we simply cannot be bothered. Hence,
# these types remain blithely unchecked.
def decorator_identity(func: CallableTypes) -> CallableTypes:
    '''
    Identity decorator returning the decorated callable unmodified.
    '''

    return func

# ....................{ DECORATORS                         }....................
@type_check
def deprecated(func: CallableTypes) -> CallableTypes:
    '''
    Decorate the passed **callable** (e.g., function, method) to log a
    non-fatal warning on the first and *only* first call of this callable.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.call import callables

    # True only if this callable has had a deprecation warning logged.
    func.__is_deprecation_logged = False

    # Callable returned by this decorator below.
    @wraps(func)
    def _deprecated_inner(*args, **kwargs) -> object:
        # If this callable has *NOT* had a deprecation warning logged, do so.
        if not func.__is_deprecation_logged:
            # Prevent this warning from being logged more than once.
            func.__is_deprecation_logged = True

            # Capitalized human-readable string describing this callable.
            func_name = callables.to_str_capitalized(func)

            # Log this warning.
            logs.log_warning('%s deprecated.', func_name)

        # Call this callable.
        return func(*args, **kwargs)

    # Return this callable decorated to log this deprecating warning.
    return _deprecated_inner

# ....................{ CHAINERS                           }....................
@type_check
def decorator_chain(*decorators: DecoratorTypes) -> CallableTypes:
    '''
    Decorator chaining together (i.e., functionally composing) all passed
    decorators (in the reverse order).

    This function dynamically synthesizes a closure iteratively decorating any
    passed callable with each passed decorator. To preserve the expected
    semantics, that callable is decorated with:

    * The first passed decorator last.
    * The last passed decorator first.

    Parameters
    ----------
    decorators : tuple[DecoratorTypes]
        Tuple of all decorators to be chained together into a new decorator.

    Returns
    ----------
    CallableTypes
        Output decorator chaining together these input decorators.

    See Also
    ----------
    https://stackoverflow.com/a/5409569/2809027
        StackOverflow answer from which this implementation is inspired.
    '''

    # Closure to be returned as the new decorator.
    #
    # Note that the entirety of this function is reducible to this one-liner:
    #
    #     return lambda x: reduce(lambda y, f: f(y), decs, x)
    #
    # For maistainability, readability, and sanity, this function is instead
    # implemented with a traditional nested closure. Efficiency is irrelevant in
    # this context; moreover, the above one-liner is unlikely to yield
    # demonstrable efficiency gains.
    @type_check
    def _decorator_chain(func: CallableTypes) -> object:
        '''
        Iteratively decorate the passed callable with all decorators passed to
        the parent :func:`decorator_chain` function (in the reverse order),
        returning the callable object returned by the last such decorator.

        Parameters
        ----------
        call : CallableTypes
            Callable to be decorated.

        Returns
        ----------
        object
            Callable object decorating the input callable with these decorators.
            While this object is typically a callable (i.e., instance of
            :class:`CallableTypes`), this is *not* necessarily the case. For
            example, the :class:`property` decorator converts the passed unbound
            method into a data descriptor; strictly speaking, descriptors do
            *not* define the special ``__call__`` method and hence are *not*
            explicitly callable.
        '''

        # For each passed decorator (in reverse order)...
        for decorator in reversed(decorators):
            # Callable decorating the prior such callable with this decorator.
            func = decorator(func)

        # Return this decorated callable.
        return func

    # Return this new decorator.
    return _decorator_chain
