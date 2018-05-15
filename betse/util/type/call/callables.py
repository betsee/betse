#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level **callable** (e.g., function, lambda, method, property) facilities.
'''

# ....................{ IMPORTS                            }....................
from betse.util.io.log import logs
from betse.util.type.types import (
    type_check,
    CallableTypes,
    CallablePartialType,
    BuiltinFunctionType,
    BuiltinMethodType,
    FunctionType,
    FunctionTypes,
    GeneratorType,
    MethodType,
    MethodTypes,
    MappingOrNoneTypes,
    SequenceOrNoneTypes,
    StrOrNoneTypes,
)
from functools import partial, wraps

if False: wraps  # silence contemptible IDE warnings

# ....................{ CONSTANTS                          }....................
TYPE_TO_NAME = {
    BuiltinFunctionType: 'builtin',
    BuiltinMethodType:   'method',
    FunctionType:        'function',
    GeneratorType:       'generator',
    MethodType:          'method',
}
'''
Dictionary mapping from the class of each possible callable to the lowercase
human-readable type of that class, ignoring low-level implementation details
(e.g., the distinction between built-in and user-defined methods).
'''

# ....................{ TESTERS                            }....................
@type_check
def is_function(func: CallableTypes) -> bool:
    '''
    ``True`` only if the passed callable is a function.
    '''

    return isinstance(func, FunctionTypes)


@type_check
def is_method(func: CallableTypes) -> bool:
    '''
    ``True`` only if the passed callable is a method.
    '''

    return isinstance(func, MethodTypes)

# ....................{ GETTERS                            }....................
@type_check
def get_name(func: CallableTypes) -> str:
    '''
    Machine-readable name of the passed callable.

    Parameters
    ----------
    func: CallableTypes
        Callable to get the name of.

    Returns
    ----------
    str
        Machine-readable name of this callable.

    See Also
    ----------
    :func:`to_str`
        Function getting a human-readable name of this callable.
    '''

    # The easy things in life are free.
    return func.__name__


@type_check
def get_doc_or_none(func: CallableTypes) -> StrOrNoneTypes:
    '''
    Human-readable docstring documenting the passed callable if any *or*
    ``None`` otherwise.

    Parameters
    ----------
    func: CallableTypes
        Callable to get the docstring of.

    Returns
    ----------
    StrOrNoneTypes
        Either:
        * If this callable has a docstring, this docstring.
        * Else, ``None``.
    '''

    return getattr(func, '__doc__', None)

# ....................{ SETTERS                            }....................
@type_check
def set_doc(func: CallableTypes, doc: str) -> None:
    '''
    Set the human-readable docstring documenting the passed callable to the
    passed string.

    Parameters
    ----------
    func: CallableTypes
        Callable to set the docstring of.
    doc : str
        Docstring to set on this callable.
    '''

    func.__doc__ = doc

# ....................{ CONVERTERS                         }....................
@type_check
def to_str(func: CallableTypes) -> str:
    '''
    Human-readable string describing the type and name of the passed callable
    (e.g., ``method Pala()``).

    For generality, this string's:

    * First character is guaranteed to be lowercase.
    * Last character is guaranteed to be the ``)`` delimiter.

    Parameters
    ----------
    func: CallableTypes
        Callable to describe.

    Returns
    ----------
    str
        Human-readable string describing this callable.
    '''

    # Type of this callable, defaulting to merely "callable" if unrecognized.
    func_type = TYPE_TO_NAME.get(type(func), 'callable')
    # print('func: {}, type: {}'.format(func.__name__, type(func)))

    # Return a human-readable string describing this callable.
    return '{} {}()'.format(func_type, func.__name__)


@type_check
def to_str_capitalized(func: CallableTypes) -> str:
    '''
    Human-readable string describing the type and name of the passed callable
    whose first character is capitalized (e.g., ``Method pala()``).

    See Also
    ----------
    :func:`to_str`
        Further details.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.text import strs

    # Return a capitalized human-readable string describing this callable.
    return strs.uppercase_char_first(to_str(func))

# ....................{ MAKERS                             }....................
@type_check
def make_partial(
    # Mandatory parameters.
    func: CallableTypes,

    # Optional parameters.
    args: SequenceOrNoneTypes = None,
    kwargs: MappingOrNoneTypes = None,
    doc: StrOrNoneTypes = None,
) -> CallablePartialType:
    '''
    **Partial callable** (i.e., callable whose signature is reduced by this
    function call to unconditionally pass one or more positional or keyword
    arguments to the passed callable) calling the passed callable, statically
    passed the passed positional and keyword arguments, and documented by the
    passed docstring.

    Parameters
    ----------
    func: CallableTypes
        Callable to partialize.
    args : SequenceOrNoneTypes
        Sequence of all positional arguments to unconditionally pass to the
        passed callable if any *or* ``None`` otherwise. Defaults to ``None``.
    kwargs : MappingOrNoneTypes
        Dictionary of all keyword arguments to unconditionally pass to the
        passed callable if any *or* ``None`` otherwise. Defaults to ``None``.
    doc : StrOrNoneTypes
        Docstring to document the returned callable with if any *or* ``None``
        otherwise. If ``None`` *and* the passed callable is documented, the
        returned callable will share the same docstring. Defaults to ``None``.

    Returns
    ----------
    CallablePartialType
        Partial callable synthesized from the passed callable and arguments.
    '''

    # Default all unpassed arguments to sane values.
    if args is None:
        args = ()
    if kwargs is None:
        kwargs = {}
    if doc is None:
        doc = get_doc_or_none(func)

    # Partial callable synthesized from the passed callable and arguments.
    func_partial = partial(func, *args, **kwargs)

    # If documenting this callable, do so.
    if doc is not None:
        set_doc(func=func_partial, doc=doc)

    # Return this partial callable.
    return func_partial
