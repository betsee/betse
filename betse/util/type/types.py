#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Low-level **type testers** (i.e., functions testing the types of passed
objects).
'''

#FIXME: Refactor all assertion statements performing callable type checking
#throughout this codebase into uses of the new @type_check decorator.

# ....................{ IMPORTS                            }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: To raise human-readable exceptions on missing mandatory dependencies
# *AND* avoid non-halting recursive imports when imported at the top-level
# of other modules in the `betse.util` package, this module may import *ONLY*
# from stock Python packages. (By definition, this excludes both BETSE and
# third-party packages.)
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

import inspect, re
from collections.abc import (
    Container, Hashable, Iterable, Iterator, Mapping, MutableMapping, Sequence,
    Set,
)
from enum import Enum, EnumMeta
from functools import wraps
from inspect import Parameter, Signature
from types import (
    BuiltinFunctionType,
    BuiltinMethodType,
    FunctionType,
    LambdaType,
    MethodType,
    ModuleType,
)

# ....................{ TYPES                              }....................
NoneType = type(None)
'''
Type (i.e., class object) of the singleton `None` object.

Curiously, although the type of the `None` object is a class object whose
`__name__` attribute is `NoneType`, there exists no globally accessible class
by that name. To circumvents this obvious oversight, this global globally
exposes this class.

This class is principally useful for annotating both:

* CallableTypes parameters accepting `None` as a valid value.
* Callables returning `None` as a valid value.
'''


HashableType = Hashable
'''
Abstract interface implemented by all **hashable objects** (i.e., implementing
the standard `__hash__()` method required by all dictionary keys).

This class is a synonym of the `collections.abc.Hashable` class, provided
merely as a convenience to callers preferring to avoid importing that class.
'''


IterableType = Iterable
'''
Abstract interface implemented by all **iterable objects** (i.e., implementing
the standard `__iter__()` method returning a new iterator object).

This class is a synonym of the `collections.abc.Iterable` class, provided
merely as a convenience to callers preferring to avoid importing that class.
'''


IteratorType = Iterator
'''
Abstract interface implemented by all **iterator objects** (i.e., implementing
the standard `__iter__()` and `__next__()` methods implementing iteration over
some container).

This class is a synonym of the `collections.abc.Iterator` class, provided
merely as a convenience to callers preferring to avoid importing that class.
'''


MappingType = Mapping
'''
Abstract interface implemented by all dictionary-like objects.

This class is a synonym of the `collections.abc.Mapping` class, provided merely
as a convenience to callers preferring to avoid importing that class.
'''


MutableMappingType = MutableMapping
'''
Abstract interface implemented by all mutable dictionary-like objects.

This class is a synonym of the `collections.abc.MutableMapping` class, provided
merely as a convenience to callers preferring to avoid importing that class.
'''


SetType = Set
'''
Abstract interface implemented by all set-like objects.

This class is a synonym of the `collections.abc.Set` class, provided merely as
a convenience to callers preferring to avoid importing that class.
'''

# ....................{ TUPLES                             }....................
CallableTypes = (
    BuiltinFunctionType,
    BuiltinMethodType,
    FunctionType,
    LambdaType,
    MethodType,
)
'''
Tuple of all **callable classes** (i.e., classes whose instances are callable
objects, including both built-in and user-defined functions, lambdas, and
methods).
'''


NumericTypes = (int, float)
'''
Tuple of all **numeric classes** (i.e., classes whose instances are single
scalar numbers excluding complex numbers, which typically require
special-purpose handling).

This tuple contains classes matching both integer and real number types.
'''


# This tuple is declared here for documentation purposes but initialized below.
SequenceTypes = None
'''
Tuple of all container base classes conforming to (but _not_ necessarily
subclassing) the canonical `collections.abc.Sequence` API.

Sequences are iterables supporting efficient element access via integer
indices. Most sequences implement the abstract base class
`collections.abc.Sequence`, including the concrete `str` string class. All
sequences define the special `__getitem__()` and `__len__()` methods, amongst
numerous others.

While all sequences are iterables, not all iterables are sequences. Generally
speaking, sequences correspond to the proper subset of iterables whose elements
are ordered. `dict` and `OrderedDict` are the canonical examples. `dict`
implements `collections.abc.Iterable` but _not_ `collections.abc.Sequence`, due
to _not_ supporting integer index-based lookup; `OrderedDict` implements both,
due to supporting such lookup.

For generality, this tuple contains classes matching both pure-Python sequences
_and_ non-Pythonic Fortran-based `numpy` arrays and matrices -- which fail to
subclass `collections.abc.Sequence` despite implementing the entirety of that
that API.
'''

# Add sequence types to the previously declared tuple.
#
# If Numpy is available, add both this API and the Numpy array type (which
# fails to subclass this API). Although Numpy is a mandatory dependency, this
# submodule is typically imported quite early in program startup, implying the
# importability of *ANY* dependency (mandatory or not) at the top level of this
# submodule to still be in question. Since subsequent logic in program startup
# is guaranteed to raise human-readable exceptions for missing dependencies,
# this error is silently ignored here.
try:
    from numpy import ndarray
    SequenceTypes = (Sequence, ndarray)
# Else, Numpy is unavailable. Add only this API.
except ImportError:
    SequenceTypes = (Sequence,)

# ....................{ SETS : private                     }....................
#FIXME: Type-check variadic keyword arguments as well.
_PARAMETER_KIND_IGNORED = {
    Parameter.POSITIONAL_ONLY, Parameter.VAR_KEYWORD,
}
'''
Set of all `inspect.Parameter.kind` constants to be ignored during annotation-
based type checking in the `@type_check` decorator.

This includes:

* Constants specific to variadic parameters (e.g., `*args`, `**kwargs`).
  Variadic parameters cannot be annotated and hence cannot be type checked.
* Constants specific to positional-only parameters, which apply to non-pure-
  Python callables (e.g., defined by C extensions). The `@type_check` decorator
  applies _only_ to pure-Python callables, which provide no syntactic means for
  specifying positional-only parameters.
'''


_RETURN_ANNOTATION_IGNORED = {Signature.empty, None}
'''
Set of all annotations for return values to be ignored during annotation-
based type checking in the `@type_check` decorator.

This includes:

* `Signature.empty`, signifying a callable whose return value is _not_
  annotated.
* `None`, signifying a callable returning no value. By convention, callables
  returning no value are typically annotated to return `None`. Technically,
  callables whose return values are annotated as `None` _could_ be explicitly
  checked to return `None` rather than a none-`None` value. Since return values
  are safely ignorable by callers, however, there appears to be little
  real-world utility in enforcing this constraint.
'''

# ....................{ DECORATORS                         }....................
# If the active Python interpreter is *NOT* optimized (e.g., option "-O" was
# *NOT* passed to this interpreter), enable type checking.
if __debug__:
    def type_check(func: CallableTypes) -> CallableTypes:
        '''
        Decorate the passed **callable** (e.g., function, method) to validate
        both all annotated parameters passed to this callable _and_ the
        annotated value returned by this callable if any.

        This decorator performs rudimentary type checking based on Python 3.x
        function annotations, as officially documented by PEP 484 ("Type
        Hints"). While PEP 484 supports arbitrarily complex type composition,
        this decorator requires _all_ parameter and return value annotations to
        be either:

        * Classes (e.g., `int`, `OrderedDict`).
        * Tuples of classes (e.g., `(int, OrderedDict)`).

        If optimizations are enabled by the active Python interpreter (e.g., due
        to option `-O` passed to this interpreter), this decorator is a noop.

        Raises
        ----------
        NameError
            If any parameter has the reserved name `__beartype_func`.
        TypeError
            If either:
            * Any parameter or return value annotation is neither:
              * A type.
              * A tuple of types.
            * The kind of any parameter is unrecognized. This should _never_
              happen, assuming no significant changes to Python semantics.
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
        # * Assert statements type checking parameters passed to this callable.
        # * A call to this callable.
        # * An assert statement type checking the value returned by this
        #   callable.
        #
        # While there exist numerous alternatives (e.g., appending to a list or
        # bytearray before joining the elements of that iterable into a string),
        # these alternatives are either slower (as in the case of a list, due to
        # the high up-front cost of list construction) or substantially more
        # cumbersome (as in the case of a bytearray). Since string concatenation
        # is heavily optimized by the official CPython interpreter, the simplest
        # approach is (curiously) the most ideal.
        func_body = '''
@wraps(__beartype_func)
def func_type_checked(*args, __beartype_func=__beartype_func, **kwargs):
'''

        # "inspect.Signature" instance encapsulating this callable's signature.
        func_sig = inspect.signature(func)

        # Human-readable name of this function for use in exceptions.
        func_name = func.__name__ + '()'

        # For the name of each parameter passed to this callable and the
        # "inspect.Parameter" instance encapsulating this parameter (in the
        # passed order)...
        for func_arg_index, func_arg in enumerate(func_sig.parameters.values()):
            # If this callable redefines a parameter initialized to a default
            # value by this wrapper, raise an exception. Permitting this
            # unlikely edge case would permit unsuspecting users to
            # "accidentally" override these defaults.
            if func_arg.name == '__beartype_func':
                raise NameError(
                    'Parameter {} reserved for use by @type_check.'.format(
                        func_arg.name))

            # If this parameter is both annotated and non-ignorable for purposes
            # of type checking, type check this parameter.
            if (func_arg.annotation is not Parameter.empty and
                func_arg.kind not in _PARAMETER_KIND_IGNORED):
                # Validate this annotation.
                _check_type_annotation(
                    annotation=func_arg.annotation,
                    label='{} parameter {} type'.format(
                        func_name, func_arg.name))

                # String evaluating to this parameter's annotated type.
                func_arg_type_expr = (
                    '__beartype_func.__annotations__[{!r}]'.format(
                        func_arg.name))

                # String evaluating to this parameter's current value when
                # passed as a keyword.
                func_arg_value_key_expr = 'kwargs[{!r}]'.format(func_arg.name)

                # If this parameter is actually a tuple of positional variadic
                # parameters, iteratively check all such parameters.
                if func_arg.kind is Parameter.VAR_POSITIONAL:
                    #FIXME: "{arg_name!r}" obviously isn't right. We need to
                    #find some means of obtaining the tuple or slice of all
                    #positionally passed parameters. No idea how, sadly. Let's
                    #google up "Parameter.VAR_POSITIONAL", please. Something's
                    #gotta give.
                    func_body += '''
    for __beartype_arg in args[{arg_index!r}:]:
        if not isinstance(__beartype_arg, {arg_type_expr}):
            raise TypeError(
                '{func_name} positional variadic parameter '
                '{arg_index}={{}} not a {{!r}}'.format(
                    trim(__beartype_arg), {arg_type_expr}))
'''.format(
                        func_name=func_name,
                        arg_name=func_arg.name,
                        arg_index=func_arg_index,
                        arg_type_expr=func_arg_type_expr,
                    )
                # Else if this parameter is keyword-only, check this parameter
                # only by lookup in the variadic "**kwargs" dictionary.
                elif func_arg.kind is Parameter.KEYWORD_ONLY:
                    func_body += '''
    if {arg_name!r} in kwargs and not isinstance(
        {arg_value_key_expr}, {arg_type_expr}):
        raise TypeError(
            '{func_name} keyword-only parameter '
            '{arg_name}={{}} not a {{!r}}'.format(
                trim({arg_value_key_expr}), {arg_type_expr}))
'''.format(
                        func_name=func_name,
                        arg_name=func_arg.name,
                        arg_type_expr=func_arg_type_expr,
                        arg_value_key_expr=func_arg_value_key_expr,
                    )
                # Else, this parameter may be passed either positionally or as
                # a keyword. Check this parameter both by lookup in the
                # variadic "**kwargs" dictionary *AND* by index into the
                # variadic "*args" tuple.
                else:
                    # String evaluating to this parameter's current value when
                    # passed positionally.
                    func_arg_value_pos_expr = 'args[{!r}]'.format(
                        func_arg_index)

                    func_body += '''
    if not (
        isinstance({arg_value_pos_expr}, {arg_type_expr})
        if {arg_index} < len(args) else
        isinstance({arg_value_key_expr}, {arg_type_expr})
        if {arg_name!r} in kwargs else True):
            raise TypeError(
                '{func_name} parameter {arg_name}={{}} not of {{!r}}'.format(
                trim({arg_value_pos_expr} if {arg_index} < len(args) else {arg_value_key_expr}),
                {arg_type_expr}))
'''.format(
                    func_name=func_name,
                    arg_name=func_arg.name,
                    arg_index=func_arg_index,
                    arg_type_expr=func_arg_type_expr,
                    arg_value_key_expr=func_arg_value_key_expr,
                    arg_value_pos_expr=func_arg_value_pos_expr,
                )

        # If this callable's return value is both annotated and non-ignorable
        # for purposes of type checking, type check this value.
        if func_sig.return_annotation not in _RETURN_ANNOTATION_IGNORED:
            # Validate this annotation.
            _check_type_annotation(
                annotation=func_sig.return_annotation,
                label='{} return type'.format(func_name))

            # Strings evaluating to this parameter's annotated type and
            # currently passed value, as above.
            func_return_type_expr = (
                "__beartype_func.__annotations__['return']")

            # Call this callable, type check the returned value, and return this
            # value from this wrapper.
            func_body += '''
    return_value = __beartype_func(*args, **kwargs)
    if not isinstance(return_value, {return_type}):
        raise TypeError(
            '{func_name} return value {{}} not of {{!r}}'.format(
                trim(return_value), {return_type}))
    return return_value
'''.format(func_name=func_name, return_type=func_return_type_expr)
        # Else, call this callable and return this value from this wrapper.
        else:
            func_body += '''
    return __beartype_func(*args, **kwargs)
'''

        # Dictionary mapping from local attribute name to value. For efficiency,
        # only those local attributes explicitly required in the body of this
        # wrapper are copied from the current namespace. (See below.)
        local_attrs = {'__beartype_func': func}

        # Dynamically define this wrapper as a closure of this decorator. For
        # obscure and presumably uninteresting reasons, Python fails to locally
        # declare this closure when the locals() dictionary is passed; to
        # capture this closure, a local dictionary must be passed instead.
        # print('\n{} wrapper: {}'.format(func_name, func_body))
        exec(func_body, globals(), local_attrs)

        # Return this wrapper.
        return local_attrs['func_type_checked']

# Else, the active Python interpreter is optimized. In this case, disable type
# checking by reducing this decorator to the identity decorator.
else:
    def type_check(func: callable) -> callable:
        return func


def _check_type_annotation(annotation: object, label: str) -> None:
    '''
    Validate the passed annotation to be a valid type supported by the
    `@type_check` decorator.

    Parameters
    ----------
    annotation : object
        Annotation to be validated.
    label : str
        Human-readable label describing this annotation, interpolated into
        exceptions raised by this function.

    Raises
    ----------
    TypeError
        If this annotation is neither a new-style class nor a tuple of
        new-style classes.
    '''

    # If this annotation is a tuple, raise an exception if any member of this
    # tuple is *NOT* a new-style class.
    if isinstance(annotation, tuple):
        for member in annotation:
            if not is_class_new(member):
                raise TypeError(
                    '{} tuple member {} not a new-style class'.format(
                        label, member))
    # Else if this annotation is *NOT* a new-style class, raise an exception.
    elif not is_class_new(annotation):
        raise TypeError(
            '{} {} neither a new-style class nor '
            'tuple of such classes'.format(label, annotation))

# ....................{ FORMATTERS                         }....................
# This string-centric function is defined in this module rather than in the
# arguably more appropriate "strs" module to drastically simplify assertion
# handlers defined below.
def trim(obj: object) -> str:
    '''
    Convert the passed object to a terse human-readable string suitable for use
    in end-user messages.
    '''

    # Maximum length of the string to be returned, defined to be the customary
    # line length of 80 characters minus:
    #
    # * Default output indentation of 4 characters.
    MAX_LEN = 76

    # Uncompiled regular expression grouping zero or more non-newline leading
    # characters preceding this maximum length *AND* zero or more trailing
    # delimiters.
    PRE_MAX_CHARS_LAST_DELIMITERS_REGEX = (
        r'^([^\n]{0,' + str(MAX_LEN) + r'}).*?([\])}>\'"]*)$')

    # String describing the passed object. For debuggability, the verbose
    # (albeit less human-readable) output of repr() is preferred to the terse
    # (albeit more human-readable) output of str().
    obj_synopsis = repr(obj)

    # If this synopsis either exceeds this maximum length *OR* contains a
    # newline, replace the substring of this synopsis from whichever of the
    # first character following this maximum length or the first newline occurs
    # first to the string end (excluding any # optional trailing delimiters)
    # with a single ellipses.
    if len(obj_synopsis) > MAX_LEN or '\n' in obj_synopsis:
        obj_synopsis = re.sub(
            PRE_MAX_CHARS_LAST_DELIMITERS_REGEX,
            r'\1...\2',
            obj_synopsis,
            flags=re.DOTALL
        )

    # Return this synopsis.
    return obj_synopsis

# ....................{ TESTERS                            }....................
def is_bool(obj: object) -> bool:
    '''
    `True` only if the passed object is **boolean** (i.e., either `True` or
    `False`).
    '''
    return isinstance(obj, bool)


def is_char(obj: object) -> bool:
    '''
    `True` only if the passed object is a **character** (i.e., string of length
    1).
    '''
    return is_str(obj) and len(obj) == 1


def is_nonnone(obj: object) -> bool:
    '''
    `True` only if the passed object is _not_ `None`.
    '''
    return obj is not None

# ....................{ TESTERS ~ callable                 }....................
def is_callable(obj: object) -> bool:
    '''
    `True` only if the passed object is **callable** (e.g., function, method,
    class defining the special `__call__()` method).
    '''

    return callable(obj)

# ....................{ TESTERS ~ class                    }....................
def is_class(obj: object) -> bool:
    '''
    `True` only if the passed object is a class.
    '''

    return isinstance(obj, type)


def is_class_new(obj: object) -> bool:
    '''
    `True` only if the passed object is a new-style class.
    '''

    # The "__name__" attribute is defined only by new-style classes and hence
    # serves as a useful means of distinguishing new- from old-style classes.
    return is_class(obj) and hasattr(obj, '__name__')

# ....................{ TESTERS ~ collection               }....................
def is_mapping(obj: object) -> bool:
    '''
    `True` only if the passed object is a **mapping** (i.e., indexable by
    strings).

    Customary mappings include `dict` and `OrderedDict` instances.
    '''
    return isinstance(obj, Mapping)

# ....................{ TESTERS ~ collection               }....................
def is_container(obj: object) -> bool:
    '''
    `True` only if the passed object is a **container** (i.e., implements the
    `__contains__` special method returning `True` only if that container
    contains the passed element).

    Most collections of interest (e.g., `dict`, `list`, `set`) are containers.
    '''
    return isinstance(obj, Container)


def is_container_nonstr(obj: object) -> bool:
    '''
    `True` only if the passed object is a **non-string container** (i.e.,
    implements the `__contains__` special method _and_ is not a string).
    '''
    return is_container(obj) and not is_str(obj)

# ....................{ TESTERS ~ collection : iterable    }....................
def is_iterable(obj: object) -> bool:
    '''
    `True` only if the passed object is an **iterable**.

    Iterables are objects capable of returning their members one at a time.
    Equivalently, iterables implement the abstract base class
    `collections.Iterable` and hence define the `__iter__()` method.
    '''
    return isinstance(obj, Iterable)


def is_iterable_nonstr(obj: object) -> bool:
    '''
    `True` only if the passed object is a **non-string iterable** (i.e.,
    implements the abstract base class `collections.Iterable` _and_ is not a
    string).
    '''
    return is_iterable(obj) and not is_str(obj)

# ....................{ TESTERS ~ sequence                 }....................
def is_sequence(obj: object) -> bool:
    '''
    `True` only if the passed object is a **sequence**.

    Sequences are iterables supporting efficient element access via integer
    indices. Equivalently, sequences implement the abstract base class
    `collections.abc.Sequence` and hence define the `__getitem__()` and
    `__len__()` methods (among numerous others).

    While all sequences are iterables, not all iterables are sequences.
    Generally speaking, sequences correspond to the proper subset of iterables
    whose elements are ordered. `dict` and `OrderedDict` are the canonical
    examples. `dict` implements `collections.Iterable` but _not_
    `collections.abc.Sequence`, due to _not_ supporting integer index-based
    lookup; `OrderedDict` implements both, due to supporting such lookup.
    '''

    return isinstance(obj, Sequence)


def is_sequence_nonstr(obj: object) -> bool:
    '''
    `True` only if the passed object is a **non-string sequence** (i.e.,
    implements the abstract base class `collections.abc.Sequence` _and_ is not
    a string).

    For generality, this functions returns `True` for both pure-Python
    non-string sequences _and_ non-Pythonic Fortran-based `numpy` arrays and
    matrices (which fail to subclass the `collections.abc.Sequence` API despite
    implementing all methods defined by that subclass).
    '''

    return (
        # Is this a pure-Python non-string sequence?
        (is_sequence(obj) and not is_str(obj)) or
        # Is this a non-Pythonic Fortran-based numpy array or matrix?
        is_numpy_array(obj)
    )


def is_sequence_nonstr_nonempty(obj: object) -> bool:
    '''
    `True` only if the passed object is a **nonempty non-string sequence**
    (i.e., implements the abstract base class `collections.abc.Sequence`, is
    not a string, _and_ contains at least one element).

    '''
    return is_sequence_nonstr(obj) and len(obj)

# ....................{ TESTERS ~ enum                     }....................
def is_enum(obj: object) -> bool:
    '''
    `True` only if the passed object is an enumeration.
    '''
    return isinstance(obj, EnumMeta)


def is_in_enum(obj: object, enum: Enum) -> bool:
    '''
    `True` only if the passed object is in the passed enumeration.

    While trivial, this tester is provided for orthogonality with the
    `assert_not_in_enum()` function.
    '''

    assert is_enum(enum), assert_not_enum(enum)
    return obj in enum

# ....................{ TESTERS ~ error                    }....................
def is_exception(obj: object) -> bool:
    '''
    `True` only if the passed object is an **exception**.
    '''

    return isinstance(obj, Exception)

# ....................{ TESTERS ~ betse : cli              }....................
def is_cli_subcommand(obj: object) -> bool:
    '''
    `True` only if the passed object is an instance of the BETSE-specific
    `CLISubcommand` class.
    '''

    # Avoid circular import dependencies.
    from betse.cli.subcommand import CLISubcommand
    return isinstance(obj, CLISubcommand)

# ....................{ TESTERS ~ betse : science          }....................
def is_cells(obj: object) -> bool:
    '''
    `True` only if the passed object is an instance of the BETSE-specific
    `Cells` class.
    '''

    # Avoid circular import dependencies.
    from betse.science.cells import Cells
    return isinstance(obj, Cells)


def is_parameters(obj: object) -> bool:
    '''
    `True` only if the passed object is an instance of the BETSE-specific
    `Parameters` class.
    '''

    # Avoid circular import dependencies.
    from betse.science.parameters import Parameters
    return isinstance(obj, Parameters)


def is_simulator(obj: object) -> bool:
    '''
    `True` only if the passed object is an instance of the BETSE-specific
    `Simulator` class.
    '''

    # Avoid circular import dependencies.
    from betse.science.sim import Simulator
    return isinstance(obj, Simulator)

def is_simrunner(obj: object) -> bool:
    '''
    `True` only if the passed object is an instance of the BETSE-specific
    `SimRunner` class.
    '''
    from betse.science.simrunner import SimRunner
    return isinstance(obj, SimRunner)


def is_tissue_picker(obj: object) -> bool:
    '''
    `True` only if the passed object is an instance of the BETSE-specific
    `TissuePicker` class.
    '''

    # Avoid circular import dependencies.
    from betse.science.tissue.picker import TissuePicker
    return isinstance(obj, TissuePicker)

# ....................{ TESTERS ~ lib : argparse           }....................
def is_arg_parser(obj: object) -> bool:
    '''
    `True` only if the passed object is an argument parser.
    '''

    # Avoid importing third-party packages at the top level.
    from argparse import ArgumentParser
    return isinstance(obj, ArgumentParser)

# ....................{ TESTERS ~ lib : matplotlib         }....................
#FIXME: Rename "_matplotlib_" everywhere below to merely "_mpl_".
def is_matplotlib_collection(obj: object) -> bool:
    '''
    `True` only if the passed object is a Matplotlib collection.
    '''

    # Avoid importing third-party packages at the top level.
    from matplotlib.collections import Collection
    return isinstance(obj, Collection)


def is_matplotlib_colormap(obj: object) -> bool:
    '''
    `True` only if the passed object is a Matplotlib colormap.
    '''

    # Avoid importing third-party packages at the top level.
    from matplotlib.colors import Colormap
    return isinstance(obj, Colormap)


def is_matplotlib_mappable(obj: object) -> bool:
    '''
    `True` only if the passed object is a Matplotlib mappable.
    '''

    # Avoid importing third-party packages at the top level.
    from matplotlib.cm import ScalarMappable
    return isinstance(obj, ScalarMappable)


def is_matplotlib_polycollection(obj: object) -> bool:
    '''
    `True` only if the passed object is a Matplotlib **polygon collection**
    (i.e., an instance of the `PolyCollection` class).
    '''

    # Avoid importing third-party packages at the top level.
    from matplotlib.collections import PolyCollection
    return isinstance(obj, PolyCollection)


def is_matplotlib_streamplot(obj: object) -> bool:
    '''
    `True` only if the passed object is a Matplotlib **streamplot** (i.e., an
    object returned by the `matplotlib.plt.streamplot()` function).
    '''

    # Avoid importing third-party packages at the top level.
    from matplotlib.streamplot import StreamplotSet
    return isinstance(obj, StreamplotSet)


def is_matplotlib_trimesh(obj: object) -> bool:
    '''
    `True` only if the passed object is a Matplotlib **triangle mesh** (i.e., an
    instance of the `matplotlib.collections.TriMesh` class).
    '''

    # Avoid importing third-party packages at the top level.
    from matplotlib.collections import TriMesh
    return isinstance(obj, TriMesh)

# ....................{ TESTERS ~ lib : numpy              }....................
#FIXME: Rename to is_np_array().
def is_numpy_array(obj: object) -> bool:
    '''
    `True` only if the passed object is a **Numpy array or matrix** (i.e.,
    instance of the `numpy.ndarray` superclass).
    '''

    # Avoid importing third-party packages at the top level.
    from numpy import ndarray
    return isinstance(obj, ndarray)

# ....................{ TESTERS ~ lib : py                 }....................
def is_py_path_local(obj: object) -> bool:
    '''
    `True` only if the passed object is a `py.path.local` instance.
    '''

    # Avoid importing third-party packages at the top level.
    from py._path.local import LocalPath
    return isinstance(obj, LocalPath)

# ....................{ TESTERS ~ lib : numpy              }....................
def is_pytest_fixture(obj: object) -> bool:
    '''
    `True` only if the passed object is a `py.test` fixture.
    '''

    # Avoid importing third-party packages at the top level.
    from _pytest.python import FixtureFunctionMarker
    return isinstance(obj, FixtureFunctionMarker)

# ....................{ TESTERS ~ numeric                  }....................
def is_int(obj: object) -> bool:
    '''
    `True` only if the passed object is an integer.
    '''
    return isinstance(obj, int)


def is_int_ge(obj: object, ge: int) -> bool:
    '''
    `True` only if the passed object is an integer greater than or equal to the
    second passed integer.
    '''
    assert is_int(ge), assert_not_int(ge)
    return is_int(obj) and obj >= ge


def is_int_gt(obj: object, gt: int) -> bool:
    '''
    `True` only if the passed object is an integer strictly greater than the
    second passed integer.
    '''
    assert is_int(gt), assert_not_int(gt)
    return is_int(obj) and obj > gt


def is_int_positive(obj: object) -> bool:
    '''
    `True` only if the passed object is a **positive integer** (i.e., is
    strictly greater than 0).
    '''
    return is_int_gt(obj, 0)


def is_numeric(obj: object) -> bool:
    '''
    `True` only if the passed object is **numeric** (i.e., instance of either
    the `int` or `float` types).
    '''
    return isinstance(obj, (int, float))

# ....................{ TESTERS ~ str                      }....................
def is_str(obj: object) -> bool:
    '''
    `True` only if the passed object is a **string** (i.e., instance of the
    `str` class).
    '''

    return isinstance(obj, str)


def is_str_or_none(obj: object) -> bool:
    '''
    `True` only if the passed object is either a string _or_ `None`.
    '''

    return isinstance(obj, str) or obj is None


def is_str_nonempty(obj: object) -> bool:
    '''
    `True` only if the passed object is a **nonempty string* (i.e., string
    comprising one or more characters and hence _not_ the empty string).
    '''

    return is_str(obj) and len(obj)


def is_str_nonempty_or_none(obj: object) -> bool:
    '''
    `True` only if the passed object is either a nonempty string _or_ `None`.
    '''

    return is_str_nonempty(obj) or obj is None

# ....................{ ASSERTERS                          }....................
def assert_not_bool(obj: object) -> str:
    '''
    String asserting the passed object to _not_ be boolean.
    '''
    return '"{}" not boolean (i.e., neither "True" nor "False").'.format(
        trim(obj))


def assert_not_char(obj: object) -> str:
    '''
    String asserting the passed object to _not_ be character.
    '''
    return '"{}" not a character (i.e., string of length 1).'.format(trim(obj))


def assert_not_class(obj: object) -> str:
    '''
    String asserting the passed object to _not_ be a class.
    '''
    return '"{}" not a class.'.format(trim(obj))


def assert_not_nonnone(label: str) -> str:
    '''
    String asserting an arbitrary object labeled by the passed human-readable
    name to be `None`.
    '''
    return '{} is "None".'.format(label.capitalize())

# ....................{ TESTERS ~ callable                 }....................
def assert_not_callable(obj: object) -> bool:
    '''
    String asserting the passed object to _not_ be callable.
    '''
    return '"{}" not callable.'.format(trim(obj))

# ....................{ ASSERTERS ~ collection             }....................
def assert_not_mapping(obj: object) -> str:
    '''
    String asserting the passed object to _not_ be a mapping.
    '''
    return '"{}" not a mapping (e.g., "dict", "OrderedDict").'.format(trim(obj))

# ....................{ ASSERTERS ~ collection : container }....................
def assert_not_container_nonstr(obj: object) -> str:
    '''
    String asserting the passed object to _not_ be a non-string container.
    '''
    return '"{}" not a non-string container (e.g., "dict", "set").'.format(
        trim(obj))

# ....................{ ASSERTERS ~ collection : iterable  }....................
def assert_not_iterable_nonstr(obj: object) -> str:
    '''
    String asserting the passed object to _not_ be a non-string iterable.
    '''
    return '"{}" not a non-string iterable (e.g., "dict", "list").'.format(
        trim(obj))


def assert_not_iterable_nonstr_nonempty(obj: object, label: str) -> str:
    '''
    String asserting the passed object labeled by the passed human-readable
    name to _not_ be a nonempty non-string iterable.
    '''
    return assert_not_iterable_nonstr(obj) if not is_iterable_nonstr(
        obj) else '{} empty.'.format(label.capitalize())

# ....................{ ASSERTERS ~ collection : sequence  }....................
def assert_not_sequence_nonstr(obj: object) -> str:
    '''
    String asserting the passed object to _not_ be a non-string sequence.
    '''
    return '"{}" not a non-string sequence (e.g., "list").'.format(trim(obj))


def assert_not_sequence_nonstr_nonempty(obj: object, label: str) -> str:
    '''
    String asserting the passed object labeled by the passed human-readable
    name to _not_ be a nonempty non-string sequence.
    '''
    return assert_not_sequence_nonstr(obj) if not is_sequence_nonstr(
        obj) else '{} empty.'.format(label.capitalize())

# ....................{ ASSERTERS ~ enum                   }....................
def assert_not_enum(obj: object) -> str:
    '''
    String asserting the passed object to _not_ be an enumeration.
    '''
    return '"{}" not an enumeration "{}".'.format(trim(obj))


def assert_not_in_enum(obj: object, enum: Enum) -> str:
    '''
    String asserting the passed object to _not_ be in the passed enumeration.
    '''
    return '"{}" not in the enumeration "{}".'.format(trim(obj), trim(enum))

# ....................{ ASSERTERS ~ error                  }....................
def assert_not_exception(obj: object) -> str:
    '''
    String asserting the passed object to _not_ be an exception.
    '''
    return '"{}" not an exception.'.format(trim(obj))

# ....................{ ASSERTERS ~ betse : cli            }....................
def assert_not_cli_subcommand(obj: object) -> bool:
    '''
    String asserting the passed object is an instance of the BETSE-specific
    `CLISubcommand` class.
    '''
    return '"{}" not a CLI subcommand.'.format(trim(obj))

# ....................{ ASSERTERS ~ betse : science        }....................
def assert_not_cells(obj: object) -> str:
    '''
    String asserting the passed object to _not_ be an instance of the BETSE-
    specific `Cells` class.
    '''
    return '"{}" not a "Cells" instance.'.format(trim(obj))


def assert_not_parameters(obj: object) -> str:
    '''
    String asserting the passed object to _not_ be an instance of the BETSE-
    specific `Parameters` class.
    '''
    return '"{}" not a "Parameters" instance.'.format(trim(obj))


def assert_not_simulator(obj: object) -> str:
    '''
    String asserting the passed object to _not_ be an instance of the BETSE-
    specific `Simulator` class.
    '''
    return '"{}" not a "Simulator" instance.'.format(trim(obj))


def assert_not_tissue_picker(obj: object) -> str:
    '''
    String asserting the passed object to _not_ be an instance of the BETSE-
    specific `TissuePicker` class.
    '''
    return '"{}" not a "TissuePicker" instance.'.format(trim(obj))

# ....................{ ASSERTERS ~ lib : argparse         }....................
def assert_not_arg_parser(obj: object) -> bool:
    '''
    String asserting the passed object to _not_ be an argument parser.
    '''
    return '"{}" not an argument parser.'.format(trim(obj))

# ....................{ ASSERTERS ~ lib : matplotlib       }....................
def assert_not_matplotlib_collection(obj: object) -> bool:
    '''
    String asserting the passed object to _not_ be a Matplotlib collection.
    '''
    return '"{}" not a Matplotlib collection.'.format(trim(obj))


def assert_not_matplotlib_colormap(obj: object) -> bool:
    '''
    String asserting the passed object to _not_ be a Matplotlib colormap.
    '''
    return '"{}" not a Matplotlib colormap.'.format(trim(obj))


def assert_not_matplotlib_mappable(obj: object) -> bool:
    '''
    String asserting the passed object to _not_ be a Matplotlib mappable.
    '''
    return '"{}" not a Matplotlib mappable.'.format(trim(obj))


def assert_not_matplotlib_polycollection(obj: object) -> bool:
    '''
    String asserting the passed object to _not_ be a Matplotlib polygon
    collection.
    '''
    return '"{}" not a Matplotlib polygon collection.'.format(trim(obj))


def assert_not_matplotlib_streamplot(obj: object) -> bool:
    '''
    String asserting the passed object to _not_ be a Matplotlib streamplot.
    '''
    return '"{}" not a Matplotlib streamplot.'.format(trim(obj))


def assert_not_matplotlib_trimesh(obj: object) -> bool:
    '''
    String asserting the passed object to _not_ be a Matplotlib triangle mesh.
    '''
    return '"{}" not a Matplotlib triangle mesh.'.format(trim(obj))

# ....................{ ASSERTERS ~ lib : numpy            }....................
def assert_not_numpy_array(obj: object) -> bool:
    '''
    String asserting the passed object to _not_ be a Numpy array or matrix.
    '''
    return '"{}" not a NumPy array or matrix.'.format(trim(obj))

# ....................{ ASSERTERS ~ lib : py               }....................
def assert_not_py_path_local(obj: object) -> bool:
    '''
    String asserting the passed object to _not_ be a `py.path.local` instance.
    '''
    return '"{}" not a "py.path.local" instance.'.format(trim(obj))

# ....................{ ASSERTERS ~ lib : pytest           }....................
def assert_not_pytest_fixture(obj: object) -> bool:
    '''
    String asserting the passed object to _not_ be a `py.test` fixture.
    '''
    return '"{}" not a py.test fixture.'.format(trim(obj))

# ....................{ ASSERTERS ~ numeric                }....................
def assert_not_int(obj: object) -> str:
    '''
    String asserting the passed object to _not_ be an integer.
    '''
    return '"{}" not an integer.'.format(trim(obj))


def assert_not_int_ge(obj: object, ge: int) -> str:
    '''
    String asserting the passed object to _not_ be an integer greater than or
    equal to the second passed integer.
    '''
    return '"{}" not an integer or not >= {}.'.format(trim(obj), ge)


def assert_not_numeric(obj: object) -> str:
    '''
    String asserting the passed object to _not_ be numeric.
    '''
    return '"{}" not numeric (i.e., neither an integer nor float).'.format(
        trim(obj))

# ....................{ ASSERTERS ~ str                    }....................
def assert_not_str(obj: object) -> str:
    '''
    String asserting the passed object to _not_ be a string.
    '''

    return '"{}" not a string.'.format(trim(obj))


def assert_not_str_or_none(obj: object) -> str:
    '''
    String asserting the passed object to be neither a string _nor_ `None`.
    '''

    return '"{}" not a string or "None".'.format(trim(obj))


def assert_not_str_nonempty(obj: object, label: str) -> str:
    '''
    String asserting the passed object labeled by the passed human-readable
    name to _not_ be a nonempty string.
    '''

    if not is_str(obj):
        return assert_not_str(obj)
    else:
        return '{} empty.'.format(label.capitalize())


def assert_not_str_nonempty_or_none(obj: object, label: str) -> str:
    '''
    String asserting the passed object labeled by the passed human-readable
    name to be neither a nonempty string _nor_ `None`.
    '''

    if not is_str_or_none(obj):
        return assert_not_str_or_none(obj)
    else:
        return assert_not_str_nonempty(obj, label)
