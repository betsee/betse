#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

r'''
Core **type** (i.e., class) functionality, enumerating a variety of core types
and :func:`instanceof`\ -friendly tuples of such types *and* the pivotal
:func:`type_check` decorator validating callable parameters to be of such types.
'''

# ....................{ IMPORTS                            }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: To raise human-readable exceptions on missing mandatory dependencies
# *AND* avoid non-halting recursive imports when imported at the top-level
# of other modules in the "betse.util" package, this module may import *ONLY*
# from stock Python packages. (By definition, this excludes both application and
# third-party packages.)
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

import inspect, re
from argparse import ArgumentParser, _SubParsersAction
from collections.abc import (
    Container, Hashable, Iterable, Iterator, Mapping, MutableMapping, Sequence,
    Set, Sized,
)
from enum import Enum, EnumMeta
from functools import partial, wraps
from io import IOBase
from inspect import Parameter, Signature
from weakref import CallableProxyType, ProxyType

# Import the following types as is into the namespace of this submodule,
# permitting callers to reference these types conveniently. Since the
# nomenclature of these types is already consistent with that used by types
# declared below (namely, both camelcase and suffixed by "Type"), these types
# are used as is rather than aliased to synonymous types below.
#
# Note that the "LambdaType" is intentionally *NOT* imported. Why? Because that
# type is exactly synonymous with "FunctionType", implying lambdas are
# indistinguishable from functions. To curtail confusion elsewhere and, in
# particular, to prevent functions from being misidentified as lambdas, all
# lambdas are currently misidentified as functions. This is the lesser of
# multiple evils, we're afraid.
from types import (
    BuiltinFunctionType,
    BuiltinMethodType,
    FunctionType,
    GeneratorType,
    MethodType,
    ModuleType,
)

# Silence IDE warnings concerning locally unused attributes. Move along, folks.
if False: (wraps, GeneratorType, ModuleType,)

# ....................{ TYPES                              }....................
ClassType = type
'''
Type of all types.
'''


FileType = IOBase
'''
Abstract base class implemented by *all* **file-like objects** (i.e., objects
implementing the standard ``read()`` and ``write()`` methods).

This class is a synonym of the `io.IOBase` class, provided merely as a
convenience to callers preferring to avoid importing that class.
'''


NoneType = type(None)
'''
Type of the singleton ``None`` object.

Curiously, although the type of the ``None`` object is a class object whose
``__name__`` attribute is ``NoneType``, there exists no globally accessible
class by that name. To circumvents this obvious oversight, this global globally
exposes this class.

This class is principally useful for annotating both:

* Callable parameters accepting ``None`` as a valid value.
* Callables returning ``None`` as a valid value.

Note that, for obscure and uninteresting reasons, the standard :mod:`types`
module defined the same type with the same name under Python 2.x but _not_ 3.x.
Depressingly, this type must now be manually redefined everywhere.
'''

# ....................{ TYPES ~ arg                        }....................
ArgParserType = ArgumentParser
'''
Type of argument parsers parsing all command-line arguments for either top-level
commands *or* subcommands of those commands.

This class is a synonym of the :class:`argparse.ArgumentParser` class,
permitting callers to avoid importing that class.
'''


ArgSubparsersType = _SubParsersAction
'''
Type of argument subparser containers parsing subcommands for parent argument
parsers parsing either top-level commands *or* subcommands of those commands.

This class is a synonym of the :class:`argparse._SubParsersAction` class,
permitting callers to avoid importing that private class.
'''

# ....................{ TYPES ~ callable                   }....................
CallablePartialType = partial
'''
Type of all **partial callables** (i.e., callables dynamically produced by the
function-like :class:`partial` class).
'''


PropertyType = property
'''
Type of all **property methods** (i.e., methods decorated by the builtin
:class:`property` class decorator).
'''

# ....................{ TYPES ~ container                  }....................
ContainerType = Container
'''
Abstract interface implemented by all **containers** (i.e., objects
implementing the standard `__contains__()` method internally called by the
`in` operator).

This class is a synonym of the :class:`collections.abc.Container` class,
permitting callers to avoid importing that class.
'''


IteratorType = Iterator
'''
Abstract interface implemented by all **iterators** (i.e., objects implementing
the standard `__iter__()` and `__next__()` methods, typically iterating over an
associated container).

This class is a synonym of the `collections.abc.Iterator` class, provided
merely as a convenience to callers preferring to avoid importing that class.
'''


SetType = Set
'''
Abstract interface implemented by all set-like objects.

This class is a synonym of the :class:`collections.abc.Set`, provided merely as
a convenience to callers preferring to avoid importing that class.
'''


SizedType = Sized
'''
Abstract interface implemented by all containers defining the special
`__len__()` method internally called by the :func:`len` builtin.

This class is a synonym of the `collections.abc.Sized` class, provided merely
as a convenience to callers preferring to avoid importing that class.
'''

# ....................{ TYPES ~ container : mapping        }....................
HashableType = Hashable
'''
Abstract interface implemented by all **hashables** (i.e., objects implementing
the standard ``__hash__()`` method required by all dictionary keys).

This class is a synonym of the `collections.abc.Hashable` class, provided
merely as a convenience to callers preferring to avoid importing that class.
'''


MappingType = Mapping
'''
Abstract interface implemented by all dictionary-like objects, both mutable and
immutable.

This class is a synonym of the `collections.abc.Mapping` class, provided merely
as a convenience to callers preferring to avoid importing that class.
'''


MappingMutableType = MutableMapping
'''
Abstract interface implemented by all mutable dictionary-like objects.

This class is a synonym of the `collections.abc.MutableMapping` class, provided
merely as a convenience to callers preferring to avoid importing that class.
'''

# ....................{ TYPES ~ enum                       }....................
# Enumeration types sufficiently obscure to warrant formalization here.

EnumType = EnumMeta
'''
Metaclass of all **enumeration types** (i.e., classes containing all enumeration
members comprising those enumerations).

This class is a synonym of the :class:`enum.EnumMeta` class, permitting callers
to avoid importing that class.

Motivation
----------
This type is widely used throughout the codebase to validate callable parameters
to be enumerations. In recognition of its popularity, this type is intentionally
named ``EnumType`` rather than ``EnumMetaType``. While the latter *would*
technically be less ambiguous, the former has the advantage of inviting
correctness throughout the codebase -- a less abundant resource.

Why? Because *all* enumeration types are instances of this type rather than the
:class:`Enum` class despite being superficially defined as instances of the
:class:`Enum` class. Thanks to metaclass abuse, enumeration types do *not*
adhere to standard Pythonic semantics. Notably, the following non-standard
invariants hold across *all* enumerations:

    >>> from betse.util.type.types import EnumType, EnumClassType, ClassType
    >>> enum_type = EnumClassType(
    ...     'Gyre', ('The', 'falcon', 'cannot', 'hear', 'the', 'falconer'))
    >>> isinstance(enum_type, EnumClassType)
    False
    >>> issubclass(enum_type, EnumClassType)
    True
    >>> isinstance(enum_type, EnumType)
    True
    >>> isinstance(enum_type, ClassType)
    True
'''


EnumClassType = Enum
'''
Abstract base class of all **enumeration types** (i.e., classes containing all
enumeration members comprising those enumerations).

This class is a synonym of the :class:`enum.Enum` class, permitting callers to
avoid importing that class.

See Also
----------
:class:`EnumType`
    Further details.
'''


EnumMemberType = Enum
'''
Abstract base class implemented by all **enumeration members** (i.e.,
alternative choices comprising their parent enumerations).

This class is a synonym of the :class:`enum.Enum` class, permitting callers
to avoid importing that class.

Caveats
----------
When type checking callable parameters, this class should *only* be referenced
where the callable permissively accepts any enumeration member type rather than
a specific enumeration member type. In the latter case, that type is simply that
enumeration's type and should be directly referenced as such: e.g.,

    >>> from betse.util.type.enums import make_enum
    >>> from betse.util.type.types import type_check
    >>> EndymionType = make_enum(
    ...     class_name='EndymionType', member_names=('BEAUTY', 'JOY',))
    >>> @type_check
    ... def our_feet_were_soft_in_flowers(superlative: EndymionType) -> str:
    ...     return str(superlative).lower()
'''

# ....................{ TUPLES                             }....................
CallableTypes = (
    BuiltinFunctionType,
    BuiltinMethodType,
    FunctionType,
    MethodType,
)
'''
Tuple of all **callable classes** (i.e., classes whose instances are callable
objects, including both built-in and user-defined functions, lambdas, and
methods).
'''


CallableOrStrTypes = CallableTypes + (str,)
'''
Tuple of all callable classes *and* the string type.
'''


FunctionTypes = (BuiltinFunctionType, FunctionType,)
'''
Tuple of all **function classes** (i.e., classes whose instances are either
built-in or user-defined functions).
'''


MethodTypes = (BuiltinMethodType, MethodType,)
'''
Tuple of all **method classes** (i.e., classes whose instances are either
built-in or user-defined methods).
'''


ModuleOrStrTypes = (str, ModuleType)
'''
Tuple of both the module *and* string type.
'''


NumericTypes = (float, int)
'''
Tuple of all **numeric classes** (i.e., classes whose instances are single
scalar numbers excluding complex numbers, which typically require
special-purpose handling).

This tuple contains classes matching both integer and real number types.

Caveats
----------
For obscure reasons, this tuple intentionally lists the :class:`float` class
*BEFORE* the :class:`int` class. (Downstream BETSEE requirements coerce
GUI-based numeric string values into numbers by casting these strings into
instances of the first item of this tuple. Reversing the order of these items in
this tuple would adversely strip the decimal portion from real number strings.)
'''


TestableTypes = (ClassType, tuple)
'''
Tuple of all **testable types** (i.e., types suitable for use as the second
parameter passed to the :func:`isinstance` and :func:`issubclass` builtins).
'''


WeakRefProxyTypes = (CallableProxyType, ProxyType)
'''
Tuple of all **weak reference proxy classes** (i.e., classes whose instances
are weak references to other instances masquerading as those instances).

This tuple contains classes matching both callable and uncallable weak
reference proxies.
'''

# ....................{ TUPLES : init                      }....................
# These tuples are declared here for documentation purposes but defined below.

NumpyArrayType = None
'''
Type of Numpy arrays if :class:`numpy` is importable *or* ``None`` otherwise.

This class is a synonym of the :class:`numpy.ndarray` class, permitting callers
to avoid importing that class.
'''


IterableTypes = None
'''
Tuple of all container base classes conforming to (but *not* necessarily
subclassing) the canonical :class:`collections.abc.Iterable` API.

See Also
----------
:class:`SequenceTypes`
    Further details.
'''


SequenceTypes = None
'''
Tuple of all container base classes conforming to (but *not* necessarily
subclassing) the canonical :class:`collections.abc.Sequence` API.

Sequences are iterables supporting efficient element access via integer indices.
Most sequences implement the :class:`collections.abc.Sequence` abstract base
class, including the concrete :class:`str` string class. All sequences define
the special ``__getitem__()`` and ``__len__()`` methods, amongst numerous
others.

While all sequences are iterables, not all iterables are sequences. Generally
speaking, sequences correspond to the proper subset of iterables whose elements
are ordered. :class:`dict` and :class:`OrderedDict` are the canonical examples.
:class:`dict` implements :class:`collections.abc.Iterable` but *not*
:class:`collections.abc.Sequence`, due to failing to support integer index-based
lookup; :class:`OrderedDict` implements both, due to supporting such lookup.

For generality, this tuple contains classes matching both pure-Python sequences
*and* non-Pythonic Fortran-based Numpy arrays and matrices -- which fail to
subclass :class:`collections.abc.Sequence` despite implementing the entirety of
that that API.
'''


# Conditionally add sequence types to previously declared tuples.
#
# If Numpy is available, add both core APIs and the Numpy array type (which
# fails to subclass these APIs). Although Numpy is a mandatory dependency, this
# submodule is typically imported quite early in program startup, implying the
# importability of *ANY* dependency (mandatory or not) at the top level of this
# submodule to still be in question. Since subsequent logic in program startup
# is guaranteed to raise human-readable exceptions for missing dependencies,
# this error is silently ignored here.
try:
    from numpy import ndarray

    NumpyArrayType = ndarray
    IterableTypes = (Iterable, NumpyArrayType)
    SequenceTypes = (Sequence, NumpyArrayType)
# Else, Numpy is unavailable. Add only core APIs.
except:
    NumpyArrayType = None
    IterableTypes = (Iterable,)
    SequenceTypes = (Sequence,)

# ....................{ TUPLES : post-init                 }....................
# Tuples of types assuming the above initialization to have been performed.

MappingOrSequenceTypes = (MappingType,) + SequenceTypes
'''
Tuple of all container base classes conforming to (but *not* necessarily
subclassing) the canonical :class:`Mapping` *or* :class:`Sequence` APIs.
'''


NumericOrIterableTypes = NumericTypes + IterableTypes
'''
Tuple of all numeric types *and* all container base classes conforming to (but
*not* necessarily subclassing) the canonical :class:`Iterable` API.
'''


NumericOrSequenceTypes = NumericTypes + SequenceTypes
'''
Tuple of all numeric types *and* all container base classes conforming to (but
*not* necessarily subclassing) the canonical :class:`Sequence` API.
'''

# ....................{ TUPLES : none                      }....................
# Tuples of types containing at least the type of the singleton "None" object.

NoneTypes = (NoneType,)
'''
Tuple of only the type of the singleton ``None`` object.

This tuple is principally intended for use in efficiently constructing other
tuples of types containing this type.
'''


BoolOrNoneTypes = (bool, NoneType)
'''
Tuple of both the boolean type *and* that of the singleton ``None`` object.
'''


CallableOrNoneTypes = CallableTypes + NoneTypes
'''
Tuple of all callable classes *and* the type of the singleton ``None`` object.
'''


ClassOrNoneTypes = (ClassType, NoneType)
'''
Tuple of the type of all types *and* that of the singleton ``None`` object.
'''


IntOrNoneTypes = (int, NoneType)
'''
Tuple of both the integer type *and* that of the singleton ``None`` object.
'''


IterableOrNoneTypes = IterableTypes + NoneTypes
'''
Tuple of all container base classes conforming to (but _not_ necessarily
subclassing) the canonical :class:`Iterable` API as well as the type of the
singleton ``None`` object.
'''


MappingOrNoneTypes = (MappingType,) + NoneTypes
'''
Tuple of all container base classes conforming to (but *not* necessarily
subclassing) the canonical :class:`Mapping` API as well as the type of the
singleton ``None`` object.
'''


MappingOrSequenceOrNoneTypes = MappingOrSequenceTypes + NoneTypes
'''
Tuple of all container base classes conforming to (but *not* necessarily
subclassing) the canonical :class:`Mapping` *or* :class:`Sequence` APIs as well
as the type of the singleton ``None`` object.
'''


NumericOrSequenceOrNoneTypes = NumericOrSequenceTypes + NoneTypes
'''
Tuple of all numeric types, all container base classes conforming to (but *not*
necessarily subclassing) the canonical :class:`int`, :class:`float`, *or*
:class:`Sequence` APIs as well as the type of the singletone ``None`` object.
'''


SequenceOrNoneTypes = SequenceTypes + NoneTypes
'''
Tuple of all container base classes conforming to (but *not* necessarily
subclassing) the canonical :class:`Sequence` API as well as the type of the
singleton ``None`` object.
'''


SetOrNoneTypes = (SetType, NoneType)
'''
Tuple of both the set type *and* the type of the singleton ``None`` object.
'''


NumericOrNoneTypes = NumericTypes + NoneTypes
'''
Tuple of all numeric types *and* the type of the singleton `None` object.
'''


StrOrNoneTypes = (str, NoneType)
'''
Tuple of both the string type *and* the type of the singleton ``None`` object.
'''


TestableOrNoneTypes = TestableTypes + NoneTypes
'''
Tuple of all testable types *and* the type of the singleton ``None`` object.
'''

# ....................{ TUPLES ~ regex                     }....................
# Yes, this is the only reliable means of obtaining the type of compiled regular
# expressions. For unknown reasons presumably concerning the archaic nature of
# Python's regular expression support, this type is *NOT* publicly exposed.
# While the private "re._pattern_type" attribute does technically provide this
# type, it does so in a private and hence inherently non-portable manner.
RegexCompiledType = type(re.compile(''))
'''
Type of all compiled regular expressions.
'''


# Yes, this type is required for type validation et the module scope elsewhere.
# Yes, this is the most time-efficient means of obtaining this type. No, this
# type is *NOT* directly importable. Although this type's classname is
# published to be "_sre.SRE_Match", the "_sre" C extension provides no such
# type for pure-Python importation.
RegexMatchType = type(re.match(r'', ''))
'''
Type of all **regular expression match objects** (i.e., objects returned by the
:func:`re.match` function).
'''


RegexTypes = (str, RegexCompiledType)
'''
Tuple of all **regular expression-like types** (i.e., types either defining
regular expressions or losslessly convertible to such types, typically accepted
by functions in the :mod:`betse.util.type.regexes` submodule).
'''


RegexMatchOrNoneTypes = (RegexMatchType, NoneType)
'''
Tuple of both the regular expression match object type *and* the type of the
singleton ``None`` object.
'''

# ....................{ SETS : private                     }....................
#FIXME: Type-check variadic keyword arguments as well.
_PARAMETER_KIND_IGNORED = {
    Parameter.POSITIONAL_ONLY, Parameter.VAR_KEYWORD,
}
'''
Set of all :attr:`inspect.Parameter.kind` constants to be ignored during
annotation-based type checking in the :func:`type_check` decorator.

This includes:

* Constants specific to variadic parameters (e.g., `*args`, `**kwargs`).
  Variadic parameters cannot be annotated and hence cannot be type checked.
* Constants specific to positional-only parameters, which apply to non-pure-
  Python callables (e.g., defined by C extensions). The :func:`type_check`
  decorator applies _only_ to pure-Python callables, which provide no syntactic
  means for specifying positional-only parameters.
'''


_RETURN_ANNOTATION_IGNORED = {Signature.empty, None}
'''
Set of all annotations for return values to be ignored during annotation-
based type checking in the :func:`type_check` decorator.

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
        both all annotated parameters passed to this callable *and* the
        annotated value returned by this callable if any.

        This decorator performs rudimentary type checking based on Python 3.x
        function annotations, as officially documented by PEP 484 ("Type
        Hints"). While PEP 484 supports arbitrarily complex type composition,
        this decorator requires _all_ parameter and return value annotations to
        be either:

        * Classes (e.g., `int`, `OrderedDict`).
        * Tuples of classes (e.g., `(int, OrderedDict)`).

        If optimizations are enabled by the active Python interpreter (e.g.,
        due to option `-O` passed to this interpreter), this decorator reduces
        to a noop.

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
        for func_arg_index, func_arg in enumerate(
            func_sig.parameters.values()):
            # If this callable redefines a parameter initialized to a default
            # value by this wrapper, raise an exception. Permitting this
            # unlikely edge case would permit unsuspecting users to
            # "accidentally" override these defaults.
            if func_arg.name == '__beartype_func':
                raise NameError(
                    'Parameter {} reserved for use by @type_check.'.format(
                        func_arg.name))

            # Value of this annotation.
            func_arg_annotation = func_arg.annotation

            # If this parameter is both annotated and non-ignorable for
            # purposes of type checking, type check this parameter.
            if (func_arg_annotation is not Parameter.empty and
                func_arg.kind not in _PARAMETER_KIND_IGNORED):
                # Human-readable label describing this annotation.
                func_arg_annotation_label = (
                    '{} parameter "{}" type annotation'.format(
                        func_name, func_arg.name))

                # Validate this annotation.
                _check_type_annotation(
                    annotation=func_arg_annotation,
                    annotation_label=func_arg_annotation_label,)

                # String evaluating to this parameter's annotated type.
                func_arg_type_expr = (
                    '__beartype_func.__annotations__[{!r}]'.format(
                        func_arg.name))

                # String evaluating to this parameter's current value when
                # passed as a keyword.
                func_arg_value_key_expr = 'kwargs[{!r}]'.format(func_arg.name)

                # If this annonation's value is a string, replace this value
                # with the value of the module attribute whose fully-qualified
                # name is this string.
                if isinstance(func_arg_annotation, str):
                    func_body += (
                        _get_code_setting_annotation_to_module_attr(
                            annotation_expr=func_arg_type_expr,
                            annotation_label=func_arg_annotation_label,
                            module_attr_name=func_arg_annotation,))

                # If this parameter is actually a tuple of positional variadic
                # parameters, iteratively check all such parameters.
                if func_arg.kind is Parameter.VAR_POSITIONAL:
                    func_body += '''
    for __beartype_arg in args[{arg_index!r}:]:
        if not isinstance(__beartype_arg, {arg_type_expr}):
            raise TypeError(
                '{func_name} positional variadic parameter '
                '{arg_index} {{}} not a {{!r}}'.format(
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

        # Value of the annotation for this callable's return value.
        func_return_annotation = func_sig.return_annotation

        # If this callable's return value is both annotated and non-ignorable
        # for purposes of type checking, type check this value.
        if func_return_annotation not in _RETURN_ANNOTATION_IGNORED:
            # Human-readable label describing this annotation.
            func_return_annotation_label = (
                '{} return type annotation'.format(func_name))

            # Validate this annotation.
            _check_type_annotation(
                annotation=func_return_annotation,
                annotation_label=func_return_annotation_label,)

            # Strings evaluating to this parameter's annotated type and
            # currently passed value, as above.
            func_return_type_expr = (
                "__beartype_func.__annotations__['return']")

            # If this annonation's value is a string, replace this value with
            # the value of the module attribute whose fully-qualified name is
            # this string.
            if isinstance(func_return_annotation, str):
                func_body += (
                    _get_code_setting_annotation_to_module_attr(
                        annotation_expr=func_return_type_expr,
                        annotation_label=func_return_annotation_label,
                        module_attr_name=func_return_annotation,))

            # Call this callable, type check the returned value, and return
            # this value from this wrapper.
            func_body += '''
    __beartype_return_value = __beartype_func(*args, **kwargs)
    if not isinstance(__beartype_return_value, {return_type}):
        raise TypeError(
            '{func_name} return value {{}} not of {{!r}}'.format(
                trim(__beartype_return_value), {return_type}))
    return __beartype_return_value
'''.format(func_name=func_name, return_type=func_return_type_expr)
        # Else, call this callable and return this value from this wrapper.
        else:
            func_body += '''
    return __beartype_func(*args, **kwargs)
'''

        # Dictionary mapping from local attribute name to value. For
        # efficiency, only attributes required by the body of this wrapper are
        # copied from the current namespace. (See below.)
        local_attrs = {'__beartype_func': func}

        # Dynamically define this wrapper as a closure of this decorator. For
        # obscure and presumably uninteresting reasons, Python fails to locally
        # declare this closure when the locals() dictionary is passed; to
        # capture this closure, a local dictionary must be passed instead.
        #
        # Note that the same result may also be achieved via the compile()
        # builtin and "types.FunctionType" class: e.g.,
        #
        #     func_code = compile(func_body, "<string>", "exec").co_consts[0]
        #     return types.FunctionType(
        #         code=func_code,
        #         globals=globals(),
        #         argdefs=('__beartype_func', func)
        #     )
        #
        # Since doing so is both more verbose and obfuscatory for no tangible
        # gain, the current circumspect approach is preferred.
        # print('\n{} wrapper: {}'.format(func_name, func_body))
        exec(func_body, globals(), local_attrs)

        # Return this wrapper.
        return local_attrs['func_type_checked']

# Else, the active Python interpreter is optimized. In this case, disable type
# checking by reducing this decorator to the identity decorator.
else:
    def type_check(func: CallableTypes) -> CallableTypes:
        return func


# Due to chicken-and-egg issues, this implementation duplicates rather than
# defers to similar functions in the "betse.util.py.modules" submodule.
def _get_code_setting_annotation_to_module_attr(
    annotation_expr: str,
    annotation_label: str,
    module_attr_name: str,
) -> object:
    '''
    Block of valid Python code dynamically replacing the value of the function
    annotation with the first passed name with the value of the module
    attribute with the second passed name.

    Parameters
    ----------
    annotation_expr : str
        Python expression evaluating to the annotation to be replaced.
    annotation_label : str
        Human-readable label describing this annotation, interpolated into
        exceptions raised by this function.
    module_attr_name : str
        Fully-qualified `.`-delimited name of the module attribute to retrieve
        (e.g., `betse.science.sim.Simulator`). This name _must_ contain at
        least one `.`. The substring of this name:
        * Preceding the last `.` signifies the fully-qualified name of the
          module containing this attribute (e.g., `betse.science.sim`).
        * Following the last `.` signifies the unqualified basename of this
          attribute in this module (e.g., `Simulator`).
    '''
    assert isinstance(annotation_expr, str), (
        '"{}" not a string.'.format(annotation_expr))
    assert isinstance(annotation_label, str), (
        '"{}" not a string.'.format(annotation_label))
    assert isinstance(module_attr_name, str), (
        '"{}" not a string.'.format(module_attr_name))

    # Fully-qualified module name and unqualified attribute basename split
    # from the passed fully-qualified module attribute name. (It is good.)
    #
    # Note that, if the passed string erroneously contains no "." delimiters
    # and hence is *NOT* a fully-qualified module attribute name, the
    # str.rsplit() method silently returns a list whose single element is this
    # string as is rather than raising an exception. While such behaviour would
    # typically be bad, the previously called _check_type_annotation() function
    # has already guaranteed this string to contain at least one "." delimiter.
    # Ergo, we're guaranteed good here.
    module_name, attr_name = module_attr_name.rsplit(sep='.', maxsplit=1)

    # Return the desired block of Python code.
    return '''
    # If this annotation's value is still a string, this value has yet to be
    # replaced with the value of the module attribute whose fully-qualified
    # name is this attribute, implying this to be the first call to this
    # callable. In this, perform this replacement, ensuring that no subsequent
    # calls to this callable are required to do so.
    if isinstance({annotation_expr}, str):
        # Attempt to import this attribute from this module, implicitly
        # raising a human-readable "ImportError" exception on failure.
        from {module_name} import {attr_name}

        # Validate this attribute to be either a new-style class or tuple of
        # such classes, preventing this attribute from being yet another
        # string. (The recursion ends here, folks.)
        _check_type_annotation(
            annotation={attr_name},
            annotation_label={annotation_label!r},
            is_str_valid=False,
        )

        # Replace the external copy of this annotation stored with this
        # function's signature by the value of the module attribute whose
        # fully-qualified name is the value of this annotation -- guaranteeing
        # that subsequent access of this annotation via
        # "__beartype_func.__annotations__" accesses the expected type rather
        # than this type's name.
        {annotation_expr} = {attr_name}
    '''.format(
        annotation_expr=annotation_expr,
        annotation_label=annotation_label,
        module_name=module_name,
        attr_name=attr_name,
    )


def _check_type_annotation(
    annotation: object,
    annotation_label: str,
    is_str_valid: bool = True,
) -> None:
    '''
    Validate the passed annotation to be a valid type supported by the
    `@type_check` decorator.

    Parameters
    ----------
    annotation : object
        Annotation to be validated.
    annotation_label : str
        Human-readable label describing this annotation, interpolated into
        exceptions raised by this function.
    is_str_valid : optional[bool]
        ``True`` only if this function accepts string annotations as valid.
        Defaults to `True`. If this boolean is:
        * `True`, this annotation is valid if this annotation's value is either
          a new-style class, tuple of such classes, or string containing one or
          more `.` characters.
        * `False`, this annotation is valid if this annotation's value is
          either a new-style class or tuple of such classes.

    Raises
    ----------
    TypeError
        If this annotation is none of the following:
        * A new-style class.
        * A tuple of new-style classes.
        * A **fully-qualified module attribute name** (i.e., string containing
          at least one `.` delimiter, such that the substring of this string:
          * Preceding the last `.` signifies the fully-qualified name of the
            module containing this attribute (e.g., `betse.science.sim`).
          * Following the last `.` signifies the unqualified basename of this
            attribute in this module (e.g., `Simulator`).
    '''

    # If this annotation is a new-style class, no further validation is needed.
    if is_class_new(annotation):
        pass
    # If this annotation is a tuple, raise an exception unless all members of
    # this tuple are new-style classes.
    elif isinstance(annotation, tuple):
        for member in annotation:
            if not is_class_new(member):
                raise TypeError(
                    '{} tuple member {} not a new-style class.'.format(
                        annotation_label, member))
    # If this function accepts string annotations as valid...
    elif is_str_valid:
        # If this annotation is a string, raise an exception unless this name
        # contains at least one "." delimiter and thus superficially appears to
        # be a fully-qualified module attribute name.
        #
        # Ideally, this function would also validate this module to be
        # importable and contain this attribute. Unfortunately, string
        # annotations are only leveraged to avoid circular import dependencies
        # (i.e., edge-cases in which two modules mutually import each other,
        # usually transitively rather than directly). Validating this module to
        # be importable and contain this attribute would necessitate importing
        # this module here. Since the @type_check decorator calling this
        # function is typically invoked via the global scope of a source
        # module, however, importing this target module here would be
        # functionally equivalent to importing that target module from that
        # source module -- triggering a circular import dependency in
        # susceptible source modules. Ergo, such validation *MUST* be deferred
        # to function call time.
        if isinstance(annotation, str):
            # If this name contains no "." delimiter, raise an exception.
            if '.' not in annotation:
                raise TypeError(
                    '{} {} not a fully-qualified module attribute name '
                    '(i.e., contains no "." delimiter).'.format(
                        annotation_label, annotation))
        # Else, this annotation is of unsupported type. Raise an exception.
        else:
            raise TypeError(
                '{} {} unsupported (i.e., neither a '
                'new-style class, '
                'tuple of such classes, nor '
                'fully-qualified module attribute name).'.format(
                    annotation_label, annotation))
    # Else, this annotation is of unsupported type. Raise an exception.
    else:
        raise TypeError(
            '{} {} unsupported (i.e., neither a '
            'new-style class nor '
            'tuple of such classes).'.format(annotation_label, annotation))

# ....................{ OBSOLETE                           }....................
#FIXME: *ALL OF THE FOLLOWING FUNCTIONALITY SHOULD EVENTUALLY BE REMOVED.* The
#@type_check decorator defined above provides a substantially superior
#solution to this decidedly... unsavoury approach.

def is_bool(obj: object) -> bool:
    '''
    ``True`` only if the passed object is **boolean** (i.e., either `True` or
    `False`).
    '''
    return isinstance(obj, bool)


def is_char(obj: object) -> bool:
    '''
    ``True`` only if the passed object is a **character** (i.e., string of length
    1).
    '''
    return is_str(obj) and len(obj) == 1


def is_nonnone(obj: object) -> bool:
    '''
    ``True`` only if the passed object is _not_ `None`.
    '''
    return obj is not None

# ....................{ TESTERS ~ callable                 }....................
def is_callable(obj: object) -> bool:
    '''
    ``True`` only if the passed object is **callable** (e.g., function, method,
    class defining the special `__call__()` method).
    '''

    return callable(obj)

# ....................{ TESTERS ~ class                    }....................
def is_class(obj: object) -> bool:
    '''
    ``True`` only if the passed object is a class.
    '''

    return isinstance(obj, ClassType)


#FIXME: Nonsense! Remove. There are only new-style classes in Python 3.x.

# Note that this function is an unavoidable duplicate of the
# betse.util.type.classes.is_new() function. Why? Because this function is
# required by the type_check decorator itself required by that function. Minor
# code duplication is a negligible price to pay for type check validation.
def is_class_new(obj: object) -> bool:
    '''
    ``True`` only if the passed class is a new- rather than old-style class.
    '''

    return is_class(obj)

# ....................{ TESTERS ~ collection               }....................
def is_mapping(obj: object) -> bool:
    '''
    ``True`` only if the passed object is a **mapping** (i.e., indexable by
    strings).

    Customary mappings include `dict` and `OrderedDict` instances.
    '''
    return isinstance(obj, Mapping)

# ....................{ TESTERS ~ collection               }....................
def is_container(obj: object) -> bool:
    '''
    ``True`` only if the passed object is a **container** (i.e., implements the
    `__contains__` special method returning ``True`` only if that container
    contains the passed element).

    Most collections of interest (e.g., `dict`, `list`, `set`) are containers.
    '''
    return isinstance(obj, Container)


def is_container_nonstr(obj: object) -> bool:
    '''
    ``True`` only if the passed object is a **non-string container** (i.e.,
    implements the `__contains__` special method *and* is not a string).
    '''
    return is_container(obj) and not is_str(obj)

# ....................{ TESTERS ~ collection : iterable    }....................
def is_iterable(obj: object) -> bool:
    '''
    ``True`` only if the passed object is an **iterable**.

    Iterables are objects capable of returning their members one at a time.
    Equivalently, iterables implement the abstract base class
    `collections.Iterable` and hence define the `__iter__()` method.
    '''
    return isinstance(obj, Iterable)


def is_iterable_nonstr(obj: object) -> bool:
    '''
    ``True`` only if the passed object is a **non-string iterable** (i.e.,
    implements the abstract base class `collections.Iterable` *and* is not a
    string).
    '''
    return is_iterable(obj) and not is_str(obj)

# ....................{ TESTERS ~ sequence                 }....................
def is_sequence(obj: object) -> bool:
    '''
    ``True`` only if the passed object is a **sequence**.

    Sequences are iterables supporting efficient element access via integer
    indices. Equivalently, sequences implement the abstract base class
    :class:`collections.abc.Sequence` and hence define the ``__getitem__()`` and
    ``__len__()`` methods (among numerous others).

    While all sequences are iterables, not all iterables are sequences.
    Generally speaking, sequences correspond to the proper subset of iterables
    whose elements are ordered. :class:`dict` and :class:`OrderedDict` are the
    canonical examples. :class:`dict` implements :class:`collections.Iterable`
    but *not* :class:`collections.abc.Sequence`, due to *not* supporting integer
    index-based lookup; :class:`OrderedDict` implements both, due to supporting
    such lookup.
    '''

    return isinstance(obj, SequenceTypes)


def is_sequence_nonstr(obj: object) -> bool:
    '''
    ``True`` only if the passed object is a **non-string sequence** (i.e.,
    implements the abstract :class:`collections.abc.Sequence` base class *and*
    is not a string).

    For generality, this functions returns ``True`` for both pure-Python
    non-string sequences *and* non-Pythonic Fortran-based `numpy` arrays and
    matrices (which fail to subclass the :class:`collections.abc.Sequence` API
    despite implementing all methods defined by that subclass).
    '''

    return is_sequence(obj) and not is_str(obj)


def is_sequence_nonstr_nonempty(obj: object) -> bool:
    '''
    ``True`` only if the passed object is a **nonempty non-string sequence**
    (i.e., implements the abstract :class:`collections.abc.Sequence` base class,
    is not a string, *and* contains at least one element).

    '''
    return is_sequence_nonstr(obj) and len(obj)

# ....................{ TESTERS ~ error                    }....................
def is_exception(obj: object) -> bool:
    '''
    ``True`` only if the passed object is an **exception**.
    '''

    return isinstance(obj, Exception)

# ....................{ TESTERS ~ betse : science          }....................
def is_cells(obj: object) -> bool:
    '''
    ``True`` only if the passed object is an instance of the BETSE-specific
    `Cells` class.
    '''

    # Avoid circular import dependencies.
    from betse.science.cells import Cells
    return isinstance(obj, Cells)


def is_parameters(obj: object) -> bool:
    '''
    ``True`` only if the passed object is an instance of the BETSE-specific
    `Parameters` class.
    '''

    # Avoid circular import dependencies.
    from betse.science.parameters import Parameters
    return isinstance(obj, Parameters)


def is_simulator(obj: object) -> bool:
    '''
    ``True`` only if the passed object is an instance of the BETSE-specific
    `Simulator` class.
    '''

    # Avoid circular import dependencies.
    from betse.science.sim import Simulator
    return isinstance(obj, Simulator)

def is_simrunner(obj: object) -> bool:
    '''
    ``True`` only if the passed object is an instance of the BETSE-specific
    `SimRunner` class.
    '''
    from betse.science.simrunner import SimRunner
    return isinstance(obj, SimRunner)

# ....................{ TESTERS ~ lib : matplotlib         }....................
def is_matplotlib_mappable(obj: object) -> bool:
    '''
    ``True`` only if the passed object is a Matplotlib mappable.
    '''

    # Avoid importing third-party packages at the top level.
    from matplotlib.cm import ScalarMappable
    return isinstance(obj, ScalarMappable)


def is_matplotlib_polycollection(obj: object) -> bool:
    '''
    ``True`` only if the passed object is a Matplotlib **polygon collection**
    (i.e., an instance of the `PolyCollection` class).
    '''

    # Avoid importing third-party packages at the top level.
    from matplotlib.collections import PolyCollection
    return isinstance(obj, PolyCollection)


def is_matplotlib_trimesh(obj: object) -> bool:
    '''
    ``True`` only if the passed object is a Matplotlib **triangle mesh** (i.e., an
    instance of the `matplotlib.collections.TriMesh` class).
    '''

    # Avoid importing third-party packages at the top level.
    from matplotlib.collections import TriMesh
    return isinstance(obj, TriMesh)

# ....................{ TESTERS ~ lib : py                 }....................
def is_py_path_local(obj: object) -> bool:
    '''
    ``True`` only if the passed object is a `py.path.local` instance.
    '''

    # Avoid importing third-party packages at the top level.
    from py._path.local import LocalPath
    return isinstance(obj, LocalPath)

# ....................{ TESTERS ~ lib : numpy              }....................
def is_pytest_fixture(obj: object) -> bool:
    '''
    ``True`` only if the passed object is a `py.test` fixture.
    '''

    # Avoid importing third-party packages at the top level.
    from _pytest.python import FixtureFunctionMarker
    return isinstance(obj, FixtureFunctionMarker)

# ....................{ TESTERS ~ numeric                  }....................
def is_int(obj: object) -> bool:
    '''
    ``True`` only if the passed object is an integer.
    '''

    return isinstance(obj, int)


def is_int_ge(obj: object, ge: int) -> bool:
    '''
    ``True`` only if the passed object is an integer greater than or equal to the
    second passed integer.
    '''

    assert is_int(ge), assert_not_int(ge)
    return is_int(obj) and obj >= ge


def is_int_gt(obj: object, gt: int) -> bool:
    '''
    ``True`` only if the passed object is an integer strictly greater than the
    second passed integer.
    '''

    assert is_int(gt), assert_not_int(gt)
    return is_int(obj) and obj > gt


def is_int_positive(obj: object) -> bool:
    '''
    ``True`` only if the passed object is a **positive integer** (i.e., is
    strictly greater than 0).
    '''

    return is_int_gt(obj, 0)


def is_numeric(obj: object) -> bool:
    '''
    ``True`` only if the passed object is **numeric** (i.e., instance of either
    the `int` or `float` types).
    '''

    return isinstance(obj, (int, float))

# ....................{ TESTERS ~ str                      }....................
def is_str(obj: object) -> bool:
    '''
    ``True`` only if the passed object is a **string** (i.e., instance of the
    `str` class).
    '''

    return isinstance(obj, str)


def is_str_or_none(obj: object) -> bool:
    '''
    ``True`` only if the passed object is either a string *or* `None`.
    '''

    return isinstance(obj, str) or obj is None


def is_str_nonempty(obj: object) -> bool:
    '''
    ``True`` only if the passed object is a **nonempty string* (i.e., string
    comprising one or more characters and hence _not_ the empty string).
    '''

    return is_str(obj) and len(obj)


def is_str_nonempty_or_none(obj: object) -> bool:
    '''
    ``True`` only if the passed object is either a nonempty string *or* `None`.
    '''

    return is_str_nonempty(obj) or obj is None

# ....................{ FORMATTERS                         }....................
# This string-centric function is defined in this module rather than in the
# arguably more appropriate "strs" module to drastically simplify assertion
# handlers defined below.
#
# To circumvent chicken-and-the-egg issues with the @type_check decorator, this
# function must be defined *AFTER* all testers required by this decorator.
@type_check
def trim(obj: object, max_len: int = 76) -> str:
    '''
    Convert the passed object into a terse human-readable string suitable for
    safe interpolation into end-user messages.

    Parameters
    ----------
    obj: object
        Object whose representation is to be trimmed, converted into a string
        via the canonical :func:`repr` builtin.
    max_len: optional[int]
        Maximum length of the string to be returned.  Defaults to the customary
        line length of 80 characters minus default output indentation of four
        characters.

    Returns
    ----------
    str
        Human-readable string trimmed from the string representation of the
        passed object to the passed maximum length.
    '''

    # Uncompiled regular expression grouping zero or more non-newline leading
    # characters preceding this maximum length *AND* zero or more trailing
    # delimiters.
    PRE_MAX_CHARS_LAST_DELIMITERS_REGEX = (
        r'^([^\n]{0,' + str(max_len) + r'}).*?([\])}>\'"]*)$')

    # String describing the passed object. For debuggability, the verbose
    # (albeit less human-readable) output of repr() is preferred to the terse
    # (albeit more human-readable) output of str().
    obj_synopsis = repr(obj)

    # If this synopsis either exceeds this maximum length *OR* contains a
    # newline, replace the substring of this synopsis from whichever of the
    # first character following this maximum length or the first newline occurs
    # first to the string end (excluding any # optional trailing delimiters)
    # with a single ellipses.
    if len(obj_synopsis) > max_len or '\n' in obj_synopsis:
        obj_synopsis = re.sub(
            PRE_MAX_CHARS_LAST_DELIMITERS_REGEX,
            r'\1...\2',
            obj_synopsis,
            flags=re.DOTALL
        )

    # Return this synopsis.
    return obj_synopsis

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

# ....................{ ASSERTERS ~ error                  }....................
def assert_not_exception(obj: object) -> str:
    '''
    String asserting the passed object to _not_ be an exception.
    '''
    return '"{}" not an exception.'.format(trim(obj))

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

# ....................{ ASSERTERS ~ lib : matplotlib       }....................
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


def assert_not_int_positive(obj: object) -> str:
    '''
    String asserting the passed object to _not_ be a positive integer.
    '''

    return '"{}" not a positive integer.'.format(trim(obj))


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
