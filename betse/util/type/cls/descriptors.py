#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level **descriptor** (i.e., objects satisfying the descriptor protocol
defined at class scope) facilities.
'''

# ....................{ IMPORTS                            }....................
from betse.exceptions import BetseTypeException
from betse.util.type.types import type_check, ClassType

if False: BetseTypeException   # squelch IDE warnings

# ....................{ GLOBALS                            }....................
_EXPR_ALIAS_ID = 0
'''
Unique arbitrary identifier with which to uniquify the class name of the next
:func:`expr_alias` descriptor.
'''

# ....................{ DATA DESCRIPTORS                   }....................
#FIXME: Unit test us up, please.
@type_check
def expr_alias(expr: str, cls: ClassType = None) -> object:
    '''
    Expression alias **data descriptor** (i.e., object satisfying the data
    descriptor protocol), dynamically aliasing a target variable bound to
    instances of the class instantiating this descriptor to an arbitrarily
    complex source Python expression suitable for use as both the left- and
    right-hand sides of Python assignment statements.

    This function (in order):

    #. Dynamically defines a new descriptor class specific to the passed string
       such that:
       * Getting the value of the target variable internally gets the current
         value of the source expression.
       * Setting the value of the target variable either:
         * If this value is of the expected type, internally sets the current
           value of the source expression to this value.
         * Else, raises a type exception.
       * Deleting the target variable deletes this variable from the instance to
         which this variable is bound much as for a typical instance variable.
         Deleting this variable does *not* delete the variable bound to the
         source expression.
    #. Creates an instance of this descriptor class.
    #. Returns this instance.

    The typical use case for this class is to bind an instance variable to a
    nested entry of the dictionary produced by loading a YAML-formatted file.
    Doing so implicitly propagates modifications to this variable back to this
    dictionary and thus back to this file on saving this dictionary in YAML
    format back to this file.

    Caveats
    ----------
    As with all descriptors, this function is intended to be called *only* at
    **class scope** (i.e., directly from within the body of a class) in a manner
    assigning the descriptor returned by this function to a class variable of
    arbitrary name. While feasible, calling this function from any other scope
    is highly discouraged.

    Parameters
    ----------
    expr : str
        Arbitrarily complex Python expression suitable for use as both the left-
        and right-hand sides of Python assignment statements, typically
        evaluating to the value of an arbitrarily nested dictionary entry:
        * Prefixed by the name of the instance variable to which this dictionary
          is bound in instances of the class containing this descriptor (e.g.,
          ``self._config``).
        * Suffixed by one or more `[`- and `]`-delimited key lookups into the
          same dictionary (e.g.,
          ``['variable settings']['noise']['dynamic noise']``).
    cls: optional[ClassType]
        Either:
        * A class, in which case attempting to set this variable to a value
          *not* an instance of this class raises an exception.
        * A callable, in which caseattempting to set this variable to a value
          passes this value to this callable first, which may then raise an
          exception of arbitrary type.
        Defaults to `None`, in which case this variable is permissively settable
        to any arbitrary value.

    Returns
    ----------
    object
        Expression alias data descriptor as detailed above.
    '''

    # Global variables set below.
    global _EXPR_ALIAS_ID

    # Avoid circular import dependencies.
    from betse.util.type.cls import classes

    # Name of this descriptor class. Since the human-readability of this name is
    # *NOT* required, a non-human-readable name is efficiently created from a
    # private prefix and a unique arbitrary identifier.
    expr_alias_class_name = '__ExprAlias' + str(_EXPR_ALIAS_ID)

    # Increment this identifier, preserving uniqueness.
    _EXPR_ALIAS_ID += 1

    # Dictionary eventually containing only the following three keys:
    #
    # * "__init__", whose value is the __init__() special method for this
    #   descriptor class defined below.
    # * "__get__", whose value is the __get__() special method for this
    #   descriptor class defined below.
    # * "__set__", whose value is the __set__() special method for this
    #   descriptor class defined below.
    expr_alias_class_methods = {
        # Type of all values accepted by this __set__() method if any.
        '__var_alias_cls': cls,
    }

    # Snippet of Python expr declaring these special methods.
    #
    # Note the differentiation between the "self_descriptor" parameter referring
    # to the current descriptor and the "self" parameter referring to the object
    # containing this descriptor, which may be referenced by the passed
    # expression and must thus be assigned the name "self".
    expr_alias_class_method_defs = None

    # If *NOT* validating the type of this variable, define this data descriptor
    # to permissively accept all possible values.
    if cls is None:
        expr_alias_class_method_defs = '''
def __get__(self_descriptor, self):
    return {expr}

def __set__(self_descriptor, self, value):
    {expr} = value
'''.format(expr=expr)
    # Else, define this data descriptor to validate the type of this variable.
    else:
        expr_alias_class_method_defs = '''
def __init__(self_descriptor, __var_alias_cls=__var_alias_cls):
    self_descriptor._var_alias_cls = __var_alias_cls

def __get__(self_descriptor, self):
    return {expr}

def __set__(self_descriptor, self, value):
    if not isinstance(value, self_descriptor._var_alias_cls):
        raise BetseTypeException(
            'Expression "{expr}" of type {{}} not settable to {{!r}}.'.format(
                self_descriptor._var_alias_cls.__name__, value))

    {expr} = value
'''.format(expr=expr)

    # Define these methods and this dictionary containing these methods.
    exec(expr_alias_class_method_defs, globals(), expr_alias_class_methods)

    # Prevent input parameters passed into this snippet by this exec() call from
    # polluting the attribute namespace of this class.
    del expr_alias_class_methods['__var_alias_cls']

    # Descriptor class with this name containing only these methods.
    expr_alias_class = classes.define_class(
        class_name=expr_alias_class_name,
        class_attr_name_to_value=expr_alias_class_methods,
    )

    # Instantiate and return the singleton descriptor for this class.
    return expr_alias_class()
