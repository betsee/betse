#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level **descriptor** (i.e., objects satisfying the descriptor protocol
defined at class scope) facilities.
'''

#FIXME: If optimizations are enabled (i.e., "if not __debug__:"), avoid
#embedding type validation in the data descriptor classes defined below.

# ....................{ IMPORTS                            }....................
from betse.exceptions import BetseEnumException, BetseTypeException
from betse.util.type.types import type_check, EnumType, TestableTypes

if False: BetseTypeException   # squelch IDE warnings

# ....................{ GLOBALS                            }....................
_EXPR_ALIAS_ID = 0
'''
Unique arbitrary identifier with which to uniquify the class name of the next
:func:`expr_alias` descriptor.
'''

# ....................{ DESCRIPTORS ~ alias                }....................
@type_check
def expr_alias(expr: str, cls: TestableTypes = None) -> object:
    '''
    Expression alias **data descriptor** (i.e., object satisfying the data
    descriptor protocol), dynamically aliasing a target variable of the passed
    type bound to instances of the class instantiating this descriptor to an
    arbitrarily complex source Python expression suitable for use as both the
    left- and right-hand sides of Python assignment statements.

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

    Expressions
    ----------
    The passed expression may be any arbitrarily complex source Python
    expression suitable for use as both the left- and right-hand sides of Python
    assignment statements. This expression may refer to the following local
    variables:

    * ``self``, the current instance of the class instantiating this descriptor.
      If this descriptor is retrieved:
      * As an **instance variable** (i.e., from the current instance of this
        class), this local is that instance.
      * As a **class variable** (i.e., from this class), this local is `None`.
        In this case, the ``cls`` local is guaranteed to be defined to the class
        instantiating this descriptor.
    * ``cls``, the class instantiating this descriptor. This local is only
      defined when getting this descriptor as a class rather than instance
      variable; equivalently, this local is undefined both when getting this
      descriptor as an instance variable *and* when setting this descriptor.
      Since one or the other of this and the ``self`` local is *always* defined,
      the class instantiating this descriptor is safely retrievable in
      expressions by referencing ``cls`` only if ``self`` is ``None``. For
      example, to get and set the name of this class:

      .. code:: python

         expr_alias(
             expr='(cls if self is None else self.__class__).__name__', cls=str)

    * ``self_descriptor``, the current instance of this descriptor. Since this
      descriptor only implements special methods required by the data descriptor
      protocol (e.g., ``__get__``, ``__set__``) and hence contains no data, this
      local is unlikely to be of interest to most callers.

    Caveats
    ----------
    As with all descriptors, this function is intended to be called *only* at
    **class scope** (i.e., directly from within the body of a class) in a manner
    assigning the descriptor returned by this function to a class variable of
    arbitrary name. While feasible, calling this function from any other scope
    is highly discouraged.

    Additionally, the instance variable bound by the descriptor returned by this
    function is *not* deletable via the :func:`del` builtin. Attempting to do so
    reliably raises an :class:`AttributeError` exception. Why? To preserve
    backward compatibility with pre-Python 3.6 versions, which provide no means
    of doing so.

    Parameters
    ----------
    expr : str
        Arbitrarily complex Python expression suitable for use as both the left-
        and right-hand sides of Python assignment statements, typically
        evaluating to the value of an arbitrarily nested dictionary key:
        * Prefixed by the name of the instance variable to which this dictionary
          is bound in instances of the class containing this descriptor (e.g.,
          ``self._config``).
        * Suffixed by one or more ``[``- and ``]``-delimited key lookups into
          the same dictionary (e.g.,
          ``['variable settings']['noise']['dynamic noise']``).
    cls: optional[TestableTypes]
        Either:
        * A class, in which case setting this variable to a value *not* an
          instance of this class raises an exception.
        * A tuple of classes, in which case setting this variable to a value
          *not* an instance of at least one class in this tuple raises an
          exception.
        * A callable, in which case setting this variable to a value passes this
          value to this callable first, which may then raise an exception of
          arbitrary type.
        Defaults to ``None``, in which case this variable is permissively
        settable to any arbitrary value.

    Returns
    ----------
    object
        Expression alias data descriptor as detailed above.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.cls import classes

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
        '__expr_alias_cls': cls,
    }

    # Snippet of Python expr declaring these special methods.
    expr_alias_class_method_defs = None

    # If *NOT* validating the type of this variable, define this data descriptor
    # to permissively accept all possible values.
    #
    # Note the differentiation between the "self_descriptor" parameter referring
    # to the current descriptor and the "self" parameter referring to the object
    # containing this descriptor, which may be referenced by the passed
    # expression and must thus be assigned the name "self".
    #
    # Note that the following special methods supported by the descriptor
    # protocol are intentionally left undefined:
    #
    # * "__set_name__(self, owner, name)", first introduced by Python 3.6 to
    #   notify descriptors of the variable names to which they have been bound.
    #   Since these names are irrelevant for aliasing purposes, this descriptor
    #   class leaves this special method unimplemented.
    # * "__delete__(self, instance)", deleting the variable to which this
    #   descriptor has been bound from the current instance of the class
    #   instantiating this descriptor. Sanely defining this method would require
    #   implementing the aforementioned __set_name__() method to store the name
    #   of the variable to which this descriptor is bound. That method is only
    #   available under Python 3.6, whereas this application is currently
    #   compatible with older Python 3.x versions. Since deletion of this
    #   descriptor is non-essential (and arguably undesirable), this descriptor
    #   class leaves this special method unimplemented.
    if cls is None:
        expr_alias_class_method_defs = '''
def __get__(self_descriptor, self, cls):
    return {expr}

def __set__(self_descriptor, self, value):
    {expr} = value
'''.format(expr=expr)
    # Else...
    else:
        # Snippet raising an exception unless the value to which this expression
        # evaluates is of the expected type(s). While this type validation could
        # also be performed by decorating the __get__() and __set__() methods
        # defined below by the @type_check decorator, doing so would impose
        # additional overhead for little gain.
        die_unless_value_type = '''
    if not isinstance(value, self_descriptor.__expr_alias_cls):
        raise BetseTypeException(
            'Expression alias value {{!r}} not a {{!r}}.'.format(
                value, self_descriptor.__expr_alias_cls))
'''

        # Define this data descriptor to validate the type of this variable.
        expr_alias_class_method_defs = '''
def __init__(self_descriptor, __expr_alias_cls=__expr_alias_cls):
    self_descriptor.__expr_alias_cls = __expr_alias_cls


def __get__(self_descriptor, self, cls):
    # Value to which this expression evaluates to be returned.
    value = {expr}

    # Raise an exception unless this value is of the expected type(s).
    {die_unless_value_type}

    # Return this value.
    return value


def __set__(self_descriptor, self, value):
    # Raise an exception unless this value is of the expected type(s).
    {die_unless_value_type}

    {expr} = value
'''.format(expr=expr, die_unless_value_type=die_unless_value_type)
    # print('expr_alias code:\n' + expr_alias_class_method_defs)

    # Define these methods and this dictionary containing these methods.
    exec(expr_alias_class_method_defs, globals(), expr_alias_class_methods)

    # Prevent input parameters passed into this snippet by this exec() call from
    # polluting the attribute namespace of this class.
    del expr_alias_class_methods['__expr_alias_cls']

    # Descriptor class with this name containing only these methods.
    expr_alias_class = classes.define_class(
        class_name=_get_expr_alias_class_name(),
        class_attr_name_to_value=expr_alias_class_methods,
    )

    # Instantiate and return the singleton descriptor for this class.
    return expr_alias_class()


@type_check
def expr_enum_alias(expr: str, enum_type: EnumType) -> object:
    '''
    Enumeration-specific expression alias **data descriptor** (i.e., object
    satisfying the data descriptor protocol), dynamically aliasing a target
    variable of the passed enumeration type bound to instances of the class
    instantiating this descriptor to an arbitrarily complex source Python
    expression suitable for use as both the left- and right-hand sides of Python
    assignment statements.

    Invariants
    ----------
    This data descriptor implicitly converts the value to which this expression
    evaluates to and from the corresponding member of this enumeration type in a
    case-insensitive manner. To facilitate this conversion, the caller *must*
    guarantee the following two invariants:

    #. This expression *must* evaluate to a string value. Where not the case, an
       exception will be explicitly raised.
    #. The names of all members of this enumeration type *must* be strictly
       uppercase. Where not the case (i.e., where members have either lowercase
       or mixed case names), an exception will be explicitly raised.

    Usage
    ----------
    If these invariants hold, then:

    * Getting the value of this data descriptor implicitly converts the string
      to which this expression evaluates into the member of this enumeration
      whose name is that string uppercased.
    * Setting the value of this data descriptor to a member of this enumeration
      implicitly sets the string to which this expression evaluates to the
      lowercased name of this member.

    Parameters
    ----------
    expr : str
        Arbitrarily complex Python expression suitable for use as both the left-
        and right-hand sides of Python assignment statements, typically
        evaluating to the value of an arbitrarily nested dictionary key:
        * Prefixed by the name of the instance variable to which this dictionary
          is bound in instances of the class containing this descriptor (e.g.,
          ``self._config``).
        * Suffixed by one or more ``[``- and ``]``-delimited key lookups into
          the same dictionary (e.g.,
          ``['variable settings']['noise']['dynamic noise']``).
    enum_type : EnumType
        Enumeration that the value of this variable *must* be a member of.
        Setting this variable to a value *not* a member of this enumeration will
        raise an exception.

    Returns
    ----------
    object
        Enemuration-specific expression alias data descriptor as detailed above.

    See Also
    ----------
    :func:`expr_alias`
        Further details.
    '''

    # Avoid circular import dependencies.
    from betse.util.type import strs
    from betse.util.type.cls import classes

    # Validate the names of all members of this enumeration to be uppercase.
    for enum_member in enum_type:
        if not strs.is_uppercase(enum_member.name):
            raise BetseEnumException(
                'Enumeration member name "{}" not uppercase.'.format(
                    enum_member.name))

    # Dictionary eventually containing only the following three keys:
    # "__init__", "__get__", and "__set__".
    expr_alias_class_methods = {
        # Type of all values accepted by this __set__() method.
        '__expr_enum_alias_type': enum_type,
    }

    # Snippet of Python expr declaring these special methods. See the
    # expr_alias() function for further details.
    expr_alias_class_method_defs = '''
def __init__(self_descriptor, __expr_enum_alias_type=__expr_enum_alias_type):
    self_descriptor.__expr_enum_alias_type = __expr_enum_alias_type


# Convert this variable's low-level string into a high-level enumeration member.
def __get__(self_descriptor, self, cls):
    # Name of the corresponding enumeration member to be returned.
    enum_member_name = {expr}

    # If this expression failed to evaluate to a string, raise an exception.
    if not isinstance(enum_member_name, str):
        raise BetseTypeException(
            'Expression alias value {{!r}} not a string.'.format(
                enum_member_name))

    # To ignore case, uppercase this name.
    enum_member_name = enum_member_name.upper()

    # If this uppercase name is *NOT* that of a member of this enumeration,
    # raise an exception. For both efficiency and readability, this logic
    # duplicates that of the enums.is_enum_member_name() function.
    if enum_member_name not in (
        self_descriptor.__expr_enum_alias_type.__members__):
        raise BetseEnumException(
            'Expression alias value {{!r}} unrecognized (i.e., not the name of '
            'an enumeration member of {{!r}}).'.format(
                enum_member_name.lower(),
                self_descriptor.__expr_enum_alias_type))

    # Get the enumeration member with this name.
    return self_descriptor.__expr_enum_alias_type[enum_member_name]


# Convert the passed high-level enumeration member into a low-level string.
def __set__(self_descriptor, self, enum_member):
    # If the passed value is *NOT* a member of this enumeration, raise an
    # exception. For both efficiency and readability, this logic duplicates that
    # of the enums.is_enum_member() function.
    if enum_member not in self_descriptor.__expr_enum_alias_type:
        raise BetseEnumException(
            'Expression alias value {{!r}} unrecognized (i.e., not an '
            'enumeration member of {{!r}}).'.format(
                enum_member, self_descriptor.__expr_enum_alias_type))

    # Set this expression to the lowercased name of this enumeration member.
    {expr} = enum_member.name.lower()
'''.format(expr=expr)
    # print('expr_alias code:\n' + expr_alias_class_method_defs)

    # Define these methods and this dictionary containing these methods.
    exec(expr_alias_class_method_defs, globals(), expr_alias_class_methods)

    # Prevent input parameters passed into this snippet by this exec() call from
    # polluting the attribute namespace of this class.
    del expr_alias_class_methods['__expr_enum_alias_type']

    # Descriptor class with this name containing only these methods.
    expr_alias_class = classes.define_class(
        class_name=_get_expr_alias_class_name(),
        class_attr_name_to_value=expr_alias_class_methods,
    )

    # Instantiate and return the singleton descriptor for this class.
    return expr_alias_class()

# ....................{ PRIVATE                            }....................
def _get_expr_alias_class_name() -> str:
    '''
    Name of the class of the next data descriptor created and returned by the
    next call to an expression alias function (e.g., :func:`expr_alias`),
    guaranteed to be unique across all such classes.

    Since the human-readability of this name is neither required nor desired,
    this is a non-human-readable name efficiently constructed from a private
    prefix and a unique arbitrary identifier.
    '''

    # Global variable set below.
    global _EXPR_ALIAS_ID

    # Increment this identifier, preserving uniqueness.
    _EXPR_ALIAS_ID += 1

    # Return a unique name suffixed by this identifier.
    return '__ExprAlias' + str(_EXPR_ALIAS_ID)
