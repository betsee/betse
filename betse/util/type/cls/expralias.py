#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level **expression alias** (i.e.,
objects satisfying the descriptor protocol
defined at class scope) facilities.
'''

#FIXME: If optimizations are enabled (i.e., "if not __debug__:"), avoid
#embedding type validation in the data descriptor classes defined below.

# ....................{ IMPORTS                            }....................
from betse.exceptions import (
    BetseDescriptorException, BetseEnumException, BetseTypeException)
from betse.util.type.types import (
    type_check,
    CallableOrNoneTypes,
    EnumType,
    StrOrNoneTypes,
    TestableOrNoneTypes,
)

if False: BetseTypeException   # squelch IDE warnings

# ....................{ GLOBALS                            }....................
_EXPR_ALIAS_ID = 0
'''
Unique arbitrary identifier with which to uniquify the class name of the next
:func:`expr_alias` descriptor.
'''

# ....................{ CLASSES                            }....................
class ExprAliasBound(object):
    '''
    Expression alias, dynamically referring to an arbitrarily complex source
    Python expression suitable for use as both the left- and right-hand sides of
    Python assignment statements of the passed type and/or satisfying the passed
    predicate.

    This alias requires the parent object in this expression to be passed to the
    :meth:`__init__` method and hence is effectively "bound" to that object at
    instantiation, implying that each instance of this alias may service
    exactly one object.

    Motivation
    ----------
    This alias is principally intended to circumvent the conflict between
    immutable types in Python (e.g., :class:`int`, :class:`str`) and Python's
    pass-by-reference semantics.

    Consider two similar scenarios:

    * Passing an integer variable to a function as an integer parameter. In this
      case, modifying the latter does *not* modify the former.
    * Passing an instance of this class encapsulating an integer variable to a
      function as a parameter expecting only an instance of this class. In this
      case, modifying this instance actually modifies the integer variable
      encapsulated by this instance -- thereby circumventing this conflict.

    Caveats
    ----------
    This class is intended to be instantiated *only* when the expression to be
    aliased is only dynamically known at runtime. When this expression is
    instead statically known at class definition time, the lower-level
    :func:`expr_alias` function should be called instead.

    Attributes
    ----------
    _expr_alias : object
        Expression alias data descriptor of dynamically defined type created and
        returned by the :func:`expr_alias` function.
    _obj : object
        Object to be bound to the variable in this expression whose name is the
        ``obj_name`` parameter passed to the :meth:`__init__` method of this
        class, defaulting to ``self``.

    See Also
    ----------
    :class:`ExprAliasUnbound`
        Expression alias *not* bound to its parent object at instantiation time.
    :func:`expr_alias`
        Lower-level function creating the data descriptor encapsulated by this
        higher-level class, intended to be used where declaring a descriptor at
        class scope is appropriate (e.g., where the expression to be aliased is
        known at class definition time).
    '''

    # ..................{ INITIALIZERS                       }..................
    def __init__(self, obj: object, *args, **kwargs) -> None:
        '''
        Initialize this expression alias.

        Parameters
        ----------
        obj : object
            Object to be bound to the variable in this expression whose name is
            the passed ``obj_name`` parameter, defaulting to ``self``.

        All remaining parameters are passed as is to the :func:`expr_alias`
        function.
        '''

        # Expression alias data descriptor passed all passed parameters. While
        # instantiating descriptors outside of class scope is feasible, doing so
        # circumvents Python's default object.__getattribute__() implementation
        # transforming every attribute access "b.x" where "x" is a data
        # descriptor declared at class scope  into
        # "type(b).__dict__['x'].__get__(b, type(b))". Circumventing this
        # machinery typically defeats the purpose of instantiating a descriptor.
        # but is warranted in this case.
        #
        # Why? Because this class neither requires nor desires this descriptor
        # to ever be accessed via special machinery. Moreover, this descriptor
        # is only instantiated by this class to avoid code duplication. Since
        # the __get__() and __set__() methods defined by this descriptor
        # already implement all functionality required by this class, calling
        # these methods suffices to implement this class. Nice, yeah?
        self._expr_alias = expr_alias(*args, **kwargs)

        # Classify all remaining parameters.
        self._obj = obj

    # ..................{ GETTERS                            }..................
    def get(self) -> object:
        '''
        Value this expression currently evaluates to.
        '''

        # Yes, the descriptor protocol explicitly supports direct method calls.
        # In particular, the "Invoking Descriptors" subsection of the official
        # "Descriptor HowTo Guide" notes:
        #
        #     A descriptor can be called directly by its method name. For
        #     example, d.__get__(obj).
        #
        # Since this data descriptor is accessed with respect to an instance
        # variable rather than a class, the second parameter passed to this
        # method call is "None".
        #
        # See: https://docs.python.org/3/howto/descriptor.html
        return self._expr_alias.__get__(self._obj, None)

    # ..................{ SETTERS                            }..................
    def set(self, value: object) -> None:
        '''
        Set this expression's current value to the passed value.
        '''

        self._expr_alias.__set__(self._obj, value)


class ExprAliasUnbound(object):
    '''
    Expression alias, dynamically referring to an arbitrarily complex source
    Python expression suitable for use as both the left- and right-hand sides of
    Python assignment statements of the passed type and/or satisfying the passed
    predicate.

    This alias requires the parent object in this expression to be passed to the
    :meth:`get` and :meth:`set` methods, implying that each instance of this
    alias may service multiple different objects.

    Attributes
    ----------
    _expr_alias : object
        Expression alias data descriptor of dynamically defined type created and
        returned by the :func:`expr_alias` function.

    See Also
    ----------
    :class:`ExprAliasBound`
        Expression alias bound to its parent object at instantiation time.
    '''

    # ..................{ INITIALIZERS                       }..................
    def __init__(self, *args, **kwargs) -> None:
        '''
        Initialize this expression alias.

        Parameters
        ----------
        All remaining parameters are passed as is to the :func:`expr_alias`
        function.
        '''

        # Expression alias data descriptor passed all passed parameters.
        self._expr_alias = expr_alias(*args, **kwargs)

    # ..................{ GETTERS                            }..................
    def get(self, obj: object) -> object:
        '''
        Value this expression currently evaluates to.

        Parameters
        ----------
        obj : optional[object]
            Object to be bound to the variable in this expression whose name was
            defined as the ``obj_name`` parameter passed to the :meth:`__init__`
            method of this class, defaulting to ``self``.

        Returns
        ----------
        object
            This expression's current value.
        '''

        return self._expr_alias.__get__(obj, None)

    # ..................{ SETTERS                            }..................
    def set(self, obj: object, value: object) -> None:
        '''
        Set this expression's current value to the passed value.

        Parameters
        ----------
        obj : object
            Object to be bound to the variable in this expression whose name was
            defined as the ``obj_name`` parameter passed to the :meth:`__init__`
            method of this class, defaulting to ``self``.
        value : object
            Value to set this expression to.
        '''

        self._expr_alias.__set__(obj, value)

# ....................{ DESCRIPTORS                        }....................
@type_check
def expr_alias(
    # Mandatory parameters.
    expr: str,

    # Optional parameters. For caller convenience, the second parameter passed
    # to this function is required to be the "cls" parameter; doing so permits
    # callers to concisely call this function with positional parameters only.
    cls: TestableOrNoneTypes = None,
    expr_settable: StrOrNoneTypes = None,
    predicate: CallableOrNoneTypes = None,
    predicate_expr: StrOrNoneTypes = None,
    predicate_label: StrOrNoneTypes = None,
    obj_name: StrOrNoneTypes = 'self',
) -> object:
    '''
    Expression alias **data descriptor** (i.e., object satisfying the data
    descriptor protocol, usually defined at class scope), dynamically aliasing a
    target variable of the passed type and/or satisfying the passed predicate
    bound to instances of the class instantiating this descriptor to an
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
      Callers preferring an alternate name for this local variable may pass the
      optional ``obj_name`` parameter. If this descriptor is retrieved:
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
    is highly discouraged. In turn, this implies this function is intended to be
    called *only* when the expression to be aliased is statically known at class
    definition time. When this is *not* the case (e.g., when this expression is
    only dynamically known at runtime), the higher-level :class:`ExprAlias`
    class should be instantiated instead.

    Additionally, the instance variable bound by the descriptor returned by this
    function is *not* deletable via the :func:`del` builtin. Attempting to do so
    reliably raises an :class:`AttributeError` exception. Why? To preserve
    backward compatibility with pre-Python 3.6 versions, which provide no means
    of doing so.

    Parameters
    ----------
    expr : str
        Arbitrarily complex Python expression suitable for use as at least the
        right-hand side of Python assignment statements, typically
        evaluating to the value of an arbitrarily nested dictionary key:
        * Prefixed by the name of the instance variable to which this dictionary
          is bound in instances of the class containing this descriptor (e.g.,
          ``self._config``).
        * Suffixed by one or more ``[``- and ``]``-delimited key lookups into
          the same dictionary (e.g.,
          ``['variable settings']['noise']['dynamic noise']``).
        If this expression is unsuitable for use as the left-hand side of Python
        assignment statements (e.g., due to expanding to a complex expression
        comprising multiple variables rather than to a single simple variable),
        the ``expr_settable`` parameter *must* be passed as well.
    expr_settable : optional[str]
        Arbitrarily complex Python expression suitable for use as the left-hand
        side of Python assignment statements, producing a data descriptor whose:
        * ``__get__()`` implementation embeds the value of the passed ``expr``
          parameter.
        * ``__set__()`` implementation embeds the value of this parameter.
        Defaults to ``None``, in which case this parameter defaults to the
        value of the ``expr`` parameter.
    cls : optional[TestableTypes]
        Either:
        * A class, in which case an exception is raised if the value of this
          expression is *not* an instance of this class.
        * A tuple of classes, in which case an exception is raised if the value
          of this expression is *not* an instance of at least one class in this
          tuple.
        Defaults to ``None``, in which case no such validation is performed.
    predicate : optional[CallableTypes]
        Callable passed the value of this expression as its first and only
        parameter and returning a boolean ``True`` only if this value satisfies
        arbitrary caller requirements and ``False`` otherwise. For uniformity,
        this callable should ideally return boolean values rather than raise
        type exceptions. If this callable returns ``False``, an exception is
        raised on behalf of this callable. If the ``predicate_label``
        parameter is *not* also passed, an exception is raised. Defaults to
        ``None``, in which case no such validation is performed.
    predicate_expr : optional[str]
        Arbitrarily complex Python expression suitable for use as the condition
        of an if statement evaluating to a boolean ``True`` only if the value of
        this expression satisfies arbitrary caller requirements and ``False``
        otherwise. If this expression evaluates to ``False``, an exception is
        raised on behalf of this expression. If the ``predicate_label``
        parameter is *not* also passed, an exception is raised. Defaults to
        ``None``, in which case no such validation is performed.
    predicate_label : optional[str]
        Human-readable adjective or adjectival phrase describing the passed
        callable predicate if any. Defaults to ``None``, in which case this
        string is synthesized from the name of this callable if needed.
    obj_name : optional[str]
        Name of the variable in this expression whose value is the current
        instance of the class instantiating this descriptor. If this name is
        already reserved for internal use by this descriptor, an exception is
        raised; reserved names include ``cls``, ``self_descriptor``, and
        ``value``. Defaults to ``self``.

    Returns
    ----------
    object
        Expression alias data descriptor as detailed above.

    See Also
    ----------
    :class:`ExprAlias`
        Higher-level class encapsulating the data descriptor returned by this
        lower-level function, intended to be used where declaring a descriptor
        at class scope is inappropriate (e.g., where the expression to be
        aliased is unknown at class definition time).
    '''

    # Avoid circular import dependencies.
    from betse.util.type.cls import classes

    # Set of the names of all parameters hard-coded into the implementations of
    # either the __get__() or __set__() methods defined below and hence reserved
    # for internal use by this descriptor.
    RESERVED_ARG_NAMES = {'cls', 'self_descriptor', 'value',}

    # If the passed object name is already reserved, raise an exception.
    if obj_name in RESERVED_ARG_NAMES:
        raise BetseDescriptorException(
            'Object name "{}" reserved for internal use.'.format(obj_name))

    # Dictionary eventually containing only the following three keys:
    #
    # * "__init__", whose value is the __init__() special method for this
    #   descriptor class defined below.
    # * "__get__", whose value is the __get__() special method for this
    #   descriptor class defined below.
    # * "__set__", whose value is the __set__() special method for this
    #   descriptor class defined below.
    class_method_name_to_func = {
        # Passed validational class if any.
        '__expr_alias_cls': cls,

        # Passed validational predicate if any.
        '__expr_alias_predicate': predicate,
    }

    # Expression to be embedded in this data descriptor's __get__() method
    expr_gettable = expr

    # Expression to be embedded in this data descriptor's __set__() method
    if expr_settable is None:
        expr_settable = expr_gettable

    # Python code snippet listing all optional arguments to be accepted by this
    # data descriptor's __init__() method.
    class_init_args = ''

    # Python code snippet implementing the body of this data descriptor's
    # __init__(), __get__(), and __set__() methods respectively.
    class_init_body = ''
    class_get_body = ''
    class_set_body = ''

    # Python code snippet validating the value this expression evaluates to.
    value_test_block = ''

    # If a validational class was passed...
    if cls is not None:
        # Pass this class to this method.
        class_init_args += ', __expr_alias_cls=__expr_alias_cls'

        # Classify this passed class in this method.
        class_init_body += '''
    self_descriptor.__expr_alias_cls = __expr_alias_cls'''

        # Raise an exception unless the value to which this expression evaluates
        # is of the expected type(s). While this type validation could also be
        # performed by decorating the __get__() and __set__() methods defined
        # below by the @type_check decorator, doing so would impose additional
        # overhead for little gain.
        value_test_block += '''
    if not isinstance(value, self_descriptor.__expr_alias_cls):
        raise BetseTypeException(
            'Expression alias value {{!r}} not a {{!r}}.'.format(
                value, self_descriptor.__expr_alias_cls))
    '''

    # If a validational predicate was passed...
    if predicate is not None:
        # If no predicate label was passed, raise an exception.
        if predicate_label is None:
            raise BetseDescriptorException(
                'Parameter "predicate_expr" but not "predicate_label" passed.')

        # Pass this predicate to this method.
        class_init_args += ', __expr_alias_predicate=__expr_alias_predicate'

        # Classify this passed predicate in this method.
        class_init_body += '''
    self_descriptor.__expr_alias_predicate = __expr_alias_predicate'''

        # Raise an exception unless the value to which this expression evaluates
        # satisfies the same predicate.
        value_test_block += '''
    if not self_descriptor.__expr_alias_predicate(value):
        raise BetseTypeException(
            'Expression alias value {{!r}} not {predicate_label}.'.format(value))
    '''.format(predicate_label=predicate_label)

    # If a validational predicate expression was passed...
    if predicate_expr is not None:
        # If no predicate label was passed, raise an exception.
        if predicate_label is None:
            raise BetseDescriptorException(
                'Parameter "predicate_expr" but not "predicate_label" passed.')

        # Raise an exception unless the value to which this expression evaluates
        # satisfies the same predicate.
        value_test_block += '''
    if not ({predicate_expr}):
        raise BetseTypeException(
            'Expression alias value {{!r}} not {predicate_label}.'.format(value))
    '''.format(predicate_expr=predicate_expr, predicate_label=predicate_label)

    # If the __init__() method's body remains empty, implement this body as a
    # noop to ensure syntactic validity.
    if not class_init_body:
        class_init_body = '''
    pass'''

    # If this expression is validated...
    if value_test_block:
        # Implement the __get__() method body to validate this expression.
        class_get_body = '''
    # Value to which this expression evaluates to be returned.
    value = {expr_gettable}

    # Validate this value.
    {value_test_block}

    # Return this value.
    return value'''.format(
        expr_gettable=expr_gettable, value_test_block=value_test_block)

        # Implement the __set__() method body to validate this expression.
        class_set_body = '''
    # Validate this value.
    {value_test_block}

    # Set this expression to this value.
    {expr_settable} = value'''.format(
        expr_settable=expr_settable, value_test_block=value_test_block)
    # Else, this expression is unvalidated. For efficiency, reduce the
    # __get__() and __set__() method bodies to the expected one-liners.
    else:
        class_get_body = '''
    return {expr_gettable}'''.format(expr_gettable=expr_gettable)
        class_set_body = '''
    {expr_settable} = value'''.format(expr_settable=expr_settable)

    #FIXME: Improve exception handling to report the name of the data descriptor
    #in question. To do so, we'll need to implement the data descriptor method
    #passing this data descriptor its name. Only available under Python 3.6,
    #we believe, but that's quite alright. (See below for method signature.)

    # Python code snippet declaring these special methods.
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
    class_body = '''
def __init__(self_descriptor{init_args}):
    {init_body}

def __get__(self_descriptor, {obj_name}, cls):
    {get_body}

def __set__(self_descriptor, {obj_name}, value):
    {set_body}
'''.format(
    init_args=class_init_args,
    init_body=class_init_body,
    get_body=class_get_body,
    set_body=class_set_body,
    obj_name=obj_name,
)

    # Define these methods and this dictionary containing these methods.
    exec(class_body, globals(), class_method_name_to_func)

    # Prevent input parameters passed into this snippet by this exec() call from
    # polluting the attribute namespace of this class.
    del class_method_name_to_func['__expr_alias_cls']
    del class_method_name_to_func['__expr_alias_predicate']

    # Descriptor class with this name containing only these methods.
    expr_alias_class = classes.define_class(
        class_name=_get_expr_alias_class_name(),
        class_attr_name_to_value=class_method_name_to_func,
    )

    # Instantiate and return the singleton descriptor for this class.
    return expr_alias_class()

# ....................{ DESCRIPTORS ~ enum                 }....................
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
    # duplicates that of the enums.is_member_name() function.
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
    # of the enums.is_member() function.
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
