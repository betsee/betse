#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level object facilities.
'''

# ....................{ IMPORTS                           }....................
from betse.exceptions import BetseAttrException, BetseCallableException
from betse.util.type.types import (
    type_check,
    CallableTypes,
    CallableOrNoneTypes,
    ClassType,
    TestableOrNoneTypes,
)

# ....................{ GETTERS : attr                    }....................
@type_check
def get_attr(obj: object, attr_name: str, **kwargs) -> object:
    '''
    Value of the attribute with the passed name bound to the passed object if
    this object defines this attribute *or* raise an exception otherwise (i.e.,
    if this object defines no such attribute), optionally validated to be of
    the passed type.

    Parameters
    ----------
    obj : object
        Object to be inspected.
    attr_name : str
        Name of the attribute to be retrieved.

    All remaining keyword arguments are passed as is to the
    :func:`get_attr_or_default` function.

    Returns
    ----------
    object
        Value of the attribute with this name bound to this object.

    Raises
    ----------
    BetseAttrException
        If no such attribute is bound to this object.
    BetseTypeException
        If the ``attr_type`` parameter is non-``None`` and the type of the
        current value of this attribute is *not* an instance of ``attr_type``.
    '''

    # Avoid circular import dependencies.
    from betse.util.py.module import pymodule
    from betse.util.type.obj.sentinels import SENTINEL

    # Value of the attribute with this name defined by this object if any *OR*
    # the sentinel otherwise.
    attr_value = get_attr_or_sentinel(obj, attr_name, **kwargs)

    # If no such attribute exists...
    if attr_value is SENTINEL:
        # Human-readable name of this object. Specifically:
        #
        # * If this object is a module, the fully-qualified name of this
        #   module (e.g., "betse.metadeps").
        # * Else, the unqualified name of the class of this object.
        #
        # Since the unqualified name of the class of module objects is simply
        # "module", differentiating these two common cases improves sanity.
        obj_name = (
            pymodule.get_name_qualified(obj) if pymodule.is_module(obj) else
            get_class_name_unqualified(obj))

        # Raise an exception embedding this object name.
        raise BetseAttrException(
            'Attribute "{}.{}" undefined.'.format(obj_name, attr_name))
    # Else, this attribute exists.

    # Return this value.
    return attr_value


@type_check
def get_attr_or_none(obj: object, attr_name: str, **kwargs) -> object:
    '''
    Attribute with the passed name bound to the passed object if this object
    defines this attribute *or* ``None`` otherwise (i.e., if this object
    defines no such attribute), optionally validated to be of
    the passed type.

    Caveats
    ----------
    For disambiguity, consider calling the :func:`get_attr_or_sentinel`
    function instead. Whereas this function fails to distinguish between
    existing attributes whose values are ``None`` and non-existing attributes
    for which this function returns ``None``, the :func:`get_attr_or_sentinel`
    function trivially disambiguates between these two common edge cases.

    Parameters
    ----------
    obj : object
        Object to be inspected.
    attr_name : str
        Name of the attribute to be obtained.

    All remaining keyword arguments are passed as is to the
    :func:`get_attr_or_default` function.

    Returns
    ----------
    object
        Attribute with this name bound to this object if any *or* ``None``
        otherwise.

    Raises
    ----------
    BetseTypeException
        If the ``attr_type`` parameter is non-``None`` and the type of the
        current value of this attribute is *not* an instance of ``attr_type``.
    '''

    # Return the current value of the attribute with this name defined by this
    # object if any *OR* "None" otherwise.
    return get_attr_or_default(
        obj=obj, attr_name=attr_name, attr_default=None, **kwargs)


@type_check
def get_attr_or_sentinel(obj: object, attr_name: str, **kwargs) -> object:
    '''
    Value of the attribute with the passed name bound to the passed object if
    this object defines this attribute *or* the sentinel singleton otherwise
    (i.e., if this object defines no such attribute), optionally validated to
    be of the passed type.

    This function enables callers to safely distinguish between non-existing
    attributes and existing attributes whose values are ``None``.

    Parameters
    ----------
    obj : object
        Object to be inspected.
    attr_name : str
        Name of the attribute to be obtained.

    All remaining keyword arguments are passed as is to the
    :func:`get_attr_or_default` function.

    Returns
    ----------
    object
        Either:

        * If this object declares this attribute, this attribute's value.
        * Else, the **sentinel singleton** (i.e.,
          :attr:`betse.util.type.obj.sentinels.SENTINEL`).

    Raises
    ----------
    BetseTypeException
        If the ``attr_type`` parameter is non-``None`` and the type of the
        current value of this attribute is *not* an instance of ``attr_type``.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.obj.sentinels import SENTINEL

    # Return the current value of the attribute with this name defined by this
    # object if any *OR* the sentinel otherwise.
    return get_attr_or_default(
        obj=obj, attr_name=attr_name, attr_default=SENTINEL, **kwargs)


@type_check
def get_attr_or_default(
    # Mandatory parameters.
    obj: object,
    attr_name: str,
    attr_default: object,

    # Optional parameters.
    attr_type: TestableOrNoneTypes = None,
) -> object:
    '''
    Value of the attribute with the passed name bound to the passed object if
    this object defines this attribute *or* the passed default value otherwise
    (i.e., if this object defines no such attribute), optionally validated to
    be of the passed type.

    Parameters
    ----------
    obj : object
        Object to be inspected.
    attr_name : str
        Name of the attribute to return the current value of.
    attr_default : object
        Default value to be returned if this object defines no such attribute.
    attr_type : TestableOrNoneTypes
        Expected type of the current value of this attribute. This function
        effectively performs the equivalent of the :meth:`type_check` decorator
        at runtime by raising an exception if all of the following apply:

        * This type is *not* ``None``.
        * This value is *not* this default value, implying this attribute to be
          defined by this object.
        * This value is *not* an instance of this type.

        Defaults to ``None``, in which case no such type checking is performed.

    Returns
    ----------
    object
        Either:

        * If this object defines this attribute, this attribute's value.
        * Else, this default value.

    Raises
    ----------
    BetseTypeException
        If the ``attr_type`` parameter is non-``None`` and the type of the
        current value of this attribute is *not* an instance of ``attr_type``.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.obj import objtest

    # Value of the attribute with this name defined by this object if any *OR*
    # this default value otherwise.
    attr_value = getattr(obj, attr_name, attr_default)

    # If this value is to be type-checked *AND* is *NOT* this default value
    # (which by definition already satisfies caller requirements regardless of
    # type), type-check this value.
    if attr_type is not None and attr_value is not attr_default:
        objtest.die_unless_instance(obj=attr_value, cls=attr_type)

    # Return this value.
    return attr_value

# ....................{ GETTERS : callable                }....................
@type_check
def get_callable(obj: object, callable_name: str) -> CallableTypes:
    '''
    Callable with the passed name defined by the passed object if this object
    defines such a callable *or* raise an exception otherwise (i.e., if this
    object defines no such callable).

    Parameters
    ----------
    obj : object
        Object to be queried.
    callable_name : str
        Name of the callable to be returned.

    Returns
    ----------
    CallableTypes
        Callable with this name defined by this object.

    Raises
    ----------
    BetseCallableException
        If this object defines no callable with this name.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.obj.sentinels import SENTINEL

    # Attribute with this name defined by this object if any *OR* the sentinel.
    # To raise human-readable exceptions, this attribute is *NOT* retrieved via
    # the higher-level get_callable_or_none() method; doing so would obscure
    # whether this attribute does not exist or does but is not callable.
    func = get_attr_or_sentinel(obj, callable_name)

    # If this attribute does *NOT* exist, raise an exception.
    if func is SENTINEL:
        raise BetseCallableException(
            'Method {}.{}() undefined.'.format(
                obj.__class__.__name__, callable_name))

    # If this attribute exists but is *NOT* callable, raise an exception.
    if not callable(func):
        raise BetseCallableException(
            'Object attribute "{}.{}" not callable: {!r}'.format(
                obj.__class__.__name__, callable_name, func))
    # Else, this attribute is callable.

    # Return this callable.
    return func


@type_check
def get_callable_or_none(obj: object, callable_name: str) -> (
    CallableOrNoneTypes):
    '''
    Callable with the passed name defined by the passed object if this object
    defines such a callable *or* ``None`` otherwise (i.e., if this object
    defines no such callable).

    Parameters
    ----------
    obj : object
        Object to be queried.
    callable_name : str
        Name of the callable to be returned.

    Returns
    ----------
    CallableOrNoneTypes
        Either:

        * If this object defines this callable, this callable.
        * Else, ``None``.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.obj.sentinels import SENTINEL

    # Attribute with this name defined by this object if any *OR* the sentinel.
    func = get_attr_or_sentinel(obj, callable_name)

    # Return this attribute if this attribute is callable *OR* "None".
    return func if func is not SENTINEL and callable(func) else None

# ....................{ GETTERS ~ class                   }....................
def get_class(obj: object) -> ClassType:
    '''
    Passed object if this object is itself a class *or* the class of this
    object otherwise (i.e., if this object is *not* a class).

    Parameters
    ----------
    obj : object
        Object to be queried for its class.

    Returns
    ----------
    ClassType
        This object if this object is a class *or* this object's class.
    '''

    # Simplicity is not a place in Simple City.
    return obj if isinstance(obj, ClassType) else type(obj)

# ....................{ GETTERS ~ class : name            }....................
def get_class_name_unqualified(obj: object) -> str:
    '''
    Unqualified name of either the passed object if this object is itself a
    class *or* the class of this object otherwise (i.e., if this object is
    *not* a class).

    Parameters
    ----------
    obj : object
        Object to be queried for its class name.

    Returns
    ----------
    str
        Unqualified name of this class.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.cls import classes

    # This object if this object is a class *OR* this object's class otherwise.
    cls = get_class(obj)

    # Return the unqualified name of this class.
    return classes.get_name_unqualified(cls)


def get_class_module_name_qualified(obj: object) -> str:
    '''
    Fully-qualified name of the module defining either the passed object if
    this object is itself a class *or* the class of this object otherwise
    (i.e., if this object is *not* a class).

    Parameters
    ----------
    obj : object
        Object to be queried for its module name.

    Returns
    ----------
    str
        Fully-qualified name of this module.

    Raises
    ----------
    BetseTypeException
        If this class has no ``__module__`` attribute, which should ideally
        *never* happen.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.cls import classes

    # This object if this object is a class *OR* this object's class otherwise.
    cls = get_class(obj)

    # Return the fully-qualified name of this class.
    return classes.get_module_name_qualified(cls)
