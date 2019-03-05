#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level **enumeration** (i.e., enumerable types created by the `Enum()` class)
facilities.
'''

# ....................{ IMPORTS                           }....................
from betse.exceptions import BetseEnumException
from betse.util.type.types import (
    type_check,
    CallableTypes,
    EnumMemberType,
    EnumType,
    EnumClassType,
    GeneratorType,
    SequenceTypes,
    StrOrNoneTypes,
)
from functools import partial

# ....................{ SUBCLASSES                        }....................
class EnumOrdered(EnumClassType):
    '''
    Enumeration whose members are comparable according to their assigned
    values.

    This :class:`Enum` subclass complies with the `Functional API`_ for
    enumerations, permitting comparable enumeration types to be defined with a
    single function call. (See the example below.)

    Oddly, the :class:`Enum` superclass does *not* support such comparisons.
    This :class:`Enum` subclass amends this oversight, implementing four of the
    six rich comparison special methods. The fifth and sixth (i.e.,
    ``__eq__()`` and ``__ne__()``) are already implemented by the :class:`Enum`
    superclass and hence need *not* be reimplemented here.

    .. _Functional API:
       https://docs.python.org/3/library/enum.html#functional-api

    See Also
    ----------
    https://docs.python.org/3/library/enum.html#orderedenum
        Standard Python documentation strongly inspiring this implementation.

    Examples
    ----------
        >>> from betse.util.type.enums import EnumOrdered
        >>> Shrike = EnumOrdered('Shrike', ('hyperion', 'endymion'))
        >>> Shrike.endymion == Shrike.endymion
        True
        >>> Shrike.hyperion != Shrike.endymion
        True
        >>> Shrike.hyperion < Shrike.endymion
        True
        >>> Shrike.hyperion <= Shrike.endymion
        True
        >>> Shrike.endymion > Shrike.hyperion
        True
        >>> Shrike.endymion >= Shrike.hyperion
        True
    '''

    # ..................{ COMPARATORS                       }..................
    # For efficiency, rich comparison special methods are expected to return
    # the "NotImplemented" constant rather than raise the "NotImplementedError"
    # exception when the current and passed objects are incomparible.

    def __ge__(self, other: object) -> bool:
        '''
        :data:`NotImplemented` if the passed object is *not* a member of the
        same enumeration as this member, ``True`` if the value of this member
        is greater than or equal to that of the passed member, or ``False``
        otherwise.
        '''

        return (
            self.value >= other.value if self.__class__ is other.__class__ else
            NotImplemented)


    def __gt__(self, other: object) -> bool:
        '''
        :data:`NotImplemented` if the passed object is *not* a member of the
        same enumeration as this member, ``True`` if the value of this member
        is greater than that of the passed member, or ``False`` otherwise.
        '''

        return (
            self.value > other.value if self.__class__ is other.__class__ else
            NotImplemented)


    def __le__(self, other: object) -> bool:
        '''
        :data:`NotImplemented` if the passed object is *not* a member of the
        same enumeration as this member, ``True`` if the value of this member
        is less than or equal to that of the passed member, or ``False``
        otherwise.
        '''

        return (
            self.value <= other.value if self.__class__ is other.__class__ else
            NotImplemented)


    def __lt__(self, other: object) -> bool:
        '''
        :data:`NotImplemented` if the passed object is *not* a member of the
        same enumeration as this member, ``True`` if the value of this member
        is less than that of the passed member, or ``False`` otherwise.
        '''

        return (
            self.value < other.value if self.__class__ is other.__class__ else
            NotImplemented)

# ....................{ MAKERS                            }....................
@type_check
def make_enum(
    # Mandatory parameters.
    class_name: str,
    member_names: SequenceTypes,

    # Optional parameters.
    is_ordered: bool = False,
    doc: StrOrNoneTypes = None,
) -> EnumType:
    '''
    **Integer-based enumeration type** (i.e., instance of the standard
    :class:`EnumClassType` type whose members are uniquely mapped to 1-based
    integers) dynamically synthesized to have the passed class name and contain
    only the members with the passed names.

    This factory function is a convenience wrapper for the
    :class:`EnumClassType` subclass, whose non-standard semantics are arguably
    more obfuscatory than helpful.

    Attributes
    ----------
    Each attribute of the returned enumeration type is an enumeration member
    whose:

    * Type is the same as that of this enumeration.
    * Value is an object uniquely identifying this member in this enumeration,
      defining the following public attributes:

      * ``name``, the Python identifier uniquely identifying this member --
        guaranteed to be the same string as that passed to this function.
      * ``value``, the 1-based integer uniquely identifying this member. The
        ``value`` of this enumeration's:

        * First member is guaranteed to be 1.
        * Last member is guaranteed to be the number of members (i.e., the
          length of this enumeration).

    Parameters
    ----------
    class_name : str
        Class name of this enumeration type, ideally unique across all
        attributes of the module calling this function.
    member_names : SequenceTypes
        Sequence of the names of all members of this enumeration type, required
        to be valid **Python identifiers** (i.e., contain only alphanumeric
        characters and underscores).
    is_ordered : bool
        ``True`` only if the members of this enumeration are ordered (and hence
        comparable) according to the ordering of member names in the passed
        ``member_names`` sequence. Defaults to ``False``, in which case the
        members of this enumeration are unordered (and hence incomparable).
    doc : StrOrNoneTypes
        Class docstring to document this type with. Defaults to ``None``, in
        which case this type remains undocumented.

    Returns
    ----------
    EnumType
        Enumeration type dynamically synthesized as defined as above.

    Examples
    ----------
        >>> from betse.util.type.enums import make_enum
        >>> CthulhicState = make_enum(
        ...     class_name='CthulhicState',
        ...     member_names=('MERCIFUL', 'INNABILITY', 'CORRELATE', 'CONTENTS',))
        >>> CthulhicState.CORRELATE.name
        'CORRELATE'
        >>> CthulhicState.CORRELATE.value
        3
    '''

    # Avoid circular import dependencies.
    from betse.util.py import pyident
    from betse.util.type.call import callers
    from betse.util.type.iterable import itertest

    # If the passed class name is syntactically invalid, raise an exception.
    pyident.die_unless_unqualified(class_name)

    # If any passed member name is *NOT* a string, raise an exception.
    itertest.die_unless_items_instance_of(iterable=member_names, cls=str)

    # Fully-qualified name of the module defining this enumeration type.
    module_name = callers.get_caller_module_name()

    # Type of enumeration superclass to be subclassed. Specifically, if the
    # caller requested that the members of this enumeration be:
    #
    # * Ordered, this is our custom ordered enumeration type.
    # * Unordered, this is the standard unordered enumeration type.
    enum_superclass = EnumOrdered if is_ordered else EnumClassType

    # Dynamically synthesize and return this enumeration type.
    enum_subclass = enum_superclass(
        value=class_name,
        names=member_names,
        module=module_name,
        qualname='{}.{}'.format(module_name, class_name),
    )

    # If passed a docstring, assign this subclass this docstring.
    if doc is not None:
        enum_subclass.__doc__ = doc

    # Return this subclass.
    return enum_subclass

# ....................{ EXCEPTIONS                        }....................
def die_unless_member(
    enum_type: EnumType, enum_member: EnumMemberType) -> None:
    '''
    Raise an exception unless the passed enumeration contains the passed
    member.

    Parameters
    ----------
    enum_type : EnumType
        Enumeration type to be inspected.
    enum_member: str
        Enumeration members to test for.
    '''

    if not is_member(enum_type, enum_member):
        raise BetseEnumException(
            'Enumeration {} member "{}" not found.'.format(
                enum_type.__name__, enum_member.name))


def die_unless_member_name(
    enum_type: EnumType, enum_member_name: str) -> None:
    '''
    Raise an exception unless the passed enumeration contains an enumeration
    member with the passed name.

    Parameters
    ----------
    enum_type : EnumType
        Enumeration type to be inspected.
    enum_member_name: str
        Name of the member to test for.
    '''

    if not is_member_name(enum_type, enum_member_name):
        raise BetseEnumException(
            'Enumeration {} member "{}" not found.'.format(
                enum_type.__name__, enum_member_name))

# ....................{ TESTERS                           }....................
@type_check
def is_member(enum_type: EnumType, enum_member: EnumMemberType) -> bool:
    '''
    ``True`` only if the passed enumeration contains the passed member.

    Parameters
    ----------
    enum_type : EnumType
        Enumeration type to be inspected.
    enum_member_name: str
        Name of the member to test for.

    Returns
    ----------
    bool
        ``True`` only if this enumeration contains this member.
    '''

    # No, thank *YOU*, Enum.__contains__().
    return enum_member in enum_type


@type_check
def is_member_name(enum_type: EnumType, enum_member_name: str) -> bool:
    '''
    ``True`` only if the passed enumeration contains an enumeration member with
    the passed name.

    Parameters
    ----------
    enum_type : EnumType
        Enumeration type to be inspected.
    enum_member_name: str
        Name of the member to test for.

    Returns
    ----------
    bool
        ``True`` only if this enumeration contains a member with this name.

    See Also
    ----------
    https://stackoverflow.com/a/29795561/2809027
        StackOverflow answer strongly inspiring this implementation.
    '''

    # This is insanity, but insanity that works. (Yes, this is the accepted
    # method for testing enumeration member existence by name. Thanks, Guido.)
    return enum_member_name in enum_type.__members__

# ....................{ GETTERS ~ enum                    }....................
@type_check
def get_converter_name_to_uppercase_enum_member(
    enum_type: EnumType) -> CallableTypes:
    '''
    Callable accepting one string parameter returning the member of the passed
    enumeration whose name is that string uppercased, raising an exception if
    this string uppercased is *not* the name of a member of this enumeration.

    The callable returned by this function is principally intended to be passed
    as the ``type`` parameter to the :meth:`ArgumentParser.add_argument`
    method, converting from lowercase command-line option string arguments to
    corresponding enumeration members.

    Parameters
    ----------
    enum_type : EnumType
        Enumeration type to be inspected.

    Returns
    ----------
    CallableTypes
        Callable whose signature is defined as above.
    '''

    # Return a closure passing the passed enumeration to an existing getter.
    return partial(get_member_from_name_uppercased, enum_type=enum_type)


@type_check
def get_member_names_lowercase(enum_type: EnumType) -> SequenceTypes:
    '''
    Sequence of the unqualified lowercased names of all members of the passed
    enumeration type in **declaration order** (i.e., the order in which these
    members were originally declared).

    Parameters
    ----------
    enum_type : EnumType
        Enumeration type to be inspected.

    Returns
    ----------
    SequenceTypes
        Sequence of the lowercased names of all members of this enumeration.
    '''

    return tuple(iter_names_lowercase(enum_type))

# ....................{ GETTERS ~ member                  }....................
@type_check
def get_enum_name_from_member(enum_member: EnumMemberType) -> str:
    '''
    Unqualified name of the enumeration class to which the passed enumeration
    member belongs (e.g., ``SolverType`` given a member ``SolverType.FULL``).

    Parameters
    ----------
    enum_member : EnumMemberType
        Enumeration member to be inspected.

    Returns
    ----------
    str
        Unqualified name of this member's enumeration class.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.obj import objects

    # Strange, but true. We don't ask questions. Neither should you.
    return objects.get_class_name_unqualified(enum_member)


@type_check
def get_class_member_name(enum_member: EnumMemberType) -> str:
    '''
    Partially-qualified name of the passed enumeration member (e.g.,
    ``SolverType.FULL``), excluding the module defining this enumeration.

    Parameters
    ----------
    enum_member : EnumMemberType
        Enumeration member to be inspected.

    Returns
    ----------
    str
        Partially-qualified name of this member.

    See Also
    ----------
    :func:`get_member_name`
        Unqualified name of this member.
    '''

    # Makes sense.
    return str(enum_member)


@type_check
def get_member_name(enum_member: EnumMemberType) -> str:
    '''
    Unqualified name of the passed enumeration member (e.g., ``FULL`` given a
    member ``SolverType.FULL``).

    Since the :attr:`enum_member.name` attribute already provides this name,
    this function principally exists for orthogonality with other getters.

    Parameters
    ----------
    enum_member : EnumMemberType
        Enumeration member to be inspected.

    Returns
    ----------
    str
        Unqualified name of this member.

    See Also
    ----------
    :func:`get_class_member_name`
        Partially-qualified name of this member.
    '''

    # Believe it.
    return enum_member.name


@type_check
def get_member_name_lowercase(enum_member: EnumMemberType) -> str:
    '''
    Unqualified lowercase name of the passed enumeration member (e.g., ``full``
    given a member ``SolverType.FULL``).

    Parameters
    ----------
    enum_member : EnumMemberType
        Enumeration member to be inspected.

    Returns
    ----------
    str
        Unqualified lowercase name of this member.
    '''

    # It cannot be, yet it is.
    return enum_member.name.lower()

# ....................{ GETTERS ~ member : name           }....................
@type_check
def get_member_from_name(
    enum_type: EnumType, enum_member_name: str) -> EnumMemberType:
    '''
    Member of this enumeration with the passed name if this member exists *or*
    raise an exception otherwise (i.e., if this member does *not* exist).

    Parameters
    ----------
    enum_type: EnumType
        Enumeration type to inspect.
    enum_member_name : EnumType
        Name of the member of this enumeration to be returned.

    Returns
    ----------
    EnumMemberType
        Member of this enumeration with this name.

    Raises
    ----------
    BetseEnumException
        If no member of this enumeration with this name exists.
    '''

    # If no member with this name exists, raise an exception.
    die_unless_member_name(enum_type, enum_member_name)

    # Return this member.
    return enum_type[enum_member_name]


@type_check
def get_member_from_name_uppercased(
    enum_type: EnumType, enum_member_name: str) -> EnumMemberType:
    '''
    Member of this enumeration with the passed name uppercased if this member
    exists *or* raise an exception otherwise (i.e., if this member does *not*
    exist).

    Parameters
    ----------
    enum_type: EnumType
        Enumeration type to inspect.
    enum_member_name : EnumType
        String that when uppercased yields the name of the member of this
        enumeration to be returned.

    Returns
    ----------
    EnumMemberType
        Member of this enumeration with this name uppercased.

    Raises
    ----------
    BetseEnumException
        If no member of this enumeration with this name uppercased exists.
    '''

    # It pays to be lazy. Admittedly, it doesn't pay much.
    return get_member_from_name(
        enum_type=enum_type, enum_member_name=enum_member_name.upper())

# ....................{ GETTERS ~ member : value          }....................
@type_check
def get_member_from_value(
    enum_type: EnumType, enum_member_value: object) -> EnumMemberType:
    '''
    Member of this enumeration associated with the passed arbitrary value if
    this member exists *or* raise an exception otherwise (i.e., if this member
    does *not* exist).

    Parameters
    ----------
    enum_type: EnumType
        Enumeration type to inspect.
    enum_member_value : object
        Arbitrary value associated with the member of this enumeration to be
        returned.

    Returns
    ----------
    EnumMemberType
        Member of this enumeration associated with this value.

    Raises
    ----------
    ValueError
        If no member of this enumeration with this value exists.
    '''

    # Return this member. Yes, for some obscure, inane, and unintelligible
    # reason, the "Enum" superclass performs name lookup by indexing syntax and
    # value lookup by... instantiating the "Enum" subclass. It makes utterly no
    # sense, but there it is. (Rationality has left the building.)
    #
    # Note that enumerations provide efficient mechanisms for determining
    # whether or not a member with a given name but *NOT* value exists. See the
    # standard enum.Enum.__new__() function for the supporting evidence. Ergo,
    # there intentionally exists no die_unless_member_value() function; while
    # such a function could be defined, doing so would effectively reduce to:
    #
    #     try:
    #         enum_type(enum_member_value)
    #     except ValueError:
    #         raise BetseEnumException(
    #             'Enumeration {} member with value {!r} not found.'.format(
    #                 enum_type.__name__, enum_member_value))
    #
    # To avoid such uselessly repetitious invocations of the enum_type()
    # pseudo-method, we simply call this pseudo-method as is and defer to the
    # mostly human-readable exception it already raises.
    return enum_type(enum_member_value)

# ....................{ ITERATORS                         }....................
@type_check
def iter_names(enum_type: EnumType) -> GeneratorType:
    '''
    Generator yielding the unmodified name of each member of the passed
    enumeration type in **declaration order** (i.e., the order in which these
    members were originally declared).

    Parameters
    ----------
    enum_type : EnumType
        Enumeration type to be inspected.

    Returns
    ----------
    GeneratorType
        Generator yielding the unmodified name of each member of this
        enumeration type.
    '''

    return (enum_member.name for enum_member in enum_type)


@type_check
def iter_names_lowercase(enum_type: EnumType) -> GeneratorType:
    '''
    Generator yielding the lowercased name of each member of the passed
    enumeration type in **declaration order** (i.e., the order in which these
    members were originally declared).

    Parameters
    ----------
    enum_type : EnumType
        Enumeration type to be inspected.

    Returns
    ----------
    GeneratorType
        Generator yielding the lowercased name of each member of this
        enumeration type.
    '''

    return (enum_member.name.lower() for enum_member in enum_type)
