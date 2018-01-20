#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level **enumeration** (i.e., enumerable types created by the `Enum()` class)
facilities.
'''

# ....................{ IMPORTS                            }....................
from betse.exceptions import BetseEnumException
from betse.util.type.types import (
    type_check,
    CallableTypes,
    EnumMemberType,
    EnumType,
    EnumClassType,
    GeneratorType,
    SequenceTypes,
)

# ....................{ SUBCLASSES                         }....................
class EnumOrdered(EnumClassType):
    '''
    Enumeration whose members are comparable according to their assigned values.

    This :class:`Enum` subclass complies with the `Functional
    API<https://docs.python.org/3/library/enum.html#functional-api>`_ for
    enumerations, permitting comparable enumeration types to be defined with a
    single function call. (See the example below.)

    Oddly, the :class:`Enum` superclass does *not* support such comparisons.
    This :class:`Enum` subclass amends this oversight, implementing four of the
    six rich comparison special methods. The fifth and sixth (i.e., ``__eq__()``
    and ``__ne__()``) are already implemented by the :class:`Enum` superclass
    and hence need *not* be reimplemented here.

    See Also
    ----------
    https://docs.python.org/3/library/enum.html#orderedenum
        `OrderedEnum` class strongly inspiring this class.

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

    # ..................{ COMPARATORS                        }..................
    # For efficiency, rich comparison special methods are expected to return the
    # "NotImplemented" constant rather than raise the "NotImplementedError"
    # exception when the current and passed objects are incomparible.

    def __ge__(self, other: object) -> bool:
        '''
        :data:`NotImplemented` if the passed object is *not* a member of the
        same enumeration as this member, ``True`` if the value of this member is
        greater than or equal to that of the passed member, or ``False``
        otherwise.
        '''

        return (
            self.value >= other.value if self.__class__ is other.__class__ else
            NotImplemented)


    def __gt__(self, other: object) -> bool:
        '''
        :data:`NotImplemented` if the passed object is *not* a member of the
        same enumeration as this member, ``True`` if the value of this member is
        greater than that of the passed member, or ``False`` otherwise.
        '''

        return (
            self.value > other.value if self.__class__ is other.__class__ else
            NotImplemented)


    def __le__(self, other: object) -> bool:
        '''
        :data:`NotImplemented` if the passed object is *not* a member of the
        same enumeration as this member, ``True`` if the value of this member is
        less than or equal to that of the passed member, or ``False`` otherwise.
        '''

        return (
            self.value <= other.value if self.__class__ is other.__class__ else
            NotImplemented)


    def __lt__(self, other: object) -> bool:
        '''
        :data:`NotImplemented` if the passed object is *not* a member of the
        same enumeration as this member, ``True`` if the value of this member is
        less than that of the passed member, or ``False`` otherwise.
        '''

        return (
            self.value < other.value if self.__class__ is other.__class__ else
            NotImplemented)

# ....................{ MAKERS                             }....................
@type_check
def make_enum(class_name: str, member_names: SequenceTypes) -> EnumType:
    '''
    **Enumeration type** (i.e., :class:`EnumType` instance or, equivalently,
    :class:`EnumClassType` subclass) dynamically synthesized to have the passed
    class name and contain exactly the members with the passed names.

    This factory function is a convenience wrapper for the
    :class:`EnumClassType` subclass, whose non-standard semantics are arguably
    more obfuscatory than helpful.

    Parameters
    ----------
    class_name : str
        Class name of this enumeration type, ideally unique across all
        attributes of the module calling this function.
    member_names : SequenceTypes
        Sequence of the names of all members of this enumeration type, required
        to be valid Python identifiers (i.e., contain only alphanumeric
        characters and the underscore)

    Returns
    ----------
    EnumType
        Enumeration type dynamically synthesized as defined as above.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.call import callers

    # Fully-qualified name of the module defining this enumeration type.
    module_name = callers.get_caller_module_name()

    # Dynamically synthesize and return this enumeration type.
    return EnumClassType(
        value=class_name,
        names=member_names,
        module=module_name,
        qualname='{}.{}'.format(module_name, class_name),
    )

# ....................{ EXCEPTIONS                         }....................
def die_unless_member(
    enum_type: EnumType, enum_member: EnumMemberType) -> None:
    '''
    Raise an exception unless the passed enumeration contains the passed member.

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

# ....................{ TESTERS                            }....................
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

# ....................{ GETTERS                            }....................
@type_check
def get_converter_name_to_uppercase_enum_member(
    enum_type: EnumType) -> CallableTypes:
    '''
    Callable accepting one string parameter returning the member of the passed
    enumeration whose name is that string uppercased, raising an exception if
    this string uppercased is *not* the name of a member of this enumeration.

    The callable returned by this function is principally intended to be passed
    as the `type` parameter to the :meth:`ArgumentParser.add_argument` method,
    converting from lowercase command-line option string arguments to
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

    # Closure specific to this enumeration type implementing this converesion.
    @type_check
    def _to_enum(enum_member_name: str) -> EnumMemberType:
        '''
        Member of this enumeration whose name is the passed string uppercased if
        such a member exists _or_ raise an exception otherwise.

        Parameters
        ----------
        enum_member_name : EnumType
            Name of the member to be returned.

        Returns
        ----------
        EnumMemberType
            Member of this enumeration with this name uppercased.
        '''

        # Uppercase name of this member.
        enum_member_name_upper = enum_member_name.upper()

        # Raise an exception unless this member exists.
        die_unless_member_name(enum_type, enum_member_name_upper)

        # Return this member.
        return enum_type[enum_member_name_upper]

    # Return this closure.
    return _to_enum


@type_check
def get_names_lowercase(enum_type: EnumType) -> SequenceTypes:
    '''
    Sequence of the lowercased names of all members of the passed enumeration
    type in **declaration order** (i.e., the order in which these members were
    originally declared).

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

# ....................{ ITERATORS                          }....................
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
