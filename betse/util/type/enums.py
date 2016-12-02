#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Low-level **enumeration** (i.e., enumerable types created by the `Enum()` class)
facilities.
'''

# ....................{ IMPORTS                            }....................
from betse.util.type.types import type_check, GeneratorType
from enum import Enum, EnumMeta

# ....................{ CLASSES                            }....................
class EnumOrdered(Enum):
    '''
    Enumeration whose members are comparable according to their assigned values.

    This :class:`Enum` subclass complies with the `Functional
    API<https://docs.python.org/3/library/enum.html#functional-api>`_ for
    enumerations, thus permitting comparable enumeration types to be defined
    with a single function call. See the example below.

    Oddly, the :class:`Enum` superclass does _not_ support such comparisons.
    This :class:`Enum` subclass amends this oversight, implementing four of the
    six rich comparison special methods. The fifth and sixth (i.e., `__eq__()`
    and `__ne__()`) are already implemented by the :class:`Enum` superclass and
    hence need _not_ be reimplemented here.

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


    # For efficiency, rich comparison special methods are expected to return the
    # "NotImplemented" constant rather than raise the "NotImplementedError"
    # exception when the current and passed objects are incomparible.
    def __ge__(self, other: object) -> bool:
        '''
        `NotImplemented` if the passed object is _not_ a member of the same
        enumeration as this member,  True` if the value of this member is
        greater than or equal to that of the passed member, or `False`
        otherwise.
        '''

        return (
            self.value >= other.value if self.__class__ is other.__class__ else
            NotImplemented)


    def __gt__(self, other: object) -> bool:
        '''
        `NotImplemented` if the passed object is _not_ a member of the same
        enumeration as this member,  True` if the value of this member is
        greater than that of the passed member, or `False` otherwise.
        '''

        return (
            self.value > other.value if self.__class__ is other.__class__ else
            NotImplemented)


    def __le__(self, other: object) -> bool:
        '''
        `NotImplemented` if the passed object is _not_ a member of the same
        enumeration as this member,  True` if the value of this member is
        less than or equal to that of the passed member, or `False` otherwise.
        '''

        return (
            self.value <= other.value if self.__class__ is other.__class__ else
            NotImplemented)


    def __lt__(self, other: object) -> bool:
        '''
        `NotImplemented` if the passed object is _not_ a member of the same
        enumeration as this member,  True` if the value of this member is
        less than that of the passed member, or `False` otherwise.
        '''

        return (
            self.value < other.value if self.__class__ is other.__class__ else
            NotImplemented)

# ....................{ ITERATORS                          }....................
@type_check
def iter_names_lowercase(enum: EnumMeta) -> GeneratorType:
    '''
    Generator yielding the lowercased name of each member of the passed
    enumeration type in **declaration order** (i.e., the order in which these
    members were originally declared).

    Parameters
    ----------
    enum : EnumMeta
        Enumeration type to be inspected.

    Returns
    ----------
    GeneratorType
        Generator yielding the lowercased name of each member of this
        enumeration type.
    '''

    return (enum_member.name.lower() for enum_member in enum)
