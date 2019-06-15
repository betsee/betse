#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level **tuple** (i.e., read-only sequence *always* supporting indexation
and optionally supporting named lookup) functionality.
'''

# ....................{ IMPORTS                           }....................
from betse.util.type.types import (
    type_check, ClassType, SequenceTypes, StrOrNoneTypes)
from collections import namedtuple

# ....................{ MAKERS                            }....................
@type_check
def make_named_subclass(
    # Mandatory parameters.
    class_name: str,
    item_names: SequenceTypes,

    # Optional parameters.
    doc: StrOrNoneTypes = None,
) -> ClassType:
    '''
    **Named tuple subclass** (i.e., subclass of the standard :class:`tuple`
    type whose items are accessible by both numeric indexation *and* by named
    attribute lookup ala instance and class variables) dynamically synthesized
    to have the passed class name and contain only the instance variables with
    the passed names.

    This factory function is a convenience wrapper for the :class:`namedtuple`
    function, whose non-standard semantics are arguably more obfuscatory than
    helpful.

    Parameters
    ----------
    class_name : str
        Class name of this tuple type, ideally unique across all attributes of
        the module calling this function.
    item_names : SequenceTypes
        Sequence of the names of all items of this type, required to be valid
        **Python identifiers** (i.e., contain only alphanumeric characters and
        the underscore).
    doc : StrOrNoneTypes
        Class docstring to document this type with. Defaults to ``None``, in
        which case this type remains undocumented.

    Returns
    ----------
    ClassType
        Enumeration type dynamically synthesized as defined as above.

    Examples
    ----------
        >>> from betse.util.type.iterable.tuples import make_subclass_named
        >>> OneBloodStrain = make_subclass_named(
        ...     class_name='OneBloodStrain',
        ...     item_names=('summer_goes', 'paint_leaves',))
        >>> common_man = OneBloodStrain(
        ...     summer_goes='and autumn comes',
        ...     paint_leaves='with sombre fires')
        >>> common_man.summer_goes
        'and autumn comes'
    '''

    # Avoid circular import dependencies.
    from betse.util.py import pyident
    from betse.util.type.call import callers
    from betse.util.type.iterable import itertest

    # If the passed class name is syntactically invalid, raise an exception.
    pyident.die_unless_unqualified(class_name)

    # If any passed item name is *NOT* a string, raise an exception.
    itertest.die_unless_items_instance_of(iterable=item_names, cls=str)

    # Fully-qualified name of the module defining this enumeration type.
    module_name = callers.get_caller_module_name()

    # Dynamically synthesize this named tuple subclass.
    named_subclass = namedtuple(
        typename=class_name,
        field_names=item_names,
    )

    #FIXME: Add support for passing the optional "module" parameter to the
    #namedtuple() function called above *AFTER* requiring Python >= 3.6, which
    #first introduced this parameter. Until this, we do this the hard way.

    # Set the fully-qualified name of the module declaring this subclass to
    # that of the caller. By default, the namedtuple() function sets this name
    # to that of this submodule -- fairly useless, all things considered.
    named_subclass.__module__ = module_name

    # If passed a docstring, assign this subclass this docstring.
    if doc is not None:
        named_subclass.__doc__ = doc

    # Return this subclass.
    return named_subclass
