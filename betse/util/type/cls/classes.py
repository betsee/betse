#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level object facilities.
'''

# ....................{ IMPORTS                            }....................
from betse.util.type.types import (
    type_check, ClassType, MappingType, SequenceTypes)

# ....................{ TESTERS                            }....................
@type_check
def define_class(
    class_name: str,
    class_attr_name_to_value: MappingType = {},
    base_classes: SequenceTypes = (),
) -> ClassType:
    '''
    Dynamically define a new class with the passed name subclassing all passed
    base classes and providing all passed class attributes (e.g., class
    variables, methods).

    Parameters
    ----------
    class_name : str
        Name of the class to be created.
    class_attr_name_to_value : MappingType
        Mapping from the name to the initial value of each class attribute
        (e.g., class variable, method) to declare this class to contain.
        Defaults to the empty dictionary, equivalent to declaring a class with
        the trivial body ``pass``.
    base_classes : optional[SequenceTypes]
        Sequence of all base classes to subclass this class from. Defaults to
        the empty tuple, equivalent to the 1-tuple ``(object,)`` containing only
        the root base class of all classes.

    Returns
    ----------
    ClassType
        Class dynamically defined with this name from these base classes and
        class attributes.
    '''

    # Thank you, bizarre 3-parameter variant of the type.__init__() method.
    return type(class_name, base_classes, class_attr_name_to_value)
