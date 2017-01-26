#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Expression alias data descriptors specific to simulation configurations.
'''

# ....................{ IMPORTS                            }....................
from betse.util.type.cls.descriptors import expr_alias
from betse.util.type.types import type_check, ClassType

# ....................{ DESCRIPTORS                        }....................
@type_check
def conf_alias(keys: str, cls: ClassType = None) -> object:
    '''
    Expression alias **data descriptor** (i.e., object satisfying the data
    descriptor protocol) specific to simulation configurations, dynamically
    aliasing a target variable bound to instances of the class instantiating
    this descriptor to a source Python expression performing one or more key
    lookups into the dictionary loaded from a YAML-formatted simulation
    configuration file.

    Parameters
    ----------
    keys : str
        Python expression evaluating to the value of an arbitrarily nested key
        of the dictionary loaded from the current simulation configuration,
        typically consisting of one or more ``[``- and ``]``-delimited key
        lookups into this same dictionary (e.g.,
        ``['variable settings']['noise']['dynamic noise']``).
    cls: optional[ClassType]
        Either the expected type of this variable *or* a callable validating
        the values to which this variable is to be set. Defaults to ``None``,
        in which case this variable is permissively settable to any values.

    Returns
    ----------
    object
        Expression alias data descriptor as detailed above.

    See Also
    ----------
    :func:`expr_alias`
        Further details.
    '''

    return expr_alias(expr='self._config' + keys, cls=cls)
