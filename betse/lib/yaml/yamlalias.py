#!/usr/bin/env python3
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Simulation configuration-specific **expression alias data descriptor** (i.e.,
object satisfying the data descriptor protocol aliasing a YAML-backed simulation
configuration option to a coventional instance variable) functionality.
'''

# ....................{ IMPORTS                            }....................
from abc import ABCMeta
from betse.util.type.descriptor.expralias import expr_alias, expr_enum_alias
from betse.util.type.types import type_check

# ....................{ SUPERCLASSES                       }....................
class YamlAliasABC(object, metaclass=ABCMeta):
    '''
    Abstract base class of all simulation configuration-specific **expression
    alias data descriptor** (i.e., object satisfying the data descriptor
    protocol aliasing a YAML-backed simulation configuration option to a
    standard instance variable) subclasses.

    Motivation
    ----------
    The class principally exists for type validation -- namely, to differentiate
    the data descriptors returned by alias functions defined by this submodule
    (e.g., :func:`yaml_alias`) from those created by other means.
    '''

    pass

# ....................{ GLOBALS                            }....................
_YAML_ALIAS_BASE_CLASSES = (YamlAliasABC,)
'''
Tuple of all base classes of all simulation configuration-specific expression
alias data descriptors.

See Also
----------
:class:`YamlAliasABC`
    Further details.
'''

# ....................{ DESCRIPTORS                        }....................
@type_check
def yaml_alias(keys: str, *args, **kwargs) -> YamlAliasABC:
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

    All remaining parameters are passed as is to the underlying
    :func:`expr_alias` function.

    Returns
    ----------
    YamlAliasABC
        Expression alias data descriptor as detailed above.

    See Also
    ----------
    :func:`expr_alias`
        Further details.
    '''

    return expr_alias(
        'self._conf' + keys, *args,
        base_classes=_YAML_ALIAS_BASE_CLASSES, **kwargs)


@type_check
def yaml_enum_alias(keys: str, *args, **kwargs) -> YamlAliasABC:
    '''
    Enumeration-specific expression alias **data descriptor** (i.e., object
    satisfying the data descriptor protocol) specific to simulation
    configurations, dynamically aliasing a target variable of the passed
    enumeration type bound to instances of the class instantiating this
    descriptor to an arbitrarily complex source Python expression performing one
    or more key lookups into the dictionary loaded from a YAML-formatted
    simulation configuration file.

    Parameters
    ----------
    keys : str
        Oner or more ``[``- and ``]``-delimited key lookups. See the
        :func:`yaml_alias` function for further details.

    All remaining parameters are passed as is to the underlying
    :func:`expr_enum_alias` function.

    Returns
    ----------
    YamlAliasABC
        Enumeration-specific expression alias data descriptor as detailed above.

    See Also
    ----------
    :func:`expr_enum_alias`
        Further details.
    '''

    return expr_enum_alias(
        'self._conf' + keys, *args,
        base_classes=_YAML_ALIAS_BASE_CLASSES, **kwargs)

# ....................{ DESCRIPTORS ~ predicate : float    }....................
def yaml_alias_float_positive(keys: str) -> YamlAliasABC:
    '''
    Simulation configuration expression alias data descriptor, dynamically
    aliasing a target floating point variable with values constrained to be
    **positive** (i.e., strictly greater than 0.0).

    See Also
    ----------
    :func:`yaml_alias`
        Further details.
    '''

    return yaml_alias(
        keys=keys,
        cls=float,
        predicate_expr='value > 0.0',
        predicate_label='positive',
    )


def yaml_alias_float_nonnegative(keys: str) -> YamlAliasABC:
    '''
    Simulation configuration expression alias data descriptor, dynamically
    aliasing a target floating point variable with values constrained to be
    **non-negative** (i.e., greater than or equal to 0.0).

    See Also
    ----------
    :func:`yaml_alias`
        Further details.
    '''

    return yaml_alias(
        keys=keys,
        cls=float,
        predicate_expr='value >= 0.0',
        predicate_label='non-negative',
    )


def yaml_alias_float_percent(keys: str) -> YamlAliasABC:
    '''
    Simulation configuration expression alias data descriptor, dynamically
    aliasing a target floating point variable with values constrained to be
    **percentages** (i.e., in the range ``[0.0, 100.0]``).

    See Also
    ----------
    :func:`yaml_alias`
        Further details.
    '''

    return yaml_alias(
        keys=keys,
        cls=float,
        predicate_expr='0.0 <= value <= 100.0',
        predicate_label='percentage',
    )

# ....................{ DESCRIPTORS ~ predicate : int      }....................
def yaml_alias_int_positive(keys: str) -> YamlAliasABC:
    '''
    Simulation configuration expression alias data descriptor, dynamically
    aliasing a target integer variable with values constrained to be
    **positive** (i.e., strictly greater than 0).

    See Also
    ----------
    :func:`yaml_alias`
        Further details.
    '''

    return yaml_alias(
        keys=keys,
        cls=int,
        predicate_expr='value > 0',
        predicate_label='positive',
    )
