#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Abstract base classes of all YAML-backed simulation configuration subclasses as
well as functionality pertaining to such classes.
'''

# ....................{ IMPORTS                            }....................
from abc import ABCMeta
from betse.util.type.cls.descriptors import expr_alias
from betse.util.type.types import type_check, TestableTypes, MappingType

# ....................{ DESCRIPTORS                        }....................
class SimConfABC(object, metaclass=ABCMeta):
    '''
    Abstract base class of all YAML-backed simulation subconfigurations,
    encapsulating both the configuration and writing of a single subsection
    (e.g., tissue profiles, animations) of the current YAML-formatted simulation
    configuration file.

    Attributes
    ----------
    _config : MappingType
        Dictionary describing the current simulation configuration, returned by
        the :func:`betse.lib.yaml.yamls.load` function and thereafter passed to
        the :func:`betse.lib.yaml.yamls.save` function.
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(self, config: MappingType) -> None:
        '''
        Initialize this subconfiguration.

        Attributes
        ----------
        config : MappingType
            Dictionary describing the current simulation subconfiguration,
            returned by the :func:`betse.lib.yaml.yamls.load` function and then
            received by the :func:`betse.lib.yaml.yamls.save` function.
        '''

        # Classify the passed parameters.
        self._config = config

# ....................{ DESCRIPTORS                        }....................
@type_check
def conf_alias(keys: str, cls: TestableTypes = None) -> object:
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
        Either the expected type or tuple of the expected types of this variable
        *or* a callable validating the values to which this variable is to be
        set. Defaults to ``None``, in which case this variable is permissively
        settable to *any* values.

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
