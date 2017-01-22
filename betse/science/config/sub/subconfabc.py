#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Abstract base classes of all YAML-backed simulation subconfigurations.
'''

# ....................{ IMPORTS                            }....................
from abc import ABCMeta
from betse.util.type.types import type_check, MappingType

# ....................{ SUPERCLASSES                       }....................
class SimSubconfABC(object, metaclass=ABCMeta):
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
