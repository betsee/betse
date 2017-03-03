#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
YAML-backed simulation visualization subconfigurations.
'''

# ....................{ IMPORTS                            }....................
#from abc import ABCMeta
from betse.exceptions import BetseMethodUnimplementedException
from betse.science.config.confabc import (
    SimConfABC, SimConfListableABC, conf_alias)
from betse.util.type.types import type_check, NumericTypes

# ....................{ SUPERCLASSES                       }....................
#FIXME: Non-ideal. Ideally, all networks subconfigurations should be refactored
#to leverage the YAML format specified by "SimConfVisualMixin".
class SimConfVisualABC(object):
    '''
    Abstract base class generalizing logic common to all YAML-backed simulation
    visualization subconfiguration subclasses.

    This mix-in encapsulates the configuration of a single visualization (either
    in- or post-simulation plot or animation) parsed from the current
    YAML-formatted simulation configuration file. For generality, this mix-in
    provides no support for a YAML ``type`` key or corresponding :attr:`name`
    property.

    Attributes (Colorbar)
    ----------
    color_max : NumericTypes
        Maximum color value to be displayed by this visualization's colorbar.
        Ignored if :attr:`is_color_autoscaled` is ``True``.
    color_min : NumericTypes
        Minimum color value to be displayed by this visualization's colorbar.
        Ignored if :attr:`is_color_autoscaled` is ``True``.
    is_color_autoscaled : bool
        ``True`` if dynamically setting the minimum and maximum colorbar values
        for this visualization to the minimum and maximum values flattened from
        the corresponding time series *or* ``False`` if statically setting these
        values to :attr:`color_min` and :attr:`color_max`.
    '''

    # ..................{ PROPERTIES                         }..................
    #FIXME: Docstring us up.
    @property
    def is_color_autoscaled(self) -> bool:
        raise BetseMethodUnimplementedException()

    @property
    def color_min(self) -> NumericTypes:
        raise BetseMethodUnimplementedException()

    @property
    def color_max(self) -> NumericTypes:
        raise BetseMethodUnimplementedException()


class SimConfVisualMixin(SimConfVisualABC):
    '''
    Abstract mix-in generalizing logic common to all YAML-backed simulation
    visualization subconfiguration subclasses.

    This mix-in encapsulates the configuration of a single visualization (either
    in- or post-simulation plot or animation) parsed from the current
    YAML-formatted simulation configuration file. For generality, this mix-in
    provides no support for a YAML ``type`` key or corresponding :attr:`name`
    property.

    Attributes (Colorbar)
    ----------
    color_max : NumericTypes
        Maximum color value to be displayed by this visualization's colorbar.
        Ignored if :attr:`is_color_autoscaled` is ``True``.
    color_min : NumericTypes
        Minimum color value to be displayed by this visualization's colorbar.
        Ignored if :attr:`is_color_autoscaled` is ``True``.
    is_color_autoscaled : bool
        ``True`` if dynamically setting the minimum and maximum colorbar values
        for this visualization to the minimum and maximum values flattened from
        the corresponding time series *or* ``False`` if statically setting these
        values to :attr:`color_min` and :attr:`color_max`.
    '''

    # ..................{ ALIASES ~ colorbar                 }..................
    is_color_autoscaled = conf_alias("['colorbar']['autoscale']", bool)
    color_min = conf_alias("['colorbar']['minimum']", NumericTypes)
    color_max = conf_alias("['colorbar']['maximum']", NumericTypes)


class SimConfVisualMolecule(SimConfVisualABC):
    '''
    Networks molecule-specific visualization subconfiguration subclasses.
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(
        self,
        is_color_autoscaled: bool,
        color_min: NumericTypes,
        color_max: NumericTypes,
    ) -> None:

        self._is_color_autoscaled = is_color_autoscaled
        self._color_min = color_min
        self._color_max = color_max

    # ..................{ PROPERTIES                         }..................
    @property
    def is_color_autoscaled(self) -> bool:
        return self._is_color_autoscaled

    @property
    def color_min(self) -> NumericTypes:
        return self._color_min

    @property
    def color_max(self) -> NumericTypes:
        return self._color_max

# ....................{ SUBCLASSES                         }....................
class SimConfVisualGeneric(SimConfVisualMixin, SimConfABC):
    '''
    YAML-backed simulation visualization subconfiguration, encapsulating the
    configuration of a single visualization (either in- or post-simulation plot
    or animation) with no ``type`` entry or corresponding :attr:`name` property
    parsed from the current YAML-formatted simulation configuration file.

    Attributes (Colorbar)
    ----------
    color_max : NumericTypes
        Maximum color value to be displayed by this visualization's colorbar.
        Ignored if :attr:`is_color_autoscaled` is ``True``.
    color_min : NumericTypes
        Minimum color value to be displayed by this visualization's colorbar.
        Ignored if :attr:`is_color_autoscaled` is ``True``.
    is_color_autoscaled : bool
        ``True`` if dynamically setting the minimum and maximum colorbar values
        for this visualization to the minimum and maximum values flattened from
        the corresponding time series *or* ``False`` if statically setting these
        values to :attr:`color_min` and :attr:`color_max`.
    '''

    pass


class SimConfVisualListable(SimConfVisualMixin, SimConfListableABC):
    '''
    YAML-backed simulation visualization subconfiguration, encapsulating the
    configuration of a single visualization (either in- or post-simulation plot
    or animation) parsed from the list of all such visualizations in the current
    YAML-formatted simulation configuration file.

    Attributes (General)
    ----------
    name : str
        Lowercase alphanumeric string uniquely identifying the type of this
        visualization (e.g., ``voltage_intra``, signifying a visualization of
        intracellular voltages). See the corresponding entry ``type`` of the
        default simulation configuration file for further commentary.
    '''

    # ..................{ CLASS                              }..................
    @classmethod
    def make_default(self) -> SimConfListableABC:

        # Duplicate the default animation listed first in our default YAML file.
        return SimConfVisualListable(conf={
            'type': 'voltage_intra',
            'colorbar': {
                'autoscale': True,
                'minimum': -70.0,
                'maximum':  10.0,
            },
        })

    # ..................{ ALIASES                            }..................
    name = conf_alias("['type']", str)
