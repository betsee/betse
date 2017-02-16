#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
YAML-backed simulation visualization subconfigurations.
'''

# ....................{ IMPORTS                            }....................
from betse.exceptions import BetseMethodUnimplementedException
from betse.science.config.confabc import SimConfListableABC, conf_alias
from betse.util.type.types import NumericTypes

# ....................{ SUBCLASSES                         }....................
class SimConfVisual(SimConfListableABC):
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

    # ..................{ SUPERCLASS                         }..................
    #FIXME: Actually implement this method.
    def default(self) -> None:
        raise BetseMethodUnimplementedException()

    # ..................{ ALIASES                            }..................
    name = conf_alias("['type']", str)

    # ..................{ ALIASES ~ colorbar                 }..................
    is_color_autoscaled = conf_alias("['colorbar']['autoscale']", bool)
    color_min = conf_alias("['colorbar']['minimum']", NumericTypes)
    color_max = conf_alias("['colorbar']['maximum']", NumericTypes)
