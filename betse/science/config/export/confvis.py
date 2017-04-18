#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
YAML-backed simulation visual subconfigurations.
'''

# ....................{ IMPORTS                            }....................
from abc import ABCMeta, abstractproperty
from betse.science.config.confabc import (
    SimConfABC, SimConfListableABC, SimConfListItemTypedABC, conf_alias)
from betse.util.io.log import logs
from betse.util.type.types import type_check, NumericTypes

# ....................{ SUPERCLASSES                       }....................
#FIXME: Non-ideal. Ideally, all networks subconfigurations should be refactored
#to leverage the YAML format specified by "SimConfVisualCellsYAMLMixin".
class SimConfVisualCellsABC(object, metaclass=ABCMeta):
    '''
    Abstract base class generalizing logic common to all cell cluster visual
    subconfigurations -- YAML-backed or otherwise.

    Attributes (Colorbar)
    ----------
    color_max : NumericTypes
        Maximum color value to be displayed by this visual's colorbar.
        Ignored if :attr:`is_color_autoscaled` is ``True``.
    color_min : NumericTypes
        Minimum color value to be displayed by this visual's colorbar.
        Ignored if :attr:`is_color_autoscaled` is ``True``.
    is_color_autoscaled : bool
        ``True`` if dynamically setting the minimum and maximum colorbar values
        for this visual to the minimum and maximum values flattened from
        the corresponding time series *or* ``False`` if statically setting these
        values to :attr:`color_min` and :attr:`color_max`.
    '''

    # ..................{ PROPERTIES                         }..................
    @abstractproperty
    def is_color_autoscaled(self) -> bool:
        pass

    @abstractproperty
    def color_min(self) -> NumericTypes:
        pass

    @abstractproperty
    def color_max(self) -> NumericTypes:
        pass


class SimConfVisualCellsYAMLMixin(SimConfVisualCellsABC):
    '''
    Abstract mix-in generalizing logic common to all YAML-backed cell cluster
    visual subconfigurations.

    This mix-in encapsulates the configuration of a single visual (either
    in- or post-simulation plot or animation) parsed from the current
    YAML-formatted simulation configuration file. For generality, this mix-in
    provides no support for a YAML ``type`` key or corresponding :attr:`name`
    property.
    '''

    # ..................{ ALIASES ~ colorbar                 }..................
    is_color_autoscaled = conf_alias("['colorbar']['autoscale']", bool)
    color_min = conf_alias("['colorbar']['minimum']", NumericTypes)
    color_max = conf_alias("['colorbar']['maximum']", NumericTypes)

# ....................{ SUBCLASSES                         }....................
#FIXME: Eliminate this subclass. For serializability, all configuration classes
#should be YAML-backed.
class SimConfVisualCellsNonYAML(SimConfVisualCellsABC):
    '''
    Cell cluster visual subconfiguration specific to the
    :mod:`betse.science.networks` package.

    Unlike all cell cluster visual subconfigurations, this class preserves
    backward compatibility by implementing superclass properties *without*
    calling the :func:`conf_alias` function returning YAML-backed data
    descriptors. Subconfiguration changes are *not* propagated back to disk.
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
class SimConfVisualCellsEmbedded(SimConfVisualCellsYAMLMixin, SimConfABC):
    '''
    YAML-backed cell cluster visual subconfiguration, encapsulating the
    configuration of a single visual applicable to all cells parsed from a
    dictionary configuring at least this visual in the current
    YAML-formatted simulation configuration file.

    This is the *only* visual configured by this dictionary. Hence, this
    dictionary contains no distinguishing ``type`` entry and this
    subconfiguration class no corresponding :attr:`name` property.
    '''

    pass

# ....................{ SUBCLASSES : list item             }....................
class SimConfVisualCellsListItem(
    SimConfVisualCellsYAMLMixin, SimConfListItemTypedABC):
    '''
    YAML-backed cell cluster visual subconfiguration, encapsulating the
    configuration of a single visual applicable to all cells parsed from a list
    of these visuals in the current YAML-formatted simulation configuration file.
    '''

    # ..................{ CLASS                              }..................
    @classmethod
    def make_default(self) -> SimConfListableABC:

        # Duplicate the default animation listed first in our default YAML file.
        return SimConfVisualCellsListItem(conf={
            'type': 'voltage_membrane',
            'enabled': True,
            'colorbar': {
                'autoscale': True,
                'minimum': -70.0,
                'maximum':  10.0,
            },
        })

    # ..................{ INITIALIZERS                       }..................
    def __init__(self, *args, **kwargs) -> None:

        # Initialize our superclass with all passed parameters.
        super().__init__(*args, **kwargs)

        # If newly defined configuration options are undefined...
        if 'enabled' not in self._conf:
            # Log a non-fatal warning.
            logs.log_warning(
                'Animation "%s" config setting "enabled" not found. '
                'Repairing to preserve backward compatibility. '
                'Consider upgrading to the newest config file format!',
                self.name)

            # Preserve backwards compatibility.
            self._conf['enabled'] = True

        # If this is an obselete visual type, silently update this to the
        # corresponding new visual type.
        if self.name == 'polarization':
            self.name = 'voltage_polarity'


class SimConfVisualCellListItem(SimConfListItemTypedABC):
    '''
    YAML-backed single-cell visual subconfiguration, encapsulating the
    configuration of a single visual (either in- or post-simulation
    plot or animation) specific to a single cell parsed from the list of all
    such visuals in the current YAML-formatted simulation configuration file.
    '''

    # ..................{ CLASS                              }..................
    @classmethod
    def make_default(self) -> SimConfListableABC:

        # Duplicate the default animation listed first in our default YAML file.
        return SimConfVisualCellListItem(conf={
            'type': 'voltage_membrane',
            'enabled': True,
        })
