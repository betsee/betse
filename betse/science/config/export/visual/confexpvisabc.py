#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
YAML-backed simulation subconfiguration classes for exporting visuals (e.g.,
plots, animations).
'''

# ....................{ IMPORTS                           }....................
from abc import ABCMeta
from betse.lib.yaml.yamlalias import yaml_alias
from betse.lib.yaml.abc.yamlabc import YamlABC
from betse.lib.yaml.abc.yamllistabc import YamlList, YamlListItemABC
from betse.science.config.export.confexpabc import SimConfExportABC
# from betse.util.io.log import logs
from betse.util.type.decorator.deccls import abstractproperty
from betse.util.type.types import type_check, NumericSimpleTypes

# ....................{ SUPERCLASSES                      }....................
#FIXME: Non-ideal. Ideally, all networks subconfigurations should be refactored
#to leverage the YAML format specified by "SimConfVisualCellsYAMLMixin".
class SimConfVisualCellsABC(object, metaclass=ABCMeta):
    '''
    Abstract base class of all cell cluster visual subconfiguration subclasses,
    YAML-backed or otherwise.

    Attributes (Colorbar)
    ----------
    color_max : NumericSimpleTypes
        Maximum color value to be displayed by this visual's colorbar. Ignored
        if :attr:`is_color_autoscaled` is ``True``.
    color_min : NumericSimpleTypes
        Minimum color value to be displayed by this visual's colorbar. Ignored
        if :attr:`is_color_autoscaled` is ``True``.
    is_color_autoscaled : bool
        ``True`` if dynamically setting the minimum and maximum colorbar values
        for this visual to the minimum and maximum values flattened from
        the corresponding time series *or* ``False`` if statically setting
        these values to :attr:`color_min` and :attr:`color_max`.
    '''

    # ..................{ PROPERTIES                        }..................
    @abstractproperty
    def is_color_autoscaled(self) -> bool:
        pass

    @abstractproperty
    def color_min(self) -> NumericSimpleTypes:
        pass

    @abstractproperty
    def color_max(self) -> NumericSimpleTypes:
        pass


class SimConfVisualCellsYAMLMixin(SimConfVisualCellsABC):
    '''
    Mixin of all **YAML-backed cell cluster visual subconfiguration** (i.e.,
    configuration of a single visual export parsed from the current
    YAML-formatted simulation configuration file) subclasses.

    For generality, this mixin provides no support for a top-level YAML key
    ``type`` or corresponding :attr:`kind` property.
    '''

    # ..................{ ALIASES ~ colorbar                }..................
    is_color_autoscaled = yaml_alias("['colorbar']['autoscale']", bool)
    color_min = yaml_alias("['colorbar']['minimum']", NumericSimpleTypes)
    color_max = yaml_alias("['colorbar']['maximum']", NumericSimpleTypes)

# ....................{ SUBCLASSES                        }....................
#FIXME: Eliminate this subclass. For serializability, all configuration classes
#should be YAML-backed.
class SimConfVisualCellsNonYAML(SimConfVisualCellsABC):
    '''
    Cell cluster visual subconfiguration specific to the
    :mod:`betse.science.networks` package.

    Unlike all cell cluster visual subconfigurations, this class preserves
    backward compatibility by implementing superclass properties *without*
    calling the :func:`yaml_alias` function returning YAML-backed data
    descriptors. Subconfiguration changes are *not* propagated back to disk.
    '''

    # ..................{ INITIALIZERS                      }..................
    @type_check
    def __init__(
        self,
        is_color_autoscaled: bool,
        color_min: NumericSimpleTypes,
        color_max: NumericSimpleTypes,
    ) -> None:

        self._is_color_autoscaled = is_color_autoscaled
        self._color_min = color_min
        self._color_max = color_max

    # ..................{ PROPERTIES                        }..................
    @property
    def is_color_autoscaled(self) -> bool:
        return self._is_color_autoscaled

    @property
    def color_min(self) -> NumericSimpleTypes:
        return self._color_min

    @property
    def color_max(self) -> NumericSimpleTypes:
        return self._color_max

# ....................{ SUBCLASSES                        }....................
class SimConfExportVisualCellsEmbedded(SimConfVisualCellsYAMLMixin, YamlABC):
    '''
    YAML-backed cell cluster visual subconfiguration, configuring a single
    visual applicable to all cells parsed from a dictionary specifying at least
    this visual in the current YAML-formatted simulation configuration file.

    This is the *only* visual configured by this dictionary. Hence, this
    dictionary contains no distinguishing ``type`` entry and this
    subconfiguration class no corresponding :attr:`name` property.
    '''

    pass

# ....................{ SUBCLASSES : item                 }....................
class SimConfExportVisualCells(
    SimConfVisualCellsYAMLMixin, SimConfExportABC):
    '''
    **Cell cluster visual export subconfiguration** (i.e., YAML-backed list
    item configuring the exportation of one or more graphical files applicable
    to all cells from the simulation configuration file containing this item).
    '''

    # ..................{ MAKERS                            }..................
    @classmethod
    @type_check
    def make_default(cls, yaml_list: YamlList) -> YamlListItemABC:

        # Duplicate the default animation listed first in our default YAML
        # configuration file.
        yaml_list_item = SimConfExportVisualCells()
        yaml_list_item.load(conf={
            'name': yaml_list.get_item_name_uniquified(
                'Cell Cluster Visual ({})'),
            'type': 'voltage_membrane',
            'enabled': True,
            'colorbar': {
                'autoscale': True,
                'minimum': -70.0,
                'maximum':  10.0,
            },
        })
        return yaml_list_item


class SimConfExportVisualCell(SimConfExportABC):
    '''
    **Single-cell visual export subconfiguration** (i.e., YAML-backed list
    item configuring the exportation of one or more graphical files applicable
    to a single cell from the simulation configuration file containing this
    item).
    '''

    # ..................{ MAKERS                            }..................
    @classmethod
    @type_check
    def make_default(cls, yaml_list: YamlList) -> YamlListItemABC:

        # Duplicate the default animation listed first in our default YAML
        # configuration file.
        yaml_list_item = SimConfExportVisualCell()
        yaml_list_item.load(conf={
            'name': yaml_list.get_item_name_uniquified(
                'Single-cell Visual ({})'),
            'type': 'voltage_membrane',
            'enabled': True,
        })
        return yaml_list_item
