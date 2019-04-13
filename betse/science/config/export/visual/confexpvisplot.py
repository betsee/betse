#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
YAML-backed simulation subconfigurations for exporting plots.
'''

# ....................{ IMPORTS                           }....................
from betse.lib.yaml.yamlalias import yaml_alias, yaml_alias_int_positive
from betse.lib.yaml.abc.yamlabc import YamlABC
from betse.lib.yaml.abc.yamllistabc import YamlList, YamlListItemABC
from betse.science.config.export.confexpabc import SimConfExportABC
from betse.science.config.export.visual.confexpvisabc import (
    SimConfVisualCellsYAMLMixin)
from betse.util.type.types import type_check

# ....................{ SUBCLASSES                        }....................
class SimConfExportPlots(YamlABC):
    '''
    YAML-backed subconfiguration for exporting *all* plots (both in- and
    post-simulation) enabled by the current YAML-formatted simulation
    configuration file.

    This subconfiguration saves (i.e., writes, serializes) in-memory plots to
    on-disk cache, image, and/or video files configured by this configuration.

    Attributes (After Solving)
    ----------
    is_after_sim : bool
        ``True`` only if post-simulation plots are to be displayed and/or saved.
    is_after_sim_show : bool
        ``True`` only if post-simulation plots are to be displayed.
    is_after_sim_save : bool
        ``True`` only if post-simulation plots are to be saved.

    Attributes (After Solving: Single-cell)
    ----------
    plots_cell_after_sim : YamlList
        YAML-backed list of all post-simulation single-cell plots to be
        animated. Ignored if :attr:``is_after_sim`` is ``False``.
    plots_cells_after_sim : YamlList
        YAML-backed list of all post-simulation cell cluster plots to be
        animated. Ignored if :attr:``is_after_sim`` is ``False``.

    Attributes (Image)
    ----------
    image_filetype : str
        Filetype of all image files saved by this configuration. Ignored if
        :attr:`is_after_sim_save` is ``False``.
    image_dpi : int
        Dots per inch (DPI) of all image files saved by this configuration.
        Ignored if :attr:`is_after_sim_save` is ``False``.
    '''

    # ..................{ ALIASES ~ after                   }..................
    is_after_sim_save = yaml_alias(
        "['results options']['after solving']['plots']['save']", bool)
    is_after_sim_show = yaml_alias(
        "['results options']['after solving']['plots']['show']", bool)

    # ..................{ ALIASES ~ save                    }..................
    image_filetype = yaml_alias(
        "['results options']['save']['plots']['filetype']", str)
    image_dpi = yaml_alias_int_positive(
        "['results options']['save']['plots']['dpi']")

    # ..................{ INITIALIZERS                      }..................
    def __init__(self, *args, **kwargs) -> None:

        # Initialize our superclass with all passed parameters.
        super().__init__(*args, **kwargs)

        # Encapsulate low-level lists of dictionaries with high-level wrappers.
        self.plots_cell_after_sim  = SimConfExportPlotCell .make_list()
        self.plots_cells_after_sim = SimConfExportPlotCells.make_list()

    # ..................{ LOADERS                           }..................
    def load(self, *args, **kwargs) -> None:

        # Load our superclass with all passed arguments.
        super().load(*args, **kwargs)

        # Simulation subconfigurations localized for convenience.
        plot_pipelines = self._conf[
            'results options']['after solving']['plots']

        # Load all subconfigurations of this configuration.
        self.plots_cell_after_sim.load(
            conf=plot_pipelines['single cell pipeline'])
        self.plots_cells_after_sim.load(
            conf=plot_pipelines['cell cluster pipeline'])


    def unload(self) -> None:

        # Unload our superclass.
        super().unload()

        # Unload all subconfigurations of this configuration.
        self.plots_cell_after_sim.unload()
        self.plots_cells_after_sim.unload()

    # ..................{ PROPERTIES ~ after                }..................
    @property
    def is_after_sim(self) -> bool:
        return self.is_after_sim_save or self.is_after_sim_show


    @is_after_sim.setter
    @type_check
    def is_after_sim(self, is_after_sim: bool) -> None:
        self.is_after_sim_save = is_after_sim
        self.is_after_sim_show = is_after_sim

# ....................{ SUBCLASSES ~ item                 }....................
class SimConfExportPlotCell(SimConfExportABC):
    '''
    **Exported single-cell plot subconfiguration** (i.e., YAML-backed list item
    configuring the exportation of one or more images specific to a single cell
    from the simulation configuration file containing this item).
    '''

    # ..................{ SUPERCLASS                        }..................
    @classmethod
    @type_check
    def make_default(cls, yaml_list: YamlList) -> YamlListItemABC:

        # Duplicate the first plot in our default configuration file.
        return cls._make_loaded(conf={
            'name': yaml_list.get_item_name_uniquified(
                'Single Cell Plot ({})'),
            'type': 'voltage_membrane',
            'enabled': True,
        })


class SimConfExportPlotCells(
    SimConfVisualCellsYAMLMixin, SimConfExportABC):
    '''
    **Exported cell cluster plot subconfiguration** (i.e., YAML-backed list
    item configuring the exportation of one or more images applicable to all
    cells from the simulation configuration file containing this item).
    '''

    # ..................{ SUPERCLASS                        }..................
    @classmethod
    @type_check
    def make_default(cls, yaml_list: YamlList) -> YamlListItemABC:

        # Duplicate the first plot in our default configuration file.
        return cls._make_loaded(conf={
            'name': yaml_list.get_item_name_uniquified(
                'Cell Cluster Plot ({})'),
            'type': 'voltage_membrane',
            'enabled': True,
            'colorbar': {
                'autoscale': True,
                'minimum': -70.0,
                'maximum':  10.0,
            },
        })
