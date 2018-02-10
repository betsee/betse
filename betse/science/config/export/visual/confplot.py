#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
YAML-backed simulation subconfigurations for exporting plots.
'''

#FIXME: Define saving-ordiented methods.

# ....................{ IMPORTS                            }....................
from betse.lib.yaml.yamlalias import yaml_alias, yaml_alias_int_positive
from betse.lib.yaml.abc.yamlabc import YamlABC
from betse.science.config.export.visual.confvisabc import (
    SimConfVisualCellsListItem, SimConfVisualCellListItem)
from betse.util.type.types import type_check

# ....................{ SUBCLASSES                         }....................
#FIXME: Rename to "SimConfPlots" for readability.
class SimConfPlotAll(YamlABC):
    '''
    YAML-backed subconfiguration for exporting *all* plots (both in- and
    post-simulation) enabled by the current YAML-formatted simulation
    configuration file.

    This subconfiguration saves (i.e., writes, serializes) in-memory plots to
    on-disk cache, image, and/or video files configured by this configuration.

    Attributes (After Solving)
    ----------
    is_after_sim : bool
        ``True`` only if this configuration displays and/or saves
        post-simulation plots.
    is_after_sim_show : bool
        ``True`` only if this configuration displays post-simulation plots.
    is_after_sim_save : bool
        ``True`` only if this configuration saves post-simulation plots.

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

    # ..................{ ALIASES ~ after                    }..................
    is_after_sim_save = yaml_alias(
        "['results options']['after solving']['plots']['save']", bool)
    is_after_sim_show = yaml_alias(
        "['results options']['after solving']['plots']['show']", bool)

    # ..................{ ALIASES ~ save                     }..................
    image_filetype = yaml_alias(
        "['results options']['save']['plots']['filetype']", str)
    image_dpi = yaml_alias_int_positive(
        "['results options']['save']['plots']['dpi']")

    # ..................{ INITIALIZERS                       }..................
    def __init__(self, *args, **kwargs) -> None:

        # Initialize our superclass with all passed parameters.
        super().__init__(*args, **kwargs)

        # Encapsulate low-level lists of dictionaries with high-level wrappers.
        self.plots_cell_after_sim = SimConfVisualCellListItem.make_list()
        self.plots_cells_after_sim = SimConfVisualCellsListItem.make_list()

    # ..................{ LOADERS                            }..................
    def load(self, *args, **kwargs) -> None:

        # Load our superclass with all passed arguments.
        super().load(*args, **kwargs)

        # Load all subconfigurations of this configuration.
        self.plots_cell_after_sim.load(conf=self._conf[
            'results options']['after solving'][
            'plots']['single cell pipeline'])
        self.plots_cells_after_sim.load(conf=self._conf[
            'results options']['after solving'][
            'plots']['cell cluster pipeline'])


    def unload(self) -> None:

        # Unload our superclass.
        super().unload()

        # Unload all subconfigurations of this configuration.
        self.plots_cell_after_sim.unload()
        self.plots_cells_after_sim.unload()

    # ..................{ PROPERTIES ~ after                 }..................
    @property
    def is_after_sim(self) -> bool:
        return self.is_after_sim_save or self.is_after_sim_show


    @is_after_sim.setter
    @type_check
    def is_after_sim(self, is_after_sim: bool) -> None:
        self.is_after_sim_save = is_after_sim
        self.is_after_sim_show = is_after_sim
