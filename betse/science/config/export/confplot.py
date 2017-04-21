#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
YAML-backed simulation plot subconfigurations.
'''

#FIXME: Define saving-ordiented methods.

# ....................{ IMPORTS                            }....................
from betse.science.config.confabc import SimConfABC, SimConfList, conf_alias
from betse.science.config.export.confvis import (
    SimConfVisualCellsListItem, SimConfVisualCellListItem)
from betse.util.type import ints
from betse.util.type.types import type_check

# ....................{ SUBCLASSES                         }....................
class SimConfPlotAll(SimConfABC):
    '''
    YAML-backed simulation plot subconfiguration, encapsulating the
    configuration of all plots (both in- and post-simulation) parsed from the
    current YAML-formatted simulation configuration file.

    This subconfiguration saves (i.e., writes, serializes) in-memory plots to
    on-disk cache, image, and/or video files configured by this configuration.

    Attributes (After)
    ----------
    is_after_sim : bool
        ``True`` only if this configuration displays and/or saves
        post-simulation plots.
    is_after_sim_show : bool
        ``True`` only if this configuration displays post-simulation plots.
    is_after_sim_save : bool
        ``True`` only if this configuration saves post-simulation plots.

    Attributes (After : Single-cell)
    ----------
    after_sim_pipeline_cell : SimConfList
        List of all post-simulation single-cell plots to be animated. Ignored if
        :attr:``is_after_sim`` is ``False``.
    after_sim_pipeline_cells : SimConfList
        List of all post-simulation cell cluster plots to be animated. Ignored
        if :attr:``is_after_sim`` is ``False``.

    Attributes (Image)
    ----------
    image_filetype : str
        Filetype of all image files saved by this configuration. Ignored if
        :attr:`is_after_sim_save` is ``False``.
    image_dpi : int
        Dots per inch (DPI) of all image files saved by this configuration.
        Ignored if :attr:`is_after_sim_save` is ``False``.
    '''

    # ..................{ INITIALIZERS                       }..................
    def __init__(self, *args, **kwargs) -> None:

        # Initialize our superclass with all passed parameters.
        super().__init__(*args, **kwargs)

        # Encapsulate low-level lists of dictionaries with high-level wrappers.
        self.after_sim_pipeline_cell = SimConfList(
            confs=self._conf[
                'results options']['after solving'][
                'plots']['single cell pipeline'],
            conf_type=SimConfVisualCellListItem,
        )
        self.after_sim_pipeline_cells = SimConfList(
            confs=self._conf[
                'results options']['after solving'][
                'plots']['cell cluster pipeline'],
            conf_type=SimConfVisualCellsListItem,
        )

        # Validate all configured integers to be positive.
        ints.die_unless_positive(self.image_dpi)

    # ..................{ ALIASES ~ after                    }..................
    is_after_sim_save = conf_alias(
        "['results options']['after solving']['plots']['save']", bool)
    is_after_sim_show = conf_alias(
        "['results options']['after solving']['plots']['show']", bool)

    # ..................{ ALIASES ~ save                     }..................
    image_filetype = conf_alias(
        "['results options']['save']['plots']['filetype']", str)
    image_dpi = conf_alias(
        "['results options']['save']['plots']['dpi']", int)

    # ..................{ PROPERTIES ~ after                 }..................
    @property
    def is_after_sim(self) -> bool:
        return self.is_after_sim_save or self.is_after_sim_show


    @is_after_sim.setter
    @type_check
    def is_after_sim(self, is_after_sim: bool) -> None:
        self.is_after_sim_save = is_after_sim
        self.is_after_sim_show = is_after_sim
