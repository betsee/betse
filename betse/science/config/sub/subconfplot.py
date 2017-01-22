#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
YAML-backed simulation plot subconfigurations.
'''

#FIXME: Define saving-ordiented methods.

# ....................{ IMPORTS                            }....................
from betse.science.config.sub.subconfabc import SimSubconfABC
from betse.util.type import ints
#from betse.util.type.types import type_check

# ....................{ SUBCLASSES                         }....................
class SimSubconfPlot(SimSubconfABC):
    '''
    YAML-backed simulation plot subconfiguration, encapsulating both the
    configuration and writing of all plots (both mid- and post-simulation)
    parsed from the current YAML-formatted simulation configuration file.

    This subconfiguration saves (i.e., writes, serializes) in-memory plots to
    on-disk cache, image, and/or video files configured by this configuration.

    Attributes
    ----------
    is_after_sim : bool
        `True` only if this configuration enables (but _not_ necessarily
        displays or saves) post-simulation plots.
    is_after_sim_show : bool
        `True` only if this configuration displays post-simulation plots.
        Ignored if `is_after_sim` is `False`.
    is_after_sim_save : bool
        `True` only if this configuration saves post-simulation plots.
        Ignored if `is_after_sim` is `False`.
    image_filetype : str
        Filetype of all image files saved by this configuration. Ignored if
        `is_after_sim_save` is `False`.
    image_dpi : int
        Dots per inch (DPI) of all image files saved by this configuration.
        Ignored if `is_after_sim_save` is `False`.
    '''

    # ..................{ INITIALIZERS                       }..................
    def __init__(self, *args, **kwargs) -> None:

        # Initialize our superclass with all passed parameters.
        super().__init__(*args, **kwargs)

        # For convenience, localize configuration subdictionaries.
        results = self._config['results options']
        after_solving = results['after solving']['plots']
        image =         results['save']['plots']

        # Post-simulation plots.
        self.is_after_sim =      after_solving['enabled']
        self.is_after_sim_show = after_solving['show']
        self.is_after_sim_save = after_solving['save']

        # Image saving.
        self.image_filetype = image['filetype']
        self.image_dpi =      image['dpi']

        # Validate all passed integers to be positive.
        ints.die_unless_positive(self.image_dpi)
