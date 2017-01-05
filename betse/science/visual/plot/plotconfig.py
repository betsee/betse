#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Plot configuration and serialization classes.
'''

#FIXME: Define saving-ordiented methods.

# ....................{ IMPORTS                            }....................
from betse.util.type import ints, types
from betse.util.type.types import type_check

# ....................{ SUPERCLASS                         }....................
class PlotConfig(object):
    '''
    Object encapsulating both the configuration and writing of all plots
    (both in- and post-simulation), parsed from the current configuration file.

    This object saves (i.e., writes, serializes) in-memory plots to on-disk
    cache, image, and/or video files configured by this configuration.

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

    # ..................{ ABSTRACT ~ static                  }..................
    @staticmethod
    def make(params: 'Parameters') -> "PlotConfig":
        '''
        Factory method producing an instance of this class encapsulating the
        passed simulation configuration.

        Parameters
        ----------------------------
        params : Parameters
            Current simulation configuration.

        Returns
        ----------------------------
        PlotConfig
            Instance of this class encapsulating this configuration.
        '''
        assert types.is_parameters(params), types.assert_not_parameters(params)

        # For convenience, localize configuration subdictionaries.
        results = params.config['results options']
        after_solving = results['after solving']['plots']
        image = results['save']['plots']

        # Create and return this instance.
        return PlotConfig(
            # Post-simulation plots.
            is_after_sim=     after_solving['enabled'],
            is_after_sim_show=after_solving['show'],
            is_after_sim_save=after_solving['save'],

            # Image saving.
            image_filetype=image['filetype'],
            image_dpi=     image['dpi'],
        )

    # ..................{ CONCRETE ~ public                  }..................
    @type_check
    def __init__(
        self,

        # Post-simulation plots.
        is_after_sim: bool,
        is_after_sim_show: bool,
        is_after_sim_save: bool,

        # Image saving.
        image_filetype: str,
        image_dpi: int,
    ) -> None:

        # Validate all passed integers to be positive.
        ints.die_unless_positive(image_dpi)

        #FIXME: Repetition is vile and demeaning. Design and leverage a new
        #@classify_params decorator here instead, please.

        # Classify the passed parameters.
        self.is_after_sim = is_after_sim
        self.is_after_sim_show = is_after_sim_show
        self.is_after_sim_save = is_after_sim_save
        self.image_filetype = image_filetype
        self.image_dpi = image_dpi
