#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Matplotlib-specific **colormap** (i.e., objects mapping from low-level
scientific data values to high-level RGBA colour values displaying that data)
facilities.

See Also
----------
https://matplotlib.org/examples/color/colormaps_reference.html
    Canonical colormap example visually depicting all default colormaps bundled
    with matplotlib.
'''

#FIXME: Document all application-specific colormaps registered by this submodule
#in our default "sim_config.yaml" file.

# ....................{ IMPORTS                            }....................
import numpy as np
from matplotlib import cm as colormaps
from matplotlib.colors import Colormap, LinearSegmentedColormap
from betse.exceptions import BetseSequenceException
from betse.util.io.log import logs
from betse.util.type import sequences
from betse.util.type.numeric import ints
from betse.util.type.types import type_check, SequenceTypes

# ....................{ ADDERS                             }....................
@type_check
def add_colormap(
    colormap_name: str, colors: SequenceTypes) -> Colormap:
    '''
    Create a colormap with the passed name and color scheme, register this
    colormap with matplotlib, and return this colormap.

    Parameters
    -----------
    colormap_name : str
        Name of the colormap to be created.
    colors : SequenceTypes
        Two-dimensional sequence whose:
        * First dimension indexes each color defining this colormap's gradient.
          This dimension *must* be a sequence containing two or more colors.
        * Second dimension is a 3-sequence indexing the red, green, and blue
          (RGB) values defining this color.

    Raises
    -----------
    BetseSequenceException
        If either:
        * The first dimension of ``colors`` contains less than two colors
          (i.e., is either empty *or* contains only one color).
        * Any second dimension of ``colors`` does *not* contain exactly
          three color values.
    BetseIntException
        If any color value in any second dimension of ``colors`` is *not*
        a valid color value in the range ``[0, 255]``.
    '''

    #FIXME: Raise an exception if this colormap name conflicts with that
    #of an existing colormap -- probably by defining a new die_if_colormap()
    #utility function in this submodule.

    # If this colormap defines less than two colors, raise an exception.
    if len(colors) < 2:
        raise BetseSequenceException(
            'Colormap scheme defines less than two colors: {!r}'.format(
                colors))

    # For each color defining this colormap...
    for color in colors:
        # If this color is *NOT* a 3-sequence, raise an exception.
        sequences.die_unless_len(sequence=color, sequence_len=3)

        # If any of this color's values is *NOT* a valid RGB colour (i.e.,
        # integer in the range [0, 255]), raise an exception.
        ints.die_unless_byte(*color)

    # Two-dimensional Numpy array, normalizing each of each color's values from
    # [0, 255] to [0.0, 1.0] (while preserving the order of colors).
    colors_normalized = np.array(colors) / 255

    # Log this registration attempt.
    # logs.log_debug(
    #     'Registering colormap "%s": %s', colormap_name, colors_normalized)

    # Colormap synthesized from this colormap name and colors.
    #
    # Unfortunately, as the names of the first two parameters accepted by this
    # function have changed across matplotlib versions, these parameters *MUST*
    # be passed positionally for safety.
    colormap = LinearSegmentedColormap.from_list(
        colormap_name, colors_normalized, N=256)

    # Register this colormap with matplotlib. Again, for portability, the most
    # succinct variant of this function is intentionally called.
    colormaps.register_cmap(colormap_name, colormap)

    # Return this colormap.
    return colormap

# ....................{ INITIALIZERS                       }....................
#FIXME: Call this function at matplotlib initialization time.
def init() -> None:
    '''
    Initialize this module by registering all application-specific colormaps
    with matplotlib, enabling these colormaps to be trivially retrieved with the
    standard :func:`matplotlib.cm.get_cmap` function.

    This function is intended to be called at matplotlib initialization time.
    '''

    # Log this initialization attempt.
    logs.log_info('Registering custom matplotlib colormaps...')

    # ..................{ COLOURS                            }..................
    # Primary colours.
    BLACK = (  0,   0,   0)
    GREEN = (  0, 255,   0)
    RED   = (255,   0,   0)
    BLUE  = (  0,   0, 255)

    # Secondary colours.
    CYAN    = (  9, 232, 239)
    MAGENTA = (239,  52, 236)
    # ORANGE  = (255, 164,  61)

    GREY_DARK = (51, 51, 51)

    GREEN_LIGHT  = (184, 255, 104)
    PURPLE_LIGHT = (219, 104, 255)

    BLUE_PALE = ( 56, 132, 255)
    RED_PALE  = (244,  66,  66)

    AQUA =  (53, 255, 211)
    AQUA2 = (71, 255, 218)

    GOLD   = (255, 231,  55)
    YELLOW = (255, 246, 130)

    # Salmon orange and/or red.
    SALMON  = (255, 111, 54)
    SALMON2 = (255, 117, 71)

    # ..................{ COLORMAPS                          }..................
    # Dictionary mapping from the name of each application-specific colormap to
    # be registered with matplotlib below to the two-dimensional sequence of all
    # colors defining this colormap. (See the add_colormap() function.)
    colormap_name_to_colors = {
        # Black-based colormaps.
        'betse_electric_cyan':    (BLACK, CYAN),
        'betse_electric_gold':    (BLACK, YELLOW),
        'betse_electric_green':   (BLACK, GREEN),
        'betse_electric_magenta': (BLACK, MAGENTA),
        'betse_electric_orange':  (BLACK, SALMON2),
        'betse_electric_blue':    (BLACK, BLUE_PALE),

        # Grey-based colormaps.
        'betse_blue_chalkboard':    (GREY_DARK, BLUE_PALE),
        'betse_cyan_chalkboard':    (GREY_DARK, CYAN),
        'betse_gold_chalkboard':    (GREY_DARK, GOLD),
        'betse_green_chalkboard':   (GREY_DARK, GREEN),
        'betse_magenta_chalkboard': (GREY_DARK, MAGENTA),
        'betse_orange_chalkboard':  (GREY_DARK, SALMON),

        # Salmon-based colormaps.
        'betse_alien_chalkboard': (SALMON2, GREY_DARK, AQUA2),
        'betse_alien_pale':       (SALMON2, BLACK, AQUA2),
        'betse_alien_solid':      (SALMON, BLACK, AQUA),

        # Purple-based colormaps.
        'betse_purple_green_chalkboard': (MAGENTA, GREY_DARK, GREEN),
        'betse_purple_green_pale':       (PURPLE_LIGHT, BLACK, GREEN_LIGHT),
        'betse_purple_green_solid':      (MAGENTA, BLACK, GREEN),

        # Blue-based colormaps.
        'betse_red_blue_chalkboard': (BLUE, GREY_DARK, RED),
        'betse_red_blue_pale':       (BLUE_PALE, BLACK, RED_PALE),
        'betse_red_blue_solid':      (BLUE, BLACK, RED),
    }

    # ..................{ REGISTRATION                       }..................
    # Register each such colormap with matplotlib.
    for colormap_name, colormap_name_colors in colormap_name_to_colors.items():
        add_colormap(colormap_name=colormap_name, colors=colormap_name_colors)
