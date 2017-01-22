#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level matplotlib-specific functionality for which no better home exists.
'''

# ....................{ IMPORTS                            }....................
from betse.exceptions import BetseMatplotlibException
from betse.util.io.warning import warnings_ignored
from betse.util.type.types import type_check, GeneratorType
from matplotlib import cm as colormaps
from matplotlib.cbook import MatplotlibDeprecationWarning
from matplotlib.colors import Colormap

# ....................{ WARNINGS                           }....................
def deprecations_ignored() -> GeneratorType:
    '''
    Single-shot context manager temporarily ignoring all matplotlib-specific
    deprecation warnings emitted by the :mod:`warnings` module for the duration
    of this context.

    See Also
    -----------
    :class:`warnings_ignored`
        Further details.
    '''

    return warnings_ignored(warning_class=MatplotlibDeprecationWarning)

# ....................{ GETTERS                            }....................
@type_check
def get_colormap(colormap_name: str) -> Colormap:
    '''
    Matplotlib colormap with the passed name.

    Parameters
    ----------
    colormap_name : str
        Name of the attribute in the :mod:`matplotlib.cm` module corresponding
        to the desired colormap (e.g., `Blues`, `Greens`, `jet`, `rainbow).

    Returns
    ----------
    Colormap
        Matplotlib colormap with this name.

    See Also
    ----------
    http://matplotlib.org/examples/color/colormaps_reference.html
        List of supported colormaps.
    '''

    # Colormap with the passed name if any or "None" otherwise.
    colormap = getattr(colormaps, colormap_name, None)

    # If no such colprmap exists, raise an exception.
    if not isinstance(colormap, Colormap):
        raise BetseMatplotlibException(
            'Matplotlib colormap "{}" not found.'.format(colormap_name))

    # Return this colormap.
    return colormap
