#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Matplotlib-specific utility functions.
'''

# ....................{ IMPORTS                            }....................
from matplotlib import cm as colormaps
from matplotlib.colors import Colormap

from betse.exceptions import BetseExceptionParameters
from betse.util.type import types

# ....................{ GETTERS                            }....................
def get_colormap(colormap_name: str) -> Colormap:
    '''
    Get the Matplotlib colormap with the passed name.

    Parameters
    ----------
    colormap_name : str
        Name of the attribute in the `matplotlib.cm` module corresponding to the
        desired colormap (e.g., `Blues`, `Greens`, `jet`, `rainbow).

    See Also
    ----------
    http://matplotlib.org/examples/color/colormaps_reference.html
        List of supported colormaps.
    '''
    assert types.is_str(colormap_name), types.assert_not_str(colormap_name)

    colormap = getattr(colormaps, colormap_name, None)
    if not isinstance(colormap, Colormap):
        raise BetseExceptionParameters(
            'Matplotlib colormap "{}" not found.'.format(colormap_name))
    return colormap