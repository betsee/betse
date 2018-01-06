#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
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
