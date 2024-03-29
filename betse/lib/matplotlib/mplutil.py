#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2023 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level matplotlib-specific functionality for which no better home exists.
'''

# ....................{ IMPORTS                            }....................
from betse.util.io.error.errwarning import ignoring_warnings
from betse.util.type.types import GeneratorType
from matplotlib import MatplotlibDeprecationWarning

# ....................{ WARNINGS                           }....................
def ignoring_deprecations_mpl() -> GeneratorType:
    '''
    Single-shot context manager temporarily ignoring all matplotlib-specific
    deprecation warnings emitted by the :mod:`warnings` module for the duration
    of this context.

    See Also
    -----------
    :class:`ignoring_warnings`
        Further details.
    '''

    return ignoring_warnings(MatplotlibDeprecationWarning)
