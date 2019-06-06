#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level matplotlib-specific functionality for which no better home exists.
'''

# ....................{ IMPORTS                           }....................
from betse.util.io.error.errwarning import ignoring_warnings
from betse.util.type.types import GeneratorType
from matplotlib.cbook import MatplotlibDeprecationWarning

# ....................{ WARNINGS                          }....................
#FIXME: Rename to ignoring_deprecations_mpl() for disambiguity.
def deprecations_ignored() -> GeneratorType:
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
