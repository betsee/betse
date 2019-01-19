#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level :mod:`matplotlib`-specific figure functionality.
'''

# ....................{ IMPORTS                           }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: To allow matplotlib defaults (e.g., backend, logging) to be replaced
# with application-specific preferences, the unsafe "matplotlib.pyplot"
# submodule must *NOT* be imported here at the top-level.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

from betse.exceptions import BetseMatplotlibException
from betse.util.io.log import logs
from betse.util.type.types import type_check, MatplotlibFigureType

# ....................{ EXCEPTIONS                        }....................
def die_unless_figure() -> bool:
    '''
    Raise an exception
    ``True`` only if one or more figures are currently open with the
    :mod:`matplotlib.pyplot` GCF API.

    Raises
    -----------
    BetseMatplotlibException
        If no figures are currently open with this API.

    See Also
    -----------
    :func:`is_figure`
        Further details.
    '''

    if not is_figure():
        raise BetseMatplotlibException('No matplotlib figures currently open.')

# ....................{ TESTERS                           }....................
def is_figure() -> bool:
    '''
    ``True`` only if one or more figures are currently open with the
    :mod:`matplotlib.pyplot` GCF API.
    '''

    # Defer heavyweight imports.
    from matplotlib import pyplot

    # Note that, despite the name, the pyplot.get_fignums() function returns a
    # list of the numbers of currently open figures rather than the number of
    # currently open figures. Welcome to Matplotlib Hell, where all is well.
    return bool(pyplot.get_fignums())

# ....................{ GETTERS                           }....................
def get_figure_current() -> MatplotlibFigureType:
    '''
    Figure most recently opened with the :mod:`matplotlib.pyplot` GCF API.

    Specifically, this function returns the same figure instance returned by
    the most recent call to the :func:`matplotlib.pyplot.figure` function.
    '''

    # Defer heavyweight imports.
    from matplotlib import pyplot

    # If no figure is currently open, raise an exception.
    die_unless_figure()

    # Three-letter acronyms do *NOT* constitute a human-readable API.
    return pyplot.gcf()

# ....................{ CLOSERS                           }....................
@type_check
def close_figure(figure: MatplotlibFigureType) -> None:
    '''
    **Close** (i.e., clear, delete, remove, garbage collect) the passed figure,
    guaranteeing that all resources consumed by this figure will be
    subsequently reclaimed after an indeterminate period of time.

    Specifically, this function (in order):

    #. Nullifies the contents of all non-axes artists of this figure.
    #. Nullifies the contents of all axes of this figure.
    #. Closes all interactive windows (if any) associated with this figure.
    #. Deassociates this figure from the :mod:`matplotlib.pyplot` GCF API.

    Caveats
    -----------
    **Figure closure is surprisingly non-trivial.** Failure to call this
    function on open figures usually results in some or all of the contents of
    those figures continuing to reside in memory. Ergo, this high-level
    function should *always* be called in lieu of low-level
    :mod:`matplotlib.pyplot` functions (e.g., :mod:`matplotlib.pyplot.close`)
    or equally low-level :class:`matplotlib.figure.Figure` methods (e.g.,
    :meth:`matplotlib.figure.Figure.cla`), which are all known to behave
    non-deterministically and hence unsafely.

    Parameters
    -----------
    figure : MatplotlibFigureType
        Figure to be closed.

    See Also
    -----------
    https://stackoverflow.com/a/17106460/2809027
        StackOverflow answer strongly inspiring this implementation.
    '''

    # Defer heavyweight imports.
    from matplotlib import pyplot

    # Log this closure.
    logs.log_debug('Closing matplotlib figure "%r"...', figure.number)

    # Nullify the contents of all non-axes artists of this figure.
    figure.clf()

    # For each axis of this figure, nullify the contents of this axis *AFTER*
    # that of non-axes. Why? We have utterly no idea. The matplotlib API is
    # dark voodoo that, frankly, not even matplotlib developers grok. It is
    # unsafe to ask too many questions.
    for figure_axis in figure.axes:
        figure_axis.cla()

    # Deassociate this figure from the "matplotlib.pyplot" GCF API, implicitly
    # closing all interactive windows (if any) associated with this figure. In
    # theory, doing so should also implicitly nullify the above content; in
    # practice, however, the matplotlib API is sufficiently unreliable across
    # version bumps that explicitly nullifying that content is all but
    # non-optional. Thanks for all of the bugs, matplotlib!
    pyplot.close(figure)


def close_figure_current() -> None:
    '''
    **Close** (i.e., clear, delete, remove, garbage collect) the figure most
    recently opened with the :mod:`matplotlib.pyplot` GCF API, guaranteeing
    that all resources consumed by this figure will be subsequently reclaimed
    after an indeterminate period of time.

    Specifically, this function closes the same figure instance returned by
    the most recent call to the :func:`matplotlib.pyplot.figure` function.

    See Also
    -----------
    :func:`close_figure`
        Further details.
    '''

    # Stay classy, matplotlib.
    close_figure(get_figure_current())


def close_figures_all() -> None:
    '''
    **Close** (i.e., clear, delete, remove, garbage collect) all figures
    currently opened with the :mod:`matplotlib.pyplot` GCF API, guaranteeing
    that all resources consumed by these figures will be subsequently reclaimed
    after an indeterminate period of time.

    See Also
    -----------
    https://stackoverflow.com/a/13174720/2809027
        StackOverflow answer strongly inspiring this implementation.
    :func:`close_figure`
        Further details.
    '''

    # Defer heavyweight imports.
    from matplotlib import pyplot

    # Log this closure.
    logs.log_debug('Closing all matplotlib figures...')

    # Mortifying API, thy name is the "pyplot" interface.
    pyplot.close('all')
