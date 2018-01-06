#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Matplotlib-specific **z-order** (i.e., positive integers ordering artist
drawing, such that artists with larger z-orders are drawn over artists with
smaller z-orders) facilities.

See Also
----------
http://matplotlib.org/examples/pylab_examples/zorder_demo.html
    Canonical z-order example from which most constants defined by this
    submodule were derived.
'''

# ....................{ CONSTANTS                          }....................
ZORDER_PATCH = 1
'''
Default **z-order** (i.e., positive integer ordering artist drawing, such that
artists with larger z-orders are drawn over artists with smaller z-orders) for
patch artists (e.g., :class:`Patch`, :class:`PatchCollection`).

This is the lowest default z-order, thus drawing patch artists under all other
artists by default.
'''


ZORDER_LINE = 2
'''
Default **z-order** (i.e., positive integer ordering artist drawing, such that
artists with larger z-orders are drawn over artists with smaller z-orders) for
line artists (e.g., `Line2D`, `LineCollection`, `StreamplotSet`).

This is the middle default z-order, thus drawing line artists over all patch
artists but under all text artists by default.
'''


ZORDER_TEXT = 3
'''
Default **z-order** (i.e., positive integer ordering artist drawing, such that
artists with larger z-orders are drawn over artists with smaller z-orders) for
text artists (e.g., `Text`).

This is the highest default z-order, thus drawing text artists over all other
artists by default.
'''


ZORDER_STREAM = (ZORDER_LINE + ZORDER_TEXT) / 2
'''
BETSE-specific **z-order** (i.e., positive integer ordering artist drawing,
such that artists with larger z-orders are drawn over artists with smaller
z-orders) for streamplots (e.g., `StreamplotSet`).

This magic number has been chosen such that streamplots with this z-order will
be drawn over all line and patch artists but under all text artists by default.
Streamplots are technically a variant of line artists but sufficiently non-
linear (and visually busy) to warrant separate handling.

It may also be pertinent to note that recent Matplotlib releases as of this
writing (1.40) accidentally broke backward compatibility with respect to
default streamplot z-order. Specifically, streamplots were assigned a default
z-order of `ZORDER_PATCH` rather than `ZORDER_LINE`:

    https://github.com/matplotlib/matplotlib/pull/5567
'''
