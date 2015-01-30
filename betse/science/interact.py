#!/usr/bin/env python3
# Copyright 2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

import matplotlib.pyplot as plt


def get_inds(event):
    """
    Returns the indices of a plot object (patch, line, point)
    for a mouse click event on that object and changes
    the alpha for that plot object.


    """
    index = event.ind
    patch = event.artist
    patch.set_alpha(0.1)
    print(*index)
    fig = plt.gcf()
    ax = plt.gca()
    fig.canvas.draw()
    return index
