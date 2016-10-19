#!/usr/bin/env python3
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Layer subclasses printing text over the current cell cluster.
'''

# ....................{ IMPORTS                            }....................
from betse.science.visual import visuals
from betse.science.visual.layer.layerabc import LayerCellsABC
from betse.util.type.types import type_check

# ....................{ BASE                               }....................
class LayerCellsIndex(LayerCellsABC):
    '''
    Layer subclass printing the 0-based index of each cell in the current cell
    cluster as a text label centered on that cell.

    Attributes
    ----------
    _is_layered : bool
        `True` only if the :meth:`layer` method has been called at least once
        for this layer instance.
    '''

    # ..................{ INITIALIZERS                       }..................
    def __init__(self) -> None:
        '''
        Initialize this layer.
        '''

        # Initialize our superclass.
        super().__init__()

        # Sanitize all instance attributes.
        self._is_layered = False

    # ..................{ SUPERCLASS                         }..................
    @type_check
    def layer(
        self, visual: 'betse.science.visual.visualabc.VisualCellsABC') -> None:
        '''
        Layer the 0-based index of each cell in the current cluster onto the
        figure axes of the passed plot or animation as a text label centered on
        that cell.

        Parameters
        ----------
        visual : VisualCellsABC
            Plot or animation to layer onto.
        '''

        # If this method has been called at least once (and hence already
        # layered cell indices onto this plot or animation), noop.
        if self._is_layered:
            return

        # Else, this method has yet to be called. Prevent subsequent calls to
        # this method from re-layering cell indices onto this plot or animation.
        self._is_layered = True

        # For the 0-based index and 2-tuple of X and Y coordinates of the center
        # of each cell, display this index centered at these coordinates.
        for cell_index, cell_center in enumerate(visual.cells.cell_centres):
            visual.axes.text(
                # Text to be displayed.
                s=cell_index,

                # X and Y coordinates to display this text at.
                x=visuals.upscale_cell_coordinates(cell_center[0]),
                y=visuals.upscale_cell_coordinates(cell_center[1]),

                # Alignment of this text at these coordinates.
                horizontalalignment='center',
                verticalalignment='center',
            )
