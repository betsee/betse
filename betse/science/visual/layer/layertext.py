#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Layer subclasses printing text over the current cell cluster.
'''

# ....................{ IMPORTS                            }....................
from betse.science.visual import visuals
from betse.science.visual.layer.layerabc import LayerCellsABC

# ....................{ CLASSES                            }....................
class LayerCellsIndex(LayerCellsABC):
    '''
    Layer printing the 0-based index of each cell in the cell cluster as a text
    label centered on that cell.
    '''

    # ..................{ SUPERCLASS                         }..................
    def _layer_first(self) -> None:

        # For the 0-based index and 2-tuple of X and Y coordinates of the center
        # of each cell, display this index centered at these coordinates.
        for cell_index, cell_center in enumerate(
            self._visual.cells.cell_centres):
            self._visual.axes.text(
                # Text to be displayed.
                s=cell_index,

                # X and Y coordinates to display this text at.
                x=visuals.upscale_cell_coordinates(cell_center[0]),
                y=visuals.upscale_cell_coordinates(cell_center[1]),

                # Alignment of this text at these coordinates.
                horizontalalignment='center',
                verticalalignment='center',
            )
