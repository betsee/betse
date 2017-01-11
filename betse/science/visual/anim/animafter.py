#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
**Post-simulation animation** (i.e., animation produced *after* rather than
*while* solving a simulation) subclasses.
'''

# ....................{ IMPORTS                            }....................
from betse.science.visual.anim.animabc import AnimCellsABC
from betse.util.type.types import type_check, SequenceTypes

# ....................{ SUBCLASSES                         }....................
class AnimCellsAfterSolving(AnimCellsABC):
    '''
    Abstract base class of all post-simulation animation subclasses.

    Each subclass of this class animates arbitrary cell data as a time series
    plotted over the cell cluster (e.g., cell membrane voltage as a function of
    time) *after* rather than *while* solving simulations.
    '''

    @type_check
    def __init__(
        self,
        p: 'betse.science.parameters.Parameters',
        *args, **kwargs
    ) -> None:
        '''
        Initialize this post-simulation animation.

        Parameters
        ----------
        p : Parameters
            Current simulation configuration.

        See the superclass `__init__()` method for all remaining parameters.
        '''

        # Initialize our superclass.
        super().__init__(
            # Pass this simulation configuration as is to our superclass.
            p=p,

            # Save and show this post-simulation animation only if this
            # configuration enables doing so.
            is_save=p.anim.is_after_sim_save,
            is_show=p.anim.is_after_sim_show,

            # Save all post-simulation animations to the same parent directory.
            save_dir_parent_basename='anim',

            # Pass all remaining arguments as is to our superclass.
            *args, **kwargs
        )


#FIXME: Merge this class into the "AnimCellsAfterSolving" superclass *AFTER*
#refactoring all subclasses to leverage layers.
class AnimCellsAfterSolvingLayered(AnimCellsAfterSolving):
    '''
    Post-simulation animation of arbitrary cell data as a time series (e.g.,
    cell membrane voltage as a function of time), plotted over the cell
    cluster using a predefined sequence of layers.
    '''

    # ..................{ SUPERCLASS                         }..................
    @type_check
    def __init__(self, layers: SequenceTypes, *args, **kwargs) -> None:
        '''
        Initialize this animation.

        Parameters
        ----------
        layers : SequenceTypes
            Sequence of all :class:`LayerCellsABC` instances collectively
            plotting each frame of this animation. See the
            :meth:`VisualCellsABC.__init__` method for further details.

        All remaining parameters are passed as is to the superclass
        :meth:`AnimCellsAfterSolving.__init__` method.
        '''

        # Initialize the superclass.
        super().__init__(layers=layers, *args, **kwargs)

        # Display and/or save this animation.
        self._animate()

# ....................{ SUBCLASSES ~ obsolete              }....................
#FIXME: Replace use of all the following subclasses with the
#"AnimCellsAfterSolvingLayered" subclass *AFTER* refactoring all subclasses to
#leverage layers.
class AnimField(AnimCellsAfterSolving):
    '''
    Abstract base class of all animations of electric field strength plotted on
    the current cell cluster.

    Attributes
    ----------
    _magnitude_time_series : SequenceTypes
        Electric field magnitudes as a function of time.
    _mesh_plot : matplotlib.image.AxesImage
        Meshplot of the current or prior frame's electric field magnitude.
    _stream_plot : matplotlib.streamplot.StreamplotSet
        Streamplot of the current or prior frame's electric field.
    _x_time_series : SequenceTypes
        Electric field X components as a function of time.
    _y_time_series : SequenceTypes
        Electric field Y components as a function of time.
    _unit_x_time_series : SequenceTypes
        Electric field X unit components as a function of time. The resulting
        electric field vectors are **unit vectors** (i.e., have magnitude 1).
    _unit_y_time_series : SequenceTypes
        Electric field Y unit components as a function of time. The resulting
        electric field vectors are **unit vectors** (i.e., have magnitude 1).
    '''

    @type_check
    def __init__(
        self,
        x_time_series: SequenceTypes,
        y_time_series: SequenceTypes,
        *args, **kwargs
    ) -> None:
        '''
        Initialize this animation.

        Parameters
        ----------
        x_time_series : SequenceTypes
            SequenceTypes (e.g., list, numpy array) of all electric field
            strength X components indexed by simulation time.
        y_time_series : SequenceTypes
            SequenceTypes (e.g., list, numpy array) of all electric field
            strength Y components indexed by simulation time.

        See the superclass `__init__()` method for all remaining parameters.
        '''

        # Pass all parameters *NOT* listed above to our superclass.
        super().__init__(
            # Since this class already plots a streamplot, prevent the
            # superclass from plotting another streamplot as an overlay.
            is_current_overlayable=False,
            *args, **kwargs)

        # Classify all remaining parameters.
        self._x_time_series = x_time_series
        self._y_time_series = y_time_series

        # Electric field magnitudes and X and Y unit components.
        self._magnitude_time_series = []
        self._unit_x_time_series = []
        self._unit_y_time_series = []

        # Prefer an alternative colormap.
        self._colormap = self._p.background_cm


class AnimVelocity(AnimCellsAfterSolving):
    '''
    Abstract base class of all animations of a velocity flow plotted on the
    current cell cluster.
    '''

    def __init__(self, *args, **kwargs) -> None:

        # Pass all parameters *NOT* listed above to our superclass.
        super().__init__(
            # Since this class already plots a streamplot, prevent the
            # superclass from plotting another streamplot as an overlay.
            is_current_overlayable=False,
            *args, **kwargs)

        # Prefer an alternative colormap.
        self._colormap = self._p.background_cm
