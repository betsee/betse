#!/usr/bin/env python3
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Abstract base classes of all Matplotlib-based layer subclasses.
'''

#FIXME: The current approach to implementing animation overlays is
#fundamentally flawed. We currently attempt to provide a crude form of plot
#composition (i.e., merging two or more types of plots together into a single
#plot) by adding new booleans to the "AnimCellsABC" base class (e.g.,
#"is_current_overlayable") -- a fundamentally unwieldy and ultimately
#unworkable approach. By definition, you cannot provide true composability frow
#within a single class hierarchy. Instead, we need to split the specific
#process of generating different types of artists (e.g., mesh plots, stream
#plots) from the general process of animating and saving frames and plots as
#follows:
#
#* Refactor all concrete subclasses of "AnimCellsABC" into one or more
#  subclasses of "LayerCellsABC" instead, which may then be instantiated and
#  composed together into a new "layers" list passed to
#  CellsLayerABC.__init__(). For example:
#  * Split the existing "AnimGapJuncTimeSeries" subclass into:
#    * A new "CellsLayerGapJunc" subclass plotting *ONLY* the gap junction
#      open state as a "LineCollection" overlay. This layer subclass would
#      probably only be used for this specific purpose.
#    * A new "CellsLayerTimeSeries" subclass plotting *ONLY* an arbitrary
#      time series for the cell cluster as a mesh plot underlay. Clearly, this
#      layer subclass would be extensively reused elsewhere as well.
#* Replace all current overlay functionality in "AnimCellsABC" with "layers".
#* Refactor the configuration file from the current hard-coded non-composable
#  approach to a dynamic list-based approach permitting zero or more
#  user-defined animations, each consisting of one or more stock BETSE-defined
#  layers, to be defined. Users would then be able to construct arbitrarily
#  simple or complex animations as required.
#
#So, yes. It's quite a bit of work. But it's absolutely essential as well,
#particularly for implementing a general-purpose BETSE GUI.

# ....................{ IMPORTS                            }....................
import weakref
from abc import ABCMeta, abstractmethod
from betse.util.type.types import type_check

# ....................{ CLASSES                            }....................
class LayerCellsABC(object, metaclass=ABCMeta):
    '''
    Abstract base class of all classes spatially plotting a single feature of
    the cell cluster for a parent plot or animation.

    Each subclass of this class plots the spatial distribution of a single
    modelled variable (e.g., membrane voltage) for one on more simulation time
    steps. Each instance of the higher-level
    :class:`betse.science.visual.visualabc.VisualCellsABC` abstract base class
    contains one or more instances of subclasses of this lower-level class.

    Separating low-level layer logic from high-level plot and animation logic
    (e.g., multithreaded animation frame iteration, video and image exporting)
    enables composition between otherwise unrelated types. Thanks to plotters,
    two or more types of plots or animations may be trivially composed into a
    unique third type of plot or animation with _no_ modification to existing
    plotters, plots, or animations.

    Attributes
    ----------
    _is_layered : bool
        `True` only if the :meth:`layer` method has been called at least once
        for this layer instance.
    _visual : VisualCellsABC
        Plot or animation to layer onto.
    '''

    # ..................{ INITIALIZERS                       }..................
    def __init__(self) -> None:
        '''
        Initialize this layer.

        This method intentionally accepts _no_ parameters except constants
        parametrizing this layer's behaviour. In particular, this method
        accepts _no_ reference to the parent
        :class:`betse.science.visual.visualabc.VisualCellsABC` instance
        containing this layer instance _or_ to any other instances also
        contained by that parent instance (e.g., Matplotlib figure or axes
        objects). Why? Because plotters are instantiated by callers _before_
        their parent `VisualCellsABC` instances are instantiated.

        See Also
        ----------
        :meth:`layer`
            Further details on class design.
        '''

        # Default instance attributes.
        self._is_layered = False
        self._visual = None


    @type_check
    def prep(
        self, visual: 'betse.science.visual.visualabc.VisualCellsABC') -> None:
        '''
        Prepare this layer to be layered onto the passed plot or animation.

        Parameters
        ----------
        visual : VisualCellsABC
            Plot or animation to layer onto.
        '''

        # Classify this plot or animation with a weak rather than strong (the
        # default) reference, thereby avoiding circular references and the
        # resulting complications thereof (e.g., increased memory overhead).
        # Since the parent plot or animation necessarily lives significantly
        # longer than this layer, no complications arise. Ergo, this attribute
        # *ALWAYS* yields this object (rather than non-deterministically
        # yielding "None" if this object is unexpectedly garbage-collected).
        self._visual = weakref.proxy(visual)

    # ..................{ LAYERS                             }..................
    def layer(self) -> None:
        '''
        Layer the spatial distribution of a single modelled variable (e.g., cell
        membrane voltage) for the current time step and each cell of the current
        cluster onto the figure axes of the current plot or animation.
        '''

        # If this method has yet to be called...
        if not self._is_layered:
            # Perform logic specific to this call.
            self._layer_first()

            # Prevent subsequent calls to this method from repeating this logic.
            self._is_layered = True
        # Else, this method has been called at least once.
        else:
            # Perform logic specific to all subsequent calls.
            self._layer_next()


    @abstractmethod
    def _layer_first(self) -> None:
        '''
        Layer the spatial distribution of a single modelled variable (e.g., cell
        membrane voltage) for the first time step and each cell of the current
        cluster onto the figure axes of the current plot or animation.

        Layer subclasses are required to implement this abstract method.
        '''

        pass


    def _layer_next(self) -> None:
        '''
        Layer the spatial distribution of a single modelled variable (e.g., cell
        membrane voltage) for the next time step and each cell of the current
        cluster onto the figure axes of the current plot or animation.

        Layer subclasses are recommended but _not_ required to reimplement this
        empty method.
        '''

        pass
