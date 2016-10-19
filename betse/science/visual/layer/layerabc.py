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
#* Add a new "plotters" parameter to the AnimCellsABC.__init__() constructor,
#  classified as a new "AnimCellsABC._plotters" instance variable. This
#  parameter and variable *MUST* be a list of "CellsPlotterABC" instances. The
#  order of plotters in this list defines the order in which these plotters are
#  drawn and hence overlaid onto one another (i.e., z-order).
#* Refactor AnimCellsABC.__init__() or a method called by that method to iterate
#  over "self._plotters" and initialize each such layer by calling
#  layer.init().
#* Refactor AnimCellsABC.plot_frame() or a related method to iterate over
#  "self._plotters" and draw each such layer by calling layer.draw().
#* Refactor all concrete subclasses of "AnimCellsABC" into one or more
#  subclasses of "CellsPlotterABC" instead, which may then be instantiated and
#  composed together into a new "plotters" list passed to
#  CellsPlotterABC.__init__(). For example:
#  * Split the existing "AnimGapJuncTimeSeries" subclass into:
#    * A new "CellsPlotterGapJunc" subclass plotting *ONLY* the gap junction
#      open state as a "LineCollection" overlay. This layer subclass would
#      probably only be used for this specific purpose.
#    * A new "CellsPlotterTimeSeries" subclass plotting *ONLY* an arbitrary
#      time series for the cell cluster as a mesh plot underlay. Clearly, this
#      layer subclass would be extensively reused elsewhere as well.
#* Replace all current overlay functionality in "AnimCellsABC" with "plotters".
#* Refactor the configuration file from the current hard-coded non-composable
#  approach to a dynamic list-based approach permitting zero or more
#  user-defined animations, each consisting of one or more stock BETSE-defined
#  plotters, to be defined. Users would then be able to construct arbitrarily
#  simple or complex animations as required.
#
#So, yes. It's quite a bit of work. But it's absolutely essential as well,
#particularly for implementing a general-purpose BETSE GUI.

# ....................{ IMPORTS                            }....................
from abc import ABCMeta, abstractmethod
# from betse.util.type.types import type_check

# ....................{ BASE                               }....................
class LayerCellsABC(object, metaclass=ABCMeta):
    '''
    Abstract base class of all classes spatially plotting a single feature of
    the cell cluster.

    Subclasses of this class plot the spatial distribution of a single modelled
    variable (e.g., membrane voltage) for one on more simulation time steps.

    Instances of these subclasses are contained by, and hence lower-level than,
    instances of the higher-level
    :class:`betse.science.visual.visualabc.VisualCellsABC` abstract base class.
    Architecturally speaking, each instance of that abstract base class
    contains one or more instances of subclasses of this abstract base class.
    Human-readably speaking, each high-level plot and animation object contains
    multiple low-level layer objects implementing the drawing of that plot or
    animation object.

    Separating low-level layer logic from high-level plot and animation logic
    (e.g., multithreaded animation frame iteration, video and image exporting)
    enables composition between otherwise unrelated types. Thanks to plotters,
    two or more types of plots or animations may be trivially composed into a
    unique third type of plot or animation with _no_ modification to existing
    plotters, plots, or animations.
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

        pass

    # ..................{ ABSTRACT                           }..................
    @abstractmethod
    def layer(
        self, visual: 'betse.science.visual.visualabc.VisualCellsABC') -> None:
        '''
        Layer the spatial distribution of a single modelled variable (e.g., cell
        membrane voltage) for the current time step and each cell of the current
        cluster onto the figure axes of the passed plot or animation.

        Parameters
        ----------
        plot : VisualCellsABC
            Plot or animation to layer onto, passed to this rather than the
            :meth:`__init__` method to avoid chicken-and-egg issues and hence:
            * Avoid long-lived circular references between layer, plot, and
              animation instances and the resulting memory costs.
            * Permit callers to:
              * Create layer instances _before_ plot or animation instances.
              * Cache previously created layer instances.
              * Share previously cached layer instances between two or more
                plot or animation instances, reducing memory footprint.
              * Pass previously cached layer instances to the :meth:`__init__`
                methods of plot or animation subclasses.
        '''

        pass
