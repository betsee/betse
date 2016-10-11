#!/usr/bin/env python3
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Abstract base classes of all Matplotlib-based plotter classes.
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
#* Define a new "betse.science.plot.anim.plotter" submodule.
#* Define a new "CellsPlotterABC" abstract base class in this submodule. Plotter
#  classes encapsulate the plotting of a single type of plot (e.g., Vmem-style
#  cell cluster meshplot, electrical current streamplot). While they may
#  internally cache one or more time series required for this plotting, they do
#  *NOT* contain the axes, figure, plot, or animation currently being plotted.
#  Since these Matplotlib artists are shared between multiple plotters, the
#  existing "PlotCellsABC" and "AnimCellsABC" base classes retain ownership of these
#  artists.
#* Define the following abstract methods in "CellsPlotterABC":
#  * An __init__() method accepting *NO* parameters except constants
#    parametrizing this plotter's behaviour. In particular, this method must
#    *NOT* accept the current "AnimABC" instance or Matplotlib figure or axes
#    objects. Why? Because plotters are instantiated by external callers rather
#    than that "AnimABC" instance *BEFORE* that "AnimABC" instance is
#    instantiated.
#  * A prep() method accepting as parameters all internal attributes
#    subsequently required by the plot() method. Subclasses *MUST* redefine
#    this method to initialize this plotter instance (e.g., by internally
#    caching one or more time series). This method is internally called by the
#    current "AnimABC" instance.
#  * A plot() method presumably accepting the index of the frame to be plotted.
#    Subclasses *MUST* redefine this method to plot the time series data for
#    this frame into the current plot or animation figure. To avoid circular
#    references (and hence unallocated figures), it wouldn't necessarily be a
#    bad idea to repass all objects required for plotting to each
#    CellsPlotterABC.plot() call rather than to the CellsPlotterABC.__init__()
#    constructor. Food for pleasant thought, anyway.
#  * *UHM.* Wait. There would clearly be *NO* circular references here, so
#    passing objects to CellsPlotterABC.__init__() once rather than repeatedly
#    to CellsPlotterABC.plot() probably makes the most sense. We absolutely
#    *MUST*, however, avoid passing the current "AnimCellsABC" instance to
#    CellsPlotterABC.__init__(). Passing that instance to
#    CellsPlotterABC.plot(), however, should both be safe and desirable, as
#    doing so would then permit plotters to access relevant data on the current
#    animation (e.g., "AnimCellsABC.time_step" providing the current frame number).
#* Add a new "plotters" parameter to the AnimCellsABC.__init__() constructor,
#  classified as a new "AnimCellsABC._plotters" instance variable. This parameter
#  and variable *MUST* be a list of "CellsPlotterABC" instances. The order of
#  plotters in this list defines the order in which these plotters are drawn
#  and hence overlaid onto one another (i.e., z-order).
#* Refactor AnimCellsABC.__init__() or a method called by that method to iterate
#  over "self._plotters" and initialize each such plotter by calling
#  plotter.init().
#* Refactor AnimCellsABC.plot_frame() or a related method to iterate over
#  "self._plotters" and draw each such plotter by calling plotter.draw().
#* Refactor all subclasses of "AnimCellsABC" into one or more subclasses of
#  "CellsPlotterABC" instead, which may then be instantiated and composed
#  together into a new "plotters" list passed to CellsPlotterABC.__init__(). For
#  example:
#  * Split the existing "AnimGapJuncTimeSeries" subclass into:
#    * A new "CellsPlotterGapJunc" subclass plotting *ONLY* the gap junction
#      open state as a "LineCollection" overlay. This plotter subclass would
#      probably only be used for this specific purpose.
#    * A new "CellsPlotterTimeSeries" subclass plotting *ONLY* an arbitrary
#      time series for the cell cluster as a mesh plot underlay. Clearly, this
#      plotter subclass would be extensively reused elsewhere as well.
#* Replace all current overlay functionality in "AnimCellsABC" with "plotters".
#* Refactor the configuration file from the current hard-coded non-composable
#  approach to a dynamic list-based approach permitting zero or more
#  user-defined animations, each consisting of one or more stock BETSE-defined
#  plotters, to be defined. Users would then be able to construct arbitrarily
#  simple or complex animations as required.
#
#Note that the "CellsPlotterABC" nomenclature used above is overly ambiguous
#and hence non-ideal. Non-ambiguous alternative names for this concept include:
#
#* "PlottableABC". Yes, we quite appreciate this one. Such objects are
#  certainly plottable, so this is coherent.
#* "ComposableABC".
#* "DrawableABC".
#
#So, yes. It's quite a bit of work. But it's absolutely essential as well,
#particularly for implementing a general-purpose BETSE GUI.

# ....................{ IMPORTS                            }....................
# import numpy as np
import weakref
from abc import ABCMeta  #, abstractmethod  #, abstractstaticmethod
# from betse.exceptions import BetseMethodException
# from betse.lib.matplotlib.matplotlibs import ZORDER_STREAM
# # from betse.util.io.log import logs
# from betse.util.type import types, objects
# from betse.util.type.types import (
#     type_check, NoneType, NumericTypes, SequenceTypes)
# from matplotlib import pyplot
# from matplotlib.collections import PolyCollection
# from matplotlib.colors import Colormap
# from matplotlib.patches import FancyArrowPatch

# ....................{ BASE                               }....................
class PlotterCellsABC(object, metaclass=ABCMeta):
    '''
    Abstract base class of all classes spatially plotting a single feature of
    the cell cluster.

    Subclasses of this class plot the spatial distribution of a single modelled
    variable (e.g., membrane voltage) for one on more simulation time steps.

    Instances of these subclasses are contained by, and hence lower-level than,
    instances of the higher-level
    :class:`betse.science.plot.plotabc.PlotCellsABC` abstract base class.
    Architecturally speaking, each instance of that abstract base class
    contains one or more instances of subclasses of this abstract base class.
    Human-readably speaking, each high-level plot and animation object contains
    multiple low-level plotter objects implementing the drawing of that plot or
    animation object.

    Separating low-level plotter logic from high-level plot and animation logic
    (e.g., multithreaded animation frame iteration, video and image exporting)
    enables composition between otherwise unrelated types. Thanks to plotters,
    two or more types of plots or animations may be trivially composed into a
    unique third type of plot or animation with _no_ modification to existing
    plotters, plots, or animations.

    Attributes
    ----------
    _plot : PlotCellsABC
        Parent plot or animation object uniquely containing this plotter object.
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(
        self,

        #FIXME: Ideally, we would type check the first three parameters like so:
        #
        #    sim: Simulator,
        #    cells: Cells,
        #    p: Parameters,
        #
        #Sadly, doing so introduces circular import issues. To circumvent this,
        #@type_check should be improved to support annotations as strings. When
        #an annotation is a string, this decorator should attempt to dynamically
        #import this string as the fully-qualified name of a class *AND* then
        #use this class to type check this parameter: e.g.,
        #
        #    sim: 'betse.science.sim.Simulator',
        #    cells: 'betse.science.cells.Cells',
        #    p: 'betse.science.parameters.Parameters',

    ) -> None:
        '''
        Initialize this plotter.

        Parameters
        ----------
        '''

        # Classify core parameters with weak rather than strong (the default)
        # references, thus avoiding circular references and the resulting
        # complications thereof (e.g., increased memory overhead). Since these
        # higher-level objects uniquely own this lower-level object, no
        # complications arise. Ergo, these attributes *ALWAYS* yield these
        # objects rather than non-deterministically returning "None" if these
        # objects are unexpectedly garbage-collected.
        # self.sim = weakref.proxy(sim)
        pass
