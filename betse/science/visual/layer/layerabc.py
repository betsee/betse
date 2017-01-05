#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
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
import numpy as np
from abc import ABCMeta, abstractmethod, abstractproperty
from betse.lib.numpy import arrays
from betse.util.py import references
from betse.util.type.callables import property_cached
from betse.util.type.types import (
    type_check, IterableTypes, SequenceTypes, SequenceOrNoneTypes,)
from numpy import ndarray

# ....................{ SUPERCLASS                         }....................
class LayerCellsABC(object, metaclass=ABCMeta):
    '''
    Abstract base class of all classes spatially plotting a single feature of
    the cell cluster for a parent plot or animation.

    Each subclass of this class plots the spatial distribution of a single
    modelled variable (e.g., membrane voltage) for one or more simulation time
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
        self._visual = references.proxy_weak(visual)

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


    # Layer subclasses are recommended but *NOT* required to reimplement this
    # empty method. Since simplistic layers plotting artists whose appearance is
    # constant across all time steps only require the abstract _layer_first()
    # method to be implemented, this method is concrete rather than abstract.
    def _layer_next(self) -> None:
        '''
        Layer the spatial distribution of a single modelled variable (e.g., cell
        membrane voltage) for the next simulation time step onto the figure axes
        of the current plot or animation.
        '''

        pass

    # ..................{ SUBCLASS                           }..................
    # Subclasses are required to implement the following abstract methods.

    @abstractmethod
    def _layer_first(self) -> None:
        '''
        Layer the spatial distribution of a single modelled variable (e.g., cell
        membrane voltage) for the first simulation time step onto the figure
        axes of the current plot or animation.
        '''

        pass

# ....................{ SUBCLASSES                         }....................
class LayerCellsMappableABC(LayerCellsABC):
    '''
    Abstract base class of all classes spatially plotting a single cell-specific
    modelled variable (e.g., cell membrane voltage) of the cell cluster whose
    values are mappable as colors onto the colorbar for a parent plot or
    animation.

    Attributes
    ----------
    _color_mappables : IterableTypes
        Iterable of all mappables internally cached and returned by the
        :meth:`color_mappables` property.
    '''

    # ..................{ INITIALIZERS                       }..................
    def __init__(self) -> None:
        '''
        Initialize this layer.
        '''

        # Initialize our superclass.
        super().__init__()

        # Default all instance attributes.
        self._color_mappables = None

    # ..................{ PROPERTIES ~ read-only             }..................
    # Read-only properties, preventing callers from setting these attributes.

    @property
    def color_mappables(self) -> IterableTypes:
        '''
        Iterable of all **mappables** (i.e.,
        :class:`matplotlib.cm.ScalarMappable` instances) previously added by
        this layer to the figure axes of the current plot or animation to be
        associated with a figure colorbar.

        By Matplotlib design, only the first mappable in this iterable defines
        the color range for the figure colorbar; all other mappables are
        artificially constrained onto the same range.
        '''

        # If this iterable of mappables has yet to be cached...
        if self._color_mappables is None:
            # Layer the first time step, presumably caching this iterable.
            self._layer_first()

            # Assert this to be the case.
            assert self._color_mappables is not None, (
                '_layer_first() failed to define "self._color_mappables".')

        # Map the figure colorbar to this cached iterable.
        return self._color_mappables

    # ..................{ SUPERCLASS                         }..................
    def _layer_first(self) -> None:
        '''
        Layer the spatial distribution of a single cell-specific modelled
        variable (e.g., cell membrane voltage) for the first simulation time
        step onto the figure axes of the current plot or animation.

        This method internally caches the :attr:`_color_mappables` attribute
        returned by the :meth:`color_mappables` property.
        '''

        # Iterable of mappables layered by the subclass for the first time step.
        self._color_mappables = self._layer_first_color_mappables()

    # ..................{ SUBCLASS                           }..................
    # Subclasses are required to implement the following abstract methods.

    @abstractproperty
    def color_data(self) -> SequenceOrNoneTypes:
        '''
        Sequence of arbitrary dimensions of all possible color values for all
        time steps plotted by this layer _or_ `None` if calculating these values
        is impractical (e.g., due to space or time constraints).

        If colorbar autoscaling is:

        * Disabled, this sequence is ignored.
        * Enabled _and_ this sequence is:
          * `None`, this sequence is ignored. In this case, the subclass is
            responsible for colorbar autoscaling.
          * Non-`None`, the colorbar is clipped to the minimum and maximum
            scalar values unravelled from this sequence.
        '''

        pass


    @abstractmethod
    def _layer_first_color_mappables(self) -> IterableTypes:
        '''
        Layer the spatial distribution of a single cell-specific modelled
        variable (e.g., cell membrane voltage) for the first simulation time
        step onto the figure axes of the current plot or animation.

        Returns
        ----------
        IterableTypes
            Iterable of all mappables cached into the :attr:`_color_mappables`
            attribute by the :meth:`_layer_first` method.
        '''

        pass

# ....................{ SUBCLASSES                         }....................
class LayerCellsMappableArrayABC(LayerCellsMappableABC):
    '''
    Abstract base class of all classes spatially plotting a single cell-specific
    modelled variable (e.g., cell membrane voltage) of the cell cluster whose
    values are mappable as colors onto the colorbar of a parent plot or
    animation from a two-dimensional Numpy array of these values for all time
    steps to be animated.

    Attributes
    ----------
    _times_membranes_midpoint_data : ndarray
        Two-dimensional Numpy array of all arbitrary cell membrane data for all
        time steps to be animated cached by the :meth:`__init__` method and
        returned by the :meth:`times_membranes_midpoint_data` property.
    '''

    # ..................{ INITIALIZERS                       }..................
    def __init__(self, times_membranes_midpoint_data: SequenceTypes) -> None:
        '''
        Initialize this layer.

        Parameters
        ----------
        times_membranes_midpoint_data : Sequence
            Two-dimensional sequence of all cell membrane data for a single
            cell membrane-specific modelled variable (e.g., cell membrane
            voltage) for all simulation time steps, whose:
            . First dimension indexes each simulation time step.
            . Second dimension indexes each cell membrane in the simulated cell
              cluster, such that each element is arbitrary cell membrane data
              spatially situated at the midpoint of this membrane for this time
              step.
        '''

        # Initialize our superclass.
        super().__init__()

        # For efficiency, convert the passed sequence into a Numpy array.
        self._times_membranes_midpoint_data = arrays.from_sequence(
            times_membranes_midpoint_data)

    # ..................{ PROPERTIES ~ read-only             }..................
    # Read-only properties, preventing callers from setting these attributes.

    #FIXME: Refactor the following methods to internally defer to the
    #corresponding Cells.map_*() methods.

    @property_cached
    def times_cells_centre_data(self) -> ndarray:
        '''
        Two-dimensional Numpy array of all arbitrary cell data for all
        simulation time steps, whose:

        . First dimension indexes each simulation time step.
        . Second dimension indexes each cell in the simulated cell cluster, such
          that each element is arbitrary cell data spatially situated at the
          the centre of this cell for this time step.
        '''

        return self._visual._cells.map_membranes_midpoint_to_cells_centre_data(
            self._times_membranes_midpoint_data)


    @property
    def times_membranes_midpoint_data(self) -> ndarray:
        '''
        Two-dimensional sequence of all arbitrary cell membrane data for all
        simulation time steps, whose:

        . First dimension indexes each simulation time step.
        . Second dimension indexes each cell membrane in the simulated cell
          cluster, such that each element is arbitrary cell membrane data
          spatially situated at the midpoint of this membrane for this time
          step.
        '''

        return self._times_membranes_midpoint_data


    @property_cached
    def times_membranes_vertex_data(self) -> ndarray:
        '''
        Two-dimensional Numpy array of all arbitrary cell membrane vertex data
        for all simulation time steps, whose:

        . First dimension indexes each simulation time step.
        . Second dimension indexes each cell membrane vertex in the simulated
          cell cluster, such that each element is arbitrary data spatially
          situated at this cell membrane vertex for this time step.
        '''

        #FIXME: Consider generalizing this logic into a new public "Cells"
        #method ala the Cells.map_mems_midpoint_to_cells_centre_data()
        #method called above.

        return np.dot(
            self._times_membranes_midpoint_data,
            self._visual.cells.matrixMap2Verts)


    @property_cached
    def times_regions_centre_data(self) -> ndarray:
        '''
        Two-dimensional Numpy array of all arbitrary Voronoi region centre data
        for all simulation time steps, whose:

        . First dimension indexes each simulation time step.
        . Second dimension indexes each polygonal regions in the Voronoi
          diagram for the simulated cell cluster, such that each element is
          arbitrary data spatially situated at the centre of this region for
          this time step.
        '''

        #FIXME: Consider generalizing this logic into a new public "Cells"
        #method ala the Cells.map_mems_midpoint_to_cells_centre_data()
        #method called above.

        # Initialize this array of the desired shape with zeroes.
        times_regions_centre_data = np.zeros((
            len(self.times_cells_centre_data),
            len(self._visual.cells.voronoi_centres)))

        # Map cell- to region-centred data for all time steps.
        times_regions_centre_data[
            :, self._visual.cells.cell_to_grid] = (
            self.times_cells_centre_data)

        # Return this array for subsequent caching.
        return times_regions_centre_data
