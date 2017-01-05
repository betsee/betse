#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Abstract base classes of all Matplotlib-based plot and animation subclasses.
'''

#FIXME: Refactor all procedural cell cluster-specific
#"betse.science.visual.plot.plotutil" functions into subclasses of the
#"LayerCellsABC" base class defined elsewhere.
#
#Ultimate power fights the dark deceit!

# ....................{ IMPORTS                            }....................
import numpy as np
from abc import ABCMeta  #, abstractmethod  #, abstractstaticmethod
from betse.exceptions import BetseMethodException
from betse.lib.matplotlib.mplzorder import ZORDER_STREAM
from betse.lib.numpy import arrays
from betse.science.visual import visuals
from betse.science.visual.layer.layerabc import (
    LayerCellsABC, LayerCellsMappableABC)
from betse.science.visual.layer.layertext import LayerCellsIndex
from betse.util.io.log import logs
from betse.util.py import references
from betse.util.type import iterables, types
from betse.util.type.iterables import SENTINEL
from betse.util.type.obj import objs
from betse.util.type.types import (
    type_check,
    IterableOrNoneTypes,
    NoneType,
    NumericTypes, NumericOrNoneTypes,
    SequenceTypes, SequenceOrNoneTypes,
)
from matplotlib import pyplot
from matplotlib.axes import Axes
from matplotlib.cm import ScalarMappable
from matplotlib.collections import PolyCollection
from matplotlib.colors import Colormap
from matplotlib.image import AxesImage
from matplotlib.patches import FancyArrowPatch
from matplotlib.streamplot import StreamplotSet

# ....................{ BASE                               }....................
class VisualCellsABC(object, metaclass=ABCMeta):
    '''
    Abstract base class of all subclasses spatially plotting or animating the
    currently simulated cell cluster.

    Subclasses of this class plot the spatial distribution of one or more
    modelled variables (e.g., charge distribution, membrane voltage) for
    either:

    * A single simulation time step, in which this subclass plots a still frame
      of these variables for this step.
    * One or more simulation time steps, in which this subclass animates a
      video of these variables for these steps.

    Attributes (Private)
    ----------
    _cells : Cells
        Current cell cluster.
    _p : Parameters
        Current simulation configuration.
    _sim : Simulator
        Current simulation.
    _label : str
        Basename of the subdirectory in the phase-specific results directory
        to which all files exported for this plot or animation are saved _and_
        the basename prefix of these files.
    _layers : list
        List of all :class:`LayerCellsABC` instances collectively composing
        this plot or animation.

    Attributes (Private: Figure)
    ----------
    _figure : Figure
        Matplotlib figure providing the current animation frame.
    _figure_title : str
        Text displayed above the figure itself.

    Attributes (Private: Axes)
    ----------
    _axes : Axes
        Matplotlib figure axes providing the current animation frame data.
    _axes_bounds : list
        Spacial extent of the current 2D environment as a 4-element list
        conisting of (in order):
        1. The minimum value of the figure's X axis.
        2. The maximum value of the figure's X axis.
        3. The minimum value of the figure's Y axis.
        4. The maximum value of the figure's Y axis.
    _axes_title : str
        Text displayed above the figure axes. If a non-`None` value for the
        `axes_title` parameter is passed to the `__init__()` method, this is
        that value; else, this is the value of the `figure_title` parameter
        passed to the same method.
    _axes_x_label : str
        Text displayed below the figure's X axis.
    _axes_y_label : str
        Text displayed to the left of the figure's Y axis.

    Attributes (Private: Color)
    ----------
    _color_mappables : IterableTypes
        Iterable of all plotted mappables (i.e., instances of the
        :class:`ScalarMappable` class added to the figure axes of this plot or
        animation) whose minimum and maximum values define the color range
        displayed both in the figure axes and colorbar of this plot or
        animation. Due to Matplotlib constraints, only the first mappable in
        this iterable is associated with this colorbar.
    _colorbar_title: str
        Text displayed above the figure colorbar.
    _color_min : float
        Minimum color value to be displayed by the colorbar. If colorbar
        autoscaling is enabled (i.e., `_is_color_autoscaled` is `True`), the
        subclass is responsible for redefining this value as appropriate.
    _color_max : float
        Maximum color value to be displayed by the colorbar. If colorbar
        autoscaling is enabled (i.e., `_is_color_autoscaled` is `True`), the
        subclass is responsible for redefining this value as appropriate.
    _colormap : Colormap
        Matplotlib colormap with which to create this animation's colorbar.
    _is_color_autoscaled : bool
        `True` if dynamically resetting the minimum and maximum colorbar values
        to be the corresponding minimum and maximum values for the current
        frame _or_ `False` if statically setting the minimum and maximum
        colorbar values to predetermined constants.
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(
        self,

        # Mandatory parameters.
        sim: 'betse.science.sim.Simulator',
        cells: 'betse.science.cells.Cells',
        p: 'betse.science.parameters.Parameters',
        is_save: bool,
        is_show: bool,
        label: str,
        figure_title: str,
        colorbar_title: str,
        is_color_autoscaled: bool,
        color_min: NumericTypes,
        color_max: NumericTypes,

        # Optional parameters.
        axes_title: (str, NoneType) = None,
        axes_x_label: str = 'Spatial Distance [um]',
        axes_y_label: str = 'Spatial Distance [um]',
        colormap: (Colormap, NoneType) = None,

        #FIXME: Refactor this to be mandatory instead.
        layers: SequenceOrNoneTypes = None,
    ) -> None:
        '''
        Initialize this plot or animation.

        Parameters
        ----------
        sim : Simulator
            Current simulation.
        cells : Cells
            Current cell cluster.
        p : Parameters
            Current simulation configuration.
        is_save : bool
            `True` only if non-interactively saving this plot or animation.
        is_show : bool
            `True` only if interactively displaying this plot or animation.
        label : str
            Terse machine-readable string (e.g., `Vmem`) serving as both:
            * The basename of the subdirectory of the phase-specific results
              directory containing all files saved by this plot or animation.
            * The basename prefix of these files.
        figure_title : str
            Text displayed above the figure itself.
        colorbar_title: str
            Text displayed above the figure colorbar.
        is_color_autoscaled : bool
            `True` if dynamically resetting the minimum and maximum colorbar
            values to be the corresponding minimum and maximum values for the
            current frame _or_ `False` if statically setting the minimum and
            maximum colorbar values to predetermined constants.
        color_min : NumericTypes
            Minimum colorbar value to be used if `is_color_autoscaled` is
            `False`. If `is_color_autoscaled` is `True`, this parameter is
            silently ignored.
        color_max : NumericTypes
            Maximum colorbar value to be used if `is_color_autoscaled` is
            `False`. If `is_color_autoscaled` is `True`, this parameter is
            silently ignored.
        axes_title : optional[str]
            Text displayed above the figure axes but below the figure title.
            Defaults to `None`, in which case no such text is displayed.
        axes_x_label : optional[str]
            Text displayed below this figure's X axis. Defaults to `None`, in
            which case a general-purpose string is defaulted to.
        axes_y_label : optional[str]
            Text displayed to the left of this figure's Y axis. Defaults to
            `None`, in which case a general-purpose string is defaulted to.
        colormap : optional[Colormap]
            Matplotlib colormap to be used by default for all animation artists
            (e.g., colorbar, images). Defaults to `None`, in which case the
            default colormap is used.
        layers : optional[SequenceTypes]
            Sequence of all :class:`LayerCellsABC` instances collectively
            plotting each frame of this plot or animation. **Order is extremely
            significant.** Specifically, the order of layers in this sequence
            defines the order in which these layers are plotted and hence
            overlaid onto one another (i.e., z-order). Defaults to `None`, in
            which case the subclass is responsible for manually plotting this
            plot or animation.
        '''

        # Classify core parameters with weak rather than strong (the default)
        # references, thus avoiding circular references and the resulting
        # complications thereof (e.g., increased memory overhead). Since these
        # objects necessarily live significantly longer than this plot, no
        # complications arise. Ergo, these attributes *ALWAYS* yield these
        # objects rather than non-deterministically yielding "None" if these
        # objects are unexpectedly garbage-collected.
        self._sim = references.proxy_weak(sim)
        self._p = references.proxy_weak(p)

        #FIXME: For currently unknown reasons, this object occasionally retains
        #the only remaining reference to the passed "Cells" instance -- which
        #*CANNOT* therefore be safely classified as a weak reference. This is
        #highly unexpected, however, and should thus be investigated.
        # self.cells = references.proxy_weak(cells)
        self._cells = cells

        # Default unpassed parameters.
        if colormap is None:
            colormap = p.default_cm
        if layers is None:
            layers = ()

        # Classify *AFTER* validating parameters.
        self._axes_title = axes_title
        self._axes_x_label = axes_x_label
        self._axes_y_label = axes_y_label
        self._colorbar_title = colorbar_title
        self._colormap = colormap
        self._is_color_autoscaled = is_color_autoscaled
        self._is_save = is_save
        self._is_show = is_show
        self._figure_title = figure_title
        self._label = label

        # Convert the passed sequence of layers into a new list of layers,
        # permitting this sequence to be safely modified *WITHOUT* modifying
        # the passed sequence. To validate the type of each such layer, call an
        # existing method rather than manually perform this conversion.
        self._layers = []
        self._append_layer(*layers)

        # If autoscaling colors, ignore the passed minimum and maximum.
        if is_color_autoscaled:
            self._color_max = None
            self._color_min = None
        # Else, classify the passed minimum and maximum.
        else:
            self._color_max = color_max
            self._color_min = color_min

        # Classify attributes to be subsequently defined.
        self._color_mappables = None
        self._writer_frames = None
        self._writer_video = None

        # Initialize this plot's figure.
        self._init_figure()

    # ..................{ INITIALIZERS ~ figure              }..................
    def _init_figure(self) -> None:
        '''
        Initialize this plot's figure.

        If an axes title is defined, the mandatory figure title is added to
        this figure as a "supertitle" and that axes title to this figure axes
        as a "subtitle;" else, this figure title is added to this figure axes.
        Likewise, the colorbar title is added to this figure's colorbar.
        '''

        # Figure encapsulating this animation as a weak rather than strong (the
        # default) reference, avoiding circular references and complications
        # thereof (e.g., memory overhead). Figures created by the "pyplot" API
        # are internally retained in Matplotlib's "Gcf" figure cache until
        # explicitly closed -- either non-interactively by a close() call or
        # interactively by the corresponding GUI window being closed. Hence,
        # strong figure references should typically *NOT* be retained.
        self._figure = references.proxy_weak(pyplot.figure())

        # Figure axes scaled to the extent of the current 2D environment as a
        # weak rather than strong (the default) reference, thus avoiding
        # circular references and complications thereof (e.g., memory
        # overhead). Since figures already contain their axes as a strong
        # reference, we need *NOT* do so as well here.
        self._axes = references.proxy_weak(pyplot.subplot(111))

        # If this object was initialized with both a figure and axes title,
        # display the former above the latter.
        if self._axes_title:
            # logs.log_debug('Setting supertitle!')
            self._figure.suptitle(
                self._figure_title, fontsize=14, fontweight='bold')
        # Else, display the figure title as the axes title.
        else:
            # logs.log_debug('Setting normal title!')
            self._axes_title = self._figure_title
        assert types.is_str_nonempty(self._axes_title), (
            types.assert_not_str_nonempty(self._axes_title, 'Axis title'))

        # Initialize the X and Y axes of this plot's figure *AFTER* classifying
        # the axes title above.
        self._init_figure_axes()


    def _init_figure_axes(self) -> None:
        '''
        Initialize the X and Y axes of this plot's figure.
        '''

        # Extent of the current 2D environment.
        self._axes_bounds = [
            self._cells.xmin * self._p.um,
            self._cells.xmax * self._p.um,
            self._cells.ymin * self._p.um,
            self._cells.ymax * self._p.um,
        ]

        # Bound these axes by this extent.
        self._axes.axis('equal')
        self._axes.axis(self._axes_bounds)

        # Display passed human-readable strings as axes attributes.
        self._axes.set_xlabel(self._axes_x_label)
        self._axes.set_ylabel(self._axes_y_label)
        self._axes.set_title(self._axes_title)


    def _reinit_figure_axes(self) -> None:
        '''
        Reinitialize the X and Y axes of this plot's figure, removing all
        artists previously drawn to these axes.

        This method is typically called by animation subclasses _before_
        redrawing cell data of the current animation frame onto these axes.
        '''

        # Remove all artists previously drawn to these axes.
        self._axes.clear()

        # Reinitialize these axes.
        self._init_figure_axes()

    # ..................{ PREPARERS ~ figure                 }..................
    @type_check
    def _prep_figure(
        self,
        color_data: SequenceOrNoneTypes = None,

        #FIXME: This is awful. For sanity, require callers always pass a
        #sequence of mappables: e.g.,
        #   color_mappables: IterableOrNoneTypes = None,

        color_mappables: (ScalarMappable,) + IterableOrNoneTypes = None,
    ) -> None:
        '''
        Prepare this plot _after_ having previously initialized this plot but
        _before_ subsequently displaying and/or saving this plot.

        Parameters
        ----------
        color_mappables : optional[ScalarMappable or IterableTypes]
            Iterable of all :class:`ScalarMappable` instances (e.g.,
            :class:`AxesImage`, :class:`ContourSet`) to associate with this
            plot or animation's colorbar. For convenience, this parameter may
            be either:
            * A single mappable, in which case this colorbar will be associated
              with this mappable as is.
            * A non-string iterable (e.g., :class:`list`, :class:`set`) of one
              or more mappables, in which case this colorbar will be
              arbitrarily associated with the first mappable in this iterable.
            Defaults to `None`, in which case this parameter defaults to the
            iterable of all mappables provided by the topmost and hence last
            mappable layer in the current layer sequence (i.e., the
            :meth:`LayerCellsMappableABC.color_mappables` property of the
            last instance of the :class:`LayerCellsMappableABC` subclass in
            the :attr:`_layers` attribute).
        color_data : optional[SequenceTypes]
            Multi-dimensional sequence of all color values to be plotted _or_
            `None` if calculating these values on initialization is impractical
            (e.g., due to space or time constraints). Defaults to `None`. If
            colorbar autoscaling is enabled _and_ this parameter is:
            * Non-`None`, the colorbar is clipped to the minimum and maximum
              scalar values unravelled from this array.
            * `None`, the subclass is responsible for colorbar autoscaling.
        '''

        # Prepare all layers to be layered onto this plot or animation *BEFORE*
        # autoscaling colors assuming layers to have been prepared.
        self._prep_layers()

        # Autoscale colors to the range implied by the passed color values.
        self._autoscale_colors(color_data)

        # If a single mappable rather than a sequence of mappables was passed,
        # convert the former to the latter.
        if isinstance(color_mappables, ScalarMappable):
            color_mappables = (color_mappables,)

        # Associate these mappables with this plot or animation's colorbar
        # *AFTER* preparing all layers defining these mappables if no mappables
        # were passed.
        self._automap_colors(color_mappables)

    # ..................{ DEINITIALIZERS                     }..................
    def close(self) -> None:
        '''
        Destroy this plot and deallocate all memory associated with this plot.

        To reduce matplotlib's memory overhead, this method (in order):

        . Explicitly closes this plot's figure.
        . Explicitly breaks all circular references between this plot's figure
          and related artist objects (e.g., between this figure and its axes).
        . Explicitly nullifies _all_ attributes of the current object.
        . Explicitly garbage collects.

        This method should only be called:

        * When this plot is non-blocking (e.g., being non-interactively saved
          rather than interactively displayed).
        * As the last action of this plot's subclass or caller.

        Attempting to subsequently call any other plot method _or_ access any
        plot field will reliably result in raised exceptions.
        '''

        # If this figure still exists, explicitly close it.
        if self._figure is not None:
            pyplot.close(self._figure)

        # For each name and value of a field bound to this object...
        for field_name, field_value in objs.iter_vars_simple_custom(self):
            # If this field itself contains a "figure" attribute, explicitly
            # nullify the latter to break this figure's circular references in a
            # manner ignoring "AttributeError: can't set attribute" exceptions.
            #
            # Note that this probably fails to break all such references, as
            # doing so appears to be infeasible in a general-purpose manner.
            # These references are baked into the Matplotlib API at a low level!
            try:
                if  hasattr(field_value, 'figure'):
                    setattr(field_value, 'figure', None)
            except AttributeError:
                pass

            # Explicitly nullify all attributes of this object, again ignoring
            # "AttributeError: can't set attribute" exceptions.
            try:
                setattr(self, field_name, None)
            except AttributeError:
                pass

        #FIXME: Unfortunately, explicitly garbage collecting immediately after
        #closing this animation appears to induce instabilities. Under the
        #"TkAgg" backend, for example, the following non-fatal runtime warning
        #is printed (presumably by the underlying Tcl/Tk library):
        #
        #    can't invoke "event" command:  application has been destroyed
        #        while executing
        #    "event generate $w <<ThemeChanged>>"
        #        (procedure "ttk::ThemeChanged" line 6)
        #        invoked from within
        #    "ttk::ThemeChanged"
        #
        #This suggests that if Python happens to garbage collect at this time,
        #the same issue will arise -- but hopefully much less frequently. Until
        #we resolve the underlying cause, we have no choice but to disable this.
        #Untrained glow bugs in a listless summer swelter!

        # Explicitly garbage collect.
        # gc.collect()

    # ..................{ PROPERTIES ~ read-only             }..................
    # Read-only properties, preventing callers from resetting these attributes.

    @property
    def cells(self) -> 'betse.science.cells.Cells':
        '''
        Current cell cluster.
        '''

        return self._cells


    @property
    def p(self) -> 'betse.science.parameters.Parameters':
        '''
        Current simulation configuration.
        '''

        return self._p


    @property
    def sim(self) -> 'betse.science.sim.Simulator':
        '''
        Current simulation.
        '''

        return self._sim


    @property
    def axes(self) -> Axes:
        '''
        Matplotlib axes for this plot or animation's figure.

        All modelled variables for this cell cluster are spatially plotted onto
        this axes at each time step of this plot or animation.
        '''

        return self._axes


    @property
    def color_min(self) -> float:
        '''
        Minimum color value displayed on this plot or animation's colorbar.
        '''

        return self._color_min


    @property
    def color_max(self) -> float:
        '''
        Maximum color value displayed on this plot or animation's colorbar.
        '''

        return self._color_max


    @property
    def colormap(self) -> Colormap:
        '''
        Matplotlib colormap, mapping all numeric data for one modelled variable
        for this cell cluster into color values displayed on this plot or
        animation's colorbar.
        '''

        return self._colormap

    # ..................{ COLORS                             }..................
    #FIXME: Improve this method to internally ignore the first time step of
    #the passed array (i.e., "color_data[0]") if this array is non-empty. Why?
    #Because this step typically contains spurious outlier data resulting in
    #the more "normal" data plotted for the remaining time steps appear to
    #exhibit no changes in color. Ignoring such outlier data should improve
    #this lamentable situation.
    @type_check
    def _autoscale_colors(self, color_data: SequenceOrNoneTypes) -> None:
        '''
        Autoscale the colorbar for this plot or animation's figure to the
        minimum and maximum scalar values unravelled from the passed sequence
        by setting the :attr:`_color_min` and :attr:`_color_max` attributes to
        such values if colorbar autoscaling is both enabled and has not already
        been performed (i.e., these attributes are `None`) _or_ noop otherwise.

        Parameters
        ----------
        color_data : SequenceOrNoneTypes
            Multi-dimensional sequence of all color values to be plotted _or_
            `None` if calculating these values on initialization is
            impractical. See the :meth:`_prep_figure` method for further
            details.
        '''

        # If colorbar autoscaling is disabled, noop.
        if not self._is_color_autoscaled:
            return

        # If no color values are passed...
        if color_data is None:
            #FIXME: Eliminate code duplication between this and the
            #_automap_colors() method. Ideally, the last mappable layer
            #instance should only be searched for once. This code is currently
            #non-trivial to unify, due to the need for a mappable layer to be
            #optional here but *NOT* in the _automap_colors() method. Ideally,
            #a mappable layer should be mandatory in both methods. Once this is
            #the case across all animations, unify the retrieval of this layer
            #into a single location.

            # Last mappable layer in this layer sequence if any or the sentinel
            # placeholder constant otherwise.
            mappable_layer = iterables.get_item_last_instance_of_or_sentinel(
                iterable=self._layers,
                cls=LayerCellsMappableABC,
            )

            # If this layer exists, defer to its color data.
            if mappable_layer is not SENTINEL:
                color_data = mappable_layer.color_data

        # If no color values are available, silently noop.
        if color_data is None:
            return

        # Flatten this multi-dimensional array to a one-dimensional array,
        # permitting efficient retrieval of minimum and maximum values.
        time_series_flat = np.ravel(color_data)

        # Set the current minimum and maximum color values.
        self._color_min = np.ma.min(time_series_flat)
        self._color_max = np.ma.max(time_series_flat)

        # Log these values.
        logs.log_debug('Autoscaling animation colors to [%d, %d].', self._color_min, self._color_max)


    @type_check
    def _automap_colors(self, color_mappables: IterableOrNoneTypes) -> None:
        '''
        Create a figure colorbar for this plot or animation associated with the
        passed mappables if any or the mappables provided by the last mappable
        layer in the current layer sequence otherwise.

        Parameters
        ----------
        color_mappables : optional[IterableTypes]
            Iterable of all :class:`ScalarMappable` instances (e.g.,
            :class:`AxesImage`, :class:`ContourSet`) to associate with this
            plot or animation's colorbar. By Matplotlib design, only the first
            mappable in this iterable is arbitrarily associated with this
            colorbar; all other mappables are ignored. Defaults to `None`, in
            which case the iterable of all mappables provided by the topmost
            and hence last mappable layer in the current layer sequence is
            defaulted to (i.e., the value of the
            :meth:`LayerCellsMappableABC.color_mappables` property of the last
            instance of the :class:`LayerCellsMappableABC` subclass in the
            :attr:`_layers` attribute).
        '''

        # If no color mappables are passed...
        if color_mappables is None:
            # Last mappable layer in this layer sequence if any or raise an
            # exception with this message otherwise.
            mappable_layer = iterables.get_item_last_instance_of(
                iterable=self._layers,
                cls=LayerCellsMappableABC,
                exception_message=(
                    'Visual "{}" mappable layer not found.'.format(
                        self._label)),
            )

            # Sequence of mappables provided by this layer.
            color_mappables = mappable_layer.color_mappables

        # Classify this sequence of mappables.
        self._color_mappables = color_mappables

        # Scale these mappables by previously established color values *AFTER*
        # classifying this sequence.
        self._rescale_color_mappables()

        # First mappable safely retrieved from this iterable of mappables.
        color_mappable_first = iterables.get_item_first(self._color_mappables)

        # Create a colorbar associated with this mappable.
        colorbar = self._figure.colorbar(color_mappable_first)
        colorbar.set_label(self._colorbar_title)


    def _rescale_color_mappables(self) -> None:
        '''
        Rescale all mappables associated with this plot or animation's colorbar
        to the current minimum and maximum color values.

        This method must be called _after_ the :meth:`_prep_figure` method has
        been called. Failing to do so will raise exceptions.
        '''

        # If the _prep_figure() method has not been called, raise an exception.
        if self._color_mappables is None:
            raise BetseMethodException(
                '{class_name}._rescale_color_mappables() called before '
                '{class_name}._prep_figure().'.format(class_name=type(self)))

        # If these values are identical, coerce them to differ. Failing to do
        # so produces spurious visual artifacts in both the axes and colorbar.
        # if self._color_min == self._color_max:
        #     self._color_min = self._color_min - 1
        #     self._color_max = self._color_max + 1

        # For each color mappable, clip that mappable to the minimum and
        # maximum values discovered above. Note this also has the beneficial
        # side-effect of establishing the colorbar's range.
        for color_mappable in self._color_mappables:
            assert types.is_matplotlib_mappable(color_mappable), (
                types.assert_not_matplotlib_mappable(color_mappable))
            color_mappable.set_clim(self._color_min, self._color_max)

    # ..................{ LAYERS                             }..................
    @type_check
    def _append_layer(self, *layers: LayerCellsABC) -> None:
        '''
        Append all passed layers to the current sequence of layers, ensuring
        the subsequently called :meth:`_visualize_layers` method will visualize
        these layers after (and hence above) all previously appended layers.

        Parameters
        ----------
        layers : Tuple[LayerCellsABC]
            Tuple of all layers to be appended to this sequence of layers.
        '''

        self._layers.extend(layers)


    def _prep_layers(self) -> None:
        '''
        Iteratively prepare all layers to be subsequently layered onto this
        plot or animation.

        If this simulation configuration requests that cells be labelled by
        their 0-based indices, a layer doing so is appended to the current
        layer sequence. Since this method is called _after_ all user-defined
        layers have been appended, these labels are guaranteed to be layered
        over rather than under all plotted cell data.
        '''

        # If labelling each cell with its 0-based index, append a layer doing
        # so *AFTER* all lower layers (e.g., cell data) have been appended,
        if self._p.enumerate_cells:
            self._append_layer(LayerCellsIndex())

        # Prepare each layer *AFTER* appending all layers above.
        for layer in self._layers:
            layer.prep(self)


    def _plot_layers(self) -> None:
        '''
        Iteratively plot all layers onto this plot or animation for the
        current time step of this simulation.

        If this is:

        * The first time step, each such layer adds one or more Matplotlib
          artists to the figure for this plot or animation.
        * Any time step _except_ the first, each such layer either:
          * Updates the contents of all artists previously added by that layer.
            This is the most efficient and hence ideal approach, but infeasible
            in cases in which updating artists is infeasible (e.g., the
            streamlines of a streamplot).
          * Removes all artists previously added by that layer _and_ adds new
            artists whose contents reflect this time step. This is the least
            efficient and hence least ideal approach.
        '''

        for layer in self._layers:
            layer.layer()

    # ..................{ PLOTTERS                           }..................
    @type_check
    def _plot_image(
        self,

        # Mandatory parameters.
        pixel_data: SequenceTypes,

        # Optional parameters.
        colormap : Colormap = None,
    ) -> AxesImage:
        '''
        Plot and return an image of the passed pixel data onto the current
        figure's axes.

        Parameters
        ----------
        pixel_data: np.ndarray
            Array of pixel data defining the image to be plotted, whose first
            two dimensions index the X and Y components of this axes. If this
            array is:
            * Two-dimensional, a greyscale image mapped onto the passed
              colormap will be plotted.
            * Three-dimensional, an RGB image _not_ mapped onto the passed
              colormap will be plotted.
            * Four-dimensional, an RGBa image _not_ mapped onto the passed
              colormap will be plotted.
        colormap : matplotlib.cm.Colormap
            Optional colormap with which to map the passed pixel data when
            greyscale (and ignored otherwise) _or_ `None`, in which case the
            default colormap will be used.

        Returns
        ----------
        matplotlib.image.AxesImage
            Image produced by plotting the passed pixel data.

        See Also
        ----------
        http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.imshow
            Further details on Matplotlib-based axes image plotting.
        '''

        # Default unpassed parameters.
        if colormap is None:
            colormap = self._colormap

        # Plot and return this image.
        return self._axes.imshow(
            pixel_data,
            origin='lower',
            extent=self._axes_bounds,
            cmap=colormap,
        )


    #FIXME: Replace entirely by the appropriate "LayerCellsStream" subclass.
    @type_check
    def _plot_stream(
        self,

        # Mandatory arguments.
        x: SequenceTypes,
        y: SequenceTypes,
        magnitude: SequenceTypes,

        # Optional arguments.
        grid_x: SequenceOrNoneTypes = None,
        grid_y: SequenceOrNoneTypes = None,
        magnitude_max: NumericOrNoneTypes = None,
        old_stream_plot: (StreamplotSet, NoneType) = None,
    ) -> StreamplotSet:
        '''
        Plot and return a streamplot of all streamlines of the passed vector
        flow data onto the current figure's axes.

        Parameters
        ----------
        x : SequenceTypes
            Two-dimensional sequence of X components of vector flow velocity.
        y : SequenceTypes
            Two-dimensional sequence of Y components of vector flow velocity.
        magnitude: SequenceTypes
            One-dimensional sequence of vector flow magnitudes.
        grid_x : SequenceTypes or NoneType
            Optional scaled X components of the cell cluster grid. Defaults to
            `None`, in which case the default `self.cells.X` array is used.
        grid_y : SequenceTypes or NoneType
            Optional scaled Y components of the cell cluster grid. Defaults to
            `None`, in which case the default `self.cells.Y` array is used.
        magnitude_max: NumericTypes or NoneType
            Optional maximum magnitude in the passed `magnitude` array.
            Defaults to `None`, in which case this array is implicitly searched
            for this value.
        old_stream_plot: StreamplotSet or NoneType
            Optional streamplot returned by a prior call to this method
            (typically for a prior frame) _or_ `None` if this is the first call
            to this method for this animation. If non-`None`, this streamplot
            will be cleared in preparation for re-streamplotting by this call.

        Returns
        ----------
        StreamplotSet
            Streamplot produced by plotting the passed vector flow data.

        See Also
        ----------
        http://matplotlib.org/api/axes_api.html#matplotlib.axes.Axes.streamplot
            Further details on Matplotlib-based axes streamplotting.
        '''

        # Default all unpassed optional arguments.
        if magnitude_max is None:
            magnitude_max = np.max(magnitude)
        if grid_x is None:
            grid_x = self._cells.X * self._p.um
        if grid_y is None:
            grid_y = self._cells.Y * self._p.um

        # If a prior streamplot to be erased was passed, do so.
        if old_stream_plot is not None:
            # Erase this streamplot's streamlines before replotting.
            old_stream_plot.lines.remove()

            # If this version of Matplotlib supports erasing the set of patches
            # corresponding to this streamplot's arrow heads, do so.
            try:
                old_stream_plot.arrows.remove()
            # Else, manually erase them by iterating over all patch objects and
            # preserving all non-arrow head patches. This will also remove all
            # arrow head patches of streamplots plotted by this subclass for
            # this frame, which is non-ideal. Matplotlib leaves us no choice.
            except NotImplementedError:
                self._axes.patches = [
                    patch
                    for patch in self._axes.patches
                    if not isinstance(patch, FancyArrowPatch)
                ]

        # Plot and return this streamplot.
        return self._axes.streamplot(
            grid_x, grid_y, x, y,
            density=self._p.stream_density,
            linewidth=(3.0*magnitude/magnitude_max) + 0.5,
            color=self._p.vcolor,
            cmap=self._colormap,
            arrowsize=3.0,

            # Draw this streamplot over all patch and line artists, by default.
            # See the "ZORDER_STREAM" docstring for further commentary.
            zorder=ZORDER_STREAM,
        )

    # ..................{ PLOTTERS ~ cell                    }..................
    #FIXME: Pretty intense, and obviously better refactored two distinct
    #"LayerCellsABC" subclasses. This will probably prove pivotal to
    #implementing deformations sanely.
    #FIXME: After doing so, excise *ALL* of the methods below.

    def _plot_cells_sans_ecm(self, *args, **kwargs) -> 'Collection':
        '''
        Plot and return an intracellular plot of all cells with colours
        corresponding to the passed vector of arbitrary cell data (e.g.,
        transmembrane cell voltages for the current time step) onto the current
        figure's axes.

        The type of plot returned is defined by this simulation's configuration.
        Specifically:

        * If this configuration requests that individual cells be plotted (i.e.,
          `show cells` is `True`), a mosaic plot of this simulation's cell
          cluster will be returned. See the `_plot_cell_mosaic()` method.
        * Else, a mesh plot interpolating the individual cells of this
          simulation's cell cluster will be returned. See the
          `_plot_cell_mesh()` method.

        Parameters
        -----------
        cell_data : np.ndarray
            Arbitrary cell data defined on an environmental grid to be plotted.

        All other passed parameters will be passed as is to the underlying plot
        method called by this method (e.g., `_plot_cell_mosaic()`).

        Returns
        --------
        Collection
            Plot produced by plotting the passed cell data.
        '''

        if self._p.showCells is True:
            return self._plot_cell_mosaic(*args, **kwargs)
        else:
            return self._plot_cell_mesh(*args, **kwargs)


    def _update_cell_plot_sans_ecm(
        self,
        cell_plot: 'Collection',
        cell_data: np.ndarray,
        *args, **kwargs
    ) -> None:
        '''
        Update _without_ recreating the passed intracellular cell plot with the
        passed vector of arbitrary cell data (e.g., transmembrane cell voltages
        for the current time step) onto the current figure's axes.

        Parameters
        -----------
        cell_plot : Collection
            Cell plot previously returned by either the `_plot_cells_sans_ecm()`
            _or_ `_revive_cell_plot_sans_ecm()` method.
        cell_data : np.ndarray
            Arbitrary cell data defined on an environmental grid to be plotted.

        All other passed parameters will be passed as is to the underlying plot
        method called by this method (e.g., `_plot_cell_mosaic()`).
        '''
        assert types.is_sequence_nonstr(cell_data), (
            types.assert_not_sequence_nonstr(cell_data))

        # If plotting individuals cells, update this plot with this data as is.
        if self._p.showCells:
            assert types.is_matplotlib_polycollection(cell_plot), (
                types.assert_not_matplotlib_polycollection(cell_plot))
            cell_plot.set_array(cell_data)
        # Else, the cell cluster is being plotted as an unstructured triangular
        # grid, requiring this data be reshaped onto this grid.
        else:
            assert types.is_matplotlib_trimesh(cell_plot), (
                types.assert_not_matplotlib_trimesh(cell_plot))

            # membranes_midpoint_data = np.zeros(len(self.cells.voronoi_centres))
            # membranes_midpoint_data[self.cells.cell_to_grid] = membranes_midpoint_data

            # Update this plot with this gridded data.
            cell_plot.set_array(cell_data)


    def _revive_cell_plots_sans_ecm(
        self,
        cell_plot: 'Collection',
        cell_data: np.ndarray,
        *args, **kwargs
    ) -> 'Collection':
        '''
        Recreate the passed intracellular plot and return a similar plot of all
        cells with colours corresponding to the passed vector of arbitrary cell
        data (e.g., transmembrane cell voltages for the current time step) onto
        the current figure's axes.

        This method is typically called to replot an animation of a cell cluster
        subject to physical cell changes (e.g., cutting, deformation). The type
        of plot returned is defined by this simulation's configuration.
        Specifically:

        * If this configuration requests that individual cells be plotted (i.e.,
          `show cells` is `True`), the passed plot _must_ be a mosaic plot
          previously returned by the `_plot_cell_mesh()` method. This plot will
          be updated in-place and returned as is without creating a new plot.
        * Else, the passed plot _must_ be a mesh plot previously returned by the
          `_plot_cell_mesh()` method -- specifically, a `TriMesh` object. Since
          the `TriMesh` API currently provides no public means of updating such
          plots in-place, this method (in order):
          . Removes the passed mesh plot from this figure's axes.
          . Creates a new mesh plot and adds that plot to this figure's axes.
          . Returns that plot.

        Parameters
        -----------
        cell_plot : Collection
            Cell plot previously returned by either the `_plot_cells_sans_ecm()` _or_
            `_revive_cell_plots_sans_ecm()` method.
        cell_data : np.ndarray
            Arbitrary cell data defined on an environmental grid to be plotted.

        All other passed parameters will be passed as is to the underlying plot
        method called by this method (e.g., `_plot_cell_mosaic()`).

        Returns
        --------
        Collection
            Plot produced by replotting the passed cell data.
        '''
        assert types.is_sequence_nonstr(cell_data), (
            types.assert_not_sequence_nonstr(cell_data))

        # If plotting individual cells, the passed cell plot *MUST* be a polygon
        # collection previously returned by the _plot_cell_mosaic() method.
        if self._p.showCells is True:
            assert types.is_matplotlib_polycollection(cell_plot), (
                types.assert_not_matplotlib_polycollection(cell_plot))

            # Update this plot in-place.
            cell_plot.set_array(cell_data)
            cell_plot.set_verts(
                arrays.from_sequence(self._cells.cell_verts) * self._p.um)

            # Return the same plot.
            return cell_plot
        # Else, the passed cell plot *MUST* be a triangle mesh previously
        # returned by the _plot_cell_mesh() method.
        else:
            assert types.is_matplotlib_trimesh(cell_plot), (
                types.assert_not_matplotlib_trimesh(cell_plot))

            cell_plot.remove()
            return self._plot_cell_mesh(cell_data=cell_data, *args, **kwargs)


    def _plot_cell_mosaic(self, cell_data: SequenceTypes) -> PolyCollection:
        '''
        Plot and return a mosaic plot of all cells with colours corresponding
        to the passed vector of arbitrary cell data (e.g., transmembrane
        voltages for all cells for the current time step) onto the current
        figure's axes.

        The returned plot will be a polygon collection such that each polygon
        signifies a cell in this simulation's cell cluster.

        Parameters
        -----------
        cell_data : SequenceTypes
            Arbitrary cell data defined on an environmental grid to be plotted.

        Returns
        --------
        PolyCollection
            Mosaic plot produced by plotting the passed cell data.
        '''

        # Cell vertices plotted as polygons.
        mosaic_plot = PolyCollection(
            verts=visuals.upscale_cell_coordinates(self._cells.cell_verts),
            cmap=self._colormap,
            edgecolors='none',
        )

        # Associate this plot with the passed cell data. Curiously, the
        # Matplotlib API provides no means of passing this data to the
        # PolyCollection.__init__() constructor. (Yes, we checked. Twice.)
        mosaic_plot.set_array(cell_data)

        # Add this plot to this figure's axes.
        self._axes.add_collection(mosaic_plot)
        return mosaic_plot


    #FIXME: This plots somewhat similarly to the presumably superior
    #"LayerCellsShadeDiscrete" subclass. Generalize this method into a new
    #"LayerCellsGouraudContinuous" subclass of the same submodule.
    def _plot_cell_mesh(self, cell_data: np.ndarray) -> 'TriMesh':
        '''
        Plot and return a mesh plot of all cells with colours corresponding to
        the passed vector of arbitrary cell data (e.g., transmembrane voltages
        for all cells for the current time step) onto the current figure's axes.

        The returned plot will be an unstructured triangular grid interpolating
        each cell of this simulation's cell cluster into a smooth continuum.

        Parameters
        -----------
        cell_data : np.ndarray
            Arbitrary cell data defined on an environmental grid to be plotted.

        Returns
        --------
        TriMesh
            Mesh plot produced by plotting the passed cell data.
        '''
        assert types.is_sequence_nonstr(cell_data), (
            types.assert_not_sequence_nonstr(cell_data))

        # If the passed cell data is defined on membrane midpoints, average that
        # to correspond to cell centres instead.
        if len(cell_data) == len(self._cells.mem_i):
            cell_data = np.dot(
                self._cells.M_sum_mems, cell_data) / self._cells.num_mems

        # Unstructured triangular grid assigned the passed cell data.
        triangular_grid = np.zeros(len(self._cells.voronoi_centres))
        triangular_grid[self._cells.cell_to_grid] = cell_data

        # cmin = membranes_midpoint_data.min()
        # cmax = membranes_midpoint_data.max()

        # Create and add this plot to this figure's axes and return this plot.
        # Unfortunately, the tripcolor() function requires:
        #
        # * The "x", "y", and "triangles" parameters to be positional.
        # * All other parameters to be keyword.
        #
        # Behold! The ultimate examplar of nonsensical API design.
        return self._axes.tripcolor(
            self._p.um * self._cells.cell_centres[:, 0],
            self._p.um * self._cells.cell_centres[:, 1],
            cell_data,
            shading='gouraud',
            cmap=self._colormap
        )
