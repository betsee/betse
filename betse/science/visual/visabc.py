#!/usr/bin/env python3
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Abstract base classes of all Matplotlib-based plot and animation subclasses.
'''

#FIXME: Refactor all procedural cell cluster-specific
#"betse.science.visual.plot.plotutil" functions into subclasses of the
#"LayerCellsABC" base class defined elsewhere.
#Ultimate power fights the dark deceit!

# ....................{ IMPORTS                            }....................
import numpy as np
from abc import ABCMeta  #, abstractmethod  #, abstractstaticmethod
from betse.exceptions import BetseSimVisualException
from betse.lib.matplotlib import mplutil
from betse.lib.matplotlib.matplotlibs import mpl_config
from betse.lib.matplotlib.mplzorder import ZORDER_PATCH, ZORDER_STREAM
from betse.lib.numpy import nparray
from betse.science.config.export.visual.confvisabc import SimConfVisualCellsABC
from betse.science.export import expmath
from betse.science.phase.phasecls import SimPhase
from betse.science.visual.layer.lyrabc import LayerCellsABC
from betse.science.visual.layer.lyrtext import LayerCellsIndex
from betse.util.io.log import logs
from betse.util.py import pyref
from betse.util.type import iterables, types
from betse.util.type.obj import objects
from betse.util.type.types import (
    type_check,
    IterableTypes,
    NoneType,
    NumericOrNoneTypes,
    SequenceTypes,
    SequenceOrNoneTypes,
    StrOrNoneTypes,
)
from matplotlib import pyplot
from matplotlib.axes import Axes
from matplotlib.cm import ScalarMappable
from matplotlib.collections import PolyCollection
from matplotlib.colors import Colormap
from matplotlib.figure import Figure
from matplotlib.image import AxesImage
from matplotlib.patches import FancyArrowPatch
from matplotlib.streamplot import StreamplotSet

# ....................{ SUPERCLASSES                       }....................
class VisualCellsABC(object, metaclass=ABCMeta):
    '''
    Abstract base class of all **cell cluster visual** (i.e., plot or animation
    spatially visualizing the current cell cluster) subclasses.

    Subclasses of this class plot the spatial distribution of one or more
    modelled variables (e.g., charge distribution, membrane voltage) for
    either:

    * A single simulation time step, in which this subclass plots a still frame
      of these variables for this step.
    * One or more simulation time steps, in which this subclass animates a
      video of these variables for these steps.

    Attributes (Private)
    ----------
    _conf : SimConfVisualCellsABC
        Configuration for this visual, synchronized with the user-defined
        YAML-backed simulation configuration file for this phase.
    _label : str
        Basename of the subdirectory in the phase-specific results directory
        to which all files exported for this visual are saved _and_
        the basename prefix of these files.
    _layers : list
        List of all :class:`LayerCellsABC` instances collectively composing
        this visual.
    _phase : SimPhase
        Current simulation phase.

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
    _axes_title : str
        Text displayed above the figure axes. If a non-``None`` value for the
        ``axes_title`` parameter is passed to the :meth:`__init__` method, this
        is that value; else, this is the value of the ``figure_title`` parameter
        passed to that method.
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
        animation. Due to matplotlib constraints, only the first mappable in
        this iterable is associated with this colorbar.
    _colorbar_title: StrOrNoneTypes
        Text displayed above the figure colorbar if any *or* ``None`` otherwise.
    _color_max : NumericSimpleTypes
        Maximum color value to be displayed by the colorbar.
    _color_min : NumericSimpleTypes
        Minimum color value to be displayed by the colorbar.
    _colormap : Colormap
        Matplotlib colormap with which to create this animation's colorbar.

    Attributes (Private: Saving)
    ----------
    _writer_savefig_kwargs : dict
        Dictionary of all keyword arguments to be passed to the
        :meth:`Figure.savefig` method called to save each visual frame for
        both still images and animated video. By default, this dictionary:
        * Enables the ``transparent`` argument, painting frame backgrounds
          transparently rather than as a solid color (usually, white).

    Attributes (Private: Time)
    ----------
    _is_time_step_first : bool
        ``True`` only if the first frame either has yet to be plotter *or* is
        currently being plotted.
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(
        self,

        # Mandatory parameters.
        phase: SimPhase,
        conf: SimConfVisualCellsABC,
        is_save: bool,
        is_show: bool,

        #FIXME: For orthogonality with nomenclature elsewhere -- including the
        #the "conf.name" attribute -- rename both this parameter and the
        #corresponding "_label" attribute to "name" and "_name" respectively.
        label: str,

        figure_title: str,

        # Optional parameters.
        axes_title: StrOrNoneTypes = None,
        axes_x_label: str = 'Spatial Distance [um]',
        axes_y_label: str = 'Spatial Distance [um]',
        colorbar_title: StrOrNoneTypes = None,
        colormap: (Colormap, NoneType) = None,

        #FIXME: Refactor this to be mandatory instead.
        layers: SequenceOrNoneTypes = None,
    ) -> None:
        '''
        Initialize this visual.

        Parameters
        ----------
        phase : SimPhase
            Current simulation phase.
        conf : SimConfVisualCellsABC
            Configuration for this visual, synchronized with the user-defined
            YAML-backed simulation configuration file for this phase.
        is_save : bool
            ``True`` only if non-interactively saving this visual.
        is_show : bool
            ``True`` only if interactively displaying this visual.
        label : str
            Terse machine-readable string (e.g., ``Vmem``) serving as both:
            * The basename of the subdirectory of the phase-specific results
              directory containing all files saved by this visual.
            * The basename prefix of these files.
        figure_title : str
            Text displayed above the figure itself.
        axes_title : optional[str]
            Text displayed above the figure axes but below the figure title.
            Defaults to `None`, in which case no such text is displayed.
        axes_x_label : optional[str]
            Text displayed below this figure's X axis. Defaults to `None`, in
            which case a general-purpose string is defaulted to.
        axes_y_label : optional[str]
            Text displayed to the left of this figure's Y axis. Defaults to
            `None`, in which case a general-purpose string is defaulted to.
        colorbar_title: optional[str]
            Text displayed above the figure colorbar if any *or* ``None``
            otherwise. If this parameter is ``None`` *and* at least one layer in
            the passed ``layers`` parameter is an instance of the
            :class:`LayerCellsColorfulABC` base class, an exception is raised.
        colormap : optional[Colormap]
            Matplotlib colormap to be used by default for all animation artists
            (e.g., colorbar, images). Defaults to ``None``, in which case the
            default colormap is used.
        layers : optional[SequenceTypes]
            Sequence of all :class:`LayerCellsABC` instances collectively
            plotting each frame of this visual. **Order is extremely
            significant.** Specifically, the order of layers in this sequence
            defines the order in which these layers are plotted and hence
            overlaid onto one another (i.e., z-order). Defaults to `None`, in
            which case the subclass is responsible for manually plotting this
            visual.
        '''

        # Default unpassed parameters.
        if colormap is None:
            colormap = phase.p.default_cm
        if layers is None:
            layers = ()

        # Classify *AFTER* validating parameters.
        self._phase = phase
        self._conf = conf
        self._axes_title = axes_title
        self._axes_x_label = axes_x_label
        self._axes_y_label = axes_y_label
        self._colorbar_title = colorbar_title
        self._colormap = colormap
        self._is_save = is_save
        self._is_show = is_show
        self._figure_title = figure_title
        self._label = label

        # Convert the passed sequence of layers into a new list of layers,
        # permitting this sequence to be safely modified *WITHOUT* modifying
        # the passed sequence. To validate the type of each such layer, call an
        # existing method rather than manually performing this conversion.
        self._layers = []
        self._append_layer(*layers)

        # If autoscaling colors, ignore the passed minimum and maximum. To avoid
        # propagating these changes back to the underlying configuration file,
        # these changes are isolated to new instance variables; the passed
        # "conf.color_max" and "conf.color_min" variables remain unchanged.
        if self._conf.is_color_autoscaled:
            self._color_max = None
            self._color_min = None
        # Else, classify the passed minimum and maximum.
        else:
            self._color_max = conf.color_max
            self._color_min = conf.color_min

        # Default all attributes to be subsequently defined.
        self._color_mappables = None
        self._writer_frames = None
        self._writer_video = None

        # Default all remaining attributes.
        self._is_time_step_first = True

        # Dictionary of keyword arguments to be passed to Figure.savefig().
        self._writer_savefig_kwargs = {
            # Visualize frame backgrounds transparently rather than solid.
            'transparent': True,
        }

        # Initialize this plot's figure *AFTER* defining all attributes.
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
        self._figure = pyref.proxy_weak(pyplot.figure())

        # Figure axes scaled to the extent of the current 2D environment as a
        # weak rather than strong (the default) reference, thus avoiding
        # circular references and complications thereof (e.g., memory
        # overhead). Since figures already contain their axes as a strong
        # reference, we need *NOT* do so as well here.
        self._axes = pyref.proxy_weak(pyplot.subplot(111))

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

        # Bound these axes by this extent.
        self._axes.axis('equal')
        self._axes.axis(self._phase.cache.upscaled.extent)

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
    #FIXME: This method should no longer accept any parameters. All parameters
    #currently accepted by this method are obsolete and hence remain
    #undocumented. They should *ALL* be removed as soon as feasible.
    def _prep_figure(self, *args, **kwargs) -> None:
        '''
        Prepare this visual *after* having previously initialized this
        visual but *before* subsequently displaying and/or saving this
        visual.
        '''

        # Prepare all layers to be layered onto this visual *BEFORE*
        # autoscaling colors assuming layers to have been prepared.
        self._prep_layers()

        #FIXME: Obsolete. Remove after adopting layers everywhere.
        if 'color_data' in kwargs:
            self._prep_colors_obsolete(*args, **kwargs)

    # ..................{ DEINITIALIZERS                     }..................
    def close(self) -> None:
        '''
        Deallocate all resources associated with this visual.

        To reduce resource overhead (namely memory), this method (in order):

        . Explicitly closes this visual's figure.
        . Explicitly breaks all circular references between this figure and
          related artist objects (e.g., between this figure and its axes).
        . Explicitly nullifies *all* attributes of the current object.
        . Explicitly garbage collects.

        Caveats
        ----------
        This method should only be called:

        * If this visual is blocking (e.g., saved non-interactively,
          displayed interactively in a blocking manner). If this plot or
          animation is non-blocking, the resources assigned this plot or
          animation may be safely deallocated only *after* the end user closes
          the corresponding GUI widget or window.
        * As the last action of either:
          * The current caller.
          * The visual subclass.

        Attempting to subsequently call any other method *or* access any other
        variable bound to this object will reliably raise an exception.
        '''

        # If this figure still exists, explicitly close it.
        if self._figure is not None:
            pyplot.close(self._figure)

        # For each name and value of a field bound to this object...
        for field_name, field_value in objects.iter_vars_custom_simple(self):
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
    def phase(self) -> SimPhase:
        '''
        Current simulation phase.
        '''

        return self._phase


    @property
    def conf(self) -> SimConfVisualCellsABC:
        '''
        Configuration for this visual, synchronized with the user-defined
        YAML-backed simulation configuration file for this phase.
        '''

        return self._conf


    @property
    def name(self) -> str:
        '''
        Non-human-readable name identifying this type of visual (e.g.,
        ``voltage_extra``).
        '''

        return self._label

    # ..................{ PROPERTIES ~ read-only : mpl       }..................
    @property
    def figure(self) -> Figure:
        '''
        Matplotlib figure for this visual, containing *all* matplotlib artists
        displayed by this visual (e.g., axes, colorbar).
        '''

        return self._axes


    @property
    def axes(self) -> Axes:
        '''
        Matplotlib axes for this visual's figure.

        All modelled variables for this cell cluster are spatially plotted onto
        this axes at each time step of this visual.
        '''

        return self._axes

    # ..................{ PROPERTIES ~ read-only : color     }..................
    #FIXME: Obsolete. Remove after adopting layers everywhere.
    @property
    def color_min(self) -> float:
        '''
        Minimum color value displayed on this visual's colorbar if
        this .
        '''

        return self._color_min


    #FIXME: Obsolete. Remove after adopting layers everywhere.
    @property
    def color_max(self) -> float:
        '''
        Maximum color value displayed on this visual's colorbar.
        '''

        return self._color_max


    #FIXME: Obsolete. Remove after adopting layers everywhere.
    @property
    def colormap(self) -> Colormap:
        '''
        Matplotlib colormap, mapping all numeric data for one modelled variable
        for this cell cluster into color values displayed on this plot or
        animation's colorbar.
        '''

        return self._colormap

    # ..................{ MAKERS                             }..................
    @type_check
    def make_colorbar(self, color_mappable: ScalarMappable) -> None:
        '''
        Add a colorbar to this visual's figure associated with the passed
        matplotlib artist satisfying the color mappable API.

        This method is typically called by the **colorful layer** (i.e., layer
        displaying data intended to be coloured) defined for this visual.
        '''

        # If no colorbar title was passed to this visual, raise an exception.
        if not self._colorbar_title:
            raise BetseSimVisualException(
                'Animation "{}" parameter "colorbar_title" '
                'undefined or empty.'.format(self._label))

        # Create and configure this colorbar.
        colorbar = self._figure.colorbar(color_mappable)
        colorbar.set_label(self._colorbar_title)

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
        Prepare all layers to be subsequently layered onto this visual.

        If this simulation configuration requests that cells be labelled by
        their 0-based indices, a layer doing so is appended to the current
        layer sequence. Since this method is called _after_ all user-defined
        layers have been appended, these labels are guaranteed to be layered
        over rather than under all plotted cell data.
        '''

        # If labelling each cell with its 0-based index, append a layer doing
        # so *AFTER* all lower layers (e.g., cell data) have been appended,
        if self._phase.p.enumerate_cells:
            self._append_layer(LayerCellsIndex())

        # For the 0-based index of each layer and layer in the layer sequence
        # *AFTER* finalizing this sequence by appending final layers above...
        for layer_index, layer in enumerate(self._layers):
            # Z-order of all artists plotted by this layer. To prevent patch
            # artists defaulting to the default "ZORDER_PATCH" from being
            # plotted over artists plotted by the first layer, this z-order is
            # guaranteed to be at least as large as that default.
            layer_zorder = layer_index + ZORDER_PATCH

            # Prepare this layer.
            layer.prep(visual=self, zorder=layer_zorder)


    def _plot_layers(self) -> None:
        '''
        Iteratively plot all layers onto this visual for the
        current time step of this simulation.

        If this is:

        * The first time step, each such layer adds one or more Matplotlib
          artists to the figure for this visual.
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

    # ..................{ COLORS                             }..................
    #FIXME: All methods in this subsection including this method are obsolete.
    #Remove after adopting layers everywhere. All other "FIXME:" comments below
    #are now ignorable, as this functionality should all go away.
    #FIXME: After doing so, also remove most color-specific attributes,
    #including:
    #
    #* "_color_mappables".
    #* "_color_min".
    #* "_color_max".
    #
    #The two notable exceptions are "_colormap" and "_colorbar_title", which
    #remain suitable for use in this base class.
    @type_check
    def _prep_colors_obsolete(
        self,

        # Obsolete optional parameters, which should no longer be passed.
        color_data: SequenceTypes,

        #FIXME: This is awful. For sanity, require callers always pass a
        #sequence of mappables: e.g.,
        #   color_mappables: IterableOrNoneTypes = None,
        #FIXME: Actually, this entire method is awful and should be refactored
        #away once layers are used everywhere. For backward compatilibility with
        #old-style animations, this is currently preserved as is.

        color_mappables: (ScalarMappable,) + IterableTypes,
    ) -> None:
        '''
        Prepare all color-specific data if this visual is colorful *or* noop
        otherwise (i.e., if this visual is colorless).

        If this visual is colorful, this method autoscale colors to the
        range implied by the passed color values *and* associate all mappables
        provided by colorful layers with this visuazilation's colorbar.

        Parameters
        ----------
        color_mappables : optional[ScalarMappable or IterableTypes]
            Iterable of all :class:`ScalarMappable` instances (e.g.,
            :class:`AxesImage`, :class:`ContourSet`) to associate with this
            visual's colorbar. For convenience, this parameter may
            be either:
            * A single mappable, in which case this colorbar will be associated
              with this mappable as is.
            * A non-string iterable (e.g., :class:`list`, :class:`set`) of one
              or more mappables, in which case this colorbar will be
              arbitrarily associated with the first mappable in this iterable.
            Defaults to `None`, in which case this parameter defaults to the
            iterable of all mappables provided by the topmost and hence last
            mappable layer in the current layer sequence (i.e., the
            :meth:`LayerCellsColorfulABC.color_mappables` property of the
            last instance of the :class:`LayerCellsColorfulABC` subclass in
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

        # If a single mappable rather than a sequence of mappables was passed,
        # convert the former to the latter.
        if isinstance(color_mappables, ScalarMappable):
            color_mappables = (color_mappables,)

        # Autoscale colors to the range implied by the passed color values.
        self._autoscale_colors(color_data)

        # Associate these mappables with this visuazilation's colorbar *AFTER*
        # preparing all layers defining these mappables if no mappables were
        # passed.
        self._automap_colors(color_mappables)


    @type_check
    def _autoscale_colors(self, color_data: SequenceTypes) -> None:
        '''
        Autoscale the colorbar for this visual's figure to the
        minimum and maximum scalar values unravelled from the passed sequence
        by setting the :attr:`_color_min` and :attr:`_color_max` attributes to
        such values if colorbar autoscaling is both enabled and has not already
        been performed (i.e., these attributes are ``None``) *or* noop
        otherwise.

        Parameters
        ----------
        color_data : SequenceOrNoneTypes
            Multi-dimensional sequence of all color values to be plotted *or*
            `None` if calculating these values on initialization is
            impractical. See the :meth:`_prep_figure` method for further
            details.
        '''

        # If colorbar autoscaling is disabled, noop.
        if not self._conf.is_color_autoscaled:
            return

        # Flatten this multi-dimensional array to a one-dimensional array,
        # permitting efficient retrieval of minimum and maximum values.
        time_series_flat = np.ravel(color_data)

        # Set the current minimum and maximum color values.
        self._color_min = np.ma.min(time_series_flat)
        self._color_max = np.ma.max(time_series_flat)

        # Log these values.
        logs.log_debug(
            'Autoscaling animation colors to [%d, %d].',
            self._color_min, self._color_max)


    @type_check
    def _automap_colors(self, color_mappables: IterableTypes) -> None:
        '''
        Create a figure colorbar for this visual associated with the
        passed mappables if any or the mappables provided by the last mappable
        layer in the current layer sequence otherwise.

        Parameters
        ----------
        color_mappables : optional[IterableTypes]
            Iterable of all :class:`ScalarMappable` instances (e.g.,
            :class:`AxesImage`, :class:`ContourSet`) to associate with this
            visual's colorbar. By Matplotlib design, only the first
            mappable in this iterable is arbitrarily associated with this
            colorbar; all other mappables are ignored. Defaults to `None`, in
            which case the iterable of all mappables provided by the topmost
            and hence last mappable layer in the current layer sequence is
            defaulted to (i.e., the value of the
            :meth:`LayerCellsColorfulABC.color_mappables` property of the last
            instance of the :class:`LayerCellsColorfulABC` subclass in the
            :attr:`_layers` attribute).
        '''

        # Classify this sequence of mappables.
        assert types.is_iterable_nonstr(color_mappables), (
            types.assert_not_iterable_nonstr(color_mappables))
        self._color_mappables = color_mappables

        self._rescale_color_mappables()

        # First mappable safely retrieved from this iterable of mappables.
        color_mappable_first = iterables.get_item_first(self._color_mappables)
        assert types.is_matplotlib_mappable(color_mappable_first), (
            types.assert_not_matplotlib_mappable(color_mappable_first))

        # Create a colorbar associated with this mappable.
        self.make_colorbar(color_mappable_first)


    def _rescale_color_mappables(self) -> None:

        # For each color mappable, clip that mappable to the minimum and
        # maximum values discovered above. Note this also has the beneficial
        # side-effect of establishing the colorbar's range.
        for color_mappable in self._color_mappables:
            assert types.is_matplotlib_mappable(color_mappable), (
                types.assert_not_matplotlib_mappable(color_mappable))
            color_mappable.set_clim(self._color_min, self._color_max)

    # ..................{ PLOTTERS                           }..................
    #FIXME: For generality, rename this method to visualize_time_step().
    #FIXME: For "PlotAfterSolving"-style plots, the first frame is uselessly
    #plotted twice. Investigate up this insanity, please. We might consider
    #simply ignoring the first call to this method if we can't get to the
    #ultimate bottom of this. Clearly, a matplotlib issue is implicated.

    @type_check
    def plot_frame(self, time_step: int) -> None:
        '''
        Display and/or save the frame corresponding to the passed sampled
        simulation time step for this visual.

        This method is iteratively called by matplotlib's
        :class:`FuncAnimation` class instantiated by our :meth:`_animate`
        method. The subclass implementation of this abstract method is expected
        to modify this animation's figure and axes so as to display the next
        frame. It is *not*, however, expected to save that figure to disk;
        frame saving is already implemented by this base class in a
        general-purpose manner.

        Specifically, this method (in order):

        . Calls the subclass :meth:`_plot_frame_figure` method to update the
          current figure with this frame's data.
        . Updates the current figure's axes title with the current time.
        . Optionally writes this frame to disk if desired.

        Parameters
        ----------
        time_step : int
            0-based index of the frame to be plotted _or_ -1 if the most recent
            frame is to be plotted.
        '''

        # Absolute 0-based index of the frame to be plotted.
        #
        # If the passed index is -1 and hence relative rather than absolute,
        # this index is assumed to be the last index of the current
        # simulation's array of time steps.
        if time_step == -1:
            time_step_absolute = len(self._phase.sim.time) - 1
        # Else, the passed index is already absolute and hence used as is.
        else:
            time_step_absolute = time_step

        # Log this animation frame.
        logs.log_debug(
            'Exporting "%s" frame %d / %d...',
            self._label, time_step_absolute, self._time_step_last)

        # Classify this time step for subsequent access by subclasses.
        self._time_step = time_step

        # Plot this frame's title *BEFORE* this frame, allowing axes changes
        # performed by the subclass implementation of the _plot_frame_figure()
        # method called below to override the default title.
        self._plot_frame_axes_title()

        # Plot this frame via external layers *BEFORE* plotting this frame
        # via subclass logic, allowing the latter to override the former.
        self._plot_layers()

        # Plot this frame via subclass logic *AFTER* performing all
        # superclass-specific plotting.
        self._plot_frame_figure()

        # Interactively display this frame if requested.
        self._show_frame(time_step_absolute)

        # Non-interactively write this frame to disk if requested.
        self._save_frame(time_step_absolute)

        # Inform subclasses that subsequent frames are no longer the first.
        self._is_time_step_first = False


    def _plot_frame_axes_title(self) -> None:
        '''
        Plot the current frame title of this animation onto this animation's
        axes.

        By default, this title interpolates the current time step and must thus
        be replotted for each animation frame.
        '''

        #FIXME: Shift into a new "betse.util.time.times" submodule.

        # Number of seconds in a minute.
        SECONDS_PER_MINUTE = 60

        # Number of seconds in an hour.
        SECONDS_PER_HOUR = SECONDS_PER_MINUTE * 60

        #FIXME: For efficiency, classify the following local variables as
        #object attributes in the __init__() method rather than recompute these
        #variables each animation frame here.

        # Human-readable suffix of the units that simulation times are reported
        # in (e.g., "ms" for milliseconds).
        time_unit_suffix = None

        # Factor by which low-level simulation times are multiplied to yield
        # human-readable simulation times in the same units of
        # "time_unit_suffix".
        time_unit_factor = None

        # Duration in seconds of the current simulation phase (e.g., "init",
        # "run"), accelerated by the current gap junction acceleration factor.
        time_len = self._phase.p.total_time_accelerated

        # If this phase runs for less than or equal to 100ms, report
        # simulation time in milliseconds (i.e., units of 0.001s).
        if time_len <= 0.1:
            time_unit_suffix = 'ms'
            time_unit_factor = 1e3
        # Else if this phase runs for less than or equal to one minute, report
        # simulation time in seconds (i.e., units of 1s).
        elif time_len <= SECONDS_PER_MINUTE:
            time_unit_suffix = 's'
            time_unit_factor = 1
        # Else if this phase runs for less than or equal to one hour, report
        # simulation time in minutes (i.e., units of 60s).
        elif time_len <= SECONDS_PER_HOUR:
            time_unit_suffix = ' minutes'
            time_unit_factor = 1/SECONDS_PER_MINUTE
        # Else, this phase is assumed to run for less than or equal to one day.
        # In this case, simulation time is reported in hours (i.e., units of
        # 60*60s).
        else:
            time_unit_suffix = ' hours'
            time_unit_factor = 1/SECONDS_PER_HOUR

        # Current time adjusted for long/short simulation.
        time_accelerated = (
            time_unit_factor *
            self._phase.sim.time[self._time_step])

        # Update this figure with this time, rounded to one decimal place.
        self._axes.set_title('{} (time: {:.1f}{})'.format(
            self._axes_title, time_accelerated, time_unit_suffix,))


    def _show_frame(self, time_step_absolute: int) -> None:
        '''
        Display the current frame of this animation if requested by the current
        simulation configuration *or* silently noop otherwise.

        Parameters
        ----------
        time_step_absolute : int
            0-based index of the current frame to be shown.
        '''

        # If *NOT* showing this frame, return immediately.
        if not self._is_show:
            return

        # Temporarily yield the time slice for the smallest amount of time
        # required by the current matplotlib backend to handle queued events in
        # the GUI-specific event loop of the current process. Failing to do so
        # produces subtle issues with plot and animation display. Technically,
        # doing so only appears to be required on one or more of the following
        # edge-case conditions:
        #
        # * The current platform is Windows, whose POSIX-incompatible process
        #   model has special needs.
        # * The current visual is non-blocking *AND* the current
        #   backend is poorly implemented or supported (e.g., "Qt4Agg").
        #
        # Since reliably detecting all such conditions is non-trivial and
        # certainly more time consuming than unconditionally yielding the time
        # slice for a negligible amount of time, we prefer the latter approach.
        #
        # In theory, the time slice should *NEVER* need to be manually yielded.
        # The "matplotlib.animation.FuncAnimation" class should handle such
        # low-level details, as confirmed by the pyplot.pause() docstring:
        #
        #  pause(interval)
        #
        #     Pause for *interval* seconds.
        #
        #     If there is an active figure it will be updated and displayed,
        #     and the gui event loop will run during the pause.
        #
        #     If there is no active figure, or if a non-interactive backend
        #     is in use, this executes time.sleep(interval).
        #
        #     This can be used for crude animation. For more complex
        #     animation, see :mod:`matplotlib.animation`.
        #
        # Note the last paragraph, implying the low-level pyplot.pause() method
        # to be an alternative (rather than addition) to the the higher-level
        # "matplotlib.animation" submodule. Unfortunately, this does *NOT*
        # appear to be the case; the pyplot.pause() method *MUST* be manually
        # called even when leveraging "matplotlib.animation" functionality.
        #
        # For further details, see also:
        #
        #     https://gitlab.com/betse/betse/issues/9
        #     https://github.com/matplotlib/matplotlib/issues/2134/
        #
        # When yielding the time slice, the following non-fatal deprecation
        # warning emitted by the pyplot.pause() function is temporarily
        # suppressed:
        #
        #     [betse] /usr/lib64/python3.4/site-packages/matplotlib/backend_bases.py:2437:
        #     MatplotlibDeprecationWarning: Using default event loop
        #     until function specific to this GUI is implemented
        #       warnings.warn(str, mplDeprecation)
        #
        # Since this warning is overly verbose and safely ignorable, we do so.
        # While mildly inefficient, ignoring this warning is preferable to the
        # alternatives of enraging, berating, and frightening end users.
        with mplutil.deprecations_ignored():
            pyplot.pause(0.0001)


    def _save_frame(self, time_step_absolute: int) -> None:
        '''
        Persist the current frame of this animation to disk if requested by the
        current simulation configuration *or* silently noop otherwise.

        Parameters
        ----------
        time_step_absolute : int
            0-based index of the current frame to be saved.
        '''

        # If *NOT* saving this frame, return immediately.
        if not self._is_save:
            return

        # If saving animation frames as images, save this frame as such.
        if self._writer_images is not None:
            self._writer_images.grab_frame(**self._writer_savefig_kwargs)

        # If saving animation frames as video, save this frame as such.
        if self._writer_video is not None:
            # For debuggability, temporarily escalate the
            # matplotlib-specific verbosity level.
            with mpl_config.reducing_log_level_to_debug_if_info():
                self._writer_video.grab_frame(
                    **self._writer_savefig_kwargs)

        # If this is the last frame to be plotted, finalize all writers
        # *AFTER* instructing these writers to write this frame.
        if time_step_absolute == self._time_step_last:
            logs.log_debug(
                'Finalizing animation "%s" saving...', self._label)
            self._close_writers()

    # ..................{ SUBCLASS                           }..................
    def _plot_frame_figure(self) -> None:
        '''
        Update this visual's figure (and typically axes) content to reflect the
        current simulation time step.

        This method defaults to a noop. Subclasses may optionally redefine this
        method with subclass-specific logic but are strongly encouraged to
        append layers containing such logic instead.
        '''

        pass

    # ..................{ PLOTTERS ~ image                   }..................
    #FIXME: Merge the useful docstring for this method into the
    #"betse.science.visual.layer.field.layerfieldsurface.LayerCellsFieldSurface"
    #class; then, excise this method.
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
            * Three-dimensional, an RGB image *not* mapped onto the passed
              colormap will be plotted.
            * Four-dimensional, an RGBa image *not* mapped onto the passed
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
            extent=self._phase.cache.upscaled.extent,
            cmap=colormap,
        )


    #FIXME: Replace entirely by the appropriate "LayerCellsFieldStream" subclass.
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
        magnitude_max: NumericSimpleTypes or NoneType
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
            grid_x = self._phase.cells.X * self._phase.p.um
        if grid_y is None:
            grid_y = self._phase.cells.Y * self._phase.p.um

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
            density=self._phase.p.stream_density,
            linewidth=(3.0*magnitude/magnitude_max) + 0.5,
            color=self._phase.p.vcolor,
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

        if self._phase.p.showCells is True:
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
        if self._phase.p.showCells:
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
        if self._phase.p.showCells:
            assert types.is_matplotlib_polycollection(cell_plot), (
                types.assert_not_matplotlib_polycollection(cell_plot))

            # Update this plot in-place.
            cell_plot.set_array(cell_data)
            cell_plot.set_verts(
                nparray.from_iterable(self._phase.cells.cell_verts) *
                self._phase.p.um)

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
            verts=expmath.upscale_coordinates(self._phase.cells.cell_verts),
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
    #"LayerCellsVectorSurfaceDiscrete" subclass. Generalize this method into a
    #new "LayerCellsGouraudContinuous" subclass of the same submodule.

    # FIXME: Your Layers API is totally incomprehensible. I can't figure out wtf is going on!
    # FIXME: I like this version better than the current one that animates things like Vmem when p.showCells is False
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
        if len(cell_data) == len(self._phase.cells.mem_i):
            # cell_data = np.dot(
            #     self._phase.cells.M_sum_mems, cell_data) / self._phase.cells.num_mems
            xi = self._phase.p.um*self._phase.cells.mem_mids_flat[:,0]
            yi = self._phase.p.um*self._phase.cells.mem_mids_flat[:, 1]

        else:

            xi = self._phase.p.um*self._phase.cells.cell_centres[:,0]
            yi = self._phase.p.um*self._phase.cells.cell_centres[:, 1]

        # Unstructured triangular grid assigned the passed cell data.
        # triangular_grid = np.zeros(len(self._phase.cells.voronoi_centres))
        # triangular_grid[self._phase.cells.cell_to_grid] = cell_data

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
            xi,
            yi,
            cell_data,
            shading='gouraud',
            cmap=self._colormap
        )
