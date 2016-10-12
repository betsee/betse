#!/usr/bin/env python3
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Abstract base classes of all Matplotlib-based plotting classes.
'''

#FIXME: Refactor all procedural cell cluster-specific "betse.plot.plot"
#functions into subclasses of the "PlotCellsABC" base class defined below.
#Ultimate power fights the dark deceit!

# ....................{ IMPORTS                            }....................
import numpy as np
import weakref
from abc import ABCMeta  #, abstractmethod  #, abstractstaticmethod
from betse.exceptions import BetseMethodException
from betse.lib.matplotlib.matplotlibs import ZORDER_STREAM
# from betse.util.io.log import logs
from betse.util.type import types, objects
from betse.util.type.types import (
    type_check, NoneType,
    NumericTypes, NumericOrNoneTypes,
    SequenceTypes, SequenceOrNoneTypes,
)
from matplotlib import pyplot
from matplotlib.collections import PolyCollection
from matplotlib.colors import Colormap
from matplotlib.image import AxesImage
from matplotlib.patches import FancyArrowPatch
from matplotlib.streamplot import StreamplotSet

# ....................{ BASE                               }....................
class PlotCellsABC(object, metaclass=ABCMeta):
    '''
    Abstract base class of all classes spatially plotting the cell cluster.

    Subclasses of this class plot the spatial distribution of one or more
    modelled variables (e.g., charge distribution, membrane voltage) for either:

    * A single simulation time step, in which this subclass plots a still frame
      of these variables for this step.
    * One or more simulation time steps, in which this subclass animates a video
      of these variables for these steps.

    Attributes (Public)
    ----------
    cells : Cells
        Current cell cluster.
    p : Parameters
        Current simulation configuration.
    sim : Simulation
        Current simulation.

    Attributes (Private)
    ----------
    _label : str
        Basename of the subdirectory in the phase-specific results directory
        to which all files exported for this plot or animation are saved _and_
        the basename prefix of these files.

    Attributes (Private: Figure)
    ----------
    _figure : Figure
        Matplotlib figure providing the current animation frame.
    _figure_title : str
        Text displayed above the figure itself.

    Attributes (Private: Axes)
    ----------
    _axes : FigureAxes
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
    _color_mappings : list
        List of all plotted mappables (i.e., instances of the `ScalarMappable`
        Matplotlib class added to this plot) whose minimum and maximum values
        define the range of colors displayed by this plot. Due to Matplotlib
        constraints, only the first mappable in this list will be associated
        with this plot's colorbar.
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
        # scaling_series: np.nadarray = None,
    ) -> None:
        '''
        Initialize this plot.

        Parameters
        ----------
        sim : Simulator
            Current simulation.
        cells : Cells
            Current cell cluster.
        p : Parameters
            Current simulation configuration.
        label : str
            Basename of the subdirectory in the phase-specific results directory
            to which all animation files will be saved _and_ the basename prefix
            of these files.
        figure_title : str
            Text displayed above the figure itself.
        colorbar_title: str
            Text displayed above the figure colorbar.
        axes_title : str
            Optional text displayed above the figure axes but below the figure
            title _or_ `None` if no text is to be displayed. Defaults to `None`.
        axes_x_label : str
            Text displayed below the figure's X axis. Defaults to a general-
            purpose string.
        axes_y_label : str
            Text displayed to the left of the figure's Y axis. Defaults to a
            general-purpose string.
        is_color_autoscaled : bool
            `True` if dynamically resetting the minimum and maximum colorbar
            values to be the corresponding minimum and maximum values for the
            current frame _or_ `False` if statically setting the minimum and
            maximum colorbar values to predetermined constants.
        color_min : NumericTypes
            Minimum colorbar value to be used if `clrAutoscale` is `False`.
        color_max : NumericTypes
            Maximum colorbar value to be used if `clrAutoscale` is `False`.
        colormap : matplotlib.cm.Colormap
            Matplotlib colormap to be used by default for all animation artists
            (e.g., colorbar, images) _or_ `None`, in which case the default
            colormap will be used.
        '''

        # Classify core parameters with weak rather than strong (the default)
        # references, thus avoiding circular references and the resulting
        # complications thereof (e.g., increased memory overhead). Since these
        # objects necessarily live significantly longer than this plot, no
        # complications arise. Ergo, these attributes *ALWAYS* yield these
        # objects rather than non-deterministically returning "None" if these
        # objects are unexpectedly garbage-collected.
        self.sim = weakref.proxy(sim)
        self.p = weakref.proxy(p)

        #FIXME: For currently unknown reasons, this object occasionally retains
        #the only remaining reference to the passed "Cells" instance -- which
        #*CANNOT* therefore be safely classified as a weak reference. This is
        #highly unexpected, however, and should thus be investigated.
        # self.cells = weakref.proxy(cells)
        self.cells = cells

        # Default unpassed parameters.
        if colormap is None:
            colormap = p.default_cm

        # Classify *AFTER* validating parameters.
        self._label = label
        self._figure_title = figure_title
        self._axes_title = axes_title
        self._axes_x_label = axes_x_label
        self._axes_y_label = axes_y_label
        self._colorbar_title = colorbar_title
        self._is_color_autoscaled = is_color_autoscaled
        self._color_min = color_min
        self._color_max = color_max
        self._colormap = colormap
        # self._scaling_series = scaling_series

        # Classify attributes to be subsequently defined.
        self._color_mappings = color_min
        self._writer_frames = None
        self._writer_video = None

        # Initialize this plot's figure.
        self._init_figure()

    # ..................{ INITIALIZERS ~ figure              }..................
    def _init_figure(self) -> None:
        '''
        Initialize this plot's figure.
        '''

        # Figure encapsulating this animation as a weak rather than strong (the
        # default) reference, avoiding circular references and complications
        # thereof (e.g., memory overhead). Figures created by the "pyplot" API
        # are internally retained in Matplotlib's "Gcf" figure cache until
        # explicitly closed -- either non-interactively by a close() call or
        # interactively by the corresponding GUI window being closed. Hence,
        # strong figure references should typically *NOT* be retained.
        self._figure = weakref.proxy(pyplot.figure())

        # Figure axes scaled to the extent of the current 2D environment as a
        # weak rather than strong (the default) reference, thus avoiding
        # circular references and complications thereof (e.g., memory overhead).
        # Since figures already contain their axes as a strong reference, we
        # need *NOT* do so as well here.
        self._axes = weakref.proxy(pyplot.subplot(111))

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
            self.cells.xmin * self.p.um,
            self.cells.xmax * self.p.um,
            self.cells.ymin * self.p.um,
            self.cells.ymax * self.p.um,
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

        # Mandatory parameters.
        color_mapping: object,

        # Optional parameters.
        color_data: SequenceTypes = None,
    ) -> None:
        '''
        Prepare this plot _after_ having previously initialized this plot but
        _before_ subsequently displaying and/or saving this plot.

        Specifically, this method:

        * If the current simulation configuration requests that each plotted
          cell by labelled with that cell's unique 0-based index, do so. To
          ensure that these labels are plotted over rather than under the
          contents of their corresponding cells, we do so only _after_ all
          subclass plotting has been performed by delaying such labelling to
          this method rather than the above `__init__()` method.
        * If the optional `axes_title` parameter was passed to `__init__()`:
          * Adds the current `_figure_title` to this figure as a "super title."
          * Adds that `axes_title` to this figure's axes as a "subtitle."
        * Else, add the current `_figure_title` to this figure's axes.
        * If the optional `colorbar_values` parameter is passed _and_
          `_is_color_autoscaled` is `True`, clip the colorbar to the minimum and
          maximum values in the `color_data` array.
        * Else, clip the colorbar to the current `clrMin` and `clrMax` values.
        * Add a colorbar whose:
          * Title is the current `_colorbar_title` string.
          * Mapping is the passed `_colorbar_mapping` parameter.

        Parameters
        ----------
        color_mapping : mpl.cm.ScalarMappable, list
            One or more color mappables (e.g., `Image`, `ContourSet`) to
            associate with this animation's colorbar. This may be either:
            * A single mappable, in which case this colorbar will be associated
              with this mappable as is.
            * A non-string sequence (e.g., `list`) of one or more mappables, in
              which case this colorbar will be associated with the **first**
              mappable in this sequence.
        color_data : optional[SequenceTypes]
            Optional multi-dimensional Numpy array containing all data values
            to be animated _or_ `None` if calculating this data during the
            animation initialization is infeasible or impractical (e.g., due to
            space and time constraints). If the colorbar is being autoscaled
            _and_ this parameter is:
            * Non-`None`, the colorbar will be clipped to the minimum and
              maximum scalar values unravelled from this array.
            * `None`, the subclass will be responsible for colorbar autoscaling.
            Defaults to `None`.
        '''

        #FIXME: Why is this logic deferred to here rather than performed in
        #__init__()? There's presumably a reasonable reason for this... but we
        #failed to document whatever that might have been.

        # If labelling each plotted cell with that cell's unique 0-based index,
        # do so.
        if self.p.enumerate_cells is True:
            for cell_index, cell_centre in enumerate(self.cells.cell_centres):
                self._axes.text(
                    self.p.um * cell_centre[0],
                    self.p.um * cell_centre[1],
                    cell_index,
                    va='center',
                    ha='center',
                )

        # If a time series is passed *AND* colorbar autoscaling is requested,
        # clip the colorbar to the minimum and maximum values of this series.
        if color_data is not None and self._is_color_autoscaled:
            #FIXME: This appears to be failing when "color_data" is
            #"self._current_density_magnitude_time_series" (e.g.,
            #"self.sim.I_gj_x_time").

            # Flatten this two-dimensional matrix to a one-dimensional array,
            # providing efficient retrieval of minimum and maximum values.
            time_series_flat = np.ravel(color_data)

            # Overwrite the current minimum and maximum color values.
            self._color_min = np.ma.min(time_series_flat)
            self._color_max = np.ma.max(time_series_flat)

        # If a single mappable rather than a list of mappables was passed,
        # convert the former to the latter.
        if not types.is_sequence_nonstr(color_mapping):
            color_mapping = (color_mapping,)

        # Classify this list of mappables.
        self._color_mappings = color_mapping

        # Scale these mappables by previously established color values.
        self._rescale_colors()

        # Create a colorbar associated with the first such mappable.
        colorbar = self._figure.colorbar(color_mapping[0])
        colorbar.set_label(self._colorbar_title)

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
        for field_name, field_value in objects.iter_fields_nonbuiltin(self):
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

    # ..................{ PROPERTIES                         }..................
    #FIXME: Convert all of the following properties back into normal instance
    #variables initialized from parameters passed to this class' __init__()
    #method. Doing so requires subclass __init__() methods to explicitly accept
    #the "p" parameter, permitting such methods to initialize these paremeters
    #in the expected way. Doing so is trivial. The current approach is both
    #obscene design overkill *AND* inefficient. This saddens us.

    # The following testers are intended to be overridden by subclasses.
    #
    # The corresponding attributes (e.g., "_is_showing" for _is_showing())
    # *MUST* be defined as dynamic methods rather than static attributes passed
    # to this class' __init__() method (e.g., as an "is_saving" parameter.) Why?
    # Because chicken-and-the-egg constraints. Specifically, the latter approach
    # prevents subclasses from passing a value dependent on the current
    # "Parameters" object to __init__(), as that object has yet to be classified
    # as the "p" attribute yet. (Ugh.)

    @property
    def _is_showing(self) -> bool:
        '''
        `True` only if interactively displaying this plot or animation.

        This property is orthogonal to the `_is_saving` property, whose value
        may concurrently also be `True`.
        '''

        return not self.p.turn_all_plots_off


    @property
    def _is_saving(self) -> bool:
        '''
        `True` only if non-interactively saving this plot or animation.

        This property is orthogonal to the `_is_showing` property, whose value
        may concurrently also be `True`.
        '''

        return self.p.autosave

    # ..................{ COLORS                             }..................
    def _rescale_colors(self):
        '''
        Rescale all color mappables (i.e., the `color_mapping` parameter passed
        to the `_prep_figure()` method) to the current minimum and maximum color
        values for this plot.

        This method may only be called _after_ the `_prep_figure()` method has
        been called. Failing to do so will result in exceptions.
        '''

        # If the _prep_figure() method has yet to be called, raise an exception.
        if self._color_mappings is None:
            raise BetseMethodException(
                '{class_name}._rescale_colors() called before '
                '{class_name}._prep_figure().'.format(class_name=type(self)))

        # If these values are identical, coerce them to differ. Failing to do so
        # produces spurious visual artifacts in both the axes and colorbar.
        # if self._color_min == self._color_max:
        #     self._color_min = self._color_min - 1
        #     self._color_max = self._color_max + 1

        # For each color mappable, clip that mappable to the minimum and maximum
        # values discovered above. Note this also has the beneficial side-effect
        # of establishing the colorbar's range.
        for color_mapping in self._color_mappings:
            assert types.is_matplotlib_mappable(color_mapping), (
                types.assert_not_matplotlib_mappable(color_mapping))
            color_mapping.set_clim(self._color_min, self._color_max)

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
            Optional maximum magnitude in the passed `magnitude` array. Defaults
            to `None`, in which case this array is searched for this value.
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
            grid_x = self.cells.X * self.p.um
        if grid_y is None:
            grid_y = self.cells.Y * self.p.um

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
            density=self.p.stream_density,
            linewidth=(3.0*magnitude/magnitude_max) + 0.5,
            color=self.p.vcolor,
            cmap=self._colormap,
            arrowsize=1.5,

            # Draw this streamplot over all patch and line artists, by default.
            # See the "ZORDER_STREAM" docstring for further commentary.
            zorder=ZORDER_STREAM,
        )

    # ..................{ PLOTTERS ~ cell                    }..................
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

        if self.p.showCells is True:
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
        if self.p.showCells:
            assert types.is_matplotlib_polycollection(cell_plot), (
                types.assert_not_matplotlib_polycollection(cell_plot))
            cell_plot.set_array(cell_data)
        # Else, the cell cluster is being plotted as an unstructured triangular
        # grid, requiring this data be reshaped onto this grid.
        else:
            assert types.is_matplotlib_trimesh(cell_plot), (
                types.assert_not_matplotlib_trimesh(cell_plot))

            # cell_data = np.zeros(len(self.cells.voronoi_centres))
            # cell_data[self.cells.cell_to_grid] = cell_data

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
        if self.p.showCells is True:
            assert types.is_matplotlib_polycollection(cell_plot), (
                types.assert_not_matplotlib_polycollection(cell_plot))

            # Update this plot in-place.
            cell_plot.set_array(cell_data)
            cell_plot.set_verts(np.asarray(self.cells.cell_verts) * self.p.um)

            # Return the same plot.
            return cell_plot
        # Else, the passed cell plot *MUST* be a triangle mesh previously
        # returned by the _plot_cell_mesh() method.
        else:
            assert types.is_matplotlib_trimesh(cell_plot), (
                types.assert_not_matplotlib_trimesh(cell_plot))


            cell_plot.remove()
            return self._plot_cell_mesh(cell_data=cell_data, *args, **kwargs)


    def _plot_cell_mosaic(self, cell_data: np.ndarray) -> PolyCollection:
        '''
        Plot and return a mosaic plot of all cells with colours corresponding to
        the passed vector of arbitrary cell data (e.g., transmembrane voltages
        for all cells for the current time step) onto the current figure's axes.

        The returned plot will be a polygon collection such that each polygon
        signifies a cell in this simulation's cell cluster.

        Parameters
        -----------
        cell_data : np.ndarray
            Arbitrary cell data defined on an environmental grid to be plotted.

        Returns
        --------
        PolyCollection
            Mosaic plot produced by plotting the passed cell data.
        '''
        assert types.is_sequence_nonstr(cell_data), (
            types.assert_not_sequence_nonstr(cell_data))

        # Cell vertices plotted as polygons.
        mosaic_plot = PolyCollection(
            verts=np.asarray(self.cells.cell_verts) * self.p.um,
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
        if len(cell_data) == len(self.cells.mem_i):
            cell_data = np.dot(
                self.cells.M_sum_mems, cell_data) / self.cells.num_mems

        # Unstructured triangular grid assigned the passed cell data.
        triangular_grid = np.zeros(len(self.cells.voronoi_centres))
        triangular_grid[self.cells.cell_to_grid] = cell_data

        # cmin = cell_data.min()
        # cmax = cell_data.max()

        # Create and add this plot to this figure's axes and return this plot.
        # Unfortunately, the tripcolor() function requires:
        #
        # * The "x", "y", and "triangles" parameters to be positional.
        # * All other parameters to be keyword.
        #
        # Behold! The ultimate examplar of nonsensical API design.
        return self._axes.tripcolor(
            self.p.um * self.cells.cell_centres[:, 0],
            self.p.um * self.cells.cell_centres[:, 1],
            cell_data,
            shading='gouraud',
            cmap=self._colormap
        )
