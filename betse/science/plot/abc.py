#!/usr/bin/env python3
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Abstract base classes of all Matplotlib-based plotting classes.
'''

#FIXME: Refactor all procedural cell cluster-specific "betse.plot.plot"
#functions into subclasses of the "PlotCells" base class defined below.
#Ultimate power fights the dark deceit!

# ....................{ IMPORTS                            }....................
import gc, weakref
import numpy as np
from abc import ABCMeta  #, abstractmethod  #, abstractstaticmethod
from betse.exceptions import BetseExceptionMethod
from betse.lib.matplotlib.matplotlibs import ZORDER_STREAM
from betse.util.py import objects
from betse.util.type import types
from matplotlib import pyplot
from matplotlib.collections import PolyCollection
from matplotlib.patches import FancyArrowPatch

# ....................{ BASE                               }....................
class PlotCells(object, metaclass=ABCMeta):
    '''
    Abstract base class of all classes spatially plotting the cell cluster.

    Instances of this class plot the spatial distribution of modelled variables
    (e.g., Vmem) over some time step(s) of the simulation.

    Attributes
    ----------
    _sim : Simulation
        Current simulation.
    _cells : Cells
        Current cell cluster.
    _p : Parameters
        Current simulation configuration.
    _axes : FigureAxes
        Matplotlib figure axes providing the current animation frame data.
    _axes_bounds : list
        Spacial extent of the current 2D environment as a 4-element list
        conisting of (in order):
        1. Minimum value of the figure's X axis.
        2. Maximum value of the figure's X axis.
        3. Minimum value of the figure's Y axis.
        4. Maximum value of the figure's Y axis.
    _axes_title : str
        Text displayed above the figure axes. If a non-`None` value for the
        `axes_title` parameter is passed to the `__init__()` method, this is
        that value; else, this is the value of the `figure_title` parameter
        passed to the same method.
    _axes_x_label : str
        Text displayed below the figure's X axis.
    _axes_y_label : str
        Text displayed to the left of the figure's Y axis.
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
    _figure : Figure
        Matplotlib figure providing the current animation frame.
    _figure_title : str
        Text displayed above the figure itself.
    _is_color_autoscaled : bool
        `True` if dynamically resetting the minimum and maximum colorbar values
        to be the corresponding minimum and maximum values for the current
        frame _or_ `False` if statically setting the minimum and maximum
        colorbar values to predetermined constants.
    _type : str
        Basename of the subdirectory in the phase-specific results directory
        to which all animation files will be saved _and_ the basename prefix of
        these files.
    '''

    # ..................{ LIFECYCLE                          }..................
    def __init__(
        self,

        # Mandatory parameters.
        sim: 'Simulator',
        cells: 'Cells',
        p: 'Parameters',
        type: str,
        figure_title: str,
        colorbar_title: str,
        is_color_autoscaled: bool,
        color_min: float,
        color_max: float,

        # Optional parameters.
        axes_title: str = None,
        axes_x_label: str = 'Spatial Distance [um]',
        axes_y_label: str = 'Spatial Distance [um]',
        colormap: 'Colormap' = None,
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
        type : str
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
        color_min : float
            Minimum colorbar value to be used if `clrAutoscale` is `False`.
        color_max : float
            Maximum colorbar value to be used if `clrAutoscale` is `False`.
        colormap : matplotlib.cm.Colormap
            Matplotlib colormap to be used by default for all animation artists
            (e.g., colorbar, images) _or_ `None`, in which case the default
            colormap will be used.
        '''
        # Validate core parameters.
        assert types.is_simulator(sim), types.assert_not_simulator(sim)
        assert types.is_cells(cells), types.assert_not_parameters(cells)
        assert types.is_parameters(p), types.assert_not_parameters(p)

        # Classify core parameters with weak rather than strong (the default)
        # references, thus avoiding circular references and all resulting
        # complications thereof (e.g., increased memory overhead). Since these
        # objects necessarily live significantly longer than this plot, no
        # complications arise. These attributes *ALWAYS* provide the expected
        # objects rather than non-deterministically returning "None".
        self._sim = weakref.proxy(sim)
        self._cells = weakref.proxy(cells)
        self._p = weakref.proxy(p)

        # Default unpassed parameters.
        if colormap is None:
            colormap = p.default_cm

        # Validate all remaining parameters *AFTER* defaulting parameters. Note
        # that the "axes_title" parameter is subsequently validated by the
        # _animate() method.
        assert types.is_str_nonempty(type), (
            types.assert_not_str_nonempty(type, 'Animation type'))
        assert types.is_str_nonempty(figure_title), (
            types.assert_not_str_nonempty(figure_title, 'Figure title'))
        assert types.is_str_nonempty(colorbar_title), (
            types.assert_not_str_nonempty(colorbar_title, 'Colorbar title'))
        assert types.is_str_nonempty(axes_x_label), (
            types.assert_not_str_nonempty(axes_x_label, 'X axis label'))
        assert types.is_str_nonempty(axes_y_label), (
            types.assert_not_str_nonempty(axes_y_label, 'Y axis label'))
        assert types.is_bool(is_color_autoscaled), types.assert_not_bool(is_color_autoscaled)
        assert types.is_numeric(color_min), types.assert_not_numeric(color_min)
        assert types.is_numeric(color_max), types.assert_not_numeric(color_max)
        assert types.is_matplotlib_colormap(colormap), (
            types.assert_not_matplotlib_colormap(colormap))

        # Classify *AFTER* validating parameters.
        self._type = type
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

        # Figure encapsulating this animation as a weak rather than strong (the
        # default) reference, avoiding circular references and complications
        # thereof (e.g., memory overhead). Figures created by the "pyplot" API
        # are internally retained in Matplotlib's "Gcf" figure cache until
        # explicitly closed -- either non-interactively by a close() call or
        # interactively by the corresponding GUI window being closed. Hence,
        # strong figure references should typically *NOT* be retained.
        self._figure = weakref.proxy(pyplot.figure())

        # Extent of the current 2D environment.
        self._axes_bounds = [
            self._cells.xmin * self._p.um,
            self._cells.xmax * self._p.um,
            self._cells.ymin * self._p.um,
            self._cells.ymax * self._p.um,
        ]

        # Figure axes scaled to the extent of the current 2D environment as a
        # weak rather than strong (the default) reference, thus avoiding
        # circular references and complications thereof (e.g., memory overhead).
        # Since figures contain their axes as a strong reference, we needn't.
        self._axes = weakref.proxy(pyplot.subplot(111))
        self._axes.axis('equal')
        self._axes.axis(self._axes_bounds)
        self._axes.set_xlabel(self._axes_x_label)
        self._axes.set_ylabel(self._axes_y_label)


    def _prep_figure(
        self,

        # Mandatory parameters.
        color_mapping: object,

        #FIXME: Rename to merely "color_data".
        #FIXME: Make this parameter mandatory. There's no justifiable reason for
        #subclasses to omit this anymore.

        # Optional parameters.
        color_series: np.ndarray = None,
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
          * Adds the passed `axes_title` to this figure's axes as a "subtitle."
        * Else, add the current `_figure_title` to this figure's axes.
        * If the optional `colorbar_values` parameter is passed _and_
          `_is_color_autoscaled` is `True`, clip the colorbar to the minimum and
          maximum values in the `color_series` array.
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
        color_series : optional[np.ndarray]
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

        # If labelling each plotted cell with that cell's unique 0-based index,
        # do so.
        if self._p.enumerate_cells is True:
            for cell_index, cell_centre in enumerate(self._cells.cell_centres):
                self._axes.text(
                    self._p.um * cell_centre[0],
                    self._p.um * cell_centre[1],
                    cell_index,
                    va='center',
                    ha='center',
                )

        # If both a figure and axes title are defined, display the figure title
        # as such above the axes title.
        if self._axes_title:
            self._figure.suptitle(
                self._figure_title, fontsize=14, fontweight='bold')
        # Else, display the figure title as the axes title.
        else:
            self._axes_title = self._figure_title

        # Add the desired axes title.
        assert types.is_str_nonempty(self._axes_title), (
            types.assert_not_str_nonempty(self._axes_title, 'Axis title'))
        self._axes.set_title(self._axes_title)

        # If a time series is passed *AND* colorbar autoscaling is requested,
        # clip the colorbar to the minimum and maximum values of this series.
        if color_series is not None and self._is_color_autoscaled:
            assert types.is_sequence_nonstr(color_series), (
                types.assert_not_sequence_nonstr(color_series))

            #FIXME: This appears to be failing when "color_series" is
            #"self._current_density_magnitude_time_series" (e.g.,
            #"self._sim.I_gj_x_time").

            # Flatten this two-dimensional matrix to a one-dimensional array,
            # providing efficient retrieval of minimum and maximum values.
            time_series_flat = np.ravel(color_series)

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


    def close(self) -> None:
        '''
        Destroy this plot and deallocate all memory associated with this plot.

        To reduce Matplotlib's memory overhead, this method (in order):

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
            raise BetseExceptionMethod(
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
    def _plot_image(
        self,

        # Mandatory parameters.
        pixel_data: np.ndarray,

        # Optional parameters.
        colormap : 'matplotlib.cm.Colormap' = None,
    ) -> 'matplotlib.image.AxesImage':
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
        assert types.is_sequence_nonstr(pixel_data), (
            types.assert_not_sequence_nonstr(pixel_data))

        # Default unpassed parameters.
        if colormap is None:
            colormap = self._colormap
        assert types.is_matplotlib_colormap(colormap), (
            types.assert_not_matplotlib_colormap(colormap))

        # Plot and return this image.
        return self._axes.imshow(
            pixel_data,
            origin='lower',
            extent=self._axes_bounds,
            cmap=colormap,
        )


    def _plot_stream(
        self,

        # Mandatory arguments.
        x: np.ndarray,
        y: np.ndarray,
        magnitude: np.ndarray,

        # Optional arguments.
        grid_x: np.ndarray = None,
        grid_y: np.ndarray = None,
        magnitude_max: float = None,
        old_stream_plot: 'matplotlib.streamplot.StreamplotSet' = None,
    ) -> 'matplotlib.streamplot.StreamplotSet':
        '''
        Plot and return a streamplot of all streamlines of the passed vector
        flow data onto the current figure's axes.

        Parameters
        ----------
        x : np.ndarray
            Two-dimensional X components of the vector flow velocity.
        y : np.ndarray
            Two-dimensional Y components of the vector flow velocity.
        magnitude: np.ndarray
            One-dimensional vector flow magnitudes.
        grid_x : np.ndarray
            Optional scaled X components of the cell cluster grid. Defaults to
            `None`, in which case the default `self._cells.X` array is used.
        grid_y : np.ndarray
            Optional scaled Y components of the cell cluster grid. Defaults to
            `None`, in which case the default `self._cells.Y` array is used.
        magnitude_max: float
            Optional maximum magnitude in the passed `magnitude` array. Defaults
            to `None`, in which case this array is searched for this value.
        old_stream_plot: StreamplotSet
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
        assert types.is_sequence_nonstr(x), (
            types.assert_not_sequence_nonstr(x))
        assert types.is_sequence_nonstr(y), (
            types.assert_not_sequence_nonstr(y))
        assert types.is_sequence_nonstr(magnitude), (
            types.assert_not_sequence_nonstr(magnitude))

        # Default all unpassed optional arguments.
        if magnitude_max is None:
            magnitude_max = np.max(magnitude)
        if grid_x is None:
            grid_x = self._cells.X * self._p.um
        if grid_y is None:
            grid_y = self._cells.Y * self._p.um
        assert types.is_numeric(magnitude_max), (
            types.assert_not_numeric(magnitude_max))
        assert types.is_sequence_nonstr(grid_x), (
            types.assert_not_sequence_nonstr(grid_x))
        assert types.is_sequence_nonstr(grid_y), (
            types.assert_not_sequence_nonstr(grid_y))

        # If a prior streamplot to be erased was passed, do so.
        if old_stream_plot is not None:
            assert types.is_matplotlib_streamplot(old_stream_plot), (
                types.assert_not_matplotlib_streamplot(old_stream_plot))

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

            #FIXME: Duplicated from below.
            # Reshape this data onto this grid.
            cell_data = np.zeros(len(self._cells.voronoi_centres))
            cell_data[self._cells.cell_to_grid] = cell_data

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
            cell_plot.set_verts(np.asarray(self._cells.cell_verts) * self._p.um)

            # Return the same plot.
            return cell_plot
        # Else, the passed cell plot *MUST* be a triangle mesh previously
        # returned by the _plot_cell_mesh() method.
        else:
            assert types.is_matplotlib_trimesh(cell_plot), (
                types.assert_not_matplotlib_trimesh(cell_plot))

            #FIXME: We're fairly certain that this isn't actually necessary and
            #that, consequently, the following two expensive statements are
            #equivalent to this noop:
            #
            #    return cell_plot
            #
            #The reason why is that the _plot_cell_mesh() method triangulates
            #from the Voronoi centres rather than edges of each cell -- and the
            #centres don *NOT* appear to actually ever change. Hence, the same
            #exact plot with possibly different cell data (which is trivially
            #and efficiently settable elsewhere by calling
            #"cell_plot.set_array(cell_data)") is returned. Tasty fudge sundaes!

            # Remove this plot and create and return a new plot. Sadly, triangle
            # meshes do *NOT* currently support in-place update. Note that we
            # could technically break privacy encapsulation to update this plot
            # in-place with the following code:
            #
            #     # ...where "triangular_grid" is as defined locally by the
            #     # _plot_cell_mesh() method.
            #     cell_plot._triangulation = triangular_grid
            #     cell_plot._paths = None
            #
            # Doing so subverts the "TriMesh" API, however, and hence is likely
            # to break on Matplotlib updates. While inefficient, the current
            # approach is vastly more robust.
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
            verts=np.asarray(self._cells.cell_verts) * self._p.um,
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


    #FIXME: Replace all calls to cell_mesh() with calls this method; then,
    #excise cell_mesh() entirely.
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

        # Create and add this plot to this figure's axes and return this plot.
        # Unfortunately, the tripcolor() function requires:
        #
        # * The "x", "y", and "triangles" parameters to be positional.
        # * All other parameters to be keyword.
        #
        # Behold! The ultimate examplar of nonsensical API design.
        return self._axes.tripcolor(
            self._p.um*self._cells.voronoi_centres[:,0],
            self._p.um*self._cells.voronoi_centres[:,1],
            triangular_grid,
            shading='gouraud',
            cmap=self._colormap,
        )
