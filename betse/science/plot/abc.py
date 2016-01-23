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
import numpy as np
from abc import ABCMeta  #, abstractmethod  #, abstractstaticmethod
from betse.util.type import types
from matplotlib import pyplot
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
    _colorbar_title: str
        Text displayed above the figure colorbar.
    _is_color_autoscaled : bool
        `True` if dynamically resetting the minimum and maximum colorbar values
        to be the corresponding minimum and maximum values for the current
        frame _or_ `False` if statically setting the minimum and maximum
        colorbar values to predetermined constants.
    _color_min : float
        Minimum colorbar value to be used. If `clrAutoscale` is `True`, the
        subclass is responsible for redefining this value as appropriate.
    _color_max : float
        Maximum colorbar value to be used. If `clrAutoscale` is `True`, the
        subclass is responsible for redefining this value as appropriate.
    _colormap : Colormap
        Matplotlib colormap with which to create this animation's colorbar.
    _figure : Figure
        Matplotlib figure providing the current animation frame.
    _figure_title : str
        Text displayed above the figure itself.
    _type : str
        Basename of the subdirectory in the phase-specific results directory
        to which all animation files will be saved _and_ the basename prefix of
        these files.
    '''

    # ..................{ PRIVATE ~ init                     }..................
    def __init__(
        self,

        # Mandatory parameters.
        sim: 'Simulator',
        cells: 'Cells',
        p: 'Parameters',
        type: str,
        figure_title: str,
        colorbar_title: str,
        axes_x_label: str,
        axes_y_label: str,
        is_color_autoscaled: bool,
        color_min: float,
        color_max: float,

        # Optional parameters.
        axes_title: str = None,
        colormap: 'Colormap' = None,
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
            Text displayed below the figure's X axis.
        axes_y_label : str
            Text displayed to the left of the figure's Y axis.
        is_color_autoscaled : bool
            `True` if dynamically resetting the minimum and maximum colorbar
            values to be the corresponding minimum and maximum values for the
            current frame _or_ `False` if statically setting the minimum and
            maximum colorbar values to predetermined constants.
        color_min : float
            Minimum colorbar value to be used if `clrAutoscale` is `False`.
        color_max : float
            Maximum colorbar value to be used if `clrAutoscale` is `False`.
        colormap : Colormap
            Matplotlib colormap to be used in this animation's colorbar or
            `None`, in which case the default colormap will be used.
        '''
        # Validate core parameters.
        assert types.is_simulator(sim), types.assert_not_simulator(sim)
        assert types.is_cells(cells), types.assert_not_parameters(cells)
        assert types.is_parameters(p), types.assert_not_parameters(p)

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
        self._sim = sim
        self._cells = cells
        self._p = p
        self._type = type
        self._figure_title = figure_title
        self._axes_title = axes_title
        self._colorbar_title = colorbar_title
        self._is_color_autoscaled = is_color_autoscaled
        self._color_min = color_min
        self._color_max = color_max
        self._colormap = colormap

        # Classify attributes to be subsequently defined.
        self._writer_frames = None
        self._writer_video = None

        # Figure encapsulating this animation.
        self._figure = pyplot.figure()

        # Extent of the current 2D environment.
        self._axes_bounds = [
            self._cells.xmin * self._p.um,
            self._cells.xmax * self._p.um,
            self._cells.ymin * self._p.um,
            self._cells.ymax * self._p.um,
        ]

        # Figure axes scaled to the extent of the current 2D environment.
        self._axes = pyplot.subplot(111)
        self._axes.axis('equal')
        self._axes.axis(self._axes_bounds)
        self._axes.set_xlabel(axes_x_label)
        self._axes.set_ylabel(axes_y_label)

    # ..................{ PRIVATE ~ plot                     }..................
    def _prep(
        self,

        # Mandatory parameters.
        color_mapping: object,

        # Optional parameters.
        color_series: np.ndarray = None,
    ) -> None:
        '''
        Prepare this plot for subsequent display and/or saving.

        Specifically (in order):

        . If the current simulation configuration requests that each plotted
          cell by labelled with that cell's unique 0-based index, do so. To
          ensure that these labels are plotted over rather than under the
          contents of their corresponding cells, we do so only _after_ all
          subclass plotting has been performed by delaying such labelling to
          this method rather than the above `__init__()` method.
        . If the optional `axes_title` parameter was passed to `__init__()`:
          * Add the current `_figure_title` to this figure as a "super title."
          * Add the passed `axes_title` to this figure's axes as a "subtitle."
        . Else, add the current `_figure_title` to this figure's axes.
        . If the optional `colorbar_values` parameter is passed _and_ the
          current `clrAutoscale` boolean is `True`, clip the colorbar to the
          minimum and maximum values in the `colorbar_values` array.
        . Else, clip the colorbar to the current `clrMin` and `clrMax` values.
        . Add a colorbar whose:
          * Title is the current `_colorbar_title` string.
          * Mapping is the passed `_colorbar_mapping` parameter.

        Parameters
        ----------
        color_mapping : object
            Mandatory Matplotlib mapping (e.g., `Image`, `ContourSet`) to which
            this colorbar applies.
        color_series : np.ndarray
            Optional multi-dimensional Numpy array containing all data values
            to be animated _or_ `None` if calculating this data during the
            animation initialization is infeasible or impractical (e.g., due to
            space and time constraints). If non-`None` _and_ colorbar
            autoscaling is requested (i.e., the initialization-time
            `clrAutoscale` parameter was `True`), the colorbar will be clipped
            to the maximum and minimum value in this matrix; else, the subclass
            is responsible for colorbar autoscaling. Defaults to `None`.
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
        if color_series is not None and self._is_color_autoscaled is True:
            assert types.is_sequence_nonstr(color_series), (
                types.assert_not_sequence_nonstr(color_series))

            # Flatten this two-dimensional matrix to a one-dimensional array,
            # providing efficient retrieval of minimum and maximum values.
            time_series_flat = np.ravel(color_series)

            # Minimum and maximum values.
            self._color_min = np.min(time_series_flat)
            self._color_max = np.max(time_series_flat)

            # If these values are identical, coerce them to differ. This ensures
            # that the colorbar will be displayed in a sane manner.
            if self._color_min == self._color_max:
                self._color_min = self._color_min - 1
                self._color_max = self._color_max + 1

        # Set the colorbar range.
        color_mapping.set_clim(self._color_min, self._color_max)

        # Display the colorbar.
        colorbar = self._figure.colorbar(color_mapping)
        colorbar.set_label(self._colorbar_title)

    # ..................{ PRIVATE ~ plotter                  }..................
    def _plot_image(
        self, pixel_data: np.ndarray) -> 'matplotlib.image.AxesImage':
        '''
        Plot the passed pixel data onto the current frame's figure axes and
        return the resulting image.

        Parameters
        ----------
        pixel_data: np.ndarray
            Array of pixel data defining the image to be plotted, whose first
            two dimensions index the X and Y components of this axes. If this
            array is:
            * Two-dimensional, a greyscale image mapped onto the current
              colormap will be plotted.
            * Three-dimensional, an RGB image _not_ mapped onto the current
              colormap will be plotted.
            * Four-dimensional, an RGBa image _not_ mapped onto the current
              colormap will be plotted.

        Returns
        ----------
        matplotlib.image.AxesImage
            Image produced by plotting the passed pixel data.

        See Also
        ----------
        http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.imshow
            Further details on Matplotlib-based axes image plotting.
        '''
        return self._axes.imshow(
            pixel_data,
            origin='lower',
            extent=self._axes_bounds,
            cmap=self._p.background_cm,
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
        Plot all streamlines of the passed vector flow data onto the current
        frame's figure axes _and_ return the resulting streamplot.

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
            color='k',
            cmap=self._colormap,
            arrowsize=1.5,
        )
