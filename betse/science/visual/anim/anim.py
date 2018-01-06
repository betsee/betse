#!/usr/bin/env python3
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Matplotlib-based animation classes.
'''

#FIXME: To avoid memory leaks, I'm fairly certain that everywhere we currently
#call the ".lines.remove()" method of a Matplotlib streamplot object, that we
#instead need to simply call remove(): e.g.,
#
#    # ...instead of this:
#    self._stream_plot.lines.remove()
#
#    # ...do this:
#    self._stream_plot.remove()
#
#The latter removes the stream plot itself from the parent figure and axes, in
#preparation for recreating the stream plot. The former merely removes the set
#of all stream lines from the stream plot. Failing to remove the stream plot
#itself could potentially cause old stream plots to be retained in memory. In
#any case, give self._stream_plot.remove() a go; if that works, we'll want to
#use that regardless of whether there are any memory leaks here or not.
#FIXME: Sadly, the self._stream_plot.remove() method has yet to be fully
#implemented (...at least, under older Matplotlib versions). Is there no better
#approach than replacing the entire streamplot each animation frame? We're
#fairly certain that we've seen Matplotlib example code that updates an
#existing streamplot each animation frame instead. Investigate the aged pandas!

# ....................{ IMPORTS                            }....................
import numpy as np
from betse.lib.numpy import nparray
from betse.science.export import expmath
from betse.science.visual.anim.animafter import (
    AnimCellsAfterSolving, AnimVelocity)
from betse.science.visual.plot.plotutil import cell_mosaic, cell_mesh
from betse.util.type.types import type_check, SequenceTypes
from matplotlib.collections import LineCollection, PolyCollection
from scipy import interpolate

# ....................{ CLASSES ~ after                    }....................
#FIXME: This class should probably no longer be used, now that the Gouraud
#shading performed by the "AnimCellsMembranesData" class has been optimized.
class AnimFlatCellsTimeSeries(AnimCellsAfterSolving):
    '''
    Animation of an arbitrary cell-centric time series (e.g., average cell
    voltage as a function of time).

    Attributes
    ----------
    _cell_plot : ???
        Artists signifying cell data for the prior or current frame.
    _cell_time_series : list
        Arbitrary cell data as a function of time to be underlayed.
    _gapjunc_plot : LineCollection
        Lines signifying gap junction state for the prior or current frame.
    _gapjunc_time_series : list
        Arbitrary gap junction data as a function of time to be overlayed.
    '''

    @type_check
    def __init__(
        self,
        time_series: SequenceTypes,
        *args, **kwargs
    ) -> None:
        '''
        Initialize this animation.

        Parameters
        ----------
        cell_time_series : SequenceTypes
            Arbitrary cell data as a function of time

        See the superclass `__init__()` method for all remaining parameters.
        '''

        # Initialize the superclass.
        super().__init__(*args, time_step_count=len(time_series), **kwargs)

        # Classify all remaining parameters.
        self._cell_time_series = time_series

        # Cell data series for the first frame.
        data_set = self._cell_time_series[0]

        # Add a collection of cell polygons with animated voltage data.
        if self._phase.p.showCells is True:
            self._cell_plot, self._axes = cell_mosaic(
                data_set, self._axes, self._phase.cells, self._phase.p, self._colormap)
        else:
            self._cell_plot, self._axes = cell_mesh(
                data_set, self._axes, self._phase.cells, self._phase.p, self._colormap)

        # Display and/or save this animation.
        self._animate(
            color_mappables=self._cell_plot,
            color_data=self._cell_time_series,
        )


    @type_check
    def _plot_frame_figure(self):

        # Cell data series for this frame.
        zv = self._cell_time_series[self._time_step]
        if self._phase.p.showCells is True:
            zz_grid = zv
        else:
            zz_grid = zv
            # zz_grid = np.zeros(len(self._phase.cells.voronoi_centres))
            # zz_grid[self._phase.cells.cell_to_grid] = zv

        # Update the cell plot for this frame.
        self._cell_plot.set_array(zz_grid)


class AnimEnvTimeSeries(AnimCellsAfterSolving):
    '''
    Animation of an arbitrary cell-agnostic time series (e.g., environmental
    voltage as a function of time), plotted over the cell cluster.

    Attributes
    ----------
    _time_series : list
        Arbitrary environmental data as a function of time to be plotted.
    _mesh_plot : matplotlib.image.AxesImage
        Meshplot of the current or prior frame's environmental data.
    '''

    @type_check
    def __init__(
        self,
        time_series: SequenceTypes,
        *args, **kwargs
    ) -> None:
        '''
        Initialize this animation.

        Parameters
        ----------
        time_series : SequenceTypes
            Arbitrary environmental data as a function of time to be plotted.

        See the superclass `__init__()` method for all remaining parameters.
        '''

        # Initialize the superclass.
        super().__init__(*args, **kwargs)

        # Classify parameters required by the _plot_frame_figure() method.
        self._time_series = time_series

        # Environmental data meshplot for the first frame.
        self._mesh_plot = self._plot_image(pixel_data=self._time_series[0])

        # Display and/or save this animation.
        self._animate(
            color_mappables=self._mesh_plot,
            color_data=self._time_series,
        )


    def _plot_frame_figure(self) -> None:

        # Environmental data meshplot for this frame.
        self._mesh_plot.set_data(self._time_series[self._time_step])


#FIXME: Disable display of the cell-centric time series, leaving only the
#gap junction-centric time series displayed. Alternately, might it be possible
#to filter the cell-centric time series with some sort of alpha-based fading
#effect -- reducing the prominance of that series but not entirely eliminating
#that series? Contemplate us up, anyways.

class AnimGapJuncTimeSeries(AnimCellsAfterSolving):
    '''
    Animation of an arbitrary gap junction-centric time series (e.g., the gap
    junction open state as a function of time) overlayed an arbitrary cell-
    centric time series (e.g., cell voltage as a function of time) on the cell
    cluster.

    Attributes
    ----------
    _cell_plot : ???
        Artists signifying cell data for the prior or current frame.
    _cell_time_series : list
        Arbitrary cell data as a function of time to be underlayed.
    _gapjunc_plot : LineCollection
        Lines signifying gap junction state for the prior or current frame.
    _gapjunc_time_series : list
        Arbitrary gap junction data as a function of time to be overlayed.
    '''

    @type_check
    def __init__(
        self,
        time_series: SequenceTypes,
        *args, **kwargs
    ) -> None:
        '''
        Initialize this animation.

        Parameters
        ----------
        cell_time_series : SequenceTypes
            Arbitrary cell data as a function of time to be underlayed.
        gapjunc_time_series : SequenceTypes
            Arbitrary gap junction data as a function of time to be overlayed.

        See the superclass `__init__()` method for all remaining parameters.
        '''

        # Initialize the superclass.
        super().__init__(*args, time_step_count=len(time_series), **kwargs)

        # Classify all remaining parameters.
        # self._cell_time_series = cell_time_series
        self._time_series = time_series

        # Gap junction data series for the first frame plotted as lines.
        self._gapjunc_plot = LineCollection(
            expmath.upscale_coordinates(self._phase.cells.nn_edges),
            array=self._time_series[0],
            cmap=self._phase.p.gj_cm,
            linewidths=2.0,
            zorder=10,
        )
        # self._gapjunc_plot.set_clim(0.0, 1.0)
        self._axes.add_collection(self._gapjunc_plot)

        # Cell data series for the first frame.
        # data_set = self._cell_time_series[0]

        # # Add a collection of cell polygons with animated voltage data.
        # if self.p.showCells is True:
        #     self._cell_plot, self._axes = cell_mosaic(
        #         data_set, self._axes, self.cells, self.p, self._colormap)
        # else:
        #     self._cell_plot, self._axes = cell_mesh(
        #         data_set, self._axes, self.cells, self.p, self._colormap)

        # Display and/or save this animation.
        self._animate(
            color_mappables=self._gapjunc_plot,
            color_data=self._time_series,
        )

    #
    # @type_check
    def _plot_frame_figure(self):

        # Update the gap junction plot for this frame.
        self._gapjunc_plot.set_array(
            self._time_series[self._time_step])

        # # Cell data series for this frame.
        # zv = self._cell_time_series[self._time_step]
        # if self.p.showCells is True:
        #     zz_grid = zv
        # else:
        #     zz_grid = np.zeros(len(self.cells.voronoi_centres))
        #     zz_grid[self.cells.cell_to_grid] = zv

        # Update the cell plot for this frame.
        # self._cell_plot.set_array(zz_grid)


class AnimMembraneTimeSeries(AnimCellsAfterSolving):
    '''
    Animation of an arbitrary cell membrane-specific time series (e.g.,
    membrane channel or pump density factor as a function of time), plotted
    over the cell cluster.

    This factor changes in response to changes in electroosmotic and
    electrophoretic movements, produced by self-generated fields and flows in
    the cell cluster.

    Attributes
    ----------
    _mem_edges : LineCollection
        Membrane edges coloured for the current or prior frame.
    _time_series : SequenceTypes
        Arbitrary cell membrane data as a function of time to be plotted.
    '''


    @type_check
    def __init__(
        self,
        time_series: SequenceTypes,
        *args, **kwargs
    ) -> None:
        '''
        Initialize this animation.

        Parameters
        ----------
        time_series : Sequence
            Arbitrary cell membrane data as a function of time to be plotted.

        See the superclass `__init__()` method for all remaining parameters.
        '''

        # Initialize the superclass.
        super().__init__(*args, time_step_count=len(time_series), **kwargs)

        # Classify parameters required by the _plot_frame_figure() method.
        self._time_series = time_series

        # Membrane edges coloured for the first frame.
        self._mem_edges = LineCollection(
            self._phase.cells.mem_edges_flat * self._phase.p.um,
            array=self._time_series[0],
            cmap=self._colormap,
            linewidths=4.0,
        )
        self._axes.add_collection(self._mem_edges)

        # Display and/or save this animation.
        self._animate(
            color_mappables=self._mem_edges,
            color_data=self._time_series,
        )


    def _plot_frame_figure(self) -> None:

        # Update membrane edges colours for this frame.
        self._mem_edges.set_array(
            self._time_series[self._time_step])


#FIXME: This animation class no longer appears to be used. Consider excising.
#FIXME: Before we eliminate this class entirely, it would be quite useful to
#ensure that the "PolyCollection" logic performed below exists as a
#"LayerCellsABC" subclass.
class AnimMorphogenTimeSeries(AnimCellsAfterSolving):
    '''
    Animation of the concentration of an arbitrary morphogen in both cells and
    the environment as a function of time, plotted over the cell cluster.

    Parameters
    ----------
    _cell_time_series : Sequence
        Morphogen concentration in cells as a function of time.
    _env_time_series : Sequence
        Morphogen concentration in the environment as a function of time.
    '''


    @type_check
    def __init__(
        self,
        cell_time_series: SequenceTypes,
        env_time_series: SequenceTypes,
        *args, **kwargs
    ) -> None:
        '''
        Initialize this animation.

        Parameters
        ----------
        cell_time_series : Sequence
            Morphogen concentration in cells as a function of time.
        env_time_series : Sequence
            Morphogen concentration in the environment as a function of time.

        See the superclass `__init__()` method for all remaining parameters.
        '''

        # Initialize the superclass.
        super().__init__(*args, time_step_count=len(cell_time_series), **kwargs)

        # Classify the passed parameters.
        self._cell_time_series = cell_time_series
        self._env_time_series = env_time_series

        #FIXME: Rename:
        #
        #* "bkgPlot" to "_"... we have no idea. Animate first. Decide later.
        #* "collection" to "_mesh_plot".

        self.bkgPlot = self._plot_image(
            pixel_data=self._env_time_series[0].reshape(self._phase.cells.X.shape))

        #FIXME: Try reducing to: self.cells.cell_verts * self.p.um

        # Polygon collection based on individual cell polygons.
        points = np.multiply(self._phase.cells.cell_verts, self._phase.p.um)
        self.collection = PolyCollection(
            points, cmap=self._colormap, edgecolors='none')
        self.collection.set_array(self._cell_time_series[0])
        self._axes.add_collection(self.collection)

        # Display and/or save this animation.
        self._animate(
            color_mappables=(self.collection, self.bkgPlot),

            # If colorbar autoscaling is requested, clip the colorbar to the
            # minimum and maximum morphogen concentrations -- regardless of
            # whether that morphogen resides in cells or the environment.
            color_data=(
                self._cell_time_series, self._env_time_series),
        )


    def _plot_frame_figure(self) -> None:

        self.collection.set_array(
            self._cell_time_series[self._time_step])
        self.bkgPlot.set_data(
            self._env_time_series[self._time_step].reshape(
                self._phase.cells.X.shape))

# ....................{ SUBCLASSES ~ velocity              }....................
class AnimVelocityIntracellular(AnimVelocity):
    '''
    Animation of fluid velocity over all intracellular gap junctions plotted on
    the cell cluster.

    Attributes
    -----------
    _mesh_plot : matplotlib.image.AxesImage
        Meshplot of the current or prior frame's velocity field magnitude.
    _stream_plot : matplotlib.streamplot.StreamplotSet
        Streamplot of the current or prior frame's velocity field _or_ `None`
        if such field has yet to be streamplotted.
    '''


    def __init__(self, *args, **kwargs) -> None:

        # Initialize the superclass.
        super().__init__(*args, **kwargs)

        # Define this attribute *BEFORE* streamplotting, which assumes this
        # attribute to exist.
        self._stream_plot = None

        #FIXME: Inefficient. This streamplot will be recreated for the first
        #time step in the exact same manner; so, it's unclear that we need to
        #do so here.

        # Streamplot the first frame's velocity field.
        vfield, vnorm = self._plot_stream_velocity_field(time_step=0)

        # Meshplot the first frame's velocity field magnitude.
        self._mesh_plot = self._plot_image(
            pixel_data=vfield,
            colormap=self._phase.p.background_cm,
        )

        #FIXME: How expensive would caching these calculations be? Oh, just do
        #it already. We have to calculate these values *ANYWAY*, so there's no
        #incentive at all in delaying the matter.

        # Display and/or save this animation. Since recalculating "vfield" for
        # each animation frame is non-trivial, this call avoids passing the
        # "time_series" parameter. Instead, the _get_velocity_field() method
        # manually rescales the colorbar on each frame according to the minimum
        # and maximum velocity field magnitude. While non-ideal, every
        # alternative is currently worse.
        self._animate(color_mappables=self._mesh_plot)


    def _plot_frame_figure(self) -> None:

        # Streamplot this frame's velocity field.
        vfield, vnorm = self._plot_stream_velocity_field(self._time_step)

        # Meshplot this frame's velocity field.
        self._mesh_plot.set_data(vfield)

        # Rescale the colorbar if needed.
        self._mesh_plot.set_clim(self._color_min, self._color_max)


    @type_check
    def _plot_stream_velocity_field(self, time_step: int) -> tuple:
        '''
        Streamplot the current velocity field for the passed frame and return a
        2-tuple describing this field.

        Returns
        ----------
        (Sequence, float)
            2-element tuple `(velocity_field, velocity_field_magnitude_max)`
            whose:
            * First element is the list of all velocity field magnitudes for
              this frame.
            * Second element is the maximum such magnitude.
        '''

        cell_centres = (
            self._phase.cells.cell_centres[:, 0],
            self._phase.cells.cell_centres[:, 1])
        cell_grid = (self._phase.cells.X, self._phase.cells.Y)

        #FIXME: Ugh. Duplicate code already performed by the superclass
        #AnimCellsABC._init_current_density() method. We clearly need a
        #general-purpose interpolation utility method. Hawkish doves in a cove!
        u_gj_x = self._phase.cells.maskECM * interpolate.griddata(
            cell_centres,
            self._phase.sim.u_cells_x_time[time_step],
            cell_grid,
            fill_value=0,
            method=self._phase.p.interp_type,
        )
        u_gj_y = self._phase.cells.maskECM * interpolate.griddata(
            cell_centres,
            self._phase.sim.u_cells_y_time[time_step],
            cell_grid,
            fill_value=0,
            method=self._phase.p.interp_type,
        )

        # Current velocity field magnitudes and the maximum such magnitude.
        vfield = np.sqrt(u_gj_x**2 + u_gj_y**2) * 1e9
        vnorm = np.max(vfield)

        # Streamplot the current velocity field for this frame.
        self._stream_plot = self._plot_stream(
            old_stream_plot=self._stream_plot,
            x=u_gj_x / vnorm,
            y=u_gj_y / vnorm,
            magnitude=vfield,
            magnitude_max=vnorm,
        )
        # self.streamV.set_UVC(u_gj_x/vnorm,u_gj_y/vnorm)

        # Rescale the colorbar range if desired.
        if self._conf.is_color_autoscaled:
            self._color_min = np.min(vfield)
            self._color_max = vnorm

        return (vfield, vnorm)


class AnimVelocityExtracellular(AnimVelocity):
    '''
    Animation of fluid velocity over all extracellular spaces plotted on the
    cell cluster.

    Attributes
    ----------
    _mesh_plot : matplotlib.image.AxesImage
        Meshplot of the current or prior frame's velocity field magnitude.
    _stream_plot : matplotlib.streamplot.StreamplotSet
        Streamplot of the current or prior frame's velocity field.
    _magnitude_time_series : np.ndarray
        Time series of all fluid velocity magnitudes.
    '''

    def __init__(self, *args, **kwargs) -> None:

        # Initialize the superclass.
        super().__init__(*args, is_ecm_required=True, **kwargs)

        # Time series of all velocity magnitudes.
        self._magnitude_time_series = np.sqrt(
            nparray.from_iterable(self._phase.sim.u_env_x_time) ** 2 +
            nparray.from_iterable(self._phase.sim.u_env_y_time) ** 2) * 1e6

        # Velocity field and maximum velocity field value for the first frame.
        vfield = self._magnitude_time_series[0]
        vnorm = np.max(vfield)

        # Velocity field meshplot for the first frame.
        self._mesh_plot = self._plot_image(
            pixel_data=vfield,
            colormap=self._phase.p.background_cm,
        )

        #FIXME: Doesn't this streamplot the last frame instead?

        # Velocity field streamplot for the first frame.
        self._stream_plot = self._axes.quiver(
            self._phase.cells.xypts[:, 0] * self._phase.p.um,
            self._phase.cells.xypts[:, 1] * self._phase.p.um,
            self._phase.sim.u_env_x_time[-1].ravel() / vnorm,
            self._phase.sim.u_env_y_time[-1].ravel() / vnorm,
        )

        # Display and/or save this animation.
        self._animate(
            color_mappables=self._mesh_plot,
            color_data=self._magnitude_time_series,
        )


    def _plot_frame_figure(self) -> None:

        # Velocity field and maximum velocity field value for this frame.
        vfield = self._magnitude_time_series[self._time_step]
        vnorm = np.max(vfield)

        # Update the current velocity meshplot.
        self._mesh_plot.set_data(vfield)

        # Update the current velocity streamplot.
        self._stream_plot.set_UVC(
            self._phase.sim.u_env_x_time[self._time_step] / vnorm,
            self._phase.sim.u_env_y_time[self._time_step] / vnorm,
        )
