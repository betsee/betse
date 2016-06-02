#!/usr/bin/env python3
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
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
from enum import Enum

import numpy as np
from matplotlib import animation
from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection, PolyCollection
from numpy import ma as ma
from scipy import interpolate

from betse.exceptions import BetseExceptionParameters
from betse.lib.matplotlib.anim import FileFrameWriter
from betse.science.plot import plot
from betse.science.plot.anim.abc import (
    AnimCells, AnimField, AnimVelocity)
from betse.util.io.log import logs
from betse.util.path import dirs, paths
from betse.util.type import types

#FIXME: Shift functions called only by this module either to a new
#"betse.science.plot.animation.helper" module or possibly as private
#methods of the "Anim" superclass. Right. After investigation, absolutely the
#latter approach. This should permit us to avoid passing *ANY* parameters to
#these methods, which is rather nice.
from betse.science.plot.plot import (
    _setup_file_saving, env_mesh, cell_mosaic, cell_mesh,
    env_quiver, cell_quiver, cell_stream, pretty_patch_plot,
)

# ....................{ CLASSES ~ while                    }....................
#FIXME: Shift this subclass into a new "while.py" submodule of this package.
#FIXME: There appears to be a largely ignorable aesthetic issue with
#deformations. The first time step centers the cell cluster; all subsequent
#time steps shift the cell cluster towards the top of the figure window. While
#nothing clips across edges, the movement is noticeably... awkward. Since this
#occurs in the older "PlotWhileSolving" implementation as well, we ignore this
#for the moment.

#FIXME: Rename "_cell_data_plot" to "_cell_body_plot".
#FIXME: Rename "_cell_edges_plot" to "_cell_edge_plot".

class AnimCellsWhileSolving(AnimCells):
    '''
    In-place animation of an arbitrary membrane-centric time series (e.g., cell
    Vmem as a function of time), plotted over the cell cluster during rather
    than after simulation modelling.

    Attributes
    ----------
    _cell_edges_plot : LineCollection
        Plot of the current or prior frame's cell edges.
    _cell_data_plot : Collection
        Plot of the current or prior frame's cell contents.
    _cell_verts_id : int
        Unique identifier for the array of cell vertices (i.e.,
        `_cells.cell_verts`) when plotting the current or prior frame.
        Retaining this identifier permits the `_plot_frame_figure()` method to
        efficiently detect and respond to physical changes (e.g., deformation
        forces, cutting events) in the fundamental structure of the previously
        plotted cell cluster.
    _is_colorbar_autoscaling_telescoped : bool
        `True` if colorbar autoscaling is permitted to increase but _not_
        decrease the colorbar range _or_ `False` otherwise (i.e., if
        colorbar autoscaling is permitted to both increase and decrease the
        colorbar range). Such telescoping assists in emphasizing stable
        long-term patterns in cell data at a cost of deemphasizing unstable
        short-term patterns. If colorbar autoscaling is disabled (i.e.,
        `is_color_autoscaled` is `False`), this will be ignored.
    _is_time_step_first : bool
        `True` only if the current frame being animated is the first.
    '''

    def __init__(
        self,

        #FIXME: Permit this option to be configured via a new configuration
        #option in the YAML file.
        is_colorbar_autoscaling_telescoped: bool = False,
        *args, **kwargs
    ) -> None:
        '''
        Initialize this animation.

        Parameters
        ----------
        is_colorbar_autoscaling_telescoped : bool
            `True` if colorbar autoscaling is permitted to increase but _not_
            decrease the colorbar range _or_ `False` otherwise (i.e., if
            colorbar autoscaling is permitted to both increase and decrease the
            colorbar range). Such telescoping assists in emphasizing stable
            long-term patterns in cell data at a cost of deemphasizing unstable
            short-term patterns. If colorbar autoscaling is disabled (i.e.,
            `is_color_autoscaled` is `False`), this will be ignored. Defaults
            to `False`.

        See the superclass `__init__()` method for all remaining parameters.
        '''
        # assert types.is_sequence_nonstr(cell_time_series), (
        #     types.assert_not_sequence_nonstr(cell_time_series))

        # Pass all parameters *NOT* listed above to our superclass.
        super().__init__(
            # Save in-place animation to a different parent directory than that
            # to which out-of-place animations are saved.
            save_dir_parent_basename='anim_while_solving',

            # Prevent the superclass from plotting electric current or
            # concentration flux. Although this class does *NOT* plot a
            # streamplot, the time series required for the current overlay is
            # unavailable until *AFTER* simulation modelling completes.
            is_current_overlayable=False,
            *args, **kwargs
        )

        # Classify all remaining parameters.
        self._is_colorbar_autoscaling_telescoped = (
            is_colorbar_autoscaling_telescoped)

        # Unique identifier for the array of cell vertices. (See docstring.)
        self._cell_verts_id = id(self._cells.cell_verts)

        # "True" only if the current frame being animated is the first.
        self._is_time_step_first = True

        # average the voltage to the cell centre
        # FIXME this is a temp change until we get this right
        vm_o = np.dot(self._cells.M_sum_mems,self._sim.vm)/self._cells.num_mems

        # self._cell_time_series = self._sim.vm_time
        self._cell_time_series = self._sim.vm_ave_time

        # cell_data_current = self._sim.vm
        cell_data_current = vm_o

        # Upscaled cell data for the first frame.
        cell_data = plot.upscale_data(cell_data_current)

        # Collection of cell polygons with animated voltage data.
        #
        # If *NOT* simulating extracellular spaces, only animate intracellular
        # spaces.
        self._cell_data_plot = self._plot_cells_sans_ecm(
            cell_data=cell_data)


        cell_data_vm = cell_data

        # Perform all superclass plotting preparation immediately *BEFORE*
        # plotting this animation's first frame.
        self._prep_figure(
            color_mapping=self._cell_data_plot,
            color_series=cell_data_vm,
        )

        # Plot this animation's first frame in a non-blocking manner.
        plt.show(block=False)

    # ..................{ PROPERTIES                         }..................
    # In-place animation is governed by different configuration file booleans
    # than are post-simulation animations.
    @property
    def _is_showing(self) -> bool:
        return self._p.plot_while_solving


    @property
    def _is_saving(self) -> bool:
        return self._p.save_solving_plot

    # ..................{ PLOTTERS                           }..................
    def _plot_frame_figure(self, frame_number: int) -> None:

        # Upscaled cell data for the current time step.
        cell_data = plot.upscale_data(self._cell_time_series[self._time_step])

        #FIXME: Duplicated from above. What we probably want to do is define a
        #new _get_cell_data() method returning this array in a centralized
        #manner callable both here and above. Rays of deluded beaming sunspray!


        # If the unique identifier for the array of cell vertices has *NOT*
        # changed, the cell cluster has *NOT* fundamentally changed and need
        # only be updated with this time step's cell data.
        if self._cell_verts_id == id(self._cells.cell_verts):
            # loggers.log_info(
            #     'Updating animation "{}" cell plots...'.format(self._type))
            self._update_cell_plots(cell_data=cell_data)
        # Else, the cell cluster has fundamentally changed (e.g., due to
        # physical deformations or cutting events) and must be recreated.
        else:
            # loggers.log_info(
            #     'Reviving animation "{}" cell plots...'.format(self._type))

            # Prevent subsequent calls to this method from erroneously
            # recreating the cell cluster again.
            self._cell_verts_id = id(self._cells.cell_verts)

            # Recreate the cell cluster.
            self._revive_cell_plots(cell_data=cell_data)

        # Update the color bar with the content of the cell body plot *AFTER*
        # possibly recreating this plot above.
        if self._is_color_autoscaled:

            cell_data_vm = cell_data

            # If autoscaling this colorbar in a telescoping manner and this is
            # *NOT* the first time step, do so.
            #
            # If this is the first time step, the previously calculated minimum
            # and maximum colors are garbage and thus *MUST* be ignored.
            if (self._is_colorbar_autoscaling_telescoped and
                not self._is_time_step_first):
                self._color_min = min(self._color_min, np.min(cell_data_vm))
                self._color_max = max(self._color_max, np.max(cell_data_vm))
            # Else, autoscale this colorbar in an unrestricted manner.
            else:
                self._color_min = np.min(cell_data_vm)
                self._color_max = np.max(cell_data_vm)

                # If autoscaling this colorbar in a telescoping manner, do so
                # for subsequent time steps.
                self._is_time_step_first = False

            # Autoscale the colorbar to these colors.
            self._rescale_colors()

        # If displaying this frame, do so.
        if self._is_showing:
            self._figure.canvas.draw()


    def _update_cell_plots(self, cell_data: np.ndarray) -> None:
        '''
        Update _without_ recreating all cell plots for this time step with the
        passed array of arbitrary cell data.

        This method is intended to be called _unless_ physical changes
        (e.g., deformation forces, cutting events) in the underlying structure
        of the cell cluster have occurred for this simulation time step.

        Parameters
        -----------
        cell_data : np.ndarray
            Arbitrary cell data defined on an environmental grid to be plotted.
        '''
        assert types.is_sequence_nonstr(cell_data), (
            types.assert_not_sequence_nonstr(cell_data))

        self._update_cell_plot_sans_ecm(
            cell_plot=self._cell_data_plot,
            cell_data=cell_data)


    def _revive_cell_plots(self, cell_data: np.ndarray) -> None:
        '''
        Recreate all cell plots for this time step with the passed array of
        arbitrary cell data.

        This method is intended to be called in response to physical changes
        (e.g., deformation forces, cutting events) in the underlying structure
        of the cell cluster for this simulation time step. This method is both
        inefficient and destructive, and should be called only when needed.

        Parameters
        -----------
        cell_data : np.ndarray
            Arbitrary cell data defined on an environmental grid to be plotted.
        '''
        assert types.is_sequence_nonstr(cell_data), (
            types.assert_not_sequence_nonstr(cell_data))

        self._cell_data_plot = self._revive_cell_plots_sans_ecm(
            cell_plot=self._cell_data_plot,
            cell_data=cell_data)

# ....................{ CLASSES ~ after                    }....................
class AnimCellsTimeSeries(AnimCells):
    '''
    Post-simulation animation of an arbitrary cell-centric time series (e.g.,
    cell voltage as a function of time), plotted over the cell cluster.

    Attributes
    ----------
    _is_ecm_ignored : bool
        `True` if ignoring extracellular spaces _or_ `False` otherwise.
    _time_series : numpy.ndarray
        Arbitrary cell data as a function of time to be plotted.
    '''

    def __init__(
        self,
        time_series: np.ndarray,
        is_ecm_ignored: bool = True,
        scaling_series: np.ndarray = None,
        *args, **kwargs
    ) -> None:
        '''
        Initialize this animation.

        Parameters
        ----------
        time_series : np.ndarray
            Arbitrary cell data as a function of time to be plotted.
        is_ecm_ignored : bool
            `True` if ignoring extracellular spaces _or_ `False` otherwise.
            Defaults to `True`.

        See the superclass `__init__()` method for all remaining parameters.
        '''
        assert types.is_sequence_nonstr(time_series), (
            types.assert_not_sequence_nonstr(time_series))
        assert types.is_bool(is_ecm_ignored), (
            types.assert_not_bool(is_ecm_ignored))

        # Pass all parameters *NOT* listed above to our superclass.
        super().__init__(
            # Since this class does *NOT* plot a streamplot, request that the
            # superclass do so for electric current or concentration flux.
            is_current_overlayable=True,
            *args, **kwargs
        )

        # Classify parameters required by the _plot_frame_figure() method.
        self._time_series = time_series
        self._is_ecm_ignored = is_ecm_ignored

        # Cell data for the first frame.
        data_points = self._time_series[0]
        data_verts = np.dot(self._cells.matrixMap2Verts, data_points)

        #FIXME: Rename "self.collection" to something more descriptive.

        if self._p.showCells is True:
            self.collection, self._axes = pretty_patch_plot(
                data_verts, self._axes, self._cells, self._p, self._colormap,
                cmin=self._color_min, cmax=self._color_max)

        # Else, plot a smooth continuum approximating the cell cluster.
        else:
            self.collection, self._axes = cell_mesh(
                data_points, self._axes, self._cells, self._p, self._colormap)


        # Display and/or save this animation.
        self._animate(
            frame_count=len(self._sim.time),
            color_mapping=self.collection,
            color_series=scaling_series if scaling_series else self._time_series,
        )

    def _plot_frame_figure(self, frame_number: int):
        assert types.is_int(frame_number), types.assert_not_int(frame_number)

        zz = self._time_series[frame_number]
        data_verts = np.dot(self._cells.matrixMap2Verts, zz)

        if self._p.showCells is True:
            self._axes.cla()
            self._axes.axis('equal')
            self._axes.axis(self._axes_bounds)
            self._axes.set_xlabel('Spatial distance [um]')
            self._axes.set_ylabel('Spatial distance [um]')

            # self.collection.set_array(zz)
            self.collection, self._axes = pretty_patch_plot(
                data_verts, self._axes, self._cells, self._p, self._colormap,
                cmin=self._color_min, cmax=self._color_max)

        else:
            zz_grid = np.zeros(len(self._cells.voronoi_centres))
            zz_grid[self._cells.cell_to_grid] = zz
            self.collection.set_array(zz_grid)

class AnimEnvTimeSeries(AnimCells):
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

    def __init__(
        self,
        time_series: np.ndarray,
        *args, **kwargs
    ) -> None:
        '''
        Initialize this animation.

        Parameters
        ----------
        time_series : np.ndarray
            Arbitrary environmental data as a function of time to be plotted.

        See the superclass `__init__()` method for all remaining parameters.
        '''
        assert types.is_sequence_nonstr(time_series), (
            types.assert_not_sequence_nonstr(time_series))

        # Pass all parameters *NOT* listed above to our superclass.
        super().__init__(
            # Since this class does *NOT* plot a streamplot, request that the
            # superclass do so for electric current or concentration flux.
            is_current_overlayable=True,
            *args, **kwargs
        )

        # Classify parameters required by the _plot_frame_figure() method.
        self._time_series = time_series

        # Environmental data meshplot for the first frame.
        self._mesh_plot = self._plot_image(pixel_data=self._time_series[0])

        # Display and/or save this animation.
        self._animate(
            frame_count=len(self._sim.time),
            color_mapping=self._mesh_plot,
            color_series=self._time_series,
        )


    def _plot_frame_figure(self, frame_number: int):
        assert types.is_int(frame_number), types.assert_not_int(frame_number)

        # Environmental data meshplot for this frame.
        self._mesh_plot.set_data(self._time_series[frame_number])

class AnimGapJuncTimeSeries(AnimCells):
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

    def __init__(
        self,
        cell_time_series: list,
        gapjunc_time_series: list,
        *args, **kwargs
    ) -> None:
        '''
        Initialize this animation.

        Parameters
        ----------
        cell_time_series : list
            Arbitrary cell data as a function of time to be underlayed.
        gapjunc_time_series : list
            Arbitrary gap junction data as a function of time to be overlayed.

        See the superclass `__init__()` method for all remaining parameters.
        '''
        assert types.is_sequence_nonstr(cell_time_series), (
            types.assert_not_sequence_nonstr(cell_time_series))
        assert types.is_sequence_nonstr(gapjunc_time_series), (
            types.assert_not_sequence_nonstr(gapjunc_time_series))

        # Pass all parameters *NOT* listed above to our superclass.
        super().__init__(*args, **kwargs)

        # Classify all remaining parameters.
        self._cell_time_series = cell_time_series
        self._gapjunc_time_series = gapjunc_time_series

        # Gap junction data series for the first frame plotted as lines.
        self._gapjunc_plot = LineCollection(
            np.asarray(self._cells.nn_edges) * self._p.um,
            array=self._gapjunc_time_series[0],
            cmap=self._p.gj_cm,
            linewidths=2.0,

            #FIXME: This isn't the best. Text should still be displayed
            #above this. Ideally, this should be something like
            #"(ZORDER_LINE + ZORDER_STREAM) / 2" instead.
            zorder=10,
        )
        self._gapjunc_plot.set_clim(0.0, 0.5)
        self._axes.add_collection(self._gapjunc_plot)

        # Cell data series for the first frame.
        data_set = self._cell_time_series[0]

        # Add a collection of cell polygons with animated voltage data.
        if self._p.showCells is True:
            self._cell_plot, self._axes = cell_mosaic(
                data_set, self._axes, self._cells, self._p, self._colormap)
        else:
            self._cell_plot, self._axes = cell_mesh(
                data_set, self._axes, self._cells, self._p, self._colormap)

        # Display and/or save this animation.
        self._animate(
            frame_count=len(self._gapjunc_time_series),
            color_mapping=self._cell_plot,
            color_series=self._cell_time_series,
        )


    def _plot_frame_figure(self, frame_number: int):
        assert types.is_int(frame_number), types.assert_not_int(frame_number)

        # Update the gap junction plot for this frame.
        self._gapjunc_plot.set_array(self._gapjunc_time_series[frame_number])

        # Cell data series for this frame.
        zv = self._cell_time_series[frame_number]
        if self._p.showCells is True:
            zz_grid = zv
        else:
            zz_grid = np.zeros(len(self._cells.voronoi_centres))
            zz_grid[self._cells.cell_to_grid] = zv

        # Update the cell plot for this frame.
        self._cell_plot.set_array(zz_grid)

class AnimMembraneTimeSeries(AnimCells):
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
    _time_series : list
        Arbitrary cell membrane data as a function of time to be plotted.
    '''

    def __init__(
        self,
        time_series: np.ndarray,
        *args, **kwargs
    ) -> None:
        '''
        Initialize this animation.

        Parameters
        ----------
        time_series : list
            Arbitrary cell membrane data as a function of time to be plotted.

        See the superclass `__init__()` method for all remaining parameters.
        '''
        assert types.is_sequence_nonstr(time_series), (
            types.assert_not_sequence_nonstr(time_series))

        # Pass all parameters *NOT* listed above to our superclass.
        super().__init__(
            # Since this class does *NOT* plot a streamplot, request that the
            # superclass do so for electric current or concentration flux.
            is_current_overlayable=True,
            *args, **kwargs
        )

        # Classify parameters required by the _plot_frame_figure() method.
        self._time_series = time_series

        # Membrane edges coloured for the first frame.
        self._mem_edges = LineCollection(
            self._cells.mem_edges_flat * self._p.um,
            array=self._time_series[0],
            cmap=self._colormap,
            linewidths=4.0,
        )
        self._axes.add_collection(self._mem_edges)

        # Display and/or save this animation.
        self._animate(
            frame_count=len(self._time_series),
            color_mapping=self._mem_edges,
            color_series=self._time_series,
        )


    def _plot_frame_figure(self, frame_number: int):
        assert types.is_int(frame_number), types.assert_not_int(frame_number)

        # Update membrane edges colours for this frame.
        self._mem_edges.set_array(self._time_series[frame_number])

class AnimMorphogenTimeSeries(AnimCells):
    '''
    Animation of the concentration of an arbitrary morphogen in both cells and
    the environment as a function of time, plotted over the cell cluster.

    Parameters
    ----------
    _cell_time_series : np.ndarray
        Morphogen concentration in cells as a function of time.
    _env_time_series : np.ndarray
        Morphogen concentration in the environment as a function of time.
    '''

    def __init__(
        self,
        cell_time_series: np.ndarray,
        env_time_series: np.ndarray,
        *args, **kwargs
    ) -> None:
        '''
        Initialize this animation.

        Parameters
        ----------
        cell_time_series : np.ndarray
            Morphogen concentration in cells as a function of time.
        env_time_series : np.ndarray
            Morphogen concentration in the environment as a function of time.

        See the superclass `__init__()` method for all remaining parameters.
        '''
        assert types.is_sequence_nonstr(cell_time_series), (
            types.assert_not_sequence_nonstr(cell_time_series))
        assert types.is_sequence_nonstr(env_time_series), (
            types.assert_not_sequence_nonstr(env_time_series))

        # Pass all parameters *NOT* listed above to our superclass.
        super().__init__(
            # Since this subclass plots no streamplot, request that the
            # superclass do so.
            is_current_overlayable=True,
            *args, **kwargs
        )

        # Classify parameters required by the _plot_frame_figure() method.
        self._cell_time_series = env_time_series
        self._env_time_series = env_time_series

        #FIXME: Rename:
        #
        #* "bkgPlot" to "_"... we have no idea. Animate first. Decide later.
        #* "collection" to "_mesh_plot".

        self.bkgPlot = self._plot_image(
            pixel_data=self._env_time_series[0].reshape(self._cells.X.shape))

        #FIXME: Try reducing to: self._cells.cell_verts * self._p.um
        # Polygon collection based on individual cell polygons.
        points = np.multiply(self._cells.cell_verts, self._p.um)
        self.collection = PolyCollection(
            points, cmap=self._colormap, edgecolors='none')
        self.collection.set_array(self._cell_time_series[0])
        self._axes.add_collection(self.collection)

        # Display and/or save this animation.
        self._animate(
            frame_count=len(self._cell_time_series),
            color_mapping=(self.collection, self.bkgPlot),

            # If colorbar autoscaling is requested, clip the colorbar to the
            # minimum and maximum morphogen concentrations -- regardless of
            # whether that morphogen resides in cells or the environment.
            color_series=(
                self._cell_time_series, self._env_time_series),
        )


    def _plot_frame_figure(self, frame_number: int):
        assert types.is_int(frame_number), types.assert_not_int(frame_number)

        self.collection.set_array(self._cell_time_series[frame_number])
        self.bkgPlot.set_data(
            self._env_time_series[frame_number].reshape(self._cells.X.shape))

# ....................{ SUBCLASSES ~ field                 }....................
class AnimFieldIntracellular(AnimField):
    '''
    Animation of the electric field over all intracellular gap junctions
    plotted on the cell cluster.

    Attributes
    ----------
    _mesh_plot : matplotlib.image.AxesImage
        Meshplot of the current or prior frame's electric field magnitude.
    _stream_plot : matplotlib.streamplot.StreamplotSet
        Streamplot of the current or prior frame's electric field.
    '''

    def __init__(self, *args, **kwargs) -> None:

        # Pass all parameters *NOT* listed above to our superclass.
        super().__init__(*args, **kwargs)

        # Number of frames to be animated.
        frame_count = len(self._sim.time)

        # Define electric field arrays for all animation frames.
        for frame_number in range(0, frame_count):
            # Electric field X and Y unit components for this frame.
            field_x = self._x_time_series[frame_number]
            field_y = self._y_time_series[frame_number]

            #FIXME: What's this about then? Buttercups and bitter nightingales!
            if len(field_x) != len(self._cells.cell_i):
                field_x = (
                    np.dot(self._cells.M_sum_mems, field_x) /
                    self._cells.num_mems)
                field_y = (
                    np.dot(self._cells.M_sum_mems, field_y) /
                    self._cells.num_mems)

            # Electric field magnitudes for this frame.
            field_magnitude = np.sqrt(field_x**2 + field_y**2)

            # If all such magnitudes are non-zero and hence safely divisible,
            # reduce all electric field X and Y components to unit vectors. To
            # avoid modifying the original arrays, use the "/" rather than
            # "/=" operator. The former produces a new array, whereas the
            # latter modifies the existing array in-place.
            if field_magnitude.all() != 0.0:
                field_x = field_x / field_magnitude
                field_y = field_y / field_magnitude

            # Add all such quantities to the corresponding lists.
            self._magnitude_time_series.append(field_magnitude)
            self._unit_x_time_series.append(field_x)
            self._unit_y_time_series.append(field_y)

        #FIXME: Why the last rather than first time step for the first frame?
        #(This seems increasingly wrong the more I bleakly stare at it.)

        # Electric field streamplot for the first frame.
        self._stream_plot, self._axes = cell_quiver(
            self._unit_x_time_series[-1], self._unit_y_time_series[-1],
            self._axes, self._cells, self._p)

        # Electric field magnitude meshplot for the first frame.
        self._mesh_plot, self._axes = cell_mesh(
            self._magnitude_time_series[-1],
            self._axes, self._cells, self._p, self._colormap)

        # Display and/or save this animation.
        self._animate(
            frame_count=frame_count,
            color_mapping=self._mesh_plot,
            color_series=self._magnitude_time_series,
        )


    def _plot_frame_figure(self, frame_number: int):
        assert types.is_int(frame_number), types.assert_not_int(frame_number)

        #FIXME: This is probably code copied from the helpers called above
        #(e.g., cell_mesh()). It'd be great to centralize this code somewhere.
        #Skinny trees fronded with blue ribbons!
        emag_grid = np.zeros(len(self._cells.voronoi_centres))
        emag_grid[self._cells.cell_to_grid] = (
            self._magnitude_time_series[frame_number])

        # Electric field streamplot for this frame.
        self._stream_plot.set_UVC(
            self._unit_x_time_series[frame_number],
            self._unit_y_time_series[frame_number])

        # Electric field magnitude meshplot for this frame.
        self._mesh_plot.set_array(emag_grid)

class AnimFieldExtracellular(AnimField):
    '''
    Animation of the electric field over all extracellular spaces plotted on
    the cell cluster.
    '''

    def __init__(self, *args, **kwargs) -> None:

        # Pass all parameters to our superclass.
        super().__init__(
            is_ecm_required=True,
            *args, **kwargs)

        # Electric field magnitude.
        efield_mag = np.sqrt(
            self._x_time_series[-1] ** 2 + self._y_time_series[-1] ** 2)

        self.msh, self._axes = env_mesh(
            efield_mag, self._axes, self._cells, self._p, self._colormap,
            ignore_showCells=True)
        self.streamE, self._axes = env_quiver(
            self._x_time_series[-1],
            self._y_time_series[-1], self._axes, self._cells, self._p)

        # Autoscale the colorbar range if desired.
        if self._is_color_autoscaled is True:
            self._color_min = np.min(efield_mag)
            self._color_max = np.max(efield_mag)

        #FIXME: How expensive would caching these calculations be? Use
        #attributes defined by our superclass for doing so.

        # Display and/or save this animation. Since recalculating "efield_mag"
        # for each animation frame is non-trivial, this call avoids passing the
        # "time_series" parameter. Instead, the __plot_frame_figure() method
        # manually rescales the colorbar on each frame according to the minimum
        # and maximum electric field magnitude. While non-ideal, every
        # alternative is currently worse.
        self._animate(
            frame_count=len(self._sim.time),
            color_mapping=self.msh,
        )


    def _plot_frame_figure(self, frame_number: int):
        assert types.is_int(frame_number), types.assert_not_int(frame_number)

        E_x = self._x_time_series[frame_number]
        E_y = self._y_time_series[frame_number]

        efield_mag = np.sqrt(E_x**2 + E_y**2)
        self.msh.set_data(efield_mag)

        if efield_mag.max() != 0.0:
            E_x = E_x/efield_mag.max()
            E_y = E_y/efield_mag.max()

        self.streamE.set_UVC(E_x, E_y)

        # Rescale the colorbar range if desired.
        if self._is_color_autoscaled is True:
            self._color_min = np.min(efield_mag)
            self._color_max = np.max(efield_mag)

            #FIXME: Make this go away. A coven of unicycles droven to the edge!

            # Set the colorbar range.
            self.msh.set_clim(self._color_min, self._color_max)

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

        # Pass all parameters *NOT* listed above to our superclass.
        super().__init__(*args, **kwargs)

        # Define this attribute *BEFORE* streamplotting, which assumes this
        # attribute to exist.
        self._stream_plot = None

        # Streamplot the first frame's velocity field.
        vfield, vnorm = self._plot_stream_velocity_field(frame_number=0)

        # Meshplot the first frame's velocity field magnitude.
        self._mesh_plot = self._plot_image(
            pixel_data=vfield,
            colormap=self._p.background_cm,
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
        self._animate(
            frame_count=len(self._sim.time),
            color_mapping=self._mesh_plot,
        )


    def _plot_frame_figure(self, frame_number: int):
        assert types.is_int(frame_number), types.assert_not_int(frame_number)

        # Streamplot this frame's velocity field.
        vfield, vnorm = self._plot_stream_velocity_field(frame_number)

        # Meshplot this frame's velocity field.
        self._mesh_plot.set_data(vfield)

        # Rescale the colorbar if needed.
        self._mesh_plot.set_clim(self._color_min, self._color_max)


    def _plot_stream_velocity_field(self, frame_number: int) -> (
        np.ndarray, float):
        '''
        Streamplot the current velocity field for the passed frame and return a
        2-tuple describing this field.

        Returns
        ----------
        `(velocity_field, velocity_field_magnitude_max)`
            2-element tuple whose:
            * First element is all velocity field magnitudes for this frame.
            * Second element is the maximum such magnitude.
        '''
        assert types.is_int(frame_number), types.assert_not_int(frame_number)

        cell_centres = (
            self._cells.cell_centres[:, 0], self._cells.cell_centres[:, 1])
        cell_grid = (self._cells.X, self._cells.Y)

        u_gj_x = interpolate.griddata(
            cell_centres,
            self._sim.u_cells_x_time[frame_number],
            cell_grid,
            fill_value=0,
            method=self._p.interp_type,
        )
        u_gj_y = interpolate.griddata(
            cell_centres,
            self._sim.u_cells_y_time[frame_number],
            cell_grid,
            fill_value=0,
            method=self._p.interp_type,
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
        if self._is_color_autoscaled is True:
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

        # Pass all parameters *NOT* listed above to our superclass.
        super().__init__(
            is_ecm_required=True,
            *args, **kwargs)

        # Time series of all velocity magnitudes.
        self._magnitude_time_series = np.sqrt(
            np.asarray(self._sim.u_env_x_time) ** 2 +
            np.asarray(self._sim.u_env_y_time) ** 2) * 1e9

        # Velocity field and maximum velocity field value for the first frame.
        vfield = self._magnitude_time_series[0]
        vnorm = np.max(vfield)

        # Velocity field meshplot for the first frame.
        self._mesh_plot = self._plot_image(
            pixel_data=vfield,
            colormap=self._p.background_cm,
        )

        #FIXME: Doesn't this streamplot the last frame instead?

        # Velocity field streamplot for the first frame.
        self._stream_plot = self._axes.quiver(
            self._cells.xypts[:,0] * self._p.um,
            self._cells.xypts[:,1] * self._p.um,
            self._sim.u_env_x_time[-1].ravel() / vnorm,
            self._sim.u_env_y_time[-1].ravel() / vnorm,
        )

        # Display and/or save this animation.
        self._animate(
            frame_count=len(self._sim.time),
            color_mapping=self._mesh_plot,
            color_series=self._magnitude_time_series,
        )


    def _plot_frame_figure(self, frame_number: int):
        assert types.is_int(frame_number), types.assert_not_int(frame_number)

        # Velocity field and maximum velocity field value for this frame.
        vfield = self._magnitude_time_series[frame_number]
        vnorm = np.max(vfield)

        # Update the current velocity meshplot.
        self._mesh_plot.set_data(vfield)

        # Update the current velocity streamplot.
        self._stream_plot.set_UVC(
                self._sim.u_env_x_time[frame_number] / vnorm,
                self._sim.u_env_y_time[frame_number] / vnorm)

# ....................{ SUBCLASSES ~ other                 }....................
class AnimCurrent(AnimCells):
    '''
    Animation of current density plotted on the cell cluster.
    '''

    def __init__(self, *args, **kwargs) -> None:

        # Pass all parameters *NOT* listed above to our superclass.
        super().__init__(*args, **kwargs)

        # Prefer an alternative colormap *BEFORE* plotting below.
        self._colormap = self._p.background_cm

        # Initialize all attributes pertaining to current density.
        self._init_current_density()

        # Current density magnitudes for the first frame in uA/cm2.
        Jmag_M = self._current_density_magnitude_time_series[0]

        # Streamplot the first frame's current density, classified to permit
        # erasure of this streamplot by the _replot_current_density() method.
        self._current_density_stream_plot = self._plot_stream(
            x=self._current_density_x_time_series[0] / Jmag_M,
            y=self._current_density_y_time_series[0] / Jmag_M,
            magnitude=Jmag_M,
        )

        # Meshplot the first frame's current density magnitude.
        self._mesh_plot = self._plot_image(
            pixel_data=Jmag_M,
            colormap=self._p.background_cm,
        )

        # Display and/or save this animation.
        self._animate(
            frame_count=len(self._sim.time),
            color_mapping=self._mesh_plot,
            color_series=self._current_density_magnitude_time_series,
        )


    def _plot_frame_figure(self, frame_number: int):
        assert types.is_int(frame_number), types.assert_not_int(frame_number)

        # Current density magnitudes for this frame in uA/cm2.
        Jmag_M = self._current_density_magnitude_time_series[frame_number]

        # Streamplot this frame's current density.
        self._replot_current_density(frame_number)

        # Meshplot this frame's current density magnitude.
        self._mesh_plot.set_data(Jmag_M)


#FIXME: Use below in lieu of string constants.
AnimDeformStyle = Enum('AnimDeformStyle', ('STREAMLINE', 'VECTOR'))

#FIXME: Reenable after we deduce why the "AnimDeform" class defined
#below no longer animates physical displacement. Since the
#AnimCellsWhileSolving.resetData() method *DOES* appear to animate physical
#displacement, that's probably the first place to start. The key appears to be
#completely recreating the entire plot -- but then, don't we do that here?
#FIXME: Split into two subclasses: one handling physical deformations and the
#other voltage deformations. There exists very little post-subclassing code
#sharred in common between the two conditional branches handling this below.

class AnimDeformTimeSeries(AnimCells):
    '''
    Animation of physical cell deformation overlayed an arbitrary cell-centric
    time series (e.g., cell voltage as a function of time) on the cell cluster.

    Attributes
    ----------
    _cell_time_series : list
        Arbitrary cell data as a function of time to be underlayed.
    '''

    def __init__(
        self,
        cell_time_series: list,
        *args, **kwargs
    ) -> None:
        '''
        Initialize this animation.

        Parameters
        ----------
        cell_time_series : list
            Arbitrary cell data as a function of time to be underlayed.

        See the superclass `__init__()` method for all remaining parameters.
        '''
        assert types.is_sequence_nonstr(cell_time_series), (
            types.assert_not_sequence_nonstr(cell_time_series))

        # Pass all parameters *NOT* listed above to our superclass.
        super().__init__(*args, **kwargs)

        # Classify all remaining parameters.
        self._cell_time_series = cell_time_series

        # Cell displacement magnitudes for the first frame.
        dd = self._cell_time_series[0]

        # Cell displacement X and Y components for the first frame.
        dx = self._sim.dx_cell_time[0]
        dy = self._sim.dy_cell_time[0]

        if self._p.showCells is True:
            dd_collection, self._axes = cell_mosaic(
                dd, self._axes, self._cells, self._p, self._colormap)
        else:
            dd_collection, self._axes = cell_mesh(
                dd, self._axes, self._cells, self._p, self._colormap)

        if self._p.ani_Deformation_style == 'vector':
            self._quiver_plot, self._axes = cell_quiver(
                dx, dy, self._axes, self._cells, self._p)
        elif self._p.ani_Deformation_style == 'streamline':
            self._stream_plot, self._axes = cell_stream(
                dx, dy, self._axes, self._cells, self._p,
                showing_cells=False)
        elif self._p.ani_Deformation_style != 'None':
            raise BetseExceptionParameters(
                'Deformation animation style "{}" not '
                '"vector", "streamline", or "None".'.format(
                    self._p.ani_Deformation_style))

        # Display and/or save this animation.
        self._animate(
            frame_count=len(self._sim.time),
            color_mapping=dd_collection,
            color_series=self._cell_time_series,
        )


    #FIXME: Quite a bit of code duplication. Generalize us up the bomb.
    def _plot_frame_figure(self, frame_number: int):
        assert types.is_int(frame_number), types.assert_not_int(frame_number)

        # Erase the prior frame before plotting this frame.
        # self.ax.patches = []

        #FIXME: Improve to avoid clearing the current axes.

        # we need to have changing cells, so we have to clear the plot and redo
        # it...
        # self.fig.clf()
        # self.ax = plt.subplot(111)
        self._axes.cla()
        self._axes.axis('equal')
        self._axes.axis(self._axes_bounds)
        self._axes.set_xlabel('Spatial distance [um]')
        self._axes.set_ylabel('Spatial distance [um]')

        # Arrays of all cell deformation X and Y components for this frame.
        dx = self._sim.dx_cell_time[frame_number]
        dy = self._sim.dy_cell_time[frame_number]

        # Array of all cell Vmem values for this frame.
        if self._p.ani_Deformation_type == 'Vmem':
            if self._p.sim_ECM is False:
                dd = self._sim.vm_time[frame_number] * 1e3
            else:
                dd = self._sim.vcell_time[frame_number] * 1e3
        # Array of all cell deformation magnitudes for this frame.
        elif self._p.ani_Deformation_type == 'Displacement':
            dd = 1e6 * np.sqrt(dx**2 + dy**2)

        # Reset the superclass colorbar mapping to this newly created mapping,
        # permitting the superclass plot_frame() method to clip this mapping.
        # dd_collection.remove()
        # self.ax.collections = []
        if self._p.showCells is True:
            dd_collection, self._axes = cell_mosaic(
                dd, self._axes, self._cells, self._p, self._colormap)
            # points = np.multiply(self.cells.cell_verts, self.p.um)
            # dd_collection = PolyCollection(
            #     points, cmap=self.colormap, edgecolors='none')
            # dd_collection.set_array(dd)
        else:
            dd_collection, self._axes = cell_mesh(
                dd, self._axes, self._cells, self._p, self._colormap)

        dd_collection.set_clim(self._color_min, self._color_max)
        # cb = self.fig.colorbar(dd_collection)
        # cb.set_label(self._colorbar_title)

        if self._p.ani_Deformation_style == 'vector':
            # self._quiver_plot.remove()
            quiver_plot, self._axes = cell_quiver(
                dx, dy, self._axes, self._cells, self._p)
        elif self._p.ani_Deformation_style == 'streamline':
            # self._stream_plot.lines.remove()
            stream_plot, self._axes = cell_stream(
                dx, dy, self._axes, self._cells, self._p,
                showing_cells=self._p.showCells)


class AnimateDeformation(object):

    def __init__(
        self,
        sim,
        cells,
        p,
        ani_repeat=True,
        save=True,
        saveFolder='animation/Deformation',
        saveFile='Deformation_',
    ):

        self.fig = plt.figure()
        self.ax = plt.subplot(111)
        self.p = p
        self.sim = sim
        self.cells = cells
        self.save = save

        self.saveFolder = saveFolder
        self.saveFile = saveFile
        self.ani_repeat = ani_repeat

        if self.save is True:
            _setup_file_saving(self,p)

        dx = self.sim.dx_cell_time[0]
        dy = self.sim.dy_cell_time[0]

        if self.p.ani_Deformation_type == 'Vmem':
            if self.p.sim_ECM is False:
                dd = self.sim.vm_time[0]*1e3
            else:
                dd = self.sim.vcell_time[0]*1e3

            self.specific_cmap = p.default_cm
        elif self.p.ani_Deformation_type == 'Displacement':
            dd = p.um*np.sqrt(dx**2 + dy**2)
            self.specific_cmap = p.background_cm
        else:
            raise BetseExceptionParameters(
                "Definition of 'data type' in deformation animation\n"
                "must be either 'Vmem' or 'Displacement'.")

        dd_collection, self.ax = cell_mesh(
            dd,self.ax,cells,p,self.specific_cmap)

        if p.ani_Deformation_style == 'vector':
            vplot, self.ax = cell_quiver(dx,dy,self.ax,cells,p)
        elif p.ani_Deformation_style == 'streamline':
            vplot, self.ax = cell_stream(
                dx,dy,self.ax,cells,p, showing_cells=p.showCells)
        else:
            raise BetseExceptionParameters(
                "Definition of 'style' in deformation animation\n"
                "must be either 'vector' or 'streamline'.")

        self.ax.axis('equal')

        xmin = cells.xmin*p.um
        xmax = cells.xmax*p.um
        ymin = cells.ymin*p.um
        ymax = cells.ymax*p.um

        self.ax.axis([xmin,xmax,ymin,ymax])

        if p.autoscale_Deformation_ani is True:
            if p.ani_Deformation_type == 'Displacement':
                # first flatten the data (needed in case cells were cut)
                all_z = []
                for xarray, yarray in zip(sim.dx_cell_time,sim.dy_cell_time):
                    zarray = np.sqrt(xarray**2 + yarray**2)
                    for val in zarray:
                        all_z.append(val*p.um)

            elif p.ani_Deformation_type == 'Vmem':
                all_z = []

                for zarray in sim.vm_time:
                    for val in zarray:
                        all_z.append(val*1e3)

            self.cmin = np.min(all_z)
            self.cmax = np.max(all_z)

            dd_collection.set_clim(self.cmin,self.cmax)

        elif p.autoscale_Deformation_ani is False:
            dd_collection.set_clim(
                p.Deformation_ani_min_clr, p.Deformation_ani_max_clr)

        cb = self.fig.colorbar(dd_collection)

        self.tit = "Displacement Field and Deformation"
        self.ax.set_title(self.tit)
        self.ax.set_xlabel('Spatial distance [um]')
        self.ax.set_ylabel('Spatial distance [um]')

        if self.p.ani_Deformation_type == 'Displacement':
            cb.set_label('Displacement [um]')

        elif self.p.ani_Deformation_type == 'Vmem':
            cb.set_label('Voltage [mV]')

        #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! BEGIN !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        # NOTE: This is the only new code that has been added to this class.

        # If animation saving is enabled, prepare to do so.
        if self.p.saveAnimations is True:
            self._type = 'Deformation'

            # Ensure that the passed directory and file basenames are actually
            # basenames and hence contain no directory separators.
            paths.die_unless_basename(self._type)

            # Path of the phase-specific parent directory of the subdirectory to
            # which these files will be saved.
            phase_dirname = None
            if self.p.plot_type == 'sim':
                phase_dirname = self.p.sim_results
            elif self.p.plot_type == 'init':
                phase_dirname = self.p.init_results
            else:
                raise BetseExceptionParameters(
                    'Anim saving unsupported during the "{}" phase.'.format(
                        self.p.plot_type))

            # Path of the subdirectory to which these files will be saved, creating
            # this subdirectory and all parents thereof if needed.
            save_dirname = paths.join(
                phase_dirname, 'animation', self._type)
            save_dirname = dirs.canonicalize_and_make_unless_dir(save_dirname)
            save_frame_filetype = 'png'

            # Template yielding the basenames of frame image files to be saved.
            # The "{{"- and "}}"-delimited substring will reduce to a "{"- and "}"-
            # delimited substring after formatting, which subsequent formatting
            # elsewhere (e.g., in the "FileFrameWriter" class) will expand with the
            # 0-based index of the current frame number.
            save_frame_template_basename = '{}_{{:07d}}.{}'.format(
                self._type, save_frame_filetype)

            # Template yielding the absolute paths of frame image files to be saved.
            self._save_frame_template = paths.join(
                save_dirname, save_frame_template_basename)

            # Object writing frames from this animation to image files.
            self._writer_frames = FileFrameWriter()
        #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! END   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        self.frames = len(sim.time)
        ani = animation.FuncAnimation(self.fig, self.aniFunc,
            frames=self.frames, interval=100, repeat=self.ani_repeat)

        #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! BEGIN !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        # NOTE: This is the only new code that has been added to this class.
        try:
            if p.turn_all_plots_off is False:
                plt.show()
            # Else if saving animation frames, do so.
            elif self.p.saveAnimations is True:
                logs.log_info(
                    'Saving animation "{}" frames...'.format(self._type))
                ani.save(
                    filename=self._save_frame_template,
                    writer=self._writer_frames)
        # plt.show() unreliably raises exceptions on window close resembling:
        #     AttributeError: 'NoneType' object has no attribute 'tk'
        # This error appears to ignorable and hence is caught and squelched.
        except AttributeError as exc:
            # If this is that exception, mercilessly squelch it.
            if str(exc) == "'NoneType' object has no attribute 'tk'":
                pass
            # Else, reraise this exception.
            else:
                raise
        #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! END   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    def aniFunc(self,i):

        # we need to have changing cells, so we have to clear the plot and redo it...
        self.fig.clf()
        self.ax = plt.subplot(111)

        dx = self.sim.dx_cell_time[i]
        dy = self.sim.dy_cell_time[i]

        if self.p.ani_Deformation_type == 'Vmem':
            if self.p.sim_ECM is False:
                dd = self.sim.vm_time[i]*1e3
            else:
                dd = self.sim.vcell_time[i]*1e3
        elif self.p.ani_Deformation_type == 'Displacement':
            dd = 1e6*np.sqrt(dx**2 + dy**2)

        dd_collection, self.ax = cell_mesh(
            dd, self.ax, self.cells, self.p, self.specific_cmap)

        if self.p.ani_Deformation_style == 'vector':
            vplot, self.ax = cell_quiver(dx,dy,self.ax,self.cells,self.p)
        elif self.p.ani_Deformation_style == 'streamline':
            vplot, self.ax = cell_stream(
                dx, dy, self.ax, self.cells, self.p,
                showing_cells=self.p.showCells)

        self.ax.axis('equal')

        xmin = self.cells.xmin*self.p.um
        xmax = self.cells.xmax*self.p.um
        ymin = self.cells.ymin*self.p.um
        ymax = self.cells.ymax*self.p.um

        self.ax.axis([xmin,xmax,ymin,ymax])

        if self.p.autoscale_Deformation_ani is False:
            dd_collection.set_clim(
                self.p.Deformation_ani_min_clr,self.p.Deformation_ani_max_clr)
        else:
            dd_collection.set_clim(self.cmin,self.cmax)

        titani = self.tit + ' (simulation time' + ' ' + str(round(self.sim.time[i],3)) + ' ' + ' s)'
        self.ax.set_title(titani)
        self.ax.set_xlabel('Spatial distance [um]')
        self.ax.set_ylabel('Spatial distance [um]')

        cb = self.fig.colorbar(dd_collection)

        if self.p.ani_Deformation_type == 'Displacement':
            cb.set_label('Displacement [um]')
        elif self.p.ani_Deformation_type == 'Vmem':
            cb.set_label('Voltage [mV]')

        if self.save is True:
            self.fig.canvas.draw()
            savename = self.savedAni + str(i) + '.png'
            plt.savefig(savename,format='png')
