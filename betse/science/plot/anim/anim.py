#!/usr/bin/env python3
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Matplotlib-based animation classes.
'''

#FIXME: Replace all instances of the "current_overlay" and "is_current_overlay"
#booleans in this class with the corresponding superclass functionality.

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

# ....................{ IMPORTS                            }....................
import os
import numpy as np
from betse.exceptions import BetseExceptionParameters
from betse.lib.matplotlib import mpl
from betse.lib.matplotlib.anim import FileFrameWriter
from betse.science.plot.anim.abc import (
    AnimCells, AnimCellsField, AnimCellsVelocity)
from betse.util.io import loggers
from betse.util.path import dirs, paths
from betse.util.type import types
from enum import Enum
from matplotlib import pyplot as plt
from matplotlib import animation
from matplotlib.collections import LineCollection, PolyCollection
from numpy import ma as ma
from scipy import interpolate

#FIXME: Shift functions called only by this module either to a new
#"betse.science.plot.animation.helper" module or possibly as private
#methods of the "Anim" superclass. Right. After investigation, absolutely the
#latter approach. This should permit us to avoid passing *ANY* parameters to
#these methods, which is rather nice.
from betse.science.plot.plot import (
    _setup_file_saving, _handle_plot, env_mesh, cell_mosaic, cell_mesh,
    I_overlay_setup, I_overlay_update, env_quiver, cell_quiver, cell_stream
)

# ....................{ CLASSES                            }....................
#FIXME: The "is_current_overlay" boolean is handled identically here as it is
#in the "AnimateMem" class. This suggests that the existing
#plot.I_overlay_setup() and plot.I_overlay_update() functions could be shifted
#into a new superclass of this class handling optional current overlays.

class AnimateCellData(AnimCells):
    '''
    Animation of arbitrary colour data plotted on the current cell cluster.

    Attributes
    ----------
    _is_ecm_ignored : bool
        `True` if ignoring extracellular spaces _or_ `False` otherwise.
    _time_series : list
        List of all colour data to be plotted, indexed by simulation time.
    '''

    def __init__(
        self,
        time_series: list,
        is_ecm_ignored: bool = True,
        *args, **kwargs
    ) -> None:
        '''
        Initialize this animation.

        Parameters
        ----------
        time_series : list
            List of all colour data to be plotted, indexed by simulation time.
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
            axes_x_label='Spatial x [um]',
            axes_y_label='Spatial y [um]',

            # Since this class does *NOT* plot a streamplot, request that the
            # superclass do so for electric current or concentration flux.
            is_plotting_current_overlay=True,
            *args, **kwargs
        )

        # Classify parameters required by the _plot_next_frame() method.
        self._time_series = time_series
        self._is_ecm_ignored = is_ecm_ignored

        data_points = self._time_series[0]

        if self._p.sim_ECM is True and is_ecm_ignored is False:
            self.collection, self._axes = env_mesh(
                data_points, self._axes, self._cells, self._p, self._colormap)
        elif self._p.showCells is True:
            self.collection, self._axes = cell_mosaic(
                data_points, self._axes, self._cells, self._p, self._colormap)
        else:
            self.collection, self._axes = cell_mesh(
                data_points, self._axes, self._cells, self._p, self._colormap)

        # Display and/or save this animation.
        self._animate(
            frame_count=len(self._sim.time),
            color_mapping=self.collection,
            color_series=self._time_series,
        )


    def _plot_frame_figure(self, frame_number: int):
        assert types.is_int(frame_number), types.assert_not_int(frame_number)

        zz = self._time_series[frame_number]

        if self._p.sim_ECM is True and self._is_ecm_ignored is False:
            if self._p.plotMask is True:
                dat_grid = ma.masked_array(
                    zz, np.logical_not(self._cells.maskM))

            self.collection.set_data(dat_grid)
        else:
            if self._p.showCells is True:
                self.collection.set_array(zz)
            else:
                zz_grid = np.zeros(len(self._cells.voronoi_centres))
                zz_grid[self._cells.cell_to_grid] = zz
                self.collection.set_array(zz_grid)


class AnimateCurrent(AnimCells):
    '''
    Animation of current plotted on the current cell cluster.

    Attributes
    ----------
    _current_density_magnitude_time_series : ndarray
        Time series of all current density magnitudes (i.e., `Jmag_M`).
    _current_density_x_time_series : list
        Time series of all current density X components.
    _current_density_y_time_series : list
        Time series of all current density Y components.
    '''

    def __init__(
        self,
        is_gj_current_only: bool,
        *args, **kwargs
    ):
        '''
        Initialize this animation.

        Parameters
        ----------
        is_gj_current_only : bool
            `True` if animating only the gap junction current _or_ `False` if
            animating all current.

        See the superclass `__init__()` method for all remaining parameters.
        '''
        assert types.is_bool(is_gj_current_only), (
            types.assert_not_bool(is_gj_current_only))

        # Pass all parameters *NOT* listed above to our superclass.
        super().__init__(
            axes_x_label='Spatial x [um]',
            axes_y_label='Spatial y [um]',
            *args, **kwargs)

        self._colormap = self._p.background_cm

        # Time series of all current density X and Y components.
        if is_gj_current_only is True:
            self._current_density_x_time_series = self._sim.I_gj_x_time
            self._current_density_y_time_series = self._sim.I_gj_y_time
        else:
            self._current_density_x_time_series = self._sim.I_tot_x_time
            self._current_density_y_time_series = self._sim.I_tot_y_time

        # Time series of all current density magnitudes (i.e., `Jmag_M`).
        self._current_density_magnitude_time_series = np.sqrt(
            np.array(self._current_density_x_time_series) ** 2 +
            np.array(self._current_density_y_time_series) ** 2) + 1e-30

        # Stream- and meshplot the first frame's current density magnitude.
        Jmag_M = self._plot_stream_current_density_magnitude(frame_number=0)
        self._mesh_plot = self._plot_image(Jmag_M)

        # Display and/or save this animation.
        self._animate(
            frame_count=len(self._sim.time),
            color_mapping=self._mesh_plot,
            color_series=self._current_density_magnitude_time_series,
        )


    def _plot_frame_figure(self, frame_number: int):
        assert types.is_int(frame_number), types.assert_not_int(frame_number)

        # Erase the prior frame's streamplot before streamplotting this frame.
        self._axes.patches = []
        self._stream_plot.lines.remove()

        # Stream- and meshplot the current density magnitude for this frame.
        Jmag_M = self._plot_stream_current_density_magnitude(frame_number)
        self._mesh_plot.set_data(Jmag_M)


    def _plot_stream_current_density_magnitude(
        self, frame_number: int) -> np.ndarray:
        '''
        Streamplot and return all current density magnitudes (i.e., `Jmag_M`)
        for the passed frame.
        '''
        assert types.is_int(frame_number), types.assert_not_int(frame_number)

        # All current density magnitudes for the current frame.
        Jmag_M = self._current_density_magnitude_time_series[frame_number]

        #FIXME: Are streamplots always passed unit vectors? Sparkling cider!

        # Classify this streamplot, thus permitting the _plot_frame_figure()
        # method to subsequently erase this streamplot's lines.
        self._stream_plot = self._plot_stream(
            x=self._current_density_x_time_series[frame_number] / Jmag_M,
            y=self._current_density_y_time_series[frame_number] / Jmag_M,
            magnitude=Jmag_M,
        )

        return Jmag_M


#FIXME: Use below in lieu of string constants.
AnimDeformationStyle = Enum('AnimDeformationStyle', ('streamline', 'vector'))

#FIXME: Reenable after we deduce why the "AnimateDeformation" class defined
#below no longer animates physical displacement.
#FIXME: Split into two subclasses: one handling physical deformations and the
#other voltage deformations. There exists very little post-subclassing code
#sharred in common between the two conditional branches handling this below.

class AnimateDeformationUnused(AnimCells):
    '''
    Animation of physical cell deformation plotted on the current cell cluster.
    '''

    def __init__(self, *args, **kwargs) -> None:

        # Pass all parameters *NOT* listed above to our superclass.
        super().__init__(
            axes_x_label='Spatial distance [um]',
            axes_y_label='Spatial distance [um]',
            *args, **kwargs)

        dx = self._sim.dx_cell_time[0]
        dy = self._sim.dy_cell_time[0]

        if self._p.ani_Deformation_type == 'Vmem':
            self._colormap = self._p.default_cm

            if self._p.sim_ECM is False:
                dd = self._sim.vm_time[0] * 1e3
            else:
                dd = self._sim.vcell_time[0] * 1e3

        elif self._p.ani_Deformation_type == 'Displacement':
            self._colormap = self._p.background_cm
            dd = self._p.um * np.sqrt(dx ** 2 + dy ** 2)

        else:
            raise BetseExceptionParameters(
                'Deformation animation type "{}" not '
                '"Vmem" or "Displacement".'.format(
                    self._p.ani_Deformation_type))

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

        # Sequence of all deformation values for use in colorbar autoscaling.
        colorbar_time_series = None

        #FIXME: Classify colorbar_time_series as a more appropriately named
        #attribute (e.g., "self._deform_time_series"); then, use that attribute
        #in _plot_frame_figure() rather than recompute such values.
        if self._is_color_autoscaled is True:
            if self._p.ani_Deformation_type == 'Displacement':
                #FIXME: Optimize. Since our superclass already ravels, this
                #should reduce to just:
                #colorbar_time_series = np.array([
                #     np.sqrt(cell_dx**2 + cell_dy**2) * self.p.um
                #     for cell_dx, cell_dy in zip(
                #        self.sim.dx_cell_time, self.sim.dy_cell_time)])

                # first flatten the data (needed in case cells were cut)
                colorbar_time_series = []

                # For sequence of X and Y cell deformation components for each
                # time step, ...
                for cell_dx, cell_dy in zip(
                    self._sim.dx_cell_time, self._sim.dy_cell_time):
                    zarray = np.sqrt(cell_dx**2 + cell_dy**2)
                    for val in zarray:
                        colorbar_time_series.append(val * self._p.um)

            elif self._p.ani_Deformation_type == 'Vmem':
                #FIXME: Optimize. Since our superclass already ravels, this
                #should reduce to just:
                #    colorbar_time_series = self.sim.vm_time * 1e3

                colorbar_time_series = []

                for zarray in self._sim.vm_time:
                    for val in zarray:
                        colorbar_time_series.append(val*1e3)

        # Display and/or save this animation.
        self._animate(
            frame_count=len(self._sim.time),
            color_mapping=dd_collection,
            color_series=colorbar_time_series,
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
        # permitting the superclass _plot_frame() method to clip this mapping.
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
                loggers.log_info(
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


class AnimateFieldIntracellular(AnimCellsField):
    '''
    Animation of the electric field over all intracellular gap junctions
    plotted on the current cell cluster.
    '''

    def __init__(self, *args, **kwargs) -> None:

        # Pass all parameters *NOT* listed above to our superclass.
        super().__init__(*args, **kwargs)

        # Electric field magnitude.
        efield_mag = np.sqrt(self._Fx_time[-1]**2 + self._Fy_time[-1]**2)

        self.msh, self._axes = cell_mesh(
            efield_mag, self._axes, self._cells, self._p, self._colormap)
        self.streamE, self._axes = cell_quiver(
            self._Fx_time[-1], self._Fy_time[-1], self._axes, self._cells, self._p)

        # Autoscale the colorbar range if desired.
        if self._is_color_autoscaled is True:
            self._color_min = np.min(efield_mag)
            self._color_max = np.max(efield_mag)

        #FIXME: How expensive would caching these calculations be?

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

        E_gj_x = self._Fx_time[frame_number]
        E_gj_y = self._Fy_time[frame_number]

        if len(E_gj_x) != len(self._cells.cell_i):
            E_gj_x = (
                np.dot(self._cells.M_sum_mems, E_gj_x) / self._cells.num_mems)
            E_gj_y = (
                np.dot(self._cells.M_sum_mems, E_gj_y) / self._cells.num_mems)

        efield_mag = np.sqrt(E_gj_x**2 + E_gj_y**2)
        emag_grid = np.zeros(len(self._cells.voronoi_centres))
        emag_grid[self._cells.cell_to_grid] = efield_mag
        self.msh.set_array(emag_grid)

        if efield_mag.all() != 0.0:
            E_gj_x = E_gj_x/efield_mag
            E_gj_y = E_gj_y/efield_mag

        self.streamE.set_UVC(E_gj_x, E_gj_y)

        # Rescale the colorbar range if desired.
        if self._is_color_autoscaled is True:
            self._color_min = np.min(efield_mag)
            self._color_max = np.max(efield_mag)

            #FIXME: Make this go away. A coven of unicycles droven to the edge!

            # Set the colorbar range.
            self.msh.set_clim(self._color_min, self._color_max)


class AnimateFieldExtracellular(AnimCellsField):
    '''
    Animation of the electric field over all extracellular spaces plotted on
    the current cell cluster.
    '''

    def __init__(self, *args, **kwargs) -> None:

        # Pass all parameters to our superclass.
        super().__init__(*args, **kwargs)

        # Validate sanity.
        if self._p.sim_ECM is False:
            raise BetseExceptionParameters(
                'Electric field animation "{}" plotted over '
                'extracellular spaces, but '
                'extracellular spaces are disabled by the '
                'current simulation configuration.'.format(
                self._type))

        # Electric field magnitude.
        efield_mag = np.sqrt(self._Fx_time[-1]**2 + self._Fy_time[-1]**2)

        self.msh, self._axes = env_mesh(
            efield_mag, self._axes, self._cells, self._p, self._colormap,
            ignore_showCells=True)
        self.streamE, self._axes = env_quiver(
            self._Fx_time[-1], self._Fy_time[-1], self._axes, self._cells, self._p)

        # Autoscale the colorbar range if desired.
        if self._is_color_autoscaled is True:
            self._color_min = np.min(efield_mag)
            self._color_max = np.max(efield_mag)

        #FIXME: How expensive would caching these calculations be?

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

        E_x = self._Fx_time[frame_number]
        E_y = self._Fy_time[frame_number]

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


class AnimateGJData(AnimCells):
    '''
    Animation of the gap junction open state as a function of time.

    Attributes
    ----------
    _cell_time_series : list
        Time series used for cell coloring.
    _gj_time_series : list
        Time series used for gap junction coloring.
    '''

    def __init__(
        self,
        gj_time_series: list,
        cell_time_series: list,
        *args, **kwargs
    ) -> None:
        '''
        Initialize this animation.

        Parameters
        ----------
        gj_time_series : list
            Time series used for gap junction coloring.
        cell_time_series : list
            Time series used for cell coloring.

        See the superclass `__init__()` method for all remaining parameters.
        '''
        assert types.is_sequence_nonstr(gj_time_series), (
            types.assert_not_sequence_nonstr(gj_time_series))
        assert types.is_sequence_nonstr(cell_time_series), (
            types.assert_not_sequence_nonstr(cell_time_series))

        # Pass all parameters *NOT* listed above to our superclass.
        super().__init__(
            axes_x_label='Spatial x [um]',
            axes_y_label='Spatial y [um]',
            *args, **kwargs)

        # Classify all remaining parameters.
        self._gj_time_series = gj_time_series
        self._cell_time_series = cell_time_series

        connects = np.asarray(self._cells.nn_edges) * self._p.um
        self.collection = LineCollection(
            connects,
            array=self._gj_time_series[0],
            cmap=self._p.gj_cm,
            linewidths=2.0,
            zorder=10,
        )
        self.collection.set_clim(0.0, 1.0)
        self._axes.add_collection(self.collection)

        # Add a collection of cell polygons with animated voltage data.
        if self._p.sim_ECM is False:
            data_set = self._cell_time_series[0]
        else:
            data_set = self._sim.vcell_time[0] * 1000

        if self._p.showCells is True:
            self.coll2, self._axes = cell_mosaic(
                data_set, self._axes, self._cells, self._p, self._colormap)
        else:
            self.coll2, self._axes = cell_mesh(
                data_set, self._axes, self._cells, self._p, self._colormap)

        # Display and/or save this animation.
        self._animate(
            frame_count=len(self._gj_time_series),
            color_mapping=self.coll2,

            #FIXME: If modelling extracellular spaces, this doesn't seem quite
            #right. In that case, shouldn't this be something resembling:
            #    color_series=self.sim.vcell_time*1000,
            #Tug boats in the muggy bastions of the gentle night!

            color_series=self._cell_time_series,
        )


    def _plot_frame_figure(self, frame_number: int):
        assert types.is_int(frame_number), types.assert_not_int(frame_number)

        zz = self._gj_time_series[frame_number]
        self.collection.set_array(zz)
        self.collection.set_clim(0.0, 1.0)

        if self._p.sim_ECM is False:
            zv = self._cell_time_series[frame_number]
        else:
            zv = self._sim.vcell_time[frame_number] * 1000

        if self._p.showCells is True:
            zz_grid = zv
        else:
            zz_grid = np.zeros(len(self._cells.voronoi_centres))
            zz_grid[self._cells.cell_to_grid] = zv

        self.coll2.set_array(zz_grid)


class AnimateVelocityIntracellular(AnimCellsVelocity):
    '''
    Animation of fluid velocity over all intracellular gap junctions plotted on
    the current cell cluster.
    '''

    def __init__(self, *args, **kwargs) -> None:

        # Pass all parameters *NOT* listed above to our superclass.
        super().__init__(*args, **kwargs)

        # Velocity field and maximum velocity field value for the first frame.
        vfield, vnorm = self._get_velocity_field(frame_number=0)

        # Vector field meshplot for the first frame.
        self._mesh_plot = self._plot_image(vfield)

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

        # Erase the prior frame's streamplot before streamplotting this frame.
        self._axes.patches = []
        self._stream_plot.lines.remove()

        # Velocity field and maximum velocity field value for this frame.
        vfield, vnorm = self._get_velocity_field(frame_number)

        # Update the current velocity field mesh.
        self._mesh_plot.set_data(vfield)

        # Rescale the colorbar if required.
        self._mesh_plot.set_clim(self._color_min, self._color_max)


    def _get_velocity_field(self, frame_number: int) -> tuple:
        '''
        Get a 2-element tuple whose first element is the velocity field and
        second element is the maximum value of this field for the current
        frame.

        For convenience, the streamplot of this field will also be redefined.

        Returns
        ----------
        `(velocity_field, velocity_field_max_value)`
            2-element tuple as described above.
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

        vfield = np.sqrt(u_gj_x**2 + u_gj_y**2)*1e9
        vnorm = np.max(vfield)

        # Replot the streamplot.
        self._stream_plot = self._plot_stream(
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


class AnimateVelocityExtracellular(AnimCellsVelocity):
    '''
    Animation of fluid velocity over all extracellular spaces plotted on the
    current cell cluster.

    Attributes
    ----------
    _velocity_magnitude_time_series : np.ndarray
        Time series of all fluid velocity magnitudes.
    '''

    def __init__(self, *args, **kwargs) -> None:

        # Pass all parameters *NOT* listed above to our superclass.
        super().__init__(*args, **kwargs)

        # Validate sanity.
        if self._p.sim_ECM is False:
            raise BetseExceptionParameters(
                'Extracellular spaces fluid velocity animation "{}" '
                'requested, but extracellular spaces are disabled by the '
                'current simulation configuration.'.format(
                self._type))

        # Time series of all velocity magnitudes.
        self._velocity_magnitude_time_series = np.sqrt(
            np.array(self._sim.u_env_x_time) ** 2 +
            np.array(self._sim.u_env_y_time) ** 2) * 1e9

        # Velocity field and maximum velocity field value for the first frame.
        vfield = self._velocity_magnitude_time_series[0]
        vnorm = np.max(vfield)

        # Velocity field meshplot for the first frame.
        self._mesh_plot = self._plot_image(vfield)

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
            color_series=self._velocity_magnitude_time_series,
        )


    def _plot_frame_figure(self, frame_number: int):
        assert types.is_int(frame_number), types.assert_not_int(frame_number)

        # Velocity field and maximum velocity field value for this frame.
        vfield = self._velocity_magnitude_time_series[frame_number]
        vnorm = np.max(vfield)

        # Update the current velocity meshplot.
        self._mesh_plot.set_data(vfield)

        # Update the current velocity streamplot.
        self._stream_plot.set_UVC(
                self._sim.u_env_x_time[frame_number] / vnorm,
                self._sim.u_env_y_time[frame_number] / vnorm)


class AnimateEnv(object):
# class AnimateEnv(Anim):

    def __init__(
        self,
        sim,
        cells,
        time,
        p,
        save=True,
        ani_repeat=False,
        clrAutoscale=True,
        clrMin=None,
        clrMax=None,
        clrmap=mpl.get_colormap('rainbow'),
        number_cells=False,
        saveFolder='animation/Venv',
        saveFile='venv_',
    ):

        self.clrmap = clrmap
        self.time = time
        self.save = save

        self.sim = sim

        self.sim_ECM = p.sim_ECM
        self.IecmPlot = p.IecmPlot

        self.cells = cells
        self.p = p

        self.fig = plt.figure()       # define figure
        self.ax = plt.subplot(111)    # define axes

        self.ax.axis('equal')

        self.saveFolder = saveFolder
        self.saveFile = saveFile
        self.ani_repeat = ani_repeat

        xmin = cells.xmin*p.um
        xmax = cells.xmax*p.um
        ymin = cells.ymin*p.um
        ymax = cells.ymax*p.um

        self.ax.axis([xmin,xmax,ymin,ymax])

        if self.save is True:
            _setup_file_saving(self,p)

        if clrAutoscale is False:
            self.cmin = clrMin
            self.cmax = clrMax

        self.meshplot = plt.imshow(
            sim.venv_time[0].reshape(cells.X.shape)*1000,
            origin='lower',
            extent=[xmin,xmax,ymin,ymax],
            cmap=p.default_cm,
        )

        if clrAutoscale is False:
            self.meshplot.set_clim(self.cmin,self.cmax)

        self.cb = self.fig.colorbar(self.meshplot)   # define colorbar for figure
        self.cb.set_label('Voltage [V]')

        self.ax.set_xlabel('Spatial x [um]')
        self.ax.set_ylabel('Spatial y [um]')
        self.ax.set_title('Environmental Voltage')

        self.frames = len(sim.time)
        ani = animation.FuncAnimation(self.fig, self.aniFunc,
            frames=self.frames, interval=100, repeat=self.ani_repeat)

        _handle_plot(p)


    def aniFunc(self,i):

        titani = 'Environmental Voltage' + ' (simulation time' + ' ' + str(round(self.sim.time[i],3)) + ' ' + ' s)'
        self.ax.set_title(titani)

        self.meshplot.set_data(self.sim.venv_time[i].reshape(self.cells.X.shape)*1000)

        if self.save is True:
            self.fig.canvas.draw()
            savename = self.savedAni + str(i) + '.png'
            plt.savefig(savename,format='png')


class AnimateMem(object):
# class AnimateMem(object):
    '''
    Animates the channel or pump density factor (`sim.rho_channel` or
    `sim.rho_pump` respectively) which changes due to electroosmotic and
    electrophoretic movements produced by self-generated fields and flows in
    the cluster.
    '''

    def __init__(
        self,
        sim,
        cells,
        time,
        p,
        save=False,
        ani_repeat=False,
        current_overlay=False,
        clrAutoscale=True,
        clrMin=None,
        clrMax=None,
        number_cells=False,
        saveFolder='animation/pump_electroosmo',
        saveFile='rhoPump_',
    ):

        self.colormap = p.default_cm
        self.time = time
        self.save = save

        self.cbtit = 'mol fraction/m2'
        self.cells = cells
        self.p = p

        self.fig = plt.figure()       # define figure
        self.ax = plt.subplot(111)    # define axes

        self.sim = sim
        self.current_overlay = current_overlay

        #FIXME: Is this a typo? "self.colormap" is already defined above to the
        #same value. Unicorns frolicking in the grape-eyed forest!
        self.clrmap = p.default_cm

        self.sim_ECM = p.sim_ECM
        self.IecmPlot = p.IecmPlot
        self.density = p.stream_density

        self.saveFolder = saveFolder
        self.saveFile = saveFile
        self.ani_repeat = ani_repeat

        self.ax.axis('equal')

        xmin = cells.xmin*p.um
        xmax = cells.xmax*p.um
        ymin = cells.ymin*p.um
        ymax = cells.ymax*p.um

        self.ax.axis([xmin,xmax,ymin,ymax])

        if self.save is True:
            _setup_file_saving(self,p)

        self.bkgBool = False

        cell_edges_flat = cells.um*cells.mem_edges_flat

        self.coll = LineCollection(
            cell_edges_flat,
            array=sim.rho_pump_time[0],
            cmap=self.clrmap,
            linewidths=4.0,
        )

        self.ax.add_collection(self.coll)
        self.ax.axis('equal')
        self.ax.axis([xmin,xmax,ymin,ymax])

        if self.current_overlay is True:
            self.streams, self.ax, self.tit_extra = I_overlay_setup(
                sim, self.ax, cells, p)
        else:
            self.tit_extra = ' '

        # Set range of the colormap.
        if clrAutoscale is True:
            # first flatten the data (needed in case cells were cut)
            all_z = []
            for zarray in sim.rho_pump_time:
                for val in zarray:
                    all_z.append(val)

            self.cmin = np.min(all_z)
            self.cmax = np.max(all_z)

        else:
            self.cmin = clrMin
            self.cmax = clrMax

        self.coll.set_clim(self.cmin,self.cmax)
        self.cb = self.fig.colorbar(self.coll)   # define colorbar for figure
        self.cb.set_label(self.cbtit)

        self.tit = 'Pump Density Factor'

        if number_cells is True:
            for i,cll in enumerate(cells.cell_centres):
                self.ax.text(p.um*cll[0],p.um*cll[1],i,va='center',ha='center')

        self.ax.set_xlabel('Spatial x [um]')
        self.ax.set_ylabel('Spatial y [um]')
        self.fig.suptitle(self.tit,fontsize=14, fontweight='bold')
        self.ax.set_title(self.tit_extra)

        self.frames = len(sim.rho_pump_time)
        ani = animation.FuncAnimation(self.fig, self.aniFunc,
            frames=self.frames, interval=100, repeat=self.ani_repeat)

        _handle_plot(p)


    def aniFunc(self,i):

        zz = self.sim.rho_pump_time[i]
        self.coll.set_array(zz)

        if self.current_overlay is True:
            self.streams, self.ax = I_overlay_update(
                i, self.sim, self.streams, self.ax, self.cells, self.p)

        titani = self.tit_extra + ' (sim time' + ' ' + str(round(self.time[i],3)) + ' ' + ' s)'
        self.ax.set_title(titani)

        if self.save is True:
            self.fig.canvas.draw()
            savename = self.savedAni + str(i) + '.png'
            plt.savefig(savename,format='png')


#FIXME: Excise the unused "current_overlay" attribute and parameter. Mush-room!
class AnimateDyeData(object):
# class AnimateDyeData(Anim):
    '''
    Animate morphogen concentration data in cell and environment as a function
    of time.
    '''

    def __init__(
        self,
        sim,
        cells,
        p,
        save=False,
        ani_repeat=False,
        current_overlay=False,
        clrAutoscale=True,
        clrMin=None,
        clrMax=None,
        clrmap=mpl.get_colormap('rainbow'),
        number_cells=False,
        saveFolder='animation',
        saveFile='sim_',
    ):

        self.zdata_t = np.multiply(np.asarray(sim.cDye_time[:]),1e3)
        self.zenv_t = np.multiply(np.asarray(sim.cDye_env_time[:]),1e3)

        self.colormap = clrmap
        self.time = sim.time
        self.save = save

        self.cells = cells
        self.p = p

        self.fig = plt.figure()       # define figure
        self.ax = plt.subplot(111)    # define axes

        self.sim = sim
        self.current_overlay = current_overlay
        self.clrmap = clrmap

        self.sim_ECM = p.sim_ECM
        self.IecmPlot = p.IecmPlot
        self.density = p.stream_density

        self.saveFolder = saveFolder
        self.saveFile = saveFile
        self.ani_repeat = ani_repeat

        self.ax.axis('equal')

        xmin = cells.xmin*p.um
        xmax = cells.xmax*p.um
        ymin = cells.ymin*p.um
        ymax = cells.ymax*p.um

        self.ax.axis([xmin,xmax,ymin,ymax])

        if self.save is True:
            _setup_file_saving(self,p)

        self.bkgPlot = self.ax.imshow(
            self.zenv_t[0].reshape(cells.X.shape),
            origin='lower',
            extent=[xmin,xmax,ymin,ymax],
            cmap=clrmap,
        )

        # define a polygon collection based on individual cell polygons
        self.points = np.multiply(cells.cell_verts, p.um)
        self.collection =  PolyCollection(
            self.points, cmap=self.colormap, edgecolors='none')
        self.collection.set_array(self.zdata_t[0])
        self.ax.add_collection(self.collection)

        # set range of the colormap

        if clrAutoscale is True:
            # first flatten the data (needed in case cells were cut)
            all_z = []
            for zarray in self.zdata_t:
                for val in zarray:
                    all_z.append(val)

            cmina = np.min(all_z)
            cmaxa = np.max(all_z)

            cminb = np.min(self.zenv_t)
            cmaxb = np.max(self.zenv_t)

            #FIXME: Consider using Python's built-in min() and max() functions
            #instead. Penultimate zenith of the zodiac arise!
            if cmaxa > cmaxb:
                self.cmax = cmaxa
            else:
                self.cmax = cmaxb

            if cmina < cminb:
                self.cmin = cmina
            else:
                self.cmin = cminb

        else:
            self.cmin = clrMin
            self.cmax = clrMax

        self.collection.set_clim(self.cmin,self.cmax)
        self.bkgPlot.set_clim(self.cmin,self.cmax)

        self.cb = self.fig.colorbar(self.collection)   # define colorbar for figure
        self.cb.set_label('Morphogen concentration [umol/L]')

        self.tit = 'Morphogen concentration in cell and environment'

        if number_cells is True:
            for i,cll in enumerate(cells.cell_centres):
                self.ax.text(p.um*cll[0],p.um*cll[1],i,va='center',ha='center')

        self.ax.set_xlabel('Spatial x [um]')
        self.ax.set_ylabel('Spatial y [um]')
        self.fig.suptitle(self.tit,fontsize=14, fontweight='bold')

        self.frames = len(self.zdata_t)
        ani = animation.FuncAnimation(self.fig, self.aniFunc,
            frames=self.frames, interval=100, repeat=self.ani_repeat)

        _handle_plot(p)


    def aniFunc(self,i):

        zz = self.zdata_t[i]
        zenv = self.zenv_t[i]

        self.collection.set_array(zz)
        self.bkgPlot.set_data(zenv.reshape(self.cells.X.shape))

        titani = 'sim time' + ' ' + str(round(self.time[i],3)) + ' ' + ' s'
        self.ax.set_title(titani)

        if self.save is True:
            self.fig.canvas.draw()
            savename = self.savedAni + str(i) + '.png'
            plt.savefig(savename,format='png')


class PlotWhileSolving(object):

    def __init__(
        self,
        cells,
        sim,
        p,
        number_cells=False,
        clrAutoscale=True,
        clrMin=None,
        clrMax=None,
    ):

        vdata = np.multiply(sim.vm,1000)   # data array for cell coloring

        self.colormap = p.default_cm
        self.fig = plt.figure()       # define figure
        self.ax = plt.subplot(111)    # define axes

        self.tit = 'Vmem check while solving'

        self.clrAutoscale = clrAutoscale

        self.cells = cells
        self.p = p

        self.number_cells = number_cells
        self.clrMin = clrMin
        self.clrMax = clrMax

        self.ax.axis('equal')

        xmin = cells.xmin*p.um
        xmax = cells.xmax*p.um
        ymin = cells.ymin*p.um
        ymax = cells.ymax*p.um

        self.ax.axis([xmin,xmax,ymin,ymax])

        if clrAutoscale is True:
            self.cmean = np.mean(vdata)
            self.cmin = round(np.min(vdata),1)
            self.cmax = round(np.max(vdata),1)
            clrCheck = self.cmax - self.cmin

            if clrCheck == 0:
                self.cmin = self.cmin - 0.1
                self.cmax = self.cmax + 0.1

        else:
            self.cmin = clrMin
            self.cmax = clrMax

        if p.sim_ECM is False:
            if p.showCells is True:
                # Add a collection of cell polygons, with animated voltage data
                self.coll2, self.ax = cell_mosaic(
                    vdata,self.ax,cells,p,p.default_cm)
            else:
                self.coll2,self.ax = cell_mesh(
                    vdata,self.ax,cells,p,p.default_cm)

        elif p.sim_ECM is True:
            dat_grid = sim.vm_Matrix[0]*1000

            if p.plotMask is True:
                dat_grid = ma.masked_array(
                    sim.vm_Matrix[0]*1000, np.logical_not(cells.maskM))

            self.coll2 = plt.imshow(
                dat_grid,
                origin='lower',
                extent=[xmin,xmax,ymin,ymax],
                cmap=self.colormap,
            )

            if p.showCells is True:
                # cell_edges_flat, _ , _= tb.flatten(cells.mem_edges)
                cell_edges_flat = cells.um*cells.mem_edges_flat
                coll = LineCollection(cell_edges_flat,colors='k')
                coll.set_alpha(0.5)
                self.ax.add_collection(coll)

         # set range of the colormap
        self.coll2.set_clim(self.cmin,self.cmax)
        self.cb = self.fig.colorbar(self.coll2)   # define colorbar for figure

        if number_cells is True and p.showCells is True:
            for i,cll in enumerate(cells.cell_centres):
                self.ax.text(p.um*cll[0],p.um*cll[1],i,va='center',ha='center')

        self.cb.set_label('Voltage [mV]')
        self.ax.set_xlabel('Spatial x [um]')
        self.ax.set_ylabel('Spatial y [um]')
        self.ax.set_title(self.tit)

        if p.save_solving_plot is True:
            if p.run_sim is True:
                # Make the BETSE-specific cache directory if not found.
                images_path = os.path.join(p.sim_results, 'plotWhileSolving')
            else:
                images_path = os.path.join(p.init_results, 'plotWhileSolving')

            betse_cache_dir = os.path.expanduser(images_path)
            os.makedirs(betse_cache_dir, exist_ok=True)
            self.savedAni = os.path.join(betse_cache_dir, 'vm_')

            self.i = 0   # an index used for saving plot filename

        # Keep the plt.show(block=False) statement as this animation is
        # different and is closed by software.
        plt.show(block=False)


    def updatePlot(self,sim,p):

        if p.sim_ECM is False:
            if self.p.showCells is True:
                zz_grid = sim.vm_time[-1]*1000
            else:
                zz_grid = np.zeros(len(self.cells.voronoi_centres))
                zz_grid[self.cells.cell_to_grid] = sim.vm_time[-1]*1000

            self.coll2.set_array(zz_grid)

        else:
            zambie = 'nulled'

            if zambie == 'tri':
                self.coll2.set_array(sim.vm*1000)
            else:
                if p.plotMask is False:
                    zv = sim.vm_Matrix[-1]*1000
                else:
                    zv = ma.masked_array(sim.vm_Matrix[-1]*1000, np.logical_not(self.cells.maskM))

                self.coll2.set_data(zv)

        if self.clrAutoscale is True:
            cmin = 1000*np.min(sim.vm_time[-1])
            cmax = 1000*np.max(sim.vm_time[-1])
            self.coll2.set_clim(cmin,cmax)

        time = sim.time[-1]

        titani = self.tit + ' ' + '(simulation time' + ' ' + str(round(time,3)) + ' ' + 's)'
        self.ax.set_title(titani)

        self.fig.canvas.draw()

        if p.save_solving_plot is True:
            self.i = self.i + 1
            savename = self.savedAni + str(self.i) + '.png'
            plt.savefig(savename,dpi=96,format='png')


    def resetData(self,cells,sim,p):

        vdata = np.multiply(sim.vm,1000)   # data array for cell coloring

        self.cells = cells
        self.p = p

        self.fig.clf()
        self.ax = plt.subplot(111)

        xmin = p.um*cells.xmin
        xmax = p.um*cells.xmax
        ymin = p.um*cells.ymin
        ymax = p.um*cells.ymax

        self.ax.axis([xmin,xmax,ymin,ymax])

        if self.clrAutoscale is True:
            self.cmin = np.min(vdata)
            self.cmax = np.max(vdata)

        elif self.clrAutoscale is False:
            self.cmin = self.clrMin
            self.cmax = self.clrMax

        if p.sim_ECM is False:
            if p.showCells is True:
                # Add a collection of cell polygons, with animated voltage data
                self.coll2, self.ax = cell_mosaic(
                    vdata,self.ax,cells,p,p.default_cm)
            else:
                self.coll2,self.ax = cell_mesh(
                    vdata,self.ax,cells,p,p.default_cm)

        elif p.sim_ECM is True:
            dat_grid = sim.vm_Matrix[0]*1000

            if p.plotMask is True:
                dat_grid = ma.masked_array(sim.vm_Matrix[0]*1000, np.logical_not(cells.maskM))

            self.coll2 = plt.imshow(
                dat_grid,
                origin='lower',
                extent=[xmin,xmax,ymin,ymax],
                cmap=self.colormap,
            )

            if p.showCells is True:
                # cell_edges_flat, _ , _= tb.flatten(cells.mem_edges)
                cell_edges_flat = cells.um*cells.mem_edges_flat
                coll = LineCollection(cell_edges_flat,colors='k')
                coll.set_alpha(0.5)
                self.ax.add_collection(coll)

            # If the "apply external voltage" event occurred and is to be
            # plotted, plot this event.
            if p.scheduled_options['extV'] is not None and p.extVPlot is True:
                boundv = sim.v_env*1e3
                self.vext_plot = self.ax.scatter(
                    p.um*cells.env_points[:,0],
                    p.um*cells.env_points[:,1],
                    cmap=self.colormap, c=boundv, zorder=10)
                self.vext_plot.set_clim(self.cmin, self.cmax)

        # set range of the colormap
        self.coll2.set_clim(self.cmin,self.cmax)
        self.cb = self.fig.colorbar(self.coll2)   # define colorbar for figure

        if self.number_cells is True and p.showCells is True:
            for i,cll in enumerate(cells.cell_centres):
                self.ax.text(p.um*cll[0],p.um*cll[1],i,va='center',ha='center')

        # self.cb.set_label('Voltage [mV]')
        self.cb.set_label('Voltage [mV]')
        self.ax.set_xlabel('Spatial x [um]')
        self.ax.set_ylabel('Spatial y [um]')
        self.ax.set_title(self.tit)
