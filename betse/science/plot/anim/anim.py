#!/usr/bin/env python3
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Matplotlib-based animation classes.
'''

#FIXME: Consider generalizing the "figure_title" and "colorbar_title"
#parameters to the "abc.Anim" base class. Quite alot of boilerplate below.

#FIXME: Unfortunately, use of pyplot and hence the pylab figure manager is a
#bit problematic for long-lived applications like BETSE. Why? Dumb memory
#leaks. Every time plt.figure() is called, a new figure is added to the pylab
#figure manager. That's a problem, as it means that figures are never
#implicitly released from memory -- even if the GUI window displaying that
#figure has long been closed. The only way to release that sort of figure from
#memory is to call its figure.close() method. Unfortunately, we can't do that
#either! Why? Because figures are displayed in a non-blocking manner, which
#means that we don't actually know when that method should be called.
#
#The solution, of course, is to stop using the pylab figure manager altogether
#and to instead instead directly instantiate figures and canvases via
#Matplotlib's object-oriented API. See:
#
#    https://stackoverflow.com/questions/16334588/create-a-figure-that-is-reference-counted/16337909#16337909
#
#Abundant stags in the furry forest!
# ....................{ IMPORTS                            }....................
import os
import numpy as np
from betse.exceptions import BetseExceptionParameters
from betse.lib.matplotlib import mpl
from betse.science.plot.anim.abc import Anim
from betse.util.type import types
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
#FIXME: Let's document a few of these initialization parameters. Hot dogs and
#warm afternoons in the lazy summertime!
#FIXME: Privatize all attributes, both here and in our base class, *AFTER*
#generalizing this refactoring to every class below.

class AnimateCellData(Anim):
    '''
    Animation of arbitrary colour data plotted on the current cell cluster.

    Attributes
    ----------
    _is_ecm_ignored : bool
        `True` if ignoring extracellular spaces _or_ `False` otherwise.
    _is_current_overlay : bool
        `True` if overlaying electric currents or concentration flux
        streamlines on 2D plotted data _or_ `False` otherwise.
    _zdata_time : list
        List of all colour data to be plotted, indexed by simulation time.
    '''

    def __init__(
        self,
        zdata_time: list,
        figure_title: str,
        colorbar_title: str,
        is_ecm_ignored: bool = True,
        is_current_overlay: bool = False,
        *args, **kwargs
    ) -> None:
        '''
        Initialize this animation.

        Parameters
        ----------
        zdata_time : list
            List of all colour data to be plotted, indexed by simulation time.
        figure_title : str
            Text displayed above the figure itself.
        colorbar_title: str
            Text displayed above the figure colorbar.
        is_ecm_ignored : bool
            `True` if ignoring extracellular spaces _or_ `False` otherwise.
            Defaults to `True`.
        is_current_overlay : bool
            `True` if overlaying electric currents or concentration flux
            streamlines on 2D plotted data _or_ `False` otherwise. Defaults to
            `False`.

        See the superclass `__init__()` method for all remaining parameters.
        '''
        # assert types.is_sequence_nonstr(zdata_time), (
        #     types.assert_not_sequence_nonstr(zdata_time))
        assert types.is_bool(is_ecm_ignored), (
            types.assert_not_bool(is_ecm_ignored))
        assert types.is_bool(is_current_overlay), (
            types.assert_not_bool(is_current_overlay))

        # Pass all parameters *NOT* listed above to our superclass.
        super().__init__(*args, **kwargs)

        # Classify parameters required by the _plot_next_frame() method.
        self._zdata_time = zdata_time
        self._is_ecm_ignored = is_ecm_ignored
        self._is_current_overlay = is_current_overlay

        data_points = self._zdata_time[0]

        if self.p.sim_ECM is True and is_ecm_ignored is False:
            self.collection, self.ax = env_mesh(
                data_points, self.ax, self.cells, self.p, self.colormap)
        elif self.p.showCells is True:
            self.collection, self.ax = cell_mosaic(
                data_points, self.ax, self.cells, self.p, self.colormap)
        else:
            self.collection, self.ax = cell_mesh(
                data_points, self.ax, self.cells, self.p, self.colormap)

        if self._is_current_overlay is True:
            self.streams, self.ax, axes_title = I_overlay_setup(
                self.sim, self.ax, self.cells, self.p)
        else:
            axes_title = None

        # Autoscale the colorbar range if desired.
        if self.clrAutoscale is True:
            #FIXME: This (probably) reduces to this list comprehension both
            #here and everywhere below:
            #all_z = [val for zarray in zdata_t for val in zarray]

            # first flatten the data (needed in case cells were cut)
            all_z = []
            for zarray in zdata_time:
                for val in zarray:
                    all_z.append(val)

            self.clrMin = np.min(all_z)
            self.clrMax = np.max(all_z)

        # Display and/or save this animation.
        self._animate(
            frame_count=len(self.sim.time),
            figure_title=figure_title,
            colorbar_title=colorbar_title,
            colorbar_mapping=self.collection,
            axes_title=axes_title,
            axes_x_label='Spatial x [um]',
            axes_y_label='Spatial y [um]',
        )


    def _plot_frame_figure(self, frame_number):
        assert types.is_int(frame_number), types.assert_not_int(frame_number)

        zz = self._zdata_time[frame_number]

        if self.p.sim_ECM is True and self._is_ecm_ignored is False:
            if self.p.plotMask is True:
                dat_grid = ma.masked_array(
                    zz, np.logical_not(self.cells.maskM))

            self.collection.set_data(dat_grid)
        else:
            if self.p.showCells is True:
                self.collection.set_array(zz)
            else:
                zz_grid = np.zeros(len(self.cells.voronoi_centres))
                zz_grid[self.cells.cell_to_grid] = zz
                self.collection.set_array(zz_grid)

        if self._is_current_overlay is True:
            self.streams, self.ax = I_overlay_update(
                frame_number,
                self.sim, self.streams, self.ax, self.cells, self.p)


class AnimateCurrent(Anim):
    '''
    Animation of current plotted on the current cell cluster.

    Attributes
    ----------
    is_gj_current_only : bool
        `True` if animating only the gap junction current _or_ `False` if
        animating all current.
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
        super().__init__(*args, **kwargs)

        self.is_gj_current_only = is_gj_current_only
        self.colormap = self.p.background_cm

        if self.is_gj_current_only is True:
            figure_title = 'Gap Junction Current'
        else:
            figure_title = 'Total Current'

        # Streamplot and meshplot the Jmag_M data series for the first frame.
        Jmag_M = self._streamplot_jmag_m(frame_number=0)

        #FIXME: We repeat this same logic in other classes. Shift into a new
        #private superclass method, please.
        #FIXME: Die, pyplot! Die!
        self.meshplot = plt.imshow(
            Jmag_M,
            origin='lower',
            extent=self._axes_bounds,
            cmap=self.colormap,
        )

        # Display and/or save this animation.
        self._animate(
            frame_count=len(self.sim.time),
            figure_title=figure_title,
            colorbar_title='Current Density [A/m2]',
            colorbar_mapping=self.meshplot,
            axes_x_label='Spatial x [um]',
            axes_y_label='Spatial y [um]',
        )


    def _plot_frame_figure(self, frame_number):
        assert types.is_int(frame_number), types.assert_not_int(frame_number)

        # Erase the prior frame's streamplot before streamplotting this frame.
        self.ax.patches = []
        self.streamplot.lines.remove()

        # Streamplot and meshplot the Jmag_M data series for this frame.
        Jmag_M = self._streamplot_jmag_m(frame_number)
        self.meshplot.set_data(Jmag_M)


    def _streamplot_jmag_m(self, frame_number) -> np.ndarray:
        '''
        Streamplot and return the magnitude of the current density (Jmag_M)
        data series for the current frame.
        '''
        assert types.is_int(frame_number), types.assert_not_int(frame_number)

        if self.is_gj_current_only is True:
            Jmag_M = np.sqrt(
                self.sim.I_gj_x_time[frame_number]**2 +
                self.sim.I_gj_y_time[frame_number]**2) + 1e-30
            J_x = self.sim.I_gj_x_time[frame_number]/Jmag_M
            J_y = self.sim.I_gj_y_time[frame_number]/Jmag_M
        else:
            Jmag_M = np.sqrt(
                self.sim.I_tot_x_time[frame_number]**2 +
                self.sim.I_tot_y_time[frame_number]**2) + 1e-30
            J_x = self.sim.I_tot_x_time[frame_number]/Jmag_M
            J_y = self.sim.I_tot_y_time[frame_number]/Jmag_M

        # Classify this streamplot, thus permitting the _plot_frame_figure()
        # method to subsequently erase this streamplot's lines.
        self.streamplot = self.ax.streamplot(
            self.cells.Xgrid * self.p.um,
            self.cells.Ygrid * self.p.um,
            J_x, J_y,
            density=self.p.stream_density,
            linewidth=(3.0*Jmag_M/Jmag_M.max()) + 0.5,
            color='k',
            cmap=self.colormap,
            arrowsize=1.5,
        )

        # Rescale the colorbar range if desired.
        if self.clrAutoscale is True:
            self.clrMin = np.min(Jmag_M)
            self.clrMax = np.max(Jmag_M)

        return Jmag_M


class AnimateField(Anim):
    '''
    Animation of the electric field over either intra- or extracellular spaces
    plotted on the current cell cluster.

    Attributes
    ----------
    _Fx_time : list
        List of all electric field strength X components indexed by
        simulation time.
    _Fy_time : list
        List of all electric field strength Y components indexed by
        simulation time.
    _is_ecm : bool
        `True` if plotting the electric field over extracellular spaces;
        `True` if plotting the electric field over intracellular spaces.
    '''

    def __init__(
        self,
        Fx_time: list,
        Fy_time: list,
        figure_title: str,
        colorbar_title: str,
        is_ecm: bool,
        *args, **kwargs
    ) -> None:
        '''
        Initialize this animation.

        Parameters
        ----------
        Fx_time : list
            List of all electric field strength X components indexed by
            simulation time.
        Fy_time : list
            List of all electric field strength Y components indexed by
            simulation time.
        is_ecm : bool
            `True` if plotting the electric field over extracellular spaces;
            `True` if plotting the electric field over intracellular spaces.

        See the superclass `__init__()` method for all remaining parameters.
        '''
        assert types.is_sequence_nonstr(Fx_time), (
            types.assert_not_sequence_nonstr(Fx_time))
        assert types.is_sequence_nonstr(Fy_time), (
            types.assert_not_sequence_nonstr(Fy_time))
        assert types.is_bool(is_ecm), types.assert_not_bool(is_ecm)

        # Pass all parameters *NOT* listed above to our superclass.
        super().__init__(*args, **kwargs)

        # Classify all remaining parameters.
        self._Fx_time = Fx_time
        self._Fy_time = Fy_time
        self._is_ecm = is_ecm

        self.colormap = self.p.background_cm
        efield_mag = np.sqrt(Fx_time[-1]**2 + Fy_time[-1]**2)

        if is_ecm is True:
            if self.p.sim_ECM is False:
                raise BetseExceptionParameters(
                    'Electric field animation "{}" plotted over '
                    'extracellular spaces, but '
                    'extracellular spaces are disabled by the '
                    'current simulation configuration.'.format(
                    self._type))

            figure_title_scope = 'Extra'
            self.msh, self.ax = env_mesh(
                efield_mag, self.ax, self.cells, self.p, self.colormap,
                ignore_showCells=True)
            self.streamE, self.ax = env_quiver(
                Fx_time[-1], Fy_time[-1], self.ax, self.cells, self.p)
        else:
            figure_title_scope = 'Intra'
            self.msh, self.ax = cell_mesh(
                efield_mag, self.ax, self.cells, self.p, self.colormap)
            self.streamE, self.ax = cell_quiver(
                Fx_time[-1], Fy_time[-1], self.ax, self.cells, self.p)

        # Autoscale the colorbar range if desired.
        if self.clrAutoscale is True:
            self.clrMin = np.min(efield_mag)
            self.clrMax = np.max(efield_mag)

        # Display and/or save this animation.
        self._animate(
            frame_count=len(self.sim.time),
            figure_title='{}cellular Spaces {}'.format(
                figure_title_scope, figure_title),
            colorbar_title=colorbar_title,
            colorbar_mapping=self.msh,
            axes_x_label='Spatial distance [um]',
            axes_y_label='Spatial distance [um]',
        )


    def _plot_frame_figure(self, frame_number):
        assert types.is_int(frame_number), types.assert_not_int(frame_number)

        if self._is_ecm is True:
            E_x = self._Fx_time[frame_number]
            E_y = self._Fy_time[frame_number]

            efield_mag = np.sqrt(E_x**2 + E_y**2)
            self.msh.set_data(efield_mag)

            if efield_mag.max() != 0.0:
                E_x = E_x/efield_mag.max()
                E_y = E_y/efield_mag.max()

            self.streamE.set_UVC(E_x, E_y)

        else:
            E_gj_x = self._Fx_time[frame_number]
            E_gj_y = self._Fy_time[frame_number]

            if len(E_gj_x) != len(self.cells.cell_i):
                E_gj_x = (
                    np.dot(self.cells.M_sum_mems, E_gj_x)/self.cells.num_mems)
                E_gj_y = (
                    np.dot(self.cells.M_sum_mems, E_gj_y)/self.cells.num_mems)

            efield_mag = np.sqrt(E_gj_x**2 + E_gj_y**2)
            emag_grid = np.zeros(len(self.cells.voronoi_centres))
            emag_grid[self.cells.cell_to_grid] = efield_mag
            self.msh.set_array(emag_grid)

            if efield_mag.all() != 0.0:
                E_gj_x = E_gj_x/efield_mag
                E_gj_y = E_gj_y/efield_mag

            self.streamE.set_UVC(E_gj_x, E_gj_y)

        # Rescale the colorbar range if desired.
        if self.clrAutoscale is True:
            self.clrMin = np.min(efield_mag)
            self.clrMax = np.max(efield_mag)


class AnimateGJData(Anim):
    '''
    Animation of the gap junction open state as a function of time.

    Attributes
    ----------
    zdata_t : list
        Data array for gap junction coloring.
    '''

    def __init__(
        self,
        figure_title: str,
        colorbar_title: str,
        *args, **kwargs
    ) -> None:
        '''
        Initialize this animation.

        Parameters
        ----------
        figure_title : str
            Text displayed above the figure itself.
        colorbar_title : str
            Text displayed above the figure colorbar.

        See the superclass `__init__()` method for all remaining parameters.
        '''

        # Pass all parameters *NOT* listed above to our superclass.
        super().__init__(*args, **kwargs)

        # Data array for gap junction coloring.
        self.zdata_t = self.sim.gjopen_time

        # Data array for cell coloring.
        self.vdata_t = [1000*arr for arr in self.sim.vm_time]

        connects = self.p.um * np.asarray(self.cells.nn_edges)
        self.collection = LineCollection(
            connects,
            array=self.zdata_t[0],
            cmap=self.p.gj_cm,
            linewidths=2.0,
            zorder=10,
        )
        self.collection.set_clim(0.0, 1.0)
        self.ax.add_collection(self.collection)

        # Add a collection of cell polygons with animated voltage data.
        if self.p.sim_ECM is False:
            data_set = self.vdata_t[0]
        else:
            data_set = self.sim.vcell_time[0]*1000

        if self.p.showCells is True:
            self.coll2, self.ax = cell_mosaic(
                data_set, self.ax, self.cells, self.p, self.colormap)
        else:
            self.coll2, self.ax = cell_mesh(
                data_set, self.ax, self.cells, self.p, self.colormap)

        # Autoscale the colorbar range if desired.
        if self.clrAutoscale is True:
            #FIXME: See above.
            # first flatten the data (needed in case cells were cut)
            all_z = []
            for zarray in self.vdata_t:
                for val in zarray:
                    all_z.append(val)

            self.clrMin = round(np.min(all_z), 1)
            self.clrMax = round(np.max(all_z), 1)

            if self.clrMax - self.clrMin == 0:
                self.clrMin = self.clrMin - 1
                self.clrMax = self.clrMax + 1

        # Display and/or save this animation. # FIXME: There's possibly a scaling issue that's climming to Vdata rather than gj
        self._animate(
            frame_count=len(self.zdata_t),
            figure_title=figure_title,
            colorbar_title=colorbar_title,
            colorbar_mapping=self.coll2,
            axes_x_label='Spatial x [um]',
            axes_y_label='Spatial y [um]',
        )


    def _plot_frame_figure(self, frame_number):
        assert types.is_int(frame_number), types.assert_not_int(frame_number)

        zz = self.zdata_t[frame_number]
        self.collection.set_array(zz)
        self.collection.set_clim(0.0, 1.0)

        if self.p.sim_ECM is False:
            zv = self.vdata_t[frame_number]
        else:
            zv = self.sim.vcell_time[frame_number]*1000

        if self.p.showCells is True:
            zz_grid = zv
        else:
            zz_grid = np.zeros(len(self.cells.voronoi_centres))
            zz_grid[self.cells.cell_to_grid] = zv

        self.coll2.set_array(zz_grid)


#FIXME: We probably want to use enum objects in place of "vtype". See also:
#    https://docs.python.org/3/library/enum.html#functional-api
#FIXME: We probably want to use a similar approach elsewhere. Grep
#above for "is_ecm".

class AnimateVelocity(Anim):
    '''
    Animation of fluid velocity over either intra- or extracellular spaces
    plotted on the current cell cluster.

    Attributes
    ----------
    _is_ecm : bool
        `True` if plotting the fluid velocity over extracellular spaces;
        `True` if plotting the fluid velocity over intracellular spaces.
    '''

    def __init__(
        self,
        vtype: str,
        figure_title: str,
        colorbar_title: str,
        *args, **kwargs
    ) -> None:
        '''
        Initialize this animation.

        Parameters
        ----------
        vtype : str
            Type of fluid velocity to be plotted. Valid values include:
            * `ECM`, plotting extracellular spaces.
            * `GJ`, plotting intracellular gap junctions.
        figure_title : str
            Text displayed above the figure itself.
        colorbar_title: str
            Text displayed above the figure colorbar.

        See the superclass `__init__()` method for all remaining parameters.
        '''
        assert types.is_str(vtype), types.assert_not_str(vtype)

        # Pass all parameters *NOT* listed above to our superclass.
        super().__init__(*args, **kwargs)

        # Classify all remaining parameters.
        self.vtype = vtype

        # Validate the passed fluid type.
        if vtype == 'ECM' and self.p.sim_ECM is False:
            raise BetseExceptionParameters(
                'Fluid velocity animation "{}" plotted over '
                'extracellular spaces, but '
                'extracellular spaces are disabled by the '
                'current simulation configuration.'.format(
                self._type))

        # Velocity field and maximum velocity field value for the first frame.
        vfield, vnorm = self._get_velocity_field(frame_number=0)

        if vtype == 'ECM':
            figure_title_scope = 'Extra'
            self.streamV = self.ax.quiver(
                self.cells.xypts[:,0] * self.p.um,
                self.cells.xypts[:,1] * self.p.um,
                self.sim.u_env_x_time[-1].ravel()/vnorm,
                self.sim.u_env_y_time[-1].ravel()/vnorm,
            )
        elif vtype == 'GJ':
            figure_title_scope = 'Intra'
        else:
            raise BetseExceptionParameters(
                'Fluid velocity animation type "{}" unrecognized.'.format(
                vtype))

        # Vector field mesh for the first frame.
        self.msh = self.ax.imshow(
            vfield,
            origin='lower',
            extent=self._axes_bounds,
            cmap=self.colormap,
        )

        # Display and/or save this animation.
        self._animate(
            frame_count=len(self.sim.time),
            figure_title='{}cellular Spaces {}'.format(
                figure_title_scope, figure_title),
            colorbar_title=colorbar_title,
            colorbar_mapping=self.msh,
            axes_x_label='Spatial distance [um]',
            axes_y_label='Spatial distance [um]',
        )


    def _plot_frame_figure(self, frame_number):
        assert types.is_int(frame_number), types.assert_not_int(frame_number)

        # If the current fluid type is intracellular gap junctions, erase the
        # prior frame's streamplot before streamplotting this frame.
        if self.vtype == 'GJ':
            self.ax.patches = []
            self.streamV.lines.remove()

        # Velocity field and maximum velocity field value for this frame.
        vfield, vnorm = self._get_velocity_field(frame_number)

        # Update the current velocity field mesh.
        self.msh.set_data(vfield)

        # If the current fluid type is extracellular spaces, streamplot this
        # frame.
        if self.vtype == 'ECM':
            self.streamV.set_UVC(
                self.sim.u_env_x_time[frame_number]/vnorm,
                self.sim.u_env_y_time[frame_number]/vnorm)


    #FIXME: Duplicate code as in the "AnimateCurrent" class abounds.
    def _get_velocity_field(self, frame_number) -> tuple:
        '''
        Get a 2-element tuple whose first element is the velocity field and
        second element is the maximum value of that field for the current
        frame.

        If the current fluid type is intracellular gap junctions, the
        streamplot of the velocity field will also be redefined.

        Returns
        ----------
        `(velocity_field, velocity_field_max_value)`
            2-element tuple as described above.
        '''
        assert types.is_int(frame_number), types.assert_not_int(frame_number)

        if self.vtype == 'ECM':
            vfield = np.sqrt(
                self.sim.u_env_x_time[frame_number]**2 +
                self.sim.u_env_y_time[frame_number]**2) * 1e9
            vnorm = np.max(vfield)

        elif self.vtype == 'GJ':
            #FIXME: Can at least "cell_grid" be cached between calls to this
            #method? We suspect "Maybe."
            cell_centres = (
                self.cells.cell_centres[:,0], self.cells.cell_centres[:,1])
            cell_grid = (self.cells.Xgrid, self.cells.Ygrid)

            u_gj_x = interpolate.griddata(
                cell_centres,
                self.sim.u_cells_x_time[frame_number],
                cell_grid,
                fill_value=0,
                method=self.p.interp_type,
            )
            u_gj_y = interpolate.griddata(
                cell_centres,
                self.sim.u_cells_y_time[frame_number],
                cell_grid,
                fill_value=0,
                method=self.p.interp_type,
            )

            u_gj_x = u_gj_x * self.cells.maskM
            u_gj_y = u_gj_y * self.cells.maskM

            vfield = np.sqrt(u_gj_x**2 + u_gj_y**2)*1e9
            vnorm = np.max(vfield)

            self.streamV = self.ax.streamplot(
                self.cells.Xgrid * self.p.um,
                self.cells.Ygrid * self.p.um,
                u_gj_x/vnorm,
                u_gj_y/vnorm,
                density=self.p.stream_density,
                linewidth=(3.0*vfield/vnorm) + 0.5,
                color='k',
                cmap=self.colormap,
                arrowsize=1.5,
            )
            # self.streamV.set_UVC(u_gj_x/vnorm,u_gj_y/vnorm)

        # Rescale the colorbar range if desired.
        if self.clrAutoscale is True:
            self.clrMin = np.min(vfield)
            self.clrMax = vnorm

        return (vfield, vnorm)


class AnimateDeformation(object):
# class AnimateDeformation(Anim):

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

        self.frames = len(sim.time)
        ani = animation.FuncAnimation(self.fig, self.aniFunc,
            frames=self.frames, interval=100, repeat=self.ani_repeat)

        _handle_plot(p)


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

        # keep the plt.show(block=False) statement as this animation is different and is closed by software
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

            # if p.showCells is True:
            #     # Add a collection of cell polygons, with animated voltage data
            #     points = np.multiply(cells.cell_verts, p.um)
            #     self.coll2 =  PolyCollection(points, array=vdata, edgecolors='none', cmap=self.colormap)
            #     self.coll2.set_alpha(1.0)
            #     self.ax.add_collection(self.coll2)
            #
            # else:
            #
            #     dat_grid = interpolate.griddata((cells.cell_centres[:, 0],cells.cell_centres[:, 1]),vdata,
            #         (cells.Xgrid,cells.Ygrid),fill_value=0,method=p.interp_type)
            #
            #     # dat_grid = np.multiply(dat_grid,cells.maskM)
            #     #
            #     if p.plotMask is True:
            #         dat_grid = ma.masked_array(dat_grid, np.logical_not(cells.maskM))
            #     #
            #     self.coll2 = plt.pcolormesh(p.um*cells.Xgrid, p.um*cells.Ygrid,dat_grid,shading='gouraud',
            #         cmap=self.colormap)

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

            # dat_grid = sim.vm_Matrix[0]*1000
            #
            # if p.plotMask is True:
            #     dat_grid = ma.masked_array(sim.vm_Matrix[0]*1000, np.logical_not(cells.maskM))
            #
            # self.coll2 = plt.imshow(dat_grid,origin='lower',extent=[xmin,xmax,ymin,ymax],cmap=self.colormap)

            # If the "apply external voltage" event occurred and is to be
            # plotted, plot this event.
            if p.scheduled_options['extV'] is not None and p.extVPlot is True:
                boundv = sim.v_env*1e3
                self.vext_plot = self.ax.scatter(
                    p.um*cells.env_points[:,0],
                    p.um*cells.env_points[:,1],
                    cmap=self.colormap, c=boundv, zorder=10)
                self.vext_plot.set_clim(self.cmin, self.cmax)

            # if p.showCells is True:
            #     # cell_edges_flat, _ , _= tb.flatten(cells.mem_edges)
            #     cell_edges_flat = cells.um*cells.mem_edges_flat
            #     coll = LineCollection(cell_edges_flat,colors='k')
            #     coll.set_alpha(0.5)
            #     self.ax.add_collection(coll)

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
