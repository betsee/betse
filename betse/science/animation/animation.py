#!/usr/bin/env python3
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Matplotlib-based animation classes.
'''

# ....................{ IMPORTS                            }....................
import os
import numpy as np
from betse.exceptions import BetseExceptionParameters
from betse.lib.matplotlib import mpl
from betse.science.animation.abc import Animation
# from betse.util.type import types
from matplotlib import pyplot as plt
from matplotlib import animation
from matplotlib.collections import LineCollection, PolyCollection
from numpy import ma as ma
from scipy import interpolate

#FIXME: Shift functions called only by this module here -- possibly as private
#methods of the "Animation" superclass.
from betse.science.visualize import (
    _setup_file_saving, _handle_plot, env_mesh, cell_mosaic, cell_mesh,
    I_overlay_setup, I_overlay_update, env_quiver, cell_quiver, cell_stream
)

# ....................{ IMPORTS                            }....................
#FIXME: Let's document a few of these initialization parameters. Hot dogs and
#warm afternoons in the lazy summertime!
#FIXME: Rename the aptly named "tit" attribute to "_title_plot".
#FIXME: Rename the "cbtit" attribute to "_title_colorbar".
#FIXME: Privatize all attributes, both here and in our base class, *AFTER*
#generalizing this refactoring to every class below.

class AnimateCellData(Animation):
    '''
    Animate color data on a plot of cells.
    '''

    def __init__(
        self,
        zdata_t,
        tit: str,
        cbtit: str,
        ignore_simECM: bool,
        current_overlay: bool = False,
        *args, **kwargs
    ) -> None:
        '''
        Initialize this animation.

        Parameters
        ----------
        current_overlay : bool
            `True` if overlaying electric currents or concentration flux
            streamlines on 2D plotted data.

        See the superclass `__init__()` method for all remaining parameters.
        '''

        # Pass all parameters *NOT* listed above to our superclass.
        super().__init__(*args, **kwargs)

        self.zdata_t = zdata_t
        self.ignore_simECm = ignore_simECM
        self.cbtit = cbtit
        self.current_overlay = current_overlay

        #FIXME: Unfortunately, use of pyplot and hence the pylab figure manager
        #is a bit problematic for long-lived applications like BETSE. Why?
        #Dumb memory leaks. Every time plt.figure() is called, a new figure is
        #added to the pylab figure manager. That's a problem, as it means that
        #figures are never implicitly released from memory -- even if the GUI
        #window displaying that figure has long been closed. The only way to
        #release that sort of figure from memory is to call its figure.close()
        #method. Unfortunately, we can't do that either! Why? Because figures
        #are displayed in a non-blocking manner, which means that we don't
        #actually know when that method should be called.
        #
        #The solution, of course, is to stop using the pylab figure manager
        #altogether and to instead instead directly instantiate figures and
        #canvases via Matplotlib's object-oriented API. See also:
        #    https://stackoverflow.com/questions/16334588/create-a-figure-that-is-reference-counted/16337909#16337909
        #Abundant stags in the furry forest!
        self.fig = plt.figure()       # define figure
        self.ax = plt.subplot(111)    # define axes

        self.sim_ECM = self.p.sim_ECM
        self.IecmPlot = self.p.IecmPlot
        self.density = self.p.stream_density

        self.ax.axis('equal')

        xmin = self.cells.xmin * self.p.um
        xmax = self.cells.xmax * self.p.um
        ymin = self.cells.ymin * self.p.um
        ymax = self.cells.ymax * self.p.um

        self.ax.axis([xmin, xmax, ymin, ymax])
        data_points = self.zdata_t[0]

        if self.p.sim_ECM is True and ignore_simECM is False:
            self.collection, self.ax = env_mesh(
                data_points, self.ax, self.cells, self.p, self.colormap)
        elif self.p.showCells is True:
            self.collection, self.ax = cell_mosaic(
                data_points, self.ax, self.cells, self.p, self.colormap)
        else:
            self.collection, self.ax = cell_mesh(
                data_points, self.ax, self.cells, self.p, self.colormap)

        if self.current_overlay is True:
            self.streams, self.ax, self.tit_extra = I_overlay_setup(
                self.sim, self.ax, self.cells, self.p)
        else:
            self.tit_extra = ' '

        # set range of the colormap
        if self.clrAutoscale is True:
            # first flatten the data (needed in case cells were cut)
            all_z = []
            for zarray in zdata_t:
                for val in zarray:
                    all_z.append(val)

            self.cmean = np.mean(all_z)
            self.cmin = np.min(all_z)
            self.cmax = np.max(all_z)
        else:
            self.cmin = self.clrMin
            self.cmax = self.clrMax

        # Define the figure colorbar.
        self.collection.set_clim(self.cmin,self.cmax)
        self.cb = self.fig.colorbar(self.collection)
        self.cb.set_label(self.cbtit)
        self.tit = tit

        if self.p.enumerate_cells is True:
            for i,cll in enumerate(self.cells.cell_centres):
                self.ax.text(
                    self.p.um*cll[0],
                    self.p.um*cll[1], i, va='center', ha='center')

        self.ax.set_xlabel('Spatial x [um]')
        self.ax.set_ylabel('Spatial y [um]')
        self.fig.suptitle(self.tit, fontsize=14, fontweight='bold')
        self.ax.set_title(self.tit_extra)

        # self.frames = len(self.zdata_t)
        # self.frames = len(self.sim.time)

        #FIXME: For efficiency, we should probably be passing "blit=True," to
        #FuncAnimation() both here and everywhere below. Lemon grass and dill!
        #FIXME: There appear to be a number of approaches to saving without
        #displaying animation frames, including:
        #
        #* Backend-based. This has the advantage of not requiring a movie file
        #  to also be saved, but the disadvantage of being labelled
        #  "experimental" and hence unsupported by Matplotlib:
        #  1. Switch to a non-interactive Matplotlib backend: e.g.,
        #     matplotlib.pyplot.switch_backend('Agg')
        #  2. Replace the call to plt.savefig() with a call to something else.
        #     I have no idea what. The following example suggests
        #     "plt.close(ani._fig)":
        #     http://yt-project.org/doc/cookbook/embedded_webm_animation.html
        #* Movie-based. This has the advantage of being officially supported by
        #  Matplotlib, but the disadvantage of requiring a movie file to
        #  *ALWAYS* be saved when saving frames. (Of course, the movie file
        #  could just be automatically deleted afterwards... but still). In
        #  this case:
        #  1. Implement the related FIXME below, which we want to do anyway.
        #  2. That's it. Pure profit all the way to the poor house!
        #* File writer-based. This may be the best of all possible approaches,
        #  as this approach requires no movie file to be saved to save frames.
        #  It does, however, require that we write our own frame writer class
        #  (and ideally contribute that class back to Matplotlib). Trivial!
        #  There appear to be two sort of movie writer classes in the
        #  "matplotlib.animation" module: those that (A) write frames to files
        #  and then run an external encoder command consuming those frames and
        #  (B) directly run an external encoder command without writing frames.
        #  Using the first type of writer, it looks like we should be able to
        #  use one of the existing "FileMovieWriter" subclasses or perhaps
        #  define our own. The FileMovieWriter.grab_frame() function appears to
        #  do everything we want, suggesting we:
        #  1. Define a new "FileFrameWriter" subclass of the "FileMovieWriter"
        #     base class. See "FFMpegFileWriter" for the boiler-plate.
        #  2. Move our savefig() related logic into a reimplementation of that
        #     subclass' grab_frame() method. Maybe? Seems about right.
        #  3. Redefine FileFrameWriter._run(self) to just be a noop.
        #  4. Set in FileFrameWriter.__init__():
        #         from collections import namedtuple
        #
        #         # Notify our superclass that the _run() method succeeded by
        #         # constructing a pseudo-process class providing the desired
        #         # attribute subsequently tested by FileMovieWriter.finish().
        #         FakeProc = namedtuple('FakeProc', 'returncode')
        #         self._proc = FakeProc(0)
        #  5. Use that writer instead of the currently configured
        #     video_encoder_name (e.g., "ffmpeg") when only writing frames and
        #     not a movie.
        #  6. Everything else should be the same. In particular, we'll still
        #     need to call the ani.save() method as detailed below.
        #
        #Too bad the Matplotlib documentation itself doesn't cover this. We
        #have gotten what we have paid for. I demand a refund!

        self._animate(frame_count=len(self.sim.time))


    def _plot_next_frame(self, frame_number):

        zz = self.zdata_t[frame_number]

        if self.p.sim_ECM is True and self.ignore_simECm is False:
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

        if self.current_overlay is True:
            self.streams, self.ax = I_overlay_update(
                frame_number,
                self.sim, self.streams, self.ax, self.cells, self.p)

        #FIXME: The str.format() should typically be called in favor of string
        #concatenation. Underworld of candy corn open to me!
        titani = self.tit_extra + ' (sim time' + ' ' + str(
            round(self.sim.time[frame_number], 3)) + ' ' + ' s)'
        self.ax.set_title(titani)

        # Save this frame to disk *AFTER* completing this frame.
        self._save_frame(frame_number)


class AnimateGJData(object):
# class AnimateGJData(Animation):
    '''
    Animate the gap junction open state as a function of time.
    '''

    def __init__(
        self,
        cells,
        sim,
        p,
        tit=' ',
        save=False,
        clrAutoscale=True,
        clrMin=None,
        clrMax=None,
        saveFolder='animation',
        saveFile='sim_',
        ani_repeat=False,
        number_cells=False,
    ):

        self.sim = sim
        self.cells = cells
        self.p = p

        self.zdata_t = sim.gjopen_time  # data array for gap junction coloring

        if p.gj_flux_sensitive is True:
            max_zdata = p.max_gj_enhancement
        else:
            max_zdata = 1.0

        self.vdata_t = [1000*arr for arr in sim.vm_time]   # data array for cell coloring
        self.colormap = p.default_cm
        self.time = sim.time

        self.gjI_t_x = sim.I_gj_x_time
        self.gjI_t_x = sim.I_gj_y_time
        self.gjvects_x = cells.nn_tx
        self.gjvects_y = cells.nn_ty

        self.fig = plt.figure()       # define figure
        self.ax = plt.subplot(111)    # define axes

        self.tit = tit

        self.save = save
        self.saveFolder = saveFolder
        self.saveFile = saveFile
        self.ani_repeat = ani_repeat

        if self.save is True:
            _setup_file_saving(self,p)

        con_segs = cells.nn_edges
        connects = p.um*np.asarray(con_segs)
        self.collection = LineCollection(
            connects,
            array=self.zdata_t[0], cmap=p.gj_cm, linewidths=2.0, zorder=10)
        self.collection.set_clim(0.0,max_zdata)
        self.ax.add_collection(self.collection)

        # Next add a collection of cell polygons, with animated voltage data

        if p.sim_ECM is False:
            data_set = self.vdata_t[0]
        else:
            data_set = sim.vcell_time[0]*1000

        if p.showCells is True:
            self.coll2, self.ax = cell_mosaic(data_set,self.ax,cells,p,p.default_cm)
        else:
            self.coll2, self.ax = cell_mesh(data_set,self.ax,cells,p,p.default_cm)

        # set range of the colormap
        if clrAutoscale is True:
             # first flatten the data (needed in case cells were cut)
            all_z = []
            for zarray in self.vdata_t:
                for val in zarray:
                    all_z.append(val)

            self.cmean = np.mean(all_z)
            self.cmin = round(np.min(all_z),1)
            self.cmax = round(np.max(all_z),1)

            clrCheck = self.cmax - self.cmin

            if clrCheck == 0:
                self.cmin = self.cmin - 1
                self.cmax = self.cmax + 1

        else:
            self.cmin = clrMin
            self.cmax = clrMax

        self.coll2.set_clim(self.cmin,self.cmax)
        self.cb = self.fig.colorbar(self.coll2)   # define colorbar for figure

        if number_cells is True:
            for i,cll in enumerate(cells.cell_centres):
                self.ax.text(p.um*cll[0],p.um*cll[1],i, va='center',ha='center')

        self.cb.set_label('Voltage [mV]')
        self.ax.set_xlabel('Spatial x [um]')
        self.ax.set_ylabel('Spatial y [um]')
        self.ax.set_title(self.tit)

        self.ax.axis('equal')

        xmin = cells.xmin*p.um
        xmax = cells.xmax*p.um
        ymin = cells.ymin*p.um
        ymax = cells.ymax*p.um

        self.ax.axis([xmin,xmax,ymin,ymax])

        self.frames = len(self.zdata_t)
        ani = animation.FuncAnimation(self.fig, self.aniFunc,
            frames=self.frames, interval=100, repeat=self.ani_repeat)

        _handle_plot(p)


    def aniFunc(self,i):

        zz = self.zdata_t[i]

        if self.p.sim_ECM is False:
            zv = self.vdata_t[i]
        else:
            zv = self.sim.vcell_time[i]*1000

        self.collection.set_array(zz)

        if self.p.showCells is True:
            zz_grid = zv
        else:
            zz_grid = np.zeros(len(self.cells.voronoi_centres))
            zz_grid[self.cells.cell_to_grid] = zv

        self.coll2.set_array(zz_grid)

        titani = self.tit + ' ' + '(simulation time' + ' ' + str(round(self.time[i],3)) + ' ' + 's)'
        self.ax.set_title(titani)

        if self.save is True:
            self.fig.canvas.draw()
            savename = self.savedAni + str(i) + '.png'
            plt.savefig(savename,dpi=96,format='png')


class AnimateCurrent(object):
# class AnimateCurrent(Animation):

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
        gj_current=True,
        clrMin=None,
        clrMax=None,
        clrmap=mpl.get_colormap('rainbow'),
        number_cells=False,
        saveFolder='animation',
        saveFile='sim_',
    ):

        self.clrmap = clrmap
        self.time = time
        self.save = save

        self.sim = sim
        self.current_overlay = current_overlay

        self.sim_ECM = p.sim_ECM
        self.IecmPlot = p.IecmPlot

        self.density = p.stream_density
        self.cells = cells
        self.p = p

        self.saveFolder = saveFolder
        self.saveFile = saveFile
        self.ani_repeat = ani_repeat

        self.gj_current = gj_current

        self.fig = plt.figure()       # define figure
        self.ax = plt.subplot(111)    # define axes

        self.ax.axis('equal')

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

        if gj_current is True:
            Jmag_M = np.sqrt(sim.I_gj_x_time[0]**2 + sim.I_gj_y_time[0]**2) + 1e-30

            J_x = sim.I_gj_x_time[0]/Jmag_M
            J_y = sim.I_gj_y_time[0]/Jmag_M

            lw = (3.0*Jmag_M/Jmag_M.max()) + 0.5

            self.meshplot = plt.imshow(Jmag_M, origin='lower',extent=[xmin,xmax,ymin,ymax], cmap=clrmap)

            self.streamplot = self.ax.streamplot(cells.Xgrid*p.um,cells.Ygrid*p.um,J_x,J_y,density=p.stream_density,
                linewidth=lw,color='k',cmap=clrmap,arrowsize=1.5)

            self.tit = 'Gap junction current'

            # set range of the colormap
            if clrAutoscale is True:

                self.cmin = np.min(Jmag_M)
                self.cmax = np.max(Jmag_M)

        else:

            Jmag_M = np.sqrt(sim.I_tot_x_time[1]**2 + sim.I_tot_y_time[1]**2) + 1e-30

            J_x = sim.I_tot_x_time[1]/Jmag_M
            J_y = sim.I_tot_y_time[1]/Jmag_M

            lw = (3.0*Jmag_M/Jmag_M.max()) + 0.5

            self.meshplot = plt.imshow(Jmag_M, origin='lower',extent=[xmin,xmax,ymin,ymax], cmap=clrmap)

            self.streamplot = self.ax.streamplot(cells.Xgrid*p.um,cells.Ygrid*p.um,J_x,J_y,density=p.stream_density,
                linewidth=lw,color='k',cmap=clrmap,arrowsize=1.5)

            self.tit = 'Total current'

            # # set range of the colormap
            if clrAutoscale is True:
                self.cmin = np.min(Jmag_M)
                self.cmax = np.max(Jmag_M)

        if clrAutoscale is False:
            self.meshplot.set_clim(self.cmin,self.cmax)

        self.cb = self.fig.colorbar(self.meshplot)   # define colorbar for figure
        self.cb.set_label('Current Density [A/m2]')

        self.ax.set_xlabel('Spatial x [um]')
        self.ax.set_ylabel('Spatial y [um]')
        self.ax.set_title(self.tit)

        self.frames = len(sim.time)
        ani = animation.FuncAnimation(self.fig, self.aniFunc,
            frames=self.frames, interval=100, repeat=self.ani_repeat)

        _handle_plot(p)

    def aniFunc(self,i):

        titani = self.tit + ' (simulation time' + ' ' + str(round(self.sim.time[i],3)) + ' ' + ' s)'
        self.ax.set_title(titani)

        if self.gj_current is True:

            Jmag_M = np.sqrt(self.sim.I_gj_x_time[i]**2 + self.sim.I_gj_y_time[i]**2) + 1e-30

            J_x = self.sim.I_gj_x_time[i]/Jmag_M
            J_y = self.sim.I_gj_y_time[i]/Jmag_M

            lw = (3.0*Jmag_M/Jmag_M.max()) + 0.5

            self.meshplot.set_data(Jmag_M)

            self.streamplot.lines.remove()
            self.ax.patches = []

            self.streamplot = self.ax.streamplot(self.cells.Xgrid*self.p.um,self.cells.Ygrid*self.p.um,J_x,J_y,
                density=self.p.stream_density, linewidth=lw,color='k',cmap=self.clrmap,arrowsize=1.5)

        else:

            Jmag_M = np.sqrt(self.sim.I_tot_x_time[i+1]**2 + self.sim.I_tot_y_time[i+1]**2) + 1e-30

            J_x = self.sim.I_tot_x_time[i+1]/Jmag_M
            J_y = self.sim.I_tot_y_time[i+1]/Jmag_M

            lw = (3.0*Jmag_M/Jmag_M.max()) + 0.5

            self.meshplot.set_data(Jmag_M)

            self.streamplot.lines.remove()
            self.ax.patches = []

            self.streamplot = self.ax.streamplot(self.cells.Xgrid*self.p.um,self.cells.Ygrid*self.p.um,J_x,J_y,
                density=self.p.stream_density,linewidth=lw,color='k',cmap=self.clrmap,arrowsize=1.5)

        cmax = np.max(Jmag_M)

        self.meshplot.set_clim(0,cmax)

        if self.save is True:
            self.fig.canvas.draw()
            savename = self.savedAni + str(i) + '.png'
            plt.savefig(savename,format='png')


class AnimateField(object):
# class AnimateField(Animation):

    def __init__(
        self,
        Fx,
        Fy,
        sim,
        cells,
        p,
        ani_repeat=True,
        save=True,
        saveFolder='animation/',
        saveFile='Field_',
        plot_ecm=False,
        title="Final Field in ",
        cb_title="Force [N]",
        colorAutoscale=True,
        colorMin=None,
        colorMax=None,
    ):

        self.fig = plt.figure()
        self.ax = plt.subplot(111)
        self.p = p
        self.Fx_time = Fx
        self.Fy_time = Fy
        self.sim = sim
        self.cells = cells
        self.save = save
        self.plot_ecm = plot_ecm

        self.title_piece = title
        self.cb_label = cb_title

        self.colorAutoscale = colorAutoscale
        self.colorMin = colorMin
        self.colorMax = colorMax

        self.saveFolder = saveFolder
        self.saveFile = saveFile
        self.ani_repeat = ani_repeat

        if self.save is True:
            _setup_file_saving(self,p)

        if p.sim_ECM is True and plot_ecm is True:
            efield_mag = np.sqrt(Fx[-1]**2 + Fy[-1]**2)

            self.msh, self.ax = env_mesh(efield_mag,self.ax,cells,p,p.background_cm, ignore_showCells=True)
            self.streamE, self.ax = env_quiver(Fx[-1],Fy[-1],self.ax,cells,p)

            tit_extra = 'Extracellular'

        elif plot_ecm is False:
            efield_mag = np.sqrt(Fx[-1]**2 + Fy[-1]**2)

            self.msh, self.ax = cell_mesh(
                efield_mag,self.ax,cells,p,p.background_cm)

            self.streamE, self.ax = cell_quiver(Fx[-1],Fy[-1],self.ax,cells,p)

            tit_extra = 'Intracellular'

        self.ax.axis('equal')

        xmin = cells.xmin*p.um
        xmax = cells.xmax*p.um
        ymin = cells.ymin*p.um
        ymax = cells.ymax*p.um

        self.ax.axis([xmin,xmax,ymin,ymax])

        if colorAutoscale is False:
            self.msh.set_clim(colorMin,colorMax)

        cb = self.fig.colorbar(self.msh)

        self.tit = self.title_piece + ' ' + tit_extra + ' Spaces'
        self.ax.set_title(self.tit)
        self.ax.set_xlabel('Spatial distance [um]')
        self.ax.set_ylabel('Spatial distance [um]')
        cb.set_label(self.cb_label)

        self.frames = len(sim.time)
        ani = animation.FuncAnimation(self.fig, self.aniFunc,
            frames=self.frames, interval=100, repeat=self.ani_repeat)

        _handle_plot(p)


    def aniFunc(self,i):

        titani = self.tit + ' (simulation time' + ' ' + str(round(self.sim.time[i],3)) + ' ' + ' s)'
        self.ax.set_title(titani)

        if self.p.sim_ECM is True and self.plot_ecm is True:
            E_x = self.Fx_time[i]
            E_y = self.Fy_time[i]

            efield = np.sqrt(E_x**2 + E_y**2)

            self.msh.set_data(efield)

            if efield.max() != 0.0:
                E_x = E_x/efield.max()
                E_y = E_y/efield.max()

            self.streamE.set_UVC(E_x,E_y)

        elif self.plot_ecm is False:
            E_gj_x = self.Fx_time[i]
            E_gj_y = self.Fy_time[i]

            if len(E_gj_x) != len(self.cells.cell_i):
                E_gj_x = np.dot(self.cells.M_sum_mems,E_gj_x)/self.cells.num_mems
                E_gj_y = np.dot(self.cells.M_sum_mems,E_gj_y)/self.cells.num_mems

            efield = np.sqrt(E_gj_x**2 + E_gj_y**2)

            emag_grid = np.zeros(len(self.cells.voronoi_centres))
            emag_grid[self.cells.cell_to_grid] = efield

            self.msh.set_array(emag_grid)

            if efield.all() != 0.0:

                E_gj_x = E_gj_x/efield
                E_gj_y = E_gj_y/efield

            self.streamE.set_UVC(E_gj_x,E_gj_y)

        cmax = np.max(efield)

        if self.colorAutoscale is True:
            self.msh.set_clim(0,cmax)

        if self.save is True:
            self.fig.canvas.draw()
            savename = self.savedAni + str(i) + '.png'
            plt.savefig(savename,format='png')


class AnimateVelocity(object):
# class AnimateVelocity(Animation):

    def __init__(
        self,
        sim,
        cells,
        p,
        ani_repeat=True,
        save=True,
        saveFolder='animation/Velocity',
        saveFile='Velocity_',
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

        if p.sim_ECM is True and p.ani_Velocity_type == 'ECM':
            vfield = np.sqrt(sim.u_env_x_time[0]**2 + sim.u_env_y_time[0]**2)*1e9

            self.msh = self.ax.imshow(
                vfield,
                origin='lower',
                extent=[cells.xmin*p.um, cells.xmax*p.um, cells.ymin*p.um, cells.ymax*p.um],
                cmap=p.default_cm,
            )

            vnorm = np.max(vfield)

            self.streamV = self.ax.quiver(
                p.um*cells.xypts[:,0],
                p.um*cells.xypts[:,1],
                sim.u_env_x_time[-1].ravel()/vnorm,
                sim.u_env_y_time[-1].ravel()/vnorm,
            )

            tit_extra = 'Extracellular'

        elif p.ani_Velocity_type == 'GJ' or p.sim_ECM is True:
            ugjx = sim.u_cells_x_time[0]
            ugjy = sim.u_cells_y_time[0]

            v_gj_x = interpolate.griddata(
                (cells.cell_centres[:,0], cells.cell_centres[:,1]),
                ugjx,
                (cells.Xgrid, cells.Ygrid),
                fill_value=0, method=p.interp_type)
            v_gj_x = v_gj_x*cells.maskM

            v_gj_y = interpolate.griddata(
                (cells.cell_centres[:,0], cells.cell_centres[:,1]),
                ugjy,
                (cells.Xgrid, cells.Ygrid),
                fill_value=0, method=p.interp_type)
            v_gj_y = v_gj_y*cells.maskM

            vfield = np.sqrt(v_gj_x**2 + v_gj_y**2)*1e9

            self.msh = self.ax.imshow(
                vfield,
                origin='lower',
                extent=[cells.xmin*p.um, cells.xmax*p.um, cells.ymin*p.um, cells.ymax*p.um],
                cmap=p.default_cm,
            )

            vnorm = np.max(vfield)
            lw = (3.0*vfield/vnorm) + 0.5

            self.streamV = self.ax.streamplot(
                cells.Xgrid*p.um,
                cells.Ygrid*p.um,
                v_gj_x/vnorm,
                v_gj_y/vnorm,
                density=p.stream_density,
                linewidth=lw,
                color='k',
                cmap=p.default_cm,
                arrowsize=1.5,
            )

            tit_extra = 'Intracellular'

        self.ax.axis('equal')

        xmin = cells.xmin*p.um
        xmax = cells.xmax*p.um
        ymin = cells.ymin*p.um
        ymax = cells.ymax*p.um

        self.ax.axis([xmin,xmax,ymin,ymax])

        if p.autoscale_Velocity_ani is False:
            self.msh.set_clim(p.Velocity_ani_min_clr,p.Velocity_ani_max_clr)

        cb = self.fig.colorbar(self.msh)

        self.tit = "Velocity in " + tit_extra + ' Spaces'
        self.ax.set_title(self.tit)
        self.ax.set_xlabel('Spatial distance [um]')
        self.ax.set_ylabel('Spatial distance [um]')
        cb.set_label('Velocity [nm/s]')

        self.frames = len(sim.time)
        ani = animation.FuncAnimation(self.fig, self.aniFunc,
            frames=self.frames, interval=100, repeat=self.ani_repeat)

        _handle_plot(p)


    def aniFunc(self,i):

        titani = self.tit + ' (simulation time' + ' ' + str(round(self.sim.time[i],3)) + ' ' + ' s)'
        self.ax.set_title(titani)

        if self.p.sim_ECM is True and self.p.ani_Velocity_type == 'ECM':
            vfield = np.sqrt(self.sim.u_env_x_time[i]**2 + self.sim.u_env_y_time[i]**2)*1e9

            self.msh.set_data(vfield)

            vnorm = np.max(vfield)

            self.streamV.set_UVC(self.sim.u_env_x_time[i]/vnorm,self.sim.u_env_y_time[i]/vnorm)

        elif self.p.ani_Velocity_type == 'GJ' or self.p.sim_ECM is False:

            ugjx = self.sim.u_cells_x_time[i]
            ugjy = self.sim.u_cells_y_time[i]

            u_gj_x = interpolate.griddata(
                (self.cells.cell_centres[:,0], self.cells.cell_centres[:,1]),
                ugjx,
                (self.cells.Xgrid, self.cells.Ygrid),
                fill_value=0, method=self.p.interp_type)
            u_gj_x = u_gj_x*self.cells.maskM

            u_gj_y = interpolate.griddata(
                (self.cells.cell_centres[:,0], self.cells.cell_centres[:,1]),
                ugjy,
                (self.cells.Xgrid, self.cells.Ygrid),
                fill_value=0, method=self.p.interp_type)
            u_gj_y = u_gj_y*self.cells.maskM

            vfield = np.sqrt(u_gj_x**2 + u_gj_y**2)*1e9

            self.msh.set_data(vfield)

            vnorm = np.max(vfield)

            self.streamV.lines.remove()
            self.ax.patches = []

            lw = (3.0*vfield/vnorm) + 0.5

            self.streamV = self.ax.streamplot(
                self.cells.Xgrid*self.p.um,
                self.cells.Ygrid*self.p.um,
                u_gj_x/vnorm,
                u_gj_y/vnorm,
                density=self.p.stream_density,
                linewidth=lw,
                color='k',
                cmap=self.p.default_cm,
                arrowsize=1.5,
            )

            # self.streamV.set_UVC(u_gj_x/vnorm,u_gj_y/vnorm)

        cmax = np.max(vfield)

        if self.p.autoscale_Velocity_ani is True:
            self.msh.set_clim(0,cmax)

        if self.save is True:
            self.fig.canvas.draw()
            savename = self.savedAni + str(i) + '.png'
            plt.savefig(savename,format='png')


class AnimateDeformation(object):
# class AnimateDeformation(Animation):

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
# class AnimateEnv(Animation):

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
# class AnimateDyeData(Animation):
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
