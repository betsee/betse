#!/usr/bin/env python3
# Copyright 2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

# FIXME saving animations as video files directly doesn't work

import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import os, os.path
from betse.util.lib import matplotlibs
from betse.util.type import types
from matplotlib import animation
from matplotlib.collections import LineCollection, PolyCollection
from scipy import interpolate

#FIXME: Let's document a few of these initialization parameters. Hot dogs and
#warm afternoons in the lazy summertime!
class AnimateCellData(object):
    '''
    Animate color data on a plot of cells.
    '''

    def __init__(
        self,
        sim: 'Simulation',
        cells: 'Cells',
        zdata_t,
        time,
        p: 'Parameters',
        tit: str = ' ',
        cbtit: str = ' ',
        save: bool = False,
        ani_repeat: bool = False,
        current_overlay: bool = False,
        clrAutoscale: bool = True,
        clrMin = None,
        clrMax = None,
        clrmap: 'Colormap' = matplotlibs.get_colormap('rainbow'),
        number_cells: bool = False,
        saveFolder = 'animation',
        saveFile = 'sim_',
        ignore_simECM: bool = False,
    ) -> None:

        self.zdata_t = zdata_t
        self.colormap = clrmap
        self.time = time
        self.save = save
        self.ignore_simECm = ignore_simECM
        self.cbtit = cbtit
        self.cells = cells
        self.p = p
        self.sim = sim
        self.current_overlay = current_overlay
        self.clrmap = clrmap

        self.saveFolder = saveFolder
        self.saveFile = saveFile
        self.ani_repeat = ani_repeat

        self.fig = plt.figure()       # define figure
        self.ax = plt.subplot(111)    # define axes

        self.sim_ECM = p.sim_ECM
        self.IecmPlot = p.IecmPlot
        self.density = p.stream_density

        self.ax.axis('equal')

        xmin = cells.xmin*p.um
        xmax = cells.xmax*p.um
        ymin = cells.ymin*p.um
        ymax = cells.ymax*p.um

        self.ax.axis([xmin,xmax,ymin,ymax])

        if self.save is True:

            set_up_filesave(self,p)

        data_points = self.zdata_t[0]

        if p.sim_ECM is True and ignore_simECM is False:

            self.collection, self.ax = env_mesh(data_points,self.ax,cells,p,p.default_cm)

        elif p.sim_ECM is False or ignore_simECM is True:

            if p.showCells is True:
                self.collection, self.ax = cell_mosaic(data_points,self.ax,cells,p,p.default_cm)

            else:
                self.collection,self.ax = cell_mesh(data_points,self.ax,cells,p,p.default_cm)

        if self.current_overlay is True:

            self.streams, self.ax, self.tit_extra = I_overlay_setup(sim,self.ax,cells,p)

        else:
            self.tit_extra = ' '

        # set range of the colormap
        if clrAutoscale is True:
            # first flatten the data (needed in case cells were cut)
            all_z = []
            for zarray in zdata_t:
                for val in zarray:
                    all_z.append(val)

            self.cmean = np.mean(all_z)
            self.cmin = np.min(all_z)
            self.cmax = np.max(all_z)

        else:
            self.cmin = clrMin
            self.cmax = clrMax

        self.collection.set_clim(self.cmin,self.cmax)

        self.cb = self.fig.colorbar(self.collection)   # define colorbar for figure
        self.cb.set_label(self.cbtit)

        self.tit = tit

        if number_cells is True:
            for i,cll in enumerate(cells.cell_centres):
                self.ax.text(p.um*cll[0],p.um*cll[1],i,va='center',ha='center')

        self.ax.set_xlabel('Spatial x [um]')
        self.ax.set_ylabel('Spatial y [um]')
        self.fig.suptitle(self.tit, fontsize=14, fontweight='bold')
        self.ax.set_title(self.tit_extra)

        # self.frames = len(self.zdata_t)
        self.frames = len(sim.time)

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

        # The animation function needs to be assigned to a local variable, otherwise the
        # animation will not take place:
        ani = animation.FuncAnimation(self.fig, self.aniFunc,
            frames=self.frames, interval=100, repeat=self.ani_repeat)

        #FIXME: It appears that movies can be saved at this exact point via the
        #following lines:
        #
        #    video_filename = 'my_filename.mp4'
        #    video_encoder_name = 'ffmpeg'
        #    ani.save(filename=video_filename, writer=video_encoder_name)
        #
        #Both the "video_filename" and "video_encoder_name" variables used
        #above should be initialized from YAML configuration items. The latter
        #should probably be configured as a list of preferred encoders: e.g.,
        #
        #    # List of Matplotlib-specific names of preferred video encoders.
        #    # The first encoder available on the current $PATH will be used.
        #    video encoders: ['ffmpeg', 'avconv', 'mencoder']
        #
        #Given that, we would then need to iteratively search this list until
        #finding the first encoder available on the current $PATH. Doesn't look
        #too hard, actually. Grandest joy in the cornucopia of easy winter!
        #FIXME: Also note that movie encoders are configurable as follows,
        #which is probably a better idea than just passing a string above:
        #
        #    Writer = animation.writers[video_encoder_name]
        #    writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
        #    ani.save(filename=video_filename, writer=writer)
        #
        #Oh, wait. No, that's overkill. All of the above parameters are also
        #accepted by the ani.save() function itself, so string name it is!

        show_plot(p)

    def aniFunc(self,i):

        zz = self.zdata_t[i]

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

            self.streams, self.ax = I_overlay_update(i,self.sim,self.streams,self.ax,self.cells,self.p)

        titani = self.tit_extra + ' (sim time' + ' ' + str(
            round(self.time[i],3)) + ' ' + ' s)'
        self.ax.set_title(titani)

        if self.save is True:
            self.fig.canvas.draw()
            savename = self.savedAni + str(i) + '.png'

            # #FIXME: Remove this debug statement later.
            # print('Saving animated frame: {}'.format(savename))
            # plt.savefig(savename, format='png')

class AnimateGJData(object):
    """
    Animate the gap junction open state as a function of time.
    """

    def __init__(self,cells,sim,p,tit=' ', save=False,saveFolder = 'animation',
        clrAutoscale = True, clrMin = None, clrMax = None,
        saveFile = 'sim_',ani_repeat=False,number_cells=False):

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

        self.cells = cells
        self.p = p
        self.sim = sim

        self.fig = plt.figure()       # define figure
        self.ax = plt.subplot(111)    # define axes

        self.tit = tit

        self.save = save
        self.saveFolder = saveFolder
        self.saveFile = saveFile
        self.ani_repeat = ani_repeat

        if self.save is True:

            set_up_filesave(self,p)

        con_segs = cells.nn_edges
        connects = p.um*np.asarray(con_segs)
        self.collection = LineCollection(connects, array=self.zdata_t[0], cmap= p.gj_cm, linewidths=1.0, zorder=10)
        self.collection.set_clim(0.0,max_zdata)
        self.ax.add_collection(self.collection)

        # Next add a collection of cell polygons, with animated voltage data

        if p.sim_ECM is False:
            data_set = self.vdata_t[0]
        else:
            data_set = sim.vcell_time[0]

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

        elif clrAutoscale is False:
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

        show_plot(p)


    def aniFunc(self,i):

        zz = self.zdata_t[i]

        if self.p.sim_ECM is False:
            zv = self.vdata_t[i]

        else:
            zv = self.sim.vcell_time[i]*1000

        self.collection.set_array(zz)

        if self.p.showCells is True:
            zz_grid = zz

        else:

            zz_grid = np.zeros(len(self.cells.voronoi_centres))
            zz_grid[self.cells.cell_to_grid] = zz

        self.coll2.set_array(zz_grid)

        titani = self.tit + ' ' + '(simulation time' + ' ' + str(round(self.time[i],3)) + ' ' + 's)'
        self.ax.set_title(titani)

        if self.save is True:
            self.fig.canvas.draw()
            savename = self.savedAni + str(i) + '.png'
            plt.savefig(savename,dpi=96,format='png')

class PlotWhileSolving(object):

    def __init__(self,cells,sim,p,number_cells=False,clrAutoscale = True, clrMin = None, clrMax = None):


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

        elif clrAutoscale is False:

            self.cmin = clrMin
            self.cmax = clrMax

        if p.sim_ECM is False:

            if p.showCells is True:
                # Add a collection of cell polygons, with animated voltage data
                self.coll2, self.ax = cell_mosaic(vdata,self.ax,cells,p,p.default_cm)

            else:

                self.coll2,self.ax = cell_mesh(vdata,self.ax,cells,p,p.default_cm)

        elif p.sim_ECM is True:

            dat_grid = sim.vm_Matrix[0]*1000

            if p.plotMask is True:
                dat_grid = ma.masked_array(sim.vm_Matrix[0]*1000, np.logical_not(cells.maskM))

            self.coll2 = plt.imshow(dat_grid,origin='lower',extent=[xmin,xmax,ymin,ymax],cmap=self.colormap)

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

        if p.save_solving_plot is True:  # FIXME do this uniquely

            if p.run_sim is True:
                # Make the BETSE-specific cache directory if not found.
                images_path = p.sim_results + '/plotWhileSolving'

            else:
                images_path = p.init_results + '/plotWhileSolving'

            betse_cache_dir = os.path.expanduser(images_path)
            os.makedirs(betse_cache_dir, exist_ok=True)
            self.savedAni = os.path.join(betse_cache_dir, 'vm_')

            self.i = 0   # an index used for saving plot filename

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
                points = np.multiply(cells.cell_verts, p.um)
                self.coll2 =  PolyCollection(points, array=vdata, edgecolors='none', cmap=self.colormap)
                self.coll2.set_alpha(1.0)
                self.ax.add_collection(self.coll2)

            else:

                dat_grid = interpolate.griddata((cells.cell_centres[:, 0],cells.cell_centres[:, 1]),vdata,
                    (cells.Xgrid,cells.Ygrid),fill_value=0,method=p.interp_type)

                # dat_grid = np.multiply(dat_grid,cells.maskM)
                #
                if p.plotMask is True:
                    dat_grid = ma.masked_array(dat_grid, np.logical_not(cells.maskM))
                #
                self.coll2 = plt.pcolormesh(p.um*cells.Xgrid, p.um*cells.Ygrid,dat_grid,shading='gouraud',
                    cmap=self.colormap)

        elif p.sim_ECM is True:

            dat_grid = sim.vm_Matrix[0]*1000

            if p.plotMask is True:
                dat_grid = ma.masked_array(sim.vm_Matrix[0]*1000, np.logical_not(cells.maskM))

            self.coll2 = plt.imshow(dat_grid,origin='lower',extent=[xmin,xmax,ymin,ymax],cmap=self.colormap)

            # If the "apply external voltage" event occurred and is to be
            # plotted, plot this event.
            if p.scheduled_options['extV'] is not None and p.extVPlot is True:
                boundv = sim.v_env*1e3
                self.vext_plot = self.ax.scatter(
                    p.um*cells.env_points[:,0],
                    p.um*cells.env_points[:,1],
                    cmap=self.colormap, c=boundv, zorder=10)
                self.vext_plot.set_clim(self.cmin, self.cmax)

            if p.showCells is True:
                # cell_edges_flat, _ , _= tb.flatten(cells.mem_edges)
                cell_edges_flat = cells.um*cells.mem_edges_flat
                coll = LineCollection(cell_edges_flat,colors='k')
                coll.set_alpha(0.5)
                self.ax.add_collection(coll)

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

class AnimateCurrent(object):

    def __init__(self,sim,cells,time,p,save=False,ani_repeat=False,current_overlay=False,clrAutoscale=True, gj_current = True,
    clrMin = None,clrMax = None,clrmap = cm.rainbow, number_cells=False,saveFolder = 'animation',saveFile = 'sim_'):

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

            set_up_filesave(self,p)

        if clrAutoscale is False:
            self.cmin = clrMin
            self.cmax = clrMax

        if gj_current is True:

            Jmag_M = np.sqrt(sim.I_gj_x_time[0]**2 + sim.I_gj_y_time[0]**2) + 1e-30

            J_x = sim.I_gj_x_time[0]/Jmag_M
            J_y = sim.I_gj_y_time[0]/Jmag_M

            lw = (3.0*Jmag_M/Jmag_M.max()) + 0.5

            self.meshplot = plt.imshow(Jmag_M, origin='lower',extent=[xmin,xmax,ymin,ymax], cmap=clrmap)

            # if p.I_overlay is True:

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

            # if p.I_overlay is True:

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

        show_plot(p)

    def aniFunc(self,i):

        titani = self.tit + ' (simulation time' + ' ' + str(round(self.sim.time[i],3)) + ' ' + ' s)'
        self.ax.set_title(titani)

        if self.gj_current is True:

            Jmag_M = np.sqrt(self.sim.I_gj_x_time[i]**2 + self.sim.I_gj_y_time[i]**2) + 1e-30

            J_x = self.sim.I_gj_x_time[i]/Jmag_M
            J_y = self.sim.I_gj_y_time[i]/Jmag_M

            lw = (3.0*Jmag_M/Jmag_M.max()) + 0.5

            self.meshplot.set_data(Jmag_M)

            # if self.p.I_overlay is True:

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

            # if self.p.I_overlay is True:

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

class AnimateEfield(object):

    def __init__(self,sim,cells,p,ani_repeat = True, save = True, saveFolder = 'animation/Efield',saveFile = 'Efield_'):

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

            set_up_filesave(self,p)

        if p.sim_ECM is True and p.ani_Efield_type == 'ECM':

            efield = np.sqrt(sim.efield_ecm_x_time[-1]**2 + sim.efield_ecm_y_time[-1]**2)

            self.msh = self.ax.imshow(efield,origin='lower', extent = [cells.xmin*p.um, cells.xmax*p.um,
                cells.ymin*p.um, cells.ymax*p.um], cmap=p.default_cm)

            if p.plot_Efield_vector is True:
                enorm = np.max(np.sqrt(sim.efield_ecm_x_time[-1]**2 + sim.efield_ecm_y_time[-1]**2))

                self.streamE = self.ax.quiver(p.um*cells.xypts[:,0], p.um*cells.xypts[:,1],
                    sim.efield_ecm_x_time[-1].ravel()/enorm,sim.efield_ecm_y_time[-1].ravel()/enorm)

            tit_extra = 'Extracellular'

        elif p.ani_Efield_type == 'GJ' or p.sim_ECM is False:

            E_gj_x = interpolate.griddata((cells.nn_mids[:,0],cells.nn_mids[:,1]),
            sim.efield_gj_x_time[-1],(cells.Xgrid,cells.Ygrid), fill_value=0,method=p.interp_type)

            E_gj_x = np.multiply(E_gj_x,cells.maskM)

            E_gj_y = interpolate.griddata((cells.nn_mids[:,0],cells.nn_mids[:,1]),
                sim.efield_gj_y_time[-1],(cells.Xgrid,cells.Ygrid), fill_value=0,method=p.interp_type)

            E_gj_y = np.multiply(E_gj_y, cells.maskM)

            efield = np.sqrt(E_gj_x**2 + E_gj_y**2)
            self.msh = self.ax.imshow(efield,origin='lower', extent = [cells.xmin*p.um, cells.xmax*p.um,
                cells.ymin*p.um, cells.ymax*p.um],cmap=p.default_cm)

            if p.ani_Efield_vector is True:

                enorm = np.max(efield)

                self.streamE = self.ax.quiver(p.um*cells.Xgrid, p.um*cells.Ygrid,
                    E_gj_x/enorm,E_gj_y/enorm,scale=10)

            tit_extra = 'Intracellular'

        self.ax.axis('equal')

        xmin = cells.xmin*p.um
        xmax = cells.xmax*p.um
        ymin = cells.ymin*p.um
        ymax = cells.ymax*p.um

        self.ax.axis([xmin,xmax,ymin,ymax])

        if p.autoscale_Efield_ani is False:
            self.msh.set_clim(p.Efield_ani_min_clr,p.Efield_ani_max_clr)

        cb = self.fig.colorbar(self.msh)

        self.tit = "Final Electric Field in " + tit_extra + ' Spaces'
        self.ax.set_title(self.tit)
        self.ax.set_xlabel('Spatial distance [um]')
        self.ax.set_ylabel('Spatial distance [um]')
        cb.set_label('Electric Field [V/m]')

        self.frames = len(sim.time)
        ani = animation.FuncAnimation(self.fig, self.aniFunc,
            frames=self.frames, interval=100, repeat=self.ani_repeat)

        show_plot(p)


    def aniFunc(self,i):

        titani = self.tit + ' (simulation time' + ' ' + str(round(self.sim.time[i],3)) + ' ' + ' s)'
        self.ax.set_title(titani)

        if self.p.sim_ECM is True and self.p.ani_Efield_type == 'ECM':

            efield = np.sqrt(self.sim.efield_ecm_x_time[i]**2 + self.sim.efield_ecm_y_time[i]**2)
            self.msh.set_data(efield)

            if self.p.ani_Efield_vector is True:

                enorm = np.max(np.sqrt(self.sim.efield_ecm_x_time[i]**2 + self.sim.efield_ecm_y_time[i]**2))
                self.streamE.set_UVC(self.sim.efield_ecm_x_time[i]/enorm,self.sim.efield_ecm_y_time[i]/enorm)

        elif self.p.ani_Efield_type == 'GJ' or self.p.sim_ECM is False:

            E_gj_x = interpolate.griddata((self.cells.nn_mids[:,0],self.cells.nn_mids[:,1]),
            self.sim.efield_gj_x_time[i],(self.cells.Xgrid,self.cells.Ygrid), fill_value=0,method=self.p.interp_type)

            E_gj_x = np.multiply(E_gj_x,self.cells.maskM)

            E_gj_y = interpolate.griddata((self.cells.nn_mids[:,0],self.cells.nn_mids[:,1]),
                self.sim.efield_gj_y_time[i],(self.cells.Xgrid,self.cells.Ygrid), fill_value=0,method=self.p.interp_type)

            E_gj_y = np.multiply(E_gj_y,self.cells.maskM)

            efield = np.sqrt(E_gj_x**2 + E_gj_y**2)

            self.msh.set_data(efield)

            if self.p.ani_Efield_vector is True:

                enorm = np.max(efield)
                self.streamE.set_UVC(E_gj_x/enorm,E_gj_y/enorm)

        cmax = np.max(efield)

        if self.p.autoscale_Efield_ani is True:
            self.msh.set_clim(0,cmax)

        if self.save is True:
            self.fig.canvas.draw()
            savename = self.savedAni + str(i) + '.png'
            plt.savefig(savename,format='png')

class AnimateField(object):  # FIXME -- should this be streamplot instead of quiver?
    """
    Animate a vector field defined on cell centres and pertaining to body forces in cell.

    """

    def __init__(self,sim,Fx_time,Fy_time,cells,p,ani_repeat = True, save = True, title = 'Force field',
        saveFolder = 'animation/Ffield',saveFile = 'Ffield_'):

        self.fig = plt.figure()
        self.ax = plt.subplot(111)
        self.p = p
        self.sim = sim
        self.cells = cells
        self.save = save
        self.Fx_time = Fx_time
        self.Fy_time = Fy_time
        self.tit = title

        self.saveFolder = saveFolder
        self.saveFile = saveFile
        self.ani_repeat = ani_repeat

        xmin = cells.xmin*p.um
        xmax = cells.xmax*p.um
        ymin = cells.ymin*p.um
        ymax = cells.ymax*p.um

        if self.save is True:

            set_up_filesave(self,p)

        Fx = (1/p.um)*self.Fx_time[0]
        Fy = (1/p.um)*self.Fy_time[0]

        F = np.sqrt(Fx**2 + Fy**2)

        cmin = F.min()
        cmax = F.max()

        if F.all() != 0.0:

            Fx = Fx/F
            Fy = Fy/F

        if p.showCells is True:

            # define a polygon collection based on individual cell polygons
            self.points = np.multiply(cells.cell_verts, p.um)
            self.fmesh =  PolyCollection(self.points, cmap=p.default_cm, edgecolors='none')
            self.fmesh.set_array(F)
            self.ax.add_collection(self.fmesh)

        else:
            # interpolate the data to the grid:

            dat_grid = interpolate.griddata((cells.cell_centres[:, 0],cells.cell_centres[:, 1]),F,
                                        (cells.Xgrid,cells.Ygrid),method=p.interp_type, fill_value=0)

            dat_grid = np.multiply(dat_grid,cells.maskM)

            if p.plotMask is True:
                dat_grid = ma.masked_array(dat_grid, np.logical_not(cells.maskM))

            self.fmesh = plt.imshow(dat_grid,origin='lower',extent=[xmin,xmax,ymin,ymax],cmap=p.default_cm)

        # add a vector plot of force components:
        self.vplot = self.ax.quiver(p.um*cells.cell_centres[:,0],p.um*cells.cell_centres[:,1],
            Fx,Fy,zorder = 10)

        self.ax.axis('equal')

        self.ax.axis([xmin,xmax,ymin,ymax])

        if p.autoscale_force_ani is False:
            self.fmesh.set_clim(p.force_ani_min_clr,p.force_ani_max_clr)
        else:
            self.fmesh.set_clim(cmin,cmax)

        cb = self.fig.colorbar(self.fmesh)

        self.ax.set_title(title)
        self.ax.set_xlabel('Spatial distance [um]')
        self.ax.set_ylabel('Spatial distance [um]')
        cb.set_label('Body Force [N/cm3]')

        self.frames = len(sim.time)
        ani = animation.FuncAnimation(self.fig, self.aniFunc,
            frames=self.frames, interval=100, repeat=self.ani_repeat)

        show_plot(p)


    #FIXME: This function appears to be broken. For example, "p" is undefined
    #below, so the line "Fx = (1/p.um)*self.Fx_time[i]" is guaranteed to fail.
    #River rapids coalesce downstream, abreast of the old swimming hole!
    def aniFunc(self,i):

        titani = self.tit + ' (simulation time' + ' ' + str(round(self.sim.time[i],3)) + ' ' + ' s)'
        self.ax.set_title(titani)

        # get the appropriate data and scale it:

        Fx = (1/p.um)*self.Fx_time[i]
        Fy = (1/p.um)*self.Fy_time[i]

        F = np.sqrt(Fx**2 + Fy**2)

        if F.all() != 0.0:

            Fx = Fx/F
            Fy = Fy/F

        if self.p.showCells is True:

            self.fmesh.set_array(F)

        else:
            self.fmesh.set_data(F)

        self.fmesh.set_clim(F.min(),F.max())
        self.vplot.set_UVC(Fx,Fy)

        if self.save is True:
            self.fig.canvas.draw()
            savename = self.savedAni + str(i) + '.png'
            plt.savefig(savename,format='png')

class AnimateVelocity(object):

    def __init__(self,sim,cells,p,ani_repeat = True, save = True, saveFolder = 'animation/Velocity',
        saveFile = 'Velocity_'):

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

            set_up_filesave(self,p)

        if p.sim_ECM is True and p.ani_Velocity_type == 'ECM':

            vfield = np.sqrt(sim.u_env_x_time[0]**2 + sim.u_env_y_time[0]**2)*1e9

            self.msh = self.ax.imshow(vfield,origin='lower', extent = [cells.xmin*p.um, cells.xmax*p.um,
                cells.ymin*p.um, cells.ymax*p.um], cmap=p.default_cm)

            vnorm = np.max(vfield)

            self.streamV = self.ax.quiver(p.um*cells.xypts[:,0], p.um*cells.xypts[:,1],
                    sim.u_env_x_time[-1].ravel()/vnorm,sim.u_env_y_time[-1].ravel()/vnorm)

            tit_extra = 'Extracellular'

        elif p.ani_Velocity_type == 'GJ' or p.sim_ECM is True:

            ugjx = sim.u_cells_x_time[0]
            ugjy = sim.u_cells_y_time[0]

            v_gj_x = interpolate.griddata((cells.cell_centres[:,0],cells.cell_centres[:,1]),ugjx,(cells.Xgrid,cells.Ygrid),
                                          fill_value=0,method=p.interp_type)

            v_gj_x = v_gj_x*cells.maskM

            v_gj_y = interpolate.griddata((cells.cell_centres[:,0],cells.cell_centres[:,1]),ugjy,(cells.Xgrid,cells.Ygrid),
                                          fill_value=0,method=p.interp_type)

            v_gj_y = v_gj_y*cells.maskM

            vfield = np.sqrt(v_gj_x**2 + v_gj_y**2)*1e9

            self.msh = self.ax.imshow(vfield,origin='lower', extent = [cells.xmin*p.um, cells.xmax*p.um,
                cells.ymin*p.um, cells.ymax*p.um],cmap=p.default_cm)

            vnorm = np.max(vfield)


            lw = (3.0*vfield/vnorm) + 0.5

            self.streamV = self.ax.streamplot(cells.Xgrid*p.um,cells.Ygrid*p.um,v_gj_x/vnorm,v_gj_y/vnorm,density=p.stream_density,
                    linewidth=lw,color='k',cmap=p.default_cm,arrowsize=1.5)

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

        show_plot(p)


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

            u_gj_x = interpolate.griddata((self.cells.cell_centres[:,0],self.cells.cell_centres[:,1]),
            ugjx,(self.cells.Xgrid,self.cells.Ygrid), fill_value=0,method=self.p.interp_type)

            u_gj_x = u_gj_x*self.cells.maskM

            u_gj_y = interpolate.griddata((self.cells.cell_centres[:,0],self.cells.cell_centres[:,1]),
                ugjy,(self.cells.Xgrid,self.cells.Ygrid), fill_value=0,method=self.p.interp_type)

            u_gj_y = u_gj_y*self.cells.maskM

            vfield = np.sqrt(u_gj_x**2 + u_gj_y**2)*1e9

            self.msh.set_data(vfield)

            vnorm = np.max(vfield)

            self.streamV.lines.remove()
            self.ax.patches = []

            lw = (3.0*vfield/vnorm) + 0.5

            self.streamV = self.ax.streamplot(self.cells.Xgrid*self.p.um,self.cells.Ygrid*self.p.um,u_gj_x/vnorm,u_gj_y/vnorm,
                density=self.p.stream_density,linewidth=lw,color='k',cmap=self.p.default_cm,arrowsize=1.5)

            # self.streamV.set_UVC(u_gj_x/vnorm,u_gj_y/vnorm)

        cmax = np.max(vfield)

        if self.p.autoscale_Velocity_ani is True:
            self.msh.set_clim(0,cmax)

        if self.save is True:
            self.fig.canvas.draw()
            savename = self.savedAni + str(i) + '.png'
            plt.savefig(savename,format='png')

class AnimateDeformation(object):

    def __init__(self,sim,cells,p,ani_repeat = True, save = True, saveFolder = 'animation/Deformation',
        saveFile = 'Deformation_'):

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

            set_up_filesave(self,p)

        if self.p.ani_Deformation_type == 'Vmem':

            dx = self.sim.dx_cell_time[0]
            dy = self.sim.dy_cell_time[0]

            if self.p.sim_ECM is False:

                dd = self.sim.vm_time[0]*1e3

            else:
                dd = self.sim.vcell_time[0]*1e3

        elif self.p.ani_Deformation_type == 'Displacement':

            dx = self.sim.dx_cell_time[0]
            dy = self.sim.dy_cell_time[0]

            dd = p.um*np.sqrt(dx**2 + dy**2)

        points = np.multiply(sim.cell_verts_time[0], p.um)
        dd_collection = PolyCollection(points, array=dd, cmap=p.default_cm, edgecolors='none')
        self.ax.add_collection(dd_collection)

        self.ax.quiver(p.um*sim.cell_centres_time[0][:,0],p.um*sim.cell_centres_time[0][:,1],dx,dy)

        self.ax.axis('equal')

        xmin = cells.xmin*p.um
        xmax = cells.xmax*p.um
        ymin = cells.ymin*p.um
        ymax = cells.ymax*p.um

        self.ax.axis([xmin,xmax,ymin,ymax])


        if p.autoscale_Deformation_ani is True:

            # first flatten the data (needed in case cells were cut)
            all_z = []
            for zarray in sim.dx_cell_time:
                for val in zarray:
                    all_z.append(val)

            for zarray in sim.dy_cell_time:
                for val in zarray:
                    all_z.append(val)

            self.cmin = p.um*np.min(all_z)
            self.cmax = p.um*np.max(all_z)

            dd_collection.set_clim(self.cmin,self.cmax)


        elif p.autoscale_Deformation_ani is False:
            dd_collection.set_clim(p.Deformation_ani_min_clr,p.Deformation_ani_max_clr)

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

        show_plot(p)


    def aniFunc(self,i):

        # we need to have changing cells, so we have to clear the plot and redo it...
        self.fig.clf()
        self.ax = plt.subplot(111)

        if self.p.ani_Deformation_type == 'Vmem':

            dx = self.sim.dx_cell_time[i]
            dy = self.sim.dy_cell_time[i]

            if self.p.sim_ECM is False:

                dd = self.sim.vm_time[i]*1e3

            else:
                dd = self.sim.vcell_time[i]*1e3

        elif self.p.ani_Deformation_type == 'Displacement':

            dx = self.sim.dx_cell_time[i]
            dy = self.sim.dy_cell_time[i]

            dd = 1e6*np.sqrt(dx**2 + dy**2)

            if dd.all() != 0:

                dx = dx/dd
                dy = dy/dd



        points = np.multiply(self.sim.cell_verts_time[i], self.p.um)
        dd_collection = PolyCollection(points, array=dd, cmap=self.p.default_cm, edgecolors='none')
        self.ax.add_collection(dd_collection)

        self.ax.quiver(self.p.um*self.sim.cell_centres_time[i][:,0],self.p.um*self.sim.cell_centres_time[i][:,1],dx,dy)

        self.ax.axis('equal')

        xmin = self.cells.xmin*self.p.um
        xmax = self.cells.xmax*self.p.um
        ymin = self.cells.ymin*self.p.um
        ymax = self.cells.ymax*self.p.um

        self.ax.axis([xmin,xmax,ymin,ymax])

        if self.p.autoscale_Deformation_ani is False:

            dd_collection.set_clim(self.p.Deformation_ani_min_clr,self.p.Deformation_ani_max_clr)

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

    def __init__(self,sim,cells,time,p,save=True,ani_repeat=False,clrAutoscale=True,
    clrMin = None,clrMax = None,clrmap = cm.rainbow, number_cells=False,saveFolder = 'animation/Venv',saveFile = 'venv_'):

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

            set_up_filesave(self,p)

        if clrAutoscale is False:
            self.cmin = clrMin
            self.cmax = clrMax

        self.meshplot = plt.imshow(sim.venv_time[0].reshape(cells.X.shape)*1000,
                                   origin='lower',extent=[xmin,xmax,ymin,ymax], cmap=p.default_cm)


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

        show_plot(p)


    def aniFunc(self,i):

        titani = 'Environmental Voltage' + ' (simulation time' + ' ' + str(round(self.sim.time[i],3)) + ' ' + ' s)'
        self.ax.set_title(titani)

        self.meshplot.set_data(self.sim.venv_time[i].reshape(self.cells.X.shape)*1000)

        if self.save is True:
            self.fig.canvas.draw()
            savename = self.savedAni + str(i) + '.png'
            plt.savefig(savename,format='png')

class AnimateMem(object):
    """
    Animates the channel or pump density factor (sim.rho_channel or sim.rho_pump) which changes due to
    electroosmotic/electrophoretic movements due to self-generated fields and flows in the cluster.

    """

    def __init__(self,sim,cells,time,p,save=False,ani_repeat=False,current_overlay=False,
        clrAutoscale = True, clrMin = None, clrMax = None,
        number_cells = False, saveFolder = 'animation/pump_electroosmo', saveFile = 'rhoPump_'):

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

            set_up_filesave(self,p)

        self.bkgBool = False

        cell_edges_flat = cells.um*cells.mem_edges_flat

        self.coll = LineCollection(cell_edges_flat, array=sim.rho_pump_time[0], cmap=self.clrmap,linewidths=4.0)
        self.ax.add_collection(self.coll)

        self.ax.axis('equal')

        self.ax.axis([xmin,xmax,ymin,ymax])


        if self.current_overlay is True:

            if p.sim_ECM is False or p.IecmPlot is False:

                Jmag_M = np.sqrt(sim.I_gj_x_time[0]**2 + sim.I_gj_y_time[0]**2) + 1e-30

                J_x = sim.I_gj_x_time[0]/Jmag_M
                J_y = sim.I_gj_y_time[0]/Jmag_M

                lw = (3.0*Jmag_M/Jmag_M.max()) + 0.5

                self.streams = self.ax.streamplot(cells.Xgrid*p.um,cells.Ygrid*p.um,J_x,J_y,
                    density=self.density,linewidth=lw,color='k',
                    cmap=self.clrmap,arrowsize=1.5)

                self.tit_extra = 'Gap junction current'

            elif p.IecmPlot is True:

                Jmag_M = np.sqrt(sim.I_tot_x_time[0]**2 + sim.I_tot_y_time[0]**2) + 1e-30

                J_x = sim.I_tot_x_time[0]/Jmag_M
                J_y = sim.I_tot_y_time[0]/Jmag_M

                lw = (3.0*Jmag_M/Jmag_M.max()) + 0.5

                self.streams = self.ax.streamplot(cells.Xgrid*p.um,cells.Ygrid*p.um,J_x,J_y,density=self.density,linewidth=lw,color='k',
                    cmap=self.clrmap,arrowsize=1.5)

                self.tit_extra = 'Extracellular current overlay'

        else:

            self.tit_extra = ' '

        # set range of the colormap

        if clrAutoscale is True:
            # first flatten the data (needed in case cells were cut)
            all_z = []
            for zarray in sim.rho_pump_time:
                for val in zarray:
                    all_z.append(val)

            self.cmin = np.min(all_z)
            self.cmax = np.max(all_z)


        elif clrAutoscale is False:
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

        show_plot(p)


    def aniFunc(self,i):

        zz = self.sim.rho_pump_time[i]

        self.coll.set_array(zz)

        if self.current_overlay is True:

            if self.sim_ECM is False or self.IecmPlot is False:

                Jmag_M = np.sqrt(self.sim.I_gj_x_time[i]**2 + self.sim.I_gj_y_time[i]**2) + 1e-30

                J_x = self.sim.I_gj_x_time[i]/Jmag_M
                J_y = self.sim.I_gj_y_time[i]/Jmag_M

                lw = (3.0*Jmag_M/Jmag_M.max()) + 0.5

                self.streams.lines.remove()
                self.ax.patches = []

                self.streams = self.ax.streamplot(self.cells.Xgrid*1e6,self.cells.Ygrid*1e6,J_x,J_y,
                    density=self.density,linewidth=lw,color='k', cmap=self.colormap,arrowsize=1.5)

            elif self.IecmPlot is True:

                Jmag_M = np.sqrt(self.sim.I_tot_x_time[i]**2 + self.sim.I_tot_y_time[i]**2) + 1e-30

                J_x = self.sim.I_tot_x_time[i]/Jmag_M
                J_y = self.sim.I_tot_y_time[i]/Jmag_M

                lw = (3.0*Jmag_M/Jmag_M.max()) + 0.5

                self.streams.lines.remove()
                self.ax.patches = []

                self.streams = self.ax.streamplot(self.cells.Xgrid*1e6,self.cells.Ygrid*1e6,
                    J_x,J_y,density=self.density,linewidth=lw,color='k', cmap=self.colormap,arrowsize=1.5)

        titani = self.tit_extra + ' (sim time' + ' ' + str(round(self.time[i],3)) + ' ' + ' s)'
        self.ax.set_title(titani)

        if self.save is True:
            self.fig.canvas.draw()
            savename = self.savedAni + str(i) + '.png'
            plt.savefig(savename,format='png')

class AnimateDyeData(object):
    """
    Animate morphogen concentration data in cell and environment as a function of time.

    """

    def __init__(self,sim,cells,p, save=False,ani_repeat=False,current_overlay=False,
        clrAutoscale = True, clrMin = None, clrMax = None, clrmap = cm.rainbow,
        number_cells = False, saveFolder = 'animation', saveFile = 'sim_'):

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

            set_up_filesave(self,p)

        self.bkgPlot = self.ax.imshow(self.zenv_t[0].reshape(cells.X.shape),origin='lower',
            extent= [xmin,xmax,ymin,ymax],cmap=clrmap)

        # define a polygon collection based on individual cell polygons
        self.points = np.multiply(cells.cell_verts, p.um)
        self.collection =  PolyCollection(self.points, cmap=self.colormap, edgecolors='none')
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

            if cmaxa > cmaxb:
                self.cmax = cmaxa
            else:
                self.cmax = cmaxb

            if cmina < cminb:
                self.cmin = cmina
            else:
                self.cmin = cminb

        elif clrAutoscale is False:
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

        show_plot(p)


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

def plotSingleCellVData(sim,celli,p,fig=None,ax=None, lncolor='k'):

    tvect_data=[x[celli]*1000 for x in sim.vm_time]

    if fig is None:
        fig = plt.figure()# define the figure and axes instances
    if ax is None:
        ax = plt.subplot(111)

    ax.plot(sim.time, tvect_data,lncolor,linewidth=2.0)

    if p.GHK_calc is True:

        tvect_data_ghk = [x[p.plot_cell]*1000 for x in sim.vm_GHK_time]
        ax.plot(sim.time, tvect_data_ghk,'r',linewidth=2.0)


    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Voltage [mV]')

    return fig, ax

def plotSingleCellCData(simdata_time,simtime,ioni,celli,fig=None,ax=None,lncolor='b',ionname='ion'):


    # ccIon_cell = [arr[ioni][celli] for arr in simdata_time]
    ccIon_cell = []

    for carray in simdata_time:

        conc = carray[ioni][celli]
        ccIon_cell.append(conc)

    if fig is None:
        fig = plt.figure()# define the figure and axes instances
    if ax is None:
        ax = plt.subplot(111)
        #ax = plt.axes()

    lab = ionname

    # xmin = simtime[0]
    # xmax = simtime[-1]
    # ymin = np.min(ccIon_cell)
    # ymax = np.max(ccIon_cell)

    ax.plot(simtime, ccIon_cell,lncolor,label=lab)
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Concentration [mol/m3]')

    return fig, ax

def plotSingleCellData(simtime,simdata_time,celli,fig=None,ax=None,lncolor='b',lab='Data'):

    data_cell = [arr[celli] for arr in simdata_time]

    if fig is None:
        fig = plt.figure()# define the figure and axes instances
    if ax is None:
        ax = plt.subplot(111)

    xmin = simtime[0]
    xmax = simtime[-1]
    ymin = np.min(data_cell)
    ymax = np.max(data_cell)

    ax.plot(simtime, data_cell,lncolor,label=lab)
    ax.set_xlabel('Time [s]')
    ax.set_ylabel(lab)

    return fig, ax

def plotFFT(simtime,simdata_time,celli,fig=None,ax=None,lncolor='b',lab='Data'):
    """
    Calculates the FFT for time-series data defined on a single cell (p.plot_cell)
    and returns a plot of the spectrum in frequency-space.

    Parameters
    -----------
    simtime:            The time vector for the plot
    simdata_time:       The full data for the plot define on all cell centres of the cluster at all sample times
    celli:              The single cell index to extract data for
    fig:                A handle to an existing figure (default None; function creates the figure)
    ax:                 Handle to an existing axis (default None; function creates the axis)
    lncolor:            Line colour for the plot
    lab:                Label for the data type (e.g. "Vmem [V]" or "Body Force [N/m3]")

    Returns
    --------
    fig, ax     Handles to the figure and axis of the FFT plot
    """

    if fig is None:
        fig = plt.figure()
    if ax is None:
        ax = plt.subplot(111)

    sample_size = len(simtime)
    sample_spacing = simtime[1] - simtime[0]

    cell_data_o = [arr[celli] for arr in simdata_time]
    # cell_data = ((1/sample_size)*(cell_data_o/np.mean(cell_data_o)) )   # normalize the signal
    cell_data = (1/sample_size)*(cell_data_o - np.mean(cell_data_o))

    f_axis = np.fft.rfftfreq(sample_size, d=sample_spacing)
    fft_data_o = np.fft.rfft(cell_data)
    fft_data = np.sqrt(np.real(fft_data_o)**2 + np.imag(fft_data_o)**2)

    xmin = f_axis[0]
    xmax = f_axis[-1]
    ymin = np.min(fft_data)
    ymax = np.max(fft_data)

    ax.plot(f_axis,fft_data)
    ax.axis([xmin,xmax,ymin,ymax])

    ax.set_xlabel('Frequency [1/s]')
    ax.set_ylabel('Signal Power')

    return fig, ax

def plotHetMem(sim,cells, p, fig=None, ax=None, zdata=None,clrAutoscale = True, clrMin = None, clrMax = None,
    clrmap=None,edgeOverlay = True,pointOverlay=None, number_cells = False, number_mems = False,
    number_ecm = False, current_overlay = False,plotIecm = False):
        """
        This plotting method assigns color-data to each node in the cell cluster that has distinct
        membrane domains for each cell. Data is interpolated to generate a smooth surface plot.
        The method returns a plot instance (fig, axes)

        When using p.sim_ECM, this plotting method overrides both plotPolyData and plotCellData.

        Parameters
        ----------
        zdata                  A data array with each scalar entry corresponding to a point in
                               mem_mids_flat. If not specified the default is z=1. If 'random'
                               is specified the method creates random vales from 0 to 1..

        clrmap                 The colormap to use for plotting. Must be specified as cm.mapname. A list of
                               available mapnames is supplied at
                               http://matplotlib.org/examples/color/colormaps_reference.html

        clrAutoscale           If True, the colorbar is autoscaled to the max and min of zdata.

        clrMin                 Sets the colorbar to a user-specified minimum value.

        clrMax                 Set the colorbar to a user-specified maximum value


        edgeOverlay             This option allows the user to specify whether or not they want cell edges overlayed.
                                Default is False, set to True to use.

        pointOverlay            This option allows user to specify whether or not they want cell_centre points plotted
                                Default is False, set to True to use.

        number_cells,           Booleans that control whether or not cell, membrane, and ecm spaces are labeled
        number_ecm,             with their indices.
        number_mems


        Returns
        -------
        fig, ax                Matplotlib figure and axes instances for the plot.

        Notes
        -------
        Uses `matplotlib.pyplot` and `numpy` arrays. With `edgeOverlay` and
        `pointOverlay` equal to `None`, this is computationally fast and *is*
        recommended for plotting data on large collectives.
        """

        if fig is None:
            fig = plt.figure()# define the figure and axes instances
        if ax is None:
            ax = plt.subplot(111)

        if clrmap is None:
            clrmap = p.default_cm

        if zdata is None:
            zdata = np.ones((p.plot_grid_size,p.plot_grid_size))

        ax.axis('equal')

        xmin = cells.xmin*p.um
        xmax = cells.xmax*p.um
        ymin = cells.ymin*p.um
        ymax = cells.ymax*p.um

        ax.axis([xmin,xmax,ymin,ymax])

        if p.plotMask is True:
            zdata = ma.masked_array(zdata, np.logical_not(cells.maskM))

        meshplt = plt.imshow(zdata,origin='lower',extent=[xmin,xmax,ymin,ymax],cmap=clrmap)

        if pointOverlay is True:
            ax.scatter(
                p.um*cells.mem_mids_flat[:,0],
                p.um*cells.mem_mids_flat[:,1], c='k',)

        if edgeOverlay is True:
            # cell_edges_flat, _ , _= tb.flatten(cells.mem_edges)
            cell_edges_flat = cells.um*cells.mem_edges_flat
            coll = LineCollection(cell_edges_flat,colors='k')
            coll.set_alpha(0.5)
            ax.add_collection(coll)

        if zdata is not None:
            # Add a colorbar for the mesh plot:
            maxval = round(np.max(1000*sim.vm_time[-1]),1)
            minval = round(np.min(1000*sim.vm_time[-1]),1)
            checkval = maxval - minval

            if checkval == 0:
                minval = minval - 0.1
                maxval = maxval + 0.1

        if zdata is not None and clrAutoscale is True:
            meshplt.set_clim(minval,maxval)
            ax_cb = fig.colorbar(meshplt,ax=ax)

        elif clrAutoscale is False:

            meshplt.set_clim(clrMin,clrMax)
            ax_cb = fig.colorbar(meshplt,ax=ax)

        else:
            ax_cb = None

        if number_cells is True:

            for i,cll in enumerate(cells.cell_centres):
                ax.text(p.um*cll[0],p.um*cll[1],i,ha='center',va='center')

        if number_mems is True:

            for i,mem in enumerate(cells.mem_mids_flat):
                ax.text(p.um*mem[0],p.um*mem[1],i,ha='center',va='center')

        if current_overlay is True:

            I_overlay(sim,cells,p,ax,clrmap,plotIecm)

        return fig, ax, ax_cb

def plotPolyData(sim, cells, p, fig=None, ax=None, zdata = None, clrAutoscale = True, clrMin = None, clrMax = None,
    clrmap = None, number_cells=False, current_overlay = False,plotIecm=False):
        """
        Assigns color-data to each polygon in a cell cluster diagram and returns a plot instance (fig, axes)

        Parameters
        ----------
        cells                  Data structure holding all world information about cell geometry

        zdata                  A data array with each scalar entry corresponding to a cell's data value
                               (for instance, concentration or voltage). If zdata is not supplied, the
                               cells will be plotted with a uniform color; if zdata = random a random
                               data set will be created and plotted.

        clrAutoscale           If True, the colorbar is autoscaled to the max and min of zdata.

        clrMin                 Sets the colorbar to a user-specified minimum value.

        clrMax                 Set the colorbar to a user-specified maximum value

        clrmap                 The colormap to use for plotting. Must be specified as cm.mapname. A list of
                               available mapnames is supplied at
                               http://matplotlib.org/examples/color/colormaps_reference.html


        Returns
        -------
        fig, ax                Matplotlib figure and axes instances for the plot.

        Notes
        -------
        Uses matplotlib.collections PolyCollection, matplotlib.cm, matplotlib.pyplot and numpy arrays
        Computationally slow -- not recommended for large collectives (500 x 500 um max)
        """
        if fig is None:
            fig = plt.figure()# define the figure and axes instances
        if ax is None:
            ax = plt.subplot(111)
            #ax = plt.axes()

        if zdata is None:  # if user doesn't supply data
            z = np.ones(len(cells.cell_verts)) # create flat data for plotting

        elif zdata == 'random':  # if user doesn't supply data
            z = np.random.random(len(cells.cell_verts)) # create some random data for plotting

        else:
            z = zdata

        # Make the polygon collection and add it to the plot.
        if clrmap is None:
            #clrmap = p.default_cm
            clrmap = cm.rainbow

        points = np.multiply(cells.cell_verts, p.um)

        coll = PolyCollection(points, array=z, cmap=clrmap, edgecolors='none')
        ax.add_collection(coll)
        ax.axis('equal')

        # Add a colorbar for the PolyCollection
        if zdata is not None:
            maxval = np.max(zdata,axis=0)
            minval = np.min(zdata,axis=0)
            # checkval = maxval - minval
            #
            # if checkval == 0:
            #     minval = minval - 0.1
            #     maxval = maxval + 0.1

        if zdata is not None and clrAutoscale is True:
            coll.set_clim(minval,maxval)
            ax_cb = fig.colorbar(coll,ax=ax)

        elif clrAutoscale is False:

            coll.set_clim(clrMin,clrMax)
            ax_cb = fig.colorbar(coll,ax=ax)

        elif zdata is None:
            ax_cb = None

        if number_cells is True:
            for i,cll in enumerate(cells.cell_centres):
                ax.text(p.um*cll[0],p.um*cll[1],i,ha='center',va='center')

        if current_overlay is True:

            I_overlay(sim,cells,p,ax,clrmap,plotIecm)

        xmin = cells.xmin*p.um
        xmax = cells.xmax*p.um
        ymin = cells.ymin*p.um
        ymax = cells.ymax*p.um

        ax.axis([xmin,xmax,ymin,ymax])

        return fig,ax,ax_cb

def plotCellData(sim,cells, p, fig=None, ax=None, zdata=None,clrAutoscale = True, clrMin = None, clrMax = None,
    clrmap=None,edgeOverlay = None,pointOverlay=None,number_cells = False,current_overlay=False,plotIecm=False):
        """
        The work-horse of pre-defined plotting methods, this method assigns color-data to each node in cell_centres
        and interpolates data to generate a smooth surface plot. The method returns a plot instance (fig, axes)

        Parameters
        ----------
        zdata_t                  A data array with each scalar entry corresponding to a point in
                               cell_centres. If not specified the default is z=1. If 'random'
                               is specified the method creates random vales from 0 to 1..

        clrmap                 The colormap to use for plotting. Must be specified as cm.mapname. A list of
                               available mapnames is supplied at
                               http://matplotlib.org/examples/color/colormaps_reference.html

        clrAutoscale           If True, the colorbar is autoscaled to the max and min of zdata.

        clrMin                 Sets the colorbar to a user-specified minimum value.

        clrMax                 Set the colorbar to a user-specified maximum value


        edgeOverlay             This option allows the user to specify whether or not they want cell edges overlayed.
                                Default is False, set to True to use.

        pointOverlay            This option allows user to specify whether or not they want cell_centre points plotted
                                Default is False, set to True to use.


        Returns
        -------
        fig, ax                Matplotlib figure and axes instances for the plot.

        Notes
        -------
        Uses `matplotlib.pyplot` and `numpy` arrays. With `edgeOverlay` and
        `pointOverlay` equal to `None`, this is computationally fast and *is*
        recommended for plotting data on large collectives.
        """

        if current_overlay is None:
            current_overlay = p.I_overlay

        if fig is None:
            fig = plt.figure()# define the figure and axes instances
        if ax is None:
            ax = plt.subplot(111)
            #ax = plt.axes()

        if zdata is None:  # if user doesn't supply data
            z = np.ones(len(cells.cell_centres)) # create flat data for plotting

        elif zdata == 'random':  # if user doesn't supply data
            z = np.random.random(len(cells.cell_centres)) # create some random data for plotting

        else:
            z = zdata

        if clrmap is None:
            clrmap = p.default_cm

        ax.axis('equal')

        xmin = cells.xmin*p.um
        xmax = cells.xmax*p.um
        ymin = cells.ymin*p.um
        ymax = cells.ymax*p.um

        ax.axis([xmin,xmax,ymin,ymax])

        dat_grid = interpolate.griddata((cells.cell_centres[:, 0],cells.cell_centres[:, 1]),z,
                                        (cells.Xgrid,cells.Ygrid),method=p.interp_type)
        dat_grid = np.nan_to_num(dat_grid)
        dat_grid = np.multiply(dat_grid,cells.maskM)

        if p.plotMask is True:
            dat_grid = ma.masked_array(dat_grid, np.logical_not(cells.maskM))

        triplt = plt.imshow(dat_grid,origin='lower',extent=[xmin,xmax,ymin,ymax],cmap=clrmap)

        # ax.axis('equal')

         # Add a colorbar for the triplot:

        maxval = np.max(z,axis=0)
        minval = np.min(z,axis=0)

        if zdata is not None and clrAutoscale is True:
            triplt.set_clim(minval,maxval)
            ax_cb = fig.colorbar(triplt,ax=ax)
        elif clrAutoscale is False:
            triplt.set_clim(clrMin,clrMax)
            ax_cb = fig.colorbar(triplt,ax=ax)
        elif zdata is None:
            ax_cb = None

        if pointOverlay is True:
            ax.scatter(p.um*cells.cell_centres[:,0],p.um*cells.cell_centres[:,1], c=z,cmap=clrmap)

        if edgeOverlay is True:
            # cell_edges_flat, _ , _= tb.flatten(cells.mem_edges)
            cell_edges_flat = cells.um*cells.mem_edges_flat
            coll = LineCollection(cell_edges_flat,colors='k')
            coll.set_alpha(0.5)
            ax.add_collection(coll)

        if current_overlay is True:
            I_overlay(sim,cells,p,ax,clrmap,plotIecm)

        if number_cells is True:
            for i,cll in enumerate(cells.cell_centres):
                ax.text(p.um*cll[0],p.um*cll[1],i,ha='center',va='center')

        return fig, ax, ax_cb

def plotEfield(sim,cells,p):

    fig = plt.figure()
    ax = plt.subplot(111)

    if p.sim_ECM is True:

        fig_2 = plt.figure()
        ax_2 = plt.subplot(111)

        efield = np.sqrt(sim.efield_ecm_x_time[-1]**2 + sim.efield_ecm_y_time[-1]**2)
        msh_2 = ax_2.imshow(efield,origin='lower', extent = [cells.xmin*p.um, cells.xmax*p.um, cells.ymin*p.um,
            cells.ymax*p.um],cmap=p.background_cm)

        ax_2.quiver(p.um*cells.xypts[:,0], p.um*cells.xypts[:,1], sim.efield_ecm_x_time[-1].ravel(),
            sim.efield_ecm_y_time[-1].ravel())

        tit_extra = 'Extracellular'

        ax_2.axis('equal')

        xmin = cells.xmin*p.um
        xmax = cells.xmax*p.um
        ymin = cells.ymin*p.um
        ymax = cells.ymax*p.um

        ax_2.axis([xmin,xmax,ymin,ymax])

        if p.autoscale_Efield is False:
            msh_2.set_clim(p.Efield_min_clr,p.Efield_max_clr)

        cb_2 = fig_2.colorbar(msh_2)

        tit = "Final Electric Field in " + tit_extra + ' Spaces'
        ax_2.set_title(tit)
        ax_2.set_xlabel('Spatial distance [um]')
        ax_2.set_ylabel('Spatial distance [um]')
        cb_2.set_label('Electric Field [V/m]')

    else:
        ax_2 = None
        fig_2 = None
        cb_2 = None

    # Plot field intracellularly:

    E_gj_x = np.dot(cells.M_sum_mems,sim.efield_gj_x_time[-1])/cells.num_mems
    E_gj_y = np.dot(cells.M_sum_mems,sim.efield_gj_y_time[-1])/cells.num_mems

    efield = np.sqrt(E_gj_x**2 + E_gj_y**2)

    # normalize the vector components:
    if efield.all() != 0:
        E_gj_x = E_gj_x/efield
        E_gj_y = E_gj_y/efield

    efield_grid = np.zeros(len(cells.voronoi_centres))
    efield_grid[cells.cell_to_grid] = efield

    msh = ax.tripcolor(p.um*cells.voronoi_centres[:,0],p.um*cells.voronoi_centres[:,1],efield_grid,
        shading='gouraud',cmap=p.background_cm)

    ax.quiver(p.um*cells.cell_centres[:,0],p.um*cells.cell_centres[:,1],E_gj_x,E_gj_y,
        pivot='mid',color = p.vcolor, units='x',headwidth=5, headlength = 7)

    tit_extra = 'Intracellular'

    ax.axis('equal')

    xmin = cells.xmin*p.um
    xmax = cells.xmax*p.um
    ymin = cells.ymin*p.um
    ymax = cells.ymax*p.um

    ax.axis([xmin,xmax,ymin,ymax])

    if p.autoscale_Efield is False:
        msh.set_clim(p.Efield_min_clr,p.Efield_max_clr)

    cb = fig.colorbar(msh)

    tit = "Final Electric Field in " + tit_extra + ' Spaces'
    ax.set_title(tit)
    ax.set_xlabel('Spatial distance [um]')
    ax.set_ylabel('Spatial distance [um]')
    cb.set_label('Electric Field [V/m]')

    return fig, ax, cb, fig_2, ax_2, cb_2

def plotMemData(cells, p, fig= None, ax = None, zdata=None,clrmap=None):
        """

        Assigns color-data to edges in a 2D Voronoi diagram and returns a plot instance (fig, axes)

        Parameters
        ----------
        zdata_t                  A data array with each scalar entry corresponding to a polygon entry in
                               vor_verts. If not specified the default is z=1. If 'random'
                               is specified the method creates random vales from 0 to 1..

        clrmap                 The colormap to use for plotting. Must be specified as cm.mapname. A list of
                               available mapnames is supplied at
                               http://matplotlib.org/examples/color/colormaps_reference.html
                               Default is cm.rainbow. Good options are cm.coolwarm, cm.Blues, cm.jet


        Returns
        -------
        fig, ax                Matplotlib figure and axes instances for the plot.

        Notes
        -------
        Uses matplotlib.collections LineCollection, matplotlib.cm, matplotlib.pyplot and numpy arrays
        Computationally slow -- not recommended for large collectives (500 x 500 um max)

        """

        if fig is None:
            fig = plt.figure()# define the figure and axes instances
        if ax is None:
            ax = plt.subplot(111)

        cell_edges_flat = cells.um*cells.mem_edges_flat

        if zdata is None:
            z = np.ones(len(cell_edges_flat))
        elif zdata == 'random':
            z = np.random.random(len(cell_edges_flat))
        else:
            z = zdata

        if clrmap is None:
            clrmap = cm.rainbow

        coll = LineCollection(cell_edges_flat, array=z, cmap=clrmap,linewidths=4.0)
        ax.add_collection(coll)

        # coll.set_clim(0,3)

        ax.axis('equal')

        # Add a colorbar for the Line Collection
        if zdata is not None:
            ax_cb = fig.colorbar(coll, ax=ax)

        ax.axis('equal')

        xmin = cells.xmin*p.um
        xmax = cells.xmax*p.um
        ymin = cells.ymin*p.um
        ymax = cells.ymax*p.um

        ax.axis([xmin,xmax,ymin,ymax])

        return fig, ax, ax_cb

def plotConnectionData(cells, p, fig = None, ax=None, zdata=None,clrmap=None,colorbar = None, pickable=None):
        """
        Assigns color-data to connections between a cell and its nearest neighbours and returns plot instance

        Parameters
        ----------

        zdata_t                  A data array with each scalar entry corresponding to a polygon entry in
                               vor_verts. If not specified the default is z=1. If 'random'
                               is specified the method creates random vales from 0 to 1..

        clrmap                 The colormap to use for plotting. Must be specified as cm.mapname. A list of
                               available mapnames is supplied at
                               http://matplotlib.org/examples/color/colormaps_reference.html
                               Default is cm.rainbow. Good options are cm.coolwarm, cm.Blues, cm.jet


        Returns
        -------
        fig, ax                Matplotlib figure and axes instances for the plot.

        Notes
        -------
        Uses matplotlib.collections LineCollection, matplotlib.cm, matplotlib.pyplot and numpy arrays

        """
        if fig is None:
            fig = plt.figure()# define the figure and axes instances
        if ax is None:
            ax = plt.subplot(111)
            #ax = plt.axes()

        if zdata is None:
            z = np.ones(len(cells.gap_jun_i))

        elif zdata == 'random':
            z = np.random.random(len(cells.gap_jun_i))

        else:
            z = zdata

        if clrmap is None:
            clrmap = cm.bone_r  # default colormap

         # Make a line collection and add it to the plot.

        con_segs = cells.cell_centres[cells.gap_jun_i]

        connects = p.um*np.asarray(con_segs)

        coll = LineCollection(connects, array=z, cmap=clrmap, linewidths=4.0, zorder=0)
        coll.set_clim(vmin=0.0,vmax=1.0)
        coll.set_picker(pickable)
        ax.add_collection(coll)

        ax.axis('equal')

        # Add a colorbar for the Line Collection
        if zdata is not None and colorbar == 1:
            ax_cb = fig.colorbar(coll, ax=ax)
        else:
            ax_cb = None

        xmin = cells.xmin*p.um
        xmax = cells.xmax*p.um
        ymin = cells.ymin*p.um
        ymax = cells.ymax*p.um

        ax.axis([xmin,xmax,ymin,ymax])

        return fig, ax, ax_cb

def plotBoundCells(points_flat,bflags,cells, p, fig=None, ax=None):
        """
        Plot elements tagged on the boundary as red points.

        Parameters
        ----------
        points_flat          A flat array of points corresponding to the bflags data structure

        bflags          A nested array of boolean flags indicating boundary tagging

        Returns
        -------
        fig, ax         Matplotlib plotting objects

        Note
        ------
        This particular plot is extremely slow -- intended for cross-checking purposes only!

        """
        if fig is None:
            fig = plt.figure()# define the figure and axes instances
        if ax is None:
            ax = plt.subplot(111)
            #ax = plt.axes()

        points_flat = np.asarray(points_flat)
        bflags = np.asarray(bflags)

        bpoints = points_flat[bflags]

        ax.plot(p.um*points_flat[:,0],p.um*points_flat[:,1],'k.')

        ax.plot(p.um*bpoints[:,0],p.um*bpoints[:,1],'r.')

        # cell_edges_flat, _ , _= tb.flatten(cells.mem_edges)
        cell_edges_flat = p.um*cells.mem_edges_flat
        coll = LineCollection(cell_edges_flat,colors='k')
        coll.set_alpha(0.5)
        ax.add_collection(coll)

        ax.axis('equal')

        xmin = cells.xmin*p.um
        xmax = cells.xmax*p.um
        ymin = cells.ymin*p.um
        ymax = cells.ymax*p.um

        ax.axis([xmin,xmax,ymin,ymax])

        return fig, ax

def plotVects(cells, p, fig=None, ax=None):
        """
        This function plots all unit vectors in the tissue system as a cross-check.
        Normals to cell membranes are shown as red arrows.
        Tangents to cell membranes are black arrows.
        Tangents to ecm edges are shown as green arrows.
        Cell membrane edges are drawn as blue lines.

        To plot streamline and vector plots with data use the pyplot quiver and streamplot functions, respectively.

        """

        if fig is None:
            fig = plt.figure()# define the figure and axes instances

        if ax is None:
            ax = plt.subplot(111)
            #ax = plt.axes()

        s = p.um

        ax.quiver(s*cells.mem_vects_flat[:,0],s*cells.mem_vects_flat[:,1],s*cells.mem_vects_flat[:,4],
                  s*cells.mem_vects_flat[:,5],color='b',label='mem tang')
        ax.quiver(s*cells.mem_vects_flat[:,0],s*cells.mem_vects_flat[:,1],s*cells.mem_vects_flat[:,2],
                  s*cells.mem_vects_flat[:,3],color='g',label ='mem norm')
        # ax.quiver(s*cells.ecm_vects[:,0],s*cells.ecm_vects[:,1],s*cells.ecm_vects[:,2],s*cells.ecm_vects[:,3],color='r')

        # cell_edges_flat, _ , _= tb.flatten(cells.mem_edges)
        cell_edges_flat = cells.um*cells.mem_edges_flat
        coll = LineCollection(cell_edges_flat,colors='k')
        ax.add_collection(coll)

        ax.axis('equal')

        xmin = cells.xmin*p.um
        xmax = cells.xmax*p.um
        ymin = cells.ymin*p.um
        ymax = cells.ymax*p.um

        ax.axis([xmin,xmax,ymin,ymax])
        plt.legend()

        return fig, ax

def streamingCurrent(sim, cells,p,fig=None, ax=None, plot_Iecm = True, zdata = None,
    clrAutoscale = True, clrMin = None, clrMax = None, clrmap= cm.coolwarm,edgeOverlay = True,number_cells = False):

    if fig is None:
        fig = plt.figure()# define the figure and axes instances
    if ax is None:
        ax = plt.subplot(111)

    ax.axis('equal')

    xmin = cells.xmin*p.um
    xmax = cells.xmax*p.um
    ymin = cells.ymin*p.um
    ymax = cells.ymax*p.um

    ax.axis([xmin,xmax,ymin,ymax])

    if p.sim_ECM is False or plot_Iecm is False:

        Jmag_M = np.sqrt(sim.I_gj_x_time[-1]**2 + sim.I_gj_y_time[-1]**2) + 1e-30

        J_x = sim.I_gj_x_time[-1]/Jmag_M
        J_y = sim.I_gj_y_time[-1]/Jmag_M

        lw = (3.0*Jmag_M/Jmag_M.max()) + 0.5

        meshplot = plt.imshow(Jmag_M,origin='lower',extent=[xmin,xmax,ymin,ymax], cmap=clrmap)

        streamplot = ax.streamplot(cells.Xgrid*p.um,cells.Ygrid*p.um,J_x,J_y,density=p.stream_density,
            linewidth=lw,color='k', cmap=clrmap,arrowsize=1.5)

        ax.set_title('Final gap junction current density')

    elif plot_Iecm is True:
        Jmag_M = np.sqrt(sim.I_tot_x_time[-1]**2 + sim.I_tot_y_time[-1]**2) + 1e-30

        J_x = sim.I_tot_x_time[-1]/Jmag_M
        J_y = sim.I_tot_y_time[-1]/Jmag_M

        lw = (3.0*Jmag_M/Jmag_M.max()) + 0.5

        meshplot = plt.imshow(Jmag_M,origin='lower',extent=[xmin,xmax,ymin,ymax], cmap=clrmap)

        ax.streamplot(
            cells.Xgrid*p.um, cells.Ygrid*p.um, J_x, J_y,
            density=p.stream_density,
            linewidth=lw,
            color='k',
            cmap=clrmap,
            arrowsize=1.5,
        )

        ax.set_title('Final total currents')

    if clrAutoscale is True:
        ax_cb = fig.colorbar(meshplot,ax=ax)

    elif clrAutoscale is False:
        meshplot.set_clim(clrMin,clrMax)
        ax_cb = fig.colorbar(meshplot,ax=ax)

    # if p.showCells is True:
    #     # cell_edges_flat, _ , _= tb.flatten(cells.mem_edges)
    #     cell_edges_flat = cells.um*cells.mem_edges_flat
    #     coll = LineCollection(cell_edges_flat,colors='k')
    #     coll.set_alpha(0.2)
    #     ax.add_collection(coll)

    if number_cells is True:

        for i,cll in enumerate(cells.cell_centres):
            ax.text(p.um*cll[0],p.um*cll[1],i,ha='center',va='center')

    return fig,ax,ax_cb

def clusterPlot(p, dyna, cells, clrmap=cm.jet):

    fig = plt.figure()
    ax = plt.subplot(111)

    # profile_names = list(p.profiles.keys())

    col_dic = {}

    cb_ticks = []
    cb_tick_labels = []

    base_points = np.multiply(cells.cell_verts, p.um)

    z = np.zeros(len(base_points))
    z[:] = 0

    cb_ticks.append(0)
    cb_tick_labels.append(p.default_tissue_name)

    col_dic['base'] = PolyCollection(
        base_points, array=z, cmap=clrmap, edgecolors='none')
    ax.add_collection(col_dic['base'])

    if len(dyna.tissue_profile_names):
        for i, name in enumerate(dyna.tissue_profile_names):
            cell_inds = dyna.cell_target_inds[name]
            points = np.multiply(cells.cell_verts[cell_inds], p.um)

            z = np.zeros(len(points))
            z[:] = i + 1

            col_dic[name] = PolyCollection(
                points, array=z, cmap=clrmap, edgecolors='none')
            col_dic[name].set_clim(0, len(dyna.tissue_profile_names))

            # col_dic[name].set_alpha(0.8)
            col_dic[name].set_zorder(p.profiles[name]['z order'])
            ax.add_collection(col_dic[name])

            # Add this profile name to the colour legend.
            cb_ticks.append(i+1)
            cb_tick_labels.append(name)

    #FIXME: Refactor into a new ActionCut.plot() method. Beans and fragrance!
    if p.plot_cutlines and p.scheduled_options['cuts'] is not None:
        # For each profile cutting a subset of the cell population...
        for cut_profile_name in p.scheduled_options['cuts'].profile_names:
            cut_profile = p.profiles[cut_profile_name]

            # Indices of all cells cut by this profile.
            cut_cell_indices = cut_profile.picker.get_cell_indices(
                cells, p, ignoreECM=True)

            points = np.multiply(cells.cell_verts[cut_cell_indices], p.um)
            col_dic[cut_profile_name] = PolyCollection(
                points, color='k', cmap=clrmap, edgecolors='none')

            # col_dic[name].set_clim(0,len(dyna.tissue_profile_names) + len(names))
            # col_dic[name].set_alpha(0.8)

            col_dic[cut_profile_name].set_zorder(cut_profile.z_order)
            ax.add_collection(col_dic[cut_profile_name])

            #FIXME: Interestingly, this doesn't appear to do anything. I have
            #no idea why. The matpotlib is weak with me. Legends and old elves!

            # Add this profile name to the colour legend.
            cb_tick_next = len(cb_ticks)
            cb_ticks.append(cb_tick_next)
            cb_tick_labels.append(cut_profile_name)

    ax_cb = None
    if len(dyna.tissue_profile_names):
        ax_cb = fig.colorbar(
            col_dic[dyna.tissue_profile_names[0]], ax=ax, ticks=cb_ticks)
        ax_cb.ax.set_yticklabels(cb_tick_labels)

    if p.enumerate_cells is True:
        for i, cll in enumerate(cells.cell_centres):
            ax.text(
                p.um*cll[0], p.um*cll[1], i,
                ha='center', va='center', zorder=20)

    ax.set_xlabel('Spatial Distance [um]')
    ax.set_ylabel('Spatial Distance [um]')
    ax.set_title('Cell Cluster')

    ax.axis('equal')

    xmin = cells.xmin*p.um
    xmax = cells.xmax*p.um
    ymin = cells.ymin*p.um
    ymax = cells.ymax*p.um

    ax.axis([xmin,xmax,ymin,ymax])

    return fig, ax, ax_cb

def exportData(cells,sim,p):


    if p.plot_type is 'sim':

        results_path = p.sim_results

    elif p.plot_type is 'init':
        results_path = p.init_results

    os.makedirs(results_path, exist_ok=True)
    savedData = os.path.join(results_path, 'ExportedData.csv')
    savedData_FFT = os.path.join(results_path, 'ExportedData_FFT.csv')

    cc_cell = []

    dd_cell = []

    ci = p.plot_cell  # index of cell to get time data for

    # create the header, first entry will be time:
    headr = 'time_s'
    t = np.asarray(sim.time)

    #-----------Vmem----------------------------------------------

    # next entry will be Vm:
    headr = headr + ',' + 'Vmem_mV'

    if p.sim_ECM is False:
        vm = [arr[ci]*1000 for arr in sim.vm_time]

    else:
        vm = []
        for vm_at_mem in sim.vm_time:
            vm_t = 1000*cell_ave(cells,vm_at_mem)[ci]
            vm.append(vm_t)

    vm = np.asarray(vm)

    # golman Vmem------------------------------------------------------

    if p.GHK_calc is True:

        vm_goldman = [x[p.plot_cell]*1000 for x in sim.vm_GHK_time]

    else:
        vm_goldman = np.zeros(len(sim.time))

    vm_goldman = np.asarray(vm_goldman)
    # next entry will be Vm:
    headr = headr + ',' + 'Goldman_Vmem_mV'

    # ------Na K pump rate-----------------------------------------------

    # Na-K-pump rate:
    if p.sim_ECM is False:
        pump_rate = [pump_array[p.plot_cell] for pump_array in sim.rate_NaKATP_time]

    else:
        pump_rate = [pump_array[cells.cell_to_mems[p.plot_cell][0]] for pump_array in sim.rate_NaKATP_time]

    pump_rate = np.asarray(pump_rate)
    headr = headr + ',' + 'NaK-ATPase_Rate_mol/m2s'

    #----------cell concentrations---------------------------------------

    # create the header starting with cell concentrations
    for i in range(0,len(sim.ionlabel)):
        label = sim.ionlabel[i]
        headr = headr + ',' + 'cell_' + label + '_mmol/L'
        cc_m = [arr[i][ci] for arr in sim.cc_time]
        cc_m = np.asarray(cc_m)
        cc_cell.append(cc_m)


    cc_cell = np.asarray(cc_cell)

    #-----optional concentrations--------------------------------------

    if p.scheduled_options['IP3'] != 0 or p.Ca_dyn is True:

        IP3_time = [arr[ci] for arr in sim.cIP3_time]
        IP3_time = np.asarray(IP3_time)
        headr = headr + ',' + 'cell_cIP3_mmol/L'

    else:
        IP3_time = np.zeros(len(sim.time))
        headr = headr + ',' + 'cell_cIP3_mmol/L'

    if p.voltage_dye ==1:
        dye_time = [arr[ci] for arr in sim.cDye_time]
        dye_time = np.asarray(dye_time)
        headr = headr + ',' + 'cell_dye_mmol/L'
    else:
        dye_time = np.zeros(len(sim.time))
        headr = headr + ',' + 'cell_dye_mmol/L'

    if p.Ca_dyn == 1 and p.ions_dict['Ca']==1:
        Ca_er = [arr[0][ci] for arr in sim.cc_er_time]
        Ca_er = np.asarray(Ca_er)
        headr = headr + ',' + 'ER_Ca2+_mmol/L'
    else:
        Ca_er = np.zeros(len(sim.time))
        headr = headr + ',' + 'CaER_mmol/L'

    #----------membrane permeabilities---------------------------------------

    # create the header starting with membrane permeabilities
    if p.sim_ECM is False:
        for i in range(0,len(sim.ionlabel)):
            label = sim.ionlabel[i]
            headr = headr + ',' + 'Dm_' + label + '_m2/s'
            dd_m = [arr[i][ci] for arr in sim.dd_time]
            dd_m = np.asarray(dd_m)
            dd_cell.append(dd_m)

    else:
        for i in range(0,len(sim.ionlabel)):
            label = sim.ionlabel[i]
            headr = headr + ',' + 'Dm_' + label + '_m2/s'
            dd_m = [arr[i][cells.cell_to_mems[ci][0]] for arr in sim.dd_time]
            dd_m = np.asarray(dd_m)
            dd_cell.append(dd_m)


    dd_cell = np.asarray(dd_cell)

    #----Transmembrane currents--------------------------
    if p.sim_ECM is False:
        Imem = [memArray[p.plot_cell] for memArray in sim.I_mem_time]
    else:
        Imem = [memArray[cells.cell_to_mems[p.plot_cell][0]] for memArray in sim.I_mem_time]

    headr = headr + ',' + 'I_A/m2'

    #----Hydrostatic pressure--------
    p_hydro = [arr[p.plot_cell] for arr in sim.P_cells_time]

    headr = headr + ',' + 'HydroP_Pa'

    # ---Osmotic pressure-----------
    if p.deform_osmo is True:

        p_osmo = [arr[p.plot_cell] for arr in sim.osmo_P_delta_time]

    else:
        p_osmo = np.zeros(len(sim.time))

    headr = headr + ',' + 'OsmoP_Pa'

    # electrostatic pressure -----------------------------------
    if p.deform_electro is True:

        p_electro = [arr[p.plot_cell] for arr in sim.P_electro_time]

    else:
        p_electro = np.zeros(len(sim.time))

    headr = headr + ',' + 'ElectroP_Pa'

    # total deformation ---------------------------------------
    if p.deformation is True and sim.run_sim is True:

        # extract time-series deformation data for the plot cell:
        dx = np.asarray([arr[p.plot_cell] for arr in sim.dx_cell_time])
        dy = np.asarray([arr[p.plot_cell] for arr in sim.dy_cell_time])

        # get the total magnitude:
        disp = p.um*np.sqrt(dx**2 + dy**2)

    else:
        disp = np.zeros(len(sim.time))

    headr = headr + ',' + 'Displacement_um'


    # FFT of voltage :
    sample_size = len(sim.time)
    sample_spacing = sim.time[1] - sim.time[0]

    cell_data = (1/sample_size)*(vm - np.mean(vm))

    f_axis = np.fft.rfftfreq(sample_size, d=sample_spacing)
    fft_data_o = np.fft.rfft(cell_data)
    fft_data = np.sqrt(np.real(fft_data_o)**2 + np.imag(fft_data_o)**2)

    dataM = np.column_stack((t,vm,vm_goldman,pump_rate,cc_cell.T,IP3_time,dye_time,Ca_er,dd_cell.T,Imem,
                             p_hydro,p_osmo,p_electro,disp))

    headr2 = 'frequency_Hz'
    headr2 = headr2 + ',' + 'FFT_Vmem'

    dataFFT = np.column_stack((f_axis,fft_data))

    np.savetxt(savedData,dataM,delimiter = ',',header = headr)
    np.savetxt(savedData_FFT,dataFFT,delimiter = ',',header = headr2)

def export2dData(simdata,cells,p):

    if p.plot_type == 'sim':
        results_path = p.sim_results

    elif p.plot_type == 'init':
        results_path = p.sim_results

    os.makedirs(results_path, exist_ok=True)
    savedData_2d = os.path.join(results_path, 'Exported2DData.csv')


    dataM = simdata
    np.savetxt(savedData_2d,dataM,delimiter=',')

def I_overlay(sim,cells,p,ax,clrmap,plotIecm = False, time=-1):

    if p.sim_ECM is False or plotIecm is False:

        Jmag_M = np.sqrt(sim.I_gj_x_time[-1]**2 + sim.I_gj_y_time[-1]**2) + 1e-30

        J_x = sim.I_gj_x_time[-1]/Jmag_M
        J_y = sim.I_gj_y_time[-1]/Jmag_M

        lw = (3.0*Jmag_M/Jmag_M.max()) + 0.5

        ax.streamplot(cells.Xgrid*p.um,cells.Ygrid*p.um,J_x,J_y,density=p.stream_density,linewidth=lw,color='k',cmap=clrmap,arrowsize=1.5)

        ax.set_title('(gap junction current overlay)')

    elif plotIecm is True:

        Jmag_M = np.sqrt(sim.I_tot_x_time[-1]**2 + sim.I_tot_y_time[-1]**2) + 1e-30

        J_x = sim.I_tot_x_time[-1]/Jmag_M
        J_y = sim.I_tot_y_time[-1]/Jmag_M

        lw = (3.0*Jmag_M/Jmag_M.max()) + 0.5

        ax.streamplot(cells.Xgrid*p.um,cells.Ygrid*p.um,J_x,J_y,density=p.stream_density,linewidth=lw,color='k',cmap=clrmap,arrowsize=1.5)

        ax.set_title('(total current overlay)')

def I_overlay_setup(sim, ax, cells, p):

    if p.sim_ECM is False or p.IecmPlot is False:

        Jmag_M = np.sqrt(sim.I_gj_x_time[-1]**2 + sim.I_gj_x_time[-1]**2) + 1e-30

        J_x = sim.I_gj_x_time[-1]/Jmag_M
        J_y = sim.I_gj_x_time[-1]/Jmag_M

        lw = (3.0*Jmag_M/Jmag_M.max()) + 0.5

        streams = ax.streamplot(cells.Xgrid*p.um,cells.Ygrid*p.um,J_x,J_y,density=p.stream_density,linewidth=lw,
            color=p.vcolor, arrowsize=1.5)

        tit_extra = 'Gap junction current'

    elif p.IecmPlot is True:

        Jmag_M = np.sqrt(sim.I_tot_x_time[0]**2 + sim.I_tot_y_time[0]**2) + 1e-30

        J_x = sim.I_tot_x_time[0]/Jmag_M
        J_y = sim.I_tot_y_time[0]/Jmag_M

        lw = (3.0*Jmag_M/Jmag_M.max()) + 0.5

        streams = ax.streamplot(cells.Xgrid*p.um,cells.Ygrid*p.um,J_x,J_y,density=p.stream_density,linewidth=lw,color='k',
            arrowsize=1.5)

        tit_extra = 'Total current overlay'

    return streams, ax, tit_extra

def I_overlay_update(i,sim,streams,ax,cells,p):

    if p.sim_ECM is False or p.IecmPlot is False:
        Jmag_M = np.sqrt(sim.I_gj_x_time[i]**2 + sim.I_gj_y_time[i]**2) + 1e-30

        J_x = sim.I_gj_x_time[i]/Jmag_M
        J_y = sim.I_gj_y_time[i]/Jmag_M

    elif p.IecmPlot is True:

        Jmag_M = np.sqrt(
            sim.I_tot_x_time[i]**2 +
            sim.I_tot_y_time[i]**2) + 1e-30

        J_x = sim.I_tot_x_time[i]/Jmag_M
        J_y = sim.I_tot_y_time[i]/Jmag_M

    lw = (3.0*Jmag_M/Jmag_M.max()) + 0.5

    streams.lines.remove()
    ax.patches = []

    streams = ax.streamplot(
        cells.Xgrid*1e6,
        cells.Ygrid*1e6, J_x, J_y,
        density=p.stream_density,
        linewidth=lw,
        color='k',
        arrowsize=1.5)

    return streams, ax

def cell_ave(cells,vm_at_mem):

    """
    Averages Vmem over membrane domains to return a mean value for each cell

    Parameters
    ----------
    cells               An instance of the Cells module cells object
    vm_at_mem           Vmem at individual membrane domains


    Returns
    --------
    v_cell              Cell Vm averaged over the whole cell

    """

    v_cell = []

    for i in cells.cell_i:
        cellinds = (cells.mem_to_cells == i).nonzero()
        v_cell_array = vm_at_mem[cellinds]
        v_cell_ave = np.mean(v_cell_array)
        v_cell.append(v_cell_ave)

    v_cell = np.asarray(v_cell)

    return v_cell

# utility functions------------------------------------------------------
def cell_quiver(datax,datay,cells,p):
    pass

def env_quiver(datax,datay,ax,cells,p):
    F_mag = np.sqrt(datax**2 + datay**2)

    if F_mag.all() != 0.0:
        Fx = datax/F_mag
        Fy = datay/F_mag

    ax.quiver(p.um*cells.xypts[:,0], p.um*cells.xypts[:,1], Fx.ravel(),
        Fy.ravel())

def cell_stream(datax,datay,cells,p):
    pass

def env_stream(datax,datay,cells,p):
    pass

def cell_mesh(data,ax,cells,p,clrmap):

    data_grid = np.zeros(len(cells.voronoi_centres))
    data_grid[cells.cell_to_grid] = data

    msh = ax.tripcolor(p.um*cells.voronoi_centres[:,0],p.um*cells.voronoi_centres[:,1],data_grid,
        shading='gouraud',cmap=clrmap)

    return msh, ax

def env_mesh(data,ax,cells,p,clrmap):

        if p.plotMask is True:
            data = ma.masked_array(data, np.logical_not(cells.maskM))

        mesh_plot = ax.imshow(data,origin='lower',
                    extent=[p.um*cells.xmin,p.um*cells.xmax,p.um*cells.ymin,p.um*cells.ymax],cmap=clrmap)

        if p.showCells is True:
            cell_edges_flat = cells.um*cells.mem_edges_flat
            coll = LineCollection(cell_edges_flat,colors='k')
            coll.set_alpha(0.5)
            ax.add_collection(coll)

        return mesh_plot, ax

def cell_mosaic(data,ax,cells,p,clrmap):

        # define a polygon collection based on individual cell polygons
    points = np.multiply(cells.cell_verts, p.um)
    collection =  PolyCollection(points, cmap=clrmap, edgecolors='none')
    collection.set_array(data)
    ax.add_collection(collection)

    return collection, ax

def show_plot(params: 'Parameters', *args, **kwargs) -> None:
    '''
    Display the current plot if the passed configuration requests plots to be
    displayed or noop otherwise.

    All passed arguments following the passed configuration will be passed as
    is to the matplotlib.pyplot.show() function.
    '''
    assert types.is_parameters(params), (
        types.assert_not_parameters(params))

    try:
        if params.turn_all_plots_off is False:
            plt.show()
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

def set_up_filesave(ani_obj,p):

    if p.plot_type == 'sim':
        # Make the BETSE-specific cache directory if not found.
        images_path = os.path.join(p.sim_results, ani_obj.saveFolder)

    elif p.plot_type is 'init':
        images_path = os.path.join(p.init_results, ani_obj.saveFolder)

    betse_cache_dir = os.path.expanduser(images_path)
    os.makedirs(betse_cache_dir, exist_ok=True)
    ani_obj.savedAni = os.path.join(betse_cache_dir, ani_obj.saveFile)
    ani_obj.ani_repeat = False

    return ani_obj














