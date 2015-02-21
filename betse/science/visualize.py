#!/usr/bin/env python3
# Copyright 2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

# FIXME create a plotting method that plots individual cell data
# FIXME create a plotting method for ecm data
# FIXME work on the vector plotting -- perhaps interpolating data to grid?
# FIXME saving animations doesn't work
# FIXME do an animation for smoothed 'vert' data and gap junctions (with fluxes)

import warnings
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection, PolyCollection
import matplotlib.cm as cm
from betse.science import toolbox as tb
from betse.science.parameters import params as p
from matplotlib import animation
import os, os.path

class AnimateCellData(object):
    """
    Animate color data on a plot of cells.

    """

    def __init__(self,cells,zdata_t,time,p,tit=' ',cbtit = ' ', colormap=p.default_cm, save=False,ani_repeat=False,
        number_cells = False, saveFolder = '/animation', saveFile = 'sim_'):

        self.zdata_t = zdata_t
        self.colormap = colormap
        self.time = time
        self.save = save

        self.cbtit = cbtit

        if self.save == True:
            # Make the BETSE-specific cache directory if not found.
            images_path = p.sim_results + saveFolder
            betse_cache_dir = os.path.expanduser(images_path)
            os.makedirs(betse_cache_dir, exist_ok=True)
            self.savedAni = os.path.join(betse_cache_dir, saveFile)

        # define a polygon collection based on individual cell polygons
        self.points = np.multiply(cells.cell_verts, p.um)
        self.collection =  PolyCollection(self.points, cmap=self.colormap, edgecolors='none')
        self.collection.set_array(self.zdata_t[0])

        # set range of the colormap
        self.cmean = np.mean(self.zdata_t)
        self.cmin = np.min(self.zdata_t)
        self.cmax = np.max(self.zdata_t)

        self.fig = plt.figure()       # define figure
        self.ax = plt.subplot(111)    # define axes
        self.collection.set_clim(self.cmin,self.cmax)

        self.cb = self.fig.colorbar(self.collection)   # define colorbar for figure
        self.cb.set_label(self.cbtit)

        self.ax.add_collection(self.collection)

        self.tit = tit

        if number_cells == True:
            for i,cll in enumerate(cells.cell_centres):
                self.ax.text(p.um*cll[0],p.um*cll[1],i)

        self.ax.set_xlabel('Spatial x [um]')
        self.ax.set_ylabel('Spatial y [um')
        self.ax.set_title(self.tit)

        self.ax.autoscale_view()

        self.frames = len(self.zdata_t)

        ani = animation.FuncAnimation(self.fig, self.aniFunc,
            frames=self.frames, interval=66, repeat=ani_repeat)

        plt.show()


    def aniFunc(self,i):

        zz = self.zdata_t[i]

        self.collection.set_array(zz)
        titani = self.tit + ' (simulation time' + ' ' + str(round(self.time[i],3)) + ' ' + ' s)'
        self.ax.set_title(titani)

        if self.save == True:
            self.fig.canvas.draw()
            savename = self.savedAni + str(i)
            plt.savefig(savename,dpi=96,format='png')

class AnimateCellData_smoothed(object):

    def __init__(self,cells,zdata_t,time,p,tit=' ',cbtit = ' ', colormap=p.default_cm, save=False,ani_repeat=False,
        number_cells = False, saveFolder = '/animation', saveFile = 'sim_'):

        self.zdata_t = zdata_t
        self.colormap = colormap
        self.time = time
        self.save = save

        self.cbtit = cbtit

        if self.save == True:
            # Make the BETSE-specific cache directory if not found.
            images_path = p.sim_results + saveFolder
            betse_cache_dir = os.path.expanduser(images_path)
            os.makedirs(betse_cache_dir, exist_ok=True)
            self.savedAni = os.path.join(betse_cache_dir, saveFile)



        # set range of the colormap
        self.cmean = np.mean(self.zdata_t)
        self.cmin = np.min(self.zdata_t)
        self.cmax = np.max(self.zdata_t)

        self.fig = plt.figure()       # define figure
        self.ax = plt.subplot(111)    # define axes

         # Next add a triplot with interpolated and animated voltage data
        self.triplt = self.ax.tripcolor(p.um*cells.cell_centres[:, 0], p.um*cells.cell_centres[:, 1],
            self.zdata_t[0],shading='gouraud', cmap=self.colormap)

        self.triplt.set_clim(self.cmin,self.cmax)

        self.cb = self.fig.colorbar(self.triplt)   # define colorbar for figure
        self.cb.set_label(self.cbtit)

        self.tit = tit

        if number_cells == True:
            for i,cll in enumerate(cells.cell_centres):
                self.ax.text(p.um*cll[0],p.um*cll[1],i)

        self.ax.set_xlabel('Spatial x [um]')
        self.ax.set_ylabel('Spatial y [um')
        self.ax.set_title(self.tit)

        self.ax.autoscale_view()

        self.frames = len(self.zdata_t)

        ani = animation.FuncAnimation(self.fig, self.aniFunc,
            frames=self.frames, interval=66, repeat=ani_repeat)

        plt.show()


    def aniFunc(self,i):

        zz = self.zdata_t[i]

        self.triplt.set_array(zz)
        titani = self.tit + ' (simulation time' + ' ' + str(round(self.time[i],3)) + ' ' + ' s)'
        self.ax.set_title(titani)

        if self.save == True:
            self.fig.canvas.draw()
            savename = self.savedAni + str(i)
            plt.savefig(savename,dpi=96,format='png')

class AnimateGJData(object):
    """
    Animate the gap junction open state as a function of time.
    """

    def __init__(self,cells,sim,p,tit=' ', colormap=p.default_cm, save=False,saveFolder = '/animation',
        saveFile = 'sim_',ani_repeat=False,number_cells=False):

        #self.zdata_t = sim.gjopen_time  # data array for gap junction coloring
        gjI_t = np.asarray(sim.Igj_time)
        normI = np.max(gjI_t)
        self.zdata_t = gjI_t/normI

        self.vdata_t = np.multiply(sim.vm_time,1000)   # data array for cell coloring
        self.colormap = colormap
        self.time = sim.time

        self.gjI_t = np.sign(sim.Igj_time)
        self.gjvects = cells.gj_vects

        self.fig = plt.figure()       # define figure
        self.ax = plt.subplot(111)    # define axes

        self.tit = tit

        self.save = save

        if self.save == True:
            # Make the BETSE-specific cache directory if not found.
            images_path = p.sim_results + saveFolder
            betse_cache_dir = os.path.expanduser(images_path)
            os.makedirs(betse_cache_dir, exist_ok=True)
            self.savedAni = os.path.join(betse_cache_dir, saveFile)

        con_segs = cells.cell_centres[cells.gap_jun_i]
        connects = p.um*np.asarray(con_segs)
        self.collection = LineCollection(connects, array=self.zdata_t[0], cmap= p.gj_cm, linewidths=3.0, zorder=5)
        self.collection.set_clim(0.0,1.0)
        self.ax.add_collection(self.collection)

        # Next add a collection of cell polygons, with animated voltage data
        points = np.multiply(cells.cell_verts, p.um)
        self.coll2 =  PolyCollection(points, array=self.vdata_t[0], edgecolors='none', cmap=self.colormap)
        self.coll2.set_alpha(1.0)
         # set range of the colormap
        self.cmean = np.mean(self.vdata_t)
        self.cmin = np.min(self.vdata_t)
        self.cmax = np.max(self.vdata_t)
        self.coll2.set_clim(self.cmin,self.cmax)
        self.cb = self.fig.colorbar(self.coll2)   # define colorbar for figure
        self.ax.add_collection(self.coll2)

        # Next add in gap junction current direction
        vx = np.multiply(self.gjI_t[0],self.gjvects[:,2])
        vy = np.multiply(self.gjI_t[0],self.gjvects[:,3])

        self.Qplot = self.ax.quiver(p.um*self.gjvects[:,0],p.um*self.gjvects[:,1],
            vx,vy,self.zdata_t[0],zorder=10, cmap=p.gj_cm,clim=[0,1])

        if number_cells == True:
            for i,cll in enumerate(cells.cell_centres):
                self.ax.text(p.um*cll[0],p.um*cll[1],i)

        self.cb.set_label('Voltage [mV]')
        self.ax.set_xlabel('Spatial x [um]')
        self.ax.set_ylabel('Spatial y [um')
        self.ax.set_title(self.tit)

        self.ax.autoscale_view()

        self.frames = len(self.zdata_t)

        ani = animation.FuncAnimation(self.fig, self.aniFunc,
               frames=self.frames, interval=100, repeat=ani_repeat)

        plt.show()


    def aniFunc(self,i):

        zz = self.zdata_t[i]
        zv = self.vdata_t[i]

        vx = np.multiply(self.gjI_t[i],self.gjvects[:,2])
        vy = np.multiply(self.gjI_t[i],self.gjvects[:,3])

        self.collection.set_array(zz)

        self.coll2.set_array(zv)
        self.Qplot.set_UVC(vx,vy,zz)

        titani = self.tit + ' ' + '(simulation time' + ' ' + str(round(self.time[i],3)) + ' ' + 's)'
        self.ax.set_title(titani)

        if self.save == True:
            self.fig.canvas.draw()
            savename = self.savedAni + str(i)
            plt.savefig(savename,dpi=96,format='png')

class AnimateGJData_smoothed(object):

    def __init__(self,cells,sim,p,tit=' ', colormap=p.default_cm, save=False,saveFolder = '/animation',
        saveFile = 'sim_',ani_repeat=False,number_cells=False):

        # self.zdata_t = sim.gjopen_time  # data array for gap junction coloring
        gjI_t = np.asarray(sim.Igj_time)
        normI = np.max(gjI_t)
        self.zdata_t = gjI_t/normI

        self.vdata_t = np.multiply(sim.vm_time,1000)   # data array for cell coloring
        self.colormap = colormap
        self.time = sim.time

        self.gjI_t = np.sign(sim.Igj_time)
        self.gjvects = cells.gj_vects

        self.fig = plt.figure()       # define figure
        self.ax = plt.subplot(111)    # define axes

        self.tit = tit

        self.save = save

        if self.save == True:
            # Make the BETSE-specific cache directory if not found.
            images_path = p.sim_results + saveFolder
            betse_cache_dir = os.path.expanduser(images_path)
            os.makedirs(betse_cache_dir, exist_ok=True)
            self.savedAni = os.path.join(betse_cache_dir, saveFile)

        con_segs = cells.cell_centres[cells.gap_jun_i]
        connects = p.um*np.asarray(con_segs)
        self.collection = LineCollection(connects, array=self.zdata_t[0], cmap=p.gj_cm, linewidths=3.0, zorder=5)
        self.collection.set_clim(0.0,1.0)
        self.ax.add_collection(self.collection)

        # Next add a triplot with interpolated and animated voltage data
        self.triplt = self.ax.tripcolor(p.um*cells.cell_centres[:, 0], p.um*cells.cell_centres[:, 1],
            self.vdata_t[0],shading='gouraud', cmap=self.colormap)
         # set range of the colormap
        self.cmean = np.mean(self.vdata_t)
        self.cmin = np.min(self.vdata_t)
        self.cmax = np.max(self.vdata_t)
        self.triplt.set_clim(self.cmin,self.cmax)
        self.cb = self.fig.colorbar(self.triplt)   # define colorbar for figure

        # Next add in gap junction current direction
        vx = np.multiply(self.gjI_t[0],self.gjvects[:,2])
        vy = np.multiply(self.gjI_t[0],self.gjvects[:,3])

        self.Qplot = self.ax.quiver(p.um*self.gjvects[:,0],p.um*self.gjvects[:,1],
            vx,vy,self.zdata_t[0],zorder=10, cmap=p.gj_cm,clim=[0,1])

        if number_cells == True:
            for i,cll in enumerate(cells.cell_centres):
                self.ax.text(p.um*cll[0],p.um*cll[1],i)

        self.cb.set_label('Voltage [mV]')
        self.ax.set_xlabel('Spatial x [um]')
        self.ax.set_ylabel('Spatial y [um')
        self.ax.set_title(self.tit)

        self.ax.autoscale_view()

        self.frames = len(self.zdata_t)

        ani = animation.FuncAnimation(self.fig, self.aniFunc,
               frames=self.frames, interval=100, repeat=ani_repeat)

        plt.show()

    def aniFunc(self,i):

        zz = self.zdata_t[i]
        zv = self.vdata_t[i]

        vx = np.multiply(self.gjI_t[i],self.gjvects[:,2])
        vy = np.multiply(self.gjI_t[i],self.gjvects[:,3])

        self.collection.set_array(zz)
        self.triplt.set_array(zv)
        self.Qplot.set_UVC(vx,vy,zz)

        titani = self.tit + ' ' + '(simulation time' + ' ' + str(round(self.time[i],3)) + ' ' + 's)'
        self.ax.set_title(titani)

        if self.save == True:
            self.fig.canvas.draw()
            savename = self.savedAni + str(i)
            plt.savefig(savename,dpi=96,format='png')


def plotSingleCellVData(simdata_time,simtime,celli,fig=None,ax=None, lncolor='b'):

    tvect_data=[x[celli]*1000 for x in simdata_time]

    if fig==None:
        fig = plt.figure()# define the figure and axes instances
    if ax == None:
        ax = plt.subplot(111)
        #ax = plt.axes()

    ax.plot(simtime, tvect_data,lncolor)
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Voltage [mV]')
    #ax.axis('equal')

    return fig, ax

def plotSingleCellCData(simdata_time,simtime,ioni,celli,fig=None,ax=None,lncolor='b',ionname='ion'):

    # ccIon = [arr[ion] for arr in simdata_time]  # get all cells at all times at one ion
    # ccIon_cell = [x[0] for x in ccIon]  # get one cell at all times at one ion

    ccIon_cell = [arr[ioni][celli] for arr in simdata_time]

    if fig==None:
        fig = plt.figure()# define the figure and axes instances
    if ax == None:
        ax = plt.subplot(111)
        #ax = plt.axes()

    lab = ionname

    xmin = simtime[0]
    xmax = simtime[-1]
    ymin = np.min(ccIon_cell)
    ymax = np.max(ccIon_cell)

    ax.plot(simtime, ccIon_cell,lncolor,label=lab)
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Concentration [mol/m3]')
    # ax.axis([xmin,xmax,ymin,ymax])

    return fig, ax

def plotSingleCellData(simtime,simdata_time,celli,fig=None,ax=None,lncolor='b',lab='Data'):

    data_cell = [arr[celli] for arr in simdata_time]

    if fig==None:
        fig = plt.figure()# define the figure and axes instances
    if ax == None:
        ax = plt.subplot(111)

    xmin = simtime[0]
    xmax = simtime[-1]
    ymin = np.min(data_cell)
    ymax = np.max(data_cell)

    ax.plot(simtime, data_cell,lncolor,label=lab)
    ax.set_xlabel('Time [s]')
    ax.set_ylabel(lab)
    # ax.axis([xmin,xmax,ymin,ymax])


    return fig, ax

def plotPolyData(cells, p, fig=None, ax=None, zdata = None,clrmap = None, number_cells=False):
        """
        Assigns color-data to each polygon in a cell cluster diagram and returns a plot instance (fig, axes)

        Parameters
        ----------
        vor_verts              Nested list of [x,y] points defining each polygon. May be ecm_verts or
                               cell_verts

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
        Uses matplotlib.collections PolyCollection, matplotlib.cm, matplotlib.pyplot and numpy arrays
        Computationally slow -- not recommended for large collectives (500 x 500 um max)
        """
        if fig==None:
            fig = plt.figure()# define the figure and axes instances
        if ax == None:
            ax = plt.subplot(111)
            #ax = plt.axes()

        if zdata == None:  # if user doesn't supply data
            z = np.ones(len(cells.cell_verts)) # create flat data for plotting

        elif zdata == 'random':  # if user doesn't supply data
            z = np.random.random(len(cells.cell_verts)) # create some random data for plotting

        else:
            z = zdata

        # Make the polygon collection and add it to the plot.
        if clrmap == None:
            clrmap = p.default_cm

        points = np.multiply(cells.cell_verts, p.um)

        coll = PolyCollection(points, array=z, cmap=clrmap, edgecolors='none')
        # coll.set_picker(True)
        ax.add_collection(coll)
        ax.axis('equal')

        # Add a colorbar for the PolyCollection
        std = np.std(zdata,axis=0)
        mean = np.mean(zdata,axis=0)

        if zdata != None and std > 1e-13:
            ax_cb = fig.colorbar(coll,ax=ax)

        elif std<1e-13:

            if mean < 0:
                coll.set_clim(mean,0)
            if mean > 0:
                coll.set_clim(0,mean)

            ax_cb = fig.colorbar(coll,ax=ax)

        elif zdata == None:
            coll.set_clim(0,1)
            ax_cb = fig.colorbar(coll,ax=ax)

        if number_cells == True:
            for i,cll in enumerate(cells.cell_centres):
                ax.text(p.um*cll[0],p.um*cll[1],i,ha='center',va='center')

        ax.autoscale_view(tight=True)

        return fig,ax,ax_cb

def plotCellData(cells, p, fig=None, ax=None, zdata=None,clrmap=None,edgeOverlay = None,pointOverlay=None):
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
                               Default is cm.rainbow. Good options are cm.coolwarm, cm.Blues, cm.jet

        edgeOverlay             This option allows the user to specify whether or not they want cell edges overlayed.
                                Default is False, set to True to use.

        pointOverlay            This option allows user to specify whether or not they want cell_centre points plotted
                                Default is False, set to True to use.


        Returns
        -------
        fig, ax                Matplotlib figure and axes instances for the plot.

        Notes
        -------
        Uses matplotlib.pyplot and numpy arrays
        With edgeOverlay and pointOverlay == None, this is computationally fast and *is* recommended for plotting data
        on large collectives.


        """

        if fig==None:
            fig = plt.figure()# define the figure and axes instances
        if ax == None:
            ax = plt.subplot(111)
            #ax = plt.axes()

        if zdata == None:  # if user doesn't supply data
            z = np.ones(len(cells.cell_centres)) # create flat data for plotting

        elif zdata == 'random':  # if user doesn't supply data
            z = np.random.random(len(cells.cell_centres)) # create some random data for plotting

        else:
            z = zdata

        if clrmap == None:
            clrmap = p.default_cm


        triplt = ax.tripcolor(p.um*cells.cell_centres[:, 0], p.um*cells.cell_centres[:, 1],
            z,shading='gouraud', cmap=clrmap)

        ax.axis('equal')

        # Add a colorbar for the z-data
        # if zdata != None:
        #     ax_cb = fig.colorbar(triplt, ax=ax)

         # Add a colorbar for the triplot:
        std = np.std(zdata,axis=0)
        mean = np.mean(zdata,axis=0)

        if zdata != None and std > 1e-13:
            ax_cb = fig.colorbar(triplt,ax=ax)

        elif std<1e-13:

            if mean < 0:
                triplt.set_clim(mean,0)
            if mean > 0:
                triplt.set_clim(0,mean)

            ax_cb = fig.colorbar(triplt,ax=ax)

        if pointOverlay == True:
            ax.scatter(p.um*cells.cell_centres[:,0],p.um*cells.cell_centres[:,1], c=z,cmap=clrmap)

        if edgeOverlay == True:
            cell_edges_flat, _ , _= tb.flatten(cells.mem_edges)
            cell_edges_flat = cells.um*np.asarray(cell_edges_flat)
            coll = LineCollection(cell_edges_flat,colors='k')
            coll.set_alpha(0.5)
            ax.add_collection(coll)


        ax.autoscale_view(tight=True)


        return fig, ax, ax_cb

def plotMemData(cells, fig= None, ax = None, zdata=None,clrmap=None):
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

        if fig==None:
            fig = plt.figure()# define the figure and axes instances
        if ax == None:
            ax = plt.subplot(111)
            #ax = plt.axes()

        cell_edges_flat, _ , _= tb.flatten(cells.mem_edges)

        cell_edges_flat = cells.um*np.asarray(cell_edges_flat)

        if zdata == None:
            z = np.ones(len(cell_edges_flat))
        elif zdata == 'random':
            z = np.random.random(len(cell_edges_flat))
        else:
            z = zdata

        if clrmap == None:
            clrmap = cm.rainbow

        coll = LineCollection(cell_edges_flat, array=z, cmap=clrmap)
        ax.add_collection(coll)

        ax.axis('equal')

        # Add a colorbar for the Line Collection
        if zdata != None:
            ax_cb = fig.colorbar(coll, ax=ax)

        ax.axis('equal')
        ax.autoscale_view(tight=True)

        return fig, ax, ax_cb

def plotConnectionData(cells, fig = None, ax=None, zdata=None,clrmap=None,colorbar = None, pickable=None):
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
        if fig==None:
            fig = plt.figure()# define the figure and axes instances
        if ax == None:
            ax = plt.subplot(111)
            #ax = plt.axes()

        if zdata == None:
            z = np.ones(len(cells.gap_jun_i))

        elif zdata == 'random':
            z = np.random.random(len(cells.gap_jun_i))

        else:
            z = zdata

        if clrmap == None:
            clrmap = cm.bone_r  # default colormap

         # Make a line collection and add it to the plot.

        con_segs = cells.cell_centres[cells.gap_jun_i]

        connects = p.um*np.asarray(con_segs)

        coll = LineCollection(connects, array=z, cmap=clrmap, linewidths=4.0, zorder=0)
        coll.set_clim(vmin=0.0,vmax=1.0)
        coll.set_picker(pickable)
        ax.add_collection(coll)

        # Plot the cell centres
        # ax.plot(p.um*cells.cell_centres[:,0],p.um*cells.cell_centres[:,1],'k.')

        #ax.quiver(s*self.gj_vects[:,0],s*self.gj_vects[:,1],s*self.gj_vects[:,2],s*self.gj_vects[:,3],z,zorder=5)

        ax.axis('equal')

        # Add a colorbar for the Line Collection
        if zdata != None and colorbar == 1:
            ax_cb = fig.colorbar(coll, ax=ax)
        else:
            ax_cb = None

        ax.autoscale_view(tight=True)

        return fig, ax, ax_cb

def plotBoundCells(points_flat,bflags,cells, fig=None, ax=None):
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
        if fig==None:
            fig = plt.figure()# define the figure and axes instances
        if ax == None:
            ax = plt.subplot(111)
            #ax = plt.axes()

        points_flat = np.asarray(points_flat)
        bflags = np.asarray(bflags)

        bpoints = points_flat[bflags]

        ax.plot(p.um*points_flat[:,0],p.um*points_flat[:,1],'k.')

        ax.plot(p.um*bpoints[:,0],p.um*bpoints[:,1],'r.')

        cell_edges_flat, _ , _= tb.flatten(cells.mem_edges)
        cell_edges_flat = p.um*np.asarray(cell_edges_flat)
        coll = LineCollection(cell_edges_flat,colors='k')
        coll.set_alpha(0.5)
        ax.add_collection(coll)

        ax.axis('equal')

        ax.autoscale_view(tight=True)

        return fig, ax

def plotVects(cells, fig=None, ax=None):
        """
        This function plots all unit vectors in the tissue system as a cross-check.
        Normals to cell membranes are shown as red arrows.
        Tangents to cell membranes are black arrows.
        Tangents to ecm edges are shown as green arrows.
        Cell membrane edges are drawn as blue lines.

        To plot streamline and vector plots with data use the pyplot quiver and streamplot functions, respectively.

        """

        if fig==None:
            fig = plt.figure()# define the figure and axes instances

        if ax == None:
            ax = plt.subplot(111)
            #ax = plt.axes()

        s = p.um

        ax.quiver(s*cells.mem_vects_flat[:,0],s*cells.mem_vects_flat[:,1],s*cells.mem_vects_flat[:,4],s*cells.mem_vects_flat[:,5],color='b')
        ax.quiver(s*cells.mem_vects_flat[:,0],s*cells.mem_vects_flat[:,1],s*cells.mem_vects_flat[:,2],s*cells.mem_vects_flat[:,3],color='g')
        ax.quiver(s*cells.ecm_vects[:,0],s*cells.ecm_vects[:,1],s*cells.ecm_vects[:,2],s*cells.ecm_vects[:,3],color='r')

        cell_edges_flat, _ , _= tb.flatten(cells.mem_edges)
        cell_edges_flat = cells.um*np.asarray(cell_edges_flat)
        coll = LineCollection(cell_edges_flat,colors='k')
        ax.add_collection(coll)

        ax.axis('equal')

        ax.autoscale_view(tight=True)

        return fig, ax

def exportData(cells,sim,p):

    results_path = p.sim_results
    os.makedirs(results_path, exist_ok=True)
    savedData = os.path.join(results_path, 'ExportedData.csv')

    cc_cell = []
    cc_env = []

    ci = p.plot_cell  # index of cell to get time data for

    # create the header, first entry will be time:
    headr = 'time_s'

    # next entry will be Vm:
    headr = headr + ',' + 'Vmem_mV'

    # create the header starting with cell concentrations
    for i in range(0,len(sim.ionlabel)):
        label = sim.ionlabel[i]
        headr = headr + ',' + 'cell_' + label + '_mmol/L'
        cc_m = [arr[i][ci] for arr in sim.cc_time]
        cc_m = np.asarray(cc_m)
        cc_cell.append(cc_m)

    # create the header moving on to env concentrations
    for i in range(0,len(sim.ionlabel)):
        label = sim.ionlabel[i]
        headr = headr + ',' + 'env_' + label + '_mmol/L'
        cc_m2 = [arr[i][ci] for arr in sim.cc_env_time]
        cc_m2 = np.asarray(cc_m2)
        cc_env.append(cc_m2)

    vm = [arr[ci]*1000 for arr in sim.vm_time]
    vm = np.asarray(vm)
    t = np.asarray(sim.time)
    cc_cell = np.asarray(cc_cell)
    cc_env = np.asarray(cc_env)

    IP3_time = [arr[ci] for arr in sim.cIP3_time]
    IP3_time = np.asarray(IP3_time)
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


    dataM = np.column_stack((t,vm,cc_cell.T, cc_env.T,IP3_time,dye_time,Ca_er))

    np.savetxt(savedData,dataM,delimiter = ',',header = headr)

def export2dData(simdata,cells,p):

    results_path = p.sim_results
    os.makedirs(results_path, exist_ok=True)
    savedData_2d = os.path.join(results_path, 'Exported2DData.csv')

    dataM = simdata
    np.savetxt(savedData_2d,dataM,delimiter=',')





