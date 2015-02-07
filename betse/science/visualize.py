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

class AnimateCellData(object):
    """
    Animate color data on a plot of cells.

    """

    def __init__(self,cells,zdata_t,time,p,colormap=cm.rainbow, save=False,ani_repeat=False):

        self.zdata_t = zdata_t
        self.colormap = colormap
        self.time = time

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

        self.ax.add_collection(self.collection)

        self.ax.set_xlabel('Spatial x [um]')
        self.ax.set_ylabel('Spatial y [um')

        self.ax.autoscale_view()

        self.frames = len(self.zdata_t)

        ani = animation.FuncAnimation(self.fig, self.aniFunc,
            frames=self.frames, interval=66, repeat=ani_repeat)

        if save == True:
            # Encode such animation to disk. Naturally, this requires external
            # dependencies (e.g., "ffmpeg"). Unfortunately, the save() method
            # only prints non-fatal warnings rather than raising fatal
            # exceptions when such dependencies are not installed. We correct
            # this by temporarily converting warnings to exceptions for the
            # duration of such call. See also:
            #     https://docs.python.org/3/library/warnings.html
            with warnings.catch_warnings():
                warnings.simplefilter('error')

                #FIXME: Ideally, save() should detect which of ffmpeg and avconv
                #is installed and enable the appropriate writer. Unfortunately,
                #it currently requires we manually specify such backend as below.
                #This is bad, as Sess uses ffmpeg whereas Ally uses avconv (due
                #to differences between Linux distributions). Correct this by
                #manually detecting which of the two (if any) is in the current
                #$PATH and enabling the appropriate writer. Annoying, but trivial.

                ani.save('basic_animation.mp4', writer='avconv')

            print('Animation saved to file.')

        plt.show()


    def aniFunc(self,i):

        zz = self.zdata_t[i]

        self.collection.set_array(zz)
        tit = 'Simulation time' + ' ' + str(round(self.time[i],1)) + ' ' + 's'
        self.ax.set_title(tit)


class AnimateGJData(object):
    """
    Animate the gap junction open state as a function of time.
    """

    def __init__(self,cells,sim,p,colormap=cm.coolwarm, save=False,ani_repeat=False):

        self.zdata_t = sim.gjopen_time  # data array for gap junction coloring
        self.vdata_t = np.multiply(sim.vm_time,1000)   # data array for cell coloring
        self.colormap = colormap
        self.time = sim.time

        self.gjI_t = np.sign(sim.Igj_time)
        self.gjvects = cells.gj_vects

        self.fig = plt.figure()       # define figure
        self.ax = plt.subplot(111)    # define axes

        con_segs = cells.cell_centres[cells.gap_jun_i]
        connects = p.um*np.asarray(con_segs)
        self.collection = LineCollection(connects, array=self.zdata_t[0], cmap=cm.bone_r, linewidths=3.0, zorder=5)
        self.collection.set_clim(0.0,1.0)
        self.ax.add_collection(self.collection)

        # Next add a collection of cell polygons, with aimated voltage data
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
            vx,vy,self.zdata_t[0],zorder=10, cmap=cm.bone_r,clim=[0,1])

        self.cb.set_label('Voltage [mV]')
        self.ax.set_xlabel('Spatial x [um]')
        self.ax.set_ylabel('Spatial y [um')

        self.ax.autoscale_view()

        self.frames = len(self.zdata_t)

        ani = animation.FuncAnimation(self.fig, self.aniFunc,
               frames=self.frames, interval=100, repeat=ani_repeat)

        if save == True:

            ani.save('basic_animation.mp4')
            print('Animation saved to file.')

        plt.show()



    def aniFunc(self,i):

        zz = self.zdata_t[i]
        zv = self.vdata_t[i]

        vx = np.multiply(self.gjI_t[i],self.gjvects[:,2])
        vy = np.multiply(self.gjI_t[i],self.gjvects[:,3])

        self.collection.set_array(zz)
        self.coll2.set_array(zv)
        self.Qplot.set_UVC(vx,vy,zz)


        # self.ax.quiver(p.um*self.gjvects[:,0],p.um*self.gjvects[:,1],
        #     vx,vy,zz,zorder=10, cmap=cm.bone_r,clim=[0,1])

        tit = 'Simulation time' + ' ' + str(round(self.time[i],3)) + ' ' + 's'
        self.ax.set_title(tit)


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


def plotPolyData(cells, fig=None, ax=None, zdata = None,clrmap = None):
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
            clrmap = cm.rainbow

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


        ax.autoscale_view(tight=True)


        return fig,ax,ax_cb

def plotCellData(cells, fig=None, ax=None, zdata=None,clrmap=None,edgeOverlay = None,pointOverlay=None):
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
            z = zdata   # FIXME make an assertion to check for right data input

        if clrmap == None:
            clrmap = cm.rainbow


        triplt = ax.tripcolor(p.um*cells.cell_centres[:, 0], p.um*cells.cell_centres[:, 1],
            z,shading='gouraud', cmap=clrmap)

        ax.axis('equal')

        # Add a colorbar for the z-data
        if zdata != None:
            ax_cb = fig.colorbar(triplt, ax=ax)

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

