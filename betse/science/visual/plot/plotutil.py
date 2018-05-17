#!/usr/bin/env python3
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level utility functions specific to single-frame plots.
'''

import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
# from betse.util.io.log import logs
from matplotlib.collections import LineCollection, PolyCollection
from scipy import interpolate


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
    # membranes_midpoint_data = ((1/sample_size)*(cell_data_o/np.mean(cell_data_o)) )   # normalize the signal
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

        When using p.is_ecm, this plotting method overrides both plotPolyData and plotCellData.

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
            cell_edges_flat = p.um*cells.mem_edges_flat
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

            I_overlay(sim,cells,p,ax,plotIecm)

        return fig, ax, ax_cb


def plotPolyData(sim, cells, p, fig=None, ax=None, zdata = None, clrAutoscale = True, clrMin = None, clrMax = None,
    clrmap = None, number_cells=False, current_overlay = False,plotIecm=False):
        """
        Assigns color-data to each polygon in a cell cluster diagram and
        returns a plot instance (fig, axes).

        Parameters
        ----------
        cells : Cells
            Data structure holding all world information about cell geometry.
        zdata : optional[numpy.ndarray]
            A data array with each scalar entry corresponding to a cell's data
            value (for instance, concentration or voltage). If zdata is not
            supplied, the cells will be plotted with a uniform color; if zdata
            is the string `random`, a random data set will be created and
            plotted.
        clrAutoscale : optional[bool]
            If `True`, the colorbar is autoscaled to the max and min of zdata.
        clrMin : optional[float]
            Set the colorbar to a user-specified minimum value.
        clrMax : optional[float]
            Set the colorbar to a user-specified maximum value.
        clrmap : optional[matplotlib.cm]
            The colormap to use for plotting. Must be specified as cm.mapname.
            A list of available mapnames is supplied at:
            http://matplotlib.org/examples/color/colormaps_reference.html

        Returns
        -------
        fig, ax
            Matplotlib figure and axes instances for the plot.

        Notes
        -------
        This method Uses `matplotlib.collections.PolyCollection`,
        `matplotlib.cm`, `matplotlib.pyplot`, and numpy arrays and hence is
        computationally slow. Avoid calling this method for large collectives
        (e.g., larger than 500 x 500 um).
        """

        if fig is None:
            fig = plt.figure()# define the figure and axes instances
        if ax is None:
            ax = plt.subplot(111)
            #ax = plt.axes()

        if zdata is None:  # if user doesn't supply data
            z = np.ones(len(cells.cell_verts)) # create flat data for plotting

        #FIXME: This is a bit cumbersome. Ideally, a new "is_zdata_random"
        #boolean parameter defaulting to "False" should be tested, instead.
        #Whack-a-mole with a big-fat-pole!

        # If random data is requested, do so. To avoid erroneous and expensive
        # elementwise comparisons when "zdata" is neither None nor a string,
        # "zdata" must be guaranteed to be a string *BEFORE* testing this
        # parameter as a string. Numpy prints scary warnings otherwise: e.g.,
        #
        #     FutureWarning: elementwise comparison failed; returning scalar
        #     instead, but in the future will perform elementwise comparison
        elif isinstance(zdata, str) and zdata == 'random':
            z = np.random.random(len(cells.cell_verts)) # create some random data for plotting
        else:
            z = zdata

        # Make the polygon collection and add it to the plot.
        if clrmap is None:
            #clrmap = p.default_cm
            clrmap = cm.rainbow

        if p.showCells is True:
            coll, ax = cell_mosaic(z,ax,cells,p,p.default_cm)
        else:
            coll, ax = cell_mesh(z,ax,cells,p,p.default_cm)

        # points = np.multiply(cells.cell_verts, p.um)
        #
        # coll = PolyCollection(points, array=z, cmap=clrmap, edgecolors='none')
        # ax.add_collection(coll)
        ax.axis('equal')

        # Add a colorbar for the PolyCollection

        if zdata is not None and clrAutoscale is True:
            maxval = np.max(zdata,axis=0)
            minval = np.min(zdata,axis=0)

            coll.set_clim(minval,maxval)
            ax_cb = fig.colorbar(coll,ax=ax)

        elif clrAutoscale is False and zdata is not None:
            coll.set_clim(clrMin,clrMax)
            ax_cb = fig.colorbar(coll,ax=ax)

        elif zdata is None:
            ax_cb = None

        if number_cells is True:
            for i,cll in enumerate(cells.cell_centres):
                ax.text(p.um*cll[0],p.um*cll[1],i,ha='center',va='center')

        if current_overlay is True:
            streams, ax = I_overlay(sim,cells,p,ax,plotIecm)

        xmin = cells.xmin*p.um
        xmax = cells.xmax*p.um
        ymin = cells.ymin*p.um
        ymax = cells.ymax*p.um

        ax.axis([xmin,xmax,ymin,ymax])

        return fig,ax,ax_cb

def plotPrettyPolyData(data, sim, cells, p, fig=None, ax=None, clrAutoscale = True, clrMin = None, clrMax = None,
    clrmap = None, number_cells=False, current_overlay = False,plotIecm=False):
        """
        Assigns color-data to each polygon mem-mid, vertex and cell centre in a cell cluster
        diagram and returns a plot instance (fig, axes).

        Parameters
        ----------
        cells : Cells
            Data structure holding all world information about cell geometry.
        data : [numpy.ndarray]
            A data array with each scalar entry corresponding to a cell's data
            value (for instance, concentration or voltage) at cell membranes. If zdata is not
            supplied, the cells will be plotted with a uniform color; if zdata
            is the string `random`, a random data set will be created and
            plotted.
        clrAutoscale : optional[bool]
            If `True`, the colorbar is autoscaled to the max and min of zdata.
        clrMin : optional[float]
            Set the colorbar to a user-specified minimum value.
        clrMax : optional[float]
            Set the colorbar to a user-specified maximum value.
        clrmap : optional[matplotlib.cm]
            The colormap to use for plotting. Must be specified as cm.mapname.
            A list of available mapnames is supplied at:
            http://matplotlib.org/examples/color/colormaps_reference.html

        Returns
        -------
        fig, ax
            Matplotlib figure and axes instances for the plot.

        Notes
        -------
        This method Uses `matplotlib.collections.PolyCollection`,
        `matplotlib.cm`, `matplotlib.pyplot`, and numpy arrays and hence is
        computationally slow. Avoid calling this method for large collectives
        (e.g., larger than 500 x 500 um).
        """

        if fig is None:
            fig = plt.figure()# define the figure and axes instances
        if ax is None:
            ax = plt.subplot(111)

        # data processing -- map to verts:
        data_verts = np.dot(data, cells.matrixMap2Verts)

        # define colorbar limits for the PolyCollection

        if clrAutoscale is True:
            maxval = data_verts.max()
            minval = data_verts.min()
            # maxval = data_verts.max()
            # minval = data_verts.min()

        else:
            maxval = clrMax
            minval = clrMin


        # Make the polygon collection and add it to the plot.
        if clrmap is None:
            clrmap = p.default_cm

        if p.showCells is True:
            coll, ax = pretty_patch_plot(
                data_verts,ax,cells,p,p.default_cm, cmin=minval, cmax=maxval)
        else:
            coll, ax = cell_mesh(data,ax,cells,p,p.default_cm)

        # add a colorbar
        coll.set_clim(minval, maxval)
        ax_cb = fig.colorbar(coll, ax=ax)

        ax.axis('equal')

        if number_cells is True:
            for i,cll in enumerate(cells.cell_centres):
                ax.text(p.um*cll[0],p.um*cll[1],i,ha='center',va='center')

        if current_overlay is True:
            streams, ax = I_overlay(sim,cells,p,ax,plotIecm)

        xmin = cells.xmin*p.um
        xmax = cells.xmax*p.um
        ymin = cells.ymin*p.um
        ymax = cells.ymax*p.um

        ax.axis([xmin,xmax,ymin,ymax])

        return fig,ax,ax_cb

def plotVectField(Fx,Fy,cells,p,plot_ecm = False,title = 'Vector field',cb_title = 'Field [V/m]',
                    colorAutoscale = True, minColor = None, maxColor=None):

    fig = plt.figure()
    ax = plt.subplot(111)

    if plot_ecm is True:

        efield = np.sqrt(Fx**2 + Fy**2)

        msh = ax.imshow(efield,origin='lower', extent = [cells.xmin*p.um, cells.xmax*p.um, cells.ymin*p.um,
            cells.ymax*p.um],cmap=p.background_cm)

        vplot, ax = env_quiver(Fx,Fy,ax,cells,p)

        tit_extra = 'Extracellular'

    elif plot_ecm is False:

        efield = np.sqrt(Fx**2 + Fy**2)

        msh, ax = cell_mesh(efield,ax,cells,p,p.background_cm)

        vplot, ax = cell_quiver(Fx,Fy,ax,cells,p)

        tit_extra = 'Intracellular'

    ax.axis('equal')

    xmin = cells.xmin*p.um
    xmax = cells.xmax*p.um
    ymin = cells.ymin*p.um
    ymax = cells.ymax*p.um

    ax.axis([xmin,xmax,ymin,ymax])

    if colorAutoscale is False:
        msh.set_clim(minColor,maxColor)

    cb = fig.colorbar(msh)

    tit = title
    ax.set_title(tit)
    ax.set_xlabel('Spatial distance [um]')
    ax.set_ylabel('Spatial distance [um]')
    cb.set_label(cb_title)

    return fig, ax, cb


def plotStreamField(
    Fx,Fy,
    cells,
    p,
    plot_ecm: bool = False,
    title: str = 'Vector field',
    cb_title: str = 'Field [V/m]',
    show_cells: bool = False,
    colorAutoscale: bool = True,
    minColor = None,
    maxColor = None,
):

    fig = plt.figure()
    ax = plt.subplot(111)

    if plot_ecm:
        efield = np.sqrt(Fx**2 + Fy**2)
        # msh = ax.imshow(
        #     efield,
        #     origin='lower',
        #     extent=[cells.xmin*p.um, cells.xmax*p.um, cells.ymin*p.um, cells.ymax*p.um],
        #     cmap=p.background_cm,
        # )
        splot, ax = env_stream(Fx, Fy, ax, cells, p, cmap=p.background_cm)
        tit_extra = 'Extracellular'
    else:
        efield = np.sqrt(Fx**2 + Fy**2)

        # msh, ax = cell_mesh(efield,ax,cells,p,p.background_cm)
        splot, ax = cell_stream(
            Fx, Fy,
            ax, cells, p,
            show_cells=show_cells,
            cmap=p.background_cm,
        )
        tit_extra = 'Intracellular'

    ax.axis('equal')

    xmin = cells.xmin*p.um
    xmax = cells.xmax*p.um
    ymin = cells.ymin*p.um
    ymax = cells.ymax*p.um

    ax.axis([xmin, xmax, ymin, ymax])

    if not colorAutoscale:
        splot.lines.set_clim(minColor, maxColor)

    cb = fig.colorbar(splot.lines)

    ax.set_title(title)
    ax.set_xlabel('Spatial distance [um]')
    ax.set_ylabel('Spatial distance [um]')
    cb.set_label(cb_title)

    return fig, ax, cb


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

        cell_edges_flat = p.um*cells.mem_edges_flat

        if zdata is None:
            z = np.ones(len(cell_edges_flat))
        #FIXME: This is a bit cumbersome. Ideally, a new "is_zdata_random"
        #boolean parameter defaulting to "False" should be tested, instead.
        #Whack-a-mole with a big-fat-pole!

        # If random data is requested, do so. To avoid erroneous and expensive
        # elementwise comparisons when "zdata" is neither None nor a string,
        # "zdata" must be guaranteed to be a string *BEFORE* testing this
        # parameter as a string. Numpy prints scary warnings otherwise: e.g.,
        #
        #     FutureWarning: elementwise comparison failed; returning scalar
        #     instead, but in the future will perform elementwise comparison
        elif isinstance(zdata, str) and zdata == 'random':
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
        #FIXME: This is a bit cumbersome. Ideally, a new "is_zdata_random"
        #boolean parameter defaulting to "False" should be tested, instead.
        #Whack-a-mole with a big-fat-pole!

        # If random data is requested, do so. To avoid erroneous and expensive
        # elementwise comparisons when "zdata" is neither None nor a string,
        # "zdata" must be guaranteed to be a string *BEFORE* testing this
        # parameter as a string. Numpy prints scary warnings otherwise: e.g.,
        #
        #     FutureWarning: elementwise comparison failed; returning scalar
        #     instead, but in the future will perform elementwise comparison
        elif isinstance(zdata, str) and zdata == 'random':
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
        cell_edges_flat = p.um*cells.mem_edges_flat
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

def streamingCurrent(
    sim, cells, p,
    fig=None,
    ax=None,
    plot_Iecm=True,
    zdata=None,
    clrAutoscale=True,
    clrMin=None,
    clrMax=None,
    clrmap=cm.coolwarm,
    edgeOverlay=True,
    number_cells=False,
):

    # Define the figure and axes instances if needed.
    if fig is None:
        fig = plt.figure()
    if ax is None:
        ax = plt.subplot(111)

    ax.axis('equal')

    xmin = cells.xmin*p.um
    xmax = cells.xmax*p.um
    ymin = cells.ymin*p.um
    ymax = cells.ymax*p.um

    ax.axis([xmin,xmax,ymin,ymax])

    #FIXME: There's a fair amount of overlap between the following branches.
    #Treetops swaying in the contumely breeze!
    if p.is_ecm is False or plot_Iecm is False:

        # multiply by 100 to get units of uA/m2
        Jmag_M = 100*np.sqrt(
            sim.I_gj_x_time[-1]**2 + sim.I_gj_y_time[-1]**2) + 1e-30

        J_x = sim.I_gj_x_time[-1]/Jmag_M
        J_y = sim.I_gj_y_time[-1]/Jmag_M

        meshplot = plt.imshow(
            Jmag_M,
            origin='lower',
            extent=[xmin,xmax,ymin,ymax],
            cmap=clrmap,
        )

        ax.streamplot(
            cells.X*p.um, cells.Y*p.um, J_x, J_y,
            density=p.stream_density,
            linewidth=(3.0*Jmag_M/Jmag_M.max()) + 0.5,
            color='k',
            cmap=clrmap,
            arrowsize=5.0,
        )

        ax.set_title('Final gap junction current density')

    elif plot_Iecm is True:
        # multiply by 100 to get units of uA/m2
        Jmag_M = 100*np.sqrt(
            sim.I_tot_x_time[-1]**2 + sim.I_tot_y_time[-1]**2) + 1e-30

        J_x = sim.I_tot_x_time[-1]/Jmag_M
        J_y = sim.I_tot_y_time[-1]/Jmag_M

        meshplot = plt.imshow(
            Jmag_M,
            origin='lower',
            extent=[xmin,xmax,ymin,ymax],
            cmap=clrmap,
        )

        ax.streamplot(
            cells.X*p.um, cells.Y*p.um, J_x, J_y,
            density=p.stream_density,
            linewidth=(3.0*Jmag_M/Jmag_M.max()) + 0.5,
            color='k',
            cmap=clrmap,
            arrowsize=5.0,
        )

        ax.set_title('Final total currents')

    if clrAutoscale is True:
        ax_cb = fig.colorbar(meshplot,ax=ax)

    elif clrAutoscale is False:
        meshplot.set_clim(clrMin,clrMax)
        ax_cb = fig.colorbar(meshplot,ax=ax)

    # if p.showCells is True:
    #     # cell_edges_flat, _ , _= tb.flatten(cells.mem_edges)
    #     cell_edges_flat = p.um*cells.mem_edges_flat
    #     coll = LineCollection(cell_edges_flat,colors='k')
    #     coll.set_alpha(0.2)
    #     ax.add_collection(coll)

    if number_cells is True:

        for i,cll in enumerate(cells.cell_centres):
            ax.text(p.um*cll[0],p.um*cll[1],i,ha='center',va='center')

    return fig,ax,ax_cb


def I_overlay(sim,cells,p,ax,plotIecm = False):
    """
    Plots an overlay of simulated currents on an existing plot.

    Parameters
    -----------

    sim         Instance of sim module
    cells       Instance of cells module
    p           Instance of parameters module
    ax          Existing figure axis to plot currents on
    plotIecm    Plot total currents (True) or only those of gap junctions (False)

    Returns
    --------
    streams             Container for streamline plot
    ax                  Modified axis
    """

    if p.is_ecm is False or plotIecm is False:

        Ix = sim.I_cell_x_time[-1]
        Iy = sim.I_cell_y_time[-1]

        streams, ax = cell_stream(Ix, Iy,ax,cells,p)

        ax.set_title('(Intracellular current overlay)')

    elif plotIecm is True:

        Ix = sim.I_tot_x_time[-1]
        Iy = sim.I_tot_y_time[-1]

        streams, ax = env_stream(Ix, Iy,ax,cells,p)

        ax.set_title('(Environment current overlay)')

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


def cell_quiver(datax, datay, ax, cells, p):
    """
    Sets up a vector plot for cell-specific data on an existing axis.

    Parameters
    -----------
    datax, datay    Data defined on cell centres or membrane midpoints
    cells           Instance of cells module
    p               Instance of parameters module
    ax              Existing figure axis to plot currents on

    Returns
    --------
    vplot               Container for vector plot, plotted at cell centres
    ax                  Modified axis
    """

    if len(datax) == len(cells.mem_i):
        Fx = np.dot(cells.M_sum_mems,datax)/cells.num_mems
        Fy = np.dot(cells.M_sum_mems,datay)/cells.num_mems
    else:
        Fx = datax
        Fy = datay

    Fmag = np.sqrt(Fx**2 + Fy**2)

    # normalize the data:
    Fmag[Fmag == 0.0] = 1.0
    Fx = Fx/Fmag
    Fy = Fy/Fmag

    vplot = ax.quiver(
        p.um*cells.cell_centres[:,0],p.um*cells.cell_centres[:,1],Fx,Fy,
        pivot='mid', color=p.vcolor, units='x',
        headwidth=5, headlength=7, zorder=10)

    return vplot, ax


def env_quiver(datax,datay,ax,cells,p):
    """
    Sets up a vector plot for environmental data on an existing axis.

    Parameters
    -----------

    datax, datay    Data defined on environmental grid
    cells           Instance of cells module
    p               Instance of parameters module
    ax              Existing figure axis to plot currents on

    Returns
    --------
    vplot               Container for vector plot, plotted at environmental grid points
    ax                  Modified axis

    """
    F_mag = np.sqrt(datax**2 + datay**2)

    if F_mag.max() != 0.0:
        Fx = datax/F_mag.max()
        Fy = datay/F_mag.max()

    else:
        Fx = datax/F_mag.mean()
        Fy = datay/F_mag.mean()

    vplot = ax.quiver(p.um*cells.xypts[:,0], p.um*cells.xypts[:,1], Fx.ravel(),
        Fy.ravel(), pivot='mid',color = p.vcolor, units='x',headwidth=5, headlength = 7,zorder=10)

    return vplot, ax


def cell_stream(
    datax, datay, ax, cells, p, show_cells: bool = False, cmap = None):
    '''
    Add a streamline plot for cell-specific data to the passed axes.

    Parameters
    -----------
    datax, datay    Data defined on cell centres or membrane midpoints
    ax              Existing figure axes to plot currents on
    cells           Instance of cells module
    p               Instance of parameters module

    Returns
    --------
    streams             Container for stream plot, plotted at plot grid
    ax                  Modified axis
    '''

    if show_cells:
        cell_edges_flat = p.um*cells.mem_edges_flat
        coll = LineCollection(cell_edges_flat,colors='k')
        coll.set_alpha(0.3)
        ax.add_collection(coll)

    if datax.shape != cells.X.shape: # if the data hasn't been interpolated yet...
        Fx = interpolate.griddata(
            (cells.cell_centres[:,0], cells.cell_centres[:,1]),
            datax,
            (cells.X, cells.Y),
            fill_value=0,
            method=p.interp_type,
        )
        Fy = interpolate.griddata(
            (cells.cell_centres[:,0], cells.cell_centres[:,1]),
            datay,
            (cells.X, cells.Y),
            fill_value=0,
            method=p.interp_type,
        )

        Fx = Fx*cells.maskECM
        Fy = Fy*cells.maskECM

    else:
        Fx = datax
        Fy = datay

    # Magnitude of the passed vector field.
    Fmag = np.sqrt(Fx**2 + Fy**2)

    # Substitute all zero magnitudes by 1.0, enabling division by these
    # magnitudes without concern for division-by-zero exceptions.
    Fmag[Fmag == 0.0] = 1.0

    # Normalize this vector field to a unit vector field.
    Fx = Fx/Fmag
    Fy = Fy/Fmag

    # Streamline width.
    line_width = (3.0*Fmag/Fmag.max()) + 0.5

    # Color(s) of each streamline, either as a scalar *OR* an array of the same
    # shape as the "Fx" and "Fy" arrays.
    if cmap is None:
        stream_color = p.vcolor
    else:
        stream_color = Fmag

    streams = ax.streamplot(
        cells.X*p.um,
        cells.Y*p.um,
        Fx, Fy,
        density=p.stream_density,
        linewidth=line_width,
        arrowsize=5.0,
        color=stream_color,
        cmap=cmap,
    )

    return streams, ax


def mem_quiver(datax,datay,ax,cells,p, cmap=None):
    """
    Sets up a quiver plot for membrane vectors with tail at cell centres.

    Parameters
    -----------

    datax, datay    Data vectors defined on membrane midpoints
    cells           Instance of cells module
    p               Instance of parameters module
    ax              Existing figure axis to plot currents on

    Returns
    --------
    mvects             Container for quiver plot
    ax                  Modified axis

    """
    # scaleval = (p.um*cells.R.mean())*7.0
    scaleval = p.um*p.wsx*0.8

    mvects = ax.quiver(
        cells.cell_centres[:, 0][cells.mem_to_cells]*p.um,
        cells.cell_centres[:, 1][cells.mem_to_cells]*p.um,
        datax*cells.R[cells.mem_to_cells]*p.um,
        datay*cells.R[cells.mem_to_cells]*p.um,
        scale=scaleval,
        color=p.vcolor,
        cmap = cmap
    )

    ax.axis('equal')
    ax.axis([cells.xmin*p.um, cells.xmax*p.um, cells.ymin*p.um, cells.ymax*p.um])

    return mvects, ax


def env_stream(datax,datay,ax,cells,p, cmap=None):
    """
    Sets up a streamline plot for environmental data on an existing axis.

    Parameters
    -----------

    datax, datay    Data defined on environmental grid
    cells           Instance of cells module
    p               Instance of parameters module
    ax              Existing figure axis to plot currents on

    Returns
    --------
    streams             Container for stream plot
    ax                  Modified axis

    """

    Fmag = np.sqrt(datax**2 + datay**2) + 1e-30

    if Fmag.all() != 0.0:
        Fx = datax/Fmag
        Fy = datay/Fmag

    if Fmag.max() != 0.0:

        lw = (3.0*Fmag/Fmag.max()) + 0.5

    else:
        lw = 3.0

    # if datax.shape == cells.X.shape:

    streams = ax.streamplot(cells.X*p.um,cells.Y*p.um, Fx, Fy, arrowsize=5.0, density=p.stream_density,
            linewidth=lw,color=Fmag, cmap=cmap)

    # elif datax.shape == cells.X.shape:
    #
    #     streams = ax.streamplot(cells.X*p.um,cells.Y*p.um, Fx, Fy,density=p.stream_density,
    #         linewidth=lw,color=p.vcolor,arrowsize=1.5)

    # else:
    #     raise BetseFunctionException("Data input to env_streams function must be \n shaped as cells.X or cells.Xgrid.")

    return streams, ax


def cell_mesh(data, ax, cells, p, clrmap):

    # If the data is defined on membrane midpoints get membrane midpoint coordinates
    if len(data) == len(cells.mem_i):
        # data = np.dot(cells.M_sum_mems,data)/cells.num_mems
        xi = p.um*cells.mem_mids_flat[:,0]
        yi = p.um*cells.mem_mids_flat[:,1]

    elif len(data) == len(cells.cell_i): # otherwise

        xi = p.um*cells.cell_centres[:,0]
        yi = p.um*cells.cell_centres[:,1]


    # data_grid = np.zeros(len(cells.voronoi_centres))
    # data_grid[cells.cell_to_grid] = data

    msh = ax.tripcolor(
        xi,
        yi,
        data,
        shading='gouraud',
        cmap=clrmap,
    )

    return msh, ax

#FIXME: Obsoleted by the more general-purpose, reliable, and efficient
#"betse.science.plot.layer.vector.lyrvecabrupt.LayerCellsVectorAbruptMembranes"
#class. Replace all remaining calls to this function with usage of that class;
#then, remove this function.
def pretty_patch_plot(
    data_verts, ax, cells, p, clrmap,
    cmin=None,
    cmax=None,
    use_other_verts=None
):
    """
    Maps data on mem midpoints to vertices, and
    uses tripcolor on every cell patch to create a
    lovely gradient. Slow but beautiful!

    data:   mem midpoint data for plotting (e.g vm)
    ax:     plot axis
    cells:  cells object
    p:      parameters object
    clrmap: colormap
    cmin, cmax   clim values for the data's colormap

    """


    # data_verts = data

    # colormap clim
    if cmin is None:
        amin = data_verts.min()
        amax = data_verts.max()
    else:
        amin = cmin
        amax = cmax

    # amin = amin + 0.1 * np.abs(amin)
    # amax = amax - 0.1 * np.abs(amax)

    # collection of cell patchs at vertices:
    if use_other_verts is None:
        cell_faces = np.multiply(cells.cell_verts, p.um)
    else:
        cell_faces = np.multiply(use_other_verts, p.um)

    # Cell membrane (Vmem) plotter (slow but beautiful!)
    for i in range(len(cell_faces)):
        x = cell_faces[i][:, 0]
        y = cell_faces[i][:, 1]

        # Average color value of each cell membrane, situated at the midpoint
        # of that membrane. This parameter is referred to as "C" in both the
        # documentation and implementation of the tripcolor() function.
        dati = data_verts[cells.cell_to_mems[i]]

        # "matplotlib.collections.TriMesh" instance providing the
        # Gouraud-shaded triangulation mesh for the non-triangular vertices of
        # this cell from the Delaunay hull of these vertices.
        col_cell = ax.tripcolor(x, y, dati, shading='gouraud', cmap=clrmap)

        #FIXME: No need to manually call set_clim() here. Instead, pass the
        #"vmin=amin, vmax=amax" parameters to the above tripcolor() call.

        # Autoscale this mesh's colours as desired.
        col_cell.set_clim(amin, amax)

    return col_cell, ax


def env_mesh(data, ax, cells, p, clrmap, ignore_showCells=False):
    """
    Sets up a mesh plot for environmental data on an existing axis.

    Parameters
    -----------

    data            Data defined on environmental grid
    cells           Instance of cells module
    p               Instance of parameters module
    ax              Existing figure axis to plot currents on

    Returns
    --------
    mesh_plot           Container for mesh plot
    ax                  Modified axis

    """

    # if p.plotMask is True:
    #     data = ma.masked_array(data, np.logical_not(cells.maskM))

    mesh_plot = ax.imshow(data,origin='lower',
                extent=[p.um*cells.xmin,p.um*cells.xmax,p.um*cells.ymin,p.um*cells.ymax],cmap=clrmap)

    if p.showCells is True and ignore_showCells is False:
        cell_edges_flat = p.um*cells.mem_edges_flat
        coll = LineCollection(cell_edges_flat,colors='k')
        coll.set_alpha(0.5)
        ax.add_collection(coll)

    return mesh_plot, ax


def cell_mosaic(
    data,
    ax: 'matplotlib.axes.Axes',
    cells: 'Cells',
    p: 'Parameters',
    clrmap: 'matplotlib.colors.Colormap',
) -> (PolyCollection, 'matplotlib.axes.Axes'):
    """
    Sets up a mosaic plot for cell data on an existing axis.

    Parameters
    -----------
    data            Data defined on environmental grid
    cells           Instance of cells module
    p               Instance of parameters module
    ax              Existing figure axis to plot currents on

    Returns
    --------
    collection          Container for mosaic plot
    ax                  Modified axis
    """

    # define a polygon collection based on individual cell polygons
    points = np.multiply(cells.cell_verts, p.um)
    collection =  PolyCollection(points, cmap=clrmap, edgecolors='none')
    collection.set_array(data)
    ax.add_collection(collection)

    return collection, ax
