#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

# FIXME do something about case where data doesn't vary so colorbar isn't useless

import matplotlib.pyplot as plt
import time
from betse.gui import interact
from betse.science.world import World
from betse.science.compute import Simulator
from betse.science.parameters import Parameters
from betse.science import visualize as viz
import matplotlib.cm as cm
from betse.science import filehandling as fh
import numpy as np

class MainDriver(object):

    def __init__(self):

        self.start_time = time.time()  # get a start value for timing the simulation

        self.cells = World(vorclose='circle',worldtype='full')  # always need instance of world
        self.cells.makeWorld()     # call functions to create the world

        #fh.saveSim(self.cells.savedWorld,self.cells)   # save the world to cache

        #cells = fh.loadWorld(cells.savedWorld)   # load a previously defined world from cache

        self.p = Parameters()

        self.sim = Simulator(self.p)   # whether running from scratch or loading, instance needs to be called

        self.sim.baseInit(self.cells, self.p)   # initialize data if working from scratch

        self.sim.runInit(self.cells,self.p)     # run and save an initialization if working from scratch

        #self.sim,self.cells, self.p = fh.loadSim(self.sim.savedInit)  # load an initialization from cache

        self.sim.runSim(self.cells,self.p,save=False)   # run and save the simulation

        #sim,cells,p = fh.loadSim(sim.savedSim)  # load the simulation from cache

        vdata0 =self.sim.vm_to*1000
        vdata = self.sim.vm_time[-1]*1000
        #
        # print(vdata0)
        # print(vdata)

        figI, axI, axcbI = viz.plotPolyData(self.cells,clrmap = cm.coolwarm,zdata=vdata0)
        figI, axI, _ = viz.plotConnectionData(self.cells, fig=figI, ax = axI, zdata=self.sim.gjopen, pickable=False)
        axI.set_ylabel('Spatial y [um]')
        axI.set_xlabel('Spatial x [um]')
        axI.set_title('Cell voltage at time zero')
        if axcbI != None:
            axcbI.set_label('Voltage [mV]')
        # for i, pt in enumerate(self.cells.cell_centres):
        #     axI.text(self.p.um*pt[0],self.p.um*pt[1],i)
        # if axcbI == None:
        #     md = np.mean(vdata0,axis=0)
        #     axI.text(0,0,md)
        plt.show(block=False)

        figV, axV, axcbV = viz.plotPolyData(self.cells,clrmap = cm.coolwarm,zdata=vdata)
        figV, axV, _ = viz.plotConnectionData(self.cells, fig=figV, ax = axV, zdata=self.sim.gjopen, pickable=False)
        axV.set_ylabel('Spatial y [um]')
        axV.set_xlabel('Spatial x [um]')
        axV.set_title('Voltage in Each Discrete Cell')
        if axcbV != None:
            axcbV.set_label('Voltage [mV]')
        # if axcbV == None:
        #     md = np.mean(vdata,axis=0)
        #     axV.text(0,0,md)
        plt.show(block=False)

        # interact.PolyPicker(self.cells,self.p,self.resume_after_picking)

        print('total cell number', self.cells.cell_number)
        print('average neighbours', int(self.cells.average_nn))
        print('The simulation took', time.time() - self.start_time, 'seconds to complete')

        plt.show()

    def resume_after_picking(self, poly_picker):


        print(poly_picker.picked_indices)
        #
        # ioni = sim.iNa
        # cdata = sim.cc_time[-1][ioni]
        # ionname = sim.ionlabel[ioni]
        #
        # figNa, axNa, axcbNa = viz.plotPolyData(cells, clrmap = cm.coolwarm,zdata=cdata)
        # figNa, axNa, _ = viz.plotConnectionData(cells, fig=figNa, ax = axNa, zdata=sim.gjopen)
        # axNa.set_ylabel('Spatial y [um]')
        # axNa.set_xlabel('Spatial x [um]')
        # tit = ionname + ' ' + 'Concentration in Cells'
        # lab = ionname + ' ' + '[mol/m3]'
        # axNa.set_title(tit)
        # axcbNa.set_label(lab)
        # plt.show(block=False)
        #
        # ioni = sim.iK
        # cdata = sim.cc_time[-1][ioni]
        # ionname = sim.ionlabel[ioni]
        #
        # figK, axK, axcbK = viz.plotPolyData(cells, clrmap = cm.coolwarm,zdata=cdata)
        # figK, axK, _ = viz.plotConnectionData(cells, fig=figK, ax = axK, zdata=sim.gjopen)
        # axK.set_ylabel('Spatial y [um]')
        # axK.set_xlabel('Spatial x [um]')
        # tit = ionname + ' ' + 'Concentration in Cells'
        # lab = ionname + ' ' + '[mol/m3]'
        # axK.set_title(tit)
        # axcbK.set_label(lab)
        # plt.show(block=False)
        #
        # figGJ, axGJ, axcbGJ =viz.plotConnectionData(cells, zdata=sim.gjopen, colorbar=1, pickable =True)
        # axGJ.set_ylabel('Spatial y [um]')
        # axGJ.set_xlabel('Spatial x [um]')
        # axGJ.set_title('Gap Junction Open Fraction')
        # axcbGJ.set_label('Gap Junction Open Fraction (1.0 = open, 0.0 = closed)')
        # plt.show(block=False)

        # fig1, ax1, axcb1 = viz.plotMemData(cells, clrmap = cm.coolwarm,zdata='random')
        # ax1.set_ylabel('Spatial y [um]')
        # ax1.set_xlabel('Spatial x [um]')
        # ax1.set_title('Foo Voltage on Discrete Membrane Domains')
        # axcb1.set_label('Foo membrane voltage [V]')
        # plt.show(block=False)
        #
        # fig6, ax6 = viz.plotBoundCells(cells.mem_mids_flat,cells.bflags_mems)
        # ax6.set_ylabel('Spatial y [um]')
        # ax6.set_xlabel('Spatial x [um]')
        # ax6.set_title('Membrane Domains Flagged at Cluster Boundary (red points)')
        # plt.show(block=False)
        #
        # fig7, ax7 = viz.plotBoundCells(cells.ecm_verts_unique,cells.bflags_ecm)
        # ax7.set_ylabel('Spatial y [um]')
        # ax7.set_xlabel('Spatial x [um]')
        # plt.show(block=False)
        #
        # fig8,ax8 = viz.plotVects()
        # ax8.set_ylabel('Spatial y [um]')
        # ax8.set_xlabel('Spatial x [um]')
        # ax8.set_title('Normal and Tangent Vectors to Membrane Domains')
        # plt.show(block=False)
        #
        # fig9, ax9, axcb9 = viz.plotCellData(zdata='random',pointOverlay=False,edgeOverlay=False)
        # ax9.set_ylabel('Spatial y [um]')
        # ax9.set_xlabel('Spatial x [um]')
        # ax9.set_title('Foo Concentration in Cells Interpolated to Surface Plot')
        # axcb9.set_label('Foo concentration [mol/m3]')
        # plt.show(block=False)

        # fig, ax = plt.subplots()
        # goop = 15                           # interested in the goopth cell
        #
        # cellpt = cells.cell_centres[goop]   # get its centre point
        # ax.plot(cellpt[0],cellpt[1],'ko')
        #
        # #cellverts = np.asarray(cells.cell_verts[goop])   # get the vertex points of the cell
        # #ax.plot(cellverts[:,0],cellverts[:,1])
        #
        # memedges = np.asarray(cells.mem_edges[goop])   # get the list of membrane edges
        #
        # ecmverts = np.asarray(cells.ecm_verts[goop])  # get the vertices of ecm points surrounding the cell
        # ax.plot(ecmverts[:,0],ecmverts[:,1],'g.')
        #
        # cellgjs = cells.cell2GJ_map[goop]
        #
        # neighcells = cells.cell_centres[cells.gap_jun_i[cellgjs]]
        #
        # coll = LineCollection(neighcells, colors='r')
        # ax.add_collection(coll)
        #
        # coll2 = LineCollection(memedges,colors='b')
        # ax.add_collection(coll2)
        #
        # cellecms_inds = cells.cell2ecm_map[goop]
        #
        # cellecm_verts = cells.ecm_verts_unique[cells.ecm_edges_i[cellecms_inds]]
        #
        # coll3 = LineCollection(cellecm_verts,colors='g')
        # ax.add_collection(coll3)
        #
        # cellecm_mids = cells.ecm_mids[cellecms_inds]
        # ax.plot(cellecm_mids[:,0],cellecm_mids[:,1],'go')
        #
        # memvects = []
        # for i, mem in enumerate(memedges):
        #     flatind = cells.rindmap_mem[goop][i]
        #     memvects.append(cells.mem_vects_flat[flatind])
        #
        # ecmvects = []
        # for val in cellecms_inds:
        #     ecmvects.append(cells.ecm_vects[val])
        #
        # memvects = np.asarray(memvects)
        # ecmvects = np.asarray(ecmvects)
        #
        # ax.quiver(memvects[:,0],memvects[:,1],memvects[:,2],memvects[:,3],color='b')
        # ax.quiver(ecmvects[:,0],ecmvects[:,1],ecmvects[:,2],ecmvects[:,3],color='g')
        #
        # plt.show(block = False)

        print('total cell number', self.cells.cell_number)
        print('average neighbours', int(self.cells.average_nn))
        print('The simulation took', time.time() - self.start_time, 'seconds to complete')

        # plt.show()

if __name__ == '__main__':
    MainDriver()

# --------------------( WASTELANDS                         )--------------------