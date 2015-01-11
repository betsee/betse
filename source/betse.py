#!/usr/bin/env python3
# Copyright 2014-2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.


import numpy as np
import scipy as sp
import scipy.spatial as sps
import matplotlib.pyplot as plt
from matplotlib.path import Path
import math
import time
from matplotlib import collections as col
from science import parameters, world
import matplotlib.cm as cm


start_time = time.time()  # get a start value for timing the simulation

# Define numerical constants for world set-up and simulation

const = parameters.Parameters()  # define the object that holds main parameters
const.wsx = 100e-6  # the x-dimension of the world space [m] recommended range 50 to 1000 um
const.wsy = 100e-6  # the y-dimension of the world space [m] recommended range 50 to 1000 um
const.rc = 5e-6  # radius of single cell
const.dc = const.rc * 2  # diameter of single cell
const.nx = int(const.wsx / const.dc)  # number of lattice sites in world x index
const.ny = int(const.wsy / const.dc)  # number of lattice sites in world y index
const.ac = 1e-6  # cell-cell separation for drawing
const.dc = const.rc * 2  # cell diameter
const.nl = 0.8  # noise level for the lattice
const.wsx = const.wsx + 5 * const.nl * const.dc  # readjust the world size for noise
const.wsy = const.wsy + 5 * const.nl * const.dc
const.search_d =1.5     # distance to search for nearest neighbours (relative to cell diameter dc) min 1.0 max 5.0
const.scale_cell = 0.9          # the amount to scale cell membranes in from ecm edges (only affects drawing)
const.cell_sides = 4      # minimum number of membrane domains per cell (must be >2)
const.scale_alpha = 1.0   # the amount to scale (1/d_cell) when calculating the concave hull (boundary search)
const.cell_height = 5.0e-6  # the height of a cell in the z-direction (for volume and surface area calculations)

cells = world.World(const, vorclose='circle',worldtype='full')
cells.makeWorld()

fig2, ax2, axcb2 = cells.plotPolyData(clrmap = cm.coolwarm,zdata='random')
ax2.set_ylabel('Spatial y [m]')
ax2.set_xlabel('Spatial x [m]')
ax2.set_title('Concentration of Foo Ion in Each Discrete Cell')
axcb2.set_label('Foo concentration [mol/m3]')
plt.show(block=False)

fig3, ax3, axcb3 = cells.plotMemData(clrmap = cm.coolwarm,zdata='random')
ax3.set_ylabel('Spatial y [um]')
ax3.set_xlabel('Spatial x [um]')
ax3.set_title('Foo Voltage on Discrete Membrane Domains')
axcb3.set_label('Foo membrane voltage [V]')
plt.show(block=False)

fig4, ax4, axcb4 =cells.plotConnectionData(zdata='random', clrmap=None)
ax4.set_ylabel('Spatial y [um]')
ax4.set_xlabel('Spatial x [um]')
ax4.set_title('Foo GJ permeability on Cell-Cell Connections')
axcb4.set_label('Foo GJ permeability [m/s]')
plt.show(block=False)

fig5, ax5, axcb5 = cells.plotVertData(cells.cell_verts,zdata='random',pointOverlay=True,edgeOverlay=True)
ax5.set_ylabel('Spatial y [um]')
ax5.set_xlabel('Spatial x [um]')
ax5.set_title('Foo Voltage of Membrane Domains Interpolated to Surface Plot')
axcb5.set_label('Foo membrane voltage [V]')
plt.show(block=False)

fig6, ax6 = cells.plotBoundCells(cells.mem_mids,cells.bflags_mems)
ax6.set_ylabel('Spatial y [um]')
ax6.set_xlabel('Spatial x [um]')
ax6.set_title('Membrane Domains Flagged at Cluster Boundary (red points)')
plt.show(block=False)

# fig7, ax7 = cells.plotBoundCells(cells.ecm_verts,cells.bflags_ecm)
# ax7.set_ylabel('Spatial y [um]')
# ax7.set_xlabel('Spatial x [um]')
# plt.show(block=False)

# fig8,ax8 = cells.plotVects()
# ax8.set_ylabel('Spatial y [um]')
# ax8.set_xlabel('Spatial x [um]')
# ax8.set_title('Normal and Tangent Vectors to Membrane Domains')
# plt.show(block=False)
#
# fig9, ax9, axcb9 = cells.plotCellData(zdata='random',pointOverlay=False,edgeOverlay=False)
# ax9.set_ylabel('Spatial y [um]')
# ax9.set_xlabel('Spatial x [um]')
# ax9.set_title('Foo Concentration in Cells Interpolated to Surface Plot')
# axcb9.set_label('Foo concentration [mol/m3]')
# plt.show(block=False)

print('total cell number', cells.cell_number)

print('average neighbours', int(cells.average_nn))

print('The simulation took', time.time() - start_time, 'seconds to complete')

plt.show()