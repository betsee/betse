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
from pprint import pprint
from matplotlib import collections as col
import matplotlib.colors as colour
from science import parameters, world
import matplotlib.cm as cm


start_time = time.time()  # get a start value for timing the simulation

# Define numerical constants for world set-up and simulation

const = parameters.Parameters()  # define the object that holds main parameters
const.wsx = 200e-6  # the x-dimension of the world space [m] recommended range 100 to 1000 um
const.wsy = 200e-6  # the y-dimension of the world space [m] recommended range 100 to 1000 um
const.rc = 5e-6  # radius of single cell
const.dc = const.rc * 2  # diameter of single cell
const.nx = int(const.wsx / const.dc)  # number of lattice sites in world x index
const.ny = int(const.wsy / const.dc)  # number of lattice sites in world y index
const.ac = 1e-6  # cell-cell separation for drawing
const.dc = const.rc * 2  # cell diameter
const.nl = 0.8  # noise level for the lattice
const.wsx = const.wsx + 5 * const.nl * const.dc  # readjust the world size for noise
const.wsy = const.wsy + 5 * const.nl * const.dc
const.search_d =2.0     # distance to search for nearest neighbours (relative to cell diameter dc) min 1.0 max 5.0
const.sf = 0.6          # the amount to scale cell membranes in from ecm edges (only affects drawing)
const.cell_sides = 4      # minimum number of membrane domains per cell (must be >2)

cells = world.World(const, crop_mask=None, vorclose='circle')

# plt.figure(1)
#
# for region in cells.vor.regions:
#     polyg = cells.vor.vertices[region]
#     plt.fill(*zip(*polyg))
#
# plt.plot(cells.vor.vertices[:,0],cells.vor.vertices[:,1],'ko')
# plt.axis('equal')
# plt.show(block=False)


# fig2,ax2 = cells.plotPolyData(clrmap = cm.coolwarm,zdata='random')
# ax2.plot(cells.crop_pts[:,0],cells.crop_pts[:,1])
# plt.show(block=False)
#
# fig3, ax3 = cells.plotMemData(clrmap = cm.coolwarm,zdata='random')
# plt.show(block=False)
#
# fig4, ax4 =cells.plotConnectionData(zdata='random', clrmap=None)
# plt.show(block=False)
#
# fig5, ax5 = cells.plotVertData(cells.cell_verts,zdata='random',pointOverlay=None,edgeOverlay=None)
# plt.show(block=False)

fig6, ax6 = cells.plotBoundCells()
plt.show(block=False)


#print('The total cell number is',cells.cell_centres.shape[0])
# print(cells.mem_mids)
# print('---')
# print(cells.ecm_mids)
print('The simulation took', time.time() - start_time, 'seconds to complete')

plt.show()