#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

import os.path
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection, PolyCollection
from betse.exceptions import BetseExceptionSimulation
from betse.science import filehandling as fh
from betse.science.cells import Cells
from betse.science.parameters import Parameters
from betse.science.plot import plot as viz
from betse.science.sim import Simulator
from betse.science.tissue.handler import TissueHandler
from betse.util.io.log import logs
from betse.util.path import files, paths
from betse import ignition, pathtree

# Initialize the current application before doing anything else.
ignition.init()
# at the moment only the default config file can be used
config_filename = pathtree.CONFIG_DEFAULT_FILENAME

# Validate and localize such filename.
files.die_unless_file(config_filename)
_config_filename = config_filename
_config_basename = paths.get_basename(_config_filename)

logs.log_info(
    'Plotting cell cluster with configuration file "{}".'.format(_config_basename))

p = Parameters(config_filename=_config_filename)  # create an instance of Parameters
p.I_overlay = False  # force the current overlay to be false as there's no data for it
sim = Simulator(p)

cells = Cells(p, worldtype='basic')

if files.is_file(cells.savedWorld):
    cells, _ = fh.loadWorld(cells.savedWorld)  # load the simulation from cache
    p.sim_ECM = cells.sim_ECM
    logs.log_info('Cell cluster loaded.')
else:
    raise BetseExceptionSimulation("Ooops! No such cell cluster file found to load!")

if p.sim_ECM is False:
    sim.baseInit_all(cells, p)
    dyna = TissueHandler(sim, cells, p)
    dyna.tissueProfiles(sim, cells, p)
else:
    sim.baseInit_all(cells, p)
    dyna = TissueHandler(sim, cells, p)
    dyna.tissueProfiles(sim, cells, p)

if p.autosave is True:
    images_path = p.init_results
    image_cache_dir = os.path.expanduser(images_path)
    os.makedirs(image_cache_dir, exist_ok=True)
    savedImg = os.path.join(image_cache_dir, 'fig_')

fig_tiss, ax_tiss, cb_tiss = viz.clusterPlot(
    p, dyna, cells, clrmap=p.default_cm)

if p.autosave is True:
    savename10 = savedImg + 'cluster_mosaic' + '.png'
    plt.savefig(savename10, format='png', transparent=True)

if p.turn_all_plots_off is False:
    plt.show(block=False)

if p.sim_ECM is True:
    plt.figure()
    ax99 = plt.subplot(111)
    plt.imshow(
        np.log10(sim.D_env_weight.reshape(cells.X.shape)),
        origin='lower',
        extent=[p.um * cells.xmin, p.um * cells.xmax, p.um * cells.ymin, p.um * cells.ymax],
        cmap=p.default_cm,
    )
    plt.colorbar()

    cell_edges_flat = p.um * cells.mem_edges_flat
    coll = LineCollection(cell_edges_flat, colors='k')
    coll.set_alpha(1.0)
    ax99.add_collection(coll)

    plt.title('Logarithm of Environmental Diffusion Weight Matrix')

    if p.autosave is True:
        savename10 = savedImg + 'env_diffusion_weights' + '.png'
        plt.savefig(savename10, format='png', transparent=True)

    if p.turn_all_plots_off is False:
        plt.show(block=False)

    plt.figure()
    plt.imshow(
        cells.maskM,
        origin='lower',
        extent=[p.um * cells.xmin, p.um * cells.xmax, p.um * cells.ymin, p.um * cells.ymax],
        cmap=p.default_cm,
    )
    plt.colorbar()
    plt.title('Cluster Masking Matrix')

    if p.autosave is True:
        savename = savedImg + 'cluster_mask' + '.png'
        plt.savefig(savename, format='png', transparent=True)

    if p.turn_all_plots_off is False:
        plt.show(block=False)

# plot gj
fig_x = plt.figure()
ax_x = plt.subplot(111)

if p.showCells is True:
    base_points = np.multiply(cells.cell_verts, p.um)
    col_cells = PolyCollection(base_points, facecolors='k', edgecolors='none')
    col_cells.set_alpha(0.3)
    ax_x.add_collection(col_cells)

con_segs = cells.nn_edges
connects = p.um * np.asarray(con_segs)
collection = LineCollection(connects, linewidths=1.0, color='b')
ax_x.add_collection(collection)
plt.axis('equal')
plt.axis([cells.xmin * p.um, cells.xmax * p.um, cells.ymin * p.um, cells.ymax * p.um])

ax_x.set_xlabel('Spatial x [um]')
ax_x.set_ylabel('Spatial y [um')
ax_x.set_title('Cell Connectivity Network')

if p.autosave is True:
    savename10 = savedImg + 'gj_connectivity_network' + '.png'
    plt.savefig(savename10, format='png', transparent=True)

if p.turn_all_plots_off is False:
    plt.show(block=False)
    plt.show()

else:
    logs.log_info(
        'Plots exported to init results folder '
        'defined in configuration file "{}".'.format(_config_basename))