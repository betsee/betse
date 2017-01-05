#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.


from betse.science.cells import Cells
from betse.science.parameters import Parameters
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
    'Seeding simulation with configuration file "{}".'.format(
        _config_basename))

p = Parameters(config_filename=_config_filename)  # create an instance of Parameters
p.I_overlay = False  # force the current overlay to be null
sim = Simulator(p)  # create an instance of Simulator as it's needed by plotting objects

if p.sim_ECM is False:

    cells = Cells(p, worldtype='basic')  # create an instance of world
    cells.containsECM = False
    logs.log_info('Cell cluster is being created...')
    cells.makeWorld(p)  # call function to create the world

    # define the tissue and boundary profiles for plotting:
    logs.log_info('Defining tissue and boundary profiles...')
    sim.baseInit_all(cells, p)
    dyna = TissueHandler(sim, cells, p)
    dyna.tissueProfiles(sim, cells, p)

    cells.redo_gj(dyna, p)  # redo gap junctions to isolate different tissue types

    if p.fluid_flow is True or p.deformation is True:  # if user desires fluid flow:

        # make a laplacian and solver for discrete transfers on closed, irregular cell network
        logs.log_info('Creating cell network Poisson solver...')
        cells.graphLaplacian(p)

    cells.save_cluster(p)

    logs.log_info('Cell cluster creation complete!')


else:

    cells = Cells(p, worldtype='full')  # create an instance of world
    cells.containsECM = True
    logs.log_info('Cell cluster is being created...')
    cells.makeWorld(p)  # call function to create the world

    # define the tissue and boundary profiles for plotting:
    logs.log_info('Defining tissue and boundary profiles...')
    sim.baseInit_all(cells, p)
    dyna = TissueHandler(sim, cells, p)
    dyna.tissueProfiles(sim, cells, p)

    cells.redo_gj(dyna, p)  # redo gap junctions to isolate different tissue types

    # make a laplacian and solver for discrete transfers on closed, irregular cell network
    if p.fluid_flow is True or p.deformation is True:
        # loggers.log_info('Creating cell network Poisson solvers...')
        cells.graphLaplacian(p)

    cells.save_cluster(p)

    logs.log_info('Cell cluster creation complete!')


sim.sim_info_report(cells, p)
