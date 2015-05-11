#!/usr/bin/env python3
# Copyright 2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.


# FIXME after that you need to do electrodiffusion in the ECM spaces...
# FIXME boundary ECMs...
# FIXME currents in ECM and gj networks...

import numpy as np
import os, os.path
import copy
from random import shuffle
from betse.science import filehandling as fh
from betse.science import visualize as viz
from betse.science import toolbox as tb
import matplotlib.pyplot as plt
from betse.exceptions import BetseExceptionSimulation
from betse.util.io import loggers
import time


class Simulator(object):
    """
    Contains the main routines used in the simulation of networked cell bioelectrical activity.
    All methods are based on matrix mathematics and are implemented in Numpy for speed.

    Fields
    ------
    self.savedInit      path and filename for saving an initialization sim
    self.savedSim       path and filename for saving a sim

    self.cc_cells       nested list of each sim ion concentration for each cell
    self.cc_env         nested list of each sim ion concentration in the environment
    self.zs             list of valence state of each sim ion
    self.Dm_cells       list of membrane diffusion constant of each ion
    self.movingIons     list of electro-diffusing ions in sim
    self.ionlabel       label of electro-diffusing ions in sim
    self.gjopen_time    nested list of gap junction open state for each cell for each sampled time in sim
    self.cc_time        nested list of each sim ion concentration for each cell at each sampled time in sim
    self.vm_time        nested list of voltage for each cell at each sampled time in sim
    self.vm_to          initial voltage of cells in simulation (typically values returned from an initialization run)
    self.fgj_time       nested list of ion fluxes through gap junctions at each sampled time in sim
    self.time           sampled time for sim

    Methods
    -------
    fileInit(p)                 Prepares save paths for initialization and simulation runs.

    baseInit(cells,p)           Prepares zeroed or base-initialized data structures for a full init or sim run.

    runInit(cells,p)            Runs and saves an initialization for world specified by cells and parameters in p.

    runSim(cells,p,save=True)   Runs a simulation for world specified by cells and parameters in p.
                                Optional save for save=True (default).

    """

    def __init__(self,p):

        self.fileInit(p)

    def fileInit(self,p):

        """
        Initializes file saving and loading directory as the betse cach.
        For now, automatically assigns file names, but later, will allow
        user-specified file names.

        """

        #FIXME: Uncomment to make required directories. Yay!
        # from betse.util.path import dirs
        # dirs.make_parent_unless_found(p.saved_init_file, p.saved_run_file)

        # Make the BETSE-specific cache directory if not found.
        betse_cache_dir = os.path.expanduser(p.init_path)
        os.makedirs(betse_cache_dir, exist_ok=True)

        sim_cache_dir = os.path.expanduser(p.sim_path)
        os.makedirs(sim_cache_dir, exist_ok=True)

        # Define data paths for saving an initialization and simulation run:
        self.savedInit = os.path.join(betse_cache_dir, p.init_filename)
        self.savedSim = os.path.join(sim_cache_dir, p.sim_filename)

    def baseInit(self,cells,p):
        """
        Creates a host of initialized data matrices for the main simulation,
        including intracellular and environmental concentrations, voltages, and specific
        diffusion constants.

        """

        # Identity matrix to easily make matrices out of scalars
        self.id_cells = np.ones(len(cells.cell_i))

        self.cc_cells = []  # cell concentrations initialized
        self.cc_env = []   # environmental concentrations initialized
        self.cc_er = []   # endoplasmic reticulum ion concentrations
        self.zs = []   # ion valence state initialized
        self.z_er = []  # ion valence states of er ions
        self.Dm_cells = []              # membrane diffusion constants initialized
        self.D_free = []                 # a list of single-valued free diffusion constants for each ion
        self.Dm_er = []                  # a list of endoplasmic reticulum membrane states
        self.movingIons = []            # moving ions indices
        self.ionlabel = {}              # dictionary to hold ion label names

        self.T = p.T                # set the base temperature for the simulation

        i = -1                           # an index to track place in ion list

        # Initialize cellular concentrations of ions:
        if p.ions_dict['Na'] == 1:

            i = i+1

            self.iNa = i
            self.movingIons.append(self.iNa)
            self.ionlabel[self.iNa] = 'sodium'

            cNa_cells = np.zeros(len(cells.cell_i))
            cNa_cells[:]=p.cNa_cell

            cNa_env = np.zeros(len(cells.cell_i))
            cNa_env[:]=p.cNa_env

            DmNa = np.zeros(len(cells.cell_i))
            DmNa[:] = p.Dm_Na

            self.cc_cells.append(cNa_cells)
            self.cc_env.append(cNa_env)
            self.zs.append(p.z_Na)
            self.Dm_cells.append(DmNa)
            self.D_free.append(p.Do_Na)


        if p.ions_dict['K'] == 1:

            i= i+ 1

            self.iK = i
            self.movingIons.append(self.iK)
            self.ionlabel[self.iK] = 'potassium'

            cK_cells = np.zeros(len(cells.cell_i))
            cK_cells[:]=p.cK_cell

            cK_env = np.zeros(len(cells.cell_i))
            cK_env[:]=p.cK_env

            DmK = np.zeros(len(cells.cell_i))
            DmK[:] = p.Dm_K

            self.cc_cells.append(cK_cells)
            self.cc_env.append(cK_env)
            self.zs.append(p.z_K)
            self.Dm_cells.append(DmK)
            self.D_free.append(p.Do_K)


        if p.ions_dict['Cl'] == 1:

            i =i+1

            self.iCl = i
            self.movingIons.append(self.iCl)
            self.ionlabel[self.iCl] = 'chloride'

            cCl_cells = np.zeros(len(cells.cell_i))
            cCl_cells[:]=p.cCl_cell

            cCl_env = np.zeros(len(cells.cell_i))
            cCl_env[:]=p.cCl_env

            DmCl = np.zeros(len(cells.cell_i))
            DmCl[:] = p.Dm_Cl

            self.cc_cells.append(cCl_cells)
            self.cc_env.append(cCl_env)
            self.zs.append(p.z_Cl)
            self.Dm_cells.append(DmCl)
            self.D_free.append(p.Do_Cl)


        if p.ions_dict['Ca'] == 1:

            i =i+1

            self.iCa = i
            self.movingIons.append(self.iCa)
            self.ionlabel[self.iCa] = 'calcium'

            cCa_cells = np.zeros(len(cells.cell_i))
            cCa_cells[:]=p.cCa_cell

            cCa_env = np.zeros(len(cells.cell_i))
            cCa_env[:]=p.cCa_env

            DmCa = np.zeros(len(cells.cell_i))
            DmCa[:] = p.Dm_Ca

            self.cc_cells.append(cCa_cells)
            self.cc_env.append(cCa_env)
            self.zs.append(p.z_Ca)
            self.Dm_cells.append(DmCa)
            self.D_free.append(p.Do_Ca)

            if p.ions_dict['Ca'] ==1:
                cCa_er = np.zeros(len(cells.cell_i))
                cCa_er[:]=p.cCa_er
                self.cc_er.append(cCa_er)
                self.z_er.append(p.z_Ca)
                self.Dm_er.append(p.Dm_Ca)

        if p.ions_dict['P'] == 1:

            i =i+1

            self.iP = i
            self.ionlabel[self.iP] = 'proteins'

            cP_cells = np.zeros(len(cells.cell_i))
            cP_cells[:]=p.cP_cell

            cP_env = np.zeros(len(cells.cell_i))
            cP_env[:]=p.cP_env

            DmP = np.zeros(len(cells.cell_i))
            DmP[:] = p.Dm_P

            self.cc_cells.append(cP_cells)
            self.cc_env.append(cP_env)
            self.zs.append(p.z_P)
            self.Dm_cells.append(DmP)
            self.D_free.append(p.Do_P)

        if p.ions_dict['M'] == 1:

            i =i+1

            self.iM = i
            self.movingIons.append(self.iM)
            self.ionlabel[self.iM] = 'charge balance anion'

            cM_cells = np.zeros(len(cells.cell_i))
            cM_cells[:]=p.cM_cell

            cM_env = np.zeros(len(cells.cell_i))
            cM_env[:]=p.cM_env

            DmM = np.zeros(len(cells.cell_i))
            DmM[:] = p.Dm_M

            self.cc_cells.append(cM_cells)
            self.cc_env.append(cM_env)
            self.zs.append(p.z_M)
            self.Dm_cells.append(DmM)
            self.D_free.append(p.Do_M)

            if p.ions_dict['Ca'] ==1:
                cM_er = np.zeros(len(cells.cell_i))
                cM_er[:]=p.cM_er
                self.cc_er.append(cM_er)
                self.z_er.append(p.z_M)
                self.Dm_er.append(p.Dm_M)

        if p.ions_dict['H'] == 1:

            i =i+1

            self.iH = i

            #self.movingIons.append(self.iH)
            self.ionlabel[self.iH] = 'protons'

            # initialize the carbonic acid for the carbonate buffer
            self.cHM_cells = np.zeros(len(cells.cell_i))
            self.cHM_cells[:] = 0.03*p.CO2

            self.cHM_env = np.zeros(len(cells.cell_i))
            self.cHM_env[:] = 0.03*p.CO2

            cH_cells = np.zeros(len(cells.cell_i))
            # cH_cells[:]=p.cH_cell
            pH_cell = 6.1 + np.log10(cM_cells/self.cHM_cells)
            cH_cells = (10**(-pH_cell))*1000  # units mmol/L


            cH_env = np.zeros(len(cells.cell_i))
            # cH_env[:]=p.cH_env
            pH_env = 6.1 + np.log10(cM_env/self.cHM_env)
            cH_env = (10**(-pH_env))*1000 # units mmol/L

            DmH = np.zeros(len(cells.cell_i))
            DmH[:] = p.Dm_H

            self.cc_cells.append(cH_cells)
            self.cc_env.append(cH_env)
            self.zs.append(p.z_H)
            self.Dm_cells.append(DmH)
            self.D_free.append(p.Do_H)

        # Initialize membrane thickness:
        self.tm = np.zeros(len(cells.cell_i))
        self.tm[:] = p.tm

        # Initialize environmental volume:
        self.envV = np.zeros(len(cells.cell_i))
        self.envV[:] = p.vol_env

        flx = np.zeros(len(cells.gj_i))
        self.fluxes_gj = [flx,flx,flx,flx,flx,flx,flx]   # stores gj fluxes for each ion
        self.gjopen_time = []   # stores gj open fraction at each time
        self.fgj_time = []      # stores the gj fluxes for each ion at each time
        self.Igj =[]            # current for each gj
        self.Igj_time = []      # current for each gj at each time

        self.Dm_er = copy.deepcopy(self.Dm_cells)

        self.vm_to = np.zeros(len(cells.cell_i))
        self.v_er = np.zeros(len(cells.cell_i))

        self.NaKATP_block = np.ones(len(cells.cell_i))  # initialize NaKATP blocking vector
        self.HKATP_block = np.ones(len(cells.cell_i))  # initialize HKATP blocking vector
        self.CaATP_block = np.ones(len(cells.cell_i))  # initialize CaATP (plasm membrane) blocking vector
        self.CaER_block = np.ones(len(cells.cell_i)) # initialize CaATP (ER membrane) blocking vector
        self.VATP_block = np.ones(len(cells.cell_i)) # initialize CaATP (ER membrane) blocking vector

        self.cc_er = np.asarray(self.cc_er)   # initialize endoplasmic reticulum concentration array
        self.cc_er_to = np.copy(self.cc_er)

        self.cDye_cell = np.zeros(len(cells.cell_i))   # initialize voltage sensitive dye array for cell and env't

        self.cDye_env = np.zeros(len(cells.cell_i))
        self.cDye_env[:] = p.cDye_to

        # add channel noise to the model:
        self.channel_noise_factor = np.random.random(len(cells.cell_i))
        self.Dm_cells[self.iK] = (p.channel_noise_level*self.channel_noise_factor + 1)*self.Dm_cells[self.iK]

        # add a random walk on protein concentration to generate dynamic noise:
        self.protein_noise_factor = p.dynamic_noise_level*(np.random.random(len(cells.cell_i)) - 0.5)

        if p.ions_dict['P']==1:
            self.cc_cells[self.iP] = self.cc_cells[self.iP]*(1+ self.protein_noise_factor)

        loggers.log_info('This world contains '+ str(cells.cell_number) + ' cells.')
        loggers.log_info('Each cell has an average of ' + str(round(cells.average_nn,2)) + ' nearest-neighbours.')
        loggers.log_info('You are running the ion profile: ' + p.ion_profile)
        loggers.log_info('Ions in this simulation: ' + str(self.ionlabel))
        loggers.log_info('If you have selected features that use omitted ions they will be ignored.')

    def baseInit_ECM(self,cells,p):

        # Identity matrix to easily make matrices out of scalars
        self.id_cells = np.ones(len(cells.cell_i))

        self.cc_cells = []  # cell concentrations initialized
        # self.cc_env = []   # environmental concentrations initialized
        self.cc_er = []   # endoplasmic reticulum ion concentrations

        self.cc_ecm = []  # extracellular spaces ion concentrations

        self.zs = []   # ion valence state initialized
        self.z_er = []  # ion valence states of er ions

        self.Dm_mems = []              # membrane diffusion constants initialized

        self.D_free = []                 # a list of single-valued free diffusion constants for each ion
        self.Dm_er = []                  # a list of endoplasmic reticulum membrane states
        self.movingIons = []            # moving ions indices
        self.ionlabel = {}              # dictionary to hold ion label names

        self.T = p.T                # set the base temperature for the simulation

        i = -1                           # an index to track place in ion list

        # Initialize cellular concentrations of ions:
        if p.ions_dict['Na'] == 1:

            i = i+1

            self.iNa = i
            self.movingIons.append(self.iNa)
            self.ionlabel[self.iNa] = 'sodium'

            cNa_cells = np.zeros(len(cells.cell_i))
            cNa_cells[:]=p.cNa_cell

            cNa_ecm = np.zeros(len(cells.ecm_i))
            cNa_ecm[:]=p.cNa_env

            DmNa = np.zeros(len(cells.mem_i))
            DmNa[:] = p.Dm_Na

            self.cc_cells.append(cNa_cells)
            self.cc_ecm.append(cNa_ecm)
            self.zs.append(p.z_Na)
            self.Dm_mems.append(DmNa)
            self.D_free.append(p.Do_Na)


        if p.ions_dict['K'] == 1:

            i= i+ 1

            self.iK = i
            self.movingIons.append(self.iK)
            self.ionlabel[self.iK] = 'potassium'

            cK_cells = np.zeros(len(cells.cell_i))
            cK_cells[:]=p.cK_cell

            cK_ecm = np.zeros(len(cells.ecm_i))
            cK_ecm[:]=p.cK_env

            DmK = np.zeros(len(cells.mem_i))
            DmK[:] = p.Dm_K

            self.cc_cells.append(cK_cells)
            self.cc_ecm.append(cK_ecm)
            self.zs.append(p.z_K)
            self.Dm_mems.append(DmK)
            self.D_free.append(p.Do_K)


        if p.ions_dict['Cl'] == 1:

            i =i+1

            self.iCl = i
            self.movingIons.append(self.iCl)
            self.ionlabel[self.iCl] = 'chloride'

            cCl_cells = np.zeros(len(cells.cell_i))
            cCl_cells[:]=p.cCl_cell

            cCl_ecm = np.zeros(len(cells.ecm_i))
            cCl_ecm[:]=p.cCl_env

            DmCl = np.zeros(len(cells.mem_i))
            DmCl[:] = p.Dm_Cl

            self.cc_cells.append(cCl_cells)
            self.cc_ecm.append(cCl_ecm)
            self.zs.append(p.z_Cl)
            self.Dm_mems.append(DmCl)
            self.D_free.append(p.Do_Cl)


        if p.ions_dict['Ca'] == 1:

            i =i+1

            self.iCa = i
            self.movingIons.append(self.iCa)
            self.ionlabel[self.iCa] = 'calcium'

            cCa_cells = np.zeros(len(cells.cell_i))
            cCa_cells[:]=p.cCa_cell

            cCa_ecm = np.zeros(len(cells.ecm_i))
            cCa_ecm[:]=p.cCa_env

            DmCa = np.zeros(len(cells.mem_i))
            DmCa[:] = p.Dm_Ca

            self.cc_cells.append(cCa_cells)
            self.cc_ecm.append(cCa_ecm)
            self.zs.append(p.z_Ca)
            self.Dm_mems.append(DmCa)
            self.D_free.append(p.Do_Ca)

            if p.ions_dict['Ca'] ==1:
                cCa_er = np.zeros(len(cells.cell_i))
                cCa_er[:]=p.cCa_er
                self.cc_er.append(cCa_er)
                self.z_er.append(p.z_Ca)
                self.Dm_er.append(p.Dm_Ca)

        if p.ions_dict['P'] == 1:

            i =i+1

            self.iP = i
            self.ionlabel[self.iP] = 'proteins'

            cP_cells = np.zeros(len(cells.cell_i))
            cP_cells[:]=p.cP_cell

            cP_ecm = np.zeros(len(cells.ecm_i))
            cP_ecm[:]=p.cP_env

            DmP = np.zeros(len(cells.mem_i))
            DmP[:] = p.Dm_P

            self.cc_cells.append(cP_cells)
            self.cc_ecm.append(cP_ecm)
            self.zs.append(p.z_P)
            self.Dm_mems.append(DmP)
            self.D_free.append(p.Do_P)

        if p.ions_dict['M'] == 1:

            i =i+1

            self.iM = i
            self.movingIons.append(self.iM)
            self.ionlabel[self.iM] = 'charge balance anion'

            cM_cells = np.zeros(len(cells.cell_i))
            cM_cells[:]=p.cM_cell

            cM_ecm = np.zeros(len(cells.ecm_i))
            cM_ecm[:]=p.cM_env

            DmM = np.zeros(len(cells.mem_i))
            DmM[:] = p.Dm_M

            self.cc_cells.append(cM_cells)
            self.cc_ecm.append(cM_ecm)
            self.zs.append(p.z_M)
            self.Dm_mems.append(DmM)
            self.D_free.append(p.Do_M)

            if p.ions_dict['Ca'] ==1:
                cM_er = np.zeros(len(cells.cell_i))
                cM_er[:]=p.cM_er
                self.cc_er.append(cM_er)
                self.z_er.append(p.z_M)
                self.Dm_er.append(p.Dm_M)

        if p.ions_dict['H'] == 1:

            i =i+1

            self.iH = i

            #self.movingIons.append(self.iH)
            self.ionlabel[self.iH] = 'protons'

            # initialize the carbonic acid for the carbonate buffer
            self.cHM_cells = np.zeros(len(cells.cell_i))
            self.cHM_cells[:] = 0.03*p.CO2

            self.cHM_ecm = np.zeros(len(cells.ecm_i))
            self.cHM_ecm[:] = 0.03*p.CO2

            cH_cells = np.zeros(len(cells.cell_i))
            # cH_cells[:]=p.cH_cell
            pH_cell = 6.1 + np.log10(cM_cells/self.cHM_cells)
            cH_cells = (10**(-pH_cell))*1000  # units mmol/L


            cH_ecm = np.zeros(len(cells.ecm_i))

            pH_ecm = 6.1 + np.log10(cM_ecm/self.cHM_ecm)
            cH_ecm = (10**(-pH_ecm))*1000 # units mmol/L

            DmH = np.zeros(len(cells.mem_i))
            DmH[:] = p.Dm_H

            self.cc_cells.append(cH_cells)
            self.cc_ecm.append(cH_ecm)
            self.zs.append(p.z_H)
            self.Dm_mems.append(DmH)
            self.D_free.append(p.Do_H)

        # Initialize membrane thickness:
        self.tm = np.zeros(len(cells.cell_i))
        self.tm[:] = p.tm

        # Initialize environmental volume:
        self.envV = np.zeros(len(cells.cell_i))
        self.envV[:] = p.vol_env

        flx = np.zeros(len(cells.gj_i))
        self.fluxes_gj = [flx,flx,flx,flx,flx,flx,flx]   # stores gj fluxes for each ion
        self.gjopen_time = []   # stores gj open fraction at each time
        self.fgj_time = []      # stores the gj fluxes for each ion at each time
        self.Igj =[]            # current for each gj
        self.Igj_time = []      # current for each gj at each time

        # self.Dm_er = copy.deepcopy(self.Dm_cells)

        self.vm_to_cells = np.zeros(len(cells.cell_i))
        self.vm_to_ecm = np.zeros(len(cells.ecm_i))
        self.v_er = np.zeros(len(cells.cell_i))

        self.NaKATP_block = np.ones(len(cells.mem_i))  # initialize NaKATP blocking vector
        self.HKATP_block = np.ones(len(cells.mem_i))  # initialize HKATP blocking vector
        self.CaATP_block = np.ones(len(cells.mem_i))  # initialize CaATP (plasm membrane) blocking vector
        self.CaER_block = np.ones(len(cells.cell_i)) # initialize CaATP (ER membrane) blocking vector
        self.VATP_block = np.ones(len(cells.mem_i)) # initialize VATP (plasma membrane) blocking vector

        self.cc_er = np.asarray(self.cc_er)   # initialize endoplasmic reticulum concentration array
        self.cc_er_to = np.copy(self.cc_er)

        self.cDye_cell = np.zeros(len(cells.cell_i))   # initialize voltage sensitive dye array for cell and env't

        self.cDye_ecm = np.zeros(len(cells.ecm_i))
        self.cDye_ecm[:] = p.cDye_to

        # add channel noise to the model:
        self.channel_noise_factor = np.random.random(len(cells.mem_i))
        self.Dm_mems[self.iK] = (p.channel_noise_level*self.channel_noise_factor + 1)*self.Dm_mems[self.iK]

        # add a random walk on protein concentration to generate dynamic noise:
        self.protein_noise_factor = p.dynamic_noise_level*(np.random.random(len(cells.cell_i)) - 0.5)

        if p.ions_dict['P']==1:
            self.cc_cells[self.iP] = self.cc_cells[self.iP]*(1+ self.protein_noise_factor)

        loggers.log_info('This world contains '+ str(cells.cell_number) + ' cells.')
        loggers.log_info('Each cell has an average of ' + str(round(cells.average_nn,2)) + ' nearest-neighbours.')
        loggers.log_info('You are running the ion profile: ' + p.ion_profile)
        loggers.log_info('Ions in this simulation: ' + str(self.ionlabel))
        loggers.log_info('If you have selected features that use omitted ions they will be ignored.')

    def tissueInit(self,cells,p):
        """
        Runs initializations for all user-specified options that will be used in the main simulation.

        """

        # add channel noise to the model:
        self.Dm_cells[self.iK] = (p.channel_noise_level*self.channel_noise_factor + 1)*self.Dm_cells[self.iK]

        # Initialize an array structure that will hold user-scheduled changes to membrane permeabilities:
        Dm_cellsA = np.asarray(self.Dm_cells)
        Dm_cellsER = np.asarray(self.Dm_er)

        self.Dm_base = np.copy(Dm_cellsA) # make a copy that will serve as the unaffected values base

        self.Dm_scheduled = np.copy(Dm_cellsA)
        self.Dm_scheduled[:] = 0

        # Initialize an array structure that will hold dynamic voltage-gated channel changes to mem permeability:
        self.Dm_vg = np.copy(Dm_cellsA)
        self.Dm_vg[:] = 0

        # Initialize an array structure that will hold dynamic calcium-gated channel changes to mem perms:
        self.Dm_cag = np.copy(Dm_cellsA)
        self.Dm_cag[:] = 0

        self.Dm_er_base = np.copy(Dm_cellsER)

        self.Dm_er_CICR = np.copy(Dm_cellsER)
        self.Dm_er_CICR[:] = 0

        self.gj_block = np.ones(len(cells.gj_i))   # initialize the gap junction blocking vector to ones

        self.cIP3 = np.zeros(len(cells.cell_i))  # initialize a vector to hold IP3 concentrations
        self.cIP3[:] = p.cIP3_to                 # set the initial concentration of IP3 from params file

        self.cIP3_env = np.zeros(len(cells.cell_i))     # initialize IP3 concentration of the environment
        self.cIP3_env[:] = p.cIP3_to_env

        # Parameter assignments for all existing options
        if p.Ca_dyn_options['CICR'] != 0:

            self.stateER = np.zeros(len(cells.cell_i))   # state of ER membrane Ca permeability

            self.maxDmCaER = p.Ca_dyn_options['CICR'][0][0]
            self.topCa = p.Ca_dyn_options['CICR'][0][1]
            self.bottomCa =  p.Ca_dyn_options['CICR'][0][2]

            if len(p.Ca_dyn_options['CICR'][1])==0:
                pass

            else:
                self.midCaR = p.Ca_dyn_options['CICR'][1][0]
                self.widthCaR = p.Ca_dyn_options['CICR'][1][1]

            if len(p.Ca_dyn_options['CICR'][2])==0:
                pass

            else:
                self.KhmIP3 = p.Ca_dyn_options['CICR'][2][0]
                self.n_IP3 = p.Ca_dyn_options['CICR'][2][1]

        if p.scheduled_options['Na_mem'] != 0:

            self.t_on_Namem = p.scheduled_options['Na_mem'][0]
            self.t_off_Namem = p.scheduled_options['Na_mem'][1]
            self.t_change_Namem = p.scheduled_options['Na_mem'][2]
            self.mem_mult_Namem = p.scheduled_options['Na_mem'][3]

        if p.scheduled_options['K_mem'] != 0:

            self.t_on_Kmem = p.scheduled_options['K_mem'][0]
            self.t_off_Kmem = p.scheduled_options['K_mem'][1]
            self.t_change_Kmem = p.scheduled_options['K_mem'][2]
            self.mem_mult_Kmem = p.scheduled_options['K_mem'][3]

        if p.scheduled_options['Cl_mem'] != 0:

            self.t_on_Clmem = p.scheduled_options['Cl_mem'][0]
            self.t_off_Clmem = p.scheduled_options['Cl_mem'][1]
            self.t_change_Clmem = p.scheduled_options['Cl_mem'][2]
            self.mem_mult_Clmem = p.scheduled_options['Cl_mem'][3]


        if p.scheduled_options['Ca_mem'] != 0:

            self.t_on_Camem = p.scheduled_options['Ca_mem'][0]
            self.t_off_Camem = p.scheduled_options['Ca_mem'][1]
            self.t_change_Camem = p.scheduled_options['Ca_mem'][2]
            self.mem_mult_Camem = p.scheduled_options['Ca_mem'][3]


        if p.scheduled_options['IP3'] != 0:

            self.t_onIP3 = p.scheduled_options['IP3'][0]
            self.t_offIP3 = p.scheduled_options['IP3'][1]
            self.t_changeIP3 = p.scheduled_options['IP3'][2]
            self.rate_IP3 = p.scheduled_options['IP3'][3]


        if p.global_options['K_env'] != 0:

            self.t_on_Kenv = p.global_options['K_env'][0]
            self.t_off_Kenv = p.global_options['K_env'][1]
            self.t_change_Kenv = p.global_options['K_env'][2]
            self.mem_mult_Kenv = p.global_options['K_env'][3]


        if p.global_options['Cl_env'] != 0:

            self.t_on_Clenv = p.global_options['Cl_env'][0]
            self.t_off_Clenv = p.global_options['Cl_env'][1]
            self.t_change_Clenv = p.global_options['Cl_env'][2]
            self.mem_mult_Clenv = p.global_options['Cl_env'][3]


        if p.global_options['Na_env'] != 0:

            self.t_on_Naenv = p.global_options['Na_env'][0]
            self.t_off_Naenv = p.global_options['Na_env'][1]
            self.t_change_Naenv = p.global_options['Na_env'][2]
            self.mem_mult_Naenv = p.global_options['Na_env'][3]


        if p.global_options['T_change'] != 0:

            self.tonT = p.global_options['T_change'][0]
            self.toffT = p.global_options['T_change'][1]
            self.trampT = p.global_options['T_change'][2]
            self.multT = p.global_options['T_change'][3]


        if p.global_options['gj_block'] != 0:

            self.tonGJ = p.global_options['gj_block'][0]
            self.toffGJ = p.global_options['gj_block'][1]
            self.trampGJ = p.global_options['gj_block'][2]

        if p.global_options['NaKATP_block'] != 0:

            self.tonNK = p.global_options['NaKATP_block'][0]
            self.toffNK = p.global_options['NaKATP_block'][1]
            self.trampNK = p.global_options['NaKATP_block'][2]

        if p.global_options['HKATP_block'] != 0:

            self.tonHK = p.global_options['HKATP_block'][0]
            self.toffHK = p.global_options['HKATP_block'][1]
            self.trampHK = p.global_options['HKATP_block'][2]


        if p.vg_options['Na_vg'] != 0:

            # Initialization of logic values for voltage gated sodium channel
            self.maxDmNa = p.vg_options['Na_vg'][0]
            self.v_activate_Na = p.vg_options['Na_vg'][1]
            self.v_inactivate_Na = p.vg_options['Na_vg'][2]
            self.v_deactivate_Na = p.vg_options['Na_vg'][3]
            self.t_alive_Na = p.vg_options['Na_vg'][4]
            self.t_dead_Na = p.vg_options['Na_vg'][5]

            # Initialize matrices defining states of vgNa channels for each cell:
            self.inactivated_Na = np.zeros(len(cells.cell_i))
            self.vgNa_state = np.zeros(len(cells.cell_i))

            self.vgNa_aliveTimer = np.zeros(len(cells.cell_i)) # sim time at which vgNa starts to close if activated
            self.vgNa_deadTimer = np.zeros(len(cells.cell_i)) # sim time at which vgNa reactivates after inactivation

        if p.vg_options['K_vg'] !=0:

            # Initialization of logic values forr voltage gated potassium channel
            self.maxDmK = p.vg_options['K_vg'][0]
            self.v_on_K = p.vg_options['K_vg'][1]
            self.v_off_K = p.vg_options['K_vg'][2]
            self.t_alive_K = p.vg_options['K_vg'][3]

            # Initialize matrices defining states of vgK channels for each cell:
            self.active_K = np.zeros(len(cells.cell_i))
            self.crossed_activate_K = np.zeros(len(cells.cell_i))
            self.crossed_inactivate_K = np.zeros(len(cells.cell_i))

            # Initialize other matrices for vgK timing logic: NEW!
            self.vgK_state = np.zeros(len(cells.cell_i))   # state can be 0 = off, 1 = open
            self.vgK_OFFtime = np.zeros(len(cells.cell_i)) # sim time at which vgK starts to close


        if p.vg_options['Ca_vg'] !=0:

            # Initialization of logic values forr voltage gated potassium channel
            self.maxDmCa = p.vg_options['Ca_vg'][0]
            self.v_on_Ca = p.vg_options['Ca_vg'][1]
            self.v_off_Ca = p.vg_options['Ca_vg'][2]
            self.ca_upper_ca = p.vg_options['Ca_vg'][3]
            self.ca_lower_ca = p.vg_options['Ca_vg'][4]

            # Initialize matrices defining states of vgK channels for each cell:
            self.active_Ca = np.zeros(len(cells.cell_i))

            self.vgCa_state = np.zeros(len(cells.cell_i))   # state can be 0 = off, 1 = open

        if p.vg_options['K_cag'] != 0:

            self.maxDmKcag = p.vg_options['K_cag'][0]
            self.Kcag_halfmax = p.vg_options['K_cag'][1]
            self.Kcag_n = p.vg_options['K_cag'][2]

            # Initialize matrices defining states of cag K channels for each cell:
            self.active_Kcag = np.zeros(len(cells.cell_i))

        # Initialize target cell sets for dynamically gated channels from user options:
        if p.gated_targets == 'none':
            self.target_cells = np.zeros(len(cells.cell_i))

        elif p.gated_targets == 'all':
            self.target_cells = np.ones(len(cells.cell_i))

        elif p.gated_targets == 'random1':
            shuffle(cells.cell_i)
            trgt = cells.cell_i[0]
            self.target_cells = np.zeros(len(cells.cell_i))
            self.target_cells[trgt] = 1

        elif p.gated_targets == 'random50':
            self.target_cells = np.random.random(len(cells.cell_i))
            self.target_cells = np.rint(self.target_cells)

        elif isinstance(p.gated_targets,list):
            self.target_cells = np.zeros(len(cells.cell_i))
            self.target_cells[p.gated_targets] = 1

        # allow for option to independently schedule an intervention to cells distinct from voltage gated:
        if p.scheduled_targets == 'none':
            self.scheduled_target_inds = []

        elif p.scheduled_targets == 'all':
            self.scheduled_target_inds = cells.cell_i

        elif p.scheduled_targets =='random1':
            shuffle(cells.cell_i)
            trgt2 = cells.cell_i[0]
            self.scheduled_target_inds = [trgt2]

        elif p.scheduled_targets == 'random50':
            shuffle(cells.cell_i)
            halflength = int(len(cells.cell_i)/2)
            self.scheduled_target_inds = [cells.cell_i[x] for x in range(0,halflength)]

        elif isinstance(p.scheduled_targets, list):
            # self.scheduled_target_inds = np.zeros(len(cells.cell_i))
            # self.scheduled_target_inds[p.scheduled_targets] = 1
            self.scheduled_target_inds = p.scheduled_targets



        loggers.log_info('This world contains '+ str(cells.cell_number) + ' cells.')
        loggers.log_info('Each cell has an average of '+ str(round(cells.average_nn,2)) + ' nearest-neighbours.')
        loggers.log_info('You are running the ion profile: '+ p.ion_profile)


        loggers.log_info('Ions in this simulation: ' + str(self.ionlabel))
        loggers.log_info('If you have selected features using other ions, they will be ignored.')

    def tissueInit_ECM(self,cells,p):

        # add channel noise to the model:
        self.Dm_mems[self.iK] = (p.channel_noise_level*self.channel_noise_factor + 1)*self.Dm_mems[self.iK]

        # Initialize an array structure that will hold user-scheduled changes to membrane permeabilities:
        Dm_memsA = np.asarray(self.Dm_mems)
        Dm_cellsER = np.asarray(self.Dm_er)

        self.Dm_base = np.copy(Dm_memsA) # make a copy that will serve as the unaffected values base

        self.Dm_scheduled = np.copy(Dm_memsA)
        self.Dm_scheduled[:] = 0

        # Initialize an array structure that will hold dynamic voltage-gated channel changes to mem permeability:
        self.Dm_vg = np.copy(Dm_memsA)
        self.Dm_vg[:] = 0

        # Initialize an array structure that will hold dynamic calcium-gated channel changes to mem perms:
        self.Dm_cag = np.copy(Dm_memsA)
        self.Dm_cag[:] = 0

        self.Dm_er_base = np.copy(Dm_cellsER)

        self.Dm_er_CICR = np.copy(Dm_cellsER)
        self.Dm_er_CICR[:] = 0

        self.gj_block = np.ones(len(cells.gj_i))   # initialize the gap junction blocking vector to ones

        self.cIP3 = np.zeros(len(cells.cell_i))  # initialize a vector to hold IP3 concentrations
        self.cIP3[:] = p.cIP3_to                 # set the initial concentration of IP3 from params file

        self.cIP3_ecm = np.zeros(len(cells.ecm_i))     # initialize IP3 concentration of the environment
        self.cIP3_ecm[:] = p.cIP3_to_env

        # Parameter assignments for all existing options
        if p.Ca_dyn_options['CICR'] != 0:

            self.stateER = np.zeros(len(cells.cell_i))   # state of ER membrane Ca permeability

            self.maxDmCaER = p.Ca_dyn_options['CICR'][0][0]
            self.topCa = p.Ca_dyn_options['CICR'][0][1]
            self.bottomCa =  p.Ca_dyn_options['CICR'][0][2]

            if len(p.Ca_dyn_options['CICR'][1])==0:
                pass

            else:
                self.midCaR = p.Ca_dyn_options['CICR'][1][0]
                self.widthCaR = p.Ca_dyn_options['CICR'][1][1]

            if len(p.Ca_dyn_options['CICR'][2])==0:
                pass

            else:
                self.KhmIP3 = p.Ca_dyn_options['CICR'][2][0]
                self.n_IP3 = p.Ca_dyn_options['CICR'][2][1]

        if p.scheduled_options['Na_mem'] != 0:

            self.t_on_Namem = p.scheduled_options['Na_mem'][0]
            self.t_off_Namem = p.scheduled_options['Na_mem'][1]
            self.t_change_Namem = p.scheduled_options['Na_mem'][2]
            self.mem_mult_Namem = p.scheduled_options['Na_mem'][3]

        if p.scheduled_options['K_mem'] != 0:

            self.t_on_Kmem = p.scheduled_options['K_mem'][0]
            self.t_off_Kmem = p.scheduled_options['K_mem'][1]
            self.t_change_Kmem = p.scheduled_options['K_mem'][2]
            self.mem_mult_Kmem = p.scheduled_options['K_mem'][3]

        if p.scheduled_options['Cl_mem'] != 0:

            self.t_on_Clmem = p.scheduled_options['Cl_mem'][0]
            self.t_off_Clmem = p.scheduled_options['Cl_mem'][1]
            self.t_change_Clmem = p.scheduled_options['Cl_mem'][2]
            self.mem_mult_Clmem = p.scheduled_options['Cl_mem'][3]


        if p.scheduled_options['Ca_mem'] != 0:

            self.t_on_Camem = p.scheduled_options['Ca_mem'][0]
            self.t_off_Camem = p.scheduled_options['Ca_mem'][1]
            self.t_change_Camem = p.scheduled_options['Ca_mem'][2]
            self.mem_mult_Camem = p.scheduled_options['Ca_mem'][3]


        if p.scheduled_options['IP3'] != 0:

            self.t_onIP3 = p.scheduled_options['IP3'][0]
            self.t_offIP3 = p.scheduled_options['IP3'][1]
            self.t_changeIP3 = p.scheduled_options['IP3'][2]
            self.rate_IP3 = p.scheduled_options['IP3'][3]


        if p.global_options['K_env'] != 0:

            self.t_on_Kenv = p.global_options['K_env'][0]
            self.t_off_Kenv = p.global_options['K_env'][1]
            self.t_change_Kenv = p.global_options['K_env'][2]
            self.mem_mult_Kenv = p.global_options['K_env'][3]


        if p.global_options['Cl_env'] != 0:

            self.t_on_Clenv = p.global_options['Cl_env'][0]
            self.t_off_Clenv = p.global_options['Cl_env'][1]
            self.t_change_Clenv = p.global_options['Cl_env'][2]
            self.mem_mult_Clenv = p.global_options['Cl_env'][3]


        if p.global_options['Na_env'] != 0:

            self.t_on_Naenv = p.global_options['Na_env'][0]
            self.t_off_Naenv = p.global_options['Na_env'][1]
            self.t_change_Naenv = p.global_options['Na_env'][2]
            self.mem_mult_Naenv = p.global_options['Na_env'][3]


        if p.global_options['T_change'] != 0:

            self.tonT = p.global_options['T_change'][0]
            self.toffT = p.global_options['T_change'][1]
            self.trampT = p.global_options['T_change'][2]
            self.multT = p.global_options['T_change'][3]


        if p.global_options['gj_block'] != 0:

            self.tonGJ = p.global_options['gj_block'][0]
            self.toffGJ = p.global_options['gj_block'][1]
            self.trampGJ = p.global_options['gj_block'][2]

        if p.global_options['NaKATP_block'] != 0:

            self.tonNK = p.global_options['NaKATP_block'][0]
            self.toffNK = p.global_options['NaKATP_block'][1]
            self.trampNK = p.global_options['NaKATP_block'][2]

        if p.global_options['HKATP_block'] != 0:

            self.tonHK = p.global_options['HKATP_block'][0]
            self.toffHK = p.global_options['HKATP_block'][1]
            self.trampHK = p.global_options['HKATP_block'][2]


        if p.vg_options['Na_vg'] != 0:

            # Initialization of logic values for voltage gated sodium channel
            self.maxDmNa = p.vg_options['Na_vg'][0]
            self.v_activate_Na = p.vg_options['Na_vg'][1]
            self.v_inactivate_Na = p.vg_options['Na_vg'][2]
            self.v_deactivate_Na = p.vg_options['Na_vg'][3]
            self.t_alive_Na = p.vg_options['Na_vg'][4]
            self.t_dead_Na = p.vg_options['Na_vg'][5]

            # Initialize matrices defining states of vgNa channels for each cell membrane:
            self.inactivated_Na = np.zeros(len(cells.mem_i))
            self.vgNa_state = np.zeros(len(cells.mem_i))

            self.vgNa_aliveTimer = np.zeros(len(cells.mem_i)) # sim time at which vgNa starts to close if activated
            self.vgNa_deadTimer = np.zeros(len(cells.mem_i)) # sim time at which vgNa reactivates after inactivation

        if p.vg_options['K_vg'] !=0:

            # Initialization of logic values forr voltage gated potassium channel
            self.maxDmK = p.vg_options['K_vg'][0]
            self.v_on_K = p.vg_options['K_vg'][1]
            self.v_off_K = p.vg_options['K_vg'][2]
            self.t_alive_K = p.vg_options['K_vg'][3]

            # Initialize matrices defining states of vgK channels for each cell:
            self.active_K = np.zeros(len(cells.mem_i))
            self.crossed_activate_K = np.zeros(len(cells.mem_i))
            self.crossed_inactivate_K = np.zeros(len(cells.mem_i))

            # Initialize other matrices for vgK timing logic: NEW!
            self.vgK_state = np.zeros(len(cells.mem_i))   # state can be 0 = off, 1 = open
            self.vgK_OFFtime = np.zeros(len(cells.mem_i)) # sim time at which vgK starts to close


        if p.vg_options['Ca_vg'] !=0:

            # Initialization of logic values for voltage gated calcium channel
            self.maxDmCa = p.vg_options['Ca_vg'][0]
            self.v_on_Ca = p.vg_options['Ca_vg'][1]
            self.v_off_Ca = p.vg_options['Ca_vg'][2]
            self.ca_upper_ca = p.vg_options['Ca_vg'][3]
            self.ca_lower_ca = p.vg_options['Ca_vg'][4]

            # Initialize matrices defining states of vgK channels for each cell membrane:
            self.active_Ca = np.zeros(len(cells.mem_i))

            self.vgCa_state = np.zeros(len(cells.mem_i))   # state can be 0 = off, 1 = open

        if p.vg_options['K_cag'] != 0:

            self.maxDmKcag = p.vg_options['K_cag'][0]
            self.Kcag_halfmax = p.vg_options['K_cag'][1]
            self.Kcag_n = p.vg_options['K_cag'][2]

            # Initialize matrices defining states of cag K channels for each cell membrane:
            self.active_Kcag = np.zeros(len(cells.mem_i))

        # Initialize target cell sets for dynamically gated channels from user options:
        if p.gated_targets == 'none':
            self.target_cells = np.zeros(len(cells.cell_i))
            self.target_mems = np.zeros(len(cells.mem_i))

        elif p.gated_targets == 'all':
            self.target_cells = np.ones(len(cells.cell_i))
            self.target_mems = np.ones(len(cells.mem_i))

        elif p.gated_targets == 'random1':
            shuffle(cells.cell_i)
            trgt = cells.cell_i[0]
            self.target_cells = np.zeros(len(cells.cell_i))
            self.target_mems = np.zeros(len(cells.mem_i))

            self.target_cells[trgt] = 1

            target_mems_inds = cells.cell_to_mems[self.target_cells]
            target_mems_inds,_,_ = tb.flatten(self.target_mems_inds)

            self.target_mems[target_mems_inds] = 1

        elif p.gated_targets == 'random50':

            self.target_cells = np.random.random(len(cells.cell_i))
            self.target_cells = np.rint(self.target_cells)

            self.target_mems = np.zeros(len(cells.mem_i))

            target_mems_inds = cells.cell_to_mems[self.target_cells]
            target_mems_inds,_,_ = tb.flatten(target_mems_inds)

            self.target_mems[target_mems_inds] = 1

        elif isinstance(p.gated_targets,list):
            self.target_cells = np.zeros(len(cells.cell_i))
            self.target_cells[p.gated_targets] = 1

            self.target_mems = np.zeros(len(cells.mem_i))

            target_mems_inds = cells.cell_to_mems[self.target_cells]
            target_mems_inds,_,_ = tb.flatten(target_mems_inds)

            self.target_mems[target_mems_inds] = 1

        # allow for option to independently schedule an intervention to cells distinct from voltage gated:
        if p.scheduled_targets == 'none':
            self.scheduled_target_inds = []
            self.scheduled_target_mem_inds = []

        elif p.scheduled_targets == 'all':
            self.scheduled_target_inds = cells.cell_i
            self.scheduled_target_mem_inds = cells.mem_i

        elif p.scheduled_targets =='random1':
            shuffle(cells.cell_i)
            trgt2 = cells.cell_i[0]

            self.scheduled_target_inds = [trgt2]

            self.scheduled_target_mems_inds = cells.cell_to_mems[trgt2]
            self.scheduled_target_mems_inds,_,_ = tb.flatten(self.scheduled_target_mems_inds)

        elif p.scheduled_targets == 'random50':
            shuffle(cells.cell_i)
            halflength = int(len(cells.cell_i)/2)
            self.scheduled_target_inds = [cells.cell_i[x] for x in range(0,halflength)]

            trgt3 = self.scheduled_target_inds

            self.scheduled_target_mems_inds = cells.cell_to_mems[trgt3]
            self.scheduled_target_mems_inds,_,_ = tb.flatten(self.scheduled_target_mems_inds)


        elif isinstance(p.scheduled_targets, list):

            self.scheduled_target_inds = p.scheduled_targets

            trgt4 = self.scheduled_target_inds

            self.scheduled_target_mems_inds = cells.cell_to_mems[trgt4]
            self.scheduled_target_mems_inds,_,_ = tb.flatten(self.scheduled_target_mems_inds)


        loggers.log_info('This world contains '+ str(cells.cell_number) + ' cells.')
        loggers.log_info('Each cell has an average of '+ str(round(cells.average_nn,2)) + ' nearest-neighbours.')
        loggers.log_info('You are running the ion profile: '+ p.ion_profile)

        loggers.log_info('Ions in this simulation: ' + str(self.ionlabel))
        loggers.log_info('If you have selected features using other ions, they will be ignored.')

    def runInit(self,cells,p):
        """
        Runs an initialization simulation from the existing data state of the Simulation object,
        and saves the resulting data (including the cell world geometry) to files that can be read later.

        Parameters:
        -----------
        cells               An instance of the World class. This is required because simulation data is specific
                            to the cell world data so it needs the reference to save.
        timesteps           The number of timesteps over which the simulation is run. Note the size of the time
                            interval is specified as dt in the parameters module.

        """
        # Reinitialize data structures that hold time data
        self.cc_time = []  # data array holding the cell concentrations at time points
        self.cc_env_time = []  # data array holding the environmental concentrations at time points
        self.cc_er_time = []  # data array holding endoplasmic reticulum concentrations at time points

        self.vm_time = []  # data array holding voltage at time points
        self.cDye_time = [] # data array holding voltage-sensitive dye concentrations at time points

        self.time = []     # time values of the simulation

        tt = np.linspace(0,p.init_tsteps*p.dt,p.init_tsteps)

        i = 0 # resample the time vector to save data at specific times:
        tsamples =[]
        resample = p.t_resample
        while i < len(tt)-resample:
            i = i + resample
            tsamples.append(tt[i])
        tsamples = set(tsamples)

        # report
        loggers.log_info('Your initialization is running for ' + str(round((p.init_tsteps*p.dt)/60,2)) +
                         ' minutes of in-world time.')

        # get the initial net, unbalanced charge and voltage in each cell:
        q_cells = get_charge(self.cc_cells,self.zs,cells.cell_vol,p)
        self.vm = get_volt(q_cells,cells.cell_sa,p)

        do_once = True  # a variable to time the loop only once

        for t in tt:   # run through the loop

            if do_once == True:
                loop_measure = time.time()

            # run the Na-K-ATPase pump:
            self.cc_cells[self.iNa],self.cc_env[self.iNa],self.cc_cells[self.iK],self.cc_env[self.iK], _, _ =\
                pumpNaKATP(self.cc_cells[self.iNa],self.cc_env[self.iNa],self.cc_cells[self.iK],self.cc_env[self.iK],
                    cells.cell_sa,cells.cell_vol,self.envV,self.vm,self.T,p,self.NaKATP_block)

            # recalculate the net, unbalanced charge and voltage in each cell:
            q_cells = get_charge(self.cc_cells,self.zs,cells.cell_vol,p)
            self.vm = get_volt(q_cells,cells.cell_sa,p)

            # if calcium is present, run the Ca-ATPase pump and fill up the endoplasmic reticulum:
            if  p.ions_dict['Ca'] == 1:

                self.cc_cells[self.iCa],self.cc_env[self.iCa], _ =\
                    pumpCaATP(self.cc_cells[self.iCa],self.cc_env[self.iCa],cells.cell_sa,cells.cell_vol,
                        self.envV,self.vm,self.T,p,self.CaATP_block)

                # recalculate the net, unbalanced charge and voltage in each cell:
                q_cells = get_charge(self.cc_cells,self.zs,cells.cell_vol,p)
                self.vm = get_volt(q_cells,cells.cell_sa,p)

                if p.Ca_dyn ==1:

                    self.cc_er[0],self.cc_cells[self.iCa], _ =\
                        pumpCaER(self.cc_er[0],self.cc_cells[self.iCa],cells.cell_sa,p.ER_vol*cells.cell_vol,
                            cells.cell_vol,self.v_er,self.T,p,self.CaER_block)

                    # recalculate the net, unbalanced charge and voltage in each cell:
                    q_cells = get_charge(self.cc_cells,self.zs,cells.cell_vol,p)
                    self.vm = get_volt(q_cells,cells.cell_sa,p)

                    q_er = get_charge(self.cc_er,self.z_er,p.ER_vol*cells.cell_vol,p)
                    v_er_o = get_volt(q_er,p.ER_sa*cells.cell_sa,p)
                    self.v_er = v_er_o - self.vm

            if p.ions_dict['H'] == 1:

                # electrofuse the H+ ion between the cytoplasm and the environment
                self.cc_env[self.iH],self.cc_cells[self.iH],f_H1 = \
                    electrofuse(self.cc_env[self.iH],self.cc_cells[self.iH],self.Dm_cells[self.iH],self.tm,cells.cell_sa,
                        self.envV,cells.cell_vol,self.zs[self.iH],self.vm,self.T,p)

                # buffer what's happening with H+ flux to or from the cell and environment:
                delH_cell = (f_H1*p.dt/cells.cell_vol)    # relative change in H wrt the cell
                delH_env =  -(f_H1*p.dt/p.vol_env)    # relative change in H wrt to environment

                self.cc_cells[self.iH], self.cc_cells[self.iM], self.cHM_cells = bicarbBuffer(
                    self.cc_cells[self.iH],self.cc_cells[self.iM],self.cHM_cells,delH_cell,p)

                self.cc_env[self.iH], self.cc_env[self.iM], self.cHM_env = bicarbBuffer(
                    self.cc_env[self.iH],self.cc_env[self.iM],self.cHM_env,delH_env,p)

                # recalculate the net, unbalanced charge and voltage in each cell:
                q_cells = get_charge(self.cc_cells,self.zs,cells.cell_vol,p)
                self.vm = get_volt(q_cells,cells.cell_sa,p)


                if p.HKATPase_dyn == 1:

                    # if HKATPase pump is desired, run the H-K-ATPase pump:
                    self.cc_cells[self.iH],self.cc_env[self.iH],self.cc_cells[self.iK],self.cc_env[self.iK], f_H2, _ =\
                    pumpHKATP(self.cc_cells[self.iH],self.cc_env[self.iH],self.cc_cells[self.iK],self.cc_env[self.iK],
                        cells.cell_sa,cells.cell_vol,self.envV,self.vm,self.T,p,self.HKATP_block)

                     # buffer what's happening with H+ flux to or from the cell and environment:
                    delH_cell = (f_H2*p.dt/cells.cell_vol)    # relative change in H wrt the cell
                    delH_env = -(f_H2*p.dt/p.vol_env)    # relative change in H wrt to environment

                    self.cc_cells[self.iH], self.cc_cells[self.iM], self.cHM_cells = bicarbBuffer(
                        self.cc_cells[self.iH],self.cc_cells[self.iM],self.cHM_cells,delH_cell,p)

                    self.cc_env[self.iH], self.cc_env[self.iM], self.cHM_env = bicarbBuffer(
                        self.cc_env[self.iH],self.cc_env[self.iM],self.cHM_env,delH_env,p)

                    # recalculate the net, unbalanced charge and voltage in each cell:
                    q_cells = get_charge(self.cc_cells,self.zs,cells.cell_vol,p)
                    self.vm = get_volt(q_cells,cells.cell_sa,p)

                if p.VATPase_dyn == 1:

                     # if HKATPase pump is desired, run the H-K-ATPase pump:
                    self.cc_cells[self.iH],self.cc_env[self.iH], f_H3 =\
                    pumpVATP(self.cc_cells[self.iH],self.cc_env[self.iH],
                        cells.cell_sa,cells.cell_vol,self.envV,self.vm,self.T,p,self.VATP_block)

                     # buffer what's happening with H+ flux to or from the cell and environment:
                    delH_cell = (f_H3*p.dt/cells.cell_vol)    # relative change in H wrt the cell
                    delH_env = -(f_H3*p.dt/p.vol_env)    # relative change in H wrt to environment

                    self.cc_cells[self.iH], self.cc_cells[self.iM], self.cHM_cells = bicarbBuffer(
                        self.cc_cells[self.iH],self.cc_cells[self.iM],self.cHM_cells,delH_cell,p)

                    self.cc_env[self.iH], self.cc_env[self.iM], self.cHM_env = bicarbBuffer(
                        self.cc_env[self.iH],self.cc_env[self.iM],self.cHM_env,delH_env,p)

                    # recalculate the net, unbalanced charge and voltage in each cell:
                    q_cells = get_charge(self.cc_cells,self.zs,cells.cell_vol,p)
                    self.vm = get_volt(q_cells,cells.cell_sa,p)

            # electro-diffuse all ions (except for proteins, which don't move!) across the cell membrane:
            shuffle(self.movingIons)  # shuffle the ion indices so it's not the same order every time step

            for i in self.movingIons:

                self.cc_env[i],self.cc_cells[i],_ = \
                    electrofuse(self.cc_env[i],self.cc_cells[i],self.Dm_cells[i],self.tm,cells.cell_sa,
                        self.envV,cells.cell_vol,self.zs[i],self.vm,self.T,p)

                # recalculate the net, unbalanced charge and voltage in each cell:
                q_cells = get_charge(self.cc_cells,self.zs,cells.cell_vol,p)
                self.vm = get_volt(q_cells,cells.cell_sa,p)

            if p.Ca_dyn == 1 and p.ions_dict['Ca'] == 1:
                    # electrodiffusion of ions between cell and endoplasmic reticulum
                    # First do calciu:m
                    self.cc_cells[self.iCa],self.cc_er[0],_ = \
                    electrofuse(self.cc_cells[self.iCa],self.cc_er[0],self.Dm_er[0],self.tm,p.ER_sa*cells.cell_sa,
                        cells.cell_vol,p.ER_vol*cells.cell_vol,self.z_er[0],self.v_er,self.T,p)

                    # next do charge compensation anion:
                    self.cc_cells[self.iM],self.cc_er[1],_ = \
                    electrofuse(self.cc_cells[self.iM],self.cc_er[1],self.Dm_er[1],self.tm,p.ER_sa*cells.cell_sa,
                        cells.cell_vol,p.ER_vol*cells.cell_vol,self.z_er[1],self.v_er,self.T,p)

                    # recalculate the net, unbalanced charge and voltage in each cell:
                    q_cells = get_charge(self.cc_cells,self.zs,cells.cell_vol,p)
                    self.vm = get_volt(q_cells,cells.cell_sa,p)

                    q_er = get_charge(self.cc_er,self.z_er,p.ER_vol*cells.cell_vol,p)
                    v_er_o = get_volt(q_er,p.ER_sa*cells.cell_sa,p)
                    self.v_er = v_er_o - self.vm

            # if p.voltage_dye=1 electrodiffuse voltage sensitive dye between cell and environment
            if p.voltage_dye ==1:

                self.cDye_env,self.cDye_cell,_ = \
                        electrofuse(self.cDye_env,self.cDye_cell,p.Dm_Dye*self.id_cells,self.tm,cells.cell_sa,
                            self.envV,cells.cell_vol,p.z_Dye,self.vm,self.T,p)

            if p.dynamic_noise == 1 and p.ions_dict['P']==1:
                # add a random walk on protein concentration to generate dynamic noise:
                self.protein_noise_factor = p.dynamic_noise_level*(np.random.random(len(cells.cell_i)) - 0.5)
                self.cc_cells[self.iP] = self.cc_cells[self.iP]*(1+ self.protein_noise_factor)

                # recalculate the net, unbalanced charge and voltage in each cell:
                q_cells = get_charge(self.cc_cells,self.zs,cells.cell_vol,p)
                self.vm = get_volt(q_cells,cells.cell_sa,p)

            # check and ensure simulation stability
            check_v(self.vm)

            if t in tsamples:
                # add the new concentration and voltage data to the time-storage matrices:
                concs = copy.deepcopy(self.cc_cells)
                self.cc_time.append(concs)

                concs_env = copy.deepcopy(self.cc_env)
                self.cc_env_time.append(concs_env)

                vmm = copy.deepcopy(self.vm)
                self.vm_time.append(vmm)

                self.time.append(t)

                if p.Ca_dyn == 1 and p.ions_dict['Ca'] == 1:
                    concs_er = copy.deepcopy(self.cc_er)
                    self.cc_er_time.append(concs_er)

                if p.voltage_dye == 1:
                    cvdye = copy.deepcopy(self.cDye_cell)
                    self.cDye_time.append(cvdye)

                if p.plot_while_solving == True:
                    pass

            # get time for loop and estimate total time for simulation
            if do_once == True:
                loop_time = time.time() - loop_measure
                time_estimate = round(loop_time*p.init_tsteps,2)
                loggers.log_info("This run should take approximately " + str(time_estimate) + ' s to compute...')
                do_once = False

        celf = copy.deepcopy(self)

        datadump = [celf,cells,p]
        fh.saveSim(self.savedInit,datadump)
        message_1 = 'Initialization run saved to' + ' ' + p.init_path
        loggers.log_info(message_1)

        self.vm_to = copy.deepcopy(self.vm)

        for i in range(0,len(self.ionlabel)):
            endconc = np.round(np.mean(self.cc_time[-1][i]),6)
            label = self.ionlabel[i]
            concmess = 'Final average cytoplasmic concentration of'+ ' '+ label + ': '
            loggers.log_info(concmess + str(endconc) + ' mmol/L')

        for i in range(0,len(self.ionlabel)):
            endconc = np.round(np.mean(self.cc_env_time[-1][i]),6)
            label = self.ionlabel[i]
            concmess = 'Final environmental concentration of'+ ' '+ label + ': '
            loggers.log_info(concmess + str(endconc) + ' mmol/L')

        final_vmean = 1000*np.round(np.mean(self.vm_time[-1]),4)
        vmess = 'Final average cell Vmem of' + ': '
        loggers.log_info(vmess + str(final_vmean) + ' mV')

        if p.ions_dict['H'] == 1:
            final_pH = -np.log10(np.mean((self.cc_time[-1][self.iH])/1000))
            loggers.log_info('Final average cell pH ' + str(np.round(final_pH,2)))

            final_pH_env = -np.log10(np.mean((self.cc_env_time[-1][self.iH])/1000))
            loggers.log_info('Final environmental pH '+ str(np.round(final_pH_env,2)))

        if p.Ca_dyn == 1 and p.ions_dict['Ca']==1:

                endconc_er = np.round(np.mean(self.cc_er[0]),6)
                label = self.ionlabel[self.iCa]
                concmess = 'Final ER concentration of'+ ' '+ label + ': '
                loggers.log_info(concmess + str(endconc_er) + ' mmol/L')

        if p.voltage_dye == 1:

            dye_env_final = np.mean(self.cDye_env)
            dye_cell_final = np.mean(self.cDye_cell)
            loggers.log_info('Final dye concentration in the environment: '+ str(np.round(dye_env_final,6)) + ' mmol/L')
            loggers.log_info('Final average dye concentration in cells: '+ str(np.round(dye_cell_final,6)) + ' mmol/L')

    def runInit_ECM(self,cells,p):
        """
        Runs an initialization simulation from the existing data state of the Simulation object,
        and saves the resulting data (including the cell world geometry) to files that can be read later.

        Parameters:
        -----------
        cells               An instance of the World class. This is required because simulation data is specific
                            to the cell world data so it needs the reference to save.
        timesteps           The number of timesteps over which the simulation is run. Note the size of the time
                            interval is specified as dt in the parameters module.

        """
        # Reinitialize data structures that hold time data
        self.cc_time = []  # data array holding the cell concentrations at time points
        self.cc_ecm_time = []  # data array holding the environmental concentrations at time points
        self.cc_er_time = []  # data array holding endoplasmic reticulum concentrations at time points

        self.vm_time = []  # data array holding voltage across membrane at time points
        self.vcells_time = []  # data array holding intracellular voltages (at membranes) at time points
        self.vecm_time = []    # data array holding ecm voltages (at membrane) at time points
        self.cDye_time = [] # data array holding voltage-sensitive dye concentrations at time points

        self.time = []     # time values of the simulation

        tt = np.linspace(0,p.init_tsteps*p.dt,p.init_tsteps)   # timstep vector

        i = 0 # resample the time vector to save data at specific times:
        tsamples =[]
        resample = p.t_resample
        while i < len(tt)-resample:
            i = i + resample
            tsamples.append(tt[i])
        tsamples = set(tsamples)

        # report
        loggers.log_info('Your initialization (with extracellular spaces) is running for '
                         + str(round((p.init_tsteps*p.dt)/60,2)) +
                         ' minutes of in-world time.')

        # get the initial net, unbalanced charge and voltage in each cell and ecm space:
        self.update_V_ecm(cells,p)

        do_once = True  # a variable to time the loop only once

        for t in tt:   # run through the loop

            if do_once == True:
                loop_measure = time.time()

            # run the Na-K-ATPase pump:
            _,_,_,_, f_Na, f_K =\
                pumpNaKATP(self.cc_cells[self.iNa][cells.mem_to_cells],self.cc_ecm[self.iNa][cells.mem_to_ecm],
                    self.cc_cells[self.iK][cells.mem_to_cells],self.cc_ecm[self.iK][cells.mem_to_ecm],
                    cells.mem_sa,cells.cell_vol[cells.mem_to_cells], cells.ecm_vol[cells.mem_to_ecm],
                    self.vm,self.T,p,self.NaKATP_block)

            # update the cell and extracellular concentrations as a result of pump fluxes:
            self.update_C_ecm(self.iNa,f_Na,cells,p)
            self.update_C_ecm(self.iK,f_K,cells,p)

            # recalculate the net, unbalanced charge and voltage in each cell:
            self.update_V_ecm(cells,p)

            # if calcium is present, run the Ca-ATPase pump and fill up the endoplasmic reticulum:
            if  p.ions_dict['Ca'] == 1:

                _,_,f_Ca_ATP = pumpCaATP(self.cc_cells[self.iCa][cells.mem_to_cells],
                    self.cc_ecm[self.iCa][cells.mem_to_ecm],cells.mem_sa, cells.cell_vol[cells.mem_to_cells],
                    cells.ecm_vol[cells.mem_to_ecm],self.vm,self.T,p,self.CaATP_block)

                # update concentration of calcium in cells and ecm from pump activity:
                self.update_C_ecm(self.iCa,f_Ca_ATP,cells,p)

                # recalculate the net, unbalanced charge and voltage in each cell:
                self.update_V_ecm(cells,p)

                if p.Ca_dyn ==1:

                    self.cc_er[0],self.cc_cells[self.iCa], f_Ca_er =\
                        pumpCaER(self.cc_er[0],self.cc_cells[self.iCa],cells.cell_sa,p.ER_vol*cells.cell_vol,
                            cells.cell_vol,self.v_er,self.T,p,self.CaER_block)

                    # recalculate the net, unbalanced charge and voltage in each cell:
                    self.update_V_ecm(cells,p)

                    q_er = get_charge(self.cc_er,self.z_er,p.ER_vol*cells.cell_vol,p)
                    self.v_er = get_volt(q_er,p.ER_sa*cells.cell_sa,p)  # calculate voltage across er membrane

            if p.ions_dict['H'] == 1:

                self.Hplus_electrofuse_ecm(cells,p)

                if p.HKATPase_dyn == 1:

                    self.Hplus_HKATP_ecm(cells,p)

                if p.VATPase_dyn == 1:

                    self.Hplus_VATP_ecm(cells,p)

            # electro-diffuse all ions (except for proteins, which don't move!) across the cell membrane:
            shuffle(self.movingIons)  # shuffle the ion indices so it's not the same order every time step

            for i in self.movingIons:

                _,_,flux_ED = \
                    electrofuse(self.cc_ecm[i][cells.mem_to_ecm],self.cc_cells[i][cells.mem_to_cells],
                        self.Dm_mems[i],self.tm[cells.mem_to_cells],cells.mem_sa,
                        cells.ecm_vol[cells.mem_to_ecm],cells.cell_vol[cells.mem_to_cells],
                        self.zs[i],self.vm,self.T,p)

                # update the concentrations in the cell and ecm due to ED fluxes at membrane:
                self.update_C_ecm(i,flux_ED,cells,p)

                # recalculate the net, unbalanced charge and voltage in each cell:
                self.update_V_ecm(cells,p)

            if p.Ca_dyn == 1 and p.ions_dict['Ca'] == 1:

                    # electrodiffusion of ions between cell and endoplasmic reticulum

                    # First do calciu:m
                    self.cc_cells[self.iCa],self.cc_er[0],_ = \
                    electrofuse(self.cc_cells[self.iCa],self.cc_er[0],self.Dm_er[0],self.tm,p.ER_sa*cells.cell_sa,
                        cells.cell_vol,p.ER_vol*cells.cell_vol,self.z_er[0],self.v_er,self.T,p)

                    # next do charge compensation anion:
                    self.cc_cells[self.iM],self.cc_er[1],_ = \
                    electrofuse(self.cc_cells[self.iM],self.cc_er[1],self.Dm_er[1],self.tm,p.ER_sa*cells.cell_sa,
                        cells.cell_vol,p.ER_vol*cells.cell_vol,self.z_er[1],self.v_er,self.T,p)

                    # recalculate the net, unbalanced charge and voltage in each cell:
                    self.update_V_ecm(cells,p)

                    q_er = get_charge(self.cc_er,self.z_er,p.ER_vol*cells.cell_vol,p)
                    self.v_er = get_volt(q_er,p.ER_sa*cells.cell_sa,p)

            # if p.voltage_dye=1 electrodiffuse voltage sensitive dye between cell and environment
            if p.voltage_dye ==1:

                _,_, flux_dye = \
                        electrofuse(self.cDye_ecm[cells.mem_to_ecm],self.cDye_cell[cells.mem_to_cells],
                            p.Dm_Dye*self.id_cells[cells.mem_to_cells],self.tm[cells.mem_to_cells],
                            cells.mem_sa,cells.ecm_vol[cells.mem_to_ecm],cells.cell_vol[cells.mem_to_cells],
                            p.z_Dye,self.vm,self.T,p)

                 # update the dye concentrations in the cell and ecm due to ED fluxes at membrane
                self.cDye_cell = self.cDye_cell + \
                                    np.dot((flux_dye/cells.cell_vol[cells.mem_to_cells])*p.dt,cells.cell_UpdateMatrix)

                self.cDye_ecm = self.cDye_ecm - \
                                    np.dot((flux_dye/cells.ecm_vol[cells.mem_to_ecm])*p.dt,cells.ecm_UpdateMatrix)

            if p.dynamic_noise == 1 and p.ions_dict['P']==1:
                # add a random walk on protein concentration to generate dynamic noise:
                self.protein_noise_factor = p.dynamic_noise_level*(np.random.random(len(cells.cell_i)) - 0.5)
                self.cc_cells[self.iP] = self.cc_cells[self.iP]*(1+ self.protein_noise_factor)

                # recalculate the net, unbalanced charge and voltage in each cell:
                self.update_V_ecm(cells,p)

            # check and ensure simulation stability
            check_v(self.vm)

            if t in tsamples:
                # add the new concentration and voltage data to the time-storage matrices:
                concs = copy.deepcopy(self.cc_cells)
                self.cc_time.append(concs)

                concs_ecm = copy.deepcopy(self.cc_ecm)
                self.cc_ecm_time.append(concs_ecm)

                vmm = copy.deepcopy(self.vm)
                self.vm_time.append(vmm)

                self.time.append(t)

                if p.Ca_dyn == 1 and p.ions_dict['Ca'] == 1:
                    concs_er = copy.deepcopy(self.cc_er)
                    self.cc_er_time.append(concs_er)

                if p.voltage_dye == 1:
                    cvdye = copy.deepcopy(self.cDye_cell)
                    self.cDye_time.append(cvdye)

                if p.plot_while_solving == True:
                    pass

            # get time for loop and estimate total time for simulation
            if do_once == True:
                loop_time = time.time() - loop_measure
                time_estimate = round(loop_time*p.init_tsteps,2)
                loggers.log_info("This run should take approximately " + str(time_estimate) + ' s to compute...')
                do_once = False

        celf = copy.deepcopy(self)

        datadump = [celf,cells,p]
        fh.saveSim(self.savedInit,datadump)
        message_1 = 'Initialization run saved to' + ' ' + p.init_path
        loggers.log_info(message_1)

        self.vm_to = copy.deepcopy(self.vm)

        for i in range(0,len(self.ionlabel)):
            endconc = np.round(np.mean(self.cc_time[-1][i]),6)
            label = self.ionlabel[i]
            concmess = 'Final average cytoplasmic concentration of'+ ' '+ label + ': '
            loggers.log_info(concmess + str(endconc) + ' mmol/L')

        for i in range(0,len(self.ionlabel)):
            endconc = np.round(np.mean(self.cc_ecm_time[-1][i]),6)
            label = self.ionlabel[i]
            concmess = 'Final average extracellular concentration of'+ ' '+ label + ': '
            loggers.log_info(concmess + str(endconc) + ' mmol/L')

        final_vmean = 1000*np.round(np.mean(self.vm_time[-1]),4)
        vmess = 'Final average cell Vmem of' + ': '
        loggers.log_info(vmess + str(final_vmean) + ' mV')

        if p.ions_dict['H'] == 1:
            final_pH = -np.log10(np.mean((self.cc_time[-1][self.iH])/1000))
            loggers.log_info('Final average cell pH ' + str(np.round(final_pH,2)))

            final_pH_ecm = -np.log10(np.mean((self.cc_ecm_time[-1][self.iH])/1000))
            loggers.log_info('Final environmental pH '+ str(np.round(final_pH_ecm,2)))

        if p.Ca_dyn == 1 and p.ions_dict['Ca']==1:

                endconc_er = np.round(np.mean(self.cc_er[0]),6)
                label = self.ionlabel[self.iCa]
                concmess = 'Final ER concentration of'+ ' '+ label + ': '
                loggers.log_info(concmess + str(endconc_er) + ' mmol/L')

        if p.voltage_dye == 1:

            dye_ecm_final = np.mean(self.cDye_ecm)
            dye_cell_final = np.mean(self.cDye_cell)
            loggers.log_info('Final dye concentration in the environment: '+ str(np.round(dye_ecm_final,6)) + ' mmol/L')
            loggers.log_info('Final average dye concentration in cells: '+ str(np.round(dye_cell_final,6)) + ' mmol/L')

    def runSim(self,cells,p,save=None):
        """
        Drives the time-loop for the main simulation, including gap-junction connections and all dynamics.
        """

        self.tissueInit(cells,p)   # Initialize all structures used for gap junctions, ion channels, and other dynamics

        # Reinitialize all time-data structures
        self.cc_time = []  # data array holding the concentrations at time points
        self.cc_env_time = [] # data array holding environmental concentrations at time points

        self.vm_time = []  # data array holding voltage at time points

        self.dvm_time = []  # data array holding derivative of voltage at time points
        self.time = []     # time values of the simulation

        self.gjopen_time = []   # stores the fractional gap junction open state at each time
        self.fgj_time = []      # stores the gj fluxes for each ion at each time
        self.Igj_time = []      # current for each gj at each time

        self.cc_er_time = []   # retains er concentrations as a function of time
        self.cIP3_time = []    # retains cellular ip3 concentrations as a function of time
        self.cDye_time = []    # retains voltage-sensitive dye concentration as a function of time

        # gap junction specific arrays:
        self.id_gj = np.ones(len(cells.gj_i))  # identity array for gap junction indices...
        self.gjopen = np.ones(len(cells.gj_i))   # holds gap junction open fraction for each gj
        self.gjl = np.zeros(len(cells.gj_i))    # gj length for each gj
        self.gjl[:] = p.gjl
        self.gjsa = np.zeros(len(cells.gj_i))        # gj x-sec surface area for each gj
        self.gjsa[:] = p.gjsa

        vm_to = copy.deepcopy(self.vm)   # create a copy of the original voltage

        # create a time-steps vector:
        tt = np.linspace(0,p.sim_tsteps*p.dt,p.sim_tsteps)

        i = 0 # resample the time vector to save data at specific times:
        tsamples =[]
        resample = p.t_resample
        while i < len(tt)-resample:
            i = i + resample
            tsamples.append(tt[i])
        tsamples = set(tsamples)

        # report
        loggers.log_info('Your simulation is running from '+ str(0) + ' to '+ str(round(p.sim_tsteps*p.dt,3))
                         + ' seconds of in-world time.')

        # get the net, unbalanced charge and corresponding voltage in each cell:
        q_cells = get_charge(self.cc_cells,self.zs,cells.cell_vol,p)
        self.vm = get_volt(q_cells,cells.cell_sa,p)

        if p.plot_while_solving == True:

            checkPlot = viz.PlotWhileSolving(cells,self,p,clrAutoscale = p.autoscale_Vmem, clrMin = p.Vmem_min_clr,
                clrMax = p.Vmem_max_clr)

        do_once = True  # a variable to time the loop only once

        for t in tt:   # run through the loop

            if do_once == True:
                loop_measure = time.time()

            self.dvm = (self.vm - self.vm_to)/p.dt    # calculate the change in the voltage derivative
            self.vm_to = copy.deepcopy(self.vm)       # reassign the history-saving vm

            if p.Ca_dyn ==1 and p.ions_dict['Ca'] == 1:

                self.dcc_ER = (self.cc_er - self.cc_er_to)/p.dt
                self.cc_er_to = copy.deepcopy(self.cc_er)

            # calculate the values of scheduled and dynamic quantities (e.g. ion channel multipliers):
            self.allDynamics(t,p)  # user-scheduled (forced) interventions

            # run the Na-K-ATPase pump:
            self.cc_cells[self.iNa],self.cc_env[self.iNa],self.cc_cells[self.iK],self.cc_env[self.iK], fNa_NaK, fK_NaK =\
                pumpNaKATP(self.cc_cells[self.iNa],self.cc_env[self.iNa],self.cc_cells[self.iK],self.cc_env[self.iK],
                    cells.cell_sa,cells.cell_vol,self.envV,self.vm,self.T,p,self.NaKATP_block)

            # recalculate the net, unbalanced charge and voltage in each cell:
            q_cells = get_charge(self.cc_cells,self.zs,cells.cell_vol,p)
            self.vm = get_volt(q_cells,cells.cell_sa,p)

            if p.ions_dict['Ca'] == 1:

                self.cc_cells[self.iCa],self.cc_env[self.iCa], _ =\
                    pumpCaATP(self.cc_cells[self.iCa],self.cc_env[self.iCa],cells.cell_sa,cells.cell_vol,
                        self.envV,self.vm,self.T,p,self.CaATP_block)

                # recalculate the net, unbalanced charge and voltage in each cell:
                q_cells = get_charge(self.cc_cells,self.zs,cells.cell_vol,p)
                self.vm = get_volt(q_cells,cells.cell_sa,p)

                if p.Ca_dyn ==1:

                    self.cc_er[0],self.cc_cells[self.iCa], _ =\
                        pumpCaER(self.cc_er[0],self.cc_cells[self.iCa],cells.cell_sa,p.ER_vol*cells.cell_vol,
                            cells.cell_vol,self.v_er,self.T,p,self.CaER_block)

                    # recalculate the net, unbalanced charge and voltage in each cell:
                    q_cells = get_charge(self.cc_cells,self.zs,cells.cell_vol,p)
                    self.vm = get_volt(q_cells,cells.cell_sa,p)

                    q_er = get_charge(self.cc_er,self.z_er,p.ER_vol*cells.cell_vol,p)
                    v_er_o = get_volt(q_er,p.ER_sa*cells.cell_sa,p)
                    self.v_er = v_er_o - self.vm

            if p.ions_dict['H'] == 1:

                # electrofuse the H+ ion between the cytoplasm and the environment
                self.cc_env[self.iH],self.cc_cells[self.iH],f_H1 = \
                    electrofuse(self.cc_env[self.iH],self.cc_cells[self.iH],self.Dm_cells[self.iH],self.tm,cells.cell_sa,
                        self.envV,cells.cell_vol,self.zs[self.iH],self.vm,self.T,p)

                # buffer what's happening with H+ flux to or from the cell and environment:
                delH_cell = (f_H1*p.dt/cells.cell_vol)    # relative change in H wrt the cell
                delH_env =  -(f_H1*p.dt/p.vol_env)    # relative change in H wrt to environment

                self.cc_cells[self.iH], self.cc_cells[self.iM], self.cHM_cells = bicarbBuffer(
                    self.cc_cells[self.iH],self.cc_cells[self.iM],self.cHM_cells,delH_cell,p)

                self.cc_env[self.iH], self.cc_env[self.iM], self.cHM_env = bicarbBuffer(
                    self.cc_env[self.iH],self.cc_env[self.iM],self.cHM_env,delH_env,p)

                # recalculate the net, unbalanced charge and voltage in each cell:
                q_cells = get_charge(self.cc_cells,self.zs,cells.cell_vol,p)
                self.vm = get_volt(q_cells,cells.cell_sa,p)


                if p.HKATPase_dyn == 1:

                    # if HKATPase pump is desired, run the H-K-ATPase pump:
                    self.cc_cells[self.iH],self.cc_env[self.iH],self.cc_cells[self.iK],self.cc_env[self.iK], f_H2, _ =\
                    pumpHKATP(self.cc_cells[self.iH],self.cc_env[self.iH],self.cc_cells[self.iK],self.cc_env[self.iK],
                        cells.cell_sa,cells.cell_vol,self.envV,self.vm,self.T,p,self.HKATP_block)

                     # buffer what's happening with H+ flux to or from the cell and environment:
                    delH_cell = (f_H2*p.dt/cells.cell_vol)    # relative change in H wrt the cell
                    delH_env = -(f_H2*p.dt/p.vol_env)    # relative change in H wrt to environment

                    self.cc_cells[self.iH], self.cc_cells[self.iM], self.cHM_cells = bicarbBuffer(
                        self.cc_cells[self.iH],self.cc_cells[self.iM],self.cHM_cells,delH_cell,p)

                    self.cc_env[self.iH], self.cc_env[self.iM], self.cHM_env = bicarbBuffer(
                        self.cc_env[self.iH],self.cc_env[self.iM],self.cHM_env,delH_env,p)

                    # recalculate the net, unbalanced charge and voltage in each cell:
                    q_cells = get_charge(self.cc_cells,self.zs,cells.cell_vol,p)
                    self.vm = get_volt(q_cells,cells.cell_sa,p)

                if p.VATPase_dyn == 1:

                     # if HKATPase pump is desired, run the H-K-ATPase pump:
                    self.cc_cells[self.iH],self.cc_env[self.iH], f_H3 =\
                    pumpVATP(self.cc_cells[self.iH],self.cc_env[self.iH],
                        cells.cell_sa,cells.cell_vol,self.envV,self.vm,self.T,p,self.VATP_block)

                     # buffer what's happening with H+ flux to or from the cell and environment:
                    delH_cell = (f_H3*p.dt/cells.cell_vol)    # relative change in H wrt the cell
                    delH_env = -(f_H3*p.dt/p.vol_env)    # relative change in H wrt to environment

                    self.cc_cells[self.iH], self.cc_cells[self.iM], self.cHM_cells = bicarbBuffer(
                        self.cc_cells[self.iH],self.cc_cells[self.iM],self.cHM_cells,delH_cell,p)

                    self.cc_env[self.iH], self.cc_env[self.iM], self.cHM_env = bicarbBuffer(
                        self.cc_env[self.iH],self.cc_env[self.iM],self.cHM_env,delH_env,p)

                    # recalculate the net, unbalanced charge and voltage in each cell:
                    q_cells = get_charge(self.cc_cells,self.zs,cells.cell_vol,p)
                    self.vm = get_volt(q_cells,cells.cell_sa,p)

            # electro-diffuse all ions (except for proteins, which don't move) across the cell membrane:
            shuffle(cells.gj_i)
            shuffle(self.movingIons)

            for i in self.movingIons:

                # electrodiffusion of ion between cell and extracellular matrix
                self.cc_env[i],self.cc_cells[i],_ = \
                    electrofuse(self.cc_env[i],self.cc_cells[i],self.Dm_cells[i],self.tm,cells.cell_sa,
                        self.envV,cells.cell_vol,self.zs[i],self.vm,self.T,p)

                # recalculate the net, unbalanced charge and voltage in each cell:
                q_cells = get_charge(self.cc_cells,self.zs,cells.cell_vol,p)
                self.vm = get_volt(q_cells,cells.cell_sa,p)

                # calculate volatge difference between cells:
                vmA,vmB = self.vm[cells.gap_jun_i][:,0], self.vm[cells.gap_jun_i][:,1]
                vgj = vmB - vmA

                # determine the open state of gap junctions:
                self.gjopen = self.gj_block*((1.0 - step(abs(vgj),p.gj_vthresh,p.gj_vgrad)) +0.2)

                # determine flux through gap junctions for this ion:
                _,_,fgj = electrofuse(self.cc_cells[i][cells.gap_jun_i][:,0],self.cc_cells[i][cells.gap_jun_i][:,1],
                    self.id_gj*self.D_free[i],self.gjl,self.gjopen*self.gjsa,cells.cell_vol[cells.gap_jun_i][:,0],
                    cells.cell_vol[cells.gap_jun_i][:,1],self.zs[i],vgj,self.T,p)

                # update cell concentration due to gap junction flux:
                self.cc_cells[i] = (self.cc_cells[i]*cells.cell_vol + np.dot((fgj*p.dt), cells.gjMatrix))/cells.cell_vol

                # recalculate the net, unbalanced charge and voltage in each cell:
                q_cells = get_charge(self.cc_cells,self.zs,cells.cell_vol,p)
                self.vm = get_volt(q_cells,cells.cell_sa,p)

                self.fluxes_gj[i] = fgj  # store gap junction flux for this ion

            # determine flux through gap junctions for IP3:
            _,_,fIP3 = electrofuse(self.cIP3[cells.gap_jun_i][:,0],self.cIP3[cells.gap_jun_i][:,1],
                self.id_gj*p.Do_IP3,self.gjl,self.gjopen*self.gjsa,cells.cell_vol[cells.gap_jun_i][:,0],
                cells.cell_vol[cells.gap_jun_i][:,1],p.z_IP3,vgj,self.T,p)

            # update cell IP3 concentration due to gap junction flux:
            self.cIP3 = (self.cIP3*cells.cell_vol + np.dot((fIP3*p.dt), cells.gjMatrix))/cells.cell_vol

            # electrodiffuse IP3 between cell and environment:
            self.cIP3_env,self.cIP3,_ = \
                        electrofuse(self.cIP3_env,self.cIP3,p.Dm_IP3*self.id_cells,self.tm,cells.cell_sa,
                            self.envV,cells.cell_vol,p.z_IP3,self.vm,self.T,p)

            if p.Ca_dyn == 1 and p.ions_dict['Ca'] == 1:
                # electrodiffusion of ions between cell and endoplasmic reticulum
                # Electrodiffusio of calcium
                self.cc_cells[self.iCa],self.cc_er[0],_ = \
                electrofuse(self.cc_cells[self.iCa],self.cc_er[0],self.Dm_er[0],self.tm,p.ER_sa*cells.cell_sa,
                    cells.cell_vol,p.ER_vol*cells.cell_vol,self.z_er[0],self.v_er,self.T,p)

                # Electrodiffusion of charge compensation anion
                self.cc_cells[self.iM],self.cc_er[1],_ = \
                electrofuse(self.cc_cells[self.iM],self.cc_er[1],self.Dm_er[1],self.tm,p.ER_sa*cells.cell_sa,
                    cells.cell_vol,p.ER_vol*cells.cell_vol,self.z_er[1],self.v_er,self.T,p)

                # recalculate the net, unbalanced charge and voltage in each cell:
                q_cells = get_charge(self.cc_cells,self.zs,cells.cell_vol,p)
                self.vm = get_volt(q_cells,cells.cell_sa,p)

                q_er = get_charge(self.cc_er,self.z_er,p.ER_vol*cells.cell_vol,p)
                v_er_o = get_volt(q_er,p.ER_sa*cells.cell_sa,p)
                self.v_er = v_er_o - self.vm

            # if p.voltage_dye=1 electrodiffuse voltage sensitive dye between cell and environment
            if p.voltage_dye ==1:

                self.cDye_env,self.cDye_cell,_ = \
                        electrofuse(self.cDye_env,self.cDye_cell,p.Dm_Dye*self.id_cells,self.tm,cells.cell_sa,
                            self.envV,cells.cell_vol,p.z_Dye,self.vm,self.T,p)

                # determine flux through gap junctions for voltage dye:
                _,_,fDye = electrofuse(self.cDye_cell[cells.gap_jun_i][:,0],self.cDye_cell[cells.gap_jun_i][:,1],
                    self.id_gj*p.Do_Dye,self.gjl,self.gjopen*self.gjsa,cells.cell_vol[cells.gap_jun_i][:,0],
                    cells.cell_vol[cells.gap_jun_i][:,1],p.z_Dye,vgj,self.T,p)

                # update cell voltage-sensitive dye concentration due to gap junction flux:
                self.cDye_cell = (self.cDye_cell*cells.cell_vol + np.dot((fDye*p.dt), cells.gjMatrix))/cells.cell_vol

            if p.dynamic_noise == 1 and p.ions_dict['P']==1:
                # add a random walk on protein concentration to generate dynamic noise:
                self.protein_noise_factor = p.dynamic_noise_level*(np.random.random(len(cells.cell_i)) - 0.5)
                self.cc_cells[self.iP] = self.cc_cells[self.iP]*(1+ self.protein_noise_factor)

                # recalculate the net, unbalanced charge and voltage in each cell:
                q_cells = get_charge(self.cc_cells,self.zs,cells.cell_vol,p)
                self.vm = get_volt(q_cells,cells.cell_sa,p)

            check_v(self.vm)

            if t in tsamples:
                # add the new concentration and voltage data to the time-storage matrices:
                concs = copy.deepcopy(self.cc_cells)
                envsc = copy.deepcopy(self.cc_env)
                flxs = copy.deepcopy(self.fluxes_gj)
                vmm = copy.deepcopy(self.vm)
                dvmm = copy.deepcopy(self.dvm)
                ccIP3 = copy.deepcopy(self.cIP3)
                ggjopen = copy.deepcopy(self.gjopen)

                self.time.append(t)

                self.cc_time.append(concs)
                self.cc_env_time.append(envsc)
                self.vm_time.append(vmm)

                self.fgj_time.append(flxs)
                self.gjopen_time.append(ggjopen)
                self.dvm_time.append(dvmm)

                self.cIP3_time.append(ccIP3)

                if p.voltage_dye ==1:
                    ccDye_cells = copy.deepcopy(self.cDye_cell)
                    self.cDye_time.append(ccDye_cells)

                if p.Ca_dyn == 1 and p.ions_dict['Ca']==1:
                    ccer = copy.deepcopy(self.cc_er)
                    self.cc_er_time.append(ccer)

                if p.plot_while_solving == True:
                    checkPlot.updatePlot(self,p)

                        # get time for loop and estimate total time for simulation
            if do_once == True:
                loop_time = time.time() - loop_measure
                time_estimate = round(loop_time*p.sim_tsteps,2)
                loggers.log_info("This run should take approximately " + str(time_estimate) + ' s to compute...')
                do_once = False

        # End off by calculating the current through the gap junction network:
        self.Igj_time = []
        for tflux in self.fgj_time:
            igj=0
            for zi, flx in zip(self.zs,tflux):
                igj = igj+ zi*flx

            igj = p.F*igj
            self.Igj_time.append(igj)

        if save==True:
            celf = copy.deepcopy(self)
            datadump = [celf,cells,p]
            fh.saveSim(self.savedSim,datadump)
            message_2 = 'Simulation run saved to' + ' ' + p.sim_path
            loggers.log_info(message_2)

        for i in range(0,len(self.ionlabel)):
            endconc = np.round(np.mean(self.cc_time[-1][i]),6)
            label = self.ionlabel[i]
            concmess = 'Final average cytoplasmic concentration of'+ ' '+ label + ': '
            loggers.log_info(concmess + str(endconc) + ' mmol/L')

        for i in range(0,len(self.ionlabel)):
            endconc = np.round(np.mean(self.cc_env_time[-1][i]),6)
            label = self.ionlabel[i]
            concmess = 'Final environmental concentration of'+ ' '+ label + ': '
            loggers.log_info(concmess + str(endconc) + ' mmol/L')


        final_vmean = 1000*np.round(np.mean(self.vm_time[-1]),4)
        vmess = 'Final average cell Vmem of ' + ': '
        loggers.log_info(vmess + str(final_vmean) + ' mV')

        if p.ions_dict['H'] == 1:
            final_pH = -np.log10(np.mean((self.cc_time[-1][self.iH])/1000))
            loggers.log_info('Final average cell pH '+ str(np.round(final_pH,2)))

            final_pH_env = -np.log10(np.mean((self.cc_env_time[-1][self.iH])/1000))
            loggers.log_info('Final environmental pH '+ str(np.round(final_pH_env,2)))

        IP3_env_final = np.mean(self.cIP3_env)
        IP3_cell_final = np.mean(self.cIP3)
        loggers.log_info('Final IP3 concentration in the environment: ' + str(np.round(IP3_env_final,6)) + ' mmol/L')
        loggers.log_info('Final average IP3 concentration in cells: ' + str(np.round(IP3_cell_final,6)) + ' mmol/L')

        if p.Ca_dyn == 1 and p.ions_dict['Ca'] == 1:

            endconc_er = np.round(np.mean(self.cc_er[0]),6)
            label = self.ionlabel[self.iCa]
            concmess = 'Final average ER concentration of'+ ' '+ label + ': '
            loggers.log_info(concmess + str(endconc_er) + ' mmol/L')

        if p.voltage_dye ==1:
            dye_env_final = np.mean(self.cDye_env)
            dye_cell_final = np.mean(self.cDye_cell)
            loggers.log_info('Final dye concentration in the environment: '+ str(np.round(dye_env_final,6))
                             + ' mmol/L')
            loggers.log_info('Final average dye concentration in cells: ' +  str(np.round(dye_cell_final,6)) +
                             ' mmol/L')

        plt.close()
        loggers.log_info('Simulation completed successfully.')

    def runSim_ECM(self,cells,p,save=None):
        """
        Drives the time-loop for the main simulation, including gap-junction connections and all dynamics.
        """

        self.tissueInit_ECM(cells,p)   # Initialize all structures used for gap junctions, ion channels, and other dynamics

        # Reinitialize all time-data structures
        self.cc_time = []  # data array holding the concentrations at time points
        self.cc_ecm_time = [] # data array holding extracellular concentrations at time points

        self.vm_time = []  # data array holding voltage at time points

        self.dvm_time = []  # data array holding derivative of voltage at time points
        self.time = []     # time values of the simulation

        self.gjopen_time = []   # stores the fractional gap junction open state at each time
        self.fgj_time = []      # stores the gj fluxes for each ion at each time
        self.Igj_time = []      # current for each gj at each time

        self.cc_er_time = []   # retains er concentrations as a function of time
        self.cIP3_time = []    # retains cellular ip3 concentrations as a function of time
        self.cDye_time = []    # retains voltage-sensitive dye concentration as a function of time

        # gap junction specific arrays:
        self.id_gj = np.ones(len(cells.gj_i))  # identity array for gap junction indices...
        self.gjopen = np.ones(len(cells.gj_i))   # holds gap junction open fraction for each gj
        self.gjl = np.zeros(len(cells.gj_i))    # gj length for each gj
        self.gjl[:] = p.gjl
        self.gjsa = np.zeros(len(cells.gj_i))        # gj x-sec surface area for each gj
        self.gjsa[:] = p.gjsa

        vm_to = copy.deepcopy(self.vm)   # create a copy of the original voltage

        # create a time-steps vector:
        tt = np.linspace(0,p.sim_tsteps*p.dt,p.sim_tsteps)

        i = 0 # resample the time vector to save data at specific times:
        tsamples =[]
        resample = p.t_resample
        while i < len(tt)-resample:
            i = i + resample
            tsamples.append(tt[i])
        tsamples = set(tsamples)

        # report
        loggers.log_info('Your simulation is running from '+ str(0) + ' to '+ str(round(p.sim_tsteps*p.dt,3))
                         + ' seconds of in-world time.')

        # get the net, unbalanced charge and corresponding voltage in each cell:
        self.update_V_ecm(cells,p)

        if p.plot_while_solving == True:

            checkPlot = viz.PlotWhileSolving(cells,self,p,clrAutoscale = p.autoscale_Vmem, clrMin = p.Vmem_min_clr,
                clrMax = p.Vmem_max_clr)

        do_once = True  # a variable to time the loop only once

        for t in tt:   # run through the loop

            if do_once == True:
                loop_measure = time.time()

            self.dvm = (self.vm - self.vm_to)/p.dt    # calculate the change in the voltage derivative
            self.vm_to = copy.deepcopy(self.vm)       # reassign the history-saving vm

            if p.Ca_dyn ==1 and p.ions_dict['Ca'] == 1:

                self.dcc_ER = (self.cc_er - self.cc_er_to)/p.dt
                self.cc_er_to = copy.deepcopy(self.cc_er)

            # calculate the values of scheduled and dynamic quantities (e.g. ion channel multipliers):
            self.allDynamics(t,p)  # user-scheduled (forced) interventions

            # run the Na-K-ATPase pump:
            _,_,_,_,fNa_NaK, fK_NaK =\
                pumpNaKATP(self.cc_cells[self.iNa][cells.mem_to_cells],self.cc_ecm[self.iNa][cells.mem_to_ecm],
                    self.cc_cells[self.iK][cells.mem_to_cells],self.cc_ecm[self.iK][cells.mem_to_ecm],
                    cells.mem_sa,cells.cell_vol[cells.mem_to_cells],cells.ecm_vol[cells.mem_to_ecm],
                    self.vm,self.T,p,self.NaKATP_block)

            # update the concentrations
            self.update_C_ecm(self.iNa,fNa_NaK,cells,p)
            self.update_C_ecm(self.iK,fK_NaK,cells,p)

            # recalculate the net, unbalanced charge and voltage in each cell:
            self.update_V_ecm(cells,p)

            if p.ions_dict['Ca'] == 1:

                _,_, f_CaATP =\
                    pumpCaATP(self.cc_cells[self.iCa][cells.mem_to_cells],self.cc_ecm[self.iCa][cells.mem_to_ecm],
                        cells.mem_sa, cells.cell_vol[cells.mem_to_cells], cells.ecm_vol[cells.mem_to_ecm],
                        self.vm,self.T,p,self.CaATP_block)

                # update calcium concentrations in cell and ecm:
                self.update_C_ecm(self.iCa,f_CaATP,cells,p)

                # recalculate the net, unbalanced charge and voltage in each cell:
                self.update_V_ecm(cells,p)

                if p.Ca_dyn ==1:

                    self.cc_er[0],self.cc_cells[self.iCa], _ =\
                        pumpCaER(self.cc_er[0],self.cc_cells[self.iCa],cells.cell_sa,p.ER_vol*cells.cell_vol,
                            cells.cell_vol,self.v_er,self.T,p,self.CaER_block)

                    # recalculate the net, unbalanced charge and voltage in each cell:
                    self.update_V_ecm(cells,p)

                    q_er = get_charge(self.cc_er,self.z_er,p.ER_vol*cells.cell_vol,p)
                    self.v_er = get_volt(q_er,p.ER_sa*cells.cell_sa,p)

            if p.ions_dict['H'] == 1:

                self.Hplus_electrofuse_ecm(cells,p)

                if p.HKATPase_dyn == 1:

                    self.Hplus_HKATP_ecm(cells,p)

                if p.VATPase_dyn == 1:

                    self.Hplus_VATP_ecm(cells,p)

            # electro-diffuse all ions (except for proteins, which don't move) across the cell membrane:
            shuffle(cells.gj_i)
            shuffle(self.movingIons)

            for i in self.movingIons:

                # electrodiffusion of ion between cell and extracellular matrix
                _,_,f_ED = \
                    electrofuse(self.cc_ecm[i][cells.mem_to_ecm],self.cc_cells[i][cells.mem_to_cells],
                        self.Dm_mems[i],self.tm[cells.mem_to_cells],cells.mem_sa,
                        cells.ecm_vol[cells.mem_to_ecm],cells.cell_vol,self.zs[i],self.vm,self.T,p)

                # update ion concentrations in cell and ecm:
                self.update_C_ecm(i,f_ED,cells,p)

                # recalculate the net, unbalanced charge and voltage in each cell:
                self.update_V_ecm(cells,p)

                # calculate voltage difference between cells:
                vmA,vmB = self.vm_cell_ave[cells.gap_jun_i][:,0], self.vm_cell_ave[cells.gap_jun_i][:,1]
                vgj = vmB - vmA

                # determine the open state of gap junctions:
                self.gjopen = self.gj_block*((1.0 - step(abs(vgj),p.gj_vthresh,p.gj_vgrad)) +0.2)

                # determine flux through gap junctions for this ion:
                _,_,fgj = electrofuse(self.cc_cells[i][cells.gap_jun_i][:,0],self.cc_cells[i][cells.gap_jun_i][:,1],
                    self.id_gj*self.D_free[i],self.gjl,self.gjopen*self.gjsa,cells.cell_vol[cells.gap_jun_i][:,0],
                    cells.cell_vol[cells.gap_jun_i][:,1],self.zs[i],vgj,self.T,p)

                # update cell concentration due to gap junction flux:
                self.cc_cells[i] = (self.cc_cells[i]*cells.cell_vol + np.dot((fgj*p.dt), cells.gjMatrix))/cells.cell_vol

                # recalculate the net, unbalanced charge and voltage in each cell:
                self.update_V_ecm(cells,p)

                self.fluxes_gj[i] = fgj  # store gap junction flux for this ion

            # determine flux through gap junctions for IP3:
            _,_,fIP3 = electrofuse(self.cIP3[cells.gap_jun_i][:,0],self.cIP3[cells.gap_jun_i][:,1],
                self.id_gj*p.Do_IP3,self.gjl,self.gjopen*self.gjsa,cells.cell_vol[cells.gap_jun_i][:,0],
                cells.cell_vol[cells.gap_jun_i][:,1],p.z_IP3,vgj,self.T,p)

            # update cell IP3 concentration due to gap junction flux:
            self.cIP3 = (self.cIP3*cells.cell_vol + np.dot((fIP3*p.dt), cells.gjMatrix))/cells.cell_vol

            # electrodiffuse IP3 between cell and environment:
            _,_,flux_IP3 = \
                        electrofuse(self.cIP3_ecm[cells.mem_to_ecm],self.cIP3[cells.mem_to_cells],
                            p.Dm_IP3*self.id_cells[cells.mem_to_cells],self.tm[cells.mem_to_cells],cells.mem_sa,
                            cells.ecm_vol[cells.mem_to_ecm],cells.cell_vol[cells.mem_to_cells],p.z_IP3,self.vm,self.T,p)

            # update the IP3 concentrations in the cell and ecm due to ED fluxes at membrane
            self.cIP3 = self.cIP3 + \
                                np.dot((flux_IP3/cells.cell_vol[cells.mem_to_cells])*p.dt,cells.cell_UpdateMatrix)

            self.cIP3_ecm = self.cIP3_ecm - \
                                np.dot((flux_IP3/cells.ecm_vol[cells.mem_to_ecm])*p.dt,cells.ecm_UpdateMatrix)

            if p.Ca_dyn == 1 and p.ions_dict['Ca'] == 1:

                # electrodiffusion of ions between cell and endoplasmic reticulum
                # Electrodiffusio of calcium
                self.cc_cells[self.iCa],self.cc_er[0],_ = \
                electrofuse(self.cc_cells[self.iCa],self.cc_er[0],self.Dm_er[0],self.tm,p.ER_sa*cells.cell_sa,
                    cells.cell_vol,p.ER_vol*cells.cell_vol,self.z_er[0],self.v_er,self.T,p)

                # Electrodiffusion of charge compensation anion
                self.cc_cells[self.iM],self.cc_er[1],_ = \
                electrofuse(self.cc_cells[self.iM],self.cc_er[1],self.Dm_er[1],self.tm,p.ER_sa*cells.cell_sa,
                    cells.cell_vol,p.ER_vol*cells.cell_vol,self.z_er[1],self.v_er,self.T,p)

                # recalculate the net, unbalanced charge and voltage in each cell:
                self.update_V_ecm(cells,p)

                q_er = get_charge(self.cc_er,self.z_er,p.ER_vol*cells.cell_vol,p)
                self.v_er = get_volt(q_er,p.ER_sa*cells.cell_sa,p)

            # if p.voltage_dye=1 electrodiffuse voltage sensitive dye between cell and environment
            if p.voltage_dye ==1:

                _,_,flux_dye = \
                        electrofuse(self.cDye_ecm[cells.mem_to_ecm],self.cDye_cell[cells.mem_to_cells],
                            p.Dm_Dye*self.id_cells[cells.mem_to_cells],self.tm[cells.mem_to_cells],cells.mem_sa,
                            cells.ecm_vol[cells.mem_to_ecm],cells.cell_vol[cells.mem_to_cells],p.z_Dye,self.vm,self.T,p)

                 # update the dye concentrations in the cell and ecm due to ED fluxes at membrane
                self.cDye_cell = self.cDye_cell + \
                                    np.dot((flux_dye/cells.cell_vol[cells.mem_to_cells])*p.dt,cells.cell_UpdateMatrix)

                self.cDye_ecm = self.cDye_ecm - \
                                    np.dot((flux_dye/cells.ecm_vol[cells.mem_to_ecm])*p.dt,cells.ecm_UpdateMatrix)

                # determine flux through gap junctions for voltage dye:
                _,_,fDye = electrofuse(self.cDye_cell[cells.gap_jun_i][:,0],self.cDye_cell[cells.gap_jun_i][:,1],
                    self.id_gj*p.Do_Dye,self.gjl,self.gjopen*self.gjsa,cells.cell_vol[cells.gap_jun_i][:,0],
                    cells.cell_vol[cells.gap_jun_i][:,1],p.z_Dye,vgj,self.T,p)

                # update cell voltage-sensitive dye concentration due to gap junction flux:
                self.cDye_cell = (self.cDye_cell*cells.cell_vol + np.dot((fDye*p.dt), cells.gjMatrix))/cells.cell_vol

            if p.dynamic_noise == 1 and p.ions_dict['P']==1:
                # add a random walk on protein concentration to generate dynamic noise:
                self.protein_noise_factor = p.dynamic_noise_level*(np.random.random(len(cells.cell_i)) - 0.5)
                self.cc_cells[self.iP] = self.cc_cells[self.iP]*(1+ self.protein_noise_factor)

                # recalculate the net, unbalanced charge and voltage in each cell:
                self.update_V_ecm(cells,p)

            check_v(self.vm_cell_ave)

            if t in tsamples:
                # add the new concentration and voltage data to the time-storage matrices:
                concs = copy.deepcopy(self.cc_cells)
                ecmsc = copy.deepcopy(self.cc_ecm)
                flxs = copy.deepcopy(self.fluxes_gj)
                vmm = copy.deepcopy(self.vm)
                dvmm = copy.deepcopy(self.dvm)
                ccIP3 = copy.deepcopy(self.cIP3)
                ggjopen = copy.deepcopy(self.gjopen)

                self.time.append(t)

                self.cc_time.append(concs)
                self.cc_ecm_time.append(ecmsc)
                self.vm_time.append(vmm)

                self.fgj_time.append(flxs)
                self.gjopen_time.append(ggjopen)
                self.dvm_time.append(dvmm)

                self.cIP3_time.append(ccIP3)

                if p.voltage_dye ==1:
                    ccDye_cells = copy.deepcopy(self.cDye_cell)
                    self.cDye_time.append(ccDye_cells)

                if p.Ca_dyn == 1 and p.ions_dict['Ca']==1:
                    ccer = copy.deepcopy(self.cc_er)
                    self.cc_er_time.append(ccer)

                if p.plot_while_solving == True:
                    checkPlot.updatePlot(self,p)

                        # get time for loop and estimate total time for simulation
            if do_once == True:
                loop_time = time.time() - loop_measure
                time_estimate = round(loop_time*p.sim_tsteps,2)
                loggers.log_info("This run should take approximately " + str(time_estimate) + ' s to compute...')
                do_once = False

        # End off by calculating the current through the gap junction network:
        self.Igj_time = []
        for tflux in self.fgj_time:
            igj=0
            for zi, flx in zip(self.zs,tflux):
                igj = igj+ zi*flx

            igj = p.F*igj
            self.Igj_time.append(igj)

        if save==True:
            celf = copy.deepcopy(self)
            datadump = [celf,cells,p]
            fh.saveSim(self.savedSim,datadump)
            message_2 = 'Simulation run saved to' + ' ' + p.sim_path
            loggers.log_info(message_2)

        for i in range(0,len(self.ionlabel)):
            endconc = np.round(np.mean(self.cc_time[-1][i]),6)
            label = self.ionlabel[i]
            concmess = 'Final average cytoplasmic concentration of'+ ' '+ label + ': '
            loggers.log_info(concmess + str(endconc) + ' mmol/L')

        for i in range(0,len(self.ionlabel)):
            endconc = np.round(np.mean(self.cc_env_time[-1][i]),6)
            label = self.ionlabel[i]
            concmess = 'Final extracellular concentration of'+ ' '+ label + ': '
            loggers.log_info(concmess + str(endconc) + ' mmol/L')

        final_vmean = 1000*np.round(np.mean(self.vm_time[-1]),4)
        vmess = 'Final average cell Vmem of ' + ': '
        loggers.log_info(vmess + str(final_vmean) + ' mV')

        if p.ions_dict['H'] == 1:
            final_pH = -np.log10(np.mean((self.cc_time[-1][self.iH])/1000))
            loggers.log_info('Final average cell pH '+ str(np.round(final_pH,2)))

            final_pH_ecm = -np.log10(np.mean((self.cc_ecm_time[-1][self.iH])/1000))
            loggers.log_info('Final extracellular pH '+ str(np.round(final_pH_ecm,2)))

        IP3_ecm_final = np.mean(self.cIP3_ecm)
        IP3_cell_final = np.mean(self.cIP3)
        loggers.log_info('Final extracellular IP3 concentration: ' + str(np.round(IP3_ecm_final,6)) + ' mmol/L')
        loggers.log_info('Final average IP3 concentration in cells: ' + str(np.round(IP3_cell_final,6)) + ' mmol/L')

        if p.Ca_dyn == 1 and p.ions_dict['Ca'] == 1:

            endconc_er = np.round(np.mean(self.cc_er[0]),6)
            label = self.ionlabel[self.iCa]
            concmess = 'Final average ER concentration of'+ ' '+ label + ': '
            loggers.log_info(concmess + str(endconc_er) + ' mmol/L')

        if p.voltage_dye ==1:
            dye_ecm_final = np.mean(self.cDye_ecm)
            dye_cell_final = np.mean(self.cDye_cell)
            loggers.log_info('Final extracellular dye concentration: '+ str(np.round(dye_ecm_final,6))
                             + ' mmol/L')
            loggers.log_info('Final average dye concentration in cells: ' +  str(np.round(dye_cell_final,6)) +
                             ' mmol/L')

        plt.close()
        loggers.log_info('Simulation completed successfully.')

    def allDynamics(self,t,p):

        target_length = len(self.scheduled_target_inds)

        if p.scheduled_options['Na_mem'] != 0:

            if p.ions_dict['Na'] == 0 or target_length == 0:
                pass

            else:

                effector_Na = pulse(t,self.t_on_Namem,self.t_off_Namem,self.t_change_Namem)

                self.Dm_scheduled[self.iNa][self.scheduled_target_inds] = self.mem_mult_Namem*effector_Na*p.Dm_Na

        if p.scheduled_options['K_mem'] != 0:

            if p.ions_dict['K'] == 0 or target_length == 0:
                pass

            else:

                effector_K = pulse(t,self.t_on_Kmem,self.t_off_Kmem,self.t_change_Kmem)

                self.Dm_scheduled[self.iK][self.scheduled_target_inds] = self.mem_mult_Kmem*effector_K*p.Dm_K

        if p.scheduled_options['Cl_mem'] != 0:

            if p.ions_dict['Cl'] == 0 or target_length == 0:
                pass

            else:

                effector_Cl = pulse(t,self.t_on_Clmem,self.t_off_Clmem,self.t_change_Clmem)

                self.Dm_scheduled[self.iCl][self.scheduled_target_inds] = self.mem_mult_Clmem*effector_Cl*p.Dm_Cl

        if p.scheduled_options['Ca_mem'] != 0:

            if p.ions_dict['Ca'] == 0 or target_length == 0:
                pass

            else:

                effector_Ca = pulse(t,self.t_on_Camem,self.t_off_Camem,self.t_change_Camem)

                self.Dm_scheduled[self.iCa][self.scheduled_target_inds] = self.mem_mult_Camem*effector_Ca*p.Dm_Ca

        if p.scheduled_options['IP3'] != 0:

            self.cIP3[self.scheduled_target_inds] = self.cIP3[self.scheduled_target_inds] + self.rate_IP3*pulse(t,self.t_onIP3,
                self.t_offIP3,self.t_changeIP3)

        if p.global_options['K_env'] != 0:

            effector_Kenv = pulse(t,self.t_on_Kenv,self.t_off_Kenv,self.t_change_Kenv)

            self.cc_env[self.iK][:] = self.mem_mult_Kenv*effector_Kenv*p.cK_env + p.cK_env

        if p.global_options['Cl_env'] != 0 and p.ions_dict['Cl'] == 1:

            effector_Clenv = pulse(t,self.t_on_Clenv,self.t_off_Clenv,self.t_change_Clenv)

            self.cc_env[self.iCl][:] = self.mem_mult_Clenv*effector_Clenv*p.cCl_env + p.cCl_env

        if p.global_options['Na_env'] != 0:

            effector_Naenv = pulse(t,self.t_on_Naenv,self.t_off_Naenv,self.t_change_Naenv)

            self.cc_env[self.iNa][:] = self.mem_mult_Naenv*effector_Naenv*p.cNa_env + p.cNa_env

        if p.global_options['T_change'] != 0:

            self.T = self.multT*pulse(t,self.tonT,self.toffT,self.trampT)*p.T + p.T

        if p.global_options['gj_block'] != 0:

            self.gj_block = (1.0 - pulse(t,self.tonGJ,self.toffGJ,self.trampGJ))

        if p.global_options['NaKATP_block'] != 0:

            self.NaKATP_block = (1.0 - pulse(t,self.tonNK,self.toffNK,self.trampNK))

        if p.global_options['HKATP_block'] != 0:

            self.HKATP_block = (1.0 - pulse(t,self.tonHK,self.toffHK,self.trampHK))


        # Voltage gated channel effects ................................................................................

        dvsign = np.sign(self.dvm)

        if p.vg_options['Na_vg'] != 0:

            if p.ions_dict['Na'] == 0:
                pass

            else:

                # Logic phase 1: find out which cells have activated their vgNa channels
                truth_vmGTvon_Na = self.vm > self.v_activate_Na  # returns bools of vm that are bigger than threshhold
                #truth_depol_Na = dvsign==1  # returns bools of vm that are bigger than threshhold
                truth_not_inactivated_Na = self.inactivated_Na == 0  # return bools of vm that can activate
                truth_vgNa_Off = self.vgNa_state == 0 # hasn't been turned on yet

                # find the cell indicies that correspond to all statements of logic phase 1:
                inds_activate_Na = (truth_vmGTvon_Na*truth_not_inactivated_Na*truth_vgNa_Off*
                                    self.target_cells).nonzero()

                self.vgNa_state[inds_activate_Na] = 1 # open the channel
                self.vgNa_aliveTimer[inds_activate_Na] = t + self.t_alive_Na # set the timers for the total active state
                self.vgNa_deadTimer[inds_activate_Na] = 0  # reset any timers for an inactivated state to zero

                # Logic phase 2: find out which cells have closed their gates due to crossing inactivating voltage:
                truth_vgNa_On = self.vgNa_state == 1  # channel must be on already
                truth_vmGTvoff_Na = self.vm > self.v_inactivate_Na  # bools of cells that have vm greater than shut-off volts

                inds_inactivate_Na = (truth_vgNa_On*truth_vmGTvoff_Na*self.target_cells).nonzero()

                self.vgNa_state[inds_inactivate_Na] = 0    # close the vg sodium channels
                self.inactivated_Na[inds_inactivate_Na] = 1   # switch these so cells do not re-activate
                self.vgNa_aliveTimer[inds_inactivate_Na] = 0            # reset any alive timers to zero
                self.vgNa_deadTimer[inds_inactivate_Na] = t + self.t_dead_Na # set the timer of the inactivated state

                 # Logic phase 3: find out if cell activation state has timed out, also rendering inactivated state:

                truth_vgNa_act_timeout = self.vgNa_aliveTimer < t   # find cells that have timed out their vgNa open state
                truth_vgNa_On = self.vgNa_state == 1 # ensure the vgNa is indeed open
                inds_timeout_Na_act = (truth_vgNa_act_timeout*truth_vgNa_On*self.target_cells).nonzero()

                self.vgNa_state[inds_timeout_Na_act] = 0             # set the state to closed
                self.vgNa_aliveTimer[inds_timeout_Na_act] = 0            # reset the timers to zero
                self.inactivated_Na[inds_timeout_Na_act] = 1    # inactivate the channel so it can't reactivate
                self.vgNa_deadTimer[inds_timeout_Na_act] = t + self.t_dead_Na # set the timer of the inactivated state

                # Logic phase 4: find out if inactivation timers have timed out:
                truth_vgNa_inact_timeout = self.vgNa_deadTimer <t  # find cells that have timed out their vgNa inact state
                truth_vgNa_Off = self.vgNa_state == 0 # check to make sure these channels are indeed closed
                inds_timeout_Na_inact = (truth_vgNa_inact_timeout*truth_vgNa_Off*self.target_cells).nonzero()

                self.vgNa_deadTimer[inds_timeout_Na_inact] = 0    # reset the inactivation timer
                self.inactivated_Na[inds_timeout_Na_inact] = 0    # remove inhibition to activation

                # Logic phase 5: find out if cells have passed below threshhold to become deactivated:
                truth_vmLTvreact_Na = self.vm < self.v_deactivate_Na # voltage is lower than the deactivate voltage

                inds_deactivate_Na = (truth_vmLTvreact_Na*self.target_cells).nonzero()

                self.inactivated_Na[inds_deactivate_Na] = 0  # turn any inhibition to activation off
                self.vgNa_state[inds_deactivate_Na] = 0   # shut the Na channel off if it's on
                self.vgNa_aliveTimer[inds_deactivate_Na] = 0       # reset any alive-timers to zero
                self.vgNa_deadTimer[inds_deactivate_Na] = 0   # reset any dead-timers to zero


                # Define ultimate activity of the vgNa channel:

                self.Dm_vg[self.iNa] = self.maxDmNa*self.vgNa_state

        if p.vg_options['K_vg'] != 0:

            if p.ions_dict['K'] == 0:
                pass

            else:

                # detecting channels to turn on:

                truth_vmGTvon_K = self.vm > self.v_on_K  # bools for cells with vm greater than the on threshold for vgK
                truth_depol_K = dvsign == 1  # bools matrix for cells that are depolarizing
                truth_vgK_OFF = self.vgK_state == 0   # bools matrix for cells that are in the off state

                # cells at these indices will become activated in this time step:
                inds_activate_K = (truth_vmGTvon_K*truth_depol_K*truth_vgK_OFF*self.target_cells).nonzero()
                self.vgK_state[inds_activate_K] = 1  # set the state of these channels to "open"
                self.vgK_OFFtime[inds_activate_K] = self.t_alive_K + t  # set the time at which these channels will close

                #  detecting channels to turn off:
                truth_vgK_ON = self.vgK_state == 1  # detect cells that are in their on state
                truth_vgK_timeout = self.vgK_OFFtime < t     # detect the cells that have expired off timers
                inds_deactivate_K = (truth_vgK_ON*truth_vgK_timeout*self.target_cells).nonzero()
                self.vgK_state[inds_deactivate_K] = 0 # turn off the channels to closed
                self.vgK_OFFtime[inds_deactivate_K] = 0

                inds_open_K = (self.vgK_state == 1).nonzero()
                self.active_K[inds_open_K] = 1

                inds_closed_K =(self.vgK_state == 0).nonzero()
                self.active_K[inds_closed_K] = 0

                self.Dm_vg[self.iK] = self.maxDmK*self.active_K


        if p.vg_options['Ca_vg'] != 0:

            if p.ions_dict['Ca'] == 0:
                pass

            else:
                # detect condition to turn vg_Ca channel on:
                truth_vmGTvon_Ca = self.vm > self.v_on_Ca  # bools for cells with vm greater than the on threshold for vgK
                truth_caLTcaOff = self.cc_cells[self.iCa] < self.ca_lower_ca # check that cellular calcium is below inactivating Ca
                truth_depol_Ca = dvsign == 1  # bools matrix for cells that are depolarizing
                truth_vgCa_OFF = self.vgCa_state == 0   # bools matrix for cells that are in the off state

                # cells at these indices will become activated in this time step:
                inds_activate_Ca = (truth_vmGTvon_Ca*truth_depol_Ca*truth_caLTcaOff*truth_vgCa_OFF*self.target_cells).nonzero()
                self.vgCa_state[inds_activate_Ca] = 1  # set the state of these channels to "open"

                # detect condition to turn off vg_Ca channel:
                truth_caGTcaOff = self.cc_cells[self.iCa] > self.ca_upper_ca   # check that calcium exceeds maximum
                truth_vgCa_ON = self.vgCa_state == 1 # check that the channel is on
                inds_inactivate_Ca = (truth_caGTcaOff*truth_vgCa_ON*self.target_cells).nonzero()
                self.vgCa_state[inds_inactivate_Ca] = 0

                # additional condition to turn off vg_Ca via depolarizing voltage:
                truth_vmGTvcaOff = self.vm > self.v_off_Ca
                inds_inactivate_Ca_2 = (truth_vmGTvcaOff*self.target_cells*truth_vgCa_ON).nonzero()
                self.vgCa_state[inds_inactivate_Ca_2] = 0


                inds_open_Ca = (self.vgCa_state == 1).nonzero()
                self.active_Ca[inds_open_Ca] = 1

                inds_closed_Ca =(self.vgCa_state == 0).nonzero()
                self.active_Ca[inds_closed_Ca] = 0

                self.Dm_vg[self.iCa] = self.maxDmCa*self.active_Ca

        if p.vg_options['K_cag'] != 0:

            if p.ions_dict['Ca'] == 0:
                pass

            else:

                inds_cagK_targets = (self.target_cells).nonzero()

                self.active_Kcag[inds_cagK_targets] = hill(self.cc_cells[self.iCa][inds_cagK_targets],
                    self.Kcag_halfmax,self.Kcag_n)

                self.Dm_cag[self.iK] = self.maxDmKcag*self.active_Kcag

        # finally, add together all effects to make change on the cell membrane permeabilities:
        self.Dm_cells = self.Dm_scheduled + self.Dm_vg + self.Dm_cag + self.Dm_base    # FIXME this is the one spot where original and ECM values contradict each other

        # Calcium Dynamics options including Calicum-Induced-Calcium-Release (CICR) and IP3 mediated calcium release....

        if p.ions_dict['Ca'] ==1 and p.Ca_dyn == 1:

            if p.Ca_dyn_options['CICR'] != 0:

                dcc_CaER_sign = np.sign(self.dcc_ER[0])

                if len(p.Ca_dyn_options['CICR'][1])==0:
                    term_Ca_reg = 1.0

                else:
                    term_Ca_reg = (np.exp(-((self.cc_cells[self.iCa]-self.midCaR)**2)/((2*self.widthCaR)**2)))

                if len(p.Ca_dyn_options['CICR'][2]) == 0:
                    term_IP3_reg = 1.0

                else:
                    term_IP3_reg = hill(self.cIP3,self.KhmIP3,self.n_IP3)

                if p.FMmod == 1:
                    span = self.topCa - self.bottomCa
                    FMmod = p.ip3FM*span
                    topCa = self.topCa - FMmod*term_IP3_reg
                else:
                    topCa = self.topCa

                truth_overHighCa = self.cc_er[0] >=  topCa
                truth_increasingCa = dcc_CaER_sign == 1
                truth_alreadyClosed = self.stateER == 0.0
                inds_open_ER = (truth_overHighCa*truth_increasingCa*truth_alreadyClosed).nonzero()

                truth_underBottomCa = self.cc_er[0]< self.bottomCa
                truth_decreasingCa = dcc_CaER_sign == -1
                truth_alreadyOpen = self.stateER == 1.0
                inds_close_ER = (truth_underBottomCa*truth_alreadyOpen).nonzero()

                self.stateER[inds_open_ER] = 1.0
                self.stateER[inds_close_ER] = 0.0

                self.Dm_er_CICR[0] = self.maxDmCaER*self.stateER*term_IP3_reg*term_Ca_reg

                self.Dm_er = self.Dm_er_CICR + self.Dm_er_base

    def update_V_ecm(self,cells,p):

        q_cells = get_charge(self.cc_cells, self.zs, cells.cell_vol, p)
        q_ecm = get_charge(self.cc_ecm, self.zs, cells.ecm_vol, p)
        self.vm, v_cell_at_mem, v_ecm_at_mem = get_volt_ECM(q_cells,q_ecm,cells)   # calculate voltages
        self.vm_cell_ave = cell_ave(cells,self.vm)  # calculate average vm for each cell

    def update_C_ecm(self,ion_i,flux,cells,p):

        self.cc_cells[ion_i] = self.cc_cells[ion_i] + \
                          np.dot((flux/cells.cell_vol[cells.mem_to_cells])*p.dt,cells.cell_UpdateMatrix)

        self.cc_ecm[ion_i] = self.cc_ecm[ion_i] - \
                                np.dot((flux/cells.ecm_vol[cells.mem_to_ecm])*p.dt,cells.ecm_UpdateMatrix)

    def Hplus_electrofuse_ecm(self,cells,p):

        # electrofuse the H+ ion between the cytoplasm and the ecms
        _,_,f_H1 = \
            electrofuse(self.cc_ecm[self.iH][cells.mem_to_ecm],self.cc_cells[self.iH][cells.mem_to_cells],
                self.Dm_mems[self.iH],self.tm[cells.mem_to_cells],cells.mem_sa,
                cells.ecm_vol[cells.mem_to_ecm],cells.cell_vol[cells.mem_to_cells],self.zs[self.iH],
                self.vm,self.T,p)

        H_cell_to = copy.deepcopy(self.cc_cells[self.iH])     # keep track of original H+ in cell and ecm
        H_ecm_to = copy.deepcopy(self.cc_ecm[self.iH])

        # update the concentration in the cytoplasm and ecms
        self.update_C_ecm(self.iH,f_H1,cells,p)

        # buffer what's happening with H+ flux to or from the cell and environment:
        delH_cell = self.cc_cells[self.iH] - H_cell_to    # relative change in H wrt the cell
        delH_ecm = self.cc_ecm[self.iH] - H_ecm_to    # relative change in H wrt to environment

        self.cc_cells[self.iH], self.cc_cells[self.iM], self.cHM_cells = bicarbBuffer(
            self.cc_cells[self.iH],self.cc_cells[self.iM],self.cHM_cells,delH_cell,p)

        self.cc_ecm[self.iH], self.cc_ecm[self.iM], self.cHM_ecm = bicarbBuffer(
            self.cc_ecm[self.iH],self.cc_ecm[self.iM],self.cHM_ecm,delH_ecm,p)

        # recalculate the net, unbalanced charge and voltage in each cell:
        self.update_V_ecm(cells,p)

    def Hplus_HKATP_ecm(self,cells,p):

        # if HKATPase pump is desired, run the H-K-ATPase pump:
        _,_,_,_, f_H2, f_K2 =\
        pumpHKATP(self.cc_cells[self.iH][cells.mem_to_cells],self.cc_ecm[self.iH][cells.mem_to_ecm],
            self.cc_cells[self.iK][cells.mem_to_cells],self.cc_ecm[self.iK][cells.mem_to_ecm],
            cells.mem_sa,cells.cell_vol[cells.mem_to_cells],cells.ecm_vol[cells.mem_to_ecm],
            self.vm,self.T,p,self.HKATP_block)

        H_cell_to = copy.deepcopy(self.cc_cells[self.iH])     # keep track of original H+ in cell and ecm
        H_ecm_to = copy.deepcopy(self.cc_ecm[self.iH])

        # update the concentrations of H+ and K+ in the cell and ecm
        self.update_C_ecm(self.iH,f_H2,cells,p)
        self.update_C_ecm(self.iK,f_K2,cells,p)

         # buffer what's happening with H+ flux to or from the cell and environment:
        delH_cell = self.cc_cells[self.iH] - H_cell_to    # relative change in H wrt the cell
        delH_ecm = self.cc_ecm[self.iH] - H_ecm_to    # relative change in H wrt to environment

        self.cc_cells[self.iH], self.cc_cells[self.iM], self.cHM_cells = bicarbBuffer(
            self.cc_cells[self.iH],self.cc_cells[self.iM],self.cHM_cells,delH_cell,p)

        self.cc_ecm[self.iH], self.cc_ecm[self.iM], self.cHM_ecm = bicarbBuffer(
            self.cc_ecm[self.iH],self.cc_ecm[self.iM],self.cHM_ecm,delH_ecm,p)

        # recalculate the net, unbalanced charge and voltage in each cell:
        self.update_V_ecm(cells,p)

    def Hplus_VATP_ecm(self,cells,p):

        # if HKATPase pump is desired, run the H-K-ATPase pump:
        _, _, f_H3 =\
        pumpVATP(self.cc_cells[self.iH][cells.mem_to_cells],self.cc_ecm[self.iH][cells.mem_to_ecm],
            cells.mem_sa,cells.cell_vol[cells.mem_to_cells],cells.ecm_vol[cells.mem_to_ecm],
            self.vm,self.T,p,self.VATP_block)

        H_cell_to = copy.deepcopy(self.cc_cells[self.iH])     # keep track of original H+ in cell and ecm
        H_ecm_to = copy.deepcopy(self.cc_ecm[self.iH])

        # update the concentration in the cytoplasm and ecms
        self.update_C_ecm(self.iH,f_H3,cells,p)

        # buffer what's happening with H+ flux to or from the cell and environment:
        delH_cell = self.cc_cells[self.iH] - H_cell_to    # relative change in H wrt the cell
        delH_ecm = self.cc_ecm[self.iH] - H_ecm_to    # relative change in H wrt to environment

        self.cc_cells[self.iH], self.cc_cells[self.iM], self.cHM_cells = bicarbBuffer(
            self.cc_cells[self.iH],self.cc_cells[self.iM],self.cHM_cells,delH_cell,p)

        self.cc_ecm[self.iH], self.cc_ecm[self.iM], self.cHM_ecm = bicarbBuffer(
            self.cc_ecm[self.iH],self.cc_ecm[self.iM],self.cHM_ecm,delH_ecm,p)

        # recalculate the net, unbalanced charge and voltage in each cell:
        self.update_V_ecm(cells,p)

def diffuse(cA,cB,Dc,d,sa,vola,volb,p):
    """
    Returns updated concentrations for diffusion between two
    connected volumes.

    Parameters
    ----------
    cA          Initial concentration of c in region A [mol/m3]
    cB          Initial concentration of c in region B [mol/m3]
    Dc          Diffusion constant of c  [m2/s]
    d           Distance between region A and region B [m]
    sa          Surface area separating region A and B [m2]
    vola        volume of region A [m3]
    volb        volume of region B [m3]
    dt          time step   [s]
    method      EULER or RK4 for Euler and Runge-Kutta 4
                integration methods, respectively.

    Returns
    --------
    cA2         Updated concentration of cA in region A [mol/m3]
    cB2         Updated concentration of cB in region B [mol/m3]
    flux        Chemical flux magnitude between region A and B [mol/s]

    """


    flux = -sa*Dc*(cB - cA)/d

    if p.method == 0:

        dmol = sa*p.dt*Dc*(cB - cA)/d

        cA2 = cA + dmol/vola
        cB2 = cB - dmol/volb

    elif p.method == 1:

        k1 = sa*Dc*(cB - cA)/d

        k2 = sa*Dc*(cB - (cA + (1/2)*k1*p.dt))/d

        k3 = sa*Dc*(cB - (cA + (1/2)*k2*p.dt))/d

        k4 = sa*Dc*(cB - (cA + k3*p.dt))/d

        dmol = (p.dt/6)*(k1 + 2*k2 + 2*k3 + k4)

        cA2 = cA + dmol/vola
        cB2 = cB - dmol/volb

    return cA2, cB2, flux

def electrofuse(cA,cB,Dc,d,sa,vola,volb,zc,Vba,T,p):
    """
    Returns updated concentrations for electro-diffusion between two
    connected volumes. Note for cell work, 'b' is 'inside', 'a' is outside, with
    a positive flux moving from a to b. The voltage is defined as
    Vb - Va (Vba), which is equivalent to Vmem.

    This function defaults to regular diffusion if Vba == 0.0

    This function takes numpy matrix values as input. All inputs must be matrices of
    the same shape.

    Parameters
    ----------
    cA          Initial concentration of c in region A [mol/m3] (out)
    cB          Initial concentration of c in region B [mol/m3] (in)
    Dc          Diffusion constant of c  [m2/s]
    d           Distance between region A and region B [m]
    sa          Surface area separating region A and B [m2]
    vola        volume of region A [m3]
    volb        volume of region B [m3]
    zc          valence of ionic species c
    Vba         voltage between B and A as Vb - Va  [V]
    p           an instance of the Parameters class


    Returns
    --------
    cA2         Updated concentration of cA in region A [mol/m3]
    cB2         Updated concentration of cB in region B [mol/m3]
    flux        Chemical flux magnitude between region A and B [mol/s]

    """

    alpha = (zc*Vba*p.F)/(p.R*T)

    #volab = (vola + volb)/2
    #qualityfactor = abs((Dc/d)*(sa/volab)*p.dt*alpha)   # quality factor should be <1.0 for stable simulations

    deno = 1 - np.exp(-alpha)   # calculate the denominator for the electrodiffusion equation,..

    izero = (deno==0).nonzero()     # get the indices of the zero and non-zero elements of the denominator
    inotzero = (deno!=0).nonzero()

    # initialize data matrices to the same shape as input data
    dmol = np.zeros(deno.shape)
    cA2 = np.zeros(deno.shape)
    cB2 = np.zeros(deno.shape)
    flux = np.zeros(deno.shape)

    if p.method == 1:

        k1 = np.zeros(deno.shape)
        k2 = np.zeros(deno.shape)
        k3 = np.zeros(deno.shape)
        k4 = np.zeros(deno.shape)

    if len(deno[izero]):   # if there's anything in the izero array:
         # calculate the flux for those elements [mol/s]:
        flux[izero] = -(sa[izero]/d[izero])*Dc[izero]*(cB[izero] - cA[izero])

        if p.sim_ECM == False:  # if we're simulating extracellular spaces, just calculate the flux

            cA2 = None
            cB2 = None

        elif p.sim_ECM == True:

            if p.method == 0:

                #dmol[izero] = sa[izero]*p.dt*Dc[izero]*(cB[izero] - cA[izero])/d[izero]
                dmol[izero] = -p.dt*flux[izero]

                cA2[izero] = cA[izero] + (dmol[izero]/vola[izero])
                cB2[izero] = cB[izero] - (dmol[izero]/volb[izero])

            elif p.method == 1:

                k1[izero] = -flux[izero]

                k2[izero] = sa[izero]*Dc[izero]*(cB[izero] - (cA[izero] + (1/2)*k1[izero]*p.dt))/d[izero]

                k3[izero] = sa[izero]*Dc[izero]*(cB[izero] - (cA[izero] + (1/2)*k2[izero]*p.dt))/d[izero]

                k4[izero] = sa[izero]*Dc[izero]*(cB[izero] - (cA[izero] + k3[izero]*p.dt))/d[izero]

                dmol[izero] = (p.dt/6)*(k1 + 2*k2 + 2*k3 + k4)

                cA2[izero] = cA[izero] + dmol[izero]/vola[izero]
                cB2[izero] = cB[izero] - dmol[izero]/volb[izero]


    if len(deno[inotzero]):   # if there's any indices in the inotzero array:

        # calculate the flux for those elements:
        flux[inotzero] = -((sa[inotzero]*Dc[inotzero]*alpha[inotzero])/d[inotzero])*\
                       ((cB[inotzero] - cA[inotzero]*np.exp(-alpha[inotzero]))/deno[inotzero])

        if p.sim_ECM == True:

            cA2 = None
            cB2 = None

        elif p.sim_ECM == False:

            if p.method == 0:

                dmol[inotzero] = -flux[inotzero]*p.dt

                cA2[inotzero] = cA[inotzero] + (dmol[inotzero]/vola[inotzero])
                cB2[inotzero] = cB[inotzero] - (dmol[inotzero]/volb[inotzero])


            elif p.method == 1:

                k1[inotzero] = -flux[inotzero]

                k2[inotzero] = ((sa[inotzero]*Dc[inotzero]*alpha[inotzero])/d[inotzero])*\
                             (cB[inotzero] - (cA[inotzero] + (1/2)*k1[inotzero]*p.dt)*np.exp(-alpha[inotzero]))/deno[inotzero]

                k3[inotzero] = ((sa[inotzero]*Dc[inotzero]*alpha[inotzero])/d[inotzero])*\
                             (cB[inotzero] - (cA[inotzero] + (1/2)*k2[inotzero]*p.dt)*np.exp(-alpha[inotzero]))/deno[inotzero]

                k4[inotzero] = ((sa[inotzero]*Dc[inotzero]*alpha[inotzero])/d[inotzero])*\
                             (cB[inotzero] - (cA[inotzero] + k3[inotzero]*p.dt)*np.exp(-alpha[inotzero]))/deno[inotzero]

                dmol[inotzero] = (p.dt/6)*(k1[inotzero] + 2*k2[inotzero] + 2*k3[inotzero] + k4[inotzero])

                cA2[inotzero] = cA[inotzero] + dmol[inotzero]/vola[inotzero]
                cB2[inotzero] = cB[inotzero] - dmol[inotzero]/volb[inotzero]


    return cA2, cB2, flux

def pumpNaKATP(cNai,cNao,cKi,cKo,sa,voli,volo,Vm,T,p,block):

    """
    Parameters
    ----------
    cNai            Concentration of Na+ inside the cell
    cNao            Concentration of Na+ outside the cell
    cKi             Concentration of K+ inside the cell
    cKo             Concentration of K+ outside the cell
    voli            Volume of the cell [m3]
    volo            Volume outside the cell [m3]
    Vm              Voltage across cell membrane [V]
    p               An instance of Parameters object


    Returns
    -------
    cNai2           Updated Na+ inside cell
    cNao2           Updated Na+ outside cell
    cKi2            Updated K+ inside cell
    cKo2            Updated K+ outside cell
    f_Na            Na+ flux (into cell +)
    f_K             K+ flux (into cell +)
    """
    deltaGATP = 20*p.R*T

    delG_Na = p.R*T*np.log(cNao/cNai) - p.F*Vm
    delG_K = p.R*T*np.log(cKi/cKo) + p.F*Vm
    delG_NaKATP = deltaGATP - (3*delG_Na + 2*delG_K)
    delG_pump = (delG_NaKATP/1000)
    delG = np.absolute(delG_pump)
    signG = np.sign(delG)

    if p.backward_pumps == False:

        alpha = sa*block*p.alpha_NaK*step(delG,p.halfmax_NaK,p.slope_NaK)

        f_Na  = -alpha*cNai*cKo      #flux as [mol/s]   scaled to concentrations Na in and K out

    elif p.backward_pumps == True:

        alpha = sa*signG*block*p.alpha_NaK*step(delG,p.halfmax_NaK,p.slope_NaK)

        truth_forwards = signG == 1    # boolean array tagging forward-running pump cells
        truth_backwards = signG == -1  # boolean array tagging backwards-running pump cells

        inds_forwards = (truth_forwards).nonzero()  # indices of forward-running cells
        inds_backwards = (truth_backwards).nonzero() # indices of backward-running cells

        f_Na = np.zeros(len(cNai))

        f_Na[inds_forwards]  = -alpha*cNai*cKo      #flux as [mol/s]   scaled to concentrations Na in and K out

        f_Na[inds_backwards]  = -alpha*cNao*cKi      #flux as [mol/s]   scaled to concentrations K in and Na out

    f_K = -(2/3)*f_Na          # flux as [mol/s]

    if p.sim_ECM == True: # if extracellular spaces are included, concentrations need to be updated by matrix math

        cNai2 = None
        cNao2 = None
        cKi2 = None
        cKo2 = None

    elif p.sim_ECM == False:

        if p.method == 0:

            dmol = f_Na*p.dt

            cNai2 = cNai + dmol/voli
            cNao2 = cNao - dmol/volo

            cKi2 = cKi - (2/3)*dmol/voli
            cKo2 = cKo + (2/3)*dmol/volo

        elif p.method == 1:

            k1 = alpha*cNai*cKo

            k2 = alpha*(cNai+(1/2)*k1*p.dt)*cKo

            k3 = alpha*(cNai+(1/2)*k2*p.dt)*cKo

            k4 = alpha*(cNai+ k3*p.dt)*cKo

            dmol = (p.dt/6)*(k1 + 2*k2 + 2*k3 + k4)

            cNai2 = cNai - dmol/voli
            cNao2 = cNao + dmol/volo

            cKi2 = cKi + (2/3)*dmol/voli
            cKo2 = cKo - (2/3)*dmol/volo


    return cNai2,cNao2,cKi2,cKo2, f_Na, f_K

def pumpCaATP(cCai,cCao,sa,voli,volo,Vm,T,p,block):

    """
    Parameters
    ----------
    cCai            Concentration of Ca2+ inside the cell
    cCao            Concentration of Ca2+ outside the cell
    voli            Volume of the cell [m3]
    volo            Volume outside the cell [m3]
    Vm              Voltage across cell membrane [V]
    p               An instance of Parameters object


    Returns
    -------
    cCai2           Updated Ca2+ inside cell
    cCao2           Updated Ca2+ outside cell
    f_Ca            Ca2+ flux (into cell +)
    """

    deltaGATP = 20*p.R*T

    delG_Ca = p.R*T*np.log(cCao/cCai) - 2*p.F*Vm
    delG_CaATP = deltaGATP - (delG_Ca)
    delG_pump = (delG_CaATP/1000)
    delG = np.absolute(delG_pump)
    signG = np.sign(delG_pump)

    if p.backward_pumps == False:

        alpha = sa*block*p.alpha_Ca*step(delG,p.halfmax_Ca,p.slope_Ca)

        f_Ca  = -alpha*(cCai)      #flux as [mol/s], scaled to concentration in cell

    elif p.backward_pumps == True:

        alpha = sa*signG*block*p.alpha_Ca*step(delG,p.halfmax_Ca,p.slope_Ca)

        truth_forwards = signG == 1
        truth_backwards = signG == -1

        inds_forwards = (truth_forwards).nonzero()  # indices of forward-running cells
        inds_backwards = (truth_backwards).nonzero() # indices of backward-running cells

        f_Ca = np.zeros(len(cCai))

        f_Ca[inds_forwards]  = -alpha*(cCai)      #flux as [mol/s], scaled to concentration in cell
        f_Ca[inds_backwards]  = -alpha*(cCao)      #flux as [mol/s], scaled to concentration out of cell

    if p.sim_ECM == True: # if extracellular spaces are included, concentrations need to be updated by matrix math

        cCai2 = None
        cCao2 = None

    elif p.sim_ECM == False:

        if p.method == 0:

            dmol = f_Ca*p.dt

            cCai2 = cCai + dmol/voli
            cCao2 = cCao - dmol/volo

        elif p.method == 1:

            k1 = alpha*cCai

            k2 = alpha*(cCai+(1/2)*k1*p.dt)

            k3 = alpha*(cCai+(1/2)*k2*p.dt)

            k4 = alpha*(cCai+ k3*p.dt)

            dmol = (p.dt/6)*(k1 + 2*k2 + 2*k3 + k4)

            cCai2 = cCai - dmol/voli
            cCao2 = cCao + dmol/volo


    return cCai2, cCao2, f_Ca

def pumpCaER(cCai,cCao,sa,voli,volo,Vm,T,p,block):

    # deltaGATP = 20*p.R*T
    #
    # delG_Ca = p.R*T*np.log(cCao/cCai) - 2*p.F*Vm
    # delG_CaATP = deltaGATP - (delG_Ca)
    # delG_pump = (delG_CaATP/1000)
    # delG = np.absolute(delG_pump)
    # signG = np.sign(delG_pump)
    #
    #
    # alpha = sa*signG*block*p.alpha_Ca*step(delG,p.halfmax_Ca,p.slope_Ca)
    #
    # truth_forwards = signG == 1
    # truth_backwards = signG == -1
    #
    # inds_forwards = (truth_forwards).nonzero()  # indices of forward-running cells
    # inds_backwards = (truth_backwards).nonzero() # indices of backward-running cells
    #
    # f_Ca = np.zeros(len(cCai))
    #
    # f_Ca[inds_forwards]  = -alpha*(cCai)      #flux as [mol/s], scaled to concentration in cell
    # f_Ca[inds_backwards]  = -alpha*(cCao)      #flux as [mol/s], scaled to concentration out of cell

    alpha = sa*p.alpha_CaER

    f_Ca  = alpha*(cCao)*(1.0 - cCai)      #flux as [mol/s]

    if p.method == 0:

        dmol = f_Ca*p.dt

        cCai2 = cCai + dmol/voli
        cCao2 = cCao - dmol/volo

    elif p.method == 1:

        k1 = alpha*cCao

        k2 = alpha*(cCao+(1/2)*k1*p.dt)

        k3 = alpha*(cCao+(1/2)*k2*p.dt)

        k4 = alpha*(cCao+ k3*p.dt)

        dmol = (p.dt/6)*(k1 + 2*k2 + 2*k3 + k4)

        cCai2 = cCai + dmol/voli
        cCao2 = cCao - dmol/volo


    return cCai2, cCao2, f_Ca

def pumpHKATP(cHi,cHo,cKi,cKo,sa,voli,volo,Vm,T,p,block):

    """
    Parameters
    ----------
    cHi            Concentration of H+ inside the cell
    cHo            Concentration of H+ outside the cell
    cKi             Concentration of K+ inside the cell
    cKo             Concentration of K+ outside the cell
    voli            Volume of the cell [m3]
    volo            Volume outside the cell [m3]
    Vm              Voltage across cell membrane [V]
    p               An instance of Parameters object


    Returns
    -------
    cHi2           Updated H+ inside cell
    cHo2           Updated H+ outside cell
    cKi2            Updated K+ inside cell
    cKo2            Updated K+ outside cell
    f_Na            Na+ flux (into cell +)
    f_K             K+ flux (into cell +)
    """

    deltaGATP = 20*p.R*T

    delG_H = p.R*T*np.log(cHo/cHi) - p.F*Vm
    delG_K = p.R*T*np.log(cKi/cKo) + p.F*Vm

    delG_HKATP = deltaGATP - (delG_H + delG_K)
    delG_pump = (delG_HKATP/1000)
    delG = np.absolute(delG_pump)
    signG = np.sign(delG)

    if p.backward_pumps == False:

        alpha = sa*block*p.alpha_HK*step(delG,p.halfmax_HK,p.slope_HK)
        f_H  = -alpha*cHi*cKo      #flux as [mol/s], scaled by concentrations in and out

    elif p.backward_pumps == True:

        alpha = sa*signG*block*p.alpha_HK*step(delG,p.halfmax_HK,p.slope_HK)

        truth_forwards = signG == 1
        truth_backwards = signG == -1

        inds_forwards = (truth_forwards).nonzero()  # indices of forward-running cells
        inds_backwards = (truth_backwards).nonzero() # indices of backward-running cells

        f_H = np.zeros(len(cHi))

        f_H[inds_forwards]  = -alpha*cHi*cKo      #flux as [mol/s], scaled by concentrations in and out
        f_H[inds_backwards]  = -alpha*cHo*cKi

    f_K = -f_H          # flux as [mol/s]

    if p.sim_ECM == True:

        cHi2 = None
        cHo2 = None
        cKi2 = None
        cKo2 = None

    elif p.sim_ECM == False:

        if p.method == 0:

            dmol = f_H*p.dt

            cHi2 = cHi + dmol/voli
            cHo2 = cHo - dmol/volo

            cKi2 = cKi - dmol/voli
            cKo2 = cKo + dmol/volo

        elif p.method == 1:

            k1 = alpha*cHi*cKo

            k2 = alpha*(cHi+(1/2)*k1*p.dt)*cKo

            k3 = alpha*(cHi+(1/2)*k2*p.dt)*cKo

            k4 = alpha*(cHi+ k3*p.dt)*cKo

            dmol = (p.dt/6)*(k1 + 2*k2 + 2*k3 + k4)

            cHi2 = cHi - dmol/voli
            cHo2 = cHo + dmol/volo

            cKi2 = cKi + dmol/voli
            cKo2 = cKo - dmol/volo


    return cHi2,cHo2,cKi2,cKo2, f_H, f_K

def pumpVATP(cHi,cHo,sa,voli,volo,Vm,T,p,block):

    deltaGATP = 20*p.R*T

    delG_H = p.R*T*np.log(cHo/cHi) - p.F*Vm  # free energy to move H+ out of cell

    delG_VATP = deltaGATP - delG_H   # free energy available to move H+ out of cell
    delG_pump = (delG_VATP/1000)
    delG = np.absolute(delG_pump)
    signG = np.sign(delG)

    if p.backward_pumps == False:

        alpha = sa*block*p.alpha_V*step(delG,p.halfmax_V,p.slope_V)
        f_H  = -alpha*cHi      #flux as [mol/s], scaled by concentrations in and out

    elif p.backward_pumps == True:

        alpha = sa*signG*block*p.alpha_V*step(delG,p.halfmax_V,p.slope_V)

        truth_forwards = signG == 1
        truth_backwards = signG == -1

        inds_forwards = (truth_forwards).nonzero()  # indices of forward-running cells
        inds_backwards = (truth_backwards).nonzero() # indices of backward-running cells

        f_H = np.zeros(len(cHi))

        f_H[inds_forwards]  = -alpha*cHi      #flux as [mol/s], scaled by concentrations in and out
        f_H[inds_backwards]  = -alpha*cHo

    if p.sim_ECM == True:

        cHi2 = None
        cHo2 = None

    elif p.sim_ECM == False:

        if p.method == 0:

            dmol = f_H*p.dt

            cHi2 = cHi + dmol/voli
            cHo2 = cHo - dmol/volo

        elif p.method == 1:

            k1 = alpha*cHi

            k2 = alpha*(cHi+(1/2)*k1*p.dt)

            k3 = alpha*(cHi+(1/2)*k2*p.dt)

            k4 = alpha*(cHi+ k3*p.dt)

            dmol = (p.dt/6)*(k1 + 2*k2 + 2*k3 + k4)

            cHi2 = cHi - dmol/voli
            cHo2 = cHo + dmol/volo

    return cHi2, cHo2, f_H

def bicarbBuffer(cH,cM,cHM,delH,p):

    cM2 = cM - delH

    cHM2 = 0.03*p.CO2

    pH = 6.1 + np.log10(cM/cHM2)

    cH2 = 10**(-pH)*1000

    return cH2, cM2, cHM2

def get_volt(q,sa,p):

    """
    Calculates the voltage for a net charge on a capacitor.

    Parameters
    ----------
    q           Net electrical charge [C]
    sa          Surface area [m2]

    Returns
    -------
    V               Voltage on the capacitive space holding charge

    """

    cap = sa*p.cm
    V = (1/cap)*q
    return V

def get_volt_ECM(q_cells,q_ecm,cells):

    """
    Makes use of the cell <--> ecm Maxwell Capacitance Matrix defined in World module
    to calculate the distinct voltage of the cell interior, its extracellular space, and the
    voltage difference across the membrane. The cell interior and the extracellular space
    are taken to be two conductors (each with self-capacitance) connected by the capacitor of
    the plasma membrane.

    Parameters
    -----------------
    q_cells             An array storing the net charge inside each cell space
    q_ecm               An array storing the net charge in each extracellular space
    cells               An object storing World module properties

    Returns
    -----------------
    v_mem               The voltage across the membrane as Vcell - Vecm
    v_cells             The voltage in the intracellular space at the membrane
    v_ecm              The voltage in the extracellular space at the membrane

    """

    v_cells = cells.Cinv_a*q_cells[cells.mem_to_cells] + cells.Cinv_b*q_ecm[cells.mem_to_ecm]
    v_ecm = cells.Cinv_c*q_cells[cells.mem_to_cells] + cells.Cinv_d*q_ecm[cells.mem_to_ecm]
    v_mem = v_cells - v_ecm

    return v_mem, v_cells, v_ecm

def get_charge(concentrations,zs,vol,p):
    """
    Calculates the total charge in a space given a set of concentrations
    and their ionic charges, along with the space volume.

    Parameters
    ---------------
    concentrations       An array of arrays listing concentrations of different ions in multiple spaces [mol/m2].
    zs                   An array of ion charge states
    vol                  Volume of the spaces [m3]
    p                    Instance of Parameters module object

    Returns
    ----------------
    netcharge            The unbalanced charge in the space or spaces [C]

    """

    q = 0

    for conc,z in zip(concentrations,zs):
        q = q+ conc*z

    netcharge = p.F*q*vol

    return netcharge

def get_molarity(concentrations,p):

    q = 0

    for conc in concentrations:
        q = q+ conc

    netmolarity = q

    return netmolarity

def cell_ave(cells,vm_at_mem):

    """
    Averages Vmem over membrane domains to return a mean value for each cell

    Parameters
    ----------
    cells               An instance of the World module cells object
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

def check_v(vm):
    """
    Does a quick check on Vmem values
    and displays error warning or exception if the value
    indicates the simulation is unstable.

    """
    # if isinstance(cA,np.float64):  # if we just have a singular value
    #     if cA < 0.0:
    #         cA2 = 0.0
    #
    # elif isinstance(cA,np.ndarray): # if we have matrix data
    highvar = np.std(vm)
    isnans = np.isnan(vm)

    if highvar > 50e-3:  # if there's anything in the isubzeros matrix...
        print("Your simulation appears to have become unstable. Please try a smaller time step to validate "
              "accuracy of the results.")

    if isnans.any():  # if there's anything in the isubzeros matrix...
        raise BetseExceptionSimulation("Your simulation has become unstable. Please try a smaller time step,"
                                       "reduce gap junction radius, and/or reduce pump rate coefficients.")

def sigmoid(x,g,y_sat):
    """
    A sigmoidal function (logistic curve) allowing user
    to specify a saturation level (y_sat) and growth rate (g).

    Parameters
    ----------
    x            Input values, may be numpy array or float
    g            Growth rate
    y_sat        Level at which growth saturates

    Returns
    --------
    y            Numpy array or float of values

    """
    y = (y_sat*np.exp(g*x))/(y_sat + (np.exp(g*x)-1))
    return y

def hill(x,K,n):

    """
    The Hill equation (log-transformed sigmoid). Function ranges
    from y = 0 to +1.

    Parameters
    ----------
    x            Input values, may be numpy array or float. Note all x>0 !
    K            Value of x at which curve is 1/2 maximum (y=0.5)
    n            Hill co-efficient n<1 negative cooperativity, n>1 positive.

    Returns
    --------
    y            Numpy array or float of values

    """
    # assert x.all() > 0

    y = x**n/((K**n)+(x**n))

    return y

def step(t,t_on,t_change):
    """
    A step function (bounded by 0 and 1) based on a logistic curve
    and allowing user to specify time for step to come on (t_on) and time for
    change from zero to one to happen.

    Parameters
    ----------
    t            Input values, may be numpy array or float
    t_on         Time step turns on
    t_change     Time for change from 0 to 1 (off to on)

    Returns
    --------
    y            Numpy array or float of values

    """
    g = (1/t_change)*10
    y = 1/(1 + (np.exp(-g*(t-t_on))))
    return y

def pulse(t,t_on,t_off,t_change):
    """
    A pulse function (bounded by 0 and 1) based on logistic curves
    and allowing user to specify time for step to come on (t_on) and time for
    change from zero to one to happen, and time for step to come off (t_change).

    Parameters
    ----------
    t            Input values, may be numpy array or float
    t_on         Time step turns on
    t_off        Time step turns off
    t_change     Time for change from 0 to 1 (off to on)

    Returns
    --------
    y            Numpy array or float of values

    """
    g = (1/t_change)*10
    y1 = 1/(1 + (np.exp(-g*(t-t_on))))
    y2 = 1/(1 + (np.exp(-g*(t-t_off))))
    y = y1 - y2
    return y

def H(x):

    y = 0.5*(np.sign(x) +1)

    return y



