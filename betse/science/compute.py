#!/usr/bin/env python3
# Copyright 2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

# FIXME manage the H+ stuff...
# FIXME straighten out the ER dynamics...
# FIXME do runSim_ECM for cut cell dynamics and proper pickling...


import numpy as np
import numpy.ma as ma
from scipy import interpolate as interp
from scipy import ndimage
import os, os.path
import copy
from random import shuffle
from betse.science import filehandling as fh
from betse.science import visualize as viz
from betse.science import toolbox as tb
from betse.science.dynamics import Dynamics
from betse.science.finitediff import laplacian, gradient
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
        self.z_array = []
        self.z_er = []  # ion valence states of er ions
        self.z_array_er = []
        self.Dm_cells = []              # membrane diffusion constants initialized
        self.D_free = []                 # a list of single-valued free diffusion constants for each ion
        self.D_gj = []                 # a list of single-valued free diffusion constants for each ion
        self.Dm_er = []                  # a list of endoplasmic reticulum membrane states
        self.movingIons = []            # moving ions indices
        self.ionlabel = {}              # dictionary to hold ion label names

        self.T = p.T                # set the base temperature for the simulation

        i = -1                           # an index to track place in ion list

        self.flx_gj_i = np.zeros(len(cells.nn_i))
        self.fluxes_gj_x = []
        self.fluxes_gj_y = []

        self.I_gj_x =np.zeros(len(cells.nn_i))     # total current in the network
        self.I_gj_y =np.zeros(len(cells.nn_i))     # total current in the network

        # Membrane current data structure initialization
        self.flx_mem_i = np.zeros(len(cells.mem_i))
        self.fluxes_mem = []

        self.I_mem =np.zeros(len(cells.mem_i))     # total current across membranes
        self.I_mem_time = []                            # membrane current unit time

        ion_names = list(p.ions_dict.keys())

        #-------------------------------------------------------------------------------------------------------------

        for name in ion_names:

            if p.ions_dict[name] == 1:

                if name != 'H':

                    i = i+1

                    str1 = 'i' + name  # create the ion index

                    setattr(self,str1,i)  # dynamically add this field to the object

                    self.ionlabel[vars(self)[str1]] = p.ion_long_name[name]

                    if name != 'P':
                        self.movingIons.append(vars(self)[str1])

                    # cell concentration for the ion
                    str_cells = 'c' + name + '_cells'
                    setattr(self,str_cells,np.zeros(len(cells.cell_i)))
                    vars(self)[str_cells][:]=p.cell_concs[name]

                    # environmental concentration for the ion
                    str_env = 'c' + name + '_env'
                    setattr(self,str_env,np.zeros(len(cells.cell_i)))
                    vars(self)[str_env][:] = p.env_concs[name]

                    # base membrane permeability for each ion
                    str_Dm = 'Dm' + name

                    setattr(self, str_Dm, np.zeros(len(cells.cell_i)))
                    vars(self)[str_Dm][:] = p.mem_perms[name]

                    # gj diffusion for each ion
                    str_Dgj = 'Dgj' + name

                    setattr(self, str_Dgj, np.zeros(len(cells.nn_i)))
                    vars(self)[str_Dgj][:] = p.free_diff[name]

                    str_z = 'z' + name

                    setattr(self, str_z, np.zeros(len(cells.cell_i)))
                    vars(self)[str_z][:] = p.ion_charge[name]

                    self.cc_cells.append(vars(self)[str_cells])
                    self.cc_env.append(vars(self)[str_env])
                    self.zs.append(p.ion_charge[name])
                    self.z_array.append(vars(self)[str_z])
                    self.Dm_cells.append(vars(self)[str_Dm])
                    self.D_free.append(p.free_diff[name])
                    self.D_gj.append(vars(self)[str_Dgj])

                    self.fluxes_gj_x.append(self.flx_gj_i)
                    self.fluxes_gj_y.append(self.flx_gj_i)
                    self.fluxes_mem.append(self.flx_mem_i)

                    if name == 'Ca':
                        self.cCa_er = np.zeros(len(cells.cell_i))
                        self.cCa_er[:]=p.cCa_er

                        self.zCa_er = np.zeros(len(cells.cell_i))
                        self.zCa_er[:]=p.z_Ca

                        self.cc_er.append(self.cCa_er)
                        self.z_er.append(p.z_Ca)
                        self.z_array_er.append(self.zCa_er)
                        self.Dm_er.append(p.Dm_Ca)

                    if name == 'M' and p.ions_dict['Ca'] == 1:

                        self.cM_er = np.zeros(len(cells.cell_i))
                        self.cM_er[:]=p.cM_er

                        self.zM_er = np.zeros(len(cells.cell_i))
                        self.zM_er[:]=p.z_M

                        self.cc_er.append(self.cM_er)
                        self.z_er.append(p.z_M)
                        self.z_array_er.append(self.zM_er)
                        self.Dm_er.append(p.Dm_M)

        # Do H+ separately as it's complicated by the buffer
        # initialize the carbonic acid for the carbonate buffer

        if p.ions_dict['H'] == 1:

            i = i+ 1

            self.iH = i

            self.ionlabel[self.iH] = 'protons'

            self.movingIons.append(self.iH)

            self.cHM_cells = np.zeros(len(cells.cell_i))
            self.cHM_cells[:] = 0.03*p.CO2

            self.cHM_env = np.zeros(len(cells.cell_i))
            self.cHM_env[:] = 0.03*p.CO2

            self.cH_cells = np.zeros(len(cells.cell_i))
            # cH_cells[:]=p.cH_cell
            self.pH_cell = 6.1 + np.log10(self.cM_cells/self.cHM_cells)
            self.cH_cells = (10**(-self.pH_cell))*1000  # units mmol/L

            self.cH_env = np.zeros(len(cells.cell_i))
            # cH_env[:]=p.cH_env
            self.pH_env = 6.1 + np.log10(self.cM_env/self.cHM_env)
            self.cH_env = (10**(-self.pH_env))*1000 # units mmol/L

            DmH = np.zeros(len(cells.cell_i))
            DmH[:] = p.Dm_H

            DgjH = np.zeros(len(cells.nn_i))
            DgjH[:] = p.free_diff[name]

            self.zH = np.zeros(len(cells.cell_i))
            self.zH[:] = p.z_H

            self.cc_cells.append(self.cH_cells)
            self.cc_env.append(self.cH_env)
            self.zs.append(p.z_H)
            self.z_array.append(self.zH)
            self.Dm_cells.append(DmH)
            self.D_gj.append(DgjH)
            self.D_free.append(p.Do_H)

            self.fluxes_gj_x.append(self.flx_gj_i)
            self.fluxes_gj_y.append(self.flx_gj_i)
            self.fluxes_mem.append(self.flx_mem_i)

        #-------------------------------------------------------------------------------------------------------

        # Initialize membrane thickness:
        self.tm = np.zeros(len(cells.cell_i))
        self.tm[:] = p.tm

        # Initialize environmental volume:
        self.envV = np.zeros(len(cells.cell_i))
        self.envV[:] = p.vol_env

        self.Dm_er = self.Dm_cells[:]

        self.vm_to = np.zeros(len(cells.cell_i))
        self.v_er = np.zeros(len(cells.cell_i))

        if p.global_options['NaKATP_block'] != 0:

            self.NaKATP_block = np.ones(len(cells.cell_i))  # initialize NaKATP blocking vector

        else:
            self.NaKATP_block = 1

        if p.HKATPase_dyn == True and p.global_options['HKATP_block'] != 0:
            self.HKATP_block = np.ones(len(cells.cell_i))  # initialize HKATP blocking vector
        else:
            self.HKATP_block = 1

        # add channel noise to the model:
        self.channel_noise_factor = np.random.random(len(cells.cell_i))
        self.Dm_cells[self.iK] = (p.channel_noise_level*self.channel_noise_factor + 1)*self.Dm_cells[self.iK]

        if p.dynamic_noise == True:
            # add a random walk on protein concentration to generate dynamic noise:
            self.protein_noise_factor = p.dynamic_noise_level*(np.random.random(len(cells.cell_i)) - 0.5)

            if p.ions_dict['P']==1:
                self.cc_cells[self.iP] = self.cc_cells[self.iP]*(1+ self.protein_noise_factor)

        self.fluxes_mem = np.asarray(self.fluxes_mem)

        # Initialize Dye and IP3

        if p.scheduled_options['IP3'] != 0 or p.Ca_dyn == True:

            self.cIP3 = np.zeros(len(cells.cell_i))  # initialize a vector to hold IP3 concentrations
            self.cIP3[:] = p.cIP3_to                 # set the initial concentration of IP3 from params file

            self.cIP3_flux_gj = np.zeros(len(cells.nn_i))
            self.cIP3_flux_mem = np.zeros(len(cells.mem_i))

            self.cIP3_flux_gj_time = []
            self.cIP3_flux_mem_time = []

            if p.sim_ECM == True:
                self.cIP3_ecm = np.zeros(len(cells.ecm_i))     # initialize IP3 concentration of the environment
                self.cIP3_ecm[:] = p.cIP3_to_env
                self.cIP3_env = np.zeros(len(cells.env_i))
                self.cIP3_env[:] = p.cIP3_to_env

                self.cIP3_flux_ecm = np.zeros(len(cells.ecm_i))
                self.cIP3_flux_env = np.zeros(len(cells.env_i))
                self.cIP3_flux_ecm_time = []
                self.cIP3_flux_env_time = []

            elif p.sim_ECM == False:
                self.cIP3_env = np.zeros(len(cells.cell_i))     # initialize IP3 concentration of the environment
                self.cIP3_env[:] = p.cIP3_to_env

        if p.voltage_dye == True:

            self.cDye_cell = np.zeros(len(cells.cell_i))   # initialize voltage sensitive dye array for cell and env't
            self.cDye_cell[:] = p.cDye_to_cell

            self.Dye_flux_gj = np.zeros(len(cells.nn_i))
            self.Dye_flux_mem = np.zeros(len(cells.mem_i))


            if p.sim_ECM == True:
                self.cDye_ecm = np.zeros(len(cells.ecm_i))
                self.cDye_ecm[:] = p.cDye_to
                self.cDye_env = np.zeros(len(cells.env_i))
                self.cDye_env[:] = p.cDye_to

                self.Dye_flux_ecm = np.zeros(len(cells.ecm_i))
                self.Dye_flux_env = np.zeros(len(cells.env_i))


            else:
                self.cDye_env = np.zeros(len(cells.cell_i))     # initialize Dye concentration in the environment
                self.cDye_env[:] = p.cDye_to

    def baseInit_ECM(self,cells,p):

        self.cc_cells = []  # cell concentrations initialized
        self.cc_er = []   # endoplasmic reticulum ion concentrations in each cell
        self.cc_env = []   # environmental concentrations initialized

        self.v_env = np.zeros(len(cells.xypts))
        self.v_cell = np.zeros(len(cells.cell_i))

        self.zs = []   # ion valence state initialized
        self.z_er = []  # ion valence states of er ions
        self.z_array = []  # ion valence array matched to cell points
        self.z_array_env = []  # ion valence array matched to env points
        self.z_array_er = []
        self.Dm_cells = []              # membrane diffusion constants initialized
        self.Dm_er = []                  # a list of endoplasmic reticulum membrane state
        self.D_free = []                 # a list of single-valued free diffusion constants for each ion
        self.D_env = []                 # an array of diffusion constants for each ion defined on env grid
        self.D_gj = []                  # an array of diffusion constants for gap junctions
        self.movingIons = []            # moving ions indices
        self.ionlabel = {}              # dictionary to hold ion label names
        self.c_env_bound = []           # moving ion concentration at global boundary

        self.T = p.T                # set the base temperature for the simulation

        # Initialize membrane thickness:
        self.tm = np.zeros(len(cells.mem_i))
        self.tm[:] = p.tm

        self.flx_gj_i = np.zeros(len(cells.nn_i))  # flux matrix across gj for individual ions
        self.fluxes_gj_x = []
        self.fluxes_gj_y = []
        self.I_gj_x =np.zeros(len(cells.nn_i))     # total current in the gj network
        self.I_gj_y =np.zeros(len(cells.nn_i))     # total current in the gj network

        # Membrane current data structure initialization
        self.flx_mem_i = np.zeros(len(cells.mem_i))
        self.fluxes_mem = []
        self.I_mem =np.zeros(len(cells.mem_i))     # total current across membranes

        self.flx_env_i = np.zeros(len(cells.xypts))
        self.fluxes_env_x = []
        self.fluxes_env_y = []
        self.I_env =np.zeros(len(cells.xypts))     # total current in environment

        if p.sim_eosmosis == True:
            self.rho_channel = np.ones(len(cells.mem_i))
        else:
            self.rho_channel = 1  # else just define it as identity.

        ion_names = list(p.ions_dict.keys())

        i = -1                           # an index to track place in ion list
        for name in ion_names:

            if p.ions_dict[name] == 1:

                if name != 'H':

                    i = i+1

                    str1 = 'i' + name  # create the ion index

                    setattr(self,str1,i)  # dynamically add this field to the object

                    self.ionlabel[vars(self)[str1]] = p.ion_long_name[name]

                    if name != 'P':
                        self.movingIons.append(vars(self)[str1])

                    # cell concentration for the ion
                    str_cells = 'c' + name + '_cells'
                    setattr(self,str_cells,np.zeros(len(cells.cell_i)))
                    vars(self)[str_cells][:]=p.cell_concs[name]

                    # environmental concentration for the ion
                    str_env = 'c' + name + '_env'
                    setattr(self,str_env,np.zeros(len(cells.xypts)))
                    vars(self)[str_env][:] = p.env_concs[name]

                    # base membrane diffusion for each ion
                    str_Dm = 'Dm' + name

                    setattr(self, str_Dm, np.zeros(len(cells.mem_i)))
                    vars(self)[str_Dm][:] = p.mem_perms[name]

                    # base membrane diffusion for each ion
                    str_Dgj = 'Dgj' + name

                    setattr(self, str_Dgj, np.zeros(len(cells.nn_i)))
                    vars(self)[str_Dgj][:] = p.free_diff[name]

                    # environmental diffusion for each ion
                    str_Denv = 'D' + name

                    setattr(self, str_Denv, np.zeros(len(cells.xypts)))
                    vars(self)[str_Denv][:] = p.free_diff[name]

                    str_z = 'z' + name

                    setattr(self, str_z, np.zeros(len(cells.cell_i)))
                    vars(self)[str_z][:] = p.ion_charge[name]

                    str_z2 = 'z2' + name

                    setattr(self, str_z2, np.zeros(len(cells.xypts)))
                    vars(self)[str_z2][:] = p.ion_charge[name]

                    self.cc_cells.append(vars(self)[str_cells])
                    self.cc_env.append(vars(self)[str_env])
                    self.c_env_bound.append(p.env_concs[name])

                    self.zs.append(p.ion_charge[name])
                    self.z_array.append(vars(self)[str_z])
                    self.z_array_env.append(vars(self)[str_z2])
                    self.Dm_cells.append(vars(self)[str_Dm])
                    self.D_env.append(vars(self)[str_Denv])
                    self.D_gj.append(vars(self)[str_Dgj])
                    self.D_free.append(p.free_diff[name])

                    self.fluxes_gj.append(self.flx_gj_i)
                    self.fluxes_mem.append(self.flx_mem_i)
                    self.fluxes_env_x.append(self.flx_env_i)
                    self.fluxes_env_y.append(self.flx_env_i)

                    if name == 'Ca':
                        self.cCa_er = np.zeros(len(cells.cell_i))
                        self.cCa_er[:]=p.cCa_er

                        self.zCa_er = np.zeros(len(cells.cell_i))
                        self.zCa_er[:]=p.z_Ca

                        self.cc_er.append(self.cCa_er)
                        self.z_er.append(p.z_Ca)
                        self.Dm_er.append(p.Dm_Ca)

                    if name == 'M' and p.ions_dict['Ca'] == 1:

                        self.cM_er = np.zeros(len(cells.cell_i))
                        self.cM_er[:]=p.cM_er

                        self.zM_er = np.zeros(len(cells.cell_i))
                        self.zM_er[:]=p.z_M

                        self.cc_er.append(self.cM_er)
                        self.z_er.append(p.z_M)
                        self.Dm_er.append(p.Dm_M)

        # Do H+ separately as it's complicated by the buffer
        # initialize the carbonic acid for the carbonate buffer

        if p.ions_dict['H'] == 1:

            i = i+ 1

            self.iH = i

            self.ionlabel[self.iH] = 'protons'

            self.movingIons.append(self.iH)

            self.cHM_cells = np.zeros(len(cells.cell_i))
            self.cHM_cells[:] = 0.03*p.CO2

            self.cHM_env = np.zeros(len(cells.xypts))
            self.cHM_env[:] = 0.03*p.CO2

            self.cH_cells = np.zeros(len(cells.cell_i))

            self.pH_cell = 6.1 + np.log10(self.cM_cells/self.cHM_cells)
            self.cH_cells = (10**(-self.pH_cell))*1000  # units mmol/L

            self.cH_env = np.zeros(len(cells.xypts))
            self.pH_env = 6.1 + np.log10(self.cM_env/self.cHM_env)
            self.cH_env = (10**(-self.pH_env))*1000 # units mmol/L

            DmH = np.zeros(len(cells.mem_i))
            DmH[:] = p.Dm_H

            self.zH = np.zeros(len(cells.cell_i))
            self.zH[:] = p.z_H

            self.zH2 = np.zeros(len(cells.xypts))
            self.zH2[:] = p.z_H

            # environmental diffusion for each ion
            DenvH = np.zeros(len(cells.xypts))
            DenvH[:] = p.free_diff['H']

            DgjH = np.zeros(len(cells.nn_i))
            DgjH[:] = p.free_diff['H']

            self.cc_cells.append(self.cH_cells)
            self.cc_env.append(self.cH_env)
            self.zs.append(p.z_H)
            self.z_array.append(self.zH)
            self.z_array_env.append(self.zH2)
            self.Dm_cells.append(DmH)
            self.D_env.append(DenvH)
            self.D_gj.append(DgjH)
            self.D_free.append(p.Do_H)

            self.fluxes_gj.append(self.flx_gj_i)
            self.fluxes_mem.append(self.flx_mem_i)
            self.fluxes_env_x.append(self.flx_env_i)
            self.fluxes_env_y.append(self.flx_env_i)

        #-------------------------------------------------------------------------------------------------------
        # set up the ecm diffusion constants, which may be heterogeneous between cell cluster and free environment:

        for i, dmat in enumerate(self.D_env):
            # for all cells in the cluster, set the internal diffusion constant:
            self.D_env[i][cells.map_cell2ecm] = dmat[cells.map_cell2ecm]*p.D_ecm_mult
            # for all membranes in the cluster, set the internal diffusion constant representing ecm junctions:
            self.D_env[i][cells.map_mem2ecm] = dmat[cells.map_mem2ecm]*p.D_ecm_mult
            # undo the transform for membranes on the boundary (set them back to free-diffusion values):
            self.D_env[i][cells.ecm_bound_k] = dmat[cells.ecm_bound_k]/p.D_ecm_mult


        self.Dm_er = self.Dm_cells[:]

        self.vm_to = np.zeros(len(cells.cell_i))
        self.v_er = np.zeros(len(cells.cell_i))

        if p.global_options['NaKATP_block'] != 0:

            self.NaKATP_block = np.ones(len(cells.cell_i))  # initialize NaKATP blocking vector

        else:
            self.NaKATP_block = 1

        if p.HKATPase_dyn == True and p.global_options['HKATP_block'] != 0:
            self.HKATP_block = np.ones(len(cells.cell_i))  # initialize HKATP blocking vector
        else:
            self.HKATP_block = 1

        # add channel noise to the model:
        self.channel_noise_factor = np.random.random(len(cells.mem_i))
        self.Dm_cells[self.iK] = (p.channel_noise_level*self.channel_noise_factor + 1)*self.Dm_cells[self.iK]

        if p.dynamic_noise == True:
            # add a random walk on protein concentration to generate dynamic noise:
            self.protein_noise_factor = p.dynamic_noise_level*(np.random.random(len(cells.cell_i)) - 0.5)

            if p.ions_dict['P']==1:
                self.cc_cells[self.iP] = self.cc_cells[self.iP]*(1+ self.protein_noise_factor)

        self.cc_cells = np.asarray(self.cc_cells)
        self.cc_env = np.asarray(self.cc_env)
        self.zs = np.asarray(self.zs)
        self.z_array = np.asarray(self.z_array)
        self.z_array_env = np.asarray(self.z_array_env)
        self.Dm_cells = np.asarray(self.Dm_cells)
        self.D_env = np.asarray(self.D_env)
        self.D_free = np.asarray(self.D_free)

        self.fluxes_gj  = np.asarray(self.fluxes_gj)
        self.fluxes_mem  = np.asarray(self.fluxes_mem)
        self.fluxes_env_x = np.asarray(self.fluxes_env_x)
        self.fluxes_env_y = np.asarray(self.fluxes_env_y)

        # boundary conditions -----------------------------------------------------------------------

        # definition of boundary values -- starting vals of these go into the config file and params --
        # scheduled dynamics might vary the values
        self.bound_V = {}
        self.bound_V['T'] = 0
        self.bound_V['B'] = 0
        self.bound_V['L'] = 0
        self.bound_V['R'] = 0

    def tissueInit(self,cells,p):

        if p.sim_ECM == True:
            pass
            #  # re-initialize diffusion constants for the ecm-ecm junctions in case value changed:
            # self.id_ecm = np.ones(len(cells.ecm_nn_i))
            #
            # self.D_ecm_juncs = []
            # for D in self.D_free:
            #     DD = D*self.id_ecm*p.D_ecm_mult
            #     self.D_ecm_juncs.append(DD)
            #
            # self.D_ecm_juncs = np.asarray(self.D_ecm_juncs)
            #
            # # create a v_ecm copy that will let us change the voltage externally:
            # self.v_env_mod = np.zeros(len(cells.env_i))

        self.dyna = Dynamics(self,cells,p)   # create the tissue dynamics object
        self.dyna.tissueProfiles(self,cells,p)  # initialize all tissue profiles

        if p.sim_ECM == True:
            self.dyna.ecmBoundProfiles(self,cells,p) # initialize boundary profiles

        # add channel noise to the model:
        self.Dm_cells[self.iK] = (p.channel_noise_level*self.channel_noise_factor + 1)*self.Dm_cells[self.iK]

        # Initialize an array structure that will hold user-scheduled changes to membrane permeabilities:
        Dm_cellsA = np.asarray(self.Dm_cells)
        Dm_cellsER = np.asarray(self.Dm_er)

        # if tb.emptyDict(p.scheduled_options) == False or tb.emptyDict(p.vg_options) == False or p.Ca_dyn == True:
        self.Dm_base = np.copy(Dm_cellsA) # make a copy that will serve as the unaffected values base

        # if tb.emptyDict(p.scheduled_options) == False:
        self.Dm_scheduled = np.copy(Dm_cellsA)
        self.Dm_scheduled[:] = 0

        # if tb.emptyDict(p.vg_options) == False:
            # Initialize an array structure that will hold dynamic voltage-gated channel changes to mem permeability:
        self.Dm_vg = np.copy(Dm_cellsA)
        self.Dm_vg[:] = 0

        # if p.Ca_dyn == True:
            # Initialize an array structure that will hold dynamic calcium-gated channel changes to mem perms:
        self.Dm_cag = np.copy(Dm_cellsA)
        self.Dm_cag[:] = 0

        self.Dm_er_base = np.copy(Dm_cellsER)

        self.Dm_er_CICR = np.copy(Dm_cellsER)
        self.Dm_er_CICR[:] = 0

        self.dcc_ER = []

        if p.global_options['gj_block'] != 0:

            self.gj_block = np.ones(len(cells.nn_i))   # initialize the gap junction blocking vector to ones

        else:

            self.gj_block = 1

        if p.scheduled_options['IP3'] != 0 or p.Ca_dyn == True:

            self.cIP3 = np.zeros(len(cells.cell_i))  # initialize a vector to hold IP3 concentrations
            self.cIP3[:] = p.cIP3_to                 # set the initial concentration of IP3 from params file

            self.cIP3_flux_gj = np.zeros(len(cells.nn_i))
            self.cIP3_flux_mem = np.zeros(len(cells.mem_i))

            self.cIP3_flux_gj_time = []
            self.cIP3_flux_mem_time = []

            if p.sim_ECM == True:
                # self.cIP3_ecm = np.zeros(len(cells.ecm_i))     # initialize IP3 concentration of the environment
                # self.cIP3_ecm[:] = p.cIP3_to_env
                self.cIP3_env = np.zeros(len(cells.xypts))
                self.cIP3_env[:] = p.cIP3_to_env

                self.cIP3_flux_env_x = np.zeros(len(cells.xypts))
                self.cIP3_flux_env_y = np.zeros(len(cells.xypts))

            elif p.sim_ECM == False:
                self.cIP3_env = np.zeros(len(cells.cell_i))     # initialize IP3 concentration of the environment
                self.cIP3_env[:] = p.cIP3_to_env

        if p.voltage_dye == True:

            self.cDye_cell = np.zeros(len(cells.cell_i))   # initialize voltage sensitive dye array for cell and env't
            self.cDye_cell[:] = p.cDye_to_cell

            self.Dye_flux_gj = np.zeros(len(cells.nn_i))
            self.Dye_flux_mem = np.zeros(len(cells.mem_i))


            if p.sim_ECM == True:
                # self.cDye_ecm = np.zeros(len(cells.ecm_i))
                # self.cDye_ecm[:] = p.cDye_to
                self.cDye_env = np.zeros(len(cells.xypts))
                self.cDye_env[:] = p.cDye_to

                self.Dye_flux_env_x = np.zeros(len(cells.xypts))
                self.Dye_flux_env_y = np.zeros(len(cells.xypts))


            else:
                self.Dye_env = np.zeros(len(cells.cell_i))     # initialize Dye concentration in the environment
                self.Dye_env[:] = p.cDye_to


        self.dyna.globalInit(self,cells,p)     # initialize any global interventions
        self.dyna.scheduledInit(self,cells,p)  # initialize any scheduled interventions
        self.dyna.dynamicInit(self,cells,p)    # initialize any dynamic channel properties

        loggers.log_info('This world contains '+ str(cells.cell_number) + ' cells.')
        loggers.log_info('Each cell has an average of '+ str(round(cells.average_nn,2)) + ' nearest-neighbours.')
        loggers.log_info('You are running the ion profile: '+ p.ion_profile)

        loggers.log_info('Ions in this simulation: ' + str(self.ionlabel))
        loggers.log_info('If you have selected features using other ions, they will be ignored.')

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
        self.cc_er_time = []   # retains er concentrations as a function of time
        self.cIP3_time = []    # retains cellular ip3 concentrations as a function of time

        self.I_mem_time = []    # initialize membrane current time vector

        self.vm_Matrix = []    # initialize matrices for resampled data sets (used in smooth plotting and streamlines)
        self.I_gj_x_time = []
        self.I_gj_y_time = []

        self.efield_gj_x_time = []   # matrices storing smooth electric field in gj connected cells
        self.efield_gj_y_time = []

        if p.voltage_dye == True:

            self.cDye_time = []    # retains voltage-sensitive dye concentration as a function of time
            self.Dye_flux_gj_time = []
            self.Dye_flux_mem_time = []

        # gap junction specific arrays:
        self.id_gj = np.ones(len(cells.nn_i))  # identity array for gap junction indices...
        self.gjopen = np.ones(len(cells.nn_i))   # holds gap junction open fraction for each gj
        self.gjl = np.zeros(len(cells.nn_i))    # gj length for each gj
        self.gjl[:] = cells.nn_len

        # get the net, unbalanced charge and corresponding voltage in each cell:
        q_cells = get_charge(self.cc_cells,self.z_array,cells.cell_vol,p)
        self.vm = get_volt(q_cells,cells.cell_sa,p)

        # vm_to = copy.deepcopy(self.vm)   # create a copy of the original voltage
        self.vm_to = self.vm[:]

        if p.Ca_dyn == True:
            self.cc_er = np.asarray(self.cc_er)
            self.cc_er_to = self.cc_er[:]

        # create a time-steps vector appropriate for the simulation type:
        if p.run_sim == True:
            tt = np.linspace(0,p.sim_tsteps*p.dt,p.sim_tsteps)
        else:
            tt = np.linspace(0,p.init_tsteps*p.dt,p.init_tsteps)

        i = 0 # resample the time vector to save data at specific times:
        tsamples =[]
        resample = p.t_resample
        while i < len(tt)-resample:
            i = i + resample
            tsamples.append(tt[i])
        tsamples = set(tsamples)

        # report

        if p.run_sim == True:


            loggers.log_info('Your simulation is running from '+ str(0) + ' to '+ str(round(p.sim_tsteps*p.dt,3))
                         + ' seconds of in-world time.')

        else:
             loggers.log_info('Your initialization is running from '+ str(0) + ' to '+ str(round(p.init_tsteps*p.dt,3))
                         + ' seconds of in-world time.')


        if p.plot_while_solving == True:

            self.checkPlot = viz.PlotWhileSolving(cells,self,p,clrAutoscale = p.autoscale_Vmem, clrMin = p.Vmem_min_clr,
                clrMax = p.Vmem_max_clr)

        do_once = True  # a variable to time the loop only once

        for t in tt:   # run through the loop

            if do_once == True:
                loop_measure = time.time()

            self.fluxes_mem.fill(0)  # reinitialize flux storage device

            self.dvm = (self.vm - self.vm_to)/p.dt    # calculate the change in the voltage derivative
            # self.vm_to = copy.deepcopy(self.vm)       # reassign the history-saving vm
            self.vm_to = self.vm[:]

            if p.Ca_dyn ==1 and p.ions_dict['Ca'] == 1:

                self.dcc_ER = (self.cc_er - self.cc_er_to)/p.dt
                # self.cc_er_to = copy.deepcopy(self.cc_er)
                self.cc_er_to = self.cc_er[:]

            # calculate the values of scheduled and dynamic quantities (e.g. ion channel multipliers):
            # self.allDynamics(t,p)  # user-scheduled (forced) interventions
            if p.run_sim == True:
                self.dyna.runAllDynamics(self,cells,p,t)

            # run the Na-K-ATPase pump:
            fNa_NaK, fK_NaK = pumpNaKATP(self.cc_cells[self.iNa],self.cc_env[self.iNa],self.cc_cells[self.iK],
                self.cc_env[self.iK],self.vm,self.T,p,self.NaKATP_block)

            # update the concentration in cells (assume environment steady and constant supply of ions)
            self.cc_cells[self.iNa] = self.cc_cells[self.iNa] + fNa_NaK*(cells.cell_sa/cells.cell_vol)*p.dt
            self.cc_cells[self.iK] = self.cc_cells[self.iK] + fK_NaK*(cells.cell_sa/cells.cell_vol)*p.dt

            # store transmembrane fluxes for each ion:
            self.fluxes_mem[self.iNa] = fNa_NaK[cells.mem_to_cells]
            self.fluxes_mem[self.iK] = fK_NaK[cells.mem_to_cells]

            # recalculate the net, unbalanced charge and voltage in each cell:
            q_cells = get_charge(self.cc_cells,self.z_array,cells.cell_vol,p)
            self.vm = get_volt(q_cells,cells.cell_sa,p)

            if p.ions_dict['Ca'] == 1:
                # run the calcium ATPase membrane pump:
                fCaATP = pumpCaATP(self.cc_cells[self.iCa],self.cc_env[self.iCa],self.vm,self.T,p)

                # update concentrations in the cell:
                self.cc_cells[self.iCa] = self.cc_cells[self.iCa] + fCaATP*(cells.cell_sa/cells.cell_vol)*p.dt

                # store the transmembrane flux for this ion
                self.fluxes_mem[self.iCa] = fCaATP[cells.mem_to_cells]

                # recalculate the net, unbalanced charge and voltage in each cell:
                q_cells = get_charge(self.cc_cells,self.z_array,cells.cell_vol,p)
                self.vm = get_volt(q_cells,cells.cell_sa,p)

                if p.Ca_dyn ==1:

                    # run the calcium ATPase endoplasmic reticulum pump:
                    fCaATP_ER = pumpCaER(self.cc_er[0],self.cc_cells[self.iCa],self.v_er,self.T,p)

                    # update calcium concentrations in the ER and cell:
                    self.cc_er[0] = self.cc_er[0] + fCaATP_ER*((cells.cell_sa)/(p.ER_vol*cells.cell_vol))*p.dt
                    self.cc_cells[self.iCa] = self.cc_cells[self.iCa] - fCaATP_ER*(cells.cell_sa/cells.cell_vol)

                    # recalculate the net, unbalanced charge and voltage in each cell:
                    q_cells = get_charge(self.cc_cells,self.z_array,cells.cell_vol,p)
                    self.vm = get_volt(q_cells,cells.cell_sa,p)

                    q_er = get_charge(self.cc_er,self.z_array_er,p.ER_vol*cells.cell_vol,p)
                    v_er_o = get_volt(q_er,p.ER_sa*cells.cell_sa,p)
                    self.v_er = v_er_o - self.vm

            if p.ions_dict['H'] == 1:


                # electrofuse the H+ ion between the cytoplasm and the environment
                f_H1 = electroflux(self.cc_env[self.iH],self.cc_cells[self.iH],self.Dm_cells[self.iH],self.tm,
                    self.zs[self.iH],self.vm,self.T,p)

                self.cc_cells[self.iH] = self.cc_cells[self.iH] + f_H1*(cells.cell_sa/cells.cell_vol)*p.dt
                self.cc_env[self.iH] = self.cc_env[self.iH] - f_H1*(cells.cell_sa/p.vol_env)*p.dt

                self.fluxes_mem[self.iH] = f_H1[cells.mem_to_cells]

                # # buffer what's happening with H+ flux to or from the cell and environment:
                delH_cell = f_H1*(cells.cell_sa/cells.cell_vol)*p.dt    # relative change in H wrt the cell
                delH_env =  -f_H1*(cells.cell_sa/p.vol_env)*p.dt    # relative change in H wrt to environment

                self.cc_cells[self.iH], self.cc_cells[self.iM], self.cHM_cells = bicarbBuffer(
                    self.cc_cells[self.iH],self.cc_cells[self.iM],self.cHM_cells,delH_cell,p)

                self.cc_env[self.iH], self.cc_env[self.iM], self.cHM_env = bicarbBuffer(
                    self.cc_env[self.iH],self.cc_env[self.iM],self.cHM_env,delH_env,p)

                # recalculate the net, unbalanced charge and voltage in each cell:
                q_cells = get_charge(self.cc_cells,self.z_array,cells.cell_vol,p)
                self.vm = get_volt(q_cells,cells.cell_sa,p)

                if p.HKATPase_dyn == 1:

                    # if HKATPase pump is desired, run the H-K-ATPase pump:
                    f_H2, f_K2 = pumpHKATP(self.cc_cells[self.iH],self.cc_env[self.iH],self.cc_cells[self.iK],
                        self.cc_env[self.iK],self.vm,self.T,p,self.HKATP_block)

                    # update the concentration in cells (assume environment steady and constant supply of ions)
                    self.cc_cells[self.iH] = self.cc_cells[self.iH] + f_H2*(cells.cell_sa/cells.cell_vol)*p.dt
                    self.cc_cells[self.iK] = self.cc_cells[self.iK] + f_K2*(cells.cell_sa/cells.cell_vol)*p.dt

                    # store fluxes for this pump:
                    self.fluxes_mem[self.iH] = self.fluxes_mem[self.iH] + f_H2[cells.mem_to_cells]
                    self.fluxes_mem[self.iK] = self.fluxes_mem[self.iK] + f_K2[cells.mem_to_cells]

                     # buffer what's happening with H+ flux to or from the cell and environment:
                    delH_cell = f_H2*(cells.cell_sa/cells.cell_vol)*p.dt    # relative change in H wrt the cell
                    delH_env =  -f_H2*(cells.cell_sa/p.vol_env)*p.dt    # relative change in H wrt to environment

                    self.cc_cells[self.iH], self.cc_cells[self.iM], self.cHM_cells = bicarbBuffer(
                        self.cc_cells[self.iH],self.cc_cells[self.iM],self.cHM_cells,delH_cell,p)

                    self.cc_env[self.iH], self.cc_env[self.iM], self.cHM_env = bicarbBuffer(
                        self.cc_env[self.iH],self.cc_env[self.iM],self.cHM_env,delH_env,p)

                    # recalculate the net, unbalanced charge and voltage in each cell:
                    q_cells = get_charge(self.cc_cells,self.z_array,cells.cell_vol,p)
                    self.vm = get_volt(q_cells,cells.cell_sa,p)

                if p.VATPase_dyn == 1:

                     # if HKATPase pump is desired, run the H-K-ATPase pump:
                    f_H3 = pumpVATP(self.cc_cells[self.iH],self.cc_env[self.iH],self.vm,self.T,p)

                    self.cc_cells[self.iH] = self.cc_cells[self.iH] + f_H3*(cells.cell_sa/cells.cell_vol)*p.dt

                    self.fluxes_mem[self.iH]  = self.fluxes_mem[self.iH] + f_H3[cells.mem_to_cells]

                     # buffer what's happening with H+ flux to or from the cell and environment:
                    delH_cell = f_H3*(cells.cell_sa/cells.cell_vol)*p.dt    # relative change in H wrt the cell
                    delH_env =  -f_H3*(cells.cell_sa/p.vol_env)*p.dt    # relative change in H wrt to environment

                    self.cc_cells[self.iH], self.cc_cells[self.iM], self.cHM_cells = bicarbBuffer(
                        self.cc_cells[self.iH],self.cc_cells[self.iM],self.cHM_cells,delH_cell,p)

                    self.cc_env[self.iH], self.cc_env[self.iM], self.cHM_env = bicarbBuffer(
                        self.cc_env[self.iH],self.cc_env[self.iM],self.cHM_env,delH_env,p)

                    # recalculate the net, unbalanced charge and voltage in each cell:
                    q_cells = get_charge(self.cc_cells,self.z_array,cells.cell_vol,p)
                    self.vm = get_volt(q_cells,cells.cell_sa,p)

            # electro-diffuse all ions (except for proteins, which don't move) across the cell membrane:

            shuffle(self.movingIons)

            for i in self.movingIons:

                # electrodiffusion of ion between cell and extracellular matrix
                f_ED = electroflux(self.cc_env[i],self.cc_cells[i],self.Dm_cells[i],self.tm,self.zs[i],self.vm,self.T,p)

                # update ion due to transmembrane flux:
                self.cc_cells[i] = self.cc_cells[i] + f_ED*(cells.cell_sa/cells.cell_vol)*p.dt

                self.fluxes_mem[i] = self.fluxes_mem[i] + f_ED[cells.mem_to_cells]

                # # recalculate the net, unbalanced charge and voltage in each cell:
                q_cells = get_charge(self.cc_cells,self.z_array,cells.cell_vol,p)
                self.vm = get_volt(q_cells,cells.cell_sa,p)

                self.update_gj(cells,p,t,i)

            # caculate electric fields:

            self.get_Efield(cells, p)

            if p.scheduled_options['IP3'] != 0 or p.Ca_dyn == True:
                # determine flux through gap junctions for IP3:

                fIP3 = electroflux(self.cIP3[cells.nn_i][:,0],self.cIP3[cells.nn_i][:,1],
                    self.id_gj*p.Do_IP3,self.gjl,p.z_IP3,self.vgj,self.T,p)

                # update cell IP3 concentration due to gap junction flux:

                deltac_gj = self.gjopen*(fIP3)*p.dt
                self.cIP3 = self.cIP3 + ((cells.cell_sa*p.gj_surface)/cells.cell_vol)*np.dot(deltac_gj, cells.gjMatrix)

                # electrodiffuse IP3 between cell and environment:
                fIP3_ED = electroflux(self.cIP3_env,self.cIP3,p.Dm_IP3*self.id_cells,self.tm,p.z_IP3,self.vm,self.T,p)

                # update concentration of IP3:
                self.cIP3 = self.cIP3 + fIP3_ED*(cells.cell_sa/cells.cell_vol)*p.dt

                self.cIP3_flux_gj = fIP3
                self.cIP3_flux_mem = fIP3_ED[cells.mem_to_cells]

            if p.Ca_dyn == 1 and p.ions_dict['Ca'] == 1:
                # electrodiffusion of ions between cell and endoplasmic reticulum
                fER_ca = electroflux(self.cc_cells[self.iCa],self.cc_er[0],self.Dm_er[0],self.tm,self.z_er[0],
                    self.v_er,self.T,p)

                # update concentration of calcium in cell and ER:
                self.cc_cells[self.iCa] = self.cc_cells[self.iCa] - fER_ca*(cells.cell_sa/cells.cell_vol)*p.dt
                self.cc_er[0] = self.cc_er[0] + fER_ca*(cells.cell_sa/(cells.cell_vol*p.ER_vol))*p.dt

                # Electrodiffusion of charge compensation anion
                fER_m = electroflux(self.cc_cells[self.iM],self.cc_er[1],self.Dm_er[1],self.tm,self.z_er[1],
                    self.v_er,self.T,p)

                # update concentration of anion in cell and ER:
                self.cc_cells[self.iM] = self.cc_cells[self.iM] - fER_m*(cells.cell_sa/cells.cell_vol)*p.dt
                self.cc_er[1] = self.cc_er[1] + fER_m*(cells.cell_sa/(cells.cell_vol*p.ER_vol))*p.dt

                # recalculate the net, unbalanced charge and voltage in each cell:
                q_cells = get_charge(self.cc_cells,self.z_array,cells.cell_vol,p)
                self.vm = get_volt(q_cells,cells.cell_sa,p)

                q_er = get_charge(self.cc_er,self.z_array_er,p.ER_vol*cells.cell_vol,p)
                v_er_o = get_volt(q_er,p.ER_sa*cells.cell_sa,p)
                self.v_er = v_er_o - self.vm

            # if p.voltage_dye=1 electrodiffuse voltage sensitive dye between cell and environment
            if p.voltage_dye ==1:

                fdye_ED = electroflux(self.cDye_env,self.cDye_cell,self.id_cells*p.Dm_Dye,self.tm,p.z_Dye,self.vm,self.T,p)

                # update dye concentration
                self.cDye_cell = self.cDye_cell + fdye_ED*(cells.cell_sa/cells.cell_vol)*p.dt

                # determine flux through gap junctions for voltage dye:
                fDye = electroflux(self.cDye_cell[cells.nn_i][:,0],self.cDye_cell[cells.nn_i][:,1],
                    self.id_gj*p.Do_Dye,self.gjl,p.z_Dye,self.vgj,self.T,p)

                # update cell voltage-sensitive dye concentration due to gap junction flux:
                deltac_dye = self.gjopen*(fDye)*p.dt
                self.cDye_cell = self.cDye_cell + ((cells.cell_sa*p.gj_surface)/cells.cell_vol)*np.dot(deltac_dye, cells.gjMatrix)

                self.Dye_flux_mem = fdye_ED[cells.mem_to_cells]
                self.Dye_flux_gj = fDye

            if p.dynamic_noise == 1 and p.ions_dict['P']==1:
                # add a random walk on protein concentration to generate dynamic noise:
                self.protein_noise_factor = p.dynamic_noise_level*(np.random.random(len(cells.cell_i)) - 0.5)
                self.cc_cells[self.iP] = self.cc_cells[self.iP]*(1+ self.protein_noise_factor)

                # recalculate the net, unbalanced charge and voltage in each cell:
                q_cells = get_charge(self.cc_cells,self.z_array,cells.cell_vol,p)
                self.vm = get_volt(q_cells,cells.cell_sa,p)

            check_v(self.vm)


            if t in tsamples:

                self.get_current(cells,p)   # get the current in the gj network connection of cells

                # add the new concentration and voltage data to the time-storage matrices:
                self.efield_gj_x_time.append(self.E_gj_x[:])

                self.efield_gj_y_time.append(self.E_gj_y[:])

                concs = self.cc_cells[:]
                self.cc_time.append(concs)
                concs = None

                envsc = self.cc_env[:]
                self.cc_env_time.append(envsc)
                envsc = None

                self.I_gj_x_time.append(self.I_gj_x[:])
                self.I_gj_y_time.append(self.I_gj_y[:])

                self.I_mem_time.append(self.I_mem[:])

                self.vm_time.append(self.vm[:])

                self.dvm_time.append(self.dvm[:])

                self.gjopen_time.append(self.gjopen[:])

                self.time.append(t)

                if p.scheduled_options['IP3'] != 0 or p.Ca_dyn == True:

                    self.cIP3_time.append(self.cIP3[:])

                    self.cIP3_flux_gj_time.append(self.cIP3_flux_gj[:])

                    self.cIP3_flux_mem_time.append(self.cIP3_flux_mem[:])

                if p.voltage_dye ==1:

                    self.cDye_time.append(self.cDye_cell[:])

                    self.Dye_flux_gj_time.append(self.Dye_flux_gj[:])

                    ffmemDye = self.Dye_flux_mem[:]
                    self.Dye_flux_mem_time.append(self.Dye_flux_mem[:])

                if p.Ca_dyn == 1 and p.ions_dict['Ca']==1:

                    self.cc_er_time.append(self.cc_er[:])


                if p.plot_while_solving == True:
                    self.checkPlot.updatePlot(self,p)

                        # get time for loop and estimate total time for simulation
            if do_once == True:
                loop_time = time.time() - loop_measure
                if p.run_sim == True:
                    time_estimate = round(loop_time*p.sim_tsteps,2)
                else:
                    time_estimate = round(loop_time*p.init_tsteps,2)
                loggers.log_info("This run should take approximately " + str(time_estimate) + ' s to compute...')
                do_once = False

        # Find embeded functions that can't be pickled...
        for key, valu in vars(self.dyna).items():
            if type(valu) == interp.interp1d:
                setattr(self.dyna,key,None)

        self.checkPlot = None

        if p.run_sim == False:

            # celf = copy.deepcopy(self)

            datadump = [self,cells,p]
            fh.saveSim(self.savedInit,datadump)
            message_1 = 'Initialization run saved to' + ' ' + p.init_path
            loggers.log_info(message_1)

        elif p.run_sim == True:


            datadump = [self,cells,p]
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


        if p.scheduled_options['IP3'] != 0 or p.Ca_dyn == True:
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

        self.tissueInit(cells,p)   # Initialize all structures used for gap junctions, ion channels, and other dynamics

        # Reinitialize all time-data structures
        self.cc_time = []  # data array holding the concentrations at time points
        self.cc_env_time = [] # data array holding extracellular concentrations at time points

        self.vm_time = []  # data array holding voltage at time points
        self.vcell_time = []
        self.venv_time = []

        self.dvm_time = []  # data array holding derivative of voltage at time points
        self.time = []     # time values of the simulation

        self.gjopen_time = []   # stores the fractional gap junction open state at each time
        self.fgj_time = []      # stores the gj fluxes for each ion at each time
        self.Igj_time = []      # current for each gj at each time

        self.I_gj_x_time = []    # initialize gap junction current data storage
        self.I_gj_y_time = []    # initialize gap junction current data storage
        self.I_env_x_time = []   # initialize environmental matrix data storage
        self.I_env_y_time = []   # initialize environmental matrix data storage
        self.I_mem_time = []   # initialize membrane matrix data storage

        self.cc_er_time = []   # retains er concentrations as a function of time
        self.cIP3_time = []    # retains cellular ip3 concentrations as a function of time

        self.efield_gj_x_time = []   # matrices storing smooth electric field in gj connected cells
        self.efield_gj_y_time = []

        self.efield_ecm_x_time = []   # matrices storing smooth electric field in ecm
        self.efield_ecm_y_time = []

        self.vm_Matrix = [] # initialize matrices for resampled data sets (used in smooth plotting and streamlines)
        vm_dato = np.zeros(len(cells.mem_i))       # FIXME this will be HUGE do we need it?>>
        dat_grid_vm = vertData(vm_dato,cells,p)
        self.vm_Matrix.append(dat_grid_vm[:])

        if p.voltage_dye == True:

            self.Dye_flux_ecm_time = []
            self.Dye_flux_env_time = []
            self.Dye_flux_gj_time = []
            self.Dye_flux_mem_time = []
            self.cDye_time = []    # retains voltage-sensitive dye concentration as a function of time

        # gap junction specific arrays:
        self.id_gj = np.ones(len(cells.nn_i))  # identity array for gap junction indices...
        self.gjopen = np.ones(len(cells.nn_i))   # holds gap junction open fraction for each gj
        self.gjl = np.zeros(len(cells.nn_i))    # gj length for each gj
        self.gjl[:] = cells.nn_len

        # get the net, unbalanced charge and corresponding voltage in each cell to initialize values of voltages:
        self.update_V_ecm(cells,p,0)

        self.vm_to = self.vm[:]   # create a copy of the original voltage

         # create a time-steps vector appropriate for simulation type:
        if p.run_sim == True:
            tt = np.linspace(0,p.sim_tsteps*p.dt,p.sim_tsteps)
        else:
            tt = np.linspace(0,p.init_tsteps*p.dt,p.init_tsteps)   # timestep vector

        i = 0 # resample the time vector to save data at specific times:
        tsamples =[]
        resample = p.t_resample
        while i < len(tt)-resample:
            i = i + resample
            tsamples.append(tt[i])
        tsamples = set(tsamples)

        # report
        if p.run_sim == True:

            loggers.log_info('Your simulation (with extracellular spaces) is running from '+ str(0) + ' to '+ str(round(p.sim_tsteps*p.dt,3))
                         + ' seconds of in-world time.')

        else:
            loggers.log_info('Your initialization (with extracellular spaces) is running from '+ str(0) + ' to '+ str(round(p.init_tsteps*p.dt,3))
                         + ' seconds of in-world time.')


        if p.plot_while_solving == True:

            self.checkPlot = viz.PlotWhileSolving(cells,self,p,clrAutoscale = p.autoscale_Vmem, clrMin = p.Vmem_min_clr,
                clrMax = p.Vmem_max_clr)

        do_once = True  # a variable to time the loop only once

        for t in tt:   # run through the loop

            if do_once == True:
                loop_measure = time.time()

            self.fluxes_mem.fill(0)  # reinitialize flux storage device

            self.dvm = (self.vm - self.vm_to)/p.dt    # calculate the change in the voltage derivative
            self.vm_to = self.vm[:]       # reassign the history-saving vm

            if p.Ca_dyn ==1 and p.ions_dict['Ca'] == 1:

                self.dcc_ER = (self.cc_er - self.cc_er_to)/p.dt
                self.cc_er_to = self.cc_er[:]

            # calculate the values of scheduled and dynamic quantities (e.g. ion channel multipliers):
            if p.run_sim == True:
                self.dyna.runAllDynamics(self,cells,p,t)

            # run the Na-K-ATPase pump:
            fNa_NaK, fK_NaK =\
                pumpNaKATP(self.cc_cells[self.iNa][cells.mem_to_cells],self.cc_env[self.iNa][cells.map_mem2ecm],
                    self.cc_cells[self.iK][cells.mem_to_cells],self.cc_env[self.iK][cells.map_mem2ecm],
                    self.vm,self.T,p,self.NaKATP_block)

            # modify the resulting fluxes by the electroosmosis membrane redistribution factor (if calculated)
            fNa_NaK = self.rho_channel*fNa_NaK
            fK_NaK = self.rho_channel*fK_NaK

            self.fluxes_mem[self.iNa] = fNa_NaK
            self.fluxes_mem[self.iK] = fK_NaK

            # update the concentrations
            self.update_C_ecm(self.iNa,fNa_NaK,cells,p)
            self.update_C_ecm(self.iK,fK_NaK,cells,p)

            # recalculate the net, unbalanced charge and voltage in each cell:
            self.update_V_ecm(cells,p,t)

            if p.ions_dict['Ca'] == 1:

                f_CaATP = pumpCaATP(self.cc_cells[self.iCa][cells.mem_to_cells],self.cc_env[self.iCa][cells.map_mem2ecm],
                        self.vm,self.T,p)

                # update calcium concentrations in cell and ecm:
                self.update_C_ecm(self.iCa,f_CaATP,cells,p)

                self.fluxes_mem[self.iCa] = f_CaATP

                # recalculate the net, unbalanced charge and voltage in each cell:
                self.update_V_ecm(cells,p,t)

                if p.Ca_dyn ==1:

                    self.cc_er[0],self.cc_cells[self.iCa], _ =\
                        pumpCaER(self.cc_er[0],self.cc_cells[self.iCa],cells.cell_sa,p.ER_vol*cells.cell_vol,
                            cells.cell_vol,self.v_er,self.T,p)

                    # recalculate the net, unbalanced charge and voltage in each cell:
                    self.update_V_ecm(cells,p,t)

                    q_er = get_charge(self.cc_er,self.z_array_er,p.ER_vol*cells.cell_vol,p)
                    self.v_er = get_volt(q_er,p.ER_sa*cells.cell_sa,p)

            if p.ions_dict['H'] == 1:

                self.Hplus_electrofuse_ecm(cells,p,t)

                if p.HKATPase_dyn == 1:

                    self.Hplus_HKATP_ecm(cells,p,t)

                if p.VATPase_dyn == 1:

                    self.Hplus_VATP_ecm(cells,p,t)

            # electro-diffuse all ions (except for proteins, which don't move) across the cell membrane:

            shuffle(self.movingIons)

            for i in self.movingIons:

                f_ED = electroflux(self.cc_env[i][cells.map_mem2ecm],self.cc_cells[i][cells.mem_to_cells],
                         self.Dm_cells[i], self.tm, self.zs[i], self.vm, self.T, p,
                         rho=self.rho_channel)

                self.fluxes_mem[i] = self.fluxes_mem[i] + f_ED

                # update ion concentrations in cell and ecm:
                self.update_C_ecm(i,f_ED,cells,p)

                # update concentrations in the extracellular spaces:
                self.update_ecm(cells,p,t,i)

                # uppdate concentrations for extracellular spaces in contact with environment:
                # self.update_environment(cells,p,t,i)

                # update flux between cells due to gap junctions
                self.update_gj(cells,p,t,i)

                # # recalculate the net, unbalanced charge and voltage in each cell:
                self.update_V_ecm(cells,p,t)

            self.get_Efield(cells,p)

            if p.scheduled_options['IP3'] != 0 or p.Ca_dyn == True:

                self.update_IP3(cells,p,t)

            if p.Ca_dyn == 1 and p.ions_dict['Ca'] == 1:

               self.update_er(cells,p,t)

            # if p.voltage_dye=1 electrodiffuse voltage sensitive dye between cell and environment
            if p.voltage_dye ==1:

                self.update_dye(cells,p,t)

            if p.dynamic_noise == 1 and p.ions_dict['P']==1:
                # add a random walk on protein concentration to generate dynamic noise:
                self.protein_noise_factor = p.dynamic_noise_level*(np.random.random(len(cells.cell_i)) - 0.5)
                self.cc_cells[self.iP] = self.cc_cells[self.iP]*(1+ self.protein_noise_factor)

                # recalculate the net, unbalanced charge and voltage in each cell:
                self.update_V_ecm(cells,p,t)

            check_v(self.vm)

            # if desired, electroosmosis of membrane channels
            if p.sim_eosmosis == True:

                self.eosmosis(cells,p)

            if t in tsamples:

                # #
                self.get_current(cells,p)   # get the current in the gj network connection of cells

                # add the new concentration and voltage data to the time-storage matrices:

                self.efield_gj_x_time.append(self.E_gj_x[:])

                self.efield_gj_y_time.append(self.E_gj_y[:])

                self.efield_ecm_x_time.append(self.E_env_x[:])

                self.efield_ecm_y_time.append(self.E_env_y[:])

                concs = self.cc_cells[:]
                self.cc_time.append(concs)
                concs = None

                ecmsc = self.cc_env[:]
                self.cc_env_time.append(ecmsc)
                ecmsc = None

                self.vm_time.append(self.vm[:])

                self.dvm_time.append(self.dvm[:])

                self.gjopen_time.append(self.gjopen[:])

                self.vcell_time.append(self.v_cell[:])

                self.venv_time.append(self.v_env[:])

                self.I_gj_x_time.append(self.I_gj_x[:])
                self.I_gj_y_time.append(self.I_gj_y[:])

                self.I_env_x_time.append(self.I_env_x[:])
                self.I_env_y_time.append(self.I_env_y[:])

                self.I_mem_time.append(self.I_mem[:])

                # calculate interpolated verts and midpoint data for Vmem:
                dat_grid_vm = vertData(self.vm[:],cells,p)    # FIXME this will be huge -- reduce it!

                self.vm_Matrix.append(dat_grid_vm[:])

                self.time.append(t)

                if p.scheduled_options['IP3'] != 0 or p.Ca_dyn == True:
                    ccIP3 = self.cIP3[:]
                    self.cIP3_time.append(ccIP3)
                    ccIP3 = None

                if p.voltage_dye ==1:
                    ccDye_cells = self.cDye_cell[:]
                    self.cDye_time.append(ccDye_cells)
                    ccDye_cells = None

                if p.Ca_dyn == 1 and p.ions_dict['Ca']==1:
                    ccer = self.cc_er[:]
                    self.cc_er_time.append(ccer)
                    ccer = None

                if p.plot_while_solving == True:
                    self.checkPlot.updatePlot(self,p)

                        # get time for loop and estimate total time for simulation
            if do_once == True:
                loop_time = time.time() - loop_measure

                if p.run_sim == True:
                    time_estimate = round(loop_time*p.sim_tsteps,2)
                else:
                    time_estimate = round(loop_time*p.init_tsteps,2)
                loggers.log_info("This run should take approximately " + str(time_estimate) + ' s to compute...')
                do_once = False


        # self.dyna.scalar_Namem = None

         # Find embeded functions that can't be pickled...
        for key, valu in vars(self.dyna).items():
            if type(valu) == interp.interp1d:
                setattr(self.dyna,key,None)

        self.checkPlot = None

        if p.run_sim == False:

            # celf = copy.deepcopy(self)

            datadump = [self,cells,p]
            fh.saveSim(self.savedInit,datadump)
            message_1 = 'Initialization run saved to' + ' ' + p.init_path
            loggers.log_info(message_1)

        elif p.run_sim == True:
            # celf = copy.deepcopy(self)
            datadump = [self,cells,p]
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

            final_pH_ecm = -np.log10(np.mean((self.cc_env_time[-1][self.iH])/1000))
            loggers.log_info('Final extracellular pH '+ str(np.round(final_pH_ecm,2)))


        if p.scheduled_options['IP3'] != 0 or p.Ca_dyn == True:

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

        print(self.D_env[self.iNa][0])

        plt.close()
        loggers.log_info('Simulation completed successfully.')

    def update_V_ecm(self,cells,p,t):

        self.rho_cells = get_charge_density(self.cc_cells, self.z_array, p)
        self.rho_env = get_charge_density(self.cc_env, self.z_array_env, p)
        self.v_cell = get_Vcell(self,cells,p)
        self.v_env = get_Venv(self,cells,p)

        self.vm = self.v_cell[cells.mem_to_cells] - self.v_env[cells.map_mem2ecm]  # calculate v_mem

    def update_C_ecm(self,ion_i,flux,cells,p):

        c_cells = copy.copy(self.cc_cells[ion_i][:])
        c_env = copy.copy(self.cc_env[ion_i][:])

        # delta_c_cells = flux*(cells.mem_sa/cells.cell_vol[cells.mem_to_cells])*p.dt
        # delta_c_env = flux*(cells.mem_sa/cells.ecm_vol)*p.dt
        #
        # self.cc_cells[ion_i] = c_cells + np.dot(delta_c_cells, cells.cell_UpdateMatrix)
        #
        # self.cc_env[ion_i] = c_env - np.dot(delta_c_env, cells.ecm_UpdateMatrix)

        d_c_cells = flux*(cells.mem_sa/cells.cell_vol[cells.mem_to_cells])
        d_c_env = flux*(cells.mem_sa/cells.ecm_vol)

        delta_cells =  np.dot(d_c_cells, cells.cell_UpdateMatrix)
        delta_env = np.dot(d_c_env, cells.ecm_UpdateMatrix)

        self.cc_cells[ion_i] = rk4(c_cells,delta_cells,p)

        self.cc_env[ion_i] = rk4(c_env,-delta_env,p)


        # self.cc_env[ion_i][cells.map_mem2ecm] = c_env[cells.map_mem2ecm] - delta_c_env

    def Hplus_electrofuse_ecm(self,cells,p,t):

        # electrofuse the H+ ion between the cytoplasm and the ecms
        _,_,f_H1 = \
            electrofuse(self.cc_ecm[self.iH][cells.mem_to_ecm],self.cc_cells[self.iH][cells.mem_to_cells],
                self.Dm_cells[self.iH],self.tm[cells.mem_to_cells],cells.mem_sa,
                cells.ecm_vol[cells.mem_to_ecm],cells.cell_vol[cells.mem_to_cells],self.zs[self.iH],
                self.vm,self.T,p)

        self.fluxes_mem[self.iH] =  self.fluxes_mem[self.iH] + f_H1

        H_cell_to = self.cc_cells[self.iH][:]     # keep track of original H+ in cell and ecm
        H_ecm_to = self.cc_ecm[self.iH][:]

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
        self.update_V_ecm(cells,p,t)

    def Hplus_HKATP_ecm(self,cells,p,t):

        # if HKATPase pump is desired, run the H-K-ATPase pump:
        _,_,_,_, f_H2, f_K2 =\
        pumpHKATP(self.cc_cells[self.iH][cells.mem_to_cells],self.cc_ecm[self.iH][cells.mem_to_ecm],
            self.cc_cells[self.iK][cells.mem_to_cells],self.cc_ecm[self.iK][cells.mem_to_ecm],
            cells.mem_sa,cells.cell_vol[cells.mem_to_cells],cells.ecm_vol[cells.mem_to_ecm],
            self.vm,self.T,p,self.HKATP_block)

        self.fluxes_mem[self.iH] =  self.fluxes_mem[self.iH] + f_H2
        self.fluxes_mem[self.iK] =  self.fluxes_mem[self.iK] + f_K2

        H_cell_to = self.cc_cells[self.iH][:]     # keep track of original H+ in cell and ecm
        H_ecm_to = self.cc_ecm[self.iH][:]

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
        self.update_V_ecm(cells,p,t)

    def Hplus_VATP_ecm(self,cells,p,t):

        # if HKATPase pump is desired, run the H-K-ATPase pump:
        _, _, f_H3 =\
        pumpVATP(self.cc_cells[self.iH][cells.mem_to_cells],self.cc_ecm[self.iH][cells.mem_to_ecm],
            cells.mem_sa,cells.cell_vol[cells.mem_to_cells],cells.ecm_vol[cells.mem_to_ecm],
            self.vm,self.T,p)

        self.fluxes_mem[self.iH] =  self.fluxes_mem[self.iH] + f_H3

        H_cell_to = self.cc_cells[self.iH][:]     # keep track of original H+ in cell and ecm
        H_ecm_to = self.cc_ecm[self.iH][:]

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
        self.update_V_ecm(cells,p,t)

    def update_gj(self,cells,p,t,i):


        # calculate voltage difference (gradient*len_gj) between gj-connected cells:
        if p.sim_ECM == True:

            self.vgj = self.v_cell[cells.nn_i][:,1]- self.v_cell[cells.nn_i][:,0]

        else:
            self.vgj = self.vm[cells.nn_i][:,1]- self.vm[cells.nn_i][:,0]


        if p.v_sensitive_gj == True:
            # determine the open state of gap junctions:
            self.gjopen = self.gj_block*((1.0 - tb.step(abs(self.vgj),p.gj_vthresh,p.gj_vgrad) + 0.1))

        # voltage gradient:
        grad_vgj = self.vgj/cells.nn_len

        grad_vgj_x = grad_vgj*cells.nn_vects[:,2]
        grad_vgj_y = grad_vgj*cells.nn_vects[:,3]

        # concentration gradient for ion i:
        grad_cgj = (self.cc_cells[i][cells.nn_i][:,1] - self.cc_cells[i][cells.nn_i][:,0])/cells.nn_len

        grad_cgj_x = grad_cgj*cells.nn_vects[:,2]
        grad_cgj_y = grad_cgj*cells.nn_vects[:,3]

        # midpoint concentration:
        c = (self.cc_cells[i][cells.nn_i][:,1] + self.cc_cells[i][cells.nn_i][:,0])/2

        fgj_x,fgj_y = nernst_planck_flux(c,grad_cgj_x,grad_cgj_y,grad_vgj_x,grad_vgj_y,
            self.D_gj[i]*self.gjopen*p.gj_surface,self.zs[i],self.T,p)

        fgj = fgj_x*cells.nn_vects[:,2] + fgj_y*cells.nn_vects[:,3]

        delta_cc = np.dot(cells.gjMatrix,fgj)

        self.cc_cells[i] = self.cc_cells[i] + p.dt*delta_cc

        self.fluxes_gj_x[i] = fgj_x  # store gap junction flux for this ion
        self.fluxes_gj_y[i] = fgj_y  # store gap junction flux for this ion

         # recalculate the net, unbalanced charge and voltage in each cell:
        if p.sim_ECM == True:
            self.update_V_ecm(cells,p,t)


# old way of doing this-----------------------------------------------------------------------------------------
         # calculate voltage difference between cells:
        # if p.sim_ECM == True:
        #     vmA,vmB = self.v_cell[cells.gap_jun_i][:,0], self.v_cell[cells.gap_jun_i][:,1]
        #
        # else:
        #     vmA,vmB = self.vm[cells.gap_jun_i][:,0], self.vm[cells.gap_jun_i][:,1]
        #
        # self.vgj = vmB - vmA
        #
        # if p.v_sensitive_gj == True:
        #     # determine the open state of gap junctions:
        #     self.gjopen = self.gj_block*((1.0 - tb.step(abs(self.vgj),p.gj_vthresh,p.gj_vgrad) + 0.1))
        #
        # # determine flux through gap junctions for this ion:
        #
        # fgj = electroflux(self.cc_cells[i][cells.gap_jun_i][:,0],self.cc_cells[i][cells.gap_jun_i][:,1],
        #     self.D_gj[i],self.gjl,self.zs[i],self.vgj,self.T,p)
        #
        # # deltac_gj = self.gjopen*(fgj)*p.dt
        #
        # # update cell concentration due to gap junction flux:
        # # self.cc_cells[i] = self.cc_cells[i] + ((cells.cell_sa*p.gj_surface)/cells.cell_vol)*np.dot(deltac_gj, cells.gjMatrix)
        #
        # d_gj = self.gjopen*(fgj)
        #
        # delta_gj = ((cells.cell_sa*p.gj_surface)/cells.cell_vol)*np.dot(d_gj, cells.gjMatrix)
        #
        # self.cc_cells[i] = rk4(self.cc_cells[i],delta_gj,p)
        #
        # # recalculate the net, unbalanced charge and voltage in each cell:
        # # self.update_V_ecm(cells,p,t)
        #
        # self.fluxes_gj[i] = fgj  # store gap junction flux for this ion

    def update_ecm(self,cells,p,t,i):  # FIXME redo this whole section!!!

        #------------------------------
         # make v_env and cc_env into 2d matrices
        cenv = self.cc_env[i][:]
        denv = self.D_env[i][:]

        self.v_env = self.v_env.reshape(cells.X.shape)
        # self.v_env = ndimage.filters.gaussian_filter(self.v_env, 1.0)
        # self.rho_env = self.rho_env.reshape(cells.X.shape)

        cenv = cenv.reshape(cells.X.shape)
        # cenv = ndimage.filters.gaussian_filter(cenv, 1.0)

        denv = denv.reshape(cells.X.shape)

        flux_mult_x = cells.env_flux_mult_x.reshape(cells.X.shape)
        flux_mult_y = cells.env_flux_mult_y.reshape(cells.X.shape)

        # calculate gradients in the environment
        grad_V_env_x, grad_V_env_y = gradient(self.v_env, cells.delta)
        grad_cc_env_x, grad_cc_env_y = gradient(cenv, cells.delta)

        # calculate fluxes for electrodiffusive transport:
        f_env_x, f_env_y = nernst_planck_flux(cenv,grad_cc_env_x,grad_cc_env_y,
            grad_V_env_x, grad_V_env_y, denv,self.zs[i],self.T,p)

        # modify fluxes to create anisotropic diffusion within the cell cluster (along cell membranes):
        f_env_x = f_env_x*flux_mult_x
        f_env_y = f_env_y*flux_mult_y

        if p.closed_bound == True:

            # if boundaries are closed, set normal fluxes at boundary to zero:
            f_env_x[:,0] = 0
            f_env_x[:,-1] = 0

            f_env_y[0,:] = 0
            f_env_y[-1,:]=0

        # calculate the spatial derivatives of the negative of each flux component to simulate divergence operator:
        f_xx, _ = gradient(f_env_x, cells.delta)
        _, f_yy = gradient(f_env_y, cells.delta)

        # change in time is:
        delta_c = f_xx + f_yy

        # cenv = cenv + delta_c*p.dt
        cenv = rk4(cenv, delta_c,p)

        # if p.closed_bound == True:
        #     # Neumann boundary condition (flux at boundary)
        #     # zero flux boundaries for concentration:
        #     cenv[:,-1] = cenv[:,-2]
        #     cenv[:,0] = cenv[:,1]
        #     cenv[0,:] = cenv[1,:]
        #     cenv[-1,:] = cenv[-2,:]

        if p.closed_bound == False:
            # if the boundary is open, set the concentration at the boundary
            # open boundary
            self.cc_env[i][cells.bL_k] = self.c_env_bound[i]
            self.cc_env[i][cells.bR_k] = self.c_env_bound[i]
            self.cc_env[i][cells.bTop_k] = self.c_env_bound[i]
            self.cc_env[i][cells.bBot_k] = self.c_env_bound[i]

        # reshape the matrices into vectors
        # self.rho_env = self.rho_env.ravel()
        self.v_env = self.v_env.ravel()
        self.cc_env[i] = cenv.ravel()

         # recalculate the net, unbalanced charge and voltage in each ecm space:
        # self.rho_env = get_charge_density(self.cc_env,self.z_array_env,p)
        # self.v_env = get_Venv(self,cells,p)

        self.fluxes_env_x[i] = -f_env_x.ravel()  # store ecm junction flux for this ion
        self.fluxes_env_y[i] = -f_env_y.ravel()  # store ecm junction flux for this ion

        #-----------------------------
        #  # calculate voltage differences between ecm <---> ecm nearest neighbours:
        # self.v_ec2ec = self.v_ecm[cells.ecm_nn_i[:,1]] - self.v_ecm[cells.ecm_nn_i[:,0]]
        #
        # # electrodiffuse through ecm <---> ecm junctions
        #
        # # _,_,f_ecm = electrofuse(self.cc_ecm[i][cells.ecm_nn_i[:,0]],self.cc_ecm[i][cells.ecm_nn_i[:,1]],
        # #         self.D_ecm_juncs[i],cells.len_ecm_junc,self.ec2ec_sa,
        # #         cells.ecm_vol[cells.ecm_nn_i[:,0]],cells.ecm_vol[cells.ecm_nn_i[:,1]],
        # #         self.zs[i],self.v_ec2ec,self.T,p)
        #
        # f_ecm = electroflux(self.cc_ecm[i][cells.ecm_nn_i[:,0]],self.cc_ecm[i][cells.ecm_nn_i[:,1]],
        #         self.D_ecm_juncs[i],cells.len_ecm_junc,self.ec2ec_sa,self.zs[i],self.v_ec2ec,self.T,p)
        #
        # self.cc_ecm[i] = (self.cc_ecm[i]*cells.ecm_vol + np.dot(f_ecm*p.dt,cells.ecmMatrix))/cells.ecm_vol
        #
        #  # recalculate the net, unbalanced charge and voltage in each ecm space:
        # self.rho_env = get_charge_density(self.cc_ecm,self.z_array_ecm,p)
        # self.v_ecm = get_Vecm(self.rho_env,p)
        #
        # self.fluxes_ecm[i] = f_ecm  # store ecm junction flux for this ion

    def update_er(self,cells,p,t):

         # electrodiffusion of ions between cell and endoplasmic reticulum
        self.cc_cells[self.iCa],self.cc_er[0],_ = \
        electrofuse(self.cc_cells[self.iCa],self.cc_er[0],self.Dm_er[0],self.tm,p.ER_sa*cells.cell_sa,
            cells.cell_vol,p.ER_vol*cells.cell_vol,self.z_er[0],self.v_er,self.T,p)

        # Electrodiffusion of charge compensation anion
        self.cc_cells[self.iM],self.cc_er[1],_ = \
        electrofuse(self.cc_cells[self.iM],self.cc_er[1],self.Dm_er[1],self.tm,p.ER_sa*cells.cell_sa,
            cells.cell_vol,p.ER_vol*cells.cell_vol,self.z_er[1],self.v_er,self.T,p)

        # recalculate the net, unbalanced charge and voltage in each cell:
        self.update_V_ecm(cells,p,t)

        q_er = get_charge(self.cc_er,self.z_array_er,p.ER_vol*cells.cell_vol,p)
        self.v_er = get_volt(q_er,p.ER_sa*cells.cell_sa,p)

    def update_dye(self,cells,p,t):

        # _,_,flux_dye = \
        #                 electrofuse(self.cDye_ecm[cells.mem_to_ecm],self.cDye_cell[cells.mem_to_cells],
        #                     p.Dm_Dye*self.id_cells[cells.mem_to_cells],self.tm[cells.mem_to_cells],cells.mem_sa,
        #                     cells.ecm_vol[cells.mem_to_ecm],cells.cell_vol[cells.mem_to_cells],p.z_Dye,self.vm,self.T,p)

        flux_dye = electroflux(self.cDye_ecm[cells.mem_to_ecm],self.cDye_cell[cells.mem_to_cells],
                            p.Dm_Dye,self.tm[cells.mem_to_cells],cells.mem_sa,
                            p.z_Dye,self.vm,self.T,p)

         # update the dye concentrations in the cell and ecm due to ED fluxes at membrane
        self.cDye_cell = self.cDye_cell + \
                            np.dot((flux_dye/cells.cell_vol[cells.mem_to_cells])*p.dt,cells.cell_UpdateMatrix)

        self.cDye_ecm = self.cDye_ecm - \
                            np.dot((flux_dye/cells.ecm_vol[cells.mem_to_ecm])*p.dt,cells.ecm_UpdateMatrix)

        # determine flux through gap junctions for voltage dye:
        # _,_,fDye_gj = electrofuse(self.cDye_cell[cells.gap_jun_i][:,0],self.cDye_cell[cells.gap_jun_i][:,1],
        #     self.id_gj*p.Do_Dye,self.gjl,self.gjopen*self.gjsa,cells.cell_vol[cells.gap_jun_i][:,0],
        #     cells.cell_vol[cells.gap_jun_i][:,1],p.z_Dye,self.vgj,self.T,p)
        fDye_gj = electroflux(self.cDye_cell[cells.nn_i][:,0],self.cDye_cell[cells.nn_i][:,1],
            p.Do_Dye,self.gjl,self.gjopen,p.z_Dye,self.vgj,self.T,p)

        # update cell voltage-sensitive dye concentration due to gap junction flux:
        self.cDye_cell = (self.cDye_cell*cells.cell_vol + np.dot((fDye_gj*p.dt), cells.gjMatrix))/cells.cell_vol

        # electrodiffuse dye through ecm <---> ecm junctions
        # _,_,flux_ecm_dye = electrofuse(self.cDye_ecm[cells.ecm_nn_i[:,0]],self.cDye_ecm[cells.ecm_nn_i[:,1]],
        #         self.id_ecm*p.Do_Dye,cells.len_ecm_junc,self.ec2ec_sa,
        #         cells.ecm_vol[cells.ecm_nn_i[:,0]],cells.ecm_vol[cells.ecm_nn_i[:,1]],
        #         p.z_Dye,self.v_ec2ec,self.T,p)
        flux_ecm_dye = electroflux(self.cDye_ecm[cells.ecm_nn_i[:,0]],self.cDye_ecm[cells.ecm_nn_i[:,1]],
                p.Do_Dye,cells.len_ecm_junc,self.ec2ec_sa,p.z_Dye,self.v_ec2ec,self.T,p)
                    #
        self.cDye_ecm = (self.cDye_ecm*cells.ecm_vol + np.dot(flux_ecm_dye*p.dt,cells.ecmMatrix))/cells.ecm_vol

        # electrodiffuse dye between environmental and ecm junctions:

        self.cDye_ecm[cells.bflags_ecm],self.cDye_env,fDye_env = electroflux(self.cDye_ecm[cells.bflags_ecm],
            self.cDye_env, self.id_env*p.Do_Dye,cells.len_ecm_junc[cells.bflags_ecm],
            self.ec2ec_sa[cells.bflags_ecm], cells.ecm_vol[cells.bflags_ecm],self.env_vol,
            p.z_Dye, self.v_ec2env, self.T, p,ignoreECM=True)

    def update_IP3(self,cells,p,t):
         # determine flux through gap junctions for IP3:
        # _,_,fIP3 = electrofuse(self.cIP3[cells.gap_jun_i][:,0],self.cIP3[cells.gap_jun_i][:,1],
        #     self.id_gj*p.Do_IP3,self.gjl,self.gjopen*self.gjsa,cells.cell_vol[cells.gap_jun_i][:,0],
        #     cells.cell_vol[cells.gap_jun_i][:,1],p.z_IP3,self.vgj,self.T,p)
        fIP3 = electroflux(self.cIP3[cells.nn_i][:,0],self.cIP3[cells.nn_i][:,1], p.Do_IP3,self.gjl,
            self.gjopen,p.z_IP3,self.vgj,self.T,p)

        # update cell IP3 concentration due to gap junction flux:
        self.cIP3 = (self.cIP3*cells.cell_vol + np.dot((fIP3*p.dt), cells.gjMatrix))/cells.cell_vol

        # electrodiffuse IP3 between cell and environment:
        # _,_,flux_IP3 = \
        #             electrofuse(self.cIP3_ecm[cells.mem_to_ecm],self.cIP3[cells.mem_to_cells],
        #                 p.Dm_IP3*self.id_cells[cells.mem_to_cells],self.tm[cells.mem_to_cells],cells.mem_sa,
        #                 cells.ecm_vol[cells.mem_to_ecm],cells.cell_vol[cells.mem_to_cells],p.z_IP3,self.vm,self.T,p)

        flux_IP3 = electroflux(self.cIP3_ecm[cells.mem_to_ecm],self.cIP3[cells.mem_to_cells], p.Dm_IP3,
                        self.tm[cells.mem_to_cells],cells.mem_sa, p.z_IP3,self.vm,self.T,p)

        # update the IP3 concentrations in the cell and ecm due to ED fluxes at membrane
        self.cIP3 = self.cIP3 + \
                            np.dot((flux_IP3/cells.cell_vol[cells.mem_to_cells])*p.dt,cells.cell_UpdateMatrix)

        self.cIP3_ecm = self.cIP3_ecm - \
                            np.dot((flux_IP3/cells.ecm_vol[cells.mem_to_ecm])*p.dt,cells.ecm_UpdateMatrix)

         # electrodiffuse IP3 through ecm <---> ecm junctions
        # _,_,flux_ecm_IP3 = electrofuse(self.cIP3_ecm[cells.ecm_nn_i[:,0]],self.cIP3_ecm[cells.ecm_nn_i[:,1]],
        #         self.id_ecm*p.Do_IP3,cells.len_ecm_junc,self.ec2ec_sa,
        #         cells.ecm_vol[cells.ecm_nn_i[:,0]],cells.ecm_vol[cells.ecm_nn_i[:,1]],
        #         p.z_IP3,self.v_ec2ec,self.T,p)
        flux_ecm_IP3 = electroflux(self.cIP3_ecm[cells.ecm_nn_i[:,0]],self.cIP3_ecm[cells.ecm_nn_i[:,1]],
                p.Do_IP3,cells.len_ecm_junc,self.ec2ec_sa,p.z_IP3,self.v_ec2ec,self.T,p)

        self.cIP3_ecm = (self.cIP3_ecm*cells.ecm_vol + np.dot(flux_ecm_IP3*p.dt,cells.ecmMatrix))/cells.ecm_vol

         # electrodiffuse IP3 between environmental and ecm junctions:
        self.cIP3_ecm[cells.bflags_ecm],self.cIP3_env,fIP3_env = electrofuse(self.cIP3_ecm[cells.bflags_ecm],
            self.cIP3_env, self.id_env*p.Do_IP3,cells.len_ecm_junc[cells.bflags_ecm],
            self.ec2ec_sa[cells.bflags_ecm], cells.ecm_vol[cells.bflags_ecm],self.env_vol,
            p.z_IP3, self.v_ec2env, self.T, p,ignoreECM=True)

    def get_Efield(self,cells,p):

         # calculate voltage difference (gradient*len_gj) between gj-connected cells:
        if p.sim_ECM == True:

            Egj = - (self.v_cell[cells.nn_i][:,1]- self.v_cell[cells.nn_i][:,0])/cells.nn_len

        else:
            Egj = - (self.vm[cells.nn_i][:,1]- self.vm[cells.nn_i][:,0])/cells.nn_len

        # get x and y components of the electric field:
        self.E_gj_x = cells.nn_vects[:,2]*Egj
        self.E_gj_y = cells.nn_vects[:,3]*Egj


        if p.sim_ECM == True:
             # in the environment:
            venv = self.v_env.reshape(cells.X.shape)
            genv_x, genv_y = gradient(venv, cells.delta)

            self.E_env_x = genv_x
            self.E_env_y = genv_y


    def get_current(self,cells,p):

        # calculate current across gap junctions in x direction:
        I_gj_x = np.zeros(len(cells.nn_i))

        for flux_array, zi in zip(self.fluxes_gj_x,self.zs):

            I_i_x = flux_array*zi*p.F

            I_gj_x = I_gj_x + I_i_x

        I_gj_y = np.zeros(len(cells.nn_i))

        for flux_array, zi in zip(self.fluxes_gj_y,self.zs):

            I_i_y = flux_array*zi*p.F

            I_gj_y = I_gj_y + I_i_y

        # interpolate the gj current components to the grid:
        self.I_gj_x = interp.griddata((cells.nn_vects[:,0],cells.nn_vects[:,1]),I_gj_x,(cells.X,cells.Y), fill_value=0)

        self.I_gj_y = interp.griddata((cells.nn_vects[:,0],cells.nn_vects[:,1]),I_gj_y,(cells.X,cells.Y),fill_value=0)

        # calculate current across cell membranes:

        self.I_mem = np.zeros(len(cells.mem_i))
        for flux_array, zi in zip(self.fluxes_mem,self.zs):

            # I_i = (flux_array*zi*p.F)/(self.gjopen*self.gjsa)
            I_i = flux_array*zi*p.F

            self.I_mem = self.I_mem + I_i

        if p.sim_ECM == True:

            self.I_env_x = np.zeros(len(cells.xypts))
            self.I_env_y = np.zeros(len(cells.xypts))

            for flux_array, zi in zip(self.fluxes_env_x,self.zs):

                I_i = flux_array*zi*p.F

                self.I_env_x = self.I_env_x + I_i

            for flux_array, zi in zip(self.fluxes_env_y,self.zs):

                I_i = flux_array*zi*p.F

                self.I_env_y = self.I_env_y + I_i

            self.I_env_x = self.I_env_x.reshape(cells.X.shape)
            self.I_env_y = self.I_env_y.reshape(cells.X.shape)

    def eosmosis(self,cells,p):

        # create interpolation functions for the components of the extracellular electric field:
        E_ECM_x_interp = interp.RectBivariateSpline(cells.X_ecm[0,:], cells.Y_ecm[:,0], self.E_ECM_x)
        E_ECM_y_interp = interp.RectBivariateSpline(cells.X_ecm[0,:], cells.Y_ecm[:,0], self.E_ECM_y)

        # interpret the value of extracellular electric field components at the centre of each membrane:
        Efield_x = E_ECM_x_interp.ev(cells.mem_vects_flat[:,0],cells.mem_vects_flat[:,1])
        Efield_y = E_ECM_y_interp.ev(cells.mem_vects_flat[:,0],cells.mem_vects_flat[:,1])

        # Take the dot product of the membrane tangent vector and the ecm E-field vector (E_tang_mem):
        E_tang_x = cells.mem_vects_flat[:,4]*Efield_x
        E_tang_y = cells.mem_vects_flat[:,5]*Efield_y

        self.E_tang_mem = E_tang_x + E_tang_y

        # map the value of rho channel to the vertices
        self.rho_channel_verts = np.dot(self.rho_channel,cells.matrixMap2Verts)

        # calculate the concentration gradient of the rho factor in the membrane:
        self.grad_c_mem = (self.rho_channel_verts[cells.mem_seg_i[:,1]] -
                self.rho_channel_verts[cells.mem_seg_i[:,0]])/cells.mem_length

        # Calculate the flux due to electroosmosis
        self.f_eosmo = p.mem_const*self.E_tang_mem*self.rho_channel*(cells.mem_length/cells.mem_sa)

        # calculate the flux due to standard diffusion:
        self.f_ficks = p.D_membrane*self.grad_c_mem*(cells.mem_length/cells.mem_sa)

        # calculate the total flux
        f_net = self.f_eosmo + self.f_ficks
        # f_net = self.f_eosmo

        # calculate the change in verts data due to electroosmosis and diffusion:
        self.rho_channel_verts = self.rho_channel_verts - (np.dot(f_net*p.dt, cells.memMatrix))
        # map data back to midpoints
        self.rho_channel = np.dot(self.rho_channel_verts,cells.matrixMap2Mids)

        # make sure nothing is non-zero:
        fix_inds = (self.rho_channel < 0).nonzero()
        self.rho_channel[fix_inds] = 0

def electroflux(cA,cB,Dc,d,zc,vBA,T,p,rho=1):

    """
    Electro-diffusion between two connected volumes. Note for cell work, 'b' is 'inside', 'a' is outside, with
    a positive flux moving from a to b. The voltage is defined as
    Vb - Va (Vba), which is equivalent to Vmem.

    This function defaults to regular diffusion if Vba == 0.0

    This function takes numpy matrix values as input. All inputs must be matrices of
    the same shape.

    This is the Goldman Flux/Current Equation (not to be confused with the Goldman Equation).

    Parameters
    ----------
    cA          concentration in region A [mol/m3] (out)
    cB          concentration in region B [mol/m3] (in)
    Dc          Diffusion constant of c  [m2/s]
    d           Distance between region A and region B [m]
    zc          valence of ionic species c
    vBA         voltage difference between region B (in) and A (out) = Vmem
    p           an instance of the Parameters class


    Returns
    --------
    flux        Chemical flux magnitude between region A and B [mol/s]

    """

    # modify the diffusion constant by the membrane density
    Dc = rho*Dc

    alpha = (zc*vBA*p.F)/(p.R*T)

    #volab = (vola + volb)/2
    #qualityfactor = abs((Dc/d)*(sa/volab)*p.dt*alpha)   # quality factor should be <1.0 for stable simulations

    deno = 1 - np.exp(-alpha)   # calculate the denominator for the electrodiffusion equation,..

    izero = (deno==0).nonzero()     # get the indices of the zero and non-zero elements of the denominator
    inotzero = (deno!=0).nonzero()

    # initialize data matrices to the same shape as input data
    flux = np.zeros(deno.shape)

    if len(deno[izero]):   # if there's anything in the izero array:
         # calculate the flux for those elements as standard diffusion [mol/m2s]:
        flux[izero] = -(Dc[izero]/d[izero])*(cB[izero] - cA[izero])

    if len(deno[inotzero]):   # if there's any indices in the inotzero array:

        # calculate the flux for those elements:
        flux[inotzero] = -((Dc[inotzero]*alpha[inotzero])/d[inotzero])*((cB[inotzero] -
                        cA[inotzero]*np.exp(-alpha[inotzero]))/deno[inotzero])

    return flux

def pumpNaKATP(cNai,cNao,cKi,cKo,Vm,T,p,block):

    """
    Parameters
    ----------
    cNai            Concentration of Na+ inside the cell
    cNao            Concentration of Na+ outside the cell
    cKi             Concentration of K+ inside the cell
    cKo             Concentration of K+ outside the cell
    Vm              Voltage across cell membrane [V]
    p               An instance of Parameters object


    Returns
    -------
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

        alpha = block*p.alpha_NaK*tb.step(delG,p.halfmax_NaK,p.slope_NaK)

        f_Na  = -alpha*cNai*cKo      #flux as [mol/m2s]   scaled to concentrations Na in and K out

    elif p.backward_pumps == True:

        alpha = signG*block*p.alpha_NaK*tb.step(delG,p.halfmax_NaK,p.slope_NaK)

        truth_forwards = signG == 1    # boolean array tagging forward-running pump cells
        truth_backwards = signG == -1  # boolean array tagging backwards-running pump cells

        inds_forwards = (truth_forwards).nonzero()  # indices of forward-running cells
        inds_backwards = (truth_backwards).nonzero() # indices of backward-running cells

        f_Na = np.zeros(len(cNai))

        f_Na[inds_forwards]  = -alpha*cNai*cKo      #flux as [mol/s]   scaled to concentrations Na in and K out

        f_Na[inds_backwards]  = -alpha*cNao*cKi      #flux as [mol/s]   scaled to concentrations K in and Na out

    f_K = -(2/3)*f_Na          # flux as [mol/s]


    return f_Na, f_K

def pumpCaATP(cCai,cCao,Vm,T,p):

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

        alpha = p.alpha_Ca*tb.step(delG,p.halfmax_Ca,p.slope_Ca)

        f_Ca  = -alpha*(cCai)      #flux as [mol/s], scaled to concentration in cell

    elif p.backward_pumps == True:

        alpha = signG*p.alpha_Ca*tb.step(delG,p.halfmax_Ca,p.slope_Ca)

        truth_forwards = signG == 1
        truth_backwards = signG == -1

        inds_forwards = (truth_forwards).nonzero()  # indices of forward-running cells
        inds_backwards = (truth_backwards).nonzero() # indices of backward-running cells

        f_Ca = np.zeros(len(cCai))

        f_Ca[inds_forwards]  = -alpha*(cCai)      #flux as [mol/s], scaled to concentration in cell
        f_Ca[inds_backwards]  = -alpha*(cCao)      #flux as [mol/s], scaled to concentration out of cell

    return f_Ca

def pumpCaER(cCai,cCao,Vm,T,p):

    alpha = p.alpha_CaER

    f_Ca  = alpha*(cCao)*(1.0 - cCai)      #flux as [mol/s]

    return f_Ca

def pumpHKATP(cHi,cHo,cKi,cKo,Vm,T,p,block):

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

        alpha = block*p.alpha_HK*tb.step(delG,p.halfmax_HK,p.slope_HK)
        f_H  = -alpha*cHi*cKo      #flux as [mol/s], scaled by concentrations in and out

    elif p.backward_pumps == True:

        alpha = signG*block*p.alpha_HK*tb.step(delG,p.halfmax_HK,p.slope_HK)

        truth_forwards = signG == 1
        truth_backwards = signG == -1

        inds_forwards = (truth_forwards).nonzero()  # indices of forward-running cells
        inds_backwards = (truth_backwards).nonzero() # indices of backward-running cells

        f_H = np.zeros(len(cHi))

        f_H[inds_forwards]  = -alpha*cHi*cKo      #flux as [mol/s], scaled by concentrations in and out
        f_H[inds_backwards]  = -alpha*cHo*cKi

    f_K = -f_H          # flux as [mol/s]

    return f_H, f_K

def pumpVATP(cHi,cHo,Vm,T,p):

    deltaGATP = 20*p.R*T

    delG_H = p.R*T*np.log(cHo/cHi) - p.F*Vm  # free energy to move H+ out of cell

    delG_VATP = deltaGATP - delG_H   # free energy available to move H+ out of cell
    delG_pump = (delG_VATP/1000)
    delG = np.absolute(delG_pump)
    signG = np.sign(delG)

    if p.backward_pumps == False:

        alpha = p.alpha_V*tb.step(delG,p.halfmax_V,p.slope_V)
        f_H  = -alpha*cHi      #flux as [mol/s], scaled by concentrations in and out

    elif p.backward_pumps == True:

        alpha = signG*p.alpha_V*tb.step(delG,p.halfmax_V,p.slope_V)

        truth_forwards = signG == 1
        truth_backwards = signG == -1

        inds_forwards = (truth_forwards).nonzero()  # indices of forward-running cells
        inds_backwards = (truth_backwards).nonzero() # indices of backward-running cells

        f_H = np.zeros(len(cHi))

        f_H[inds_forwards]  = -alpha*cHi      #flux as [mol/s], scaled by concentrations in and out
        f_H[inds_backwards]  = -alpha*cHo

    return f_H

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

    # q = 0
    #
    # for conc,z in zip(concentrations,zs):
    #     q = q+ conc*z
    #
    netcharge = np.sum(p.F*vol*concentrations*zs,axis=0)

    return netcharge

def get_charge_density(concentrations,zs,p):

    """
    Calculates the charge density given ion concentrations in an array of spaces.

    Parameters
    --------------
    concentrations:  an array of array of concentration of ions in spaces [mol/m3]
    zs:              valence of each ion
    p:               Parameters object instance

    Returns
    -------------
    netcharge     the net charge density in spaces C/m3
    """

    # q = 0
    #
    # for conc,z in zip(concentrations,zs):
    #
    #     q = q + conc*z
    #
    # netcharge = p.F*q

    netcharge = np.sum(p.F*zs*concentrations, axis=0)

    return netcharge

def get_Vcell(self,cells,p):

    """
    Calculates the voltage in each cell from Poisson equation charge density

    Parameters
    --------------
    rho_cell:      an array of charge density in each cell space [C/m3]
    p:             an instance of the Parameters object

    Returns
    -------------
    v_cell          an array of voltages in each cell space  [V]

    """
    #
    # v_cell = self.rho_cells[:]*(p.rc**2)/(4*p.eo*80.0)


    #-------------------

    #   # # Poisson solver----------------------------------------------------------------

    # rho = np.zeros(len(cells.xypts))
    #
    # # map charge density in cells to cell and membrane points of the global lattice:
    # rho[cells.map_cell2ecm] = self.rho_cells
    # rho[cells.map_mem2ecm] = self.rho_cells[cells.mem_to_cells]
    # # rho = rho_env.ravel()

    rho = interp.griddata((cells.cell_centres[:,0],cells.cell_centres[:,1]),
                self.rho_cells,(cells.X,cells.Y))

    rho = np.nan_to_num(rho)

    rho = rho.ravel()

    # create a solution vector in the same shape as the source vector
    # V = np.zeros(len(cells.xypts))

    # modify the source charge distribution in line with electrostatic Poisson equation:
    fxy = -rho/(80*p.eo)

    # # modify the RHS of the equation to incorporate Dirichlet boundary conditions on the cell cluster:
    # fxy[cells.bound_pts_k] = 0
    # fxy[cells.bTop_kp] = (self.bound_V['T']/cells.delta**2)
    # fxy[cells.bL_kp] = (self.bound_V['L']/cells.delta**2)
    # fxy[cells.bR_kp] = (self.bound_V['R']/cells.delta**2)

    # Solve Poisson's electrostatic equation:
    V = np.dot(cells.Ainv,fxy)

    # if the boundary conditions set the outside of the matrix:
    V[cells.bBot_k] = self.bound_V['B']
    V[cells.bTop_k] = self.bound_V['T']
    V[cells.bL_k] = self.bound_V['L']
    V[cells.bR_k] = self.bound_V['R']

    v_cell = V[cells.map_cell2ecm]

    return v_cell

def get_Venv(self,cells,p):

    """
    Calculates the voltage in each extracellular space from Poisson equation charge density

    Parameters
    ---------------
    rho_ecm:        an array listing charge density in each ecm space [C/m3]
    p:              an instance of the Parameters object

    Returns
    ---------
    v_env           the voltage in each ecm space   [V]

    """

    # v_env = (self.rho_env*(p.cell_space**2))/(8*p.eo*80.0)
    # v_env[cells.ecm_bound_k] = 0

    # Poisson solver----------------------------------------------------------------

    # rho_env = self.rho_env[cells.core_k][:]
    # v_env = self.v_env[:]
    #
    # # modify the RHS of the equation to incorporate Dirichlet boundary conditions:
    #
    # rho_env[cells.bBot_kp] = rho_env[cells.bBot_kp] - (self.bound_V['B']/cells.delta**2)
    # rho_env[cells.bTop_kp] = rho_env[cells.bTop_kp] - (self.bound_V['T']/cells.delta**2)
    # rho_env[cells.bL_kp] = rho_env[cells.bL_kp] - (self.bound_V['L']/cells.delta**2)
    # rho_env[cells.bR_kp] = rho_env[cells.bR_kp] - (self.bound_V['R']/cells.delta**2)
    #
    # # Solve Poisson's electrostatic equation:
    # V = np.dot(cells.Ainv,-rho_env)
    #
    # # Re-pack the solution to the original vector:
    # v_env[cells.core_k] = V
    #
    # # set boundary values on the original vector
    # v_env[cells.bTop_k] = self.bound_V['T']
    # v_env[cells.bBot_k] = self.bound_V['B']
    # v_env[cells.bL_k] = self.bound_V['L']
    # v_env[cells.bR_k] = self.bound_V['R']

    # unravel the source vector
    rho = copy.copy(self.rho_env)
    # rho = rho_env.ravel()

    # create a solution vector in the same shape as the source vector
    V = np.zeros(rho.shape)

    # modify the source charge distribution in line with electrostatic Poisson equation:
    fxy = -rho/(80*p.eo)
    # fxy = -rho

    # # modify the RHS of the equation to incorporate Dirichlet boundary conditions:
    fxy[cells.bBot_kp] = (self.bound_V['B']/cells.delta**2)
    fxy[cells.bTop_kp] = (self.bound_V['T']/cells.delta**2)
    fxy[cells.bL_kp] = (self.bound_V['L']/cells.delta**2)
    fxy[cells.bR_kp] = (self.bound_V['R']/cells.delta**2)

    # Solve Poisson's electrostatic equation:
    V = np.dot(cells.Ainv,fxy)

    # if the boundary conditions set the outside of the matrix:
    V[cells.bBot_k] = self.bound_V['B']
    V[cells.bTop_k] = self.bound_V['T']
    V[cells.bL_k] = self.bound_V['L']
    V[cells.bR_k] = self.bound_V['R']

    return V

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
    isnans = np.isnan(vm)

    if isnans.any():  # if there's anything in the isubzeros matrix...
        raise BetseExceptionSimulation("Your simulation has become unstable. Please try a smaller time step,"
                                       "reduce gap junction radius, and/or reduce pump rate coefficients.")

def vertData(data, cells, p):
    """
    Interpolate data from midpoints to verts
    and sample it over a uniform grid.
    Produces data suitable for 2D mesh and
    streamline plots.

    Parameters
    -----------
    data          A numpy vector of data points on cell mids
    cells         An instance of the World object
    p             An instance of the Parameters object

    Returns
    ----------
    dat_grid      THe data sampled on a uniform grid

    """

    verts_data = np.dot(data,cells.matrixMap2Verts)
    plot_data = np.hstack((data,verts_data))

    dat_grid = interp.griddata((cells.plot_xy[:,0],cells.plot_xy[:,1]),plot_data,(cells.Xgrid,cells.Ygrid))
    dat_grid = np.nan_to_num(dat_grid)
    dat_grid = np.multiply(dat_grid,cells.cluster_mask)

    return dat_grid

def nernst_planck_flux(c, gcx, gcy, gvx, gvy,D,z,T,p):

    alpha = (D*z*p.q)/(p.kb*T)

    fx =  -D*gcx - alpha*gvx*c

    fy =  -D*gcy - alpha*gvy*c

    return fx, fy

def rk4(c,deltac,p):

    if p.method == 0:

        c2 = c + deltac*p.dt

    elif p.method == 1:

        k1 = deltac*p.dt

        k2 = (1/2)*k1*p.dt

        k3 = (1/2)*k2*p.dt

        k4 = k3*p.dt

        c2 = c + (p.dt/6)*(k1 + 2*k2 + 2*k3 + k4)

    return c2



#-----------------------------------------------------------------------------------------------------------------------
# WASTELANDS
#-----------------------------------------------------------------------------------------------------------------------

# def pumpNaKATP(cNai,cNao,cKi,cKo,sa,voli,volo,Vm,T,p,block):
#
#     """
#     Parameters
#     ----------
#     cNai            Concentration of Na+ inside the cell
#     cNao            Concentration of Na+ outside the cell
#     cKi             Concentration of K+ inside the cell
#     cKo             Concentration of K+ outside the cell
#     voli            Volume of the cell [m3]
#     volo            Volume outside the cell [m3]
#     Vm              Voltage across cell membrane [V]
#     p               An instance of Parameters object
#
#
#     Returns
#     -------
#     cNai2           Updated Na+ inside cell
#     cNao2           Updated Na+ outside cell
#     cKi2            Updated K+ inside cell
#     cKo2            Updated K+ outside cell
#     f_Na            Na+ flux (into cell +)
#     f_K             K+ flux (into cell +)
#     """
#     deltaGATP = 20*p.R*T
#
#     delG_Na = p.R*T*np.log(cNao/cNai) - p.F*Vm
#     delG_K = p.R*T*np.log(cKi/cKo) + p.F*Vm
#     delG_NaKATP = deltaGATP - (3*delG_Na + 2*delG_K)
#     delG_pump = (delG_NaKATP/1000)
#     delG = np.absolute(delG_pump)
#     signG = np.sign(delG)
#
#
#     if p.backward_pumps == False:
#
#         alpha = sa*block*p.alpha_NaK*tb.step(delG,p.halfmax_NaK,p.slope_NaK)
#
#         f_Na  = -alpha*cNai*cKo      #flux as [mol/s]   scaled to concentrations Na in and K out
#
#     elif p.backward_pumps == True:
#
#         alpha = sa*signG*block*p.alpha_NaK*tb.step(delG,p.halfmax_NaK,p.slope_NaK)
#
#         truth_forwards = signG == 1    # boolean array tagging forward-running pump cells
#         truth_backwards = signG == -1  # boolean array tagging backwards-running pump cells
#
#         inds_forwards = (truth_forwards).nonzero()  # indices of forward-running cells
#         inds_backwards = (truth_backwards).nonzero() # indices of backward-running cells
#
#         f_Na = np.zeros(len(cNai))
#
#         f_Na[inds_forwards]  = -alpha*cNai*cKo      #flux as [mol/s]   scaled to concentrations Na in and K out
#
#         f_Na[inds_backwards]  = -alpha*cNao*cKi      #flux as [mol/s]   scaled to concentrations K in and Na out
#
#     f_K = -(2/3)*f_Na          # flux as [mol/s]
#
#     if p.sim_ECM == True: # if extracellular spaces are included, concentrations need to be updated by matrix math
#
#         cNai2 = None
#         cNao2 = None
#         cKi2 = None
#         cKo2 = None
#
#     elif p.sim_ECM == False:
#
#         if p.method == 0:
#
#             dmol = f_Na*p.dt
#
#             cNai2 = cNai + dmol/voli
#             cNao2 = cNao - dmol/volo
#
#             cKi2 = cKi - (2/3)*dmol/voli
#             cKo2 = cKo + (2/3)*dmol/volo
#
#         elif p.method == 1:
#
#             k1 = alpha*cNai*cKo
#
#             k2 = alpha*(cNai+(1/2)*k1*p.dt)*cKo
#
#             k3 = alpha*(cNai+(1/2)*k2*p.dt)*cKo
#
#             k4 = alpha*(cNai+ k3*p.dt)*cKo
#
#             dmol = (p.dt/6)*(k1 + 2*k2 + 2*k3 + k4)
#
#             cNai2 = cNai - dmol/voli
#             cNao2 = cNao + dmol/volo
#
#             cKi2 = cKi + (2/3)*dmol/voli
#             cKo2 = cKo - (2/3)*dmol/volo
#
#
#     return cNai2,cNao2,cKi2,cKo2, f_Na, f_K


#---------------------------
    # Dc = rho*Dc
    #
    # grad_c = (cB - cA)/d
    #
    # c_mid = (cA + cB)/2
    #
    # grad_V = vBA/d
    #
    # flux = -Dc*grad_c - ((Dc*p.q*zc)/(p.kb*T))*c_mid*grad_V
    #
    # flux = rho*flux



#----------------------------

#def electrofuse(cA,cB,Dc,d,sa,vola,volb,zc,Vba,T,p,rho=1,ignoreECM = False):
#     """
#     Returns updated concentrations for electro-diffusion between two
#     connected volumes. Note for cell work, 'b' is 'inside', 'a' is outside, with
#     a positive flux moving from a to b. The voltage is defined as
#     Vb - Va (Vba), which is equivalent to Vmem.
#
#     This function defaults to regular diffusion if Vba == 0.0
#
#     This function takes numpy matrix values as input. All inputs must be matrices of
#     the same shape.
#
#     Parameters
#     ----------
#     cA          Initial concentration of c in region A [mol/m3] (out)
#     cB          Initial concentration of c in region B [mol/m3] (in)
#     Dc          Diffusion constant of c  [m2/s]
#     d           Distance between region A and region B [m]
#     sa          Surface area separating region A and B [m2]
#     vola        volume of region A [m3]
#     volb        volume of region B [m3]
#     zc          valence of ionic species c
#     Vba         voltage between B and A as Vb - Va  [V]
#     p           an instance of the Parameters class
#
#
#     Returns
#     --------
#     cA2         Updated concentration of cA in region A [mol/m3]
#     cB2         Updated concentration of cB in region B [mol/m3]
#     flux        Chemical flux magnitude between region A and B [mol/s]
#
#     """
#
#     # modify the diffusion constant by the membrane density
#     Dc = rho*Dc
#
#     alpha = (zc*Vba*p.F)/(p.R*T)
#
#     #volab = (vola + volb)/2
#     #qualityfactor = abs((Dc/d)*(sa/volab)*p.dt*alpha)   # quality factor should be <1.0 for stable simulations
#
#     deno = 1 - np.exp(-alpha)   # calculate the denominator for the electrodiffusion equation,..
#
#     izero = (deno==0).nonzero()     # get the indices of the zero and non-zero elements of the denominator
#     inotzero = (deno!=0).nonzero()
#
#     # initialize data matrices to the same shape as input data
#     dmol = np.zeros(deno.shape)
#     cA2 = np.zeros(deno.shape)
#     cB2 = np.zeros(deno.shape)
#     flux = np.zeros(deno.shape)
#
#     if p.method == 1:
#
#         k1 = np.zeros(deno.shape)
#         k2 = np.zeros(deno.shape)
#         k3 = np.zeros(deno.shape)
#         k4 = np.zeros(deno.shape)
#
#     if len(deno[izero]):   # if there's anything in the izero array:
#          # calculate the flux for those elements [mol/s]:
#         flux[izero] = -(sa[izero]/d[izero])*Dc[izero]*(cB[izero] - cA[izero])
#
#         if p.sim_ECM == True and ignoreECM == False:  # if we're simulating extracellular spaces, just calculate the flux
#
#             cA2 = None
#             cB2 = None
#
#         elif p.sim_ECM == False or ignoreECM == True:
#
#             if p.method == 0:
#
#                 #dmol[izero] = sa[izero]*p.dt*Dc[izero]*(cB[izero] - cA[izero])/d[izero]
#                 dmol[izero] = -p.dt*flux[izero]
#
#                 cA2[izero] = cA[izero] + (dmol[izero]/vola[izero])
#                 cB2[izero] = cB[izero] - (dmol[izero]/volb[izero])
#
#             elif p.method == 1:
#
#                 k1[izero] = -flux[izero]
#
#                 k2[izero] = sa[izero]*Dc[izero]*(cB[izero] - (cA[izero] + (1/2)*k1[izero]*p.dt))/d[izero]
#
#                 k3[izero] = sa[izero]*Dc[izero]*(cB[izero] - (cA[izero] + (1/2)*k2[izero]*p.dt))/d[izero]
#
#                 k4[izero] = sa[izero]*Dc[izero]*(cB[izero] - (cA[izero] + k3[izero]*p.dt))/d[izero]
#
#                 dmol[izero] = (p.dt/6)*(k1 + 2*k2 + 2*k3 + k4)
#
#                 cA2[izero] = cA[izero] + dmol[izero]/vola[izero]
#                 cB2[izero] = cB[izero] - dmol[izero]/volb[izero]
#
#
#     if len(deno[inotzero]):   # if there's any indices in the inotzero array:
#
#         # calculate the flux for those elements:
#         flux[inotzero] = -((sa[inotzero]*Dc[inotzero]*alpha[inotzero])/d[inotzero])*\
#                        ((cB[inotzero] - cA[inotzero]*np.exp(-alpha[inotzero]))/deno[inotzero])
#
#         if p.sim_ECM == True and ignoreECM == False:
#
#             cA2 = None
#             cB2 = None
#
#         elif p.sim_ECM == False or ignoreECM == True:
#
#             if p.method == 0:
#
#                 dmol[inotzero] = -flux[inotzero]*p.dt
#
#                 cA2[inotzero] = cA[inotzero] + (dmol[inotzero]/vola[inotzero])
#                 cB2[inotzero] = cB[inotzero] - (dmol[inotzero]/volb[inotzero])
#
#
#             elif p.method == 1:
#
#                 k1[inotzero] = -flux[inotzero]
#
#                 k2[inotzero] = ((sa[inotzero]*Dc[inotzero]*alpha[inotzero])/d[inotzero])*\
#                              (cB[inotzero] - (cA[inotzero] + (1/2)*k1[inotzero]*p.dt)*np.exp(-alpha[inotzero]))/deno[inotzero]
#
#                 k3[inotzero] = ((sa[inotzero]*Dc[inotzero]*alpha[inotzero])/d[inotzero])*\
#                              (cB[inotzero] - (cA[inotzero] + (1/2)*k2[inotzero]*p.dt)*np.exp(-alpha[inotzero]))/deno[inotzero]
#
#                 k4[inotzero] = ((sa[inotzero]*Dc[inotzero]*alpha[inotzero])/d[inotzero])*\
#                              (cB[inotzero] - (cA[inotzero] + k3[inotzero]*p.dt)*np.exp(-alpha[inotzero]))/deno[inotzero]
#
#                 dmol[inotzero] = (p.dt/6)*(k1[inotzero] + 2*k2[inotzero] + 2*k3[inotzero] + k4[inotzero])
#
#                 cA2[inotzero] = cA[inotzero] + dmol[inotzero]/vola[inotzero]
#                 cB2[inotzero] = cB[inotzero] - dmol[inotzero]/volb[inotzero]
#
#
#     return cA2, cB2, flux

#--------------------------------------------------

# def get_volt_ECM(q_cells,q_ecm,cells):
#
#     """
#     Makes use of the cell <--> ecm Maxwell Capacitance Matrix defined in World module
#     to calculate the distinct voltage of the cell interior, its extracellular space, and the
#     voltage difference across the membrane. The cell interior and the extracellular space
#     are taken to be two conductors (each with self-capacitance) connected by the capacitor of
#     the plasma membrane.
#
#     Parameters
#     -----------------
#     q_cells             An array storing the net charge inside each cell space
#     q_ecm               An array storing the net charge in each extracellular space
#     cells               An object storing World module properties
#
#     Returns
#     -----------------
#     v_mem               The voltage across the membrane as Vcell - Vecm
#     v_cells             The voltage in the intracellular space at the membrane
#     v_ecm              The voltage in the extracellular space at the membrane
#
#     """
#
#     v_cells = cells.Cinv_a*q_cells[cells.mem_to_cells] + cells.Cinv_b*q_ecm[cells.mem_to_ecm]
#     v_ecm = cells.Cinv_c*q_cells[cells.mem_to_cells] + cells.Cinv_d*q_ecm[cells.mem_to_ecm]
#     v_mem = v_cells - v_ecm
#
#     return v_mem, v_cells, v_ecm

    # def tissueInit(self,cells,p):
    #     """
    #     Runs initializations for all user-specified options that will be used in the main simulation.
    #
    #     """
    #
    #     # add channel noise to the model:
    #     self.Dm_cells[self.iK] = (p.channel_noise_level*self.channel_noise_factor + 1)*self.Dm_cells[self.iK]
    #
    #     # Initialize an array structure that will hold user-scheduled changes to membrane permeabilities:
    #     Dm_cellsA = np.asarray(self.Dm_cells)
    #     Dm_cellsER = np.asarray(self.Dm_er)
    #
    #     self.Dm_base = np.copy(Dm_cellsA) # make a copy that will serve as the unaffected values base
    #
    #     self.Dm_scheduled = np.copy(Dm_cellsA)
    #     self.Dm_scheduled[:] = 0
    #
    #     # Initialize an array structure that will hold dynamic voltage-gated channel changes to mem permeability:
    #     self.Dm_vg = np.copy(Dm_cellsA)
    #     self.Dm_vg[:] = 0
    #
    #     # Initialize an array structure that will hold dynamic calcium-gated channel changes to mem perms:
    #     self.Dm_cag = np.copy(Dm_cellsA)
    #     self.Dm_cag[:] = 0
    #
    #     self.Dm_er_base = np.copy(Dm_cellsER)
    #
    #     self.Dm_er_CICR = np.copy(Dm_cellsER)
    #     self.Dm_er_CICR[:] = 0
    #
    #     self.gj_block = np.ones(len(cells.gj_i))   # initialize the gap junction blocking vector to ones
    #
    #     self.cIP3 = np.zeros(len(cells.cell_i))  # initialize a vector to hold IP3 concentrations
    #     self.cIP3[:] = p.cIP3_to                 # set the initial concentration of IP3 from params file
    #
    #     self.cIP3_env = np.zeros(len(cells.cell_i))     # initialize IP3 concentration of the environment
    #     self.cIP3_env[:] = p.cIP3_to_env
    #
    #
    #     # Initialize target cell sets for dynamically gated channels from user options:
    #     if p.gated_targets == 'none':
    #         self.target_cells = np.zeros(len(cells.cell_i))
    #
    #     elif p.gated_targets == 'all':
    #         self.target_cells = np.ones(len(cells.cell_i))
    #
    #     elif p.gated_targets == 'random1':
    #         shuffle(cells.cell_i)
    #         trgt = cells.cell_i[0]
    #         self.target_cells = np.zeros(len(cells.cell_i))
    #         self.target_cells[trgt] = 1
    #
    #     elif p.gated_targets == 'random50':
    #         self.target_cells = np.random.random(len(cells.cell_i))
    #         self.target_cells = np.rint(self.target_cells)
    #
    #     elif isinstance(p.gated_targets,list):
    #         self.target_cells = np.zeros(len(cells.cell_i))
    #         self.target_cells[p.gated_targets] = 1
    #
    #     # allow for option to independently schedule an intervention to cells distinct from voltage gated:
    #     if p.scheduled_targets == 'none':
    #         self.scheduled_target_inds = []
    #
    #     elif p.scheduled_targets == 'all':
    #         self.scheduled_target_inds = cells.cell_i
    #
    #     elif p.scheduled_targets =='random1':
    #         shuffle(cells.cell_i)
    #         trgt2 = cells.cell_i[0]
    #         self.scheduled_target_inds = [trgt2]
    #
    #     elif p.scheduled_targets == 'random50':
    #         shuffle(cells.cell_i)
    #         halflength = int(len(cells.cell_i)/2)
    #         self.scheduled_target_inds = [cells.cell_i[x] for x in range(0,halflength)]
    #
    #     elif isinstance(p.scheduled_targets, list):
    #         # self.scheduled_target_inds = np.zeros(len(cells.cell_i))
    #         # self.scheduled_target_inds[p.scheduled_targets] = 1
    #         self.scheduled_target_inds = p.scheduled_targets
    #
    #
    #
    #     loggers.log_info('This world contains '+ str(cells.cell_number) + ' cells.')
    #     loggers.log_info('Each cell has an average of '+ str(round(cells.average_nn,2)) + ' nearest-neighbours.')
    #     loggers.log_info('You are running the ion profile: '+ p.ion_profile)
    #
    #
    #     loggers.log_info('Ions in this simulation: ' + str(self.ionlabel))
    #     loggers.log_info('If you have selected features using other ions, they will be ignored.')


    # def update_environment(self,cells,p,t,i):
    #
    #     # calculate voltage differences between ecm <---> env nearest environmental neighbours:
    #     self.v_ec2env =  self.v_env - self.v_ecm[cells.bflags_ecm]
    #
    #     # electrodiffuse between ecm <---> env and simultaneously update concentrations
    #     self.cc_ecm[i][cells.bflags_ecm],self.cc_env[i],f_env = electrofuse(self.cc_ecm[i][cells.bflags_ecm],
    #             self.cc_env[i], self.D_env_juncs[i],cells.len_ecm_junc[cells.bflags_ecm],
    #             self.ec2ec_sa[cells.bflags_ecm], cells.ecm_vol[cells.bflags_ecm],self.env_vol,
    #             self.zs[i], self.v_ec2env, self.T, p,ignoreECM=True)
    #
    #     self.rho_env = get_charge_density(self.cc_ecm, self.z_array_ecm, p)
    #     self.v_ecm = get_Vecm(self.rho_env,p)
    #
    #     self.rho_env = get_charge_density(self.cc_env, self.z_array_env, p)
    #     self.v_env = get_Vecm(self.rho_env,p)
    #
    #     self.fluxes_env[i] = f_env  # store ecm junction flux to the environment for this ion
    #
    #     # if p.ecm_bound_open == True: # if the geometry is open to the environment
    #     #     self.cc_env[i] = self.cc_env_bulk[i]  # make outer ecm spaces equal to bulk environment
    #     #     self.v_env[:] = 0
    #
    #     if p.scheduled_options['extV'] != 0 and p.run_sim == True:
    #         self.dyna.externalVoltage(self,cells,p,t,self.v_env)      # apply a voltage to the environmental connectors





