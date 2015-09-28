#!/usr/bin/env python3
# Copyright 2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

import numpy as np
import math
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
from betse.science import finitediff as fd
import matplotlib.pyplot as plt
from betse.exceptions import BetseExceptionSimulation
from betse.util.io import loggers
import time

class Simulator(object):
    """
    Contains the main routines used in the simulation of networked cell bioelectrical activity.
    All methods are based on matrix mathematics and are implemented in Numpy for speed.

    Methods
    -------
    fileInit(p)                 Prepares save paths for initialization and simulation runs

    baseInit(cells,p)           Prepares data structures for a cell-only simulation

    baseInit_ECM(cells,p)      Prepares data structures for a simulation with extracellular spaces

    tissueInit(cells,p)         Prepares data structures pertaining to tissue profiles and dynamic activity

    runSim(cells,p)            Runs and saves a simulation (init or sim) for world with cells only

    runSim_ECM(cells, p)        Runs and saves a simulation (init or sim) for world with cells and full environment

    update_V_ecm(cells,p,t)    For sims with environmental spaces, gets charge densities in cells and environment
                                and calculates respective voltages.

    update_C_ecm(ion_i,flux, cells, p)     For sims with full environment, updates concentration of ion with index
                                            ion_i in cell and environment for a flux leaving the cell.

    Hplus_electrofuse_ecm(cells,p,t)        For sims with full environment, updates H+ concentrations in cell and
                                            environment, which are further influenced by the bicarbonate buffer action.

    Hplus_HKATP_ecm(cells,p,t)              For sims with full environment, runs the HKATPase pump and updates H+ and K+
                                            concentrations under the action of the bicarbonate buffer system.

    Hplus_VATP_ecm(cells,p,t)               For sims with full environment, runs the VATPase pump and updates H+
                                            concentrations in cell and environment, which are further influenced by the
                                            bicarbonate buffer action.

    update_gj(cells,p,t,i)                  Calculates the voltage gradient between two cells, the gating character of
                                            gap junctions, and updates concentration for ion 'i' and voltage of each
                                            cell after electrodiffusion of ion 'i' between gap junction connected cells.

    update_ecm(cells,p,t,i)                 Updates the environmental spaces by calculating electrodiffusive transport
                                            of ion 'i'.

    update_er(cells,p,t)                    Updates concentrations in the endoplasmic reticulum.

    update_dye(cells,p,t)                   Updates concentration of morphogen in cells and environment by calculating
                                            electrodiffusive transport across membranes, between gj connected cells, and
                                            through environmental spaces.

    update_IP3(cells,p,t)                   Updates concentration of IP3 in cells and environment by calculating
                                            electrodiffusive transport across membranes, between gj connected cells, and
                                            through environmental spaces.

    get_Efield(cells,p)                     Calculates electric fields in cells and environment.

    get_Bfield(cells,p)                     Calculates magnetic fields in cells and environment.

    get_current(cells,p)                    Calculates currents in cells and environment.

    getFlow(cells,p)                        Calculates electroosmotic flows in cells and environment.

    eosmosis(cells,p)                       Calculates lateral movement of membrane pumps and channels via tangential
                                            forces exerted by endogenous electric fields and electroosmotic flows.

    get_ion(label)                          Supply the ion name as a string input ('Na', 'K', 'Ca', etc) and it outputs
                                            the sim index of that ion type.

    initDenv(cells,p)                       Initializes the environmental diffusion matrix and corresponding weight
                                            matrices, including tight and adherin junctions.





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

        # initialize vectors for electroosmosis in the cell collection:
        self.u_cells_x = np.zeros(len(cells.cell_i))
        self.u_cells_y = np.zeros(len(cells.cell_i))
        self.P_cells = np.zeros(len(cells.cell_i))

        if p.gj_flux_sensitive == True:

            self.gj_rho = np.zeros(len(cells.nn_i))

        else:

            self.gj_rho = 0


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
                        # self.Dm_er.append(p.Dm_Ca)

                    if name == 'M' and p.ions_dict['Ca'] == 1:

                        self.cM_er = np.zeros(len(cells.cell_i))
                        self.cM_er[:]=p.cCa_er

                        self.zM_er = np.zeros(len(cells.cell_i))
                        self.zM_er[:]=p.z_M

                        self.cc_er.append(self.cM_er)
                        self.z_er.append(p.z_M)
                        self.z_array_er.append(self.zM_er)
                        # self.Dm_er.append(p.Dm_M)

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

            # self.cH_cells = np.zeros(len(cells.cell_i))
            # cH_cells[:]=p.cH_cell
            self.pH_cell = 6.1 + np.log10(self.cM_cells/self.cHM_cells)
            self.cH_cells = (10**(-self.pH_cell))  # units mmol/L

            # self.cH_env = np.zeros(len(cells.cell_i))
            # cH_env[:]=p.cH_env
            self.pH_env = 6.1 + np.log10(self.cM_env/self.cHM_env)
            self.cH_env = (10**(-self.pH_env)) # units mmol/L

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

        self.Dm_er = np.zeros((2,len(cells.cell_i)))
        self.Dm_er[0,:] = p.Dm_Ca
        self.Dm_er[1,:] = p.Dm_M

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

        if p.VATPase_dyn == True and p.global_options['VATP_block'] != 0:
            self.VATP_block = np.ones(len(cells.cell_i))  # initialize HKATP blocking vector
        else:
            self.VATP_block = 1

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

            self.cIP3_env = np.zeros(len(cells.cell_i))     # initialize IP3 concentration of the environment
            self.cIP3_env[:] = p.cIP3_to_env

        if p.voltage_dye == True:

            self.cDye_cell = np.zeros(len(cells.cell_i))   # initialize voltage sensitive dye array for cell and env't
            self.cDye_cell[:] = p.cDye_to_cell

            self.Dye_flux_gj = np.zeros(len(cells.nn_i))
            self.Dye_flux_mem = np.zeros(len(cells.mem_i))

            self.cDye_env = np.zeros(len(cells.cell_i))     # initialize Dye concentration in the environment
            self.cDye_env[:] = p.cDye_to
        #
        self.z_array = np.asarray(self.z_array)
        self.D_gj = np.asarray(self.D_gj)

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
        self.Dtj_rel = []               # relative diffusion constants for ions across tight junctions

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

        # Electroosmosis Initialization:

        # initialize vectors for env flow (note enhanced data type!):
        self.u_env_x = np.zeros(cells.grid_obj.u_shape)
        self.u_env_y = np.zeros(cells.grid_obj.v_shape)

        self.P_env = np.zeros(cells.grid_obj.cents_shape)

        # initialize vectors for electroosmosis in the cell collection wrt each gap junction (note data type!):
        self.u_cells_x = np.zeros(len(cells.cell_i))
        self.u_cells_y = np.zeros(len(cells.cell_i))
        self.P_cells = np.zeros(len(cells.cell_i))


        if p.gj_flux_sensitive == True:

            self.gj_rho = np.zeros(len(cells.nn_i))

        else:

            self.gj_rho = 0

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
                    self.Dtj_rel.append(p.Dtj_rel[name])

                    self.fluxes_gj_x.append(self.flx_gj_i)
                    self.fluxes_gj_y.append(self.flx_gj_i)
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
                        self.z_array_er.append(self.zCa_er)
                        # self.Dm_er.append(p.Dm_Ca)

                    if name == 'M' and p.ions_dict['Ca'] == 1:

                        self.cM_er = np.zeros(len(cells.cell_i))
                        self.cM_er[:]=p.cCa_er

                        self.zM_er = np.zeros(len(cells.cell_i))
                        self.zM_er[:]=p.z_M

                        self.cc_er.append(self.cM_er)
                        self.z_er.append(p.z_M)
                        self.z_array_er.append(self.zM_er)
                        # self.Dm_er.append(p.Dm_M)

        # Do H+ separately as it's complicated by the buffer
        # initialize the carbonic acid for the carbonate buffer

        if p.ions_dict['H'] == 1:

            i = i+ 1

            self.iH = i

            self.ionlabel[self.iH] = 'protons'

            self.movingIons.append(self.iH)

            self.cHM_cells = np.zeros(len(cells.cell_i))
            self.cHM_cells[:] = 0.03*p.CO2

            # self.cHM_env = np.zeros(len(cells.xypts))
            # self.cHM_env[:] = 0.03*p.CO2
            self.cHM_env = 0.03*p.CO2

            # self.cH_cells = np.zeros(len(cells.cell_i))

            self.pH_cell = 6.1 + np.log10(self.cM_cells/self.cHM_cells)
            self.cH_cells = (10**(-self.pH_cell)) # units mmol/L

            # self.cH_env = np.zeros(len(cells.xypts))
            self.pH_env = 6.1 + np.log10(self.cM_env/self.cHM_env)
            self.cH_env = (10**(-self.pH_env)) # units mmol/L

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


            self.c_env_bound.append(p.env_concs['H'])

            self.cc_cells.append(self.cH_cells)
            self.cc_env.append(self.cH_env)

            self.zs.append(p.z_H)
            self.z_array.append(self.zH)
            self.z_array_env.append(self.zH2)
            self.Dm_cells.append(DmH)
            self.D_env.append(DenvH)
            self.D_gj.append(DgjH)
            self.D_free.append(p.Do_H)
            self.Dtj_rel.append(p.Dtj_rel['H'])

            self.fluxes_gj_x.append(self.flx_gj_i)
            self.fluxes_gj_y.append(self.flx_gj_i)
            self.fluxes_mem.append(self.flx_mem_i)
            self.fluxes_env_x.append(self.flx_env_i)
            self.fluxes_env_y.append(self.flx_env_i)

        #-------------------------------------------------------------------------------------------------------

        # Define the diffusion matrix for the endoplasmic reticulum:
        self.Dm_er = np.zeros((2,len(cells.cell_i)))
        self.Dm_er[0,:] = p.Dm_Ca
        self.Dm_er[1,:] = p.Dm_M

        self.vm_to = np.zeros(len(cells.cell_i))
        self.v_er = np.zeros(len(cells.cell_i))

        if p.global_options['NaKATP_block'] != 0:

            self.NaKATP_block = np.ones(len(cells.cell_i))  # initialize NaKATP blocking vector

        else:
            self.NaKATP_block = 1

        if p.HKATPase_dyn == True and p.global_options['HKATP_block'] != 0:
            self.HKATP_block = np.ones(len(cells.mem_i))  # initialize HKATP blocking vector
        else:
            self.HKATP_block = 1

        if p.VATPase_dyn == True and p.global_options['VATP_block'] != 0:
            self.VATP_block = np.ones(len(cells.mem_i))  # initialize HKATP blocking vector
        else:
            self.VATP_block = 1

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
        self.cc_er = np.asarray(self.cc_er)
        self.zs = np.asarray(self.zs)
        self.z_array = np.asarray(self.z_array)
        self.z_array_env = np.asarray(self.z_array_env)
        self.z_array_er = np.asarray(self.z_array_er)
        self.Dm_cells = np.asarray(self.Dm_cells)
        self.Dm_er = np.asarray(self.Dm_er)
        self.D_env = np.asarray(self.D_env)
        self.D_free = np.asarray(self.D_free)
        self.D_gj = np.asarray(self.D_gj)

        self.fluxes_gj_x  = np.asarray(self.fluxes_gj_x)
        self.fluxes_gj_y  = np.asarray(self.fluxes_gj_y)
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

        # initialize the environmental diffusion matrix:
        self.initDenv(cells,p)

    def tissueInit(self,cells,p):

        if p.sim_ECM == True:
            #  Initialize diffusion constants for the extracellular transport:
            self.initDenv(cells,p)

        self.dyna = Dynamics(self,cells,p)   # create the tissue dynamics object
        self.dyna.tissueProfiles(self,cells,p)  # initialize all tissue profiles

        if p.sim_ECM == True:
            # create a copy-base of the environmental junctions diffusion constants:
            self.D_env_base = copy.copy(self.D_env)

        # add channel noise to the model:
        self.Dm_cells[self.iK] = (p.channel_noise_level*self.channel_noise_factor + 1)*self.Dm_cells[self.iK]

        # Initialize an array structure that will hold user-scheduled changes to membrane permeabilities:
        Dm_cellsA = np.asarray(self.Dm_cells)
        Dm_cellsER = np.copy(self.Dm_er)

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

        self.Dm_morpho = np.copy(Dm_cellsA)
        self.Dm_morpho[:] = 0

        if p.global_options['gj_block'] != 0:

            self.gj_block = np.ones(len(cells.nn_i))   # initialize the gap junction blocking vector to ones

        else:

            self.gj_block = 1

        if p.scheduled_options['IP3'] != 0 or p.Ca_dyn == True:

            self.cIP3 = np.zeros(len(cells.cell_i))  # initialize a vector to hold IP3 concentrations
            self.cIP3[:] = p.cIP3_to                 # set the initial concentration of IP3 from params file

            self.IP3_flux_x_gj = np.zeros(len(cells.nn_i))
            self.IP3_flux_y_gj = np.zeros(len(cells.nn_i))
            self.IP3_flux_mem = np.zeros(len(cells.mem_i))

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

            self.Dye_flux_x_gj = np.zeros(len(cells.nn_i))
            self.Dye_flux_y_gj = np.zeros(len(cells.nn_i))
            self.Dye_flux_mem = np.zeros(len(cells.mem_i))

            self.c_dye_bound = p.cDye_to    # concentration of dye at the global boundaries

            if p.Dye_target_channel != 'None':

                # get the ion index of the target ion:
                self.dye_target = self.get_ion(p.Dye_target_channel)

                # make a copy of the appropriate ion list
                self.Dm_mod_dye = np.copy(self.Dm_cells[self.dye_target][:])

            else:

                self.dye_target = None

            if p.sim_ECM == True:

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

        # initialize time-storage vectors for electroosmosis:
        self.P_cells_time = []
        self.u_cells_x_time = []
        self.u_cells_y_time = []

        if p.voltage_dye == True:

            self.cDye_time = []    # retains voltage-sensitive dye concentration as a function of time
            self.Dye_flux_x_gj_time = []
            self.Dye_flux_y_gj_time = []
            self.Dye_flux_mem_time = []

        if p.scheduled_options['IP3'] != 0 or p.Ca_dyn == True:

            self.cIP3_time = []    # retains IP3 concentration as a function of time
            self.IP3_flux_x_gj_time = []
            self.IP3_flux_y_gj_time = []
            self.IP3_flux_mem_time = []

        # gap junction specific arrays:
        self.id_gj = np.ones(len(cells.nn_i))  # identity array for gap junction indices...
        self.gjopen = np.ones(len(cells.nn_i))   # holds gap junction open fraction for each gj
        self.gjl = np.zeros(len(cells.nn_i))    # gj length for each gj
        self.gjl[:] = cells.nn_len

        self.rho_channel = np.ones(len(cells.mem_i))

        # get the net, unbalanced charge and corresponding voltage in each cell:
        self.update_V_ecm(cells,p,0)

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
                 # reassign the history-saving vm
            self.vm_to = self.vm[:]

            if p.Ca_dyn ==1 and p.ions_dict['Ca'] == 1:

                self.dcc_ER = (self.cc_er - self.cc_er_to)/p.dt
                # self.cc_er_to = copy.deepcopy(self.cc_er)
                self.cc_er_to = self.cc_er[:]

            # calculate the values of scheduled and dynamic quantities (e.g. ion channel multipliers):
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
            self.update_V_ecm(cells,p,t)

            if p.ions_dict['Ca'] == 1:
                # run the calcium ATPase membrane pump:
                fCaATP = pumpCaATP(self.cc_cells[self.iCa],self.cc_env[self.iCa],self.vm,self.T,p)

                # update concentrations in the cell:
                self.cc_cells[self.iCa] = self.cc_cells[self.iCa] + fCaATP*(cells.cell_sa/cells.cell_vol)*p.dt

                # store the transmembrane flux for this ion
                self.fluxes_mem[self.iCa] = fCaATP[cells.mem_to_cells]

                # recalculate the net, unbalanced charge and voltage in each cell:
                self.update_V_ecm(cells,p,t)

                if p.Ca_dyn ==1:

                    # run the calcium ATPase endoplasmic reticulum pump:
                    fCaATP_ER = pumpCaER(self.cc_er[0],self.cc_cells[self.iCa],self.v_er,self.T,p)



                    # update calcium concentrations in the ER and cell:
                    self.cc_er[0] = self.cc_er[0] + fCaATP_ER*((cells.cell_sa)/(p.ER_vol*cells.cell_vol))*p.dt
                    self.cc_cells[self.iCa] = self.cc_cells[self.iCa] - fCaATP_ER*(cells.cell_sa/cells.cell_vol)*p.dt

                    # recalculate the net, unbalanced charge and voltage in each cell:
                    self.update_V_ecm(cells,p,t)

                    q_er = get_charge(self.cc_er,self.z_array_er,p.ER_vol*cells.cell_vol,p)
                    v_er_o = get_volt(q_er,p.ER_sa*cells.cell_sa,p)
                    self.v_er = v_er_o - self.vm

            if p.ions_dict['H'] == 1:

                # electrofuse the H+ ion between the cytoplasm and the environment
                f_H1 = electroflux(self.cc_env[self.iH],self.cc_cells[self.iH],self.Dm_cells[self.iH],self.tm,
                    self.zs[self.iH],self.vm,self.T,p)

                # update the anion rather than H+, assuming that the bicarbonate buffer is working:
                self.cc_cells[self.iM] = self.cc_cells[self.iM] - f_H1*(cells.cell_sa/cells.cell_vol)*p.dt
                # self.cc_env[self.iM] = self.cc_env[self.iM] + f_H1*(cells.cell_sa/p.vol_env)*p.dt

                self.fluxes_mem[self.iH] = f_H1[cells.mem_to_cells]

                # Calculate the new pH and H+ concentration:
                self.pH_cell = 6.1 + np.log10(self.cc_cells[self.iM]/self.cHM_cells)
                self.cc_cells[self.iH] = 10**(-self.pH_cell)

                # recalculate the net, unbalanced charge and voltage in each cell:
                self.update_V_ecm(cells,p,t)

                if p.HKATPase_dyn == 1:

                    # if HKATPase pump is desired, run the H-K-ATPase pump:
                    f_H2, f_K2 = pumpHKATP(self.cc_cells[self.iH],self.cc_env[self.iH],self.cc_cells[self.iK],
                        self.cc_env[self.iK],self.vm,self.T,p,self.HKATP_block)

                    # update the concentration in cells (assume environment steady and constant supply of ions)
                    self.cc_cells[self.iM] = self.cc_cells[self.iM] - f_H2*(cells.cell_sa/cells.cell_vol)*p.dt
                    self.cc_cells[self.iK] = self.cc_cells[self.iK] + f_K2*(cells.cell_sa/cells.cell_vol)*p.dt

                    # store fluxes for this pump:
                    self.fluxes_mem[self.iH] = self.fluxes_mem[self.iH] + f_H2[cells.mem_to_cells]
                    self.fluxes_mem[self.iK] = self.fluxes_mem[self.iK] + f_K2[cells.mem_to_cells]

                    # Calculate the new pH and H+ concentration:
                    self.pH_cell = 6.1 + np.log10(self.cc_cells[self.iM]/self.cHM_cells)
                    self.cc_cells[self.iH] = 10**(-self.pH_cell)

                    # recalculate the net, unbalanced charge and voltage in each cell:
                    self.update_V_ecm(cells,p,t)

                if p.VATPase_dyn == 1:

                     # if HKATPase pump is desired, run the H-K-ATPase pump:
                    f_H3 = pumpVATP(self.cc_cells[self.iH],self.cc_env[self.iH],self.vm,self.T,p,self.VATP_block)

                    self.cc_cells[self.iM] = self.cc_cells[self.iM] - f_H3*(cells.cell_sa/cells.cell_vol)*p.dt

                    self.fluxes_mem[self.iH]  = self.fluxes_mem[self.iH] + f_H3[cells.mem_to_cells]

                    # Calculate the new pH and H+ concentration:
                    self.pH_cell = 6.1 + np.log10(self.cc_cells[self.iM]/self.cHM_cells)
                    self.cc_cells[self.iH] = 10**(-self.pH_cell)

                    # recalculate the net, unbalanced charge and voltage in each cell:
                    self.update_V_ecm(cells,p,t)

            # electro-diffuse all ions (except for proteins, which don't move) across the cell membrane:

            shuffle(self.movingIons)

            for i in self.movingIons:

                # electrodiffusion of ion between cell and extracellular matrix
                f_ED = electroflux(self.cc_env[i],self.cc_cells[i],self.Dm_cells[i],self.tm,self.zs[i],self.vm,self.T,p)

                # update ion due to transmembrane flux:
                self.cc_cells[i] = self.cc_cells[i] + f_ED*(cells.cell_sa/cells.cell_vol)*p.dt

                self.fluxes_mem[i] = self.fluxes_mem[i] + f_ED[cells.mem_to_cells]

                # # recalculate the net, unbalanced charge and voltage in each cell:
                self.update_V_ecm(cells,p,t)

                self.update_gj(cells,p,t,i)

            # calculate electric fields:

            self.get_Efield(cells, p)

            # calculate electroosmotic flow:
            if p.base_eosmo == True and cells.lapGJinv != 0:
                self.getFlow(cells,p)

            if p.scheduled_options['IP3'] != 0 or p.Ca_dyn == True:
                # determine flux through gap junctions for IP3:

                self.update_IP3(cells,p,t)

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
                self.update_V_ecm(cells,p,t)

                q_er = get_charge(self.cc_er,self.z_array_er,p.ER_vol*cells.cell_vol,p)
                v_er_o = get_volt(q_er,p.ER_sa*cells.cell_sa,p)
                self.v_er = v_er_o

            # if p.voltage_dye=1 electrodiffuse voltage sensitive dye between cell and environment
            if p.voltage_dye ==1:

                self.update_dye(cells,p,t)

            if p.dynamic_noise == 1 and p.ions_dict['P']==1:
                # add a random walk on protein concentration to generate dynamic noise:
                self.protein_noise_factor = p.dynamic_noise_level*(np.random.random(len(cells.cell_i)) - 0.5)
                self.cc_cells[self.iP] = self.cc_cells[self.iP]*(1+ self.protein_noise_factor)
                self.update_V_ecm(cells,p,t)

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

                if p.base_eosmo == True:

                    self.P_cells_time.append(self.P_cells[:])
                    self.u_cells_x_time.append(self.u_cells_x[:])
                    self.u_cells_y_time.append(self.u_cells_y[:])

                self.gjopen_time.append(self.gjopen[:])

                self.time.append(t)

                if p.scheduled_options['IP3'] != 0 or p.Ca_dyn == True:

                    self.cIP3_time.append(self.cIP3[:])

                    self.IP3_flux_x_gj_time.append(self.IP3_flux_x_gj[:])
                    self.IP3_flux_y_gj_time.append(self.IP3_flux_y_gj[:])

                    self.IP3_flux_mem_time.append(self.IP3_flux_mem[:])

                if p.voltage_dye ==1:

                    self.cDye_time.append(self.cDye_cell[:])

                    self.Dye_flux_x_gj_time.append(self.Dye_flux_x_gj[:])
                    self.Dye_flux_y_gj_time.append(self.Dye_flux_y_gj[:])

                    self.Dye_flux_mem_time.append(self.Dye_flux_mem[:])

                if p.Ca_dyn == 1 and p.ions_dict['Ca']==1:

                    self.cc_er_time.append(np.copy(self.cc_er[:]))


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
            final_pH = -np.log10(np.mean((self.cc_time[-1][self.iH])))
            loggers.log_info('Final average cell pH '+ str(np.round(final_pH,2)))

            final_pH_env = -np.log10(np.mean((self.cc_env_time[-1][self.iH])))
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

        self.f_gj_x_time = []      # stores the gj fluxes for each ion at each time

        self.Igj_time = []      # current for each gj at each time

        self.I_gj_x_time = []    # initialize gap junction current data storage
        self.I_gj_y_time = []    # initialize gap junction current data storage
        self.I_env_x_time = []   # initialize environmental matrix data storage
        self.I_env_y_time = []   # initialize environmental matrix data storage
        self.I_mem_time = []   # initialize membrane matrix data storage

        self.cc_er_time = []   # retains er concentrations as a function of time
        self.cIP3_time = []    # retains cellular ip3 concentrations as a function of time
        self.cIP3_env_time = []

        self.efield_gj_x_time = []   # matrices storing smooth electric field in gj connected cells
        self.efield_gj_y_time = []

        self.efield_ecm_x_time = []   # matrices storing smooth electric field in ecm
        self.efield_ecm_y_time = []

        # initialize time-storage vectors for electroosmotic data:
        self.P_env_time = []
        self.u_env_x_time = []
        self.u_env_y_time = []

        self.u_cells_x_time = []
        self.u_cells_y_time = []
        self.P_cells_time = []

        self.vm_Matrix = [] # initialize matrices for resampled data sets (used in smooth plotting and streamlines)
        vm_dato = np.zeros(len(cells.mem_i))
        dat_grid_vm = vertData(vm_dato,cells,p)
        self.vm_Matrix.append(dat_grid_vm[:])

        if p.Ca_dyn == True:
            self.cc_er_to = np.copy(self.cc_er[:])

        if p.voltage_dye == True:

            self.Dye_flux_env_x_time = []
            self.Dye_flux_env_y_time = []
            self.Dye_flux_x_gj_time = []
            self.Dye_flux_y_gj_time = []
            self.Dye_flux_mem_time = []
            self.cDye_time = []    # retains voltage-sensitive dye concentration as a function of time
            self.cDye_env_time = []

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

            loggers.log_info('Your simulation (with extracellular spaces) is running from '+ str(0) +
                             ' to '+ str(round(p.sim_tsteps*p.dt,3))
                         + ' seconds of in-world time.')

        else:
            loggers.log_info('Your initialization (with extracellular spaces) is running from '+ str(0) +
                             ' to '+ str(round(p.init_tsteps*p.dt,3))
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

            #-----------------PUMPS-------------------------------------------------------------------------------------

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

                    f_Ca_ER = pumpCaER(self.cc_er[0],self.cc_cells[self.iCa],self.v_er,self.T,p)

                    # update calcium concentrations in the ER and cell:
                    self.cc_er[0] = self.cc_er[0] + f_Ca_ER*((cells.cell_sa)/(p.ER_vol*cells.cell_vol))*p.dt
                    self.cc_cells[self.iCa] = self.cc_cells[self.iCa] - f_Ca_ER*(cells.cell_sa/cells.cell_vol)*p.dt

                    # recalculate the net, unbalanced charge and voltage in each cell:
                    self.update_V_ecm(cells,p,t)

                    # calculate the net, unbalanced charge and voltage in the endoplasmic reticulum:
                    q_er = get_charge(self.cc_er,self.z_array_er,p.ER_vol*cells.cell_vol,p)
                    v_er_o = get_volt(q_er,p.ER_sa*cells.cell_sa,p)

                    self.v_er = v_er_o - self.v_cell

            if p.ions_dict['H'] == 1:

                self.Hplus_electrofuse_ecm(cells,p,t)

                if p.HKATPase_dyn == 1:

                    self.Hplus_HKATP_ecm(cells,p,t)

                if p.VATPase_dyn == 1:

                    self.Hplus_VATP_ecm(cells,p,t)

            #----------------ELECTRODIFFUSION---------------------------------------------------------------------------

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

                # update flux between cells due to gap junctions
                self.update_gj(cells,p,t,i)

                # # recalculate the net, unbalanced charge and voltage in each cell:
                self.update_V_ecm(cells,p,t)

            self.get_Efield(cells,p)

            if p.base_eosmo == True:

                self.getFlow(cells,p)

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
            if p.sim_eosmosis == True and cells.gradMem != None:

                self.eosmosis(cells,p)    # modify membrane pump and channel density according to Nernst-Planck

            if t in tsamples:

                # #
                self.get_current(cells,p)   # get the current in the gj network connection of cells

                # add the new concentration and voltage data to the time-storage matrices:

                self.efield_gj_x_time.append(self.E_gj_x[:])

                self.efield_gj_y_time.append(self.E_gj_y[:])

                self.efield_ecm_x_time.append(self.E_env_x[:])

                self.efield_ecm_y_time.append(self.E_env_y[:])

                concs = np.copy(self.cc_cells[:])
                concs.tolist()
                self.cc_time.append(concs)
                concs = None

                ecmsc = np.copy(self.cc_env[:])
                ecmsc.tolist()
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

                if p.base_eosmo == True:

                    self.P_env_time.append(np.float64(self.P_env[:]))
                    self.u_env_x_time.append(self.u_at_c[:])
                    self.u_env_y_time.append(self.v_at_c[:])

                    self.P_cells_time.append(np.float64(self.P_cells[:]))
                    self.u_cells_x_time.append(np.float64(self.u_cells_x))
                    self.u_cells_y_time.append(np.float64(self.u_cells_y))

                # calculate interpolated verts and midpoint data for Vmem:
                dat_grid_vm = vertData(self.vm[:],cells,p)

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

                    self.cDye_env_time.append(self.cDye_env[:])

                    self.Dye_flux_x_gj_time.append(self.Dye_flux_x_gj[:])
                    self.Dye_flux_y_gj_time.append(self.Dye_flux_y_gj[:])

                    self.Dye_flux_env_x_time.append(self.Dye_flux_env_x[:])
                    self.Dye_flux_env_y_time.append(self.Dye_flux_env_y[:])

                if p.Ca_dyn == 1 and p.ions_dict['Ca']==1:
                    self.cc_er_time.append(np.copy(self.cc_er[:]))


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
            final_pH = -np.log10(np.mean((self.cc_time[-1][self.iH])))
            loggers.log_info('Final average cell pH '+ str(np.round(final_pH,2)))

            final_pH_ecm = -np.log10(np.mean((self.cc_env_time[-1][self.iH])))
            loggers.log_info('Final extracellular pH '+ str(np.round(final_pH_ecm,2)))

        if p.scheduled_options['IP3'] != 0 or p.Ca_dyn == True:

            IP3_env_final = np.mean(self.cIP3_env)
            IP3_cell_final = np.mean(self.cIP3)
            loggers.log_info('Final extracellular IP3 concentration: ' + str(np.round(IP3_env_final,6)) + ' mmol/L')
            loggers.log_info('Final average IP3 concentration in cells: ' + str(np.round(IP3_cell_final,6)) + ' mmol/L')

        if p.Ca_dyn == 1 and p.ions_dict['Ca'] == 1:

            endconc_er = np.round(np.mean(self.cc_er[0]),6)
            label = self.ionlabel[self.iCa]
            concmess = 'Final average ER concentration of'+ ' '+ label + ': '
            loggers.log_info(concmess + str(endconc_er) + ' mmol/L')

        if p.voltage_dye ==1:
            dye_ecm_final = np.mean(self.cDye_env)
            dye_cell_final = np.mean(self.cDye_cell)
            loggers.log_info('Final extracellular morphogen concentration: '+ str(np.round(dye_ecm_final,6))
                             + ' mmol/L')
            loggers.log_info('Final average morphogen concentration in cells: ' +  str(np.round(dye_cell_final,6)) +
                             ' mmol/L')

        plt.close()
        loggers.log_info('Simulation completed successfully.')

    def update_V_ecm(self,cells,p,t):


        if p.sim_ECM == True:

            self.rho_cells = get_charge_density(self.cc_cells, self.z_array, p)
            self.rho_env = get_charge_density(self.cc_env, self.z_array_env, p)
            self.v_env = get_Venv(self,cells,p)
            self.v_cell = get_Vcell(self,cells,p)

            self.vm = self.v_cell[cells.mem_to_cells] - self.v_env[cells.map_mem2ecm]  # calculate v_mem

        else:

            self.rho_cells = get_charge_density(self.cc_cells, self.z_array, p)
            self.vm = get_Vcell(self,cells,p)

    def update_C_ecm(self,ion_i,flux,cells,p):

        c_cells = self.cc_cells[ion_i][:]
        c_env = self.cc_env[ion_i][:]

        d_c_cells = flux*(cells.mem_sa/cells.cell_vol[cells.mem_to_cells])
        d_c_env = -flux*(cells.mem_sa/cells.ecm_vol)

        delta_cells =  np.dot(d_c_cells, cells.cell_UpdateMatrix)
        delta_env = np.dot(d_c_env, cells.ecm_UpdateMatrix)

        self.cc_cells[ion_i] = rk4(c_cells,delta_cells,p)

        self.cc_env[ion_i] = rk4(c_env,delta_env,p)

    def Hplus_electrofuse_ecm(self,cells,p,t):

        # electrofuse the H+ ion between the cytoplasm and the ecms
        f_H1 = \
            electroflux(self.cc_env[self.iH][cells.map_mem2ecm],self.cc_cells[self.iH][cells.mem_to_cells],
                self.Dm_cells[self.iH],self.tm[cells.mem_to_cells],self.zs[self.iH],self.vm,self.T,p)

        self.fluxes_mem[self.iH] =  self.fluxes_mem[self.iH] + f_H1

        # Update the anion (bicarbonate) concentration instead of H+, assuming bicarb buffer holds:
        self.update_C_ecm(self.iM,-f_H1,cells,p)

        # Calculate the new pH and H+ concentration:
        self.pH_cell = 6.1 + np.log10(self.cc_cells[self.iM]/self.cHM_cells)
        self.cc_cells[self.iH] = 10**(-self.pH_cell)

        self.pH_env = 6.1 + np.log10(self.cc_env[self.iM]/self.cHM_env)
        self.cc_env[self.iH] = 10**(-self.pH_env)

        # recalculate the net, unbalanced charge and voltage in each cell:
        self.update_V_ecm(cells,p,t)

    def Hplus_HKATP_ecm(self,cells,p,t):

        # if HKATPase pump is desired, run the H-K-ATPase pump:
        f_H2, f_K2 = pumpHKATP(self.cc_cells[self.iH][cells.mem_to_cells],self.cc_env[self.iH][cells.map_mem2ecm],
            self.cc_cells[self.iK][cells.mem_to_cells],self.cc_env[self.iK][cells.map_mem2ecm],
            self.vm,self.T,p,self.HKATP_block)

        self.fluxes_mem[self.iH] =  self.fluxes_mem[self.iH] + f_H2
        self.fluxes_mem[self.iK] =  self.fluxes_mem[self.iK] + f_K2

        # calculate the update to K+ in the cell and ecm:
        self.update_C_ecm(self.iK,f_K2,cells,p)

        # Update the anion (bicarbonate) concentration instead of H+, assuming bicarb buffer holds:
        self.update_C_ecm(self.iM,-f_H2,cells,p)

        self.pH_cell = 6.1 + np.log10(self.cc_cells[self.iM]/self.cHM_cells)
        self.cc_cells[self.iH] = 10**(-self.pH_cell)

        self.pH_env = 6.1 + np.log10(self.cc_env[self.iM]/self.cHM_env)
        self.cc_env[self.iH] = 10**(-self.pH_env)

        # recalculate the net, unbalanced charge and voltage in each cell:
        self.update_V_ecm(cells,p,t)

    def Hplus_VATP_ecm(self,cells,p,t):

        # if HKATPase pump is desired, run the H-K-ATPase pump:
        f_H3 = pumpVATP(self.cc_cells[self.iH][cells.mem_to_cells],self.cc_env[self.iH][cells.map_mem2ecm],
            self.vm,self.T,p,self.VATP_block)

        self.fluxes_mem[self.iH] =  self.fluxes_mem[self.iH] + f_H3

        # Update the anion (bicarbonate) concentration instead of H+, assuming bicarb buffer holds:
        self.update_C_ecm(self.iM,-f_H3,cells,p)

        self.pH_cell = 6.1 + np.log10(self.cc_cells[self.iM]/self.cHM_cells)
        self.cc_cells[self.iH] = 10**(-self.pH_cell)

        self.pH_env = 6.1 + np.log10(self.cc_env[self.iM]/self.cHM_env)
        self.cc_env[self.iH] = 10**(-self.pH_env)

        # recalculate the net, unbalanced charge and voltage in each cell:
        self.update_V_ecm(cells,p,t)

    def update_gj(self,cells,p,t,i):

        # calculate voltage difference (gradient*len_gj) between gj-connected cells:
        if p.sim_ECM == True:

            self.vgj = self.v_cell[cells.nn_i][:,1]- self.v_cell[cells.nn_i][:,0]

        else:
            self.vgj = self.vm[cells.nn_i][:,1]- self.vm[cells.nn_i][:,0]


        if p.gj_flux_sensitive == True and cells.nnAveMatrix != None:

            if p.gj_respond_flow == True:

                # map intracellular flow velocity to gap junctions:
                ux_gj = (self.u_cells_x[cells.nn_i][:,0] + self.u_cells_x[cells.nn_i][:,1])/2
                uy_gj = (self.u_cells_y[cells.nn_i][:,0] + self.u_cells_y[cells.nn_i][:,1])/2

                # get the magnitude of flow at the gap junction:
                u_gj_s = np.sqrt(ux_gj**2 + uy_gj**2)

                # sum the flow vector lengths for individual cells
                # u_sum = np.dot(cells.gj2cellSum,u_gj_s)

                # average flow vector lengths for whole cell collection:
                u_sum = np.mean(u_gj_s)

            else: # respond to total current instead:

                 # calculate current across gap junctions in x direction:
                I_gj_x = np.zeros(len(cells.nn_i))

                for flux_array, zi in zip(self.fluxes_gj_x,self.zs):

                    I_i_x = flux_array*zi*p.F

                    I_gj_x = I_gj_x + I_i_x

                I_gj_y = np.zeros(len(cells.nn_i))

                for flux_array, zi in zip(self.fluxes_gj_y,self.zs):

                    I_i_y = flux_array*zi*p.F

                    I_gj_y = I_gj_y + I_i_y

                u_gj_s = np.sqrt(I_gj_x**2 + I_gj_y**2)

                # sum the flow vector lengths for individual cells
                # u_sum = np.dot(cells.gj2cellSum,u_gj_s)

                u_sum =np.mean(u_gj_s)


            if u_sum.any() == 0:

                pass

            else:
                # normalize cell vector length by the cell sum -- this is the input to each rho_gj growth
                # U_gj_cell_norm = u_gj_s/u_sum[cells.nn_i][:,0]
                U_gj_cell_norm = u_gj_s/u_sum

                delta_gj_rho = (U_gj_cell_norm)*(p.max_gj_enhancement - self.gj_rho) -\
                               (p.u_decay_rate*p.max_gj_enhancement)*self.gj_rho

                # average the delta (calculated on nn duplicates) between cells:
                delta_gj_rho = np.dot(cells.nnAveMatrix,delta_gj_rho)

                self.gj_rho = self.gj_rho + p.dt*delta_gj_rho*p.alpha_rho_gj


        else:

            self.gj_rho = np.zeros(len(cells.nn_i))


        if p.v_sensitive_gj == True:
            # determine the open state of gap junctions:
            self.gjopen = self.gj_rho + self.gj_block*((1.0 - tb.step(abs(self.vgj),p.gj_vthresh,p.gj_vgrad) + 0.1))

        else:
            self.gjopen = (1 + self.gj_rho)*self.gj_block


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

        # electroosmotic fluid velocity -- averaged at gap junctions:
        if p.base_eosmo == True:
            ux = (self.u_cells_x[cells.nn_i][:,0] + self.u_cells_x[cells.nn_i][:,1])/2
            uy = (self.u_cells_y[cells.nn_i][:,0] + self.u_cells_y[cells.nn_i][:,1])/2

        else:
            ux = 0
            uy =0

        fgj_x,fgj_y = nernst_planck_flux(c,grad_cgj_x,grad_cgj_y,grad_vgj_x,grad_vgj_y,ux,uy,
            self.D_gj[i],self.zs[i],self.T,p)


        fgj = fgj_x*cells.nn_vects[:,2] + fgj_y*cells.nn_vects[:,3]

        delta_cc = np.dot(cells.gjMatrix*p.gj_surface*self.gjopen,fgj)

        self.cc_cells[i] = self.cc_cells[i] + p.dt*delta_cc

        self.fluxes_gj_x[i] = fgj_x  # store gap junction flux for this ion
        self.fluxes_gj_y[i] = fgj_y  # store gap junction flux for this ion

    def update_ecm(self,cells,p,t,i):

        if p.closed_bound == True:
            btag = 'closed'

        else:
            btag = 'open'
        # make v_env and cc_env into 2d matrices
        cenv = self.cc_env[i][:]
        # denv = self.D_env[i][:]

        v_env = self.v_env[:].reshape(cells.X.shape)

        # enforce zero voltage boundary:
        v_env[:,0] = self.bound_V['L']
        v_env[:,-1] = self.bound_V['R']
        v_env[0,:] = self.bound_V['B']
        v_env[-1,:] = self.bound_V['T']

        cenv = cenv.reshape(cells.X.shape)

        # prepare concentrations and diffusion constants for MACs grid format
        # by resampling the values at the u v coordinates of the flux:
        cenv_x = np.zeros(cells.grid_obj.u_shape)
        cenv_y = np.zeros(cells.grid_obj.v_shape)

        # create the proper shape for the concentrations and state appropriate boundary conditions::
        cenv_x[:,1:] = cenv

        cenv_y[1:,:] = cenv

        if p.closed_bound == True: # insulation boundary conditions
            cenv_x[:,0] = cenv_x[:,1]
            cenv_x[:,-1] = cenv_x[:,-2]
            cenv_x[0,:] = cenv_x[1,:]
            cenv_x[-1,:] = cenv_x[-2,:]

            cenv_y[0,:] = cenv_y[1,:]
            cenv_y[-1,:] = cenv_y[-2,:]
            cenv_y[:,0] = cenv_y[:,1]
            cenv_y[:,-1] = cenv_y[:,-2]


        else:   # open and electrically grounded boundary conditions
            cenv_x[:,0] =  self.c_env_bound[i]
            cenv_x[:,-1] =  self.c_env_bound[i]
            cenv_x[0,:] =  self.c_env_bound[i]
            cenv_x[-1,:] =  self.c_env_bound[i]

            cenv_y[0,:] =  self.c_env_bound[i]
            cenv_y[-1,:] =  self.c_env_bound[i]
            cenv_y[:,0] =  self.c_env_bound[i]
            cenv_y[:,-1] =  self.c_env_bound[i]

        # calculate gradients in the environment
        grad_V_env_x, grad_V_env_y = cells.grid_obj.grid_gradient(v_env,bounds='closed')

        grad_cc_env_x, grad_cc_env_y = cells.grid_obj.grid_gradient(cenv,bounds=btag)

        # calculate fluxes for electrodiffusive transport:
        if p.base_eosmo == True:
            uenvx = self.u_env_x[:]
            uenvy = self.u_env_y[:]

        else:
            uenvx = 0
            uenvy = 0

        f_env_x, f_env_y = np_flux_special(cenv_x,cenv_y,grad_cc_env_x,grad_cc_env_y,
            grad_V_env_x, grad_V_env_y, uenvx,uenvy,self.D_env_u[i],self.D_env_v[i],
            self.zs[i],self.T,p)

        if p.closed_bound == False:
            f_env_x[:,0] = f_env_x[:,1]
            f_env_x[:,-1]= f_env_x[:,-2]
            f_env_x[0,:] = f_env_x[1,:]
            f_env_x[-1,:] = f_env_x[-2,:]

            f_env_y[:,0] = f_env_y[:,1]
            f_env_y[:,-1]= f_env_y[:,-2]
            f_env_y[0,:] = f_env_y[1,:]
            f_env_y[-1,:] = f_env_y[-2,:]

        else:

            f_env_x[:,0] = 0
            f_env_x[:,-1]= 0
            f_env_x[0,:] = 0
            f_env_x[-1,:] = 0

            f_env_y[:,0] = 0
            f_env_y[:,-1]= 0
            f_env_y[0,:] = 0
            f_env_y[-1,:] = 0

        # # # calculate the divergence of the total flux, which is equivalent to the total change per unit time
        # delta_c = fd.flux_summer(f_env_x,f_env_y,cells.X)*(1/cells.delta)

        d_fenvx = -(f_env_x[:,1:] - f_env_x[:,0:-1])/cells.delta
        d_fenvy = -(f_env_y[1:,:] - f_env_y[0:-1,:])/cells.delta

        delta_c = d_fenvx + d_fenvy

        #-----------------------
        cenv = cenv + delta_c*p.dt

        if p.closed_bound == True:
            # Neumann boundary condition (flux at boundary)
            # zero flux boundaries for concentration:
            cenv[:,-1] = cenv[:,-2]
            cenv[:,0] = cenv[:,1]
            cenv[0,:] = cenv[1,:]
            cenv[-1,:] = cenv[-2,:]

        elif p.closed_bound == False:
            # if the boundary is open, set the concentration at the boundary
            # open boundary
            cenv[:,-1] = self.c_env_bound[i]
            cenv[:,0] = self.c_env_bound[i]
            cenv[0,:] = self.c_env_bound[i]
            cenv[-1,:] = self.c_env_bound[i]

        # reshape the matrices back into vectors:
        self.cc_env[i] = cenv.ravel()

        fenvx = (f_env_x[:,1:] + f_env_x[:,0:-1])/2
        fenvy = (f_env_y[1:,:] + f_env_y[0:-1,:])/2

        self.fluxes_env_x[i] = fenvx.ravel()  # store ecm junction flux for this ion
        self.fluxes_env_y[i] = fenvy.ravel()  # store ecm junction flux for this ion

    def update_er(self,cells,p,t):

         # electrodiffusion of ions between cell and endoplasmic reticulum
        f_Ca_ER = \
        electroflux(self.cc_cells[self.iCa],self.cc_er[0],self.Dm_er[0],self.tm,self.z_er[0],self.v_er,self.T,p)

        # Electrodiffusion of charge compensation anion
        f_M_ER = \
        electroflux(self.cc_cells[self.iM],self.cc_er[1],self.Dm_er[1],self.tm,self.z_er[1],self.v_er,self.T,p)

        # update calcium concentrations in the ER and cell:
        self.cc_er[0] = self.cc_er[0] + f_Ca_ER*((cells.cell_sa)/(p.ER_vol*cells.cell_vol))*p.dt
        self.cc_cells[self.iCa] = self.cc_cells[self.iCa] - f_Ca_ER*(cells.cell_sa/cells.cell_vol)*p.dt

        # update anion concentration in the ER and cell:
        self.cc_er[1] = self.cc_er[1] + f_M_ER*((cells.cell_sa)/(p.ER_vol*cells.cell_vol))*p.dt
        self.cc_cells[self.iM] = self.cc_cells[self.iM] - f_M_ER*(cells.cell_sa/cells.cell_vol)*p.dt

        # recalculate the net, unbalanced charge and voltage in each cell:
        self.update_V_ecm(cells,p,t)

        q_er = get_charge(self.cc_er,self.z_array_er,p.ER_vol*cells.cell_vol,p)
        self.v_er = get_volt(q_er,p.ER_sa*cells.cell_sa,p) - self.v_cell

    def update_dye(self,cells,p,t):

        # Update dye concentration in the gj connected cell network:

        # voltage gradient:
        grad_vgj = self.vgj/cells.nn_len

        grad_vgj_x = grad_vgj*cells.nn_vects[:,2]
        grad_vgj_y = grad_vgj*cells.nn_vects[:,3]

        # concentration gradient for Dye:
        grad_cgj = (self.cDye_cell[cells.nn_i][:,1] - self.cDye_cell[cells.nn_i][:,0])/cells.nn_len

        grad_cgj_x = grad_cgj*cells.nn_vects[:,2]
        grad_cgj_y = grad_cgj*cells.nn_vects[:,3]

        # midpoint concentration:
        cdye = (self.cDye_cell[cells.nn_i][:,1] + self.cDye_cell[cells.nn_i][:,0])/2

        # electroosmotic fluid velocity:
        if p.base_eosmo == True:
            ux = (self.u_cells_x[cells.nn_i][:,0] + self.u_cells_x[cells.nn_i][:,1])/2
            uy = (self.u_cells_y[cells.nn_i][:,0] + self.u_cells_y[cells.nn_i][:,1])/2

        else:
            ux = 0
            uy = 0

        fgj_x_dye,fgj_y_dye = nernst_planck_flux(cdye,grad_cgj_x,grad_cgj_y,grad_vgj_x,grad_vgj_y,ux,uy,
            p.Do_Dye*self.gjopen,p.z_Dye,self.T,p)

        fgj_dye = fgj_x_dye*cells.nn_vects[:,2] + fgj_y_dye*cells.nn_vects[:,3]

        delta_cc = np.dot(cells.gjMatrix*p.gj_surface*self.gjopen,fgj_dye)

        self.cDye_cell = self.cDye_cell + p.dt*delta_cc

        self.Dye_flux_x_gj = fgj_x_dye[:]  # store gap junction flux for this ion
        self.Dye_flux_y_gj = fgj_y_dye[:]  # store gap junction flux for this ion

        if p.sim_ECM == False:

            fdye_ED = electroflux(self.cDye_env,self.cDye_cell,self.id_cells*p.Dm_Dye,self.tm,p.z_Dye,self.vm,self.T,p)

            # update dye concentration
            self.cDye_cell = self.cDye_cell + fdye_ED*(cells.cell_sa/cells.cell_vol)*p.dt

        elif p.sim_ECM == True:

            flux_dye = electroflux(self.cDye_env[cells.map_mem2ecm],self.cDye_cell[cells.mem_to_cells],
                            np.ones(len(cells.mem_i))*p.Dm_Dye,self.tm,p.z_Dye,self.vm,self.T,p)

            # update the dye concentrations in the cell and ecm due to ED fluxes at membrane
            d_c_cells = flux_dye*(cells.mem_sa/cells.cell_vol[cells.mem_to_cells])
            d_c_env = flux_dye*(cells.mem_sa/cells.ecm_vol)

            delta_cells =  np.dot(d_c_cells, cells.cell_UpdateMatrix)
            delta_env = np.dot(d_c_env, cells.ecm_UpdateMatrix)

            self.cDye_cell = rk4(self.cDye_cell,delta_cells,p)

            self.cDye_env = rk4(self.cDye_env,-delta_env,p)

            # transport dye through environment: _________________________________________________________
            if p.closed_bound == True:
                btag = 'closed'

            else:
                btag = 'open'
             # make v_env and cc_env into 2d matrices
            cenv = self.cDye_env[:]
            denv = p.Do_Dye*np.ones(len(cells.xypts))

            v_env = self.v_env.reshape(cells.X.shape)

            v_env[:,0] = self.bound_V['L']
            v_env[:,-1] = self.bound_V['R']
            v_env[0,:] = self.bound_V['B']
            v_env[-1,:] = self.bound_V['T']

            cenv = cenv.reshape(cells.X.shape)

            # prepare concentrations and diffusion constants for MACs grid format
            # by resampling the values at the u v coordinates of the flux:
            cenv_x = np.zeros(cells.grid_obj.u_shape)
            cenv_y = np.zeros(cells.grid_obj.v_shape)

            # create the proper shape for the concentrations and state appropriate boundary conditions::
            cenv_x[:,1:] = cenv[:]
            cenv_x[:,0] = cenv_x[:,1]
            cenv_y[1:,:] = cenv[:]
            cenv_y[0,:] = cenv_y[1,:]

            if p.closed_bound == True: # insulation boundary conditions
                cenv_x[:,0] = cenv_x[:,1]
                cenv_x[:,-1] = cenv_x[:,-2]
                cenv_x[0,:] = cenv_x[1,:]
                cenv_x[-1,:] = cenv_x[-2,:]

                cenv_y[0,:] = cenv_y[1,:]
                cenv_y[-1,:] = cenv_y[-2,:]
                cenv_y[:,0] = cenv_y[:,1]
                cenv_y[:,-1] = cenv_y[:,-2]

            else:   # open and electrically grounded boundary conditions
                cenv_x[:,0] =  self.c_dye_bound
                cenv_x[:,-1] =   self.c_dye_bound
                cenv_x[0,:] =   self.c_dye_bound
                cenv_x[-1,:] =   self.c_dye_bound

                cenv_y[0,:] =   self.c_dye_bound
                cenv_y[-1,:] =   self.c_dye_bound
                cenv_y[:,0] =   self.c_dye_bound
                cenv_y[:,-1] =   self.c_dye_bound

            denv = denv.reshape(cells.X.shape)

            denv_x = interp.griddata((cells.xypts[:,0],cells.xypts[:,1]),denv.ravel(),
                    (cells.grid_obj.u_X,cells.grid_obj.u_Y),method='nearest',fill_value = p.Do_Dye)

            denv_y = interp.griddata((cells.xypts[:,0],cells.xypts[:,1]),denv.ravel(),
                    (cells.grid_obj.v_X,cells.grid_obj.v_Y),method='nearest',fill_value=p.Do_Dye)

            denv_x = denv_x*self.D_env_weight_u
            denv_y = denv_y*self.D_env_weight_v

            # calculate gradients in the environment
            grad_V_env_x, grad_V_env_y = cells.grid_obj.grid_gradient(v_env,bounds='closed')

            grad_cc_env_x, grad_cc_env_y = cells.grid_obj.grid_gradient(cenv,bounds=btag)

            # calculate fluxes for electrodiffusive transport:

            if p.base_eosmo == True:
                uenvx = self.u_env_x[:]
                uenvy = self.u_env_y[:]

            else:
                uenvx = 0
                uenvy = 0

            f_env_x_dye, f_env_y_dye = np_flux_special(cenv_x,cenv_y,grad_cc_env_x,grad_cc_env_y,
                grad_V_env_x, grad_V_env_y, uenvx,uenvy,denv_x,denv_y,p.z_Dye,self.T,p)

            # calculate the divergence of the total flux, which is equivalent to the total change per unit time:
            delta_c = fd.flux_summer(f_env_x_dye,f_env_y_dye,cells.X)/cells.delta

            cenv = cenv + delta_c*p.dt

            if p.closed_bound == True:
                # Neumann boundary condition (flux at boundary)
                # zero flux boundaries for concentration:
                cenv[:,-1] = cenv[:,-2]
                cenv[:,0] = cenv[:,1]
                cenv[0,:] = cenv[1,:]
                cenv[-1,:] = cenv[-2,:]

            elif p.closed_bound == False:
                # if the boundary is open, set the concentration at the boundary
                # open boundary
                cenv[:,-1] =  self.c_dye_bound
                cenv[:,0] =  self.c_dye_bound
                cenv[0,:] =  self.c_dye_bound
                cenv[-1,:] =  self.c_dye_bound

            # reshape the matrices into vectors:
            # self.v_env = self.v_env.ravel()
            self.cDye_env = cenv.ravel()

            fenvx = (f_env_x_dye[:,1:] + f_env_x_dye[:,0:-1])/2
            fenvy = (f_env_y_dye[1:,:] + f_env_y_dye[0:-1,:])/2

            self.Dye_flux_env_x = fenvx.ravel()  # store ecm junction flux for this ion
            self.Dye_flux_env_y = fenvy.ravel()  # store ecm junction flux for this ion

    def update_IP3(self,cells,p,t):

        # Update dye concentration in the gj connected cell network:
        # voltage gradient:
        grad_vgj = self.vgj/cells.nn_len

        grad_vgj_x = grad_vgj*cells.nn_vects[:,2]
        grad_vgj_y = grad_vgj*cells.nn_vects[:,3]

        # concentration gradient for Dye:
        grad_cgj = (self.cIP3[cells.nn_i][:,1] - self.cIP3[cells.nn_i][:,0])/cells.nn_len

        grad_cgj_x = grad_cgj*cells.nn_vects[:,2]
        grad_cgj_y = grad_cgj*cells.nn_vects[:,3]

        # midpoint concentration:
        cip3 = (self.cIP3[cells.nn_i][:,1] + self.cIP3[cells.nn_i][:,0])/2

        # electroosmotic fluid velocity:
        if p.base_eosmo == True:
            ux = (self.u_cells_x[cells.nn_i][:,0] + self.u_cells_x[cells.nn_i][:,1])/2
            uy = (self.u_cells_y[cells.nn_i][:,0] + self.u_cells_y[cells.nn_i][:,1])/2

        else:
            ux = 0
            uy = 0

        fgj_x_ip3,fgj_y_ip3 = nernst_planck_flux(cip3,grad_cgj_x,grad_cgj_y,grad_vgj_x,grad_vgj_y,ux,uy,
            p.Do_IP3*self.gjopen,p.z_IP3,self.T,p)

        fgj_ip3 = fgj_x_ip3*cells.nn_vects[:,2] + fgj_y_ip3*cells.nn_vects[:,3]

        delta_cc = np.dot(cells.gjMatrix*p.gj_surface*self.gjopen,fgj_ip3)

        self.cIP3 = self.cIP3 + p.dt*delta_cc

        self.IP3_flux_x_gj = fgj_x_ip3[:]  # store gap junction flux for this ion
        self.IP3_flux_y_gj = fgj_y_ip3[:]  # store gap junction flux for this ion

        if p.sim_ECM == False:

            fip3_ED = electroflux(self.cIP3_env,self.cIP3,self.id_cells*p.Dm_IP3,self.tm,p.z_IP3,self.vm,self.T,p)

            # update dye concentration
            self.cIP3 = self.cIP3 + fip3_ED*(cells.cell_sa/cells.cell_vol)*p.dt

        elif p.sim_ECM == True:

            flux_ip3 = electroflux(self.cIP3_env[cells.map_mem2ecm],self.cIP3[cells.mem_to_cells],
                            np.ones(len(cells.mem_i))*p.Dm_IP3,self.tm,p.z_IP3,self.vm,self.T,p)

            # update the dye concentrations in the cell and ecm due to ED fluxes at membrane
            d_c_cells = flux_ip3*(cells.mem_sa/cells.cell_vol[cells.mem_to_cells])
            d_c_env = flux_ip3*(cells.mem_sa/cells.ecm_vol)

            delta_cells =  np.dot(d_c_cells, cells.cell_UpdateMatrix)
            delta_env = np.dot(d_c_env, cells.ecm_UpdateMatrix)

            self.cIP3 = rk4(self.cIP3,delta_cells,p)

            self.cIP3_env = rk4(self.cIP3_env,-delta_env,p)

            # transport dye through environment: _________________________________________________________
            if p.closed_bound == True:
                btag = 'closed'

            else:
                btag = 'open'
             # make v_env and cc_env into 2d matrices
            cenv = self.cIP3_env[:]
            denv = p.Do_IP3*np.ones(len(cells.xypts))

            v_env = self.v_env.reshape(cells.X.shape)

            v_env[:,0] = self.bound_V['L']
            v_env[:,-1] = self.bound_V['R']
            v_env[0,:] = self.bound_V['B']
            v_env[-1,:] = self.bound_V['T']

            cenv = cenv.reshape(cells.X.shape)

            # prepare concentrations and diffusion constants for MACs grid format
            # by resampling the values at the u v coordinates of the flux:
            cenv_x = np.zeros(cells.grid_obj.u_shape)
            cenv_y = np.zeros(cells.grid_obj.v_shape)

            # create the proper shape for the concentrations and state appropriate boundary conditions::
            cenv_x[:,1:] = cenv[:]
            cenv_x[:,0] = cenv_x[:,1]
            cenv_y[1:,:] = cenv[:]
            cenv_y[0,:] = cenv_y[1,:]

            if p.closed_bound == True: # insulation boundary conditions
                cenv_x[:,0] = cenv_x[:,1]
                cenv_x[:,-1] = cenv_x[:,-2]
                cenv_x[0,:] = cenv_x[1,:]
                cenv_x[-1,:] = cenv_x[-2,:]

                cenv_y[0,:] = cenv_y[1,:]
                cenv_y[-1,:] = cenv_y[-2,:]
                cenv_y[:,0] = cenv_y[:,1]
                cenv_y[:,-1] = cenv_y[:,-2]

            else:   # open and electrically grounded boundary conditions
                cenv_x[:,0] =  p.cIP3_to_env
                cenv_x[:,-1] =  p.cIP3_to_env
                cenv_x[0,:] =  p.cIP3_to_env
                cenv_x[-1,:] =  p.cIP3_to_env

                cenv_y[0,:] =  p.cIP3_to_env
                cenv_y[-1,:] =  p.cIP3_to_env
                cenv_y[:,0] =  p.cIP3_to_env
                cenv_y[:,-1] =  p.cIP3_to_env

            denv = denv.reshape(cells.X.shape)

            denv_x = np.zeros(cells.grid_obj.u_shape)
            denv_y = np.zeros(cells.grid_obj.v_shape)

            # create the proper shape for the diffusion constants and state continuous boundaries:
            denv_x[:,1:] = denv
            denv_x[:,0] = denv_x[:,1]

            denv_y[1:,:] = denv
            denv_y[0,:] = denv_y[1,:]

            # calculate gradients in the environment
            grad_V_env_x, grad_V_env_y = cells.grid_obj.grid_gradient(v_env,bounds='closed')

            grad_cc_env_x, grad_cc_env_y = cells.grid_obj.grid_gradient(cenv,bounds=btag)

            # calculate fluxes for electrodiffusive transport:

            if p.base_eosmo == True:
                uenvx = self.u_env_x[:]
                uenvy = self.u_env_y[:]

            else:
                uenvx = 0
                uenvy = 0

            f_env_x_ip3, f_env_y_ip3 = np_flux_special(cenv_x,cenv_y,grad_cc_env_x,grad_cc_env_y,
                grad_V_env_x, grad_V_env_y, uenvx,uenvy,denv_x,denv_y,p.z_IP3,self.T,p)

            # calculate the divergence of the total flux, which is equivalent to the total change per unit time:
            delta_c = fd.flux_summer(f_env_x_ip3,f_env_y_ip3,cells.X)/cells.delta

            cenv = cenv + delta_c*p.dt

            if p.closed_bound == True:
                # Neumann boundary condition (flux at boundary)
                # zero flux boundaries for concentration:
                cenv[:,-1] = cenv[:,-2]
                cenv[:,0] = cenv[:,1]
                cenv[0,:] = cenv[1,:]
                cenv[-1,:] = cenv[-2,:]

            elif p.closed_bound == False:
                # if the boundary is open, set the concentration at the boundary
                # open boundary
                cenv[:,-1] = p.cIP3_to_env
                cenv[:,0] = p.cIP3_to_env
                cenv[0,:] = p.cIP3_to_env
                cenv[-1,:] = p.cIP3_to_env

            # reshape the matrices into vectors:
            # self.v_env = self.v_env.ravel()
            self.cIP3_env = cenv.ravel()

            fenvx = (f_env_x_ip3[:,1:] + f_env_x_ip3[:,0:-1])/2
            fenvy = (f_env_y_ip3[1:,:] + f_env_y_ip3[0:-1,:])/2

            self.IP3_flux_env_x = fenvx.ravel()  # store ecm junction flux for this ion
            self.IP3_flux_env_y = fenvy.ravel()  # store ecm junction flux for this ion

    def get_Efield(self,cells,p):

         # calculate voltage difference (gradient*len_gj) between gj-connected cells:
        if p.sim_ECM == True:

            self.Egj = - (self.v_cell[cells.nn_i][:,1]- self.v_cell[cells.nn_i][:,0])/cells.nn_len

            # in the environment:
            venv = self.v_env.reshape(cells.X.shape)
            genv_x, genv_y = fd.gradient(venv, cells.delta)

            self.E_env_x = -genv_x
            self.E_env_y = -genv_y

        else:
            self.Egj = - (self.vm[cells.nn_i][:,1]- self.vm[cells.nn_i][:,0])/cells.nn_len

        # get x and y components of the electric field:
        self.E_gj_x = cells.nn_vects[:,2]*self.Egj
        self.E_gj_y = cells.nn_vects[:,3]*self.Egj

    def get_Bfield(self,cells,p):

         # calculate current across gap junctions in x an y directions:
        I_gj_x = np.zeros(len(cells.nn_i))

        for flux_array, zi in zip(self.fluxes_gj_x,self.zs):

            I_i_x = flux_array*zi*p.F

            I_gj_x = I_gj_x + I_i_x

        I_gj_y = np.zeros(len(cells.nn_i))

        for flux_array, zi in zip(self.fluxes_gj_y,self.zs):

            I_i_y = flux_array*zi*p.F

            I_gj_y = I_gj_y + I_i_y

        sourceBx = -p.mu*np.dot(cells.gj2cellMatrix,I_gj_x)
        sourceBy = -p.mu*np.dot(cells.gj2cellMatrix,I_gj_y)

        self.Bx = np.dot(cells.lapGJinv,sourceBx)
        self.By = np.dot(cells.lapGJinv,sourceBy)

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
        self.I_gj_x = interp.griddata((cells.nn_vects[:,0],cells.nn_vects[:,1]),I_gj_x,(cells.Xgrid,cells.Ygrid),
                                      method='nearest',fill_value=0)
        self.I_gj_x = np.multiply(self.I_gj_x,cells.maskM)

        self.I_gj_y = interp.griddata((cells.nn_vects[:,0],cells.nn_vects[:,1]),I_gj_y,(cells.Xgrid,cells.Ygrid),
                                      method='nearest',fill_value=0)
        self.I_gj_y = np.multiply(self.I_gj_y,cells.maskM)

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

    def getFlow(self,cells,p):
        """
        Calculate the electroosmotic-magneto-hydrodynamic fluid flow in the cell and extracellular
         networks.

        """

        if p.sim_ECM== True:

            # method 1-------------------------------------------------------------------------------------------------
            # force of gravity:
            if p.closed_bound == True:

                btag = 'closed'

            else:

                btag = 'open'

             # estimate the "pipe radius" for extracellular flow based on the diffusion constant weighting for the world:
            sa_env_x = cells.grid_obj.delta*p.cell_height*self.D_env_weight_u
            sa_env_y = cells.grid_obj.delta*p.cell_height*self.D_env_weight_v

            r_env_x = np.sqrt(sa_env_x/math.pi)          # average env radius, u grid
            r_env_y = np.sqrt(sa_env_y/math.pi)          # average env radius, v grid

            # fluid "conductivity" coefficient:
            alpha_x = (r_env_x**2)/(8*p.mu_water)
            alpha_y = (r_env_y**2)/(8*p.mu_water)

            if p.gravity == True:
                F_gravity_x = np.zeros(cells.grid_obj.u_shape)
                F_gravity_y = -np.ones(cells.grid_obj.v_shape)*9.81*1000

            else:
                F_gravity_x = np.zeros(cells.grid_obj.u_shape)
                F_gravity_y = np.zeros(cells.grid_obj.v_shape)

            # calculate the electroosmotic force (shaped to the u and v grids):
            # first calculate the electric field on the u and v junctions:
            venv = self.v_env.reshape(cells.X.shape)
            env_x, env_y = cells.grid_obj.grid_gradient(venv,bounds=btag)

            # create u and v sized grids for the charge density:
            rho_env_x = np.zeros(cells.grid_obj.u_shape)
            rho_env_y = np.zeros(cells.grid_obj.v_shape)

            # map the charge density to the grid
            rho_env_x[:,1:] = self.rho_env.reshape(cells.X.shape)/100
            rho_env_x[:,0] = rho_env_x[:,1]

            rho_env_y[1:,:] = self.rho_env.reshape(cells.X.shape)/100
            rho_env_y[0,:] = rho_env_y[1,:]

            # these are negative because the gradient of the voltage is the electric field and we just took the grad
            # above but didn't carry through the negative sign.
            Fe_x = -(rho_env_x)*env_x
            Fe_y = -(rho_env_y)*env_y

            # sum the forces:
            Fx = Fe_x + F_gravity_x
            Fy = Fe_y + F_gravity_y

            # # scale the forces by the fluid conductivity (which may change over the space):
            Sx = Fx*alpha_x
            Sy = Fy*alpha_y

            # calculate the divergence of the scaled conductivity as the sum of the two spatial derivatives:
            Sx_dx = (Sx[:,1:] - Sx[:,0:-1])/cells.grid_obj.delta
            Sy_dy = (Sy[1:,:] - Sy[0:-1,:])/cells.grid_obj.delta

            div_S = Sx_dx + Sy_dy

            source = div_S.ravel()

            # calculate the alpha-scaled pressure from the divergence of the force:
            if p.closed_bound == True:

                P = np.dot(cells.lapENV_P_inv, source)
                P = P.reshape(cells.grid_obj.cents_shape)

                # enforce zero normal gradient boundary conditions on P:
                P[:,0] = P[:,1]
                P[:,-1] = P[:,-2]
                P[0,:] = P[1,:]
                P[-1,:] = P[-2,:]

            else:
                P = np.dot(cells.lapENVinv, source)
                P = P.reshape(cells.grid_obj.cents_shape)

                # enforce zero pressure boundary conditions on P:
                P[:,0] = 0
                P[:,-1] = 0
                P[0,:] = 0
                P[-1,:] = 0

            # smooth out the pressure
            # P = fd.integrator(P)

            # Take the grid gradient of the scaled internal pressure:
            gPx, gPy = cells.grid_obj.grid_gradient(P,bounds=btag)

             # subtract the correction pressure term from the solution to yield a divergence-free flow field

            self.u_env_x = fd.integrator(alpha_x*Fx) - fd.integrator(gPx)
            self.u_env_y = fd.integrator(alpha_y*Fy) - fd.integrator(gPy)

            # reinforce boundary conditions
            if p.closed_bound == True:
                #left
                self.u_env_x[:,0] = 0
                # right
                self.u_env_x[:,-1] = 0
                # top
                self.u_env_x[-1,:] = 0
                # bottom
                self.u_env_x[0,:] = 0

                # left
                self.u_env_y[:,0] = 0
                # right
                self.u_env_y[:,-1] = 0
                # top
                self.u_env_y[-1,:] = 0
                # bottom
                self.u_env_y[0,:] = 0

            else:
                # left
                self.u_env_x[:,0] = self.u_env_x[:,1]
                # right
                self.u_env_x[:,-1] = self.u_env_x[:,-2]
                # top
                self.u_env_x[-1,:] = self.u_env_x[-2,:]
                # bottom
                self.u_env_x[0,:] = self.u_env_x[1,:]

                # left
                self.u_env_y[:,0] = self.u_env_y[:,1]
                # right
                self.u_env_y[:,-1] = self.u_env_y[:,-2]
                # top
                self.u_env_y[-1,:] = self.u_env_y[-2,:]
                # bottom
                self.u_env_y[0,:] = self.u_env_y[1,:]


            # interpolate u and v values at the centre for easy plotting:
            self.u_at_c = (self.u_env_x[:,0:-1]+self.u_env_x[:,1:])/2

            self.v_at_c = (self.u_env_y[0:-1,:]+self.u_env_y[1:,:])/2

            # resample alpha:
            alpha = (alpha_x[:,0:-1] + alpha_x[:,1:]+alpha_y[0:-1,:] + alpha_y[1:,:])/4

            alpha_zero = list(*(alpha.ravel() == 0).nonzero())

            # if len(alpha_zero) == 0:
            #
            #     self.P_env = P[:]*(1/alpha)

            self.P_env = P[:]


        #---------------Flow through gap junction connected cells-------------------------------------------------------

        # gravity force:

        if p.gravity == True:
            Fgx = np.zeros(len(cells.nn_i))
            Fgy = np.zeros(len(cells.nn_i))
            Fgy[:] = -9.81*1000

            Fgj_gravity = cells.nn_vects[:,2]*Fgx + cells.nn_vects[:,3]*Fgy

        else:
            Fgj_gravity = np.zeros(len(cells.nn_vects))

        sa_term = p.gj_surface*self.gjopen

        sagj = sa_term*cells.ave_sa_all    # average total gj surface area
        rgj = np.sqrt(sagj/math.pi)          # average gj radius

        alpha_gj = ((rgj**2)/(8*p.mu_water))

        # to get the body force at each gap junction, first map the charge density from cell to gj:
        rho_gj = (self.rho_cells[cells.nn_i][:,0] + self.rho_cells[cells.nn_i][:,1])/2

        # body force is equal to the electric field at the gap junction multiplied by the charge density there.
        F_gj = rho_gj*self.Egj

        F_source = F_gj + Fgj_gravity

        # sum the tangential body force pressure at the gap junctions for each cell:
        # scale the forces to the alpha value as 'rgj' may vary over space:
        S_source = alpha_gj*F_source

        Fgj_sum = np.dot(cells.gjMatrix,S_source)

        # # calculate the pressure in each cell required to create a divergence-free (mass conserved) flow field
        # NOTE! this represents the fluid conductivity term alpha_gj*grad_pressure
        # We needed to to it this way because alpha_gj varries over space.
        # Note also that the 'saterm' (variable GJ surface area) doesn't show up in Fgj_sum because it shows up on the
        # opposite side of the GJinv equation -- therefore it cancels out.
        self.P_cells = np.dot(cells.lapGJinv, -Fgj_sum)

        gradP = (self.P_cells[cells.nn_i][:,1] - self.P_cells[cells.nn_i][:,0])/cells.nn_len

        u_cells = (alpha_gj)*(F_gj +Fgj_gravity) - gradP

        u_cells_x = u_cells*cells.nn_vects[:,2]
        u_cells_y = u_cells*cells.nn_vects[:,3]

        # average components to the cell centres:
        self.u_cells_x = np.dot(cells.gj2cellMatrix,u_cells_x)
        self.u_cells_y = np.dot(cells.gj2cellMatrix,u_cells_y)

        # resample alpha to the cell centres:
        alpha_ave = np.dot(cells.gj2cellMatrix,alpha_gj)

        alpha_zero = list(*(alpha_ave == 0).nonzero())

        alpha_ave[alpha_zero] = 1

        self.P_cells = self.P_cells*(1/alpha_ave)

    def eosmosis(self,cells,p):

        """
        Electroosmosis of ion pumps and channels to potentially create directional fluxes in individual cells.

        This is presently simulated by calculating the Nernst-Planck concentration flux of a weighting
        agent, rho, which moves under its own concentration gradient and
        through the influence of the extracellular voltage gradient and fluid flows tangential to the membrane.

        """

        vmem = self.v_env[cells.map_mem2ecm]

        # components of fluid flow velocity at the membrane:
        if p.base_eosmo == True:
            ux_mem = self.u_at_c.ravel()[cells.map_mem2ecm]
            uy_mem = self.v_at_c.ravel()[cells.map_mem2ecm]

        else:
            ux_mem = 0
            uy_mem = 0

        # get the gradient of rho concentration around each membrane:
        grad_c = np.dot(cells.gradMem,self.rho_channel)

        # gcx = grad_c*cells.mem_vects_flat[:,4]   # FIXME no concentration gradient will be used for now
        # gcy = grad_c*cells.mem_vects_flat[:,5]
        gcx = 0
        gcy = 0

        # electric field at each membrane

        gvx = - self.E_env_x.ravel()[cells.map_mem2ecm]
        gvy = - self.E_env_y.ravel()[cells.map_mem2ecm]

        # calculate the total Nernst-Planck flux at each membrane:

        fx, fy = nernst_planck_flux(self.rho_channel, gcx, gcy, gvx, gvy,ux_mem,uy_mem,p.D_membrane,p.z_channel,
            self.T,p)

        # component of total flux in direction of membrane
        ftot = fx*cells.mem_vects_flat[:,4] + fy*cells.mem_vects_flat[:,5]

        # calculate the membrane-specific divergence of the flux which will be the "delta rho" at each membrane:
        delta_rho = np.dot(cells.gradMem,ftot)

        # update the channel property:
        self.rho_channel = self.rho_channel + p.dt*delta_rho

        # # make sure nothing is non-zero:
        fix_inds = (self.rho_channel < 0).nonzero()
        self.rho_channel[fix_inds] = 0

    def get_ion(self,label):
        """
        Given a string input, returns the simulation index of the appropriate ion.

        """

        if label == 'Na':

            ion = self.iNa

        elif label == 'K':

            ion = self.iK

        elif label == 'Ca':

            ion = self.iCa

        elif label == 'Cl':

            ion = self.iCl

        else:
            ion = []

        return ion

    def initDenv(self,cells,p):

        self.D_env_u = np.zeros((self.D_env.shape[0],cells.grid_obj.u_shape[0],cells.grid_obj.u_shape[1]))
        self.D_env_v = np.zeros((self.D_env.shape[0],cells.grid_obj.v_shape[0],cells.grid_obj.v_shape[1]))

        for i, dmat in enumerate(self.D_env):

            if p.env_type == False: # if air surrounds, first set everything to zero and add in cluster data...
                self.D_env[i][:] = 0
            # for all cells and mems in the cluster, set the internal diffusion constant for adherens junctions:
            dummyMems = np.ones(len(cells.mem_i))*self.D_free[i]*p.D_adh

            # get a list of all membranes for boundary cells:
            neigh_to_bcells,_,_ = tb.flatten(cells.cell_nn[cells.bflags_cells])

            all_bound_mem_inds = cells.cell_to_mems[cells.bflags_cells]
            interior_bound_mem_inds = cells.cell_to_mems[neigh_to_bcells]
            interior_bound_mem_inds,_,_ = tb.flatten(interior_bound_mem_inds)
            all_bound_mem_inds, _ ,_ = tb.flatten(all_bound_mem_inds)

            # set external membrane of boundary cells to the diffusion constant of tight junctions:
            dummyMems[all_bound_mem_inds] = self.D_free[i]*p.D_tj*self.Dtj_rel[i]
            dummyMems[interior_bound_mem_inds] = self.D_free[i]*p.D_tj*self.Dtj_rel[i]
            dummyMems[cells.bflags_mems] = self.D_free[i]     # FIXME note that this has been "opened up" since v0.2

            # interp the membrane data to an ecm grid, fill values correspond to environmental diffusion consts:
            if p.env_type == True:
                Denv_o = interp.griddata((cells.mem_vects_flat[:,0],cells.mem_vects_flat[:,1]),dummyMems,
                    (cells.X,cells.Y),method='nearest',fill_value=self.D_free[i])

            else:
                Denv_o = interp.griddata((cells.mem_vects_flat[:,0],cells.mem_vects_flat[:,1]),dummyMems,
                    (cells.X,cells.Y),method='nearest',fill_value=0)

            Denv_o = Denv_o.ravel()
            Denv_o[cells.inds_env] = self.D_free[i]

            # create an ecm diffusion grid filled with the environmental values
            self.D_env[i] = Denv_o

        # create a matrix that weights the relative transport efficiency in the world space:
        D_env_weight = self.D_env[self.iP]/self.D_env[self.iP].max()
        self.D_env_weight = D_env_weight.reshape(cells.X.shape)
        self.D_env_weight_base = np.copy(self.D_env_weight)

        for i, dmat in enumerate(self.D_env):

            if p.env_type == True:

                self.D_env_u[i] = interp.griddata((cells.xypts[:,0],cells.xypts[:,1]),dmat.ravel(),
                    (cells.grid_obj.u_X,cells.grid_obj.u_Y),method='nearest',fill_value = self.D_free[i])

                self.D_env_v[i] = interp.griddata((cells.xypts[:,0],cells.xypts[:,1]),dmat.ravel(),
                    (cells.grid_obj.v_X,cells.grid_obj.v_Y),method='nearest',fill_value=self.D_free[i])

            else:
                self.D_env_u[i] = interp.griddata((cells.xypts[:,0],cells.xypts[:,1]),dmat.ravel(),
                    (cells.grid_obj.u_X,cells.grid_obj.u_Y),method='nearest',fill_value = 0)

                self.D_env_v[i] = interp.griddata((cells.xypts[:,0],cells.xypts[:,1]),dmat.ravel(),
                    (cells.grid_obj.v_X,cells.grid_obj.v_Y),method='nearest',fill_value = 0)

        self.D_env_weight_u = self.D_env_u[self.iP]/self.D_env_u[self.iP].max()

        self.D_env_weight_v = self.D_env_v[self.iP]/self.D_env_v[self.iP].max()

        if p.closed_bound == True:  # set full no slip boundary condition at exterior bounds

            self.D_env_weight_u[:,0] = 0
            self.D_env_weight_u[:,-1] = 0
            self.D_env_weight_u[0,:] = 0
            self.D_env_weight_u[-1,:] = 0

            self.D_env_weight_v[:,0] = 0
            self.D_env_weight_v[:,-1] = 0
            self.D_env_weight_v[0,:] = 0
            self.D_env_weight_v[-1,:] = 0

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

    exp_alpha = np.exp(-alpha)

    deno = 1 - exp_alpha   # calculate the denominator for the electrodiffusion equation,..

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
                        cA[inotzero]*exp_alpha[inotzero])/deno[inotzero])

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

    alpha = block*p.alpha_NaK*tb.step(delG,p.halfmax_NaK,p.slope_NaK)

    f_Na  = -alpha*(cNai)*(cKo**(1/2))      #flux as [mol/m2s]   scaled to concentrations Na in and K out

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

    f_Ca  = alpha*(cCao**(1/2))*(1.0 - cCai)**(1/2)      #flux as [mol/s]

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
        f_H  = -alpha*np.sqrt(cHi)*np.sqrt(cKo)      #flux as [mol/s], scaled by concentrations in and out

    elif p.backward_pumps == True:

        alpha = signG*block*p.alpha_HK*tb.step(delG,p.halfmax_HK,p.slope_HK)

        truth_forwards = signG == 1
        truth_backwards = signG == -1

        inds_forwards = (truth_forwards).nonzero()  # indices of forward-running cells
        inds_backwards = (truth_backwards).nonzero() # indices of backward-running cells

        f_H = np.zeros(len(cHi))

        f_H[inds_forwards]  = -alpha*np.sqrt(cHi)*np.sqrt(cKo)      #flux as [mol/s], scaled by concentrations in and out
        f_H[inds_backwards]  = -alpha*np.sqrt(cHo)*np.sqrt(cKi)

    f_K = -f_H          # flux as [mol/s]

    return f_H, f_K

def pumpVATP(cHi,cHo,Vm,T,p,block):

    deltaGATP = 20*p.R*T

    delG_H = p.R*T*np.log(cHo/cHi) - p.F*Vm  # free energy to move H+ out of cell

    delG_VATP = deltaGATP - delG_H   # free energy available to move H+ out of cell
    delG_pump = (delG_VATP/1000)
    delG = np.absolute(delG_pump)
    signG = np.sign(delG)

    if p.backward_pumps == False:

        alpha = block*p.alpha_V*tb.step(delG,p.halfmax_V,p.slope_V)
        f_H  = -alpha*cHi      #flux as [mol/s], scaled by concentrations in and out

    elif p.backward_pumps == True:

        alpha = block*signG*p.alpha_V*tb.step(delG,p.halfmax_V,p.slope_V)

        truth_forwards = signG == 1
        truth_backwards = signG == -1

        inds_forwards = (truth_forwards).nonzero()  # indices of forward-running cells
        inds_backwards = (truth_backwards).nonzero() # indices of backward-running cells

        f_H = np.zeros(len(cHi))

        f_H[inds_forwards]  = -alpha*cHi      #flux as [mol/s], scaled by concentrations in and out
        f_H[inds_backwards]  = -alpha*cHo

    return f_H

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

    if p.sim_ECM == False:
        v_cell = (self.rho_cells*cells.cell_vol*p.tm)/(p.eo*80*cells.cell_sa)

    else:
        # get the value of the environmental voltage at each cell membrane:
        venv_at_mem = self.v_env[cells.map_mem2ecm]

        # sum the environmental voltage at each mem for each cell and take the average:
        cell_ave_Venv = np.dot(cells.M_sum_mems,venv_at_mem)/cells.num_mems

        # calculate the voltage in each cell:
        # v_cell = (self.rho_cells*cells.cell_vol*p.tm)/(p.eo*80*cells.cell_sa) + cell_ave_Venv
        v_cell = (self.rho_cells*cells.cell_vol*p.tm)/(p.eo*80*cells.cell_sa)

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

    # modify the source charge distribution in line with electrostatic Poisson equation:
    # note this should be divided by the electric permeability, but it produces way too high a voltage
    # in lieu of a feasible solution, the divisor is increased from 80*p.eo by self.ff
    self.ff = 1e6

    # smooth out the charge density:
    self.rho_env = self.rho_env.reshape(cells.X.shape)
    self.rho_env = fd.integrator(self.rho_env)
    self.rho_env = self.rho_env.ravel()

    fxy = -self.rho_env/(80*self.ff*p.eo)
    # fxy = -self.rho_env/(100*p.eo)

    # # modify the RHS of the equation to incorporate Dirichlet boundary conditions on Poisson voltage:
    fxy[cells.bBot_k] = (self.bound_V['B']/cells.delta**2)
    fxy[cells.bTop_k] = (self.bound_V['T']/cells.delta**2)
    fxy[cells.bL_k] = (self.bound_V['L']/cells.delta**2)
    fxy[cells.bR_k] = (self.bound_V['R']/cells.delta**2)

    # Solve Poisson's electrostatic equation:
    V = np.dot(cells.lapENVinv,fxy)

    # if the boundary conditions set the outside of the matrix:
    V[cells.bBot_k] = self.bound_V['B']
    V[cells.bTop_k] = self.bound_V['T']
    V[cells.bL_k] = self.bound_V['L']
    V[cells.bR_k] = self.bound_V['R']

    V = V.reshape(cells.X.shape)
    V = fd.integrator(V)
    V = V.ravel()

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

    dat_grid = interp.griddata((cells.plot_xy[:,0],cells.plot_xy[:,1]),plot_data,(cells.Xgrid,cells.Ygrid), method='nearest',
                               fill_value=0)
    #
    dat_grid = np.multiply(dat_grid,cells.maskM)

    return dat_grid

def nernst_planck_flux(c, gcx, gcy, gvx, gvy,ux,uy,D,z,T,p):

    alpha = (D*z*p.q)/(p.kb*T)

    fx =  -D*gcx - alpha*gvx*c + ux*c

    fy =  -D*gcy - alpha*gvy*c + uy*c

    return fx, fy

def np_flux_special(cx,cy,gcx,gcy,gvx,gvy,ux,uy,Dx,Dy,z,T,p):

    alphax = (Dx*z*p.q)/(p.kb*T)
    alphay = (Dy*z*p.q)/(p.kb*T)

    fx =  -Dx*gcx - alphax*gvx*cx + cx*ux

    fy =  -Dy*gcy - alphay*gvy*cy + cy*uy

    return fx, fy

def rk4(c,deltac,p):
    """
    This function is here for the day we figure out how to
    do a suitable RK4 solver for the complex physics we're
    simulating.

    For now it just computes forwards Euler method.
    """

    c2 = c + deltac*p.dt


    return c2



#-----------------------------------------------------------------------------------------------------------------------
# WASTELANDS
#-----------------------------------------------------------------------------------------------------------------------

             # method 2--------------------------------------------------------------------------------------------------
            # # force of gravity:
            # F_gravity_x = np.zeros(cells.grid_obj.u_X.shape)
            # F_gravity_y = -np.ones(cells.grid_obj.v_X.shape)*9.81*1000
            #
            # venv = self.v_env.reshape(cells.X.shape)
            # env_x, env_y = cells.grid_obj.grid_gradient(venv)
            #
            # rho_env_x = np.zeros(cells.grid_obj.u_X.shape)
            # rho_env_y = np.zeros(cells.grid_obj.v_X.shape)
            #
            # rho_env_x[:,0:-1] = self.rho_env.reshape(cells.X.shape)/(self.ff)
            # rho_env_y[0:-1,:] = self.rho_env.reshape(cells.X.shape)/(self.ff)
            #
            # Fe_x = (rho_env_x)*env_x
            # Fe_y = (rho_env_y)*env_y
            #
            # ts = p.dt*np.sqrt(cells.grid_obj.delta)
            # # ts = p.dt
            #
            # u = copy.copy(self.u_env_x)
            # v = copy.copy(self.u_env_y)
            #
            # # reinforce boundary conditions -- closed boundary
            #
            # if p.closed_bound == True:
            #     #left
            #     u[:,0] = 0
            #     # right
            #     u[:,-1] = 0
            #     # top
            #     u[-1,:] = 0
            #     # bottom
            #     u[0,:] = 0
            #
            #     # left
            #     v[:,0] = 0
            #     # right
            #     v[:,-1] = 0
            #     # top
            #     v[-1,:] = 0
            #     # bottom
            #     v[0,:] = 0
            #
            # else:
            #
            #     #left
            #     u[:,0] =  u[:,1]
            #     # right
            #     u[:,-1] = u[:,-2]
            #     # top
            #     u[-1,:] = u[-2,:]
            #     # bottom
            #     u[0,:] = u[1,:]
            #
            #     # left
            #     v[:,0] = v[:,1]
            #     # right
            #     v[:,-1] = v[:,-2]
            #     # top
            #     v[-1,:] = v[-2,:]
            #     # bottom
            #     v[0,:] = v[1,:]
            #
            # # scale the flow by the diffusion constant weighting matrix for the environment:
            # # u = self.D_env_weight_u*u
            # # v = self.D_env_weight_v*v
            #
            # # calculate the flow, omitting the pressure term:
            # lap_u = fd.laplacian(u,cells.grid_obj.delta)
            # lap_v = fd.laplacian(v,cells.grid_obj.delta)
            #
            # u = u + (ts/p.rho)*(p.mu_water*lap_u + fd.integrator(Fe_x))
            #
            # v = v + (ts/p.rho)*(p.mu_water*lap_v + fd.integrator(Fe_y))
            #
            # # take the divergence of the interm flow field using a forward difference
            # # that creates a matrix the same size as the pressure matrix:
            #
            # u_dx = (u[:,1:] - u[:,0:-1])/cells.grid_obj.delta
            # v_dy = (v[1:,:] - v[0:-1,:])/cells.grid_obj.delta
            #
            # div_u = u_dx + v_dy
            # div_u = fd.integrator(div_u)
            #
            # source = (p.rho/ts)*div_u.ravel()
            #
            # P = np.dot(cells.lapENV_P_inv, source)
            # P = P.reshape(cells.grid_obj.cents_shape)
            #
            # # enforce zero gradient boundary conditions on P:
            # P[:,0] = P[:,1]
            # P[:,-1] = P[:,-2]
            # P[0,:] = P[1,:]
            # P[-1,:] = P[-2,:]
            #
            # # Take the gradient of the pressue:
            # # gPx, gPy = cells.grid_obj.grid_gradient(self.P_env,bounds='open')
            # gPxo, gPyo = fd.gradient(P,cells.grid_obj.delta)
            #
            # gPx = np.zeros(cells.grid_obj.u_shape)
            # gPx[:,0:-1] = gPxo
            # gPx[:,-1] = gPxo[:,-1]
            #
            # gPy = np.zeros(cells.grid_obj.v_shape)
            # gPy[0:-1,:] = gPyo
            # gPy[-1,:] = gPyo[-1,:]
            #
            # # subtract the pressure from the solution to yeild a divergence-free flow field
            # self.u_env_x = u - fd.integrator(gPx)*(ts/p.rho)
            # self.u_env_y = v - fd.integrator(gPy)*(ts/p.rho)
            #
            # # scale the flow by the diffusion constant weighting matrix for the environment:
            # # self.u_env_x = self.D_env_weight_u*self.u_env_x
            # # self.u_env_y = self.D_env_weight_v*self.u_env_y
            #
            # # interpolate u and v values at the centre for easy plotting:
            # self.u_at_c = np.float64(self.u_env_x[:,0:-1])
            #
            # self.v_at_c = np.float64(self.u_env_y[0:-1,:])
            #
            # # reinforce boundary conditions
            #
            # if p.closed_bound == True:
            #     #left
            #     self.u_at_c[:,0] = 0
            #     # right
            #     self.u_at_c[:,-1] = 0
            #     # top
            #     self.u_at_c[-1,:] = 0
            #     # bottom
            #     self.u_at_c[0,:] = 0
            #
            #     # left
            #     self.v_at_c[:,0] = 0
            #     # right
            #     self.v_at_c[:,-1] = 0
            #     # top
            #     self.v_at_c[-1,:] = 0
            #     # bottom
            #     self.v_at_c[0,:] = 0
            #
            # else:
            #
            #     #left
            #     self.u_at_c[:,0] = self.u_at_c[:,1]
            #     # right
            #     self.u_at_c[:,-1] = self.u_at_c[:,-2]
            #     # top
            #     self.u_at_c[-1,:] = self.u_at_c[-2,:]
            #     # bottom
            #     self.u_at_c[0,:] = self.u_at_c[1,:]
            #
            #     # left
            #     self.v_at_c[:,0] = self.v_at_c[:,1]
            #     # right
            #     self.v_at_c[:,-1] = self.v_at_c[:,-2]
            #     # top
            #     self.v_at_c[-1,:] = self.v_at_c[-2,:]
            #     # bottom
            #     self.v_at_c[0,:] = self.v_at_c[1,:]
            #
            # self.P_env = np.float64(P[:])






