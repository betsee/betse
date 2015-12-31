#!/usr/bin/env python3
# Copyright 2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

#FIXME: Higher-level "betse.util.path" functions should typically be called
#rather than low-level "os.path" functions. Daggers of scintillating wonder!
import copy, math, os, time
import os.path
import matplotlib.pyplot as plt
import numpy as np
from betse.science import filehandling as fh
from betse.science import visualize as viz
from betse.science import toolbox as tb
from betse.science.tissue.handler import TissueHandler
from betse.science import finitediff as fd
from betse.exceptions import BetseExceptionSimulation
from betse.util.io import loggers
from random import shuffle
from scipy import interpolate as interp

class Simulator(object):
    """
    Contains the main routines used in the simulation of networked cell
    bioelectrical activity. For efficiency, all methods are implemented in
    terms of Numpy-based linear algebra.

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
        self.molar_mass = []


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

        self.P_cells = np.zeros(len(cells.cell_i))  # initialize pressure

        # initialize vectors for potential electroosmosis in the cell collection:
        if p.fluid_flow is True:
            self.u_cells_x = np.zeros(len(cells.cell_i))
            self.u_cells_y = np.zeros(len(cells.cell_i))

        if p.deformation is True:

            # initialize vectors for potential deformation:
            self.d_cells_x = np.zeros(len(cells.cell_i))
            self.d_cells_y = np.zeros(len(cells.cell_i))

        if p.gj_flux_sensitive is True:

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
                    self.molar_mass.append(p.molar_mass[name])
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
            self.molar_mass.append(p.M_H)
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

        if p.HKATPase_dyn is True and p.global_options['HKATP_block'] != 0:
            self.HKATP_block = np.ones(len(cells.cell_i))  # initialize HKATP blocking vector
        else:
            self.HKATP_block = 1

        if p.VATPase_dyn is True and p.global_options['VATP_block'] != 0:
            self.VATP_block = np.ones(len(cells.cell_i))  # initialize HKATP blocking vector
        else:
            self.VATP_block = 1

        # add channel noise to the model:
        self.channel_noise_factor = np.random.random(len(cells.cell_i))
        self.Dm_cells[self.iK] = (p.channel_noise_level*self.channel_noise_factor + 1)*self.Dm_cells[self.iK]

        if p.dynamic_noise is True:
            # add a random walk on protein concentration to generate dynamic noise:
            self.protein_noise_factor = p.dynamic_noise_level*(np.random.random(len(cells.cell_i)) - 0.5)

            if p.ions_dict['P']==1:
                self.cc_cells[self.iP] = self.cc_cells[self.iP]*(1+ self.protein_noise_factor)

        self.fluxes_mem = np.asarray(self.fluxes_mem)

        # Initialize Dye and IP3

        if p.scheduled_options['IP3'] != 0 or p.Ca_dyn is True:

            self.cIP3 = np.zeros(len(cells.cell_i))  # initialize a vector to hold IP3 concentrations
            self.cIP3[:] = p.cIP3_to                 # set the initial concentration of IP3 from params file

            self.cIP3_flux_gj = np.zeros(len(cells.nn_i))
            self.cIP3_flux_mem = np.zeros(len(cells.mem_i))

            self.cIP3_env = np.zeros(len(cells.cell_i))     # initialize IP3 concentration of the environment
            self.cIP3_env[:] = p.cIP3_to_env

        if p.voltage_dye is True:

            self.cDye_cell = np.zeros(len(cells.cell_i))   # initialize voltage sensitive dye array for cell and env't
            self.cDye_cell[:] = p.cDye_to_cell

            self.Dye_flux_gj = np.zeros(len(cells.nn_i))
            self.Dye_flux_mem = np.zeros(len(cells.mem_i))

            self.cDye_env = np.zeros(len(cells.cell_i))     # initialize Dye concentration in the environment
            self.cDye_env[:] = p.cDye_to
        #
        self.z_array = np.asarray(self.z_array)
        self.D_gj = np.asarray(self.D_gj)
        # self.molar_mass = np.asarray(self.molar_mass)

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

        self.molar_mass = []

        self.P_cells = np.zeros(len(cells.cell_i))  # initialize pressure in cells

        if p.fluid_flow is True:
            # Electroosmosis Initialization:
            # initialize vectors for env flow (note enhanced data type!):
            self.u_env_x = np.zeros(cells.X.shape)
            self.u_env_y = np.zeros(cells.X.shape)

            # initialize vectors for electroosmosis in the cell collection wrt each gap junction (note data type!):
            self.u_cells_x = np.zeros(len(cells.cell_i))
            self.u_cells_y = np.zeros(len(cells.cell_i))

        if p.deformation is True:
            # initialize vectors for potential deformation:
            self.d_cells_x = np.zeros(len(cells.cell_i))
            self.d_cells_y = np.zeros(len(cells.cell_i))


        if p.gj_flux_sensitive is True:

            self.gj_rho = np.zeros(len(cells.nn_i))

        else:

            self.gj_rho = 0

        if p.sim_eosmosis is True:
            self.rho_pump = np.ones(len(cells.mem_i))
            self.rho_channel = np.ones(len(cells.mem_i))
        else:
            self.rho_pump = 1  # else just define it as identity.
            self.rho_channel = 1

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
                    self.molar_mass.append(p.molar_mass[name])
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
            self.molar_mass.append(p.M_H)
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

        if p.HKATPase_dyn is True and p.global_options['HKATP_block'] != 0:
            self.HKATP_block = np.ones(len(cells.mem_i))  # initialize HKATP blocking vector
        else:
            self.HKATP_block = 1

        if p.VATPase_dyn is True and p.global_options['VATP_block'] != 0:
            self.VATP_block = np.ones(len(cells.mem_i))  # initialize HKATP blocking vector
        else:
            self.VATP_block = 1

        # add channel noise to the model:
        self.channel_noise_factor = np.random.random(len(cells.mem_i))
        self.Dm_cells[self.iK] = (p.channel_noise_level*self.channel_noise_factor + 1)*self.Dm_cells[self.iK]

        if p.dynamic_noise is True:
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
        # self.molar_mass = np.asarray(self.molar_mass)

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

        if p.sim_ECM is True:
            #  Initialize diffusion constants for the extracellular transport:
            self.initDenv(cells,p)

        self.dyna = TissueHandler(self,cells,p)   # create the tissue dynamics object
        self.dyna.tissueProfiles(self,cells,p)  # initialize all tissue profiles

        if p.sim_ECM is True:
            # create a copy-base of the environmental junctions diffusion constants:
            self.D_env_base = copy.copy(self.D_env)

        # add channel noise to the model:
        self.Dm_cells[self.iK] = (p.channel_noise_level*self.channel_noise_factor + 1)*self.Dm_cells[self.iK]

        # Initialize an array structure that will hold user-scheduled changes to membrane permeabilities:
        Dm_cellsA = np.asarray(self.Dm_cells)
        Dm_cellsER = np.copy(self.Dm_er)

        # if tb.emptyDict(p.scheduled_options) is False or tb.emptyDict(p.vg_options) is False or p.Ca_dyn is True:
        self.Dm_base = np.copy(Dm_cellsA) # make a copy that will serve as the unaffected values base

        # if tb.emptyDict(p.scheduled_options) is False:
        self.Dm_scheduled = np.copy(Dm_cellsA)
        self.Dm_scheduled[:] = 0

        # if tb.emptyDict(p.vg_options) is False:
            # Initialize an array structure that will hold dynamic voltage-gated channel changes to mem permeability:
        self.Dm_vg = np.copy(Dm_cellsA)
        self.Dm_vg[:] = 0

        # if p.Ca_dyn is True:
        # Initialize an array structure that will hold dynamic calcium-gated channel changes to mem perms:
        self.Dm_cag = np.copy(Dm_cellsA)
        self.Dm_cag[:] = 0

        self.Dm_stretch = np.copy(Dm_cellsA)   # array for stretch activated ion channels...

        self.Dm_er_base = np.copy(Dm_cellsER)

        self.Dm_er_CICR = np.copy(Dm_cellsER)
        self.Dm_er_CICR[:] = 0

        self.dcc_ER = []

        self.Dm_morpho = np.copy(Dm_cellsA)
        self.Dm_morpho[:] = 0

        self.P_mod = np.copy(self.P_cells[:])
        self.P_base = np.copy(self.P_cells[:])


        if p.global_options['gj_block'] != 0:

            self.gj_block = np.ones(len(cells.nn_i))   # initialize the gap junction blocking vector to ones

        else:

            self.gj_block = 1

        if p.scheduled_options['IP3'] != 0 or p.Ca_dyn is True:

            self.cIP3 = np.zeros(len(cells.cell_i))  # initialize a vector to hold IP3 concentrations
            self.cIP3[:] = p.cIP3_to                 # set the initial concentration of IP3 from params file

            self.IP3_flux_x_gj = np.zeros(len(cells.nn_i))
            self.IP3_flux_y_gj = np.zeros(len(cells.nn_i))
            self.IP3_flux_mem = np.zeros(len(cells.mem_i))



            if p.sim_ECM is True:
                # self.cIP3_ecm = np.zeros(len(cells.ecm_i))     # initialize IP3 concentration of the environment
                # self.cIP3_ecm[:] = p.cIP3_to_env
                self.cIP3_env = np.zeros(len(cells.xypts))
                self.cIP3_env[:] = p.cIP3_to_env

                self.cIP3_flux_env_x = np.zeros(len(cells.xypts))
                self.cIP3_flux_env_y = np.zeros(len(cells.xypts))

            elif p.sim_ECM is False:
                self.cIP3_env = np.zeros(len(cells.cell_i))     # initialize IP3 concentration of the environment
                self.cIP3_env[:] = p.cIP3_to_env

        if p.voltage_dye is True:

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

            if p.sim_ECM is True:

                self.cDye_env = np.zeros(len(cells.xypts))
                self.cDye_env[:] = p.cDye_to

                self.Dye_flux_env_x = np.zeros(len(cells.xypts))
                self.Dye_flux_env_y = np.zeros(len(cells.xypts))

            else:
                self.Dye_env = np.zeros(len(cells.cell_i))     # initialize Dye concentration in the environment
                self.Dye_env[:] = p.cDye_to

        # Initialize all user-specified interventions and dynamic channels.
        self.dyna.runAllInit(self,cells,p)

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

        self.dd_time = []  # data array holding membrane permeabilites at time points
        self.vm_time = []  # data array holding voltage at time points
        self.vm_GHK_time = [] # data array holding GHK vm estimates
        # self.dvm_time = []  # data array holding derivative of voltage at time points
        self.time = []     # time values of the simulation
        self.gjopen_time = []   # stores the fractional gap junction open state at each time
        self.cc_er_time = []   # retains er concentrations as a function of time
        self.cIP3_time = []    # retains cellular ip3 concentrations as a function of time
        self.osmo_P_delta_time = []  # osmotic pressure difference between cell interior and exterior as func of time

        self.I_mem_time = []    # initialize membrane current time vector

        self.vm_Matrix = []    # initialize matrices for resampled data sets (used in smooth plotting and streamlines)
        self.I_gj_x_time = []
        self.I_gj_y_time = []

        self.efield_gj_x_time = []   # matrices storing smooth electric field in gj connected cells
        self.efield_gj_y_time = []

        self.P_cells_time = []
        self.u_cells_x_time = []
        self.u_cells_y_time = []

        self.rho_cells_time = []

        self.I_tot_x_time = [0]
        self.I_tot_y_time = [0]

        self.F_electro_time = []
        self.F_electro_x_time = []
        self.F_electro_y_time = []
        self.F_hydro_x_time = []
        self.F_hydro_y_time = []
        self.F_hydro_time = []
        self.P_electro_time = []

        self.rate_NaKATP_time =[]

        # if p.deform_osmo is True and p.run_sim is False:
        #
        #     self.P_cells = np.zeros(len(cells.cell_i)) #initialize total pressure in cells, but only during an init

        if p.deformation is True and p.run_sim is True:

            self.ecm_verts_unique_to = cells.ecm_verts_unique[:] # make a copy of original ecm verts as disp ref point

            self.cell_centres_time = []
            self.mem_mids_time = []
            self.maskM_time = []
            self.mem_edges_time = []
            self.cell_verts_time = []

            self.dx_cell_time = []
            self.dy_cell_time = []

            self.dx_time = []
            self.dy_time = []

            self.phi = np.zeros(len(cells.cell_i))
            self.phi_time = []

        if p.voltage_dye is True:

            self.cDye_time = []    # retains voltage-sensitive dye concentration as a function of time

        if p.scheduled_options['IP3'] != 0 or p.Ca_dyn is True:

            self.cIP3_time = []    # retains IP3 concentration as a function of time

        # gap junction specific arrays:
        self.id_gj = np.ones(len(cells.nn_i))  # identity array for gap junction indices...
        self.gjopen = np.ones(len(cells.nn_i))   # holds gap junction open fraction for each gj
        self.gjl = np.zeros(len(cells.nn_i))    # gj length for each gj
        self.gjl[:] = cells.gj_len

        self.rho_pump = np.ones(len(cells.mem_i))

        # get the net, unbalanced charge and corresponding voltage in each cell:
        self.update_V_ecm(cells,p,0)

        # vm_to = copy.deepcopy(self.vm)   # create a copy of the original voltage
        self.vm_to = self.vm[:]

        if p.Ca_dyn is True:
            self.cc_er = np.asarray(self.cc_er)
            self.cc_er_to = self.cc_er[:]

        # create a time-steps vector appropriate for the simulation type:
        if p.run_sim is True:
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

        if p.run_sim is True:


            loggers.log_info('Your simulation is running from '+ str(0) + ' to '+ str(round(p.sim_tsteps*p.dt,3))
                         + ' seconds of in-world time.')

        else:
             loggers.log_info('Your initialization is running from '+ str(0) + ' to '+ str(round(p.init_tsteps*p.dt,3))
                         + ' seconds of in-world time.')


        if p.plot_while_solving is True:

            self.checkPlot = viz.PlotWhileSolving(cells,self,p,clrAutoscale = p.autoscale_Vmem, clrMin = p.Vmem_min_clr,
                clrMax = p.Vmem_max_clr)

        do_once = True  # a variable to time the loop only once

        for t in tt:   # run through the loop

            if do_once is True:
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
            if p.run_sim is True:
                self.dyna.runAllDynamics(self,cells,p,t)

            # run the Na-K-ATPase pump:
            fNa_NaK, fK_NaK, self.rate_NaKATP = pumpNaKATP(self.cc_cells[self.iNa],self.cc_env[self.iNa],self.cc_cells[self.iK],
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

            # get forces from any hydrostatic (self.P_Cells) pressure:
            self.getHydroF(cells,p)

            # calculate pressures

            if p.deform_osmo is True:

                 self.osmotic_P(cells,p)

            if p.deform_electro is True:

                self.electro_P(cells,p)

            if p.gravity is True:

                self.gravity_P(cells,p)

            # calculate fluid flow:
            if p.fluid_flow is True:

                self.getFlow(cells,p)

            if p.deformation is True and p.run_sim is True:

                self.run_sim = True

                if p.td_deform is False:

                    self.getDeformation(cells,t,p)

                elif p.td_deform is True:

                    self.timeDeform(cells,t,p)

            if p.scheduled_options['IP3'] != 0 or p.Ca_dyn is True:
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

                if p.GHK_calc is True:

                    self.ghk_calculator(cells,p)
                    self.vm_GHK_time.append(self.vm_GHK) # data array holding GHK vm estimates

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

                ddc = np.copy(self.Dm_cells[:])
                ddc.tolist()
                self.dd_time.append(ddc)
                ddc = None

                self.I_gj_x_time.append(self.I_gj_x[:])
                self.I_gj_y_time.append(self.I_gj_y[:])

                self.I_tot_x_time.append(self.I_tot_x[:])
                self.I_tot_y_time.append(self.I_tot_y[:])

                self.I_mem_time.append(self.I_mem[:])

                self.vm_time.append(self.vm[:])

                self.rho_cells_time.append(self.rho_cells[:])

                self.rate_NaKATP_time.append(self.rate_NaKATP[:])

                self.P_cells_time.append(self.P_cells[:])
                self.F_hydro_x_time.append(self.F_hydro_x[:])
                self.F_hydro_y_time.append(self.F_hydro_y[:])

                if p.deform_osmo is True:

                    self.osmo_P_delta_time.append(self.osmo_P_delta[:])


                if p.deform_electro is True:
                    self.F_electro_time.append(self.F_electro[:])
                    self.F_electro_x_time.append(self.F_electro_x[:])
                    self.F_electro_y_time.append(self.F_electro_y[:])

                    self.P_electro_time.append(self.P_electro[:])

                if p.deformation is True and p.run_sim is True:

                    self.implement_deform_timestep(cells,t,p)

                    self.cell_centres_time.append(cells.cell_centres[:])
                    self.mem_mids_time.append(cells.mem_mids_flat[:])
                    self.maskM_time.append(cells.maskM[:])
                    self.mem_edges_time.append(cells.mem_edges_flat[:])
                    self.cell_verts_time.append(cells.cell_verts[:])

                    self.dx_cell_time.append(self.d_cells_x[:])
                    self.dy_cell_time.append(self.d_cells_y[:])


                if p.fluid_flow is True:

                    self.u_cells_x_time.append(self.u_cells_x[:])
                    self.u_cells_y_time.append(self.u_cells_y[:])

                self.gjopen_time.append(self.gjopen[:])

                self.time.append(t)

                if p.scheduled_options['IP3'] != 0 or p.Ca_dyn is True:

                    self.cIP3_time.append(self.cIP3[:])

                if p.voltage_dye ==1:

                    self.cDye_time.append(self.cDye_cell[:])

                if p.Ca_dyn == 1 and p.ions_dict['Ca']==1:

                    self.cc_er_time.append(np.copy(self.cc_er[:]))

                if p.plot_while_solving is True:
                    self.checkPlot.updatePlot(self,p)

                        # get time for loop and estimate total time for simulation
            if do_once is True:
                loop_time = time.time() - loop_measure
                if p.run_sim is True:
                    time_estimate = round(loop_time*p.sim_tsteps,2)
                else:
                    time_estimate = round(loop_time*p.init_tsteps,2)
                loggers.log_info("This run should take approximately " + str(time_estimate) + ' s to compute...')
                do_once = False

        # Find embeded functions that can't be pickled...
        for key, valu in vars(self.dyna).items():
            if type(valu) == interp.interp1d or callable(valu):
                setattr(self.dyna,key,None)

        cells.points_tree = None

        self.checkPlot = None

        if p.run_sim is False:

            datadump = [self,cells,p]
            fh.saveSim(self.savedInit,datadump)
            message_1 = 'Initialization run saved to' + ' ' + p.init_path
            loggers.log_info(message_1)

        elif p.run_sim is True:

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


        final_vmean = 1000*np.round(np.mean(self.vm_time[-1]),6)
        vmess = 'Final average cell Vmem of ' + ': '
        loggers.log_info(vmess + str(final_vmean) + ' mV')

        if p.GHK_calc is True:
            final_vmean_GHK = 1000*np.round(np.mean(self.vm_GHK_time[-1]),6)
            vmess = 'Final average cell Vmem calculated using GHK: ' + ': '
            loggers.log_info(vmess + str(final_vmean_GHK) + ' mV')


        if p.ions_dict['H'] == 1:
            final_pH = -np.log10(np.mean((self.cc_time[-1][self.iH])))
            loggers.log_info('Final average cell pH '+ str(np.round(final_pH,2)))

            final_pH_env = -np.log10(np.mean((self.cc_env_time[-1][self.iH])))
            loggers.log_info('Final environmental pH '+ str(np.round(final_pH_env,2)))


        if p.scheduled_options['IP3'] != 0 or p.Ca_dyn is True:
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

        self.dd_time = []  # data array storing cell membrane permeabilities

        self.vm_time = []  # data array holding voltage at time points
        self.vcell_time = []
        self.venv_time = []

        self.vm_GHK_time = [] # data array holding GHK vm estimate

        # self.dvm_time = []  # data array holding derivative of voltage at time points
        self.time = []     # time values of the simulation

        self.gjopen_time = []   # stores the fractional gap junction open state at each time

        self.I_gj_x_time = []    # initialize gap junction current data storage
        self.I_gj_y_time = []    # initialize gap junction current data storage
        self.I_tot_x_time = []
        self.I_tot_y_time = []

        self.I_mem_time = []   # initialize membrane matrix data storage

        self.cc_er_time = []   # retains er concentrations as a function of time
        self.cIP3_time = []    # retains cellular ip3 concentrations as a function of time
        self.cIP3_env_time = []

        self.efield_gj_x_time = []   # matrices storing smooth electric field in gj connected cells
        self.efield_gj_y_time = []

        self.efield_ecm_x_time = []   # matrices storing smooth electric field in ecm
        self.efield_ecm_y_time = []

        # initialize time-storage vectors for electroosmotic data:
        self.u_env_x_time = []
        self.u_env_y_time = []

        self.rate_NaKATP_time = []

        self.u_cells_x_time = []
        self.u_cells_y_time = []
        self.P_cells_time = []

        self.osmo_P_delta_time = []  # osmotic pressure difference between cell interior and exterior as func of time

        self.rho_cells_time = []

        self.rho_pump_time = []    # store pump and channel states as function of time...
        self.rho_channel_time = []

        self.vm_Matrix = [] # initialize matrices for resampled data sets (used in smooth plotting and streamlines)
        vm_dato = np.zeros(len(cells.mem_i))
        dat_grid_vm = vertData(vm_dato,cells,p)
        self.vm_Matrix.append(dat_grid_vm[:])

        self.F_electro_time = []
        self.F_electro_x_time = []
        self.F_electro_y_time = []
        self.F_hydro_x_time = []
        self.F_hydro_y_time = []
        self.P_electro_time = []

        # if p.deform_osmo is True and p.run_sim is False:
        #
        #     self.P_cells = np.zeros(len(cells.cell_i)) #initialize total pressure in cells

        if p.deformation is True and p.run_sim is True:

            self.ecm_verts_unique_to = cells.ecm_verts_unique[:]

            self.cell_centres_time = []
            self.mem_mids_time = []
            self.maskM_time = []
            self.mem_edges_time = []
            self.cell_verts_time = []

            self.dx_cell_time = []
            self.dy_cell_time = []

            self.dx_time = []
            self.dy_time = []

        if p.Ca_dyn is True:
            self.cc_er_to = np.copy(self.cc_er[:])

        if p.voltage_dye is True:
            self.cDye_time = []    # retains voltage-sensitive dye concentration as a function of time
            self.cDye_env_time = []

        # gap junction specific arrays:
        self.id_gj = np.ones(len(cells.nn_i))  # identity array for gap junction indices...
        self.gjopen = np.ones(len(cells.nn_i))   # holds gap junction open fraction for each gj
        self.gjl = np.zeros(len(cells.nn_i))    # gj length for each gj
        self.gjl[:] = cells.gj_len

        # get the net, unbalanced charge and corresponding voltage in each cell to initialize values of voltages:
        self.update_V_ecm(cells,p,0)

        self.vm_to = self.vm[:]   # create a copy of the original voltage

         # create a time-steps vector appropriate for simulation type:
        if p.run_sim is True:
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
        if p.run_sim is True:

            loggers.log_info('Your simulation (with extracellular spaces) is running from '+ str(0) +
                             ' to '+ str(round(p.sim_tsteps*p.dt,3))
                         + ' seconds of in-world time.')

        else:
            loggers.log_info('Your initialization (with extracellular spaces) is running from '+ str(0) +
                             ' to '+ str(round(p.init_tsteps*p.dt,3))
                         + ' seconds of in-world time.')


        if p.plot_while_solving is True:

            self.checkPlot = viz.PlotWhileSolving(cells,self,p,clrAutoscale = p.autoscale_Vmem, clrMin = p.Vmem_min_clr,
                clrMax = p.Vmem_max_clr)

        do_once = True  # a variable to time the loop only once

        for t in tt:   # run through the loop

            if do_once is True:
                loop_measure = time.time()

            self.fluxes_mem.fill(0)  # reinitialize flux storage device

            self.dvm = (self.vm - self.vm_to)/p.dt    # calculate the change in the voltage derivative
            self.vm_to = self.vm[:]       # reassign the history-saving vm

            if p.Ca_dyn ==1 and p.ions_dict['Ca'] == 1:

                self.dcc_ER = (self.cc_er - self.cc_er_to)/p.dt
                self.cc_er_to = self.cc_er[:]

            # calculate the values of scheduled and dynamic quantities (e.g. ion channel multipliers):
            if p.run_sim is True:
                self.dyna.runAllDynamics(self,cells,p,t)

            #-----------------PUMPS-------------------------------------------------------------------------------------

            # run the Na-K-ATPase pump:
            fNa_NaK, fK_NaK, self.rate_NaKATP =\
                pumpNaKATP(self.cc_cells[self.iNa][cells.mem_to_cells],self.cc_env[self.iNa][cells.map_mem2ecm],
                    self.cc_cells[self.iK][cells.mem_to_cells],self.cc_env[self.iK][cells.map_mem2ecm],
                    self.vm,self.T,p,self.NaKATP_block)

            # modify the resulting fluxes by the electroosmosis membrane redistribution factor (if calculated)
            fNa_NaK = self.rho_pump*fNa_NaK
            fK_NaK = self.rho_pump*fK_NaK

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

            # get forces from any hydrostatic (self.P_Cells) pressure:
            self.getHydroF(cells,p)

            # calculate pressures:

            if p.deform_osmo is True:

                self.osmotic_P(cells,p)

            if p.deform_electro is True:

                self.electro_P(cells,p)

            if p.gravity is True:

                self.gravity_P(cells,p)

            if p.fluid_flow is True:

                self.getFlow(cells,p)

            # if desired, electroosmosis of membrane channels
            if p.sim_eosmosis is True and cells.gradMem is not None:

                self.eosmosis(cells,p)    # modify membrane pump and channel density according to Nernst-Planck


            if p.deformation is True and p.run_sim is True:

                self.run_sim = True

                if p.td_deform is False:

                    self.getDeformation(cells,t,p)

                elif p.td_deform is True:

                    self.timeDeform(cells,t,p)

            if p.scheduled_options['IP3'] != 0 or p.Ca_dyn is True:

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


            if t in tsamples:

                if p.GHK_calc is True:

                    self.ghk_calculator(cells,p)
                    self.vm_GHK_time.append(self.vm_GHK)

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

                ddc = np.copy(self.Dm_cells[:])
                ddc.tolist()
                self.dd_time.append(ddc)
                ddc = None

                self.vm_time.append(self.vm[:])

                # self.dvm_time.append(self.dvm[:])

                self.gjopen_time.append(self.gjopen[:])

                self.vcell_time.append(self.v_cell[:])

                self.venv_time.append(self.v_env[:])

                self.I_gj_x_time.append(self.I_gj_x[:])
                self.I_gj_y_time.append(self.I_gj_y[:])

                self.I_mem_time.append(self.I_mem[:])

                self.I_tot_x_time.append(self.I_tot_x)
                self.I_tot_y_time.append(self.I_tot_y)

                self.rho_cells_time.append(self.rho_cells[:])

                self.rate_NaKATP_time.append(self.rate_NaKATP[:])

                if p.fluid_flow is True:
                    self.u_env_x_time.append(self.u_env_x[:])
                    self.u_env_y_time.append(self.u_env_y[:])

                    self.u_cells_x_time.append(self.u_cells_x[:])
                    self.u_cells_y_time.append(self.u_cells_y[:])

                # calculate interpolated verts and midpoint data for Vmem:
                dat_grid_vm = vertData(self.vm[:],cells,p)

                self.vm_Matrix.append(dat_grid_vm[:])

                self.time.append(t)

                self.P_cells_time.append(self.P_cells[:])
                self.F_hydro_x_time.append(self.F_hydro_x[:])
                self.F_hydro_y_time.append(self.F_hydro_y[:])

                if p.deform_osmo is True: # if osmotic pressure is enabled
                    self.osmo_P_delta_time.append(self.osmo_P_delta[:])

                if p.deform_electro is True:
                    self.F_electro_time.append(self.F_electro[:])
                    self.F_electro_x_time.append(self.F_electro_x[:])
                    self.F_electro_y_time.append(self.F_electro_y[:])
                    self.P_electro_time.append(self.P_electro[:])

                if p.deformation is True and p.run_sim is True:

                    self.implement_deform_timestep(cells,t,p)

                    self.cell_centres_time.append(cells.cell_centres[:])
                    self.mem_mids_time.append(cells.mem_mids_flat[:])
                    self.maskM_time.append(cells.maskM[:])
                    self.mem_edges_time.append(cells.mem_edges_flat[:])
                    self.cell_verts_time.append(cells.cell_verts[:])

                    self.dx_cell_time.append(self.d_cells_x[:])
                    self.dy_cell_time.append(self.d_cells_y[:])


                if p.scheduled_options['IP3'] != 0 or p.Ca_dyn is True:
                    ccIP3 = self.cIP3[:]
                    self.cIP3_time.append(ccIP3)
                    ccIP3 = None

                if p.voltage_dye ==1:
                    ccDye_cells = self.cDye_cell[:]
                    self.cDye_time.append(ccDye_cells)
                    ccDye_cells = None

                    self.cDye_env_time.append(self.cDye_env[:])

                if p.Ca_dyn == 1 and p.ions_dict['Ca']==1:
                    self.cc_er_time.append(np.copy(self.cc_er[:]))

                if p.sim_eosmosis is True:

                    self.rho_channel_time.append(self.rho_channel[:])
                    self.rho_pump_time.append(self.rho_pump[:])

                if p.plot_while_solving is True:
                    self.checkPlot.updatePlot(self,p)

                        # get time for loop and estimate total time for simulation
            if do_once is True:
                loop_time = time.time() - loop_measure

                if p.run_sim is True:
                    time_estimate = round(loop_time*p.sim_tsteps,2)
                else:
                    time_estimate = round(loop_time*p.init_tsteps,2)
                loggers.log_info("This run should take approximately " + str(time_estimate) + ' s to compute...')
                do_once = False

         # Find embeded functions that can't be pickled...
        for key, valu in vars(self.dyna).items():
            if type(valu) == interp.interp1d or callable(valu):
                setattr(self.dyna,key,None)

        cells.points_tree = None

        self.checkPlot = None

        if p.run_sim is False:

            datadump = [self,cells,p]
            fh.saveSim(self.savedInit,datadump)
            message_1 = 'Initialization run saved to' + ' ' + p.init_path
            loggers.log_info(message_1)

        elif p.run_sim is True:
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

        final_vmean = 1000*np.round(np.mean(self.vm_time[-1]),6)
        vmess = 'Final average cell Vmem: ' + ': '
        loggers.log_info(vmess + str(final_vmean) + ' mV')

        if p.GHK_calc is True:

            final_vmean_GHK = 1000*np.round(np.mean(self.vm_GHK_time[-1]),6)
            vmess = 'Final average cell Vmem calculated using GHK equation: ' + ': '
            loggers.log_info(vmess + str(final_vmean_GHK) + ' mV')

        if p.ions_dict['H'] == 1:
            final_pH = -np.log10(np.mean((self.cc_time[-1][self.iH])))
            loggers.log_info('Final average cell pH '+ str(np.round(final_pH,2)))

            final_pH_ecm = -np.log10(np.mean((self.cc_env_time[-1][self.iH])))
            loggers.log_info('Final extracellular pH '+ str(np.round(final_pH_ecm,2)))

        if p.scheduled_options['IP3'] != 0 or p.Ca_dyn is True:

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


        if p.sim_ECM is True:

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
        d_c_env = -flux*(cells.mem_sa/cells.ecm_vol[cells.map_mem2ecm])

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
        if p.sim_ECM is True:

            vmems = self.v_cell[cells.mem_to_cells]
            # vmems = self.vm

            self.vgj = vmems[cells.nn_i]- vmems[cells.mem_i]

        else:

            vmems = self.vm[cells.mem_to_cells]

            self.vgj = vmems[cells.nn_i]- vmems[cells.mem_i]


        if p.gj_flux_sensitive is True and cells.nnAveMatrix is not None:

            if p.gj_respond_flow is True:

                # map intracellular flow velocity to gap junctions:
                ux_gj = (self.u_cells_x[cells.mem_i] + self.u_cells_x[cells.nn_i])/2
                uy_gj = (self.u_cells_y[cells.mem_i] + self.u_cells_y[cells.nn_i])/2

                # get the magnitude of flow at the gap junction:
                u_gj_s = np.sqrt(ux_gj**2 + uy_gj**2)

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


        if p.v_sensitive_gj is True:
            # determine the open state of gap junctions:
            self.gjopen = self.gj_rho + self.gj_block*((1.0 - tb.step(abs(self.vgj),p.gj_vthresh,p.gj_vgrad) + 0.1))
            # self.gjopen = self.gj_rho + self.gj_block*((1.0 - tb.hill(abs(self.vgj),p.gj_vthresh,p.gj_vgrad) + 0.1))

        else:
            self.gjopen = (1 + self.gj_rho)*self.gj_block


        # voltage gradient:
        grad_vgj = self.vgj/cells.gj_len

        grad_vgj_x = grad_vgj*cells.cell_nn_tx
        grad_vgj_y = grad_vgj*cells.cell_nn_ty

        # concentration gradient for ion i:

        conc_mem = self.cc_cells[i][cells.mem_to_cells]
        grad_cgj = (conc_mem[cells.nn_i] - conc_mem[cells.mem_i])/cells.gj_len

        grad_cgj_x = grad_cgj*cells.cell_nn_tx
        grad_cgj_y = grad_cgj*cells.cell_nn_ty

        # midpoint concentration:
        c = (conc_mem[cells.nn_i] + conc_mem[cells.mem_i])/2

        # electroosmotic fluid velocity -- averaged at gap junctions:
        if p.fluid_flow is True:
            ux = (self.u_cells_x[cells.cell_nn_i[:,0]] + self.u_cells_x[cells.cell_nn_i[:,1]])/2
            uy = (self.u_cells_y[cells.cell_nn_i[:,0]] + self.u_cells_y[cells.cell_nn_i[:,1]])/2

        else:
            ux = 0
            uy =0

        fgj_x,fgj_y = nernst_planck_flux(c,grad_cgj_x,grad_cgj_y,grad_vgj_x,grad_vgj_y,ux,uy,
            p.gj_surface*self.gjopen*self.D_gj[i],self.zs[i],self.T,p)


        fgj = fgj_x*cells.cell_nn_tx + fgj_y*cells.cell_nn_ty

        delta_cc = np.dot(cells.gjMatrix,-fgj*cells.mem_sa)/cells.cell_vol

        self.cc_cells[i] = self.cc_cells[i] + p.dt*delta_cc

        self.fluxes_gj_x[i] = fgj_x  # store gap junction flux for this ion
        self.fluxes_gj_y[i] = fgj_y  # store gap junction flux for this ion

    def update_ecm(self,cells,p,t,i):

        if p.closed_bound is True:
            btag = 'closed'

        else:
            btag = 'open'

        # make v_env and cc_env into 2d matrices
        cenv = self.cc_env[i][:]
        # denv = self.D_env[i][:]

        v_env = self.v_env[:].reshape(cells.X.shape)

        # enforce voltage at boundary:
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

        if p.closed_bound is True: # insulation boundary conditions

            pass
            # cenv_x[:,0] = cenv_x[:,1]
            # cenv_x[:,-1] = cenv_x[:,-2]
            # cenv_x[0,:] = cenv_x[1,:]
            # cenv_x[-1,:] = cenv_x[-2,:]
            #
            # cenv_y[0,:] = cenv_y[1,:]
            # cenv_y[-1,:] = cenv_y[-2,:]
            # cenv_y[:,0] = cenv_y[:,1]
            # cenv_y[:,-1] = cenv_y[:,-2]

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
        if p.fluid_flow is True:

            uenvx = np.zeros(cells.grid_obj.u_shape)
            uenvy = np.zeros(cells.grid_obj.v_shape)

            uenvx[:,1:] = self.u_env_x
            uenvy[1:,:] = self.u_env_y

            if p.closed_bound is False:

                uenvx[:,0] = uenvx[:,1]
                uenvx[:,-1]= uenvx[:,-2]
                uenvx[0,:] = uenvx[1,:]
                uenvx[-1,:] = uenvx[-2,:]

                uenvy[:,0] = uenvy[:,1]
                uenvy[:,-1]= uenvy[:,-2]
                uenvy[0,:] = uenvy[1,:]
                uenvy[-1,:] = uenvy[-2,:]

            else:

                uenvx[:,0] = 0
                uenvx[:,-1]= 0
                uenvx[0,:] = 0
                uenvx[-1,:] = 0

                uenvy[:,0] = 0
                uenvy[:,-1]= 0
                uenvy[0,:] = 0
                uenvy[-1,:] = 0


        else:
            uenvx = 0
            uenvy = 0

        f_env_x, f_env_y = np_flux_special(cenv_x,cenv_y,grad_cc_env_x,grad_cc_env_y,
            grad_V_env_x, grad_V_env_y, uenvx,uenvy,self.D_env_u[i],self.D_env_v[i],
            self.zs[i],self.T,p)

        if p.closed_bound is False:

            pass
            # f_env_x[:,0] = f_env_x[:,1]
            # f_env_x[:,-1]= f_env_x[:,-2]
            # f_env_x[0,:] = f_env_x[1,:]
            # f_env_x[-1,:] = f_env_x[-2,:]
            #
            # f_env_y[:,0] = f_env_y[:,1]
            # f_env_y[:,-1]= f_env_y[:,-2]
            # f_env_y[0,:] = f_env_y[1,:]
            # f_env_y[-1,:] = f_env_y[-2,:]

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

        # delta_c = fd.integrator(delta_c)

        #-----------------------
        cenv = cenv + delta_c*p.dt

        if p.closed_bound is True:
            # Neumann boundary condition (flux at boundary)
            # zero flux boundaries for concentration:
            pass
            # cenv[:,-1] = cenv[:,-2]
            # cenv[:,0] = cenv[:,1]
            # cenv[0,:] = cenv[1,:]
            # cenv[-1,:] = cenv[-2,:]

        elif p.closed_bound is False:
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

        fenvx = fd.integrator(fenvx)
        fenvy = fd.integrator(fenvy)

        self.fluxes_env_x[i] = fenvx.ravel()  # store ecm junction flux for this ion
        self.fluxes_env_y[i] = fenvy.ravel()  # store ecm junction flux for this ion

    def update_ecm_o(self,cells,p,t,i):

        if p.closed_bound is True:
            btag = 'closed'

        else:
            btag = 'open'


        # make v_env, d_env, and cc_env into 2d matrices
        cenv = self.cc_env[i][:].reshape(cells.X.shape)
        denv = self.D_env[i][:].reshape(cells.X.shape)

        v_env = self.v_env[:].reshape(cells.X.shape)

        # enforce boundary conditions on denv:
        if p.closed_bound is True:
            denv[:,-1] = 0
            denv[:,0] = 0
            denv[0,:] = 0
            denv[-1,:] = 0

        # enforce voltage at boundary:
        v_env[:,0] = self.bound_V['L']
        v_env[:,-1] = self.bound_V['R']
        v_env[0,:] = self.bound_V['B']
        v_env[-1,:] = self.bound_V['T']

        # calculate gradients in the environment
        grad_V_env_x, grad_V_env_y = fd.gradient(v_env,cells.delta)

        grad_cc_env_x, grad_cc_env_y = fd.gradient(cenv,cells.delta)

        # calculate fluxes for electrodiffusive transport:
        if p.fluid_flow is True:
            uenvx = self.u_env_x
            uenvy = self.u_env_y

        else:
            uenvx = 0
            uenvy = 0

        f_env_x, f_env_y = np_flux_special(cenv,cenv,grad_cc_env_x,grad_cc_env_y,
            grad_V_env_x, grad_V_env_y, uenvx,uenvy,denv,denv,
            self.zs[i],self.T,p)

        # set boundary conditions on resulting flux

        if p.closed_bound is False: # no normal gradient
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

        # smooth out the fluxes by integration:
        f_env_x = cells.grid_obj.grid_int(f_env_x, bounds='btag')
        f_env_y = cells.grid_obj.grid_int(f_env_y, bounds='btag')

        # calculate the negative divergence of the total flux (amount entering area per unit time):

        delta_c = -fd.divergence(f_env_x,f_env_y,cells.delta,cells.delta)

        # smooth out the delta_c by integration:
        delta_c = cells.grid_obj.grid_int(delta_c, bounds='btag')

        #-----------------------
        cenv = cenv + delta_c*p.dt

        if p.closed_bound is True:
            # Neumann boundary condition (flux at boundary)
            # zero flux boundaries for concentration:
            cenv[:,-1] = cenv[:,-2]
            cenv[:,0] = cenv[:,1]
            cenv[0,:] = cenv[1,:]
            cenv[-1,:] = cenv[-2,:]

        elif p.closed_bound is False:
            # if the boundary is open, set the concentration at the boundary
            # open boundary
            cenv[:,-1] = self.c_env_bound[i]
            cenv[:,0] = self.c_env_bound[i]
            cenv[0,:] = self.c_env_bound[i]
            cenv[-1,:] = self.c_env_bound[i]

        # reshape the matrices back into vectors:
        self.cc_env[i] = cenv.ravel()

        fenvx = f_env_x[:]
        fenvy = f_env_y[:]

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

        #------------------------------------------------------------
        # pump dye and update result

        if p.sim_ECM is True:

            if p.pump_Dye is True:

                if p.pump_Dye_out is True:

                    # active pumping of dye from environment and into cell
                    deltaGATP = 20*p.R*self.T

                    delG_dye = p.R*p.T*np.log(self.cDye_cell[cells.mem_to_cells]/self.cDye_env[cells.map_mem2ecm]) \
                               + p.z_Dye*p.F*self.vm

                    delG_dyeATP = deltaGATP - delG_dye
                    delG_pump = (delG_dyeATP/1000)

                    # alpha = p.pump_Dye_alpha*tb.step(delG_pump,6,3)
                    alpha = p.pump_Dye_alpha*delG_pump

                    f_Dye_pump  = self.rho_pump*alpha*(self.cDye_env[cells.map_mem2ecm])

                    d_dye_cells = f_Dye_pump*(cells.mem_sa/cells.cell_vol[cells.mem_to_cells])
                    d_dye_env = -f_Dye_pump*(cells.mem_sa/cells.ecm_vol[cells.map_mem2ecm])

                    delta_cells =  np.dot(d_dye_cells, cells.cell_UpdateMatrix)
                    delta_env = np.dot(d_dye_env, cells.ecm_UpdateMatrix)

                else:

                     # active pumping of dye from environment and into cell
                    deltaGATP = 20*p.R*self.T

                    delG_dye = p.R*p.T*np.log(self.cDye_env[cells.map_mem2ecm]/self.cDye_cell[cells.mem_to_cells]) \
                               - p.z_Dye*p.F*self.vm

                    delG_dyeATP = deltaGATP - delG_dye
                    delG_pump = (delG_dyeATP/1000)

                    # alpha = p.pump_Dye_alpha*tb.step(delG_pump,6,3)
                    alpha = p.pump_Dye_alpha*delG_pump

                    f_Dye_pump  = self.rho_pump*alpha*(self.cDye_cell[cells.mem_to_cells])

                    d_dye_cells = -f_Dye_pump*(cells.mem_sa/cells.cell_vol[cells.mem_to_cells])
                    d_dye_env = f_Dye_pump*(cells.mem_sa/cells.ecm_vol[cells.map_mem2ecm])

                    delta_cells =  np.dot(d_dye_cells, cells.cell_UpdateMatrix)
                    delta_env = np.dot(d_dye_env, cells.ecm_UpdateMatrix)

                self.cDye_cell = self.cDye_cell + delta_cells*p.dt

                self.cDye_env = self.cDye_env + delta_env*p.dt


        elif p.sim_ECM is False:

            if p.pump_Dye is True:

                if p.pump_Dye_out is True:

                    # active pumping of dye from environment and into cell
                    deltaGATP = 20*p.R*self.T

                    delG_dye = p.R*p.T*np.log(self.cDye_cell/self.cDye_env) \
                               + p.z_Dye*p.F*self.vm

                    delG_dyeATP = deltaGATP - delG_dye
                    delG_pump = (delG_dyeATP/1000)

                    # alpha = p.pump_Dye_alpha*tb.step(delG_pump,6,3)
                    alpha = p.pump_Dye_alpha*delG_pump

                    f_Dye_pump  = alpha*(self.cDye_env)

                    d_dye_cells = f_Dye_pump*(cells.cell_sa/cells.cell_vol)

                    # delta_cells =  np.dot(d_dye_cells, cells.cell_UpdateMatrix)


                else:

                     # active pumping of dye from cell and into environment
                    deltaGATP = 20*p.R*self.T

                    delG_dye = p.R*p.T*np.log(self.cDye_env/self.cDye_cell) \
                               - p.z_Dye*p.F*self.vm

                    delG_dyeATP = deltaGATP - delG_dye
                    delG_pump = (delG_dyeATP/1000)

                    # alpha = p.pump_Dye_alpha*tb.step(delG_pump,6,3)
                    alpha = p.pump_Dye_alpha*delG_pump

                    f_Dye_pump  = alpha*(self.cDye_cell)

                    d_dye_cells = -f_Dye_pump*(cells.cell_sa/cells.cell_vol)

                    # delta_cells =  np.dot(d_dye_cells, cells.cell_UpdateMatrix)

                self.cDye_cell = self.cDye_cell + d_dye_cells*p.dt

        #------------------------------------------------------------

        # Update dye concentration in the gj connected cell network:

        # voltage gradient:
        grad_vgj = self.vgj/cells.gj_len

        grad_vgj_x = grad_vgj*cells.cell_nn_tx
        grad_vgj_y = grad_vgj*cells.cell_nn_ty

        # concentration gradient for Dye:
        Dye_mems = self.cDye_cell[cells.mem_to_cells]

        grad_cgj = (Dye_mems[cells.nn_i] - Dye_mems[cells.mem_i])/cells.gj_len

        grad_cgj_x = grad_cgj*cells.cell_nn_tx
        grad_cgj_y = grad_cgj*cells.cell_nn_ty

        # midpoint concentration:
        cdye = (Dye_mems[cells.nn_i] + Dye_mems[cells.mem_i])/2

        # electroosmotic fluid velocity:
        if p.fluid_flow is True:

            ux = (self.u_cells_x[cells.cell_nn_i[:,0]] + self.u_cells_x[cells.cell_nn_i[:,1]])/2
            uy = (self.u_cells_y[cells.cell_nn_i[:,0]] + self.u_cells_y[cells.cell_nn_i[:,1]])/2

        else:
            ux = 0
            uy = 0

        fgj_x_dye,fgj_y_dye = nernst_planck_flux(cdye,grad_cgj_x,grad_cgj_y,grad_vgj_x,grad_vgj_y,ux,uy,
            p.Do_Dye*self.gjopen,p.z_Dye,self.T,p)

        fgj_dye = fgj_x_dye*cells.cell_nn_tx + fgj_y_dye*cells.cell_nn_ty

        delta_cc = np.dot(cells.gjMatrix*p.gj_surface*self.gjopen,-fgj_dye*cells.mem_sa)/cells.cell_vol

        self.cDye_cell = self.cDye_cell + p.dt*delta_cc

        self.Dye_flux_x_gj = fgj_x_dye[:]  # store gap junction flux for this ion
        self.Dye_flux_y_gj = fgj_y_dye[:]  # store gap junction flux for this ion

        if p.sim_ECM is False:

            fdye_ED = electroflux(self.cDye_env,self.cDye_cell,self.id_cells*p.Dm_Dye,self.tm,p.z_Dye,self.vm,self.T,p)

            # update dye concentration
            self.cDye_cell = self.cDye_cell + fdye_ED*(cells.cell_sa/cells.cell_vol)*p.dt

        elif p.sim_ECM is True:

            flux_dye = electroflux(self.cDye_env[cells.map_mem2ecm],self.cDye_cell[cells.mem_to_cells],
                            np.ones(len(cells.mem_i))*p.Dm_Dye,self.tm,p.z_Dye,self.vm,self.T,p)

            # update the dye concentrations in the cell and ecm due to ED fluxes at membrane
            d_c_cells = self.rho_channel*flux_dye*(cells.mem_sa/cells.cell_vol[cells.mem_to_cells])
            d_c_env = self.rho_channel*flux_dye*(cells.mem_sa/cells.ecm_vol[cells.map_mem2ecm])

            delta_cells =  np.dot(d_c_cells, cells.cell_UpdateMatrix)
            delta_env = np.dot(d_c_env, cells.ecm_UpdateMatrix)

            self.cDye_cell = rk4(self.cDye_cell,delta_cells,p)

            self.cDye_env = rk4(self.cDye_env,-delta_env,p)

            # transport dye through environment: _________________________________________________________
            if p.closed_bound is True:
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

            if p.closed_bound is True: # insulation boundary conditions
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

            if p.fluid_flow is True:
                uenvx = self.u_env_x
                uenvy = self.u_env_y

            else:
                uenvx = 0
                uenvy = 0

            f_env_x_dye, f_env_y_dye = np_flux_special(cenv_x,cenv_y,grad_cc_env_x,grad_cc_env_y,
                grad_V_env_x, grad_V_env_y, uenvx,uenvy,denv_x,denv_y,p.z_Dye,self.T,p)

            # calculate the divergence of the total flux, which is equivalent to the total change per unit time:
            d_fenvx = -(f_env_x_dye[:,1:] - f_env_x_dye[:,0:-1])/cells.delta
            d_fenvy = -(f_env_y_dye[1:,:] - f_env_y_dye[0:-1,:])/cells.delta

            delta_c = d_fenvx + d_fenvy

            cenv = cenv + delta_c*p.dt

            if p.closed_bound is True:
                # Neumann boundary condition (flux at boundary)
                # zero flux boundaries for concentration:
                cenv[:,-1] = cenv[:,-2]
                cenv[:,0] = cenv[:,1]
                cenv[0,:] = cenv[1,:]
                cenv[-1,:] = cenv[-2,:]

            elif p.closed_bound is False:
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
        grad_vgj = self.vgj/cells.gj_len

        grad_vgj_x = grad_vgj*cells.cell_nn_tx
        grad_vgj_y = grad_vgj*cells.cell_nn_ty

        # concentration gradient for Dye:

        IP3mem = self.cIP3[cells.mem_to_cells]

        grad_cgj = (IP3mem[cells.nn_i] - IP3mem[cells.mem_i])/cells.gj_len

        grad_cgj_x = grad_cgj*cells.cell_nn_tx
        grad_cgj_y = grad_cgj*cells.cell_nn_ty

        # midpoint concentration:
        cip3 = (IP3mem[cells.nn_i] + IP3mem[cells.mem_i])/2

        # electroosmotic fluid velocity:
        if p.fluid_flow is True:
            ux = (self.u_cells_x[cells.cell_nn_i[:,0]] + self.u_cells_x[cells.cell_nn_i[:,1]])/2
            uy = (self.u_cells_y[cells.cell_nn_i[:,0]] + self.u_cells_y[cells.cell_nn_i[:,1]])/2

        else:
            ux = 0
            uy = 0

        fgj_x_ip3,fgj_y_ip3 = nernst_planck_flux(cip3,grad_cgj_x,grad_cgj_y,grad_vgj_x,grad_vgj_y,ux,uy,
            p.Do_IP3*self.gjopen,p.z_IP3,self.T,p)

        fgj_ip3 = fgj_x_ip3*cells.cell_nn_tx + fgj_y_ip3*cells.cell_nn_ty

        delta_cc = np.dot(cells.gjMatrix*p.gj_surface*self.gjopen,-fgj_ip3*cells.mem_sa)/cells.cell_vol

        self.cIP3 = self.cIP3 + p.dt*delta_cc

        self.IP3_flux_x_gj = fgj_x_ip3[:]  # store gap junction flux for this ion
        self.IP3_flux_y_gj = fgj_y_ip3[:]  # store gap junction flux for this ion

        if p.sim_ECM is False:

            fip3_ED = electroflux(self.cIP3_env,self.cIP3,self.id_cells*p.Dm_IP3,self.tm,p.z_IP3,self.vm,self.T,p)

            # update dye concentration
            self.cIP3 = self.cIP3 + fip3_ED*(cells.cell_sa/cells.cell_vol)*p.dt

        elif p.sim_ECM is True:

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
            if p.closed_bound is True:
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

            if p.closed_bound is True: # insulation boundary conditions
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

            if p.fluid_flow is True:
                uenvx = self.u_env_x
                uenvy = self.u_env_y

            else:
                uenvx = 0
                uenvy = 0

            f_env_x_ip3, f_env_y_ip3 = np_flux_special(cenv_x,cenv_y,grad_cc_env_x,grad_cc_env_y,
                grad_V_env_x, grad_V_env_y, uenvx,uenvy,denv_x,denv_y,p.z_IP3,self.T,p)

            # calculate the divergence of the total flux, which is equivalent to the total change per unit time:
            d_fenvx = -(f_env_x_ip3[:,1:] - f_env_x_ip3[:,0:-1])/cells.delta
            d_fenvy = -(f_env_y_ip3[1:,:] - f_env_y_ip3[0:-1,:])/cells.delta

            delta_c = d_fenvx + d_fenvy

            cenv = cenv + delta_c*p.dt

            if p.closed_bound is True:
                # Neumann boundary condition (flux at boundary)
                # zero flux boundaries for concentration:
                cenv[:,-1] = cenv[:,-2]
                cenv[:,0] = cenv[:,1]
                cenv[0,:] = cenv[1,:]
                cenv[-1,:] = cenv[-2,:]

            elif p.closed_bound is False:
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
        if p.sim_ECM is True:

            vmem = self.v_cell[cells.mem_to_cells]

            self.Egj = - (vmem[cells.nn_i]- vmem[cells.mem_i])/cells.gj_len

            # in the environment:
            venv = self.v_env.reshape(cells.X.shape)
            genv_x, genv_y = fd.gradient(venv, cells.delta)

            self.E_env_x = -genv_x
            self.E_env_y = -genv_y

        else:

            # vmem = self.vm[cells.mem_to_cells]
            self.Egj = - (self.vm[cells.cell_nn_i[:,1]] - self.vm[cells.cell_nn_i[:,0]])/cells.gj_len

        # get x and y components of the electric field:
        self.E_gj_x = cells.cell_nn_tx*self.Egj
        self.E_gj_y = cells.cell_nn_ty*self.Egj

    def get_Bfield(self,cells,p):

        axo = self.Ax_time[-1]
        ayo = self.Ay_time[-1]

        sourceAx = -p.mu*self.I_tot_x.ravel()
        sourceAy = -p.mu*self.I_tot_y.ravel()

        self.Ax = np.dot(cells.lapENV_P_inv,sourceAx)
        self.Ay = np.dot(cells.lapENV_P_inv,sourceAy)

        self.Ax = self.Ax.reshape(cells.X.shape)
        self.Ay = self.Ay.reshape(cells.X.shape)

        # enforce zero normal gradient boundary conditions on A:
        self.Ax[:,0] = self.Ax[:,1]
        self.Ax[:,-1] = self.Ax[:,-2]
        self.Ax[0,:] = self.Ax[1,:]
        self.Ax[-1,:] = self.Ax[-2,:]

        self.Ay[:,0] = self.Ay[:,1]
        self.Ay[:,-1] = self.Ay[:,-2]
        self.Ay[0,:] = self.Ay[1,:]
        self.Ay[-1,:] = self.Ay[-2,:]

        self.Ax_time.append(self.Ax)
        self.Ay_time.append(self.Ay)

        self.dAx = (self.Ax - axo)/p.dt
        self.dAy = (self.Ay - ayo)/p.dt



        # print(np.max(self.dAx),np.max(self.dAy))

    def get_current(self,cells,p):

        # zero components of total current vector:

        self.I_tot_x = np.zeros(cells.Xgrid.shape)
        self.I_tot_y = np.zeros(cells.Xgrid.shape)

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
        self.I_gj_x = interp.griddata((cells.mem_mids_flat[:,0],cells.mem_mids_flat[:,1]),I_gj_x,(cells.Xgrid,cells.Ygrid),
                                      method=p.interp_type,fill_value=0)
        self.I_gj_x = np.multiply(self.I_gj_x,cells.maskM)

        self.I_gj_y = interp.griddata((cells.mem_mids_flat[:,0],cells.mem_mids_flat[:,1]),I_gj_y,(cells.Xgrid,cells.Ygrid),
                                      method=p.interp_type,fill_value=0)
        self.I_gj_y = np.multiply(self.I_gj_y,cells.maskM)

        self.I_tot_x = self.I_tot_x + self.I_gj_x
        self.I_tot_y = self.I_tot_y + self.I_gj_y

        # calculate current across cell membranes:

        self.I_mem = np.zeros(len(cells.mem_i))
        for flux_array, zi in zip(self.fluxes_mem,self.zs):

            # I_i = (flux_array*zi*p.F)/(self.gjopen*self.gjsa)
            I_i = flux_array*zi*p.F

            self.I_mem = self.I_mem + I_i

            I_mem_x = self.I_mem*cells.mem_vects_flat[:,2]
            I_mem_y = self.I_mem*cells.mem_vects_flat[:,3]

         # interpolate the trans-membrane current components to the grid:
        self.I_mem_x = interp.griddata((cells.mem_vects_flat[:,0],cells.mem_vects_flat[:,1]),I_mem_x,(cells.Xgrid,cells.Ygrid),
                                      method=p.interp_type,fill_value=0)
        self.I_mem_x = np.multiply(self.I_mem_x,cells.maskM)

        self.I_mem_y = interp.griddata((cells.mem_vects_flat[:,0],cells.mem_vects_flat[:,1]),I_mem_y,(cells.Xgrid,cells.Ygrid),
                                      method=p.interp_type,fill_value=0)
        self.I_mem_y = np.multiply(self.I_mem_y,cells.maskM)

        self.I_tot_x = self.I_tot_x + self.I_mem_x
        self.I_tot_y = self.I_tot_y + self.I_mem_y

        if p.sim_ECM is True:

            self.I_env_x = np.zeros(len(cells.xypts))
            self.I_env_y = np.zeros(len(cells.xypts))

            for flux_array, zi in zip(self.fluxes_env_x,self.zs):

                I_i = flux_array*zi*p.F

                self.I_env_x = self.I_env_x + I_i

            for flux_array, zi in zip(self.fluxes_env_y,self.zs):

                I_i = flux_array*zi*p.F

                self.I_env_y = self.I_env_y + I_i

            I_env_x = interp.griddata((cells.xypts[:,0],cells.xypts[:,1]),self.I_env_x,(cells.Xgrid,cells.Ygrid),
                                      method=p.interp_type,fill_value=0)

            I_env_y = interp.griddata((cells.xypts[:,0],cells.xypts[:,1]),self.I_env_y,(cells.Xgrid,cells.Ygrid),
                                      method=p.interp_type,fill_value=0)

            self.I_tot_x = self.I_tot_x + I_env_x
            self.I_tot_y = self.I_tot_y + I_env_y

    def getFlow(self,cells,p):
        """
        Calculate the electroosmotic fluid flow in the cell and extracellular
        networks.

        """

        if p.sim_ECM is True:

            # method 1-------------------------------------------------------------------------------------------------
            # force of gravity:
            if p.closed_bound is True:

                btag = 'closed'

            else:

                btag = 'open'

            # estimate the inverse viscosity for extracellular flow based on the diffusion constant weighting
            # for the world:
            alpha = (1/p.mu_water)*self.D_env_weight
            #---------------------------------------------------------

            if p.deform_electro is True:

                # calculate the electroosmotic force (shaped to the u and v grids):
                # first calculate the electric field on the u and v junctions:
                venv = self.v_env.reshape(cells.X.shape)

                # take the gradient of the venv and calculate the electric field
                env_x, env_y = fd.gradient(venv,cells.delta)

                # reshape and rescale the charge density:
                rho_env = self.rho_env.reshape(cells.X.shape)*(1/p.ff_env)

                # these are negative because the gradient of the voltage is the electric field and we just took the grad
                # above but didn't carry through the negative sign.
                Fe_x = -rho_env*env_x
                Fe_y = -rho_env*env_y

            else:

                Fe_x = np.zeros(cells.X.shape)
                Fe_y = np.zeros(cells.Y.shape)

            # sum the forces:
            Fx = Fe_x
            Fy = Fe_y

            source_x = -Fx*alpha
            source_x = cells.grid_obj.grid_int(source_x,bounds=btag) # perform finite volume integration of source

            source_y = -Fy*alpha
            source_y = cells.grid_obj.grid_int(source_y,bounds=btag) # perform finite volume integration of source

            # # calculated the fluid flow using the time-independent Stokes Flow equation:
            if p.closed_bound is True:

                # enforce closed conditions on source:
                source_x[:,0] = 0
                source_x[:,-1] = 0
                source_x[0,:] = 0
                source_x[-1,:] = 0

                source_y[:,0] = 0
                source_y[:,-1] = 0
                source_y[0,:] = 0
                source_y[-1,:] = 0

                ux_ecm_o = np.dot(cells.lapENVinv,source_x.ravel())
                uy_ecm_o = np.dot(cells.lapENVinv,source_y.ravel())

            else:
                ux_ecm_o = np.dot(cells.lapENV_P_inv,source_x.ravel())
                uy_ecm_o = np.dot(cells.lapENV_P_inv,source_y.ravel())


            # calculate the divergence of the flow field as the sum of the two spatial derivatives:
            div_uo = fd.divergence(ux_ecm_o.reshape(cells.X.shape),uy_ecm_o.reshape(cells.X.shape),
                cells.delta,cells.delta)

            # perform finite volume integration on the divergence:
            # div_uo = fd.integrator(div_uo)
            div_uo = cells.grid_obj.grid_int(div_uo,bounds=btag)

            # calculate the alpha-scaled internal pressure from the divergence of the force:
            if p.closed_bound is True:

                P = np.dot(cells.lapENV_P_inv, div_uo.ravel())
                P = P.reshape(cells.grid_obj.cents_shape)

                # enforce zero normal gradient boundary conditions on P:
                P[:,0] = P[:,1]
                P[:,-1] = P[:,-2]
                P[0,:] = P[1,:]
                P[-1,:] = P[-2,:]

            else:

                P = np.dot(cells.lapENVinv, div_uo.ravel())
                P = P.reshape(cells.grid_obj.cents_shape)

                # enforce zero pressure boundary conditions on P:
                P[:,0] = 0
                P[:,-1] = 0
                P[0,:] = 0
                P[-1,:] = 0

            # Take the grid gradient of the scaled internal pressure:
            gPx, gPy = fd.gradient(P,cells.delta)

             # subtract the pressure term from the solution to yield a divergence-free flow field
            u_env_x = ux_ecm_o.reshape(cells.X.shape) - gPx
            u_env_y = uy_ecm_o.reshape(cells.X.shape) - gPy

            # velocities at cell centres:
            self.u_env_x = u_env_x[:]
            self.u_env_y = u_env_y[:]

            # reinforce boundary conditions
            if p.closed_bound is True:
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

        #---------------Flow through gap junction connected cells-------------------------------------------------------

        # calculate the inverse viscocity for the cell collection, which is scaled by gj conductivity:
        alpha_gj_o = p.gj_surface*self.gjopen*(1/p.mu_water)
        # average to individual cells:
        alpha_gj = np.dot(cells.gjMatrix,alpha_gj_o)/cells.num_nn

        if p.deform_electro is True:

            Fe_cell_x = self.F_electro_x
            Fe_cell_y = self.F_electro_y

        else:
            Fe_cell_x = np.zeros(len(cells.cell_i))
            Fe_cell_y = np.zeros(len(cells.cell_i))

        # net force is the sum of electrostatic and hydrostatic pressure induced body forces:
        F_net_x = Fe_cell_x + self.F_hydro_x
        F_net_y = Fe_cell_y + self.F_hydro_y

        # integrate body forces:
        F_net_x = cells.integrator(F_net_x)
        F_net_y = cells.integrator(F_net_y)

        # set boundary values (applicable to u_gj_x/y solution):
        F_net_x[cells.bflags_cells] = 0
        F_net_y[cells.bflags_cells] = 0

        # Calculate flow under body forces:
        u_gj_xo = np.dot(cells.lapGJinv,-alpha_gj*F_net_x)
        u_gj_yo = np.dot(cells.lapGJinv,-alpha_gj*F_net_y)

        # calculate divergence of the flow field using general definition:
        # first interpolate flow field at membrane midpoints:
        ux_mem = interp.griddata((cells.cell_centres[:,0],cells.cell_centres[:,1]),u_gj_xo,
                         (cells.mem_mids_flat[:,0],cells.mem_mids_flat[:,1]),fill_value = 0)

        uy_mem = interp.griddata((cells.cell_centres[:,0],cells.cell_centres[:,1]),u_gj_yo,
                                 (cells.mem_mids_flat[:,0],cells.mem_mids_flat[:,1]), fill_value = 0)

        # get the component of the velocity field normal to the membranes:
        u_n = ux_mem*cells.mem_vects_flat[:,2] + uy_mem*cells.mem_vects_flat[:,3]

        # calculate divergence as the sum of this vector x each surface area, divided by cell volume:
        div_u = (np.dot(cells.M_sum_mems, u_n*cells.mem_sa)/cells.cell_vol)

        # calculate the reaction pressure required to counter-balance the flow field:
        P_react = np.dot(cells.lapGJ_P_inv,div_u)

        # calculate its gradient:
        gradP_react = (P_react[cells.cell_nn_i[:,1]] - P_react[cells.cell_nn_i[:,0]])/(cells.nn_len)

        gP_x = gradP_react*cells.cell_nn_tx
        gP_y = gradP_react*cells.cell_nn_ty

        # average the components of the reaction force field at cell centres and get boundary values:
        gPx_cell = np.dot(cells.M_sum_mems,gP_x)/cells.num_mems   # FIXME if this stops working move back to gjMatrix
        gPy_cell = np.dot(cells.M_sum_mems,gP_y)/cells.num_mems

        self.u_cells_x = u_gj_xo - gPx_cell
        self.u_cells_y = u_gj_yo - gPy_cell

        # enforce the boundary conditions:
        self.u_cells_x[cells.bflags_cells] = 0
        self.u_cells_y[cells.bflags_cells] = 0

    def eosmosis(self,cells,p):

        """
        Electroosmosis of ion pumps and channels to potentially create directional fluxes in individual cells.

        This is presently simulated by calculating the Nernst-Planck concentration flux of a weighting
        agent, rho, which moves under its own concentration gradient and
        through the influence of the extracellular voltage gradient and fluid flows tangential to the membrane.

        """

        vmem = self.v_env[cells.map_mem2ecm]

        # components of fluid flow velocity at the membrane:
        if p.fluid_flow is True:
            ux_mem = self.u_env_x.ravel()[cells.map_mem2ecm]
            uy_mem = self.u_env_y.ravel()[cells.map_mem2ecm]

        else:
            ux_mem = 0
            uy_mem = 0

        # get the gradient of rho concentration around each membrane:
        grad_c = np.dot(cells.gradMem,self.rho_pump)
        grad_c_ch = np.dot(cells.gradMem,self.rho_channel)

        # gcx = grad_c*cells.mem_vects_flat[:,4]   # FIXME no concentration gradient will be used for now
        # gcy = grad_c*cells.mem_vects_flat[:,5]
        gcx = 0
        gcy = 0

        # electric field at each membrane

        gvx = - self.E_env_x.ravel()[cells.map_mem2ecm]
        gvy = - self.E_env_y.ravel()[cells.map_mem2ecm]

        # calculate the total Nernst-Planck flux at each membrane:

        fx, fy = nernst_planck_flux(self.rho_pump, gcx, gcy, gvx, gvy,ux_mem,uy_mem,p.D_membrane,p.z_pump,
            self.T,p)

        # component of total flux in direction of membrane
        ftot = fx*cells.mem_vects_flat[:,4] + fy*cells.mem_vects_flat[:,5]

        # calculate the membrane-specific divergence of the flux which will be the "delta rho" at each membrane:
        delta_rho = np.dot(cells.gradMem,ftot)

        #------------------------------------------------
        # calculate the total Nernst-Planck flux at each membrane for rho_channel factor:

        fx2, fy2 = nernst_planck_flux(self.rho_channel, gcx, gcy, gvx, gvy,ux_mem,uy_mem,p.D_membrane,p.z_channel,
            self.T,p)

        # component of total flux in direction of membrane
        ftot2 = fx2*cells.mem_vects_flat[:,4] + fy2*cells.mem_vects_flat[:,5]

        # calculate the membrane-specific divergence of the flux which will be the "delta rho" at each membrane:
        delta_rho2 = np.dot(cells.gradMem,ftot2)

        #------------------------------------------------

        # update the pump and channel density properties:
        self.rho_pump = self.rho_pump + p.dt*delta_rho

        self.rho_channel = self.rho_channel + p.dt*delta_rho2

        # # make sure nothing is non-zero:
        fix_inds = (self.rho_pump < 0).nonzero()
        self.rho_pump[fix_inds] = 0

        fix_inds2 = (self.rho_channel < 0).nonzero()
        self.rho_channel[fix_inds2] = 0

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
            loggers.warning('Oops! Morphogen gated ion channel target not found!')
            ion = []

        return ion

    def initDenv(self,cells,p):

        self.D_env_u = np.zeros((self.D_env.shape[0],cells.grid_obj.u_shape[0],cells.grid_obj.u_shape[1]))
        self.D_env_v = np.zeros((self.D_env.shape[0],cells.grid_obj.v_shape[0],cells.grid_obj.v_shape[1]))

        for i, dmat in enumerate(self.D_env):

            if p.env_type is False: # if air surrounds, first set everything to zero and add in cluster data...
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
            dummyMems[cells.bflags_mems] = self.D_free[i]

            # interp the membrane data to an ecm grid, fill values correspond to environmental diffusion consts:
            if p.env_type is True:
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

            if p.env_type is True:

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

        if p.closed_bound is True:  # set full no slip boundary condition at exterior bounds

            self.D_env_weight_u[:,0] = 0
            self.D_env_weight_u[:,-1] = 0
            self.D_env_weight_u[0,:] = 0
            self.D_env_weight_u[-1,:] = 0

            self.D_env_weight_v[:,0] = 0
            self.D_env_weight_v[:,-1] = 0
            self.D_env_weight_v[0,:] = 0
            self.D_env_weight_v[-1,:] = 0

    def osmotic_P(self,cells,p):

        # initialize osmotic pressures in cells and env

        self.osmo_P_cell = np.zeros(len(self.cc_cells[0]))
        self.osmo_P_env = np.zeros(len(self.cc_env[0]))

        # calculate osmotic pressure in cells based on total molarity:
        for c_ion in self.cc_cells:

            self.osmo_P_cell = c_ion*p.R*self.T + self.osmo_P_cell

        # calculate osmotic pressure in environment based on total molarity:

        for c_ion_env in self.cc_env:

            self.osmo_P_env = c_ion_env*p.R*self.T + self.osmo_P_env


        if p.sim_ECM is False:

            self.osmo_P_delta = self.osmo_P_cell - self.osmo_P_env

        else:
            # smooth out the environmental osmotic pressure:
            self.osmo_P_env = self.osmo_P_env.reshape(cells.X.shape)
            # self.osmo_P_env = fd.integrator(self.osmo_P_env)
            self.osmo_P_env = self.osmo_P_env.ravel()

            self.osmo_P_delta = self.osmo_P_cell[cells.mem_to_cells] - self.osmo_P_env[cells.map_mem2ecm]

            # average the pressure to individual cells
            self.osmo_P_delta = np.dot(cells.M_sum_mems,self.osmo_P_delta)/cells.num_mems

        # Calculate the transmembrane flow of water due to osmotic pressure.
        # High, positive osmotic pressure leads to water flow into the cell. Existing pressure in the cell
        # resists the degree of osmotic influx. The effect also depends on aquaporin fraction in membrane:
        u_osmo_o = (self.osmo_P_delta - self.P_cells)*(p.aquaporins/p.mu_water)

        # map osmotic influx to membranes:
        u_osmo = u_osmo_o[cells.mem_to_cells]

        # calculate the flow due to net mass flux into the cell due to selective active/passive ion transport:
        u_mass_flux = (1/p.rho)*self.get_mass_flux(cells,p)

        u_net = u_osmo + u_mass_flux

        # obtain the divergence of the flow in a timestep, which yields the fractional volume change:
        div_u_osmo = p.dt*np.dot(cells.M_sum_mems,u_net*cells.mem_sa)/cells.cell_vol

        # pressure developing in the cell depends on how much the volume can change:
        P_react = (1 - (1/p.youngMod))*div_u_osmo

        if p.run_sim is True:

            # the inflow of mass adds to base pressure in cells:
            self.P_base = self.P_base + P_react[:]

        else:

            # the inflow of mass adds to base pressure in cells:
            self.P_cells = self.P_cells + P_react[:]

        #------------------------------------------------------------------------------------------------------------

        # actual volume change depends on the mechanical properties (young's modulus) of the
        # tissue:

        self.delta_vol = (1/p.youngMod)*div_u_osmo

        # print("Reaction pressure", P_react)
        # print("volume change", self.delta_vol)
        # print('--------------------')

        # update concentrations and volume in the cell:
        vo = cells.cell_vol[:]

        v1 = (1 + self.delta_vol)*vo

        self.cc_cells =self.cc_cells*(vo/v1)

        if p.voltage_dye is True:

            self.cDye_cell= self.cDye_cell*(vo/v1)

        if p.scheduled_options['IP3'] != 0 or p.Ca_dyn is True:

            self.cIP3 = self.cIP3*(vo/v1)

        cells.cell_vol = v1[:]

    def getHydroF(self,cells,p):
        #----Calculate body forces due to hydrostatic pressure gradients---------------------------------------------

        # determine body force due to hydrostatic pressure gradient between cells:
        gPcells = -(self.P_cells[cells.cell_nn_i[:,1]] - self.P_cells[cells.cell_nn_i[:,0]])/cells.nn_len

        F_hydro_x = gPcells*cells.cell_nn_tx
        F_hydro_y = gPcells*cells.cell_nn_ty

        # calculate a shear electrostatic body force at the cell centre:
        self.F_hydro_x = np.dot(cells.M_sum_mems, F_hydro_x) / cells.num_mems
        self.F_hydro_y = np.dot(cells.M_sum_mems, F_hydro_y) / cells.num_mems

        self.F_hydro = np.sqrt(self.F_hydro_x ** 2 + self.F_hydro_y ** 2)

    def electro_P(self,cells,p):
        """
        Calculates electrostatic pressure in collection of cells and
        an electrostatic body force.

        """

        # # average charge density at the gap junctions between cells:
        # q_cells = 10*self.rho_cells*(1/p.ff_cell)
        # make sure we're using a realistic charge density by working with the referenced cell membrane capacitance:
        if p.sim_ECM is False:

            q_cells = (p.cm*cells.cell_sa*self.vm)/cells.cell_vol

        else:
            q_cells = (p.cm*cells.cell_sa*self.v_cell)/cells.cell_vol

        # charge density interpolated to membranes

        Q_mem = interp.griddata((cells.cell_centres[:,0],cells.cell_centres[:,1]),q_cells,
                         (cells.mem_mids_flat[:,0],cells.mem_mids_flat[:,1]),fill_value = 0)

        # charge density averaged at nearest neighbour membranes:

        # Q_mem = (rho_cells[cells.cell_nn_i[:,0]] + rho_cells[cells.cell_nn_i[:,1]])/2


        if p.sim_ECM is False:

            Eab_o = -(self.vm[cells.cell_nn_i[:,1]] -
                    self.vm[cells.cell_nn_i[:,0]])/(cells.nn_len)


        else:

            Eab_o = -(self.v_cell[cells.cell_nn_i[:,1]] -
                    self.v_cell[cells.cell_nn_i[:,0]])/(cells.nn_len)

        F_x = Q_mem*Eab_o*cells.cell_nn_tx
        F_y = Q_mem*Eab_o*cells.cell_nn_ty

        # calculate a shear electrostatic body force at the cell centre:
        self.F_electro_x = np.dot(cells.M_sum_mems, F_x)/cells.num_mems
        self.F_electro_y = np.dot(cells.M_sum_mems, F_y)/cells.num_mems

        self.F_electro = np.sqrt(self.F_electro_x**2 + self.F_electro_y**2)

        P_x = (self.F_electro_x*cells.cell_vol)/cells.cell_sa
        P_y = (self.F_electro_y*cells.cell_vol)/cells.cell_sa

        self.P_electro = np.sqrt(P_x**2 + P_y**2)

    def gravity_P(self,cells,p):

        # calculate gravity hydrostatic pressure head as P = density x gravity constant x height
        self.P_gravity = np.ones(len(cells.mem_i))*9.81*1000*(cells.mem_mids_flat[:,1] - cells.ymin)

    def ghk_calculator(self,cells,p):
        """
        Uses simulation parameters in the Goldman (GHK) equation
        to calculate an alternative Vmem for validation purposes.

        """

        # begin by initializing all summation arrays for the cell network:
        sum_PmAnion_out = np.zeros(len(cells.cell_i))
        sum_PmAnion_in = np.zeros(len(cells.cell_i))
        sum_PmCation_out = np.zeros(len(cells.cell_i))
        sum_PmCation_in = np.zeros(len(cells.cell_i))

        for i, z in enumerate(self.zs):

            ion_type = np.sign(z)

            if ion_type == -1:

                if p.sim_ECM is True:

                    Dm = np.dot(cells.M_sum_mems,self.Dm_cells[i])/cells.num_mems

                    sum_PmAnion_in = sum_PmAnion_in + Dm*self.cc_cells[i]*(1/p.tm)
                    sum_PmAnion_out = sum_PmAnion_out + Dm*self.cc_env[i][cells.map_cell2ecm]*(1/p.tm)

                else:
                    sum_PmAnion_in = sum_PmAnion_in + self.Dm_cells[i]*self.cc_cells[i]*(1/p.tm)
                    sum_PmAnion_out = sum_PmAnion_out + self.Dm_cells[i]*self.cc_env[i]*(1/p.tm)

            if ion_type == 1:

                if p.sim_ECM is True:

                    Dm = np.dot(cells.M_sum_mems,self.Dm_cells[i])/cells.num_mems

                    sum_PmCation_in = sum_PmCation_in + Dm*self.cc_cells[i]*(1/p.tm)
                    sum_PmCation_out = sum_PmCation_out + Dm*self.cc_env[i][cells.map_cell2ecm]*(1/p.tm)

                else:
                    sum_PmCation_in = sum_PmCation_in + self.Dm_cells[i]*self.cc_cells[i]*(1/p.tm)
                    sum_PmCation_out = sum_PmCation_out + self.Dm_cells[i]*self.cc_env[i]*(1/p.tm)

        self.vm_GHK = ((p.R*self.T)/p.F)*np.log((sum_PmCation_out + sum_PmAnion_in)/(sum_PmCation_in + sum_PmAnion_out))

    def getDeformation(self,cells,t,p):
        """
        Calculates the deformation of the cell cluster under the action
        of intracellular forces and pressures, assuming steady-state
        (slow) changes.

        The method assumes that material is entirely incompressible.
        First, the equation of linear elastic motion is used to calculate
        deformation assuming full compressibility. Then, the divergence
        is calculated, an internal reaction pressure is calculated from
        the divergence, and the gradient of the reaction pressure is subtracted
        from the initial solution to create a divergence-free deformation field.

        """

        # Determine action forces ------------------------------------------------

        # body forces from hydrostatic pressure
        F_hydro_x = self.F_hydro_x
        F_hydro_y = self.F_hydro_y


        # first determine body force components due to electrostatics, if desired:
        if p.deform_electro is True:
            F_electro_x = self.F_electro_x
            F_electro_y = self.F_electro_y

        else:

            F_electro_x = np.zeros(len(cells.cell_i))
            F_electro_y = np.zeros(len(cells.cell_i))

        # Take the total component of pressure from all contributions:
        F_cell_x = F_electro_x + F_hydro_x
        F_cell_y = F_electro_y + F_hydro_y

        # integrate the forces, as is mandated by finite volume methods:
        F_cell_x = cells.integrator(F_cell_x)
        F_cell_y = cells.integrator(F_cell_y)

        #--calculate displacement field for incompressible medium------------------------------------------------

        # calculate the initial displacement field (not divergence free!) for the forces using the linear elasticity
        # equation:

        if p.fixed_cluster_bound is True:

            u_x_o = np.dot(cells.lapGJinv,-(1/p.lame_mu)*(F_cell_x))
            u_y_o = np.dot(cells.lapGJinv,-(1/p.lame_mu)*(F_cell_y))

        else:

            u_x_o = np.dot(cells.lapGJ_P_inv,-(1/p.lame_mu)*(F_cell_x))
            u_y_o = np.dot(cells.lapGJ_P_inv,-(1/p.lame_mu)*(F_cell_y))

         # first interpolate displacement field at membrane midpoints:
        ux_mem = interp.griddata((cells.cell_centres[:,0],cells.cell_centres[:,1]),u_x_o,
                         (cells.mem_mids_flat[:,0],cells.mem_mids_flat[:,1]),fill_value = 0)

        uy_mem = interp.griddata((cells.cell_centres[:,0],cells.cell_centres[:,1]),u_y_o,
                                 (cells.mem_mids_flat[:,0],cells.mem_mids_flat[:,1]), fill_value = 0)

        # get the component of the displacement field normal to the membranes:
        u_n = ux_mem*cells.mem_vects_flat[:,2] + uy_mem*cells.mem_vects_flat[:,3]

        # calculate divergence as the sum of this vector x each surface area, divided by cell volume:
        div_u = (np.dot(cells.M_sum_mems, u_n*cells.mem_sa)/cells.cell_vol)

        # calculate the reaction pressure required to counter-balance the flow field:

        # if p.fixed_cluster_bound is True:

        P_react = np.dot(cells.lapGJ_P_inv,div_u)

        # calculate its gradient:
        gradP_react = (P_react[cells.cell_nn_i[:,1]] - P_react[cells.cell_nn_i[:,0]])/(cells.nn_len)

        gP_x = gradP_react*cells.cell_nn_tx
        gP_y = gradP_react*cells.cell_nn_ty

        # average the components of the reaction force field at cell centres and get boundary values:
        gPx_cell = np.dot(cells.M_sum_mems,gP_x)/cells.num_mems
        gPy_cell = np.dot(cells.M_sum_mems,gP_y)/cells.num_mems

        # calculate the displacement of cell centres under the applied force under incompressible conditions:
        self.d_cells_x = u_x_o - gPx_cell
        self.d_cells_y = u_y_o - gPy_cell

        # # enforce boundary conditions:
        if p.fixed_cluster_bound is True:
            self.d_cells_x[cells.bflags_cells] = 0
            self.d_cells_y[cells.bflags_cells] = 0
            # self.d_cells_x[cells.nn_bound] = 0
            # self.d_cells_y[cells.nn_bound] = 0

    def timeDeform(self,cells,t,p):
        """
        Calculates the deformation of the cell cluster under the action
        of intracellular pressure, considering the full time-dependent
        linear elasticity equation for an incompressible medium.

        The solution method for this equation is similar to the
        steady-state method of deformation(). First the displacement
        field is calculated assuming compressibility,
        a reaction pressure is calculated from the divergence of the
        initial field, and the gradient of the internal pressure is
        subtracted from the initial field to produce a divergence
        free solution.

        This method is working much better than the timeDeform_o()
        so is presently in active use.

        """

        # Check for the adequacy of the time step:
        step_check = (p.dt/(2*p.rc))*np.sqrt(p.lame_mu/1000)

        if step_check > 1.0:

            new_ts = (0.9*2*p.rc)/(np.sqrt(p.lame_mu/1000))

            raise BetseExceptionSimulation(
                    'Time dependent deformation is tricky business, requiring a small time step! '
                    'The time step you are using is too large to bother going further with. '
                    'Please set your time step to ' + str(new_ts) + ' and try again.')

        k_const = (p.dt**2)*(p.lame_mu/1000)


        # # Determine action forces ------------------------------------------------

        # body force from hydrostatic pressure:
        F_hydro_x = self.F_hydro_x
        F_hydro_y = self.F_hydro_y

        # first determine body force components due to electrostatics, if desired:
        if p.deform_electro is True:
            F_electro_x = self.F_electro_x
            F_electro_y = self.F_electro_y

        else:

            F_electro_x = np.zeros(len(cells.cell_i))
            F_electro_y = np.zeros(len(cells.cell_i))


        # Take the total component of pressure from all contributions:
        F_cell_x = F_electro_x + F_hydro_x
        F_cell_y = F_electro_y + F_hydro_y

        # integrate the forces:
        F_cell_x = cells.integrator(F_cell_x)
        F_cell_y = cells.integrator(F_cell_y)

        # # and the solution:
        #
        # self.d_cells_x = cells.integrator(self.d_cells_x)
        # self.d_cells_y = cells.integrator(self.d_cells_y)

        #-------------------------------------------------------------------------------------------------

        self.dx_time.append(self.d_cells_x[:]) # append the solution to the time-save vector
        self.dy_time.append(self.d_cells_y[:])


        # Initial value solution--------------------------------------------------------------------------------
        if t == 0.0:

            wave_speed = np.sqrt(p.lame_mu/1000)
            wave_speed = np.float(wave_speed)
            wave_speed = np.round(wave_speed,2)

            loggers.log_info(
                'Your wave speed is approximately: ' +
                 str(wave_speed) + ' m/s '
            )

            loggers.log_info('Try a world size of at least: ' + str(round((5/3)*(wave_speed/500)*1e6))
                             + ' um for resonance.')

            if p.fixed_cluster_bound is True:

                self.d_cells_x = k_const*np.dot(cells.lapGJ,self.dx_time[-1]) + (k_const/p.lame_mu)*F_cell_x + \
                                 self.dx_time[-1]
                self.d_cells_y = k_const*np.dot(cells.lapGJ,self.dy_time[-1]) + (k_const/p.lame_mu)*F_cell_y + \
                                 self.dy_time[-1]

                # self.d_cells_x[cells.bflags_cells] = 0
                # self.d_cells_y[cells.bflags_cells] = 0
                #
                # self.d_cells_x[cells.nn_bound] = 0
                # self.d_cells_y[cells.nn_bound] = 0

            else:

                self.d_cells_x = k_const*np.dot(cells.lapGJ_P,self.dx_time[-1]) + (k_const/p.lame_mu)*F_cell_x + \
                                 self.dx_time[-1]
                self.d_cells_y = k_const*np.dot(cells.lapGJ_P,self.dy_time[-1]) + (k_const/p.lame_mu)*F_cell_y + \
                                 self.dy_time[-1]



        elif t > 0.0:

            # do the non-initial value, standard solution iteration:

            # calculate the velocity for viscous damping:
            d_ux_dt = (self.dx_time[-1] - self.dx_time[-2])/(p.dt)
            d_uy_dt = (self.dy_time[-1] - self.dy_time[-2])/(p.dt)

            gamma = ((p.dt**2)*(p.mu_tissue*p.lame_mu))/(1000*(2*p.rc))

            if p.fixed_cluster_bound is True:

                self.d_cells_x = k_const*np.dot(cells.lapGJ,self.dx_time[-1]) - gamma*d_ux_dt + \
                           (k_const/p.lame_mu)*F_cell_x + 2*self.dx_time[-1] -self.dx_time[-2]

                self.d_cells_y = k_const*np.dot(cells.lapGJ,self.dy_time[-1]) - gamma*d_uy_dt + \
                           (k_const/p.lame_mu)*F_cell_y + 2*self.dy_time[-1] -self.dy_time[-2]

            else:

                self.d_cells_x = k_const*np.dot(cells.lapGJ_P,self.dx_time[-1]) - gamma*d_ux_dt + \
                           (k_const/p.lame_mu)*F_cell_x + 2*self.dx_time[-1] -self.dx_time[-2]

                self.d_cells_y = k_const*np.dot(cells.lapGJ_P,self.dy_time[-1]) - gamma*d_uy_dt + \
                           (k_const/p.lame_mu)*F_cell_y + 2*self.dy_time[-1] -self.dy_time[-2]


        # calculate divergence of u  -----------------------------------------------------------------------

         # first interpolate displacement field at membrane midpoints:
        ux_mem = interp.griddata((cells.cell_centres[:,0],cells.cell_centres[:,1]),self.d_cells_x,
                         (cells.mem_mids_flat[:,0],cells.mem_mids_flat[:,1]),fill_value = 0)

        uy_mem = interp.griddata((cells.cell_centres[:,0],cells.cell_centres[:,1]),self.d_cells_y,
                                 (cells.mem_mids_flat[:,0],cells.mem_mids_flat[:,1]), fill_value = 0)

        # get the component of the displacement field normal to the membranes:
        u_n = ux_mem*cells.mem_vects_flat[:,2] + uy_mem*cells.mem_vects_flat[:,3]

        # calculate divergence as the sum of this vector x each surface area, divided by cell volume:
        div_u = (np.dot(cells.M_sum_mems, u_n*cells.mem_sa)/cells.cell_vol)

        # calculate the reaction pressure required to counter-balance the flow field:

        P_react = np.dot(cells.lapGJ_P_inv,div_u)

        # self.P_cells = (p.lame_mu/k_const)*P_react[:]

        # calculate its gradient:
        gradP_react = (P_react[cells.cell_nn_i[:,1]] - P_react[cells.cell_nn_i[:,0]])/(cells.nn_len)

        gP_x = gradP_react*cells.cell_nn_tx
        gP_y = gradP_react*cells.cell_nn_ty

        # average the components of the reaction force field at cell centres and get boundary values:
        gPx_cell = np.dot(cells.M_sum_mems,gP_x)/cells.num_mems   # FIXME if this stops working change back to gjMatrix
        gPy_cell = np.dot(cells.M_sum_mems,gP_y)/cells.num_mems

        # calculate the displacement of cell centres under the applied force under incompressible conditions:
        self.d_cells_x = self.d_cells_x - gPx_cell
        self.d_cells_y = self.d_cells_y - gPy_cell

        if p.fixed_cluster_bound is True: # enforce zero displacement boundary condition:

            self.d_cells_x[cells.bflags_cells] = 0
            self.d_cells_y[cells.bflags_cells] = 0

            self.d_cells_x[cells.nn_bound] = 0
            self.d_cells_y[cells.nn_bound] = 0

        # check the displacement for NANs:
        check_v(self.d_cells_x)

    def get_density(self,cells,p):

        mass_water = 1000*cells.cell_vol

        mass_sum = np.zeros(len(cells.cell_i))

        for i, concs in enumerate(self.cc_cells):

            MM = np.float(self.molar_mass[i])
            mass_sum = MM*concs*cells.cell_vol + mass_sum

        total_mass = mass_sum + mass_water
        self.density = total_mass/cells.cell_vol

    def get_mass_flux(self,cells,p):
        """
        Sum up individual trans-membrane and
        trans-gap junction ion fluxes to obtain the
        net flow of mass into a cell.

        Assumes that each ion travels with a hydration
        shell of 6 water molecules.

        """

        # FIXME add in gap junction component to mass flux
         # calculate mass flux across cell membranes:
        mass_flux = np.zeros(len(cells.mem_i))

        for flux_array, mm in zip(self.fluxes_mem,self.molar_mass):

            m_flx = flux_array*(mm + 6*18.01e-3)  # flux x molar mass of ion x 6 water molecules at 18e-3 kg/mol

            mass_flux = mass_flux + m_flx

        return mass_flux

        # # total mass change in cell
        # mass_change = self.mass_flux*p.dt*cells.mem_sa
        # # sum the change over the membranes to get the total mass change of salts:
        # self.delta_m_salts = np.dot(cells.M_sum_mems,mass_change)

    def implement_deform(self,cells,i,p):

        # map individual cell deformations to their membranes. In this case it's better than interpolation
        ux_at_mem = self.d_cells_x[cells.mem_to_cells]
        uy_at_mem = self.d_cells_y[cells.mem_to_cells]

        ux_at_ecm = np.dot(cells.M_sum_mem_to_ecm, ux_at_mem)
        uy_at_ecm = np.dot(cells.M_sum_mem_to_ecm, uy_at_mem)

        # get new ecm verts:
        new_ecm_verts_x = self.ecm_verts_unique_to[:,0] + np.dot(cells.deforM,ux_at_ecm)
        new_ecm_verts_y = self.ecm_verts_unique_to[:,1] + np.dot(cells.deforM,uy_at_ecm)

        ecm_new = np.column_stack((new_ecm_verts_x,new_ecm_verts_y))

        # set the voronoi points originally tagged to the ecm to the value of these new points
        cells.voronoi_grid[cells.map_voronoi2ecm] = ecm_new[:]

        # recreate ecm_verts_unique:
        cells.ecm_verts_unique = ecm_new[:]

        # Repackage ecm verts so that the World module can do its magic:
        ecm_new_flat = ecm_new[cells.ecmInds]  # first expand it to a flattened form (include duplictes)

        # next repackage the structure to include individual cell data
        cells.ecm_verts = [] # null the original ecm verts data structure...

        for i in range(0,len(cells.cell_to_mems)):

            ecm_nest = ecm_new_flat[cells.cell_to_mems[i]]

            ecm_nest = np.asarray(ecm_nest)      # convert region to a numpy array so it can be sorted

            cells.ecm_verts.append(ecm_nest)

        cells.ecm_verts = np.asarray(cells.ecm_verts)   # Voila! Deformed ecm_verts!

        cells.deformWorld(p)

    def implement_deform_timestep(self,cells,t,p):
        # map individual cell deformations to their membranes. In this case it's better than interpolation
        ux_at_mem = self.d_cells_x[cells.mem_to_cells]
        uy_at_mem = self.d_cells_y[cells.mem_to_cells]

        ux_at_ecm = np.dot(cells.M_sum_mem_to_ecm, ux_at_mem)
        uy_at_ecm = np.dot(cells.M_sum_mem_to_ecm, uy_at_mem)

        # get new ecm verts:
        new_ecm_verts_x = self.ecm_verts_unique_to[:,0] + np.dot(cells.deforM,ux_at_ecm)
        new_ecm_verts_y = self.ecm_verts_unique_to[:,1] + np.dot(cells.deforM,uy_at_ecm)

        ecm_new = np.column_stack((new_ecm_verts_x,new_ecm_verts_y))

        # set the voronoi points originally tagged to the ecm to the value of these new points
        cells.voronoi_grid[cells.map_voronoi2ecm] = ecm_new[:]

        # recreate ecm_verts_unique:
        cells.ecm_verts_unique = ecm_new[:]

        # Repackage ecm verts so that the World module can do its magic:
        ecm_new_flat = ecm_new[cells.ecmInds]  # first expand it to a flattened form (include duplictes)

        # next repackage the structure to include individual cell data
        cells.ecm_verts = [] # null the original ecm verts data structure...

        for i in range(0,len(cells.cell_to_mems)):

            ecm_nest = ecm_new_flat[cells.cell_to_mems[i]]

            ecm_nest = np.asarray(ecm_nest)      # convert region to a numpy array so it can be sorted

            cells.ecm_verts.append(ecm_nest)

        cells.ecm_verts = np.asarray(cells.ecm_verts)   # Voila! Deformed ecm_verts!

        cells.deformWorld(p)

        #----------------------------------------
        if p.plot_while_solving is True and t > 0:

            self.checkPlot.resetData(cells,self,p)

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
    # delG = np.absolute(delG_pump)
    # signG = np.sign(delG)

    # alpha = block*p.alpha_NaK*(tb.step(delG_pump,p.halfmax_NaK,p.slope_NaK))
    alpha = block*p.alpha_NaK*(delG_pump - p.halfmax_NaK)

    f_Na  = -alpha*(cNai**p.Na_exp)*(cKo**(p.K_exp))      #flux as [mol/m2s]   scaled to concentrations Na in and K out

    f_K = -(2/3)*f_Na          # flux as [mol/m2s]


    return f_Na, f_K, -f_Na

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
    # delG = np.absolute(delG_pump)
    # signG = np.sign(delG_pump)
    #
    # alpha = p.alpha_Ca*tb.step(delG_pump,p.halfmax_Ca,p.slope_Ca)
    alpha = p.alpha_Ca*delG_pump

    f_Ca  = -alpha*(cCai)      #flux as [mol/s], scaled to concentration in cell

    return f_Ca

def pumpCaER(cCai,cCao,Vm,T,p):
    """
    Pumps calcium out of the cell and into the endoplasmic reticulum.
    Vm is the voltage across the endoplasmic reticulum membrane.

    """

    deltaGATP = 20*p.R*T

    delG_Ca = p.R*T*np.log(cCai/cCao) + 2*p.F*Vm
    delG_CaATP = deltaGATP - (delG_Ca)
    delG_pump = (delG_CaATP/1000)

    alpha = p.alpha_CaER*tb.step(delG_pump,p.halfmax_Ca,p.slope_Ca)

    f_Ca  = alpha*(cCao**(1/2))*(cCai**(1/2))      #flux as [mol/s]

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
    # delG = np.absolute(delG_pump)
    # signG = np.sign(delG)

    # alpha = block*p.alpha_HK*tb.step(delG_pump,p.halfmax_HK,p.slope_HK)
    alpha = block*p.alpha_HK*delG_pump
    f_H  = -alpha*(cHi)*(cKo)      #flux as [mol/s], scaled by concentrations in and out

    f_K = -f_H          # flux as [mol/s]

    return f_H, f_K

def pumpVATP(cHi,cHo,Vm,T,p,block):

    deltaGATP = 20*p.R*T

    delG_H = p.R*T*np.log(cHo/cHi) - p.F*Vm  # free energy to move H+ out of cell

    delG_VATP = deltaGATP - delG_H   # free energy available to move H+ out of cell
    delG_pump = (delG_VATP/1000)
    # delG = np.absolute(delG_pump)
    # signG = np.sign(delG)

    # alpha = block*p.alpha_V*tb.step(delG_pump,p.halfmax_V,p.slope_V)
    alpha = block*p.alpha_V*delG_pump
    f_H  = -alpha*cHi      #flux as [mol/s], scaled by concentrations in and out

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

    if p.sim_ECM is False:
        # v_cell = (self.rho_cells*cells.cell_vol*p.tm)/(p.eo*80*cells.cell_sa)
        v_cell = (1/(4*math.pi*p.eo*80*cells.R*p.ff_cell))*(self.rho_cells*cells.cell_vol)


    else:

        v_cell = (1/(4*math.pi*p.eo*80*cells.R*p.ff_cell))*(self.rho_cells*cells.cell_vol)

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


    self.rho_env = self.rho_env.reshape(cells.X.shape)

    # Perform Finite Volume integration on the environmental charge
    self.rho_env = fd.integrator(self.rho_env)
    # self.rho_env = cells.grid_obj.grid_int(self.rho_env[:], bounds = 'closed')

    rho_env = self.rho_env

    # make sure charge at the global boundary is zero:
    rho_env[:,0] = 0
    rho_env[:,-1] = 0
    rho_env[-1,:] = 0
    rho_env[0,:] = 0

    rho_env = rho_env.ravel()

    fxy = -rho_env/(80*p.ff_env*p.eo)

    # # modify the RHS of the equation to incorporate Dirichlet boundary conditions on external voltage:
    fxy[cells.bBot_k] = (self.bound_V['B']/cells.delta**2)
    fxy[cells.bTop_k] = (self.bound_V['T']/cells.delta**2)
    fxy[cells.bL_k] = (self.bound_V['L']/cells.delta**2)
    fxy[cells.bR_k] = (self.bound_V['R']/cells.delta**2)

    # Solve Poisson's electrostatic equation:
    V = np.dot(cells.lapENVinv,fxy)

    # set the boundary conditions for the outside of the matrix:
    V[cells.bBot_k] = self.bound_V['B']
    V[cells.bTop_k] = self.bound_V['T']
    V[cells.bL_k] = self.bound_V['L']
    V[cells.bR_k] = self.bound_V['R']

    V = V.reshape(cells.X.shape)

    # smooth out the voltage
    V = fd.integrator(V)
    V = V.ravel()

    return V

def get_Vall(self,cells,p):

    self.ff = 1

    rho_stack = np.hstack((self.rho_cells,self.rho_env[cells.map_mem2ecm]))
    vol_stack = np.hstack((cells.cell_vol,cells.ecm_vol[cells.map_mem2ecm]))

    f_stack = (rho_stack*vol_stack)/(p.eo*80*self.ff)

    V_stack = np.dot(f_stack,cells.VMatrix_inv)
    v_cells = V_stack[cells.a:cells.b]
    v_ecms = V_stack[cells.c:cells.d]

    v_env = get_Venv(self, cells, p)

    vm = v_cells[cells.mem_to_cells] - v_env[cells.map_mem2ecm]

    # v_env[cells.map_mem2ecm] = v_ecms

    # print('******')
    # print(v_cells)
    # print('---')
    # print(v_ecms)
    # print('----')
    # print(vm)

    return v_cells, v_env, vm

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
    data          A numpy vector of data points on cell membrane mids
    cells         An instance of the Cells object
    p             An instance of the Parameters object

    Returns
    ----------
    dat_grid      THe data sampled on a uniform grid

    """

    # interpolate vmem defined on mem mids to cell vertices:
    verts_data = np.dot(data,cells.matrixMap2Verts)

    # amalgamate both mem mids and verts data into one stack:
    plot_data = np.hstack((data,verts_data))

    # interpolate the stack to the plotting grid:
    dat_grid = interp.griddata((cells.plot_xy[:,0],cells.plot_xy[:,1]),plot_data,(cells.Xgrid,cells.Ygrid),
                               method=p.interp_type,
                               fill_value=0)
    # # smooth out the data a bit:
    # dat_grid = fd.integrator(dat_grid)

    # get rid of values that bleed into the environment:
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

    # def timeDeform_o(self,cells,t,p):
    #     """
    #     Calculates the deformation of the cell cluster under the action
    #     of intracellular pressure, considering the full time dependent
    #     linear elasticity equation for an incompressible medium.
    #
    #     This version uses the "Earthquake" scientist's theory of S-waves
    #     and the reduction of the linear elasticity equation to solve
    #     for the cross product of the displacement field and force.
    #
    #     Unfortunately, this doesn't seem to produce the right wave
    #     solutions as it's defined here, so it's been taken out of
    #     use for the moment.
    #
    #     """
    #
    #     # Check for the adequacy of the time step:
    #     step_check = (p.dt/(2*p.rc))*np.sqrt(p.lame_mu/1000)
    #
    #     if step_check > 1.0:
    #
    #         new_ts = (0.9*2*p.rc)/(np.sqrt(p.lame_mu/1000))
    #
    #         raise BetseExceptionSimulation(
    #                 'Time dependent deformation is tricky business, requiring a small time step! '
    #                 'The time step you are using is too large to bother going further with. '
    #                 'Please set your time step to ' + str(new_ts) + ' and try again.')
    #
    #     k_const = (p.dt**2)*(p.lame_mu/1000)
    #
    #
    #     # # Determine action forces ------------------------------------------------
    #     # first determine body force components due to electrostatics, if desired:
    #     if p.deform_electro is True:
    #         F_electro_x = self.F_electro_x
    #         F_electro_y = self.F_electro_y
    #
    #     else:
    #
    #         F_electro_x = np.zeros(len(cells.cell_i))
    #         F_electro_y = np.zeros(len(cells.cell_i))
    #
    #
    #     # determine body force due to osmotic water flow, if desired
    #     if p.deform_osmo is True:
    #
    #         F_osmo_x = self.F_osmo_x
    #         F_osmo_y = self.F_osmo_y
    #
    #     else:
    #         F_osmo_x = np.zeros(len(cells.cell_i))
    #         F_osmo_y = np.zeros(len(cells.cell_i))
    #
    #     # Take the total component of pressure from all contributions:
    #     F_cell_x = F_electro_x + F_osmo_x
    #     F_cell_y = F_electro_y + F_osmo_y
    #
    #     self.dx_time.append(self.d_cells_x[:]) # append the solution to the time-save vector
    #     self.dy_time.append(self.d_cells_y[:])
    #
    #     self.phi_time.append(self.phi[:])
    #
    #     # calculate the curl of the force ----------------------------------------------------------------
    #
    #     _, _, curlF = cells.curl(F_cell_x,F_cell_y,0)
    #
    #
    #     # Initial value solution--------------------------------------------------------------------------------
    #     if t == 0.0:
    #
    #         wave_speed = np.sqrt(p.lame_mu/1000)
    #         wave_speed = np.float(wave_speed)
    #         wave_speed = np.round(wave_speed,2)
    #         freq = (wave_speed/p.wsx)*1e-3
    #         freq = np.round(freq,2)
    #
    #         loggers.log_info(
    #             'Your wave speed is about: ' +
    #              str(wave_speed) + ' m/s '
    #             "so with your world size, you'd be looking for resonances at:"
    #             ' ' + str(freq) + ' kHz. '
    #         )
    #
    #         loggers.log_info('Try a world size of at least: ' + str(wave_speed*1000) + ' um for resonance.')
    #
    #
    #         self.phi = k_const*np.dot(cells.lapGJ,self.phi_time[-1]) + (k_const/p.lame_mu)*curlF + self.phi_time[-1]
    #
    #         # set boundary condition for phi:
    #         self.phi[cells.bflags_cells] = 0
    #
    #
    #     elif t > 0.0:
    #
    #
    #         # do the non-initial value, standard solution iteration:
    #
    #         # calculate the velocity for viscous damping:
    #         d_phi_dt = (self.phi - self.phi_time[-2])/(p.dt)
    #         # vy = (self.d_cells_y - self.dy_time[-2])/(p.dt)
    #
    #         gamma = ((p.dt**2)*(p.mu_tissue))/(1000*(2*p.rc))
    #
    #         self.phi = k_const*np.dot(cells.lapGJ,self.phi) - gamma*d_phi_dt + \
    #                    (k_const/p.lame_mu)*curlF + 2*self.phi -self.phi_time[-2]
    #
    #          # set boundary condition for phi:
    #         self.phi[cells.bflags_cells] = 0
    #
    #     # obtain displacement from phi -----------------------------------------------------------------------
    #
    #     curl_phi_x, curl_phi_y, _ = cells.curl(0,0,self.phi)
    #
    #
    #     if p.fixed_cluster_bound is False:
    #
    #         self.d_cells_x = np.dot(cells.lapGJ_P_inv,-curl_phi_x)
    #         self.d_cells_y = np.dot(cells.lapGJ_P_inv,-curl_phi_y)
    #
    #     else:
    #
    #         self.d_cells_x = np.dot(cells.lapGJinv,-curl_phi_x)
    #         self.d_cells_y = np.dot(cells.lapGJinv,-curl_phi_y)
    #
    #
    #     if p.fixed_cluster_bound is True:
    #
    #         self.d_cells_x[cells.bflags_cells] = 0
    #         self.d_cells_y[cells.bflags_cells] = 0
    #
    #     # check the displacement for NANs:
    #     check_v(self.d_cells_x)



# if p.deform_electro is True:
        #
        #     # calculate a divergence-free material flow field
        #     P_mem_net = self.P_mem
        #
        #     F = - (P_mem_net[cells.mem_nn[:,1]] -
        #            P_mem_net[cells.mem_nn[:,0]])/cells.mem_distance
        #
        #     # components of the driving force field:
        #     Fx = F*cells.mem_tx
        #     Fy = F*cells.mem_ty
        #
        #     # divergence of the driving force field:
        #     divF_x = np.dot(cells.M_sum_mems,Fx*cells.mem_sa*cells.mem_vects_flat[:,2])*(1/cells.cell_vol)
        #     divF_y = np.dot(cells.M_sum_mems,Fy*cells.mem_sa*cells.mem_vects_flat[:,3])*(1/cells.cell_vol)
        #
        #     divF = divF_x + divF_y
        #
        #     # reaction pressure:
        #     P_react = cells.num_mems*np.dot(cells.mem_LapM_inv,divF)
        #
        #     P_react_mem = P_react[cells.mem_to_cells]
        #
        #     # reaction force field:
        #     gradP_react = -(P_react_mem[cells.mem_nn[:,1]] - P_react_mem[cells.mem_nn[:,0]])/cells.mem_distance
        #
        #     gPx = gradP_react*cells.mem_tx
        #     gPy = gradP_react*cells.mem_ty
        #
        #     # net body force field:
        #     fx_o = Fx + gPx
        #     fy_o = Fy + gPy
        #
        #     # average the force field to cell centres, multiply by volume to get total force:
        #     fx = (np.dot(cells.M_sum_mems, fx_o)/cells.num_mems)*cells.cell_vol
        #     fy = (np.dot(cells.M_sum_mems, fy_o)/cells.num_mems)*cells.cell_vol
        #
        #     # calculate the stress at the membrane by dividing by membrane surface area:
        #     sx = fx[cells.mem_to_cells]/cells.mem_sa
        #     sy = fy[cells.mem_to_cells]/cells.mem_sa
        #
        #     # Calculate net pressure
        #     P_net_mem = self.P_mem + P_react_mem
        #     self.P_cells = np.dot(cells.M_sum_mems,self.P_mem)/cells.num_mems
        #
        #
        # #--------------------------------------------------------------------------------------------
        # # calculate the strain with respect to the membrane midpoints:
        #
        # # define the strain-stress matrix at the cell membrane midpoints:
        # Y = p.youngMod   # Young's modulus (elastic modulus)
        # poi = 0.5        # Poisson ratio
        #
        # eta_x = (1/Y)*(sx)
        # eta_y = (1/Y)*(sy)
        #
        #
        # if p.fixed_cluster_bound is True:  # if the outer boundary of the cell cluster is prevented from deforming
        #
        #     # set the strain at the boundary to zero:
        #     eta_x[cells.bflags_mems] = 0
        #     eta_y[cells.bflags_mems] = 0
        #
        # # sum the strain from opposing membranes to the shared ecm wall:
        # eta_ecm_x = np.dot(cells.M_sum_mem_to_ecm,eta_x)
        # eta_ecm_y = np.dot(cells.M_sum_mem_to_ecm,eta_y)
        #
        # # calculate the displacement gradient field:
        # d_x = eta_ecm_x*cells.chord_mag
        # d_y = eta_ecm_y*cells.chord_mag
        #
        # #----------------------------------------------------------
        # # Deform key structures:
        #
        # # get new ecm mids:
        # new_ecm_mids_x = cells.ecm_mids[:,0] + d_x
        # new_ecm_mids_y = cells.ecm_mids[:,1] + d_y
        #
        # cells.ecm_mids = np.column_stack((new_ecm_mids_x,new_ecm_mids_y))
        #
        # # get new ecm verts:
        # new_ecm_verts_x = cells.ecm_verts_unique[:,0] + np.dot(cells.deforM,d_x)
        # new_ecm_verts_y = cells.ecm_verts_unique[:,1] + np.dot(cells.deforM,d_y)
        #
        # ecm_new = np.column_stack((new_ecm_verts_x,new_ecm_verts_y))
        #
        # # recreate ecm_verts_unique:
        # cells.ecm_verts_unique = ecm_new[:]
        #
        # if p.sim_ECM is True:
        #
        #     # set the voronoi points originally tagged to the ecm to the value of these new points
        #     cells.voronoi_grid[cells.map_voronoi2ecm] = ecm_new[:]
        #
        #
        # #----------------------------------
        #
        # # Repackage ecm verts so that the Cells module can do its magic:
        #
        # ecm_new_flat = ecm_new[cells.ecmInds]  # first expand it to a flattened form (include duplictes)
        #
        # # next repackage the structure to include individual cell data
        #
        # cells.ecm_verts = [] # null the original ecm verts data structure...
        #
        # for i in range(0,len(cells.cell_to_mems)):
        #
        #     ecm_nest = ecm_new_flat[cells.cell_to_mems[i]]
        #
        #     ecm_nest = np.asarray(ecm_nest)      # convert region to a numpy array so it can be sorted
        #     cent = ecm_nest.mean(axis=0)     # calculate the centre point
        #     angles = np.arctan2(ecm_nest[:,1]-cent[1], ecm_nest[:,0] - cent[0])  # calculate point angles
        #     sorted_region = ecm_nest[np.argsort(angles)]   # sort indices counter-clockwise
        #     sorted_region_b = sorted_region.tolist()
        #
        #     if p.deform_osmo is False:
        #         cells.ecm_verts.append(sorted_region_b)
        #
        #     else:
        #         cells.ecm_verts.append(ecm_nest)
        #
        # cells.ecm_verts = np.asarray(cells.ecm_verts)   # Voila! Deformed ecm_verts!
        #
        # # deform cell world:
        # cells.deformWorld(p)
        #
        # if p.sim_ECM is True:
        #     self.initDenv(cells,p)
        #
        #
        # # if p.deform_osmo is True:
        # # #
        # # #     size_eta_mem = eta_x*cells.mem_vects_flat[:,2] + eta_y*cells.mem_vects_flat[:,3]
        # # #     eta_cell = np.dot(cells.M_sum_mems,size_eta_mem)/cells.num_mems
        # # #     cells.cell_vol = cells.cell_vol*(1+eta_cell)  # increase in cell vol
        # #     self.cc_cells = self.cc_cells*(cell_vol_o/cells.cell_vol)
        #
        #
        # if p.deform_osmo is True:
        #
        #     # with osmosis, we assume cells can deform under the strain. Therefore, update concentrations
        #     # with new volumes:
        #     self.cc_cells = self.cc_cells*(cell_vol_o/cells.cell_vol)
        #
        #
        #     #  if we assume the change in volume has now happened, what is the change in density?
        #     self.density = self.get_density(cells,p)
        #
        #     # # update the pressure at the membrane after expansion:
        #     P_cll = F_net/cells.cell_sa
        #     self.P_mem = P_cll[cells.mem_to_cells]
        #
        #     self.P_cells = np.dot(cells.M_sum_mems,self.P_mem)/cells.num_mems





