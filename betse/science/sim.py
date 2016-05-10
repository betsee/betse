#!/usr/bin/env python3
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

import copy
import os
import os.path
import time
from random import shuffle

import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate as interp
from scipy.ndimage.filters import gaussian_filter

from betse.exceptions import BetseExceptionSimulation
from betse.science import filehandling as fh
from betse.science import finitediff as fd
from betse.science import toolbox as tb
from betse.science.plot.anim.anim import AnimCellsWhileSolving
from betse.science import sim_toolbox as stb
from betse.science.tissue.channels_o import Gap_Junction
from betse.science.tissue.handler import TissueHandler
from betse.util.io.log import logs

class Simulator(object):
    '''
    Contains the main routines used in the simulation of networked cell
    bioelectrical activity. For efficiency, all methods are implemented in
    terms of Numpy-based linear algebra.

    Methods
    -------
    baseInit(cells,p)           Prepares data structures for a cell-only simulation

    baseInit_ECM(cells,p)      Prepares data structures for a simulation with extracellular spaces

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

    Attributes
    ----------
    _anim_cells_while_solving : AnimCellsWhileSolving
        In-place animation of cell voltage as a function of time plotted during
        (rather than after) simulation modelling if both requested and a
        simulation is currently being modelled _or_ `None` otherwise.
    vcell_time : np.ndarray
        Voltage at the inner membrane surface of each cell as a function of
        time.
    venv_time : np.ndarray
        Voltage at the outer membrane surface of each cell as a function of
        time.
    vm : np.ndarray
        Transmembrane voltage of each cell for the current time step.
    vm_time : np.ndarray
        Transmembrane voltage of each cell as a function of time.
    vm_Matrix : np.ndarray
        Transmembrane voltage of each cell as a function of time, resampled for
        use in smooth visualization (e.g., streamplots).
    '''

    def __init__(self, p: 'Parameters'):

        #FIXME: Defer until later. To quote the "simrunner" module, which
        #explicitly calls this public method:
        #   "Reinitialize save and load directories in case params defines new
        #    ones for this sim."
        #Hence, this method should instead be called as the first statement in
        #both the run_loop_no_ecm() and run_loop_with_ecm() methods.
        self.fileInit(p)

    def fileInit(self,p):
        '''
        Initializes the pathnames of top-level files and directories comprising
        the BETSE cache for subsequent initialization and simulation runs.

        This method currently implicitly assigns file names, but will
        (hopefully) permit caller-specified pathnames at some point.
        '''

        # Make the BETSE-specific cache directory if not found.
        betse_cache_dir = os.path.expanduser(p.init_path)
        os.makedirs(betse_cache_dir, exist_ok=True)

        sim_cache_dir = os.path.expanduser(p.sim_path)
        os.makedirs(sim_cache_dir, exist_ok=True)

        # Define data paths for saving an initialization and simulation run:
        self.savedInit = os.path.join(betse_cache_dir, p.init_filename)
        self.savedSim = os.path.join(sim_cache_dir, p.sim_filename)

    def baseInit_all(self, cells, p):
        """
        Creates a host of initialized data matrices for the main simulation,
        including intracellular and environmental concentrations, voltages, and specific
        diffusion constants.

        This method is only done once per world-creation, and therefore contains crucial
        parameters, such as the type of ions included in the simulation, which can't be
        changed after running an initialization.

        """

        i = -1  # dynamic index

        # define sim_ECM and non_sim_ECM data length:
        if p.sim_ECM is True:
            self.mems_data_length = len(cells.mem_i)
            self.env_data_length = len(cells.xypts)
        else:
            self.mems_data_length = len(cells.cell_i)
            self.env_data_length = len(cells.cell_i)

        # load in the object that mathematically handles individual gap junction functionality:
        self.gj_funk = Gap_Junction(p)

        # Identity matrix to easily make matrices out of scalars
        self.id_cells = np.ones(len(cells.cell_i))

        self.cc_cells = []  # cell concentrations initialized
        self.cc_er = []  # endoplasmic reticulum ion concentrations in each cell
        self.cc_env = []  # environmental concentrations initialized

        self.zs = []  # ion valence state initialized
        self.z_er = []  # ion valence states of er ions
        self.z_array = []  # ion valence array matched to cell points

        self.z_array_er = []
        self.Dm_cells = []  # membrane diffusion constants initialized
        self.Dm_er = []  # a list of endoplasmic reticulum membrane state
        self.D_free = []  # a list of single-valued free diffusion constants for each ion
        self.D_gj = []  # an array of diffusion constants for gap junctions
        self.movingIons = []  # moving ions indices
        self.ionlabel = {}  # dictionary to hold ion label names
        self.molar_mass = []

        self.T = p.T  # set the base temperature for the simulation

        self.flx_gj_i = np.zeros(len(cells.nn_i))  # flux matrix across gj for individual ions
        self.fluxes_gj_x = []
        self.fluxes_gj_y = []

        self.I_gj_x = np.zeros(len(cells.nn_i))  # total current in the gj network
        self.I_gj_y = np.zeros(len(cells.nn_i))  # total current in the gj network

        # Membrane current data structure initialization
        self.flx_mem_i = np.zeros(len(cells.mem_i))
        self.fluxes_mem = []
        self.I_mem = np.zeros(len(cells.mem_i))  # total current across membranes

        self.P_cells = np.zeros(len(cells.cell_i))  # initialize pressure in cells

        if p.sim_ECM is True:  # special items specific to simulation of extracellular spaces only:

            # vectors storing separate cell and env voltages
            self.v_env = np.zeros(len(cells.xypts))
            self.v_cell = np.zeros(len(cells.cell_i))

            self.z_array_env = []  # ion valence array matched to env points
            self.D_env = []  # an array of diffusion constants for each ion defined on env grid
            self.c_env_bound = []  # moving ion concentration at global boundary
            self.Dtj_rel = []  # relative diffusion constants for ions across tight junctions

            # Initialize membrane thickness:
            self.tm = np.zeros(len(cells.mem_i))
            self.tm[:] = p.tm

            # initialize environmental fluxes and current data stuctures:
            self.flx_env_i = np.zeros(len(cells.xypts))
            self.fluxes_env_x = []
            self.fluxes_env_y = []
            self.I_env = np.zeros(len(cells.xypts))  # total current in environment

        else:  # items specific to simulation *without* extracellular spaces:
            # Initialize environmental volume:
            self.envV = np.zeros(len(cells.cell_i))
            self.envV[:] = p.vol_env

            # Initialize membrane thickness:
            self.tm = np.zeros(len(cells.cell_i))
            self.tm[:] = p.tm

        if p.fluid_flow is True:
            # Electroosmosis Initialization:

            # initialize vectors for electroosmosis in the cell collection wrt each gap junction (note data type!):
            self.u_cells_x = np.zeros(len(cells.cell_i))
            self.u_cells_y = np.zeros(len(cells.cell_i))

            if p.sim_ECM is True:
                # initialize vectors for env flow (note enhanced data type!):
                self.u_env_x = np.zeros(cells.X.shape)
                self.u_env_y = np.zeros(cells.X.shape)

                # FIXME: place a check here for cells.matrices involved in fluid handling and
                # create it if it doesn't exist

        if p.deformation is True:
            # initialize vectors for potential deformation:
            self.d_cells_x = np.zeros(len(cells.cell_i))
            self.d_cells_y = np.zeros(len(cells.cell_i))

            # FIXME: place a check here for cells.matrices involved in fluid handling and
            # create it if it doesn't exist

        if p.gj_flux_sensitive is True:

            self.gj_rho = np.zeros(len(cells.nn_i))

        else:

            self.gj_rho = 0

        if p.sim_eosmosis is True:  # if simulating electrodiffusive movement of membrane pumps and channels:
            self.rho_pump = np.ones(len(cells.mem_i))
            self.rho_channel = np.ones(len(cells.mem_i))
            self.rho_pump_o = np.ones(len(cells.cell_i))
            self.rho_channel_o = np.ones(len(cells.cell_i))
        else:
            self.rho_pump = 1  # else just define it as identity.
            self.rho_channel = 1
            self.rho_pump_o = 1
            self.rho_channel_o = 1

        ion_names = list(p.ions_dict.keys())

        for name in ion_names:  # go through ion list/dictionary and initialize all sim structures

            if p.ions_dict[name] == 1:

                if name != 'H':

                    i = i+1 # update the dynamic index

                    str1 = 'i' + name  # create the ion index

                    setattr(self, str1, i)  # dynamically add this field to the object

                    self.ionlabel[vars(self)[str1]] = p.ion_long_name[name]

                    if name != 'P':
                        self.movingIons.append(vars(self)[str1])

                    # cell concentration for the ion
                    str_cells = 'c' + name + '_cells'
                    setattr(self, str_cells, np.zeros(len(cells.cell_i)))
                    vars(self)[str_cells][:] = p.cell_concs[name]

                    # environmental concentration for the ion
                    str_env = 'c' + name + '_env'

                    setattr(self, str_env, np.zeros(self.env_data_length))

                    vars(self)[str_env][:] = p.env_concs[name]

                    # base transmembrane diffusion for each ion
                    str_Dm = 'Dm' + name

                    setattr(self, str_Dm, np.zeros(self.mems_data_length))

                    vars(self)[str_Dm][:] = p.mem_perms[name]

                    # base gap junction (intercellular) diffusion for each ion
                    str_Dgj = 'Dgj' + name

                    setattr(self, str_Dgj, np.zeros(len(cells.nn_i)))
                    vars(self)[str_Dgj][:] = p.free_diff[name]

                    # environmental diffusion for each ion:
                    if p.sim_ECM:
                        str_Denv = 'D' + name

                        setattr(self, str_Denv, np.zeros(len(cells.xypts)))
                        vars(self)[str_Denv][:] = p.free_diff[name]

                    # ion charge characteristic for intracellular:
                    str_z = 'z' + name

                    setattr(self, str_z, np.zeros(len(cells.cell_i)))
                    vars(self)[str_z][:] = p.ion_charge[name]

                    if p.sim_ECM:  # ion charge characteristic for extracellular:
                        str_z2 = 'z2' + name

                        setattr(self, str_z2, np.zeros(len(cells.xypts)))
                        vars(self)[str_z2][:] = p.ion_charge[name]

                    self.cc_cells.append(vars(self)[str_cells])
                    self.cc_env.append(vars(self)[str_env])

                    self.zs.append(p.ion_charge[name])
                    self.molar_mass.append(p.molar_mass[name])
                    self.z_array.append(vars(self)[str_z])
                    self.Dm_cells.append(vars(self)[str_Dm])
                    self.D_gj.append(vars(self)[str_Dgj])
                    self.D_free.append(p.free_diff[name])

                    self.fluxes_gj_x.append(self.flx_gj_i)
                    self.fluxes_gj_y.append(self.flx_gj_i)
                    self.fluxes_mem.append(self.flx_mem_i)

                    if p.sim_ECM:
                        self.c_env_bound.append(p.env_concs[name])
                        self.z_array_env.append(vars(self)[str_z2])
                        self.D_env.append(vars(self)[str_Denv])
                        self.Dtj_rel.append(p.Dtj_rel[name])
                        self.fluxes_env_x.append(self.flx_env_i)
                        self.fluxes_env_y.append(self.flx_env_i)

                    if name == 'Ca':
                        self.cCa_er = np.zeros(len(cells.cell_i))
                        self.cCa_er[:] = p.cCa_er

                        self.zCa_er = np.zeros(len(cells.cell_i))
                        self.zCa_er[:] = p.z_Ca

                        self.cc_er.append(self.cCa_er)
                        self.z_er.append(p.z_Ca)
                        self.z_array_er.append(self.zCa_er)

                    if name == 'M' and p.ions_dict['Ca'] == 1:
                        self.cM_er = np.zeros(len(cells.cell_i))
                        self.cM_er[:] = p.cCa_er

                        self.zM_er = np.zeros(len(cells.cell_i))
                        self.zM_er[:] = p.z_M

                        self.cc_er.append(self.cM_er)
                        self.z_er.append(p.z_M)
                        self.z_array_er.append(self.zM_er)

        # Do H+ separately as it's complicated by the bicarbonate buffer

        if p.ions_dict['H'] == 1:

            i = i + 1

            self.iH = i

            self.ionlabel[self.iH] = 'protons'

            self.movingIons.append(self.iH)

            # create concentrations of dissolved carbon dioxide (carbonic acid, non-dissociated):

            self.cHM_cells = np.zeros(len(cells.cell_i))
            self.cHM_cells[:] = 0.03 * p.CO2

            self.cHM_env = np.zeros(self.env_data_length)

            self.cHM_env[:] = 0.03 * p.CO2

            self.cH_cells = np.zeros(len(cells.cell_i))
            self.cH_cells[:] = p.cH_cell

            # use Henderson-Hasselbach equation to obtain pH and cH concentrations:

            self.pH_cell = 6.1 + np.log10(self.cM_cells / self.cHM_cells)
            self.cH_cells = (10 ** (-self.pH_cell)) * 1e3  # units mmol/L

            self.pH_env = 6.1 + np.log10(self.cM_env / self.cHM_env)
            self.cH_env = (10 ** (-self.pH_env)) * 1e3  # units mmol/L

            # initialize diffusion constants

            DmH = np.zeros(self.mems_data_length)

            DmH[:] = p.Dm_H

            self.zH = np.zeros(len(cells.cell_i))
            self.zH[:] = p.z_H

            # gap junction diffusion constant for H+
            DgjH = np.zeros(len(cells.nn_i))
            DgjH[:] = p.free_diff['H']

            if p.sim_ECM is True:
                self.zH2 = np.zeros(len(cells.xypts))
                self.zH2[:] = p.z_H

                # environmental diffusion for H+
                DenvH = np.zeros(len(cells.xypts))
                DenvH[:] = p.free_diff['H']

                # add fixed boundary concentration of H+
                self.c_env_bound.append(p.env_concs['H'])

            # append items to main data vectors:
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

            if p.sim_ECM is True:
                self.z_array_env.append(self.zH2)
                self.D_env.append(DenvH)
                self.Dtj_rel.append(p.Dtj_rel['H'])
                self.fluxes_env_x.append(self.flx_env_i)
                self.fluxes_env_y.append(self.flx_env_i)

        # -------------------------------------------------------------------------------------------------------

        if p.ions_dict['Ca'] == 1:  # initialize the endoplasmic reticulum
            # Define the diffusion matrix for the endoplasmic reticulum:
            self.Dm_er = np.zeros((2, len(cells.cell_i)))
            self.Dm_er[0, :] = p.Dm_Ca
            self.Dm_er[1, :] = p.Dm_M

            self.v_er = np.zeros(len(cells.cell_i))


            # initialize a time-zero vmem vector:

        self.vm_to = np.zeros(self.mems_data_length)  # FIXME is this used anywhere?

        # convert all data structures to Numpy arrays:
        self.cc_cells = np.asarray(self.cc_cells)
        self.cc_env = np.asarray(self.cc_env)

        self.zs = np.asarray(self.zs)
        self.z_array = np.asarray(self.z_array)

        self.Dm_cells = np.asarray(self.Dm_cells)

        self.D_free = np.asarray(self.D_free)
        self.D_gj = np.asarray(self.D_gj)
        self.molar_mass = np.asarray(self.molar_mass)

        self.fluxes_gj_x = np.asarray(self.fluxes_gj_x)
        self.fluxes_gj_y = np.asarray(self.fluxes_gj_y)
        self.fluxes_mem = np.asarray(self.fluxes_mem)

        if p.ions_dict['Ca'] == 1:  # items specific for Calcium dynamics
            self.cc_er = np.asarray(self.cc_er)
            self.z_array_er = np.asarray(self.z_array_er)
            self.Dm_er = np.asarray(self.Dm_er)

        if p.sim_ECM is True:  # items specific for extracellular spaces simulation:

            self.z_array_env = np.asarray(self.z_array_env)
            self.D_env = np.asarray(self.D_env)
            self.fluxes_env_x = np.asarray(self.fluxes_env_x)
            self.fluxes_env_y = np.asarray(self.fluxes_env_y)

            # boundary conditions for voltages:

            # voltage (scheduled dynamics might vary these values)
            self.bound_V = {}
            self.bound_V['T'] = 0
            self.bound_V['B'] = 0
            self.bound_V['L'] = 0
            self.bound_V['R'] = 0

            # initialize the environmental diffusion matrix:
            self.initDenv(cells, p)

        # gap junction specific arrays:
        self.id_gj = np.ones(len(cells.mem_i))  # identity array for gap junction indices...
        self.gjopen = np.ones(len(cells.mem_i))  # holds gap junction open fraction for each gj
        self.gjl = np.zeros(len(cells.mem_i))  # gj length for each gj
        self.gjl[:] = cells.gj_len

    def init_tissue(self, cells: 'Cells', p: 'Parameters') -> None:
        '''
        Prepares data structures pertaining to tissue profiles, dynamic
        activity, and optional methods such as electroosmotic fluid,
        which can be changed in between an initialization and simulation
        run.

        This method is called at the start of all simulations.
        '''

        self.gj_funk = Gap_Junction(p)

        if p.sim_ECM is True:
            #  Initialize diffusion constants for the extracellular transport:
            self.initDenv(cells,p)

        self.dyna = TissueHandler(self, cells, p)   # create the tissue dynamics object
        self.dyna.tissueProfiles(self, cells, p)  # initialize all tissue profiles

        if p.sim_ECM is True:
            # create a copy-base of the environmental junctions diffusion constants:
            self.D_env_base = copy.copy(self.D_env)

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

        # Noise Initialization ---------------------
        # add channel noise to the model:
        self.channel_noise_factor = np.random.random(self.mems_data_length)

        self.Dm_cells[self.iK] = (p.channel_noise_level * self.channel_noise_factor + 1) * self.Dm_cells[self.iK]

        if p.dynamic_noise is True:
            # add a random walk on protein concentration to generate dynamic noise:
            self.protein_noise_factor = p.dynamic_noise_level * (np.random.random(len(cells.cell_i)) - 0.5)

            if p.ions_dict['P'] == 1:
                self.cc_cells[self.iP] = self.cc_cells[self.iP] * (1 + self.protein_noise_factor)

        #--Blocks initialization--------------------

        if p.global_options['gj_block'] != 0:

            self.gj_block = np.ones(len(cells.mem_i))   # initialize the gap junction blocking vector to ones

        else:

            self.gj_block = 1

        # initialize dynamic pump blocking vectors:

        if p.global_options['NaKATP_block'] != 0:

            self.NaKATP_block = np.ones(self.mems_data_length)  # initialize NaKATP blocking vector

        else:
            self.NaKATP_block = 1

        if p.HKATPase_dyn is True and p.global_options['HKATP_block'] != 0:
            self.HKATP_block = np.ones(self.mems_data_length)  # initialize HKATP blocking vector
        else:
            self.HKATP_block = 1

        if p.VATPase_dyn is True and p.global_options['VATP_block'] != 0:
            self.VATP_block = np.ones(self.mems_data_length)  # initialize HKATP blocking vector
        else:
            self.VATP_block = 1

        # -----auxillary molecules initialization -------------------------

        if p.scheduled_options['IP3'] != 0 or p.Ca_dyn is True:

            self.cIP3 = np.zeros(len(cells.cell_i))  # initialize a vector to hold IP3 concentrations
            self.cIP3[:] = p.cIP3_to                 # set the initial concentration of IP3 from params file

            self.IP3_flux_x_gj = np.zeros(len(cells.nn_i))
            self.IP3_flux_y_gj = np.zeros(len(cells.nn_i))
            self.IP3_flux_mem = np.zeros(len(cells.mem_i))

            if p.sim_ECM is True:

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

    def run_sim_core(self, cells: 'Cells', p: 'Parameters') -> None:
        '''
        Runs and saves the current simulation phase (e.g. init or sim phases)

        This method drives the time-loop for the main simulation, including
        gap-junction connections and calls to all dynamic channels and entities (sim phase
        only).

        '''

        self.init_tissue(cells, p)  # Initialize all structures used for gap junctions, ion channels, and other dynamics

        # Reinitialize all time-data structures
        self.clear_storage(cells, p)

        # get the net, unbalanced charge and corresponding voltage in each cell to initialize values of voltages:
        self.update_V(cells, p)

        self.vm_to = np.copy(self.vm)  # create a copy of the original voltage

        if p.Ca_dyn is True:
            self.cc_er_to = np.copy(self.cc_er)

        # Display and/or save an animation during solving and calculate:
        #
        # * "tt", a time-steps vector appropriate for the current run.
        # * "tsamples", that vector resampled to save data at fewer times.
        tt, tsamples = self._plot_loop(cells, p)

        do_once = True  # a variable to time the loop only once

        for t in tt:  # run through the loop

            # start the timer to approximate time for the simulation
            if do_once is True:
                loop_measure = time.time()

            # Reinitialize flux storage device.
            self.fluxes_mem.fill(0)

            # Calculate the change in the voltage derivative.
            self.dvm = (self.vm - self.vm_to) / p.dt

            self.vm_to = np.copy(self.vm)  # reassign the history-saving vm

            if p.Ca_dyn == 1 and p.ions_dict['Ca'] == 1:
                self.dcc_ER = (self.cc_er - self.cc_er_to) / p.dt

                self.cc_er_to = np.copy(self.cc_er)

            # Calculate the values of scheduled and dynamic quantities (e.g.
            # ion channel multipliers).
            if p.run_sim is True:
                self.dyna.runAllDynamics(self, cells, p, t)

            # -----------------PUMPS-------------------------------------------------------------------------------------

            if p.sim_ECM is True:
                # run the Na-K-ATPase pump:
                fNa_NaK, fK_NaK, self.rate_NaKATP = stb.pumpNaKATP(
                    self.cc_cells[self.iNa][cells.mem_to_cells],
                    self.cc_env[self.iNa][cells.map_mem2ecm],
                    self.cc_cells[self.iK][cells.mem_to_cells],
                    self.cc_env[self.iK][cells.map_mem2ecm],
                    self.vm,
                    self.T,
                    p,
                    self.NaKATP_block,
                )

            else:

                fNa_NaK, fK_NaK, self.rate_NaKATP = stb.pumpNaKATP(
                            self.cc_cells[self.iNa],
                            self.cc_env[self.iNa],
                            self.cc_cells[self.iK],
                            self.cc_env[self.iK],
                            self.vm,
                            self.T,
                            p,
                            self.NaKATP_block,
                        )

            if p.sim_ECM is True:
                # modify the fluxes by electrodiffusive membrane redistribution factor and add fluxes to storage:
                self.fluxes_mem[self.iNa] = self.rho_pump * fNa_NaK
                self.fluxes_mem[self.iK] = self.rho_pump * fK_NaK

            else:

                self.fluxes_mem[self.iNa] = self.rho_pump*fNa_NaK[cells.mem_to_cells]
                self.fluxes_mem[self.iK] = self.rho_pump*fK_NaK[cells.mem_to_cells]

            # update the concentrations of Na and K in cells and environment:
            self.update_C(self.iNa, fNa_NaK, cells, p)
            self.update_C(self.iK, fK_NaK, cells, p)

            # recalculate the net, unbalanced charge and voltage in each cell:
            self.update_V(cells, p)

            # -------------------------------------- FIXME shift these to after main comps

            if p.ions_dict['Ca'] == 1:  # FIXME do in Ca handler!

                self.ca_handler(cells,p)



            if p.ions_dict['H'] == 1:

                self.acid_handler(cells,p)

            # ----------------ELECTRODIFFUSION---------------------------------------------------------------------------

            # electro-diffuse all ions (except for proteins, which don't move) across the cell membrane:

            shuffle(self.movingIons)

            for i in self.movingIons:

                if p.sim_eosmosis is True and i == self.iK:  # confine electroosmotic movement to K+ channels:

                    rho_i = self.rho_channel

                else:
                    rho_i = 1

                if p.sim_ECM is True:

                    f_ED = stb.electroflux(self.cc_env[i][cells.map_mem2ecm], self.cc_cells[i][cells.mem_to_cells],
                        self.Dm_cells[i], self.tm, self.zs[i], self.vm, self.T, p,
                        rho=rho_i)

                else:

                    f_ED = stb.electroflux(self.cc_env[i],self.cc_cells[i],self.Dm_cells[i],self.tm,self.zs[i],
                                    self.vm,self.T,p,rho=rho_i)


                if p.sim_ECM is True:
                    self.fluxes_mem[i] = self.fluxes_mem[i] + f_ED

                else:
                    self.fluxes_mem[i] = self.fluxes_mem[i] + f_ED[cells.mem_to_cells]

                # update ion concentrations in cell and ecm:
                self.update_C(i, f_ED, cells, p)

                # update flux between cells due to gap junctions
                self.update_gj(cells, p, t, i)

                if p.sim_ECM:
                    # update concentrations in the extracellular spaces:
                    self.update_ecm(cells, p, t, i)

                # # recalculate the net, unbalanced charge and voltage in each cell:
                self.update_V(cells, p)

            self.get_Efield(cells, p)
            # get forces from any hydrostatic (self.P_Cells) pressure:
            self.getHydroF(cells, p)

            # calculate pressures:

            if p.deform_osmo is True:
                self.osmotic_P(cells, p)

            if p.deform_electro is True:
                self.electro_P(cells, p)

            if p.fluid_flow is True and p.run_sim is True:  # FIXME fluid flow isn't plotting!

                self.run_sim = True

                self.getFlow(cells, p)

            # if desired, electroosmosis of membrane channels
            if p.sim_eosmosis is True and p.run_sim is True:
                self.run_sim = True

                self.eosmosis(cells, p)  # modify membrane pump and channel density according to Nernst-Planck

            if p.deformation is True and p.run_sim is True:

                self.run_sim = True

                if p.td_deform is False:

                    self.getDeformation(cells, t, p)

                elif p.td_deform is True:

                    self.timeDeform(cells, t, p)

            if p.scheduled_options['IP3'] != 0 or p.Ca_dyn is True:
                self.update_IP3(cells, p, t)

            if p.Ca_dyn == 1 and p.ions_dict['Ca'] == 1:
                self.update_er(cells, p, t)

            # if p.voltage_dye=1 electrodiffuse voltage sensitive dye between cell and environment
            if p.voltage_dye == 1:
                self.update_dye(cells, p, t)

            if p.dynamic_noise == 1 and p.ions_dict['P'] == 1:
                # add a random walk on protein concentration to generate dynamic noise:
                self.protein_noise_factor = p.dynamic_noise_level * (np.random.random(len(cells.cell_i)) - 0.5)
                self.cc_cells[self.iP] = self.cc_cells[self.iP] * (1 + self.protein_noise_factor)

                # recalculate the net, unbalanced charge and voltage in each cell:
                self.update_V(cells, p)

            stb.check_v(self.vm)

            # ---------time sampling and data storage---------------------------------------------------

            if t in tsamples:
                self.get_current(cells, p)  # get the current in the gj network connection of cells

                self.write2storage(t, cells, p)  # write data to time storage vectors

                # Update the currently displayed and/or saved animation.
                self._replot_loop(p)

            # Get time for loop and estimate total time for simulation.
            if do_once is True:
                loop_time = time.time() - loop_measure

                if p.run_sim is True:
                    time_estimate = round(loop_time * p.sim_tsteps, 2)
                else:
                    time_estimate = round(loop_time * p.init_tsteps, 2)
                logs.log_info("This run should take approximately " + str(time_estimate) + ' s to compute...')
                do_once = False

        # Find embedded functions that can't be pickled...
        fh.safe_pickle(self, p)

        cells.points_tree = None

        # Explicitly close the prior animation to conserve memory.
        self._deplot_loop()

        # save the init or sim and report results of potential interest to the user.
        self.save_and_report(cells, p)

        plt.close()
        logs.log_info('Simulation completed successfully.')


    #.................{  INITIALIZERS & FINALIZERS  }............................................
    def clear_storage(self, cells, p):
        """
        Re-initializes time storage vectors at the begining of a sim or init.

        """

        # clear mass flux storage vectors:
        self.fluxes_gj_x  = np.zeros(self.fluxes_gj_x.shape)
        self.fluxes_gj_y = np.zeros(self.fluxes_gj_y.shape)
        self.fluxes_mem = np.zeros(self.fluxes_mem.shape)

        self.cc_time = []  # data array holding the concentrations at time points
        self.cc_env_time = [] # data array holding environmental concentrations at time points

        self.dd_time = []  # data array holding membrane permeabilites at time points
        self.vm_time = []  # data array holding voltage at time points
        self.vm_GHK_time = [] # data array holding GHK vm estimates
        self.dvm_time = []  # data array holding derivative of voltage at time points
        self.time = []     # time values of the simulation
        self.gjopen_time = []   # stores the fractional gap junction open state at each time
        self.osmo_P_delta_time = []  # osmotic pressure difference between cell interior and exterior as func of time

        self.I_mem_time = []    # initialize membrane current time vector

        self.vm_Matrix = []    # initialize matrices for resampled data sets (used in smooth plotting and streamlines)
        self.I_gj_x_time = []
        self.I_gj_y_time = []
        self.I_tot_x_time = []
        self.I_tot_y_time = []

        self.efield_gj_x_time = []   # matrices storing smooth electric field in gj connected cells
        self.efield_gj_y_time = []

        self.P_cells_time = []
        self.u_cells_x_time = []
        self.u_cells_y_time = []

        self.rho_cells_time = []

        self.F_electro_time = []
        self.F_electro_x_time = []
        self.F_electro_y_time = []
        self.F_hydro_x_time = []
        self.F_hydro_y_time = []
        self.F_hydro_time = []
        self.P_electro_time = []

        self.rate_NaKATP_time =[]

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

        if p.sim_eosmosis is True:
            self.rho_channel_time = []
            self.rho_pump_time = []

        if p.Ca_dyn is True:
            self.cc_er_time = []
            self.cc_er_to = np.copy(self.cc_er[:])

        if p.sim_ECM is True:

            # clear flux storage vectors for environment
            self.fluxes_env_x = np.zeros(self.fluxes_env_x.shape)
            self.fluxes_env_y = np.zeros(self.fluxes_env_y.shape)


            self.vcell_time = []
            self.venv_time = []

            self.cc_er_time = []   # retains er concentrations as a function of time
            self.cIP3_time = []    # retains cellular ip3 concentrations as a function of time
            self.cIP3_env_time = []

            self.efield_ecm_x_time = []   # matrices storing smooth electric field in ecm
            self.efield_ecm_y_time = []

            # initialize time-storage vectors for electroosmotic data:
            self.u_env_x_time = []
            self.u_env_y_time = []

            self.rho_pump_time = []    # store pump and channel states as function of time...
            self.rho_channel_time = []

            vm_dato = np.zeros(len(cells.mem_i))
            dat_grid_vm = stb.vertData(vm_dato,cells,p)
            self.vm_Matrix.append(dat_grid_vm[:])

            self.cDye_env_time = []

    def write2storage(self,t,cells,p):

        if p.GHK_calc is True:
                stb.ghk_calculator(self,cells,p)
                self.vm_GHK_time.append(self.vm_GHK) # data array holding GHK vm estimates

        # add the new concentration and voltage data to the time-storage matrices:
        self.efield_gj_x_time.append(self.E_gj_x[:])

        self.efield_gj_y_time.append(self.E_gj_y[:])

        concs = np.copy(self.cc_cells[:])
        self.cc_time.append(concs)
        concs = None

        envsc = np.copy(self.cc_env[:])
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
        self.dvm_time.append(self.dvm[:])
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
            self.implement_deform_timestep(cells, t, p)

            #FIXME: Shift into implement_deform_timestep(). Magical OK!
            self.cell_centres_time.append(cells.cell_centres[:])
            self.mem_mids_time.append(cells.mem_mids_flat[:])
            self.maskM_time.append(cells.maskM[:])
            self.mem_edges_time.append(cells.mem_edges_flat[:])
            self.cell_verts_time.append(cells.cell_verts[:])

            self.dx_cell_time.append(self.d_cells_x[:])
            self.dy_cell_time.append(self.d_cells_y[:])

        if p.fluid_flow is True and p.run_sim is True:
            self.u_cells_x_time.append(self.u_cells_x[:])
            self.u_cells_y_time.append(self.u_cells_y[:])

        if p.sim_eosmosis is True and p.run_sim is True:
            self.rho_channel_time.append(self.rho_channel[:])
            self.rho_pump_time.append(self.rho_pump[:])

        self.gjopen_time.append(self.gjopen[:])
        self.time.append(t)

        if p.scheduled_options['IP3'] != 0 or p.Ca_dyn is True:
            self.cIP3_time.append(self.cIP3[:])

        if p.voltage_dye ==1:
            self.cDye_time.append(self.cDye_cell[:])

        if p.Ca_dyn == 1 and p.ions_dict['Ca']==1:
            self.cc_er_time.append(np.copy(self.cc_er[:]))

        if p.sim_ECM is True:

            self.efield_ecm_x_time.append(self.E_env_x[:])

            self.efield_ecm_y_time.append(self.E_env_y[:])

            # ecmsc = np.copy(self.cc_env[:])
            # ecmsc.tolist()
            # self.cc_env_time.append(ecmsc)
            # ecmsc = None

            self.vcell_time.append(self.v_cell[:])
            self.venv_time.append(self.v_env[:])

            if p.fluid_flow is True and p.run_sim is True:
                self.u_env_x_time.append(self.u_env_x[:])
                self.u_env_y_time.append(self.u_env_y[:])

            # calculate interpolated verts and midpoint data for Vmem:
            dat_grid_vm = stb.vertData(self.vm[:],cells,p)

            self.vm_Matrix.append(dat_grid_vm[:])

            if p.voltage_dye ==1:
                self.cDye_env_time.append(self.cDye_env[:])

    def save_and_report(self,cells,p):

        # save the init or sim:
        if p.run_sim is False:
            datadump = [self, cells, p]
            fh.saveSim(self.savedInit, datadump)
            message_1 = 'Initialization run saved to' + ' ' + p.init_path
            logs.log_info(message_1)
        else:
            datadump = [self, cells, p]
            fh.saveSim(self.savedSim, datadump)
            message_2 = 'Simulation run saved to' + ' ' + p.sim_path
            logs.log_info(message_2)

        # report final output to user:
        for i in range(0,len(self.ionlabel)):
            endconc = np.round(np.mean(self.cc_time[-1][i]),6)
            label = self.ionlabel[i]
            concmess = 'Final average cytoplasmic concentration of'+ ' '+ label + ': '
            logs.log_info(concmess + str(endconc) + ' mmol/L')

        for i in range(0,len(self.ionlabel)):
            endconc = np.round(np.mean(self.cc_env_time[-1][i]),6)
            label = self.ionlabel[i]
            concmess = 'Final environmental concentration of'+ ' '+ label + ': '
            logs.log_info(concmess + str(endconc) + ' mmol/L')

        final_vmean = 1000*np.round(np.mean(self.vm_time[-1]),6)
        vmess = 'Final average cell Vmem of ' + ': '
        logs.log_info(vmess + str(final_vmean) + ' mV')

        if p.GHK_calc is True:
            final_vmean_GHK = 1000*np.round(np.mean(self.vm_GHK_time[-1]),6)
            vmess = 'Final average cell Vmem calculated using GHK: ' + ': '
            logs.log_info(vmess + str(final_vmean_GHK) + ' mV')

        if p.ions_dict['H'] == 1:
            final_pH = -np.log10(1.0e-3*np.mean((self.cc_time[-1][self.iH])))
            logs.log_info('Final average cell pH ' + str(np.round(final_pH, 2)))

            final_pH_env = -np.log10(np.mean(1.0e-3*(self.cc_env_time[-1][self.iH])))
            logs.log_info('Final environmental pH ' + str(np.round(final_pH_env, 2)))


        if p.scheduled_options['IP3'] != 0 or p.Ca_dyn is True:
            IP3_env_final = np.mean(self.cIP3_env)
            IP3_cell_final = np.mean(self.cIP3)
            logs.log_info('Final IP3 concentration in the environment: ' + str(np.round(IP3_env_final, 6)) + ' mmol/L')
            logs.log_info('Final average IP3 concentration in cells: ' + str(np.round(IP3_cell_final, 6)) + ' mmol/L')

        if p.Ca_dyn == 1 and p.ions_dict['Ca'] == 1:

            endconc_er = np.round(np.mean(self.cc_er[0]),6)
            label = self.ionlabel[self.iCa]
            concmess = 'Final average ER concentration of'+ ' '+ label + ': '
            logs.log_info(concmess + str(endconc_er) + ' mmol/L')

        if p.voltage_dye ==1:
            dye_env_final = np.mean(self.cDye_env)
            dye_cell_final = np.mean(self.cDye_cell)
            logs.log_info('Final average morphogen concentration in the environment: ' + str(np.round(dye_env_final, 6))
                          + ' mmol/L')
            logs.log_info('Final average morphogen concentration in cells: ' + str(np.round(dye_cell_final, 6)) +
                             ' mmol/L')

    def sim_info_report(self,cells,p):

        logs.log_info('This world contains ' + str(cells.cell_number) + ' cells.')
        logs.log_info('Each cell has an average of ' + str(round(cells.average_nn, 2)) + ' nearest-neighbours.')
        logs.log_info('You are running the ion profile: ' + p.ion_profile)

        logs.log_info('Ions in this simulation: ' + str(self.ionlabel))
        logs.log_info(
            'If you have selected features using other ions, '
            'they will be ignored.')

        logs.log_info('Considering extracellular spaces: ' + str(p.sim_ECM))
        logs.log_info('Electroosmotic fluid flow: ' + str(p.fluid_flow))
        logs.log_info('Ion pump and channel electodiffusion in membrane: ' + str(p.sim_eosmosis))
        logs.log_info('Force-induced cell deformation: ' + str(p.deformation))
        logs.log_info('Osmotic pressure: ' + str(p.deform_osmo))
        logs.log_info('Electrostatic pressure: ' + str(p.deform_electro))


    # ................{  DOERS & GETTERS }...............................................
    def get_Vall(self, cells, p) -> (np.ndarray, np.ndarray, np.ndarray):
        """
        Calculates transmembrane voltage (Vmem) and voltages inside the cell
        and extracellular space, assuming the cell is a capacitor consisting of
        two concentric spheres.

        Parameters
        --------
        cells:
        p:

        Returns
        --------
        vm : np.ndarray
            Transmembrane voltage.
        v_cell : np.ndarray
            Voltage at the inner membrane surface.
        v_env : np.ndarray
            Voltage at the outer membrane surface.
        """

        if p.sim_ECM is False:
            # if we're not modelling extracellular spaces, but assuming a perfectly mixing environment,
            # assume the cell is a concentric spherical capacitor with opposite but equal magnitude
            # charge inside the cell and out
            # calculate the voltage for the system using the measured capacitance of the cell membrane:

            vm = (1/(p.cm*cells.cell_sa))*self.rho_cells*cells.cell_vol

            # in this case, the voltage in the environment is assumed to be zero:
            v_cell = vm
            v_env = 0

        else:
            # total charge in cells per unit surface area:
            Qcells = (self.rho_cells*cells.cell_vol)/cells.cell_sa

            # smooth out the environmental charge:
            # self.rho_env = gaussian_filter(self.rho_env.reshape(cells.X.shape),2)
            # self.rho_env = fd.integrator(self.rho_env.reshape(cells.X.shape))
            # self.rho_env = self.rho_env.ravel()

            # interpolate charge from environmental grid to the ecm_mids:
            rho_ecm = interp.griddata((cells.xypts[:,0],cells.xypts[:,1]),
                                      self.rho_env, (cells.ecm_mids[:,0], cells.ecm_mids[:,1]), method='nearest',
                                      fill_value = 0)

                # total charge per unit surface area in the extracellular spaces:
            Qecm = rho_ecm*(p.cell_space/2)

            # concatenate the cell and ecm charge vectors to the maxwell capacitance vector:
            Q_max_vect = np.hstack((Qcells,Qecm))

           # original solver in terms of pseudo-inverse matrix:
            v_max_vect = np.dot(cells.M_max_cap_inv, Q_max_vect)

            # separate voltages for cells and ecm spaces
            v_cell = v_max_vect[cells.cell_range_a:cells.cell_range_b]
            v_ecm = v_max_vect[cells.ecm_range_a:cells.ecm_range_b]

            #Map the environmental voltage to the regular grid:
            v_env = np.zeros(len(cells.xypts))

            # map the ecm voltage to membrane midpoints:
            v_ecm_at_mem = v_ecm[cells.mem_to_ecm_mids]
            # map it again to the environmental grid
            v_env[cells.map_mem2ecm] = v_ecm_at_mem

            # smooth out the environmental voltage:
            v_env = gaussian_filter(v_env.reshape(cells.X.shape),2)
            # v_env = fd.integrator(v_env.reshape(cells.X.shape))
            v_env = v_env.ravel()

            # # set the conditions for the global boundaries:
            v_env[cells.bBot_k] = self.bound_V['B']
            v_env[cells.bTop_k] = self.bound_V['T']
            v_env[cells.bL_k] = self.bound_V['L']
            v_env[cells.bR_k] = self.bound_V['R']

            # calculate the vm
            vm = v_cell[cells.mem_to_cells] - v_env[cells.map_mem2ecm]

            # null out charge and voltage in the environmental space:
            self.rho_env[cells.inds_env] = 0
            v_env[cells.inds_env] = 0

        return vm, v_cell, v_env

    def update_V(self,cells,p):

        if p.sim_ECM is True:
            # get the charge in cells and the environment:
            self.rho_cells = stb.get_charge_density(self.cc_cells, self.z_array, p)
            self.rho_env = stb.get_charge_density(self.cc_env, self.z_array_env, p)
            self.vm, self.v_cell, self.v_env = self.get_Vall(cells,p)
        else:
             self.rho_cells = stb.get_charge_density(self.cc_cells, self.z_array, p)
             self.vm, _, _ = self.get_Vall(cells,p)

    def update_C(self,ion_i,flux,cells,p):  # FIXME this is a perfect place to implement RK4

        c_cells = self.cc_cells[ion_i][:]
        c_env = self.cc_env[ion_i][:]

        if p.sim_ECM is True:

            d_c_cells = flux*(cells.mem_sa/cells.cell_vol[cells.mem_to_cells])
            d_c_env = -flux*(cells.mem_sa/cells.ecm_vol[cells.map_mem2ecm])

            delta_cells =  np.dot(d_c_cells, cells.cell_UpdateMatrix)
            delta_env = np.dot(d_c_env, cells.ecm_UpdateMatrix)

            self.cc_cells[ion_i] = c_cells + delta_cells*p.dt
            self.cc_env[ion_i] = c_env + delta_env*p.dt

        else:

            delta_cells = flux* (cells.cell_sa/cells.cell_vol)
            delta_env = -flux*(cells.cell_sa/p.vol_env)

            self.cc_cells[ion_i] = c_cells + delta_cells * p.dt

            c_env = c_env + delta_env * p.dt
            # assume auto-mixing of environmental concentrations:
            self.cc_env[ion_i] = c_env.mean()

        # ensure that there are no negative values in the cells or the extracellular spaces:
        for i, arr in enumerate(self.cc_cells):
            self.cc_cells[i] = stb.no_negs(arr)
        for i, arr in enumerate(self.cc_env):
            self.cc_env[i] = stb.no_negs(arr)

    def acid_handler(self,cells,p): # FIXME update for sim_ECM = True *&* False

        # electrofuse the H+ ion between the cytoplasm and the environment
        if p.sim_ECM is True:

            # Electrofuse the H+ ion between the cytoplasm and the ecms.
            f_H1 = stb.electroflux(
                self.cc_env[self.iH][cells.map_mem2ecm],
                self.cc_cells[self.iH][cells.mem_to_cells],
                self.Dm_cells[self.iH],
                self.tm[cells.mem_to_cells],
                self.zs[self.iH],
                self.vm,
                self.T,
                p,
            )

            self.fluxes_mem[self.iH] = self.fluxes_mem[self.iH] + f_H1


        else:

            f_H1 = stb.electroflux(self.cc_env[self.iH],self.cc_cells[self.iH],self.Dm_cells[self.iH],self.tm,
                self.zs[self.iH],self.vm,self.T,p)

            self.fluxes_mem[self.iH] = f_H1[cells.mem_to_cells]

        # update H+ in cells and environment, first in absence of bicarbonate buffering:
        self.update_C(self.iH, f_H1, cells, p)


        # run the bicarbonate buffer to ensure realistic concentrations and pH in cell and environment:
        self.cc_cells[self.iH], _, self.cc_cells[self.iM], self.pH_cell = stb.bicarbonate_buffer(
                                                                                        self.cc_cells[self.iH],
                                                                                        self.cHM_cells,
                                                                                        self.cc_cells[self.iM],
                                                                                        p)

        self.cc_env[self.iH], _, self.cc_env[self.iM], self.pH_env = stb.bicarbonate_buffer(
                                                                                            self.cc_env[self.iH],
                                                                                            self.cHM_env,
                                                                                            self.cc_env[self.iM],
                                                                                            p)

        # recalculate the net, unbalanced charge and voltage in each cell:
        self.update_V(cells,p)

        if p.HKATPase_dyn == 1 and p.run_sim is True: # if there's an H,K ATPase pump

            if p.sim_ECM is True:

                f_H2, f_K2 = stb.pumpHKATP(self.cc_cells[self.iH][cells.mem_to_cells],
                    self.cc_env[self.iH][cells.map_mem2ecm],
                    self.cc_cells[self.iK][cells.mem_to_cells], self.cc_env[self.iK][cells.map_mem2ecm],
                    self.vm, self.T, p, self.HKATP_block)

                # modify fluxes by any uneven redistribution of pump location:
                f_H2 = self.rho_pump * f_H2
                f_K2 = self.rho_pump * f_K2

                self.HKATPase_rate = f_H2[:]

                self.fluxes_mem[self.iH] = self.fluxes_mem[self.iH] + f_H2
                self.fluxes_mem[self.iK] = self.fluxes_mem[self.iK] + f_K2


            else:

                # if HKATPase pump is desired, run the H-K-ATPase pump:
                f_H2, f_K2 = stb.pumpHKATP(self.cc_cells[self.iH],self.cc_env[self.iH],self.cc_cells[self.iK],
                    self.cc_env[self.iK],self.vm,self.T,p,self.HKATP_block)

                # store fluxes for this pump:
                self.fluxes_mem[self.iH] = self.fluxes_mem[self.iH] + self.rho_pump * f_H2[cells.mem_to_cells]
                self.fluxes_mem[self.iK] = self.fluxes_mem[self.iK] + self.rho_pump * f_K2[cells.mem_to_cells]

            # update the concentration in cells (assume environment steady and constant supply of ions)
            # update bicarbonate instead of H+, assuming buffer action holds:

            # calculate the update to K+ in the cell and environment:
            self.update_C(self.iK, f_K2, cells, p)

            # Update the anion (bicarbonate) concentration instead of H+, assuming bicarb buffer holds:
            self.update_C(self.iM, -f_H2, cells, p)


            # Calculate the new pH and H+ concentrations:
            # run the bicarbonate buffer to ensure realistic concentrations and pH in cell and environment:
            self.cc_cells[self.iH], _, self.cc_cells[
                self.iM], self.pH_cell = stb.bicarbonate_buffer(
                self.cc_cells[self.iH],
                self.cHM_cells,
                self.cc_cells[self.iM],
                p)

            self.cc_env[self.iH], _, self.cc_env[self.iM], self.pH_env = stb.bicarbonate_buffer(
                self.cc_env[self.iH],
                self.cHM_env,
                self.cc_env[self.iM],
                p)

            # recalculate the net, unbalanced charge and voltage in each cell:
            self.update_V(cells,p)

        if p.VATPase_dyn == 1 and p.run_sim is True:  # if there's a V-ATPase pump

            if p.sim_ECM is True:

                # if HKATPase pump is desired, run the H-K-ATPase pump:
                f_H3 = stb.pumpVATP(self.cc_cells[self.iH][cells.mem_to_cells], self.cc_env[self.iH][cells.map_mem2ecm],
                    self.vm, self.T, p, self.VATP_block)

                # modify flux by any uneven redistribution of pump location:
                f_H3 = self.rho_pump * f_H3

                self.fluxes_mem[self.iH] = self.fluxes_mem[self.iH] + f_H3

            else:

                # if HKATPase pump is desired, run the H-K-ATPase pump:
                f_H3 = stb.pumpVATP(self.cc_cells[self.iH],self.cc_env[self.iH],self.vm,self.T,p,self.VATP_block)
                self.fluxes_mem[self.iH] = self.fluxes_mem[self.iH] + self.rho_pump * f_H3[cells.mem_to_cells]

            # Update the anion (bicarbonate) concentration instead of H+, assuming bicarb buffer holds:
            self.update_C(self.iM, -f_H3, cells, p)

            # Calculate the new pH and H+ concentration:
             # run the bicarbonate buffer to ensure realistic concentrations and pH in cell and environment:

            self.cc_cells[self.iH], _, self.cc_cells[self.iM], self.pH_cell = stb.bicarbonate_buffer(
                 self.cc_cells[self.iH],
                 self.cHM_cells,
                 self.cc_cells[self.iM],
                 p)

            self.cc_env[self.iH], _, self.cc_env[self.iM], self.pH_env = stb.bicarbonate_buffer(
                 self.cc_env[self.iH],
                 self.cHM_env,
                 self.cc_env[self.iM],
                 p)


            # recalculate the net, unbalanced charge and voltage in each cell:
            self.update_V(cells,p)

    def ca_handler(self,cells,p):

        if p.sim_ECM is True:

            # run Ca-ATPase

            f_CaATP = stb.pumpCaATP(self.cc_cells[self.iCa][cells.mem_to_cells],
                self.cc_env[self.iCa][cells.map_mem2ecm],
                self.vm, self.T, p)

            f_CaATP = self.rho_pump * f_CaATP

            # add Ca++ flux to storage:
            self.fluxes_mem[self.iCa] = f_CaATP

        else:

            # run Ca-ATPase

            f_CaATP = stb.pumpCaATP(self.cc_cells[self.iCa], self.cc_env[self.iCa], self.vm, self.T, p)

            # store the transmembrane flux for this ion
            self.fluxes_mem[self.iCa] = self.rho_pump*f_CaATP[cells.mem_to_cells]



        # update calcium concentrations in cell and ecm:
        self.update_C(self.iCa, f_CaATP, cells, p)


        # recalculate the net, unbalanced charge and voltage in each cell:
        self.update_V(cells, p)

        if p.Ca_dyn == 1:  # do endoplasmic reticulum handling  # FIXME not right!

            f_Ca_ER = stb.pumpCaER(
                self.cc_er[0],
                self.cc_cells[self.iCa],
                self.v_er,
                self.T,
                p,
            )

            # update calcium concentrations in the ER and cell:
            self.cc_er[0] = self.cc_er[0] + f_Ca_ER * ((cells.cell_sa) / (p.ER_vol * cells.cell_vol)) * p.dt
            self.cc_cells[self.iCa] = self.cc_cells[self.iCa] - f_Ca_ER * (
                cells.cell_sa / cells.cell_vol) * p.dt

            # recalculate the net, unbalanced charge and voltage in each cell:
            self.update_V(cells, p)

            # calculate the net, unbalanced charge and voltage in the endoplasmic reticulum:
            q_er = stb.get_charge(self.cc_er, self.z_array_er, p.ER_vol * cells.cell_vol, p)
            v_er_o = stb.get_volt(q_er, p.ER_sa * cells.cell_sa, p)

            self.v_er = v_er_o - self.v_cell

    def molecule_mover(self,cX_cell,cX_env,cells,p):

        # ------------------------------------------------------------
        # pump dye and update result

        if p.sim_ECM is True:

            if p.pump_Dye is True:

                if p.pump_Dye_out is True:

                    # active pumping of dye from environment and into cell
                    deltaGATP = 20 * p.R * self.T

                    delG_dye = p.R * p.T * np.log(self.cDye_cell[cells.mem_to_cells] / self.cDye_env[cells.map_mem2ecm]) \
                               + p.z_Dye * p.F * self.vm

                    delG_dyeATP = deltaGATP - delG_dye
                    delG_pump = (delG_dyeATP / 1000)

                    # alpha = p.pump_Dye_alpha*tb.step(delG_pump,6,3)
                    alpha = p.pump_Dye_alpha * delG_pump

                    f_Dye_pump = self.rho_pump * alpha * (self.cDye_env[cells.map_mem2ecm])

                    d_dye_cells = self.rho_pump * f_Dye_pump * (cells.mem_sa / cells.cell_vol[cells.mem_to_cells])
                    d_dye_env = -self.rho_pump * f_Dye_pump * (cells.mem_sa / cells.ecm_vol[cells.map_mem2ecm])

                    delta_cells = np.dot(d_dye_cells, cells.cell_UpdateMatrix)
                    delta_env = np.dot(d_dye_env, cells.ecm_UpdateMatrix)

                else:

                    # active pumping of dye from environment and into cell
                    deltaGATP = 20 * p.R * self.T

                    delG_dye = p.R * p.T * np.log(self.cDye_env[cells.map_mem2ecm] / self.cDye_cell[cells.mem_to_cells]) \
                               - p.z_Dye * p.F * self.vm

                    delG_dyeATP = deltaGATP - delG_dye
                    delG_pump = (delG_dyeATP / 1000)

                    # alpha = p.pump_Dye_alpha*tb.step(delG_pump,6,3)
                    alpha = p.pump_Dye_alpha * delG_pump

                    f_Dye_pump = self.rho_pump * alpha * (self.cDye_cell[cells.mem_to_cells])

                    d_dye_cells = -self.rho_pump * f_Dye_pump * (cells.mem_sa / cells.cell_vol[cells.mem_to_cells])
                    d_dye_env = self.rho_pump * f_Dye_pump * (cells.mem_sa / cells.ecm_vol[cells.map_mem2ecm])

                    delta_cells = np.dot(d_dye_cells, cells.cell_UpdateMatrix)
                    delta_env = np.dot(d_dye_env, cells.ecm_UpdateMatrix)

                self.cDye_cell = self.cDye_cell + delta_cells * p.dt

                self.cDye_env = self.cDye_env + delta_env * p.dt

                # ensure that there are no negative values
                self.cDye_cell = stb.no_negs(self.cDye_cell)
                self.cDye_env = stb.no_negs(self.cDye_env)

        elif p.sim_ECM is False:

            if p.pump_Dye is True:

                if p.pump_Dye_out is True:

                    # active pumping of dye from environment and into cell
                    deltaGATP = 20 * p.R * self.T

                    delG_dye = p.R * p.T * np.log(self.cDye_cell / self.cDye_env) \
                               + p.z_Dye * p.F * self.vm

                    delG_dyeATP = deltaGATP - delG_dye
                    delG_pump = (delG_dyeATP / 1000)

                    # alpha = p.pump_Dye_alpha*tb.step(delG_pump,6,3)
                    alpha = p.pump_Dye_alpha * delG_pump

                    f_Dye_pump = self.rho_pump * alpha * (self.cDye_env)

                    d_dye_cells = f_Dye_pump * (cells.cell_sa / cells.cell_vol)

                    # delta_cells =  np.dot(d_dye_cells, cells.cell_UpdateMatrix)


                else:

                    # active pumping of dye from cell and into environment
                    deltaGATP = 20 * p.R * self.T

                    delG_dye = p.R * p.T * np.log(self.cDye_env / self.cDye_cell) \
                               - p.z_Dye * p.F * self.vm

                    delG_dyeATP = deltaGATP - delG_dye
                    delG_pump = (delG_dyeATP / 1000)

                    # alpha = p.pump_Dye_alpha*tb.step(delG_pump,6,3)
                    alpha = p.pump_Dye_alpha * delG_pump

                    f_Dye_pump = self.rho_pump * alpha * (self.cDye_cell)

                    d_dye_cells = -f_Dye_pump * (cells.cell_sa / cells.cell_vol)

                    # delta_cells =  np.dot(d_dye_cells, cells.cell_UpdateMatrix)

                self.cDye_cell = self.cDye_cell + d_dye_cells * p.dt

        # electrodiffuse dye between cell and extracellular space--------------------------------------------------
        if p.sim_ECM is False:

            fdye_ED = stb.electroflux(self.cDye_env, self.cDye_cell, self.id_cells * p.Dm_Dye, self.tm, p.z_Dye,
                self.vm,
                self.T, p, rho=self.rho_channel_o)

            # update dye concentration
            self.cDye_cell = self.cDye_cell + fdye_ED * (cells.cell_sa / cells.cell_vol) * p.dt

        elif p.sim_ECM is True:

            flux_dye = stb.electroflux(self.cDye_env[cells.map_mem2ecm], self.cDye_cell[cells.mem_to_cells],
                np.ones(len(cells.mem_i)) * p.Dm_Dye, self.tm, p.z_Dye, self.vm, self.T, p)

            # update the dye concentrations in the cell and ecm due to ED fluxes at membrane
            d_c_cells = self.rho_channel * flux_dye * (cells.mem_sa / cells.cell_vol[cells.mem_to_cells])
            d_c_env = self.rho_channel * flux_dye * (cells.mem_sa / cells.ecm_vol[cells.map_mem2ecm])

            delta_cells = np.dot(d_c_cells, cells.cell_UpdateMatrix)
            delta_env = np.dot(d_c_env, cells.ecm_UpdateMatrix)

            self.cDye_cell = self.cDye_cell + delta_cells * p.dt

            self.cDye_env = self.cDye_env - delta_env * p.dt

            # ensure that there are no negative values
            self.cDye_cell = stb.no_negs(self.cDye_cell)
            self.cDye_env = stb.no_negs(self.cDye_env)

        # ------------------------------------------------------------

        # Update dye concentration in the gj connected cell network:

        # voltage gradient:
        grad_vgj = self.vgj / cells.gj_len

        grad_vgj_x = grad_vgj * cells.cell_nn_tx
        grad_vgj_y = grad_vgj * cells.cell_nn_ty

        # concentration gradient for Dye:
        Dye_mems = self.cDye_cell[cells.mem_to_cells]

        grad_cgj = (Dye_mems[cells.nn_i] - Dye_mems[cells.mem_i]) / cells.gj_len

        grad_cgj_x = grad_cgj * cells.cell_nn_tx
        grad_cgj_y = grad_cgj * cells.cell_nn_ty

        # midpoint concentration:
        cdye = (Dye_mems[cells.nn_i] + Dye_mems[cells.mem_i]) / 2

        # electroosmotic fluid velocity:
        if p.fluid_flow is True:

            ux = (self.u_cells_x[cells.cell_nn_i[:, 0]] + self.u_cells_x[cells.cell_nn_i[:, 1]]) / 2
            uy = (self.u_cells_y[cells.cell_nn_i[:, 0]] + self.u_cells_y[cells.cell_nn_i[:, 1]]) / 2

        else:
            ux = 0
            uy = 0

        fgj_x_dye, fgj_y_dye = stb.nernst_planck_flux(cdye, grad_cgj_x, grad_cgj_y, grad_vgj_x, grad_vgj_y, ux, uy,
            p.Do_Dye * self.gjopen, p.z_Dye, self.T, p)

        fgj_dye = fgj_x_dye * cells.cell_nn_tx + fgj_y_dye * cells.cell_nn_ty

        # divergence calculation for individual cells (finite volume expression)
        delta_cc = np.dot(cells.gjMatrix * p.gj_surface * self.gjopen, -fgj_dye * cells.mem_sa) / cells.cell_vol

        self.cDye_cell = self.cDye_cell + p.dt * delta_cc

        self.Dye_flux_x_gj = fgj_x_dye[:]  # store gap junction flux for this ion
        self.Dye_flux_y_gj = fgj_y_dye[:]  # store gap junction flux for this ion

        # transport dye through environment: _________________________________________________________
        if p.sim_ECM is True:

            if p.closed_bound is True:
                btag = 'closed'

            else:
                btag = 'open'
                # make v_env and cc_env into 2d matrices
            cenv = self.cDye_env[:]
            denv = p.Do_Dye * np.ones(len(cells.xypts))

            v_env = self.v_env.reshape(cells.X.shape)

            v_env[:, 0] = self.bound_V['L']
            v_env[:, -1] = self.bound_V['R']
            v_env[0, :] = self.bound_V['B']
            v_env[-1, :] = self.bound_V['T']

            cenv = cenv.reshape(cells.X.shape)

            # prepare concentrations and diffusion constants for MACs grid format
            # by resampling the values at the u v coordinates of the flux:
            cenv_x = np.zeros(cells.grid_obj.u_shape)
            cenv_y = np.zeros(cells.grid_obj.v_shape)

            # create the proper shape for the concentrations and state appropriate boundary conditions::
            cenv_x[:, 1:] = cenv[:]
            cenv_x[:, 0] = cenv_x[:, 1]
            cenv_y[1:, :] = cenv[:]
            cenv_y[0, :] = cenv_y[1, :]

            if p.closed_bound is True:  # insulation boundary conditions
                cenv_x[:, 0] = cenv_x[:, 1]
                cenv_x[:, -1] = cenv_x[:, -2]
                cenv_x[0, :] = cenv_x[1, :]
                cenv_x[-1, :] = cenv_x[-2, :]

                cenv_y[0, :] = cenv_y[1, :]
                cenv_y[-1, :] = cenv_y[-2, :]
                cenv_y[:, 0] = cenv_y[:, 1]
                cenv_y[:, -1] = cenv_y[:, -2]

            else:  # open and electrically grounded boundary conditions
                cenv_x[:, 0] = self.c_dye_bound
                cenv_x[:, -1] = self.c_dye_bound
                cenv_x[0, :] = self.c_dye_bound
                cenv_x[-1, :] = self.c_dye_bound

                cenv_y[0, :] = self.c_dye_bound
                cenv_y[-1, :] = self.c_dye_bound
                cenv_y[:, 0] = self.c_dye_bound
                cenv_y[:, -1] = self.c_dye_bound

            denv = denv.reshape(cells.X.shape)

            denv_x = interp.griddata((cells.xypts[:, 0], cells.xypts[:, 1]), denv.ravel(),
                (cells.grid_obj.u_X, cells.grid_obj.u_Y), method='nearest', fill_value=p.Do_Dye)

            denv_y = interp.griddata((cells.xypts[:, 0], cells.xypts[:, 1]), denv.ravel(),
                (cells.grid_obj.v_X, cells.grid_obj.v_Y), method='nearest', fill_value=p.Do_Dye)

            # denv_x = denv_x*self.D_env_weight_u
            # denv_y = denv_y*self.D_env_weight_v

            # calculate gradients in the environment
            grad_V_env_x, grad_V_env_y = cells.grid_obj.grid_gradient(v_env, bounds='closed')

            grad_cc_env_x, grad_cc_env_y = cells.grid_obj.grid_gradient(cenv, bounds=btag)

            # calculate fluxes for electrodiffusive transport:

            if p.fluid_flow is True:

                uenvx = np.zeros(cells.grid_obj.u_shape)
                uenvy = np.zeros(cells.grid_obj.v_shape)

                uenvx[:, 1:] = self.u_env_x
                uenvy[1:, :] = self.u_env_y

                if p.closed_bound is False:

                    uenvx[:, 0] = uenvx[:, 1]
                    uenvx[:, -1] = uenvx[:, -2]
                    uenvx[0, :] = uenvx[1, :]
                    uenvx[-1, :] = uenvx[-2, :]

                    uenvy[:, 0] = uenvy[:, 1]
                    uenvy[:, -1] = uenvy[:, -2]
                    uenvy[0, :] = uenvy[1, :]
                    uenvy[-1, :] = uenvy[-2, :]

                else:

                    uenvx[:, 0] = 0
                    uenvx[:, -1] = 0
                    uenvx[0, :] = 0
                    uenvx[-1, :] = 0

                    uenvy[:, 0] = 0
                    uenvy[:, -1] = 0
                    uenvy[0, :] = 0
                    uenvy[-1, :] = 0

            else:
                uenvx = 0
                uenvy = 0

            f_env_x_dye, f_env_y_dye = stb.np_flux_special(cenv_x, cenv_y, grad_cc_env_x, grad_cc_env_y,
                grad_V_env_x, grad_V_env_y, uenvx, uenvy, denv_x, denv_y, p.z_Dye, self.T, p)

            # calculate the divergence of the total (negative) flux to obtain the total change per unit time:
            d_fenvx = -(f_env_x_dye[:, 1:] - f_env_x_dye[:, 0:-1]) / cells.delta
            d_fenvy = -(f_env_y_dye[1:, :] - f_env_y_dye[0:-1, :]) / cells.delta

            delta_c = d_fenvx + d_fenvy

            cenv = cenv + delta_c * p.dt

            cenv = fd.integrator(cenv)

            if p.closed_bound is True:
                # Neumann boundary condition (flux at boundary)
                # zero flux boundaries for concentration:
                cenv[:, -1] = cenv[:, -2]
                cenv[:, 0] = cenv[:, 1]
                cenv[0, :] = cenv[1, :]
                cenv[-1, :] = cenv[-2, :]

            elif p.closed_bound is False:
                # if the boundary is open, set the concentration at the boundary
                # open boundary
                cenv[:, -1] = self.c_dye_bound
                cenv[:, 0] = self.c_dye_bound
                cenv[0, :] = self.c_dye_bound
                cenv[-1, :] = self.c_dye_bound

            # reshape the matrices into vectors:
            # self.v_env = self.v_env.ravel()
            self.cDye_env = cenv.ravel()

            # average flux at the midpoint of the MACs grid:
            fenvx = (f_env_x_dye[:, 1:] + f_env_x_dye[:, 0:-1]) / 2
            fenvy = (f_env_y_dye[1:, :] + f_env_y_dye[0:-1, :]) / 2

            self.Dye_flux_env_x = fenvx.ravel()  # store ecm junction flux for this ion
            self.Dye_flux_env_y = fenvy.ravel()  # store ecm junction flux for this ion

            # ensure that there are no negative values
            self.cDye_cell = stb.no_negs(self.cDye_cell)
            self.cDye_env = stb.no_negs(self.cDye_env)

    def update_gj(self,cells,p,t,i):

        # calculate voltage difference (gradient*len_gj) between gj-connected cells:
        if p.sim_ECM is True:

            vmems = self.v_cell
            # vmems = self.vm

            self.vgj = self.v_cell[cells.cell_nn_i[:,1]]- self.v_cell[cells.cell_nn_i[:,0]]

        else:

            self.vgj = self.vm[cells.cell_nn_i[:,1]] - self.vm[[cells.cell_nn_i[:,0]]]


        if p.v_sensitive_gj is True:
            # determine the open state of gap junctions:
            # self.gjopen = self.gj_block*((1.0 - tb.step(abs(self.vgj),p.gj_vthresh,p.gj_vgrad) + 0.1))

            # calculate the steady state value of gj conductivity for this voltage:

            # gj_inf = ((1 - 0.04)/(1 + np.exp(0.217*(self.vgj*1e3 - p.gj_vthresh)))) + 0.04

            gmin = self.gj_funk.gmin

            alpha_gj = self.gj_funk.alpha_gj
            beta_gj = self.gj_funk.beta_gj_p

            vgj = np.abs(self.vgj)


            dgjopen_dt = (1 - self.gjopen)*alpha_gj(vgj) - (self.gjopen - gmin)*beta_gj(vgj)

            self.gjopen = self.gjopen + 1e3*dgjopen_dt*p.dt

            # threshold to ensure 0 to 1 status
            inds_gj_over = (self.gjopen > 1.0).nonzero()
            self.gjopen[inds_gj_over] = 1.0

            inds_gj_under = (self.gjopen < 0.0).nonzero()
            self.gjopen[inds_gj_under] = 0.0

            self.gjopen = self.gj_block*self.gjopen

            # print(self.gjopen.min(),self.gjopen.max())
            # print('---------------------------')



        else:
            self.gjopen = self.gj_block


        # voltage gradient:
        grad_vgj = self.vgj/cells.gj_len

        grad_vgj_x = grad_vgj*cells.cell_nn_tx
        grad_vgj_y = grad_vgj*cells.cell_nn_ty

        # concentration gradient for ion i:

        conc_mem = self.cc_cells[i]
        grad_cgj = (conc_mem[cells.cell_nn_i[:,1]] - conc_mem[cells.cell_nn_i[:,0]])/cells.gj_len

        grad_cgj_x = grad_cgj*cells.cell_nn_tx
        grad_cgj_y = grad_cgj*cells.cell_nn_ty

        # midpoint concentration:
        c = (conc_mem[cells.cell_nn_i[:,1]] + conc_mem[cells.cell_nn_i[:,0]])/2

        # electroosmotic fluid velocity -- averaged at gap junctions:
        if p.fluid_flow is True:
            ux = (self.u_cells_x[cells.cell_nn_i[:,0]] + self.u_cells_x[cells.cell_nn_i[:,1]])/2
            uy = (self.u_cells_y[cells.cell_nn_i[:,0]] + self.u_cells_y[cells.cell_nn_i[:,1]])/2

        else:
            ux = 0
            uy =0

        fgj_x,fgj_y = stb.nernst_planck_flux(c,grad_cgj_x,grad_cgj_y,grad_vgj_x,grad_vgj_y,ux,uy,
            p.gj_surface*self.gjopen*self.D_gj[i],self.zs[i],self.T,p)

        # component of flux tangent to gap junctions:
        fgj = fgj_x*cells.cell_nn_tx + fgj_y*cells.cell_nn_ty

        # divergence calculation (finite volume expression)
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

            # try adding it to the second layer to get proper diffusion in...

            cenv_x[:,1] =  self.c_env_bound[i]
            cenv_x[:,-2] =  self.c_env_bound[i]
            cenv_x[1,:] =  self.c_env_bound[i]
            cenv_x[-2,:] =  self.c_env_bound[i]

            cenv_y[1,:] =  self.c_env_bound[i]
            cenv_y[-2,:] =  self.c_env_bound[i]
            cenv_y[:,1] =  self.c_env_bound[i]
            cenv_y[:,-2] =  self.c_env_bound[i]

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

        f_env_x, f_env_y = stb.np_flux_special(cenv_x,cenv_y,grad_cc_env_x,grad_cc_env_y,
            grad_V_env_x, grad_V_env_y, uenvx,uenvy,self.D_env_u[i],self.D_env_v[i],
            self.zs[i],self.T,p)

        if p.closed_bound is False:

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

        # delta_c = fd.integrator(delta_c)

        #-----------------------
        cenv = cenv + delta_c*p.dt

        cenv = fd.integrator(cenv)  # smooth out the concentration

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

            cenv[:,-2] = self.c_env_bound[i]
            cenv[:,1] = self.c_env_bound[i]
            cenv[1,:] = self.c_env_bound[i]
            cenv[-2,:] = self.c_env_bound[i]

        # reshape the matrices back into vectors:
        self.cc_env[i] = cenv.ravel()

        fenvx = (f_env_x[:,1:] + f_env_x[:,0:-1])/2
        fenvy = (f_env_y[1:,:] + f_env_y[0:-1,:])/2

        # fenvx = fd.integrator(fenvx)
        # fenvy = fd.integrator(fenvy)

        self.fluxes_env_x[i] = fenvx.ravel()  # store ecm junction flux for this ion
        self.fluxes_env_y[i] = fenvy.ravel()  # store ecm junction flux for this ion

    # FIXME get the following in to Ca handler

    def update_er(self,cells,p,t):

         # electrodiffusion of ions between cell and endoplasmic reticulum
        f_Ca_ER = \
        stb.electroflux(self.cc_cells[self.iCa],self.cc_er[0],self.Dm_er[0],self.tm,self.z_er[0],self.v_er,self.T,p)

        # Electrodiffusion of charge compensation anion
        f_M_ER = \
        stb.electroflux(self.cc_cells[self.iM],self.cc_er[1],self.Dm_er[1],self.tm,self.z_er[1],self.v_er,self.T,p)

        # update calcium concentrations in the ER and cell:
        self.cc_er[0] = self.cc_er[0] + f_Ca_ER*((cells.cell_sa)/(p.ER_vol*cells.cell_vol))*p.dt
        self.cc_cells[self.iCa] = self.cc_cells[self.iCa] - f_Ca_ER*(cells.cell_sa/cells.cell_vol)*p.dt

        # update anion concentration in the ER and cell:
        self.cc_er[1] = self.cc_er[1] + f_M_ER*((cells.cell_sa)/(p.ER_vol*cells.cell_vol))*p.dt
        self.cc_cells[self.iM] = self.cc_cells[self.iM] - f_M_ER*(cells.cell_sa/cells.cell_vol)*p.dt

        # recalculate the net, unbalanced charge and voltage in each cell:
        self.update_V(cells,p,t)

        q_er = stb.get_charge(self.cc_er,self.z_array_er,p.ER_vol*cells.cell_vol,p)
        self.v_er = stb.get_volt(q_er,p.ER_sa*cells.cell_sa,p) - self.v_cell

    # FIXME the following 2 are to be done with the new molecule mover, streamlined for +/- ecm

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

                    d_dye_cells = self.rho_pump*f_Dye_pump*(cells.mem_sa/cells.cell_vol[cells.mem_to_cells])
                    d_dye_env = -self.rho_pump*f_Dye_pump*(cells.mem_sa/cells.ecm_vol[cells.map_mem2ecm])

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

                    d_dye_cells = -self.rho_pump*f_Dye_pump*(cells.mem_sa/cells.cell_vol[cells.mem_to_cells])
                    d_dye_env = self.rho_pump*f_Dye_pump*(cells.mem_sa/cells.ecm_vol[cells.map_mem2ecm])

                    delta_cells =  np.dot(d_dye_cells, cells.cell_UpdateMatrix)
                    delta_env = np.dot(d_dye_env, cells.ecm_UpdateMatrix)

                self.cDye_cell = self.cDye_cell + delta_cells*p.dt

                self.cDye_env = self.cDye_env + delta_env*p.dt

                # ensure that there are no negative values
                self.cDye_cell = stb.no_negs(self.cDye_cell)
                self.cDye_env = stb.no_negs(self.cDye_env)

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

                    f_Dye_pump  = self.rho_pump*alpha*(self.cDye_env)

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

                    f_Dye_pump  = self.rho_pump*alpha*(self.cDye_cell)

                    d_dye_cells = -f_Dye_pump*(cells.cell_sa/cells.cell_vol)

                    # delta_cells =  np.dot(d_dye_cells, cells.cell_UpdateMatrix)

                self.cDye_cell = self.cDye_cell + d_dye_cells*p.dt

        # electrodiffuse dye between cell and extracellular space--------------------------------------------------
        if p.sim_ECM is False:

            fdye_ED = stb.electroflux(self.cDye_env,self.cDye_cell,self.id_cells*p.Dm_Dye,self.tm,p.z_Dye,self.vm,
                self.T,p,rho=self.rho_channel_o)

            # update dye concentration
            self.cDye_cell = self.cDye_cell + fdye_ED*(cells.cell_sa/cells.cell_vol)*p.dt

        elif p.sim_ECM is True:

            flux_dye = stb.electroflux(self.cDye_env[cells.map_mem2ecm],self.cDye_cell[cells.mem_to_cells],
                            np.ones(len(cells.mem_i))*p.Dm_Dye,self.tm,p.z_Dye,self.vm,self.T,p)

            # update the dye concentrations in the cell and ecm due to ED fluxes at membrane
            d_c_cells = self.rho_channel*flux_dye*(cells.mem_sa/cells.cell_vol[cells.mem_to_cells])
            d_c_env = self.rho_channel*flux_dye*(cells.mem_sa/cells.ecm_vol[cells.map_mem2ecm])

            delta_cells =  np.dot(d_c_cells, cells.cell_UpdateMatrix)
            delta_env = np.dot(d_c_env, cells.ecm_UpdateMatrix)

            self.cDye_cell = self.cDye_cell + delta_cells*p.dt

            self.cDye_env = self.cDye_env - delta_env*p.dt

            # ensure that there are no negative values
            self.cDye_cell = stb.no_negs(self.cDye_cell)
            self.cDye_env = stb.no_negs(self.cDye_env)

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

        fgj_x_dye,fgj_y_dye = stb.nernst_planck_flux(cdye,grad_cgj_x,grad_cgj_y,grad_vgj_x,grad_vgj_y,ux,uy,
            p.Do_Dye*self.gjopen,p.z_Dye,self.T,p)

        fgj_dye = fgj_x_dye*cells.cell_nn_tx + fgj_y_dye*cells.cell_nn_ty

        # divergence calculation for individual cells (finite volume expression)
        delta_cc = np.dot(cells.gjMatrix*p.gj_surface*self.gjopen,-fgj_dye*cells.mem_sa)/cells.cell_vol

        self.cDye_cell = self.cDye_cell + p.dt*delta_cc

        self.Dye_flux_x_gj = fgj_x_dye[:]  # store gap junction flux for this ion
        self.Dye_flux_y_gj = fgj_y_dye[:]  # store gap junction flux for this ion

        # transport dye through environment: _________________________________________________________
        if p.sim_ECM is True:

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

            # denv_x = denv_x*self.D_env_weight_u
            # denv_y = denv_y*self.D_env_weight_v

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

            f_env_x_dye, f_env_y_dye = stb.np_flux_special(cenv_x,cenv_y,grad_cc_env_x,grad_cc_env_y,
                grad_V_env_x, grad_V_env_y, uenvx,uenvy,denv_x,denv_y,p.z_Dye,self.T,p)

            # calculate the divergence of the total (negative) flux to obtain the total change per unit time:
            d_fenvx = -(f_env_x_dye[:,1:] - f_env_x_dye[:,0:-1])/cells.delta
            d_fenvy = -(f_env_y_dye[1:,:] - f_env_y_dye[0:-1,:])/cells.delta

            delta_c = d_fenvx + d_fenvy

            cenv = cenv + delta_c*p.dt

            cenv = fd.integrator(cenv)

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

            # average flux at the midpoint of the MACs grid:
            fenvx = (f_env_x_dye[:,1:] + f_env_x_dye[:,0:-1])/2
            fenvy = (f_env_y_dye[1:,:] + f_env_y_dye[0:-1,:])/2

            self.Dye_flux_env_x = fenvx.ravel()  # store ecm junction flux for this ion
            self.Dye_flux_env_y = fenvy.ravel()  # store ecm junction flux for this ion

            # ensure that there are no negative values
            self.cDye_cell = stb.no_negs(self.cDye_cell)
            self.cDye_env = stb.no_negs(self.cDye_env)

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

        fgj_x_ip3,fgj_y_ip3 = stb.nernst_planck_flux(cip3,grad_cgj_x,grad_cgj_y,grad_vgj_x,grad_vgj_y,ux,uy,
            p.Do_IP3*self.gjopen,p.z_IP3,self.T,p)

        fgj_ip3 = fgj_x_ip3*cells.cell_nn_tx + fgj_y_ip3*cells.cell_nn_ty

        delta_cc = np.dot(cells.gjMatrix*p.gj_surface*self.gjopen,-fgj_ip3*cells.mem_sa)/cells.cell_vol

        self.cIP3 = self.cIP3 + p.dt*delta_cc

        self.IP3_flux_x_gj = fgj_x_ip3[:]  # store gap junction flux for this ion
        self.IP3_flux_y_gj = fgj_y_ip3[:]  # store gap junction flux for this ion

        if p.sim_ECM is False:

            fip3_ED = stb.electroflux(self.cIP3_env,self.cIP3,self.id_cells*p.Dm_IP3,self.tm,p.z_IP3,self.vm,self.T,p)

            # update dye concentration
            self.cIP3 = self.cIP3 + fip3_ED*(cells.cell_sa/cells.cell_vol)*p.dt

        elif p.sim_ECM is True:

            flux_ip3 = stb.electroflux(self.cIP3_env[cells.map_mem2ecm],self.cIP3[cells.mem_to_cells],
                            np.ones(len(cells.mem_i))*p.Dm_IP3,self.tm,p.z_IP3,self.vm,self.T,p)

            # update the dye concentrations in the cell and ecm due to ED fluxes at membrane
            d_c_cells = flux_ip3*(cells.mem_sa/cells.cell_vol[cells.mem_to_cells])
            d_c_env = flux_ip3*(cells.mem_sa/cells.ecm_vol)

            delta_cells =  np.dot(d_c_cells, cells.cell_UpdateMatrix)
            delta_env = np.dot(d_c_env, cells.ecm_UpdateMatrix)

            self.cIP3 = self.cIP3 + delta_cells*p.dt

            self.cIP3_env = self.cIP3_env - delta_env*p.dt

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

            f_env_x_ip3, f_env_y_ip3 = stb.np_flux_special(cenv_x,cenv_y,grad_cc_env_x,grad_cc_env_y,
                grad_V_env_x, grad_V_env_y, uenvx,uenvy,denv_x,denv_y,p.z_IP3,self.T,p)

            # calculate the divergence of the total flux, which is equivalent to the total change per unit time:
            d_fenvx = -(f_env_x_ip3[:,1:] - f_env_x_ip3[:,0:-1])/cells.delta
            d_fenvy = -(f_env_y_ip3[1:,:] - f_env_y_ip3[0:-1,:])/cells.delta

            delta_c = d_fenvx + d_fenvy

            cenv = cenv + delta_c*p.dt

            cenv = fd.integrator(cenv)

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

            # true flux is indeed negative
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

            self.E_env_x = -genv_x.ravel()*cells.ave2ecmV
            self.E_env_y = -genv_y.ravel()*cells.ave2ecmV

            self.E_env_x = gaussian_filter(self.E_env_x.reshape(cells.X.shape),2)
            self.E_env_y = gaussian_filter(self.E_env_y.reshape(cells.X.shape),2)

        else:

            self.Egj = - (self.vm[cells.cell_nn_i[:,1]] - self.vm[cells.cell_nn_i[:,0]])/cells.gj_len

        # self.Egj = self.Egj*cells.ave2cellV # scale from cell space volume to whole cell volume

        # get x and y components of the electric field:
        self.E_gj_x = cells.cell_nn_tx*self.Egj
        self.E_gj_y = cells.cell_nn_ty*self.Egj

    def get_current(self,cells,p):

        # zero components of total current vector:

        self.I_tot_x = np.zeros(cells.X.shape)
        self.I_tot_y = np.zeros(cells.X.shape)

        # calculate current across gap junctions in x direction:
        I_gj_x = np.zeros(len(cells.mem_i))

        for flux_array, zi in zip(self.fluxes_gj_x,self.zs):

            I_i_x = flux_array*zi*p.F*cells.mem_sa

            I_gj_x = I_gj_x + I_i_x

        # calculate current across gap junctions in x direction:
        I_gj_y = np.zeros(len(cells.mem_i))

        for flux_array, zi in zip(self.fluxes_gj_y,self.zs):

            I_i_y = flux_array*zi*p.F*cells.mem_sa

            I_gj_y = I_gj_y + I_i_y

        # interpolate the gj current components to the grid:
        self.I_gj_x = interp.griddata((cells.mem_mids_flat[:,0],cells.mem_mids_flat[:,1]),I_gj_x,(cells.X,cells.Y),
                                      method=p.interp_type,fill_value=0)

        # self.I_gj_x = np.multiply(self.I_gj_x,cells.maskECM)

        self.I_gj_y = interp.griddata((cells.mem_mids_flat[:,0],cells.mem_mids_flat[:,1]),I_gj_y,(cells.X,cells.Y),
                                      method=p.interp_type,fill_value=0)

        self.I_gj_x = self.I_gj_x/(cells.delta*p.cell_height)
        self.I_gj_y = self.I_gj_y/(cells.delta*p.cell_height)

        # self.I_gj_y = np.multiply(self.I_gj_y,cells.maskECM)

        self.I_tot_x = self.I_tot_x + self.I_gj_x
        self.I_tot_y = self.I_tot_y + self.I_gj_y

        # calculate current across cell membranes:

        self.I_mem = np.zeros(len(cells.mem_i))
        for flux_array, zi in zip(self.fluxes_mem,self.zs):

            I_i = flux_array*zi*p.F*cells.mem_sa

            self.I_mem = self.I_mem + I_i

            # components are negative as transmembrane fluxes point into the cell, but mem normals point out:
            I_mem_x = -self.I_mem*cells.mem_vects_flat[:,2]
            I_mem_y = -self.I_mem*cells.mem_vects_flat[:,3]

         # interpolate the trans-membrane current components to the grid:
        self.I_mem_x = interp.griddata((cells.mem_vects_flat[:,0],cells.mem_vects_flat[:,1]),I_mem_x,(cells.X,cells.Y),
                                      method=p.interp_type,fill_value=0)
        # self.I_mem_x = np.multiply(self.I_mem_x,cells.maskM)

        self.I_mem_y = interp.griddata((cells.mem_vects_flat[:,0],cells.mem_vects_flat[:,1]),I_mem_y,(cells.X,cells.Y),
                                      method=p.interp_type,fill_value=0)
        # self.I_mem_y = np.multiply(self.I_mem_y,cells.maskM)

        self.I_mem_x = self.I_mem_x/(cells.delta*p.cell_height)
        self.I_mem_y = self.I_mem_y/(cells.delta*p.cell_height)

        # add membrane current to total current:

        self.I_tot_x = self.I_tot_x + self.I_mem_x
        self.I_tot_y = self.I_tot_y + self.I_mem_y

        if p.sim_ECM is True:

            self.I_env_x = np.zeros(len(cells.xypts))
            self.I_env_y = np.zeros(len(cells.xypts))

            for flux_array, zi in zip(self.fluxes_env_x,self.zs):

                I_i = flux_array*zi*p.F*p.cell_space*p.cell_height

                self.I_env_x = self.I_env_x + I_i

            for flux_array, zi in zip(self.fluxes_env_y,self.zs):

                I_i = flux_array*zi*p.F*p.cell_space*p.cell_height

                self.I_env_y = self.I_env_y + I_i

            I_env_x = self.I_env_x.reshape(cells.X.shape)/(cells.delta*p.cell_height)
            I_env_y = self.I_env_y.reshape(cells.X.shape)/(cells.delta*p.cell_height)

            self.I_tot_x = self.I_tot_x + I_env_x
            self.I_tot_y = self.I_tot_y + I_env_y

    # FIXME consider moving to a flow module

    def getFlow_o(self,cells,p):

        """
        Calculate the electroosmotic fluid flow in the cell and extracellular
        networks.

        """

        if p.sim_ECM is True:

            # force of gravity:
            if p.closed_bound is True:

                btag = 'closed'

            else:

                btag = 'open'

            # estimate the inverse viscosity for extracellular flow based on the diffusion constant weighting
            # for the world:
            alpha = (1/p.mu_water)*self.D_env_weight*1.0e-6

            E_ave_x = self.E_env_x
            E_ave_y = self.E_env_y

            if p.deform_electro is True:

                # determine the geometric factor to map charge density from volume to surface:
                # Qfactor = p.cell_space
                #
                # Fe_x = Qfactor*self.rho_env.reshape(cells.X.shape)*E_ave_x
                # Fe_y = Qfactor*self.rho_env.reshape(cells.X.shape)*E_ave_y

                # map charge density to rectangular grid volume:
                Qenv = self.rho_env*cells.ave2ecmV

                # Qenv = Qenv.reshape(cells.X.shape)

                Qenv = gaussian_filter(Qenv.reshape(cells.X.shape),2)

                Fe_x = Qenv*E_ave_x
                Fe_y = Qenv*E_ave_y


            else:

                Fe_x = np.zeros(cells.X.shape)
                Fe_y = np.zeros(cells.Y.shape)

            # sum the forces:
            Fx = Fe_x
            Fy = Fe_y

            source_x = -Fx*alpha
            # source_x = cells.grid_obj.grid_int(source_x,bounds=btag) # perform finite volume integration of source

            source_y = -Fy*alpha
            # source_y = cells.grid_obj.grid_int(source_y,bounds=btag) # perform finite volume integration of source

            # # calculated the fluid flow using the time-independent Stokes Flow equation:
            ux_ecm_o = np.dot(cells.lapENVinv,source_x.ravel())
            uy_ecm_o = np.dot(cells.lapENVinv,source_y.ravel())

            # calculate the divergence of the flow field as the sum of the two spatial derivatives:
            div_uo = fd.divergence(ux_ecm_o.reshape(cells.X.shape),uy_ecm_o.reshape(cells.X.shape),
                cells.delta,cells.delta)

            # perform finite volume integration on the divergence:
            # div_uo = fd.integrator(div_uo)
            # div_uo = cells.grid_obj.grid_int(div_uo,bounds=btag)

            # calculate the alpha-scaled internal pressure from the divergence of the force:
            P = np.dot(cells.lapENV_P_inv, div_uo.ravel())
            P = P.reshape(cells.X.shape)

            # enforce zero normal gradient boundary conditions on P:
            P[:,0] = P[:,1]
            P[:,-1] = P[:,-2]
            P[0,:] = P[1,:]
            P[-1,:] = P[-2,:]

            # Take the grid gradient of the scaled internal pressure:
            gPx, gPy = fd.gradient(P,cells.delta)

             # subtract the pressure term from the solution to yield a divergence-free flow field
            u_env_x = ux_ecm_o.reshape(cells.X.shape) - gPx
            u_env_y = uy_ecm_o.reshape(cells.X.shape) - gPy

            # velocities at cell centres:
            self.u_env_x = u_env_x[:]
            self.u_env_y = u_env_y[:]

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


        #---------------Flow through gap junction connected cells-------------------------------------------------------

        # calculate the inverse viscocity for the cell collection, which is scaled by gj conductivity:
        alpha_gj_o = p.gj_surface*self.gjopen*(1/p.mu_water)
        # average to individual cells:
        alpha_gj = np.dot(cells.M_sum_mems,alpha_gj_o)/cells.num_mems

        if p.deform_electro is True:

            Fe_cell_x = self.F_electro_x
            Fe_cell_y = self.F_electro_y

        else:
            Fe_cell_x = np.zeros(len(cells.cell_i))
            Fe_cell_y = np.zeros(len(cells.cell_i))

        if p.deform_osmo is True:

            F_osmo_x = self.F_hydro_x
            F_osmo_y = self.F_hydro_y

        else:

            F_osmo_x = np.zeros(len(cells.cell_i))
            F_osmo_y = np.zeros(len(cells.cell_i))

        # net force is the sum of electrostatic and hydrostatic pressure induced body forces:
        F_net_x = Fe_cell_x + F_osmo_x
        F_net_y = Fe_cell_y + F_osmo_y

        # integrate body forces:
        # F_net_x = cells.integrator(F_net_x)
        # F_net_y = cells.integrator(F_net_y)

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
        P_react = np.dot(cells.lapGJ_P_inv,2*div_u)

        # calculate its gradient:
        gradP_react = (P_react[cells.cell_nn_i[:,1]] - P_react[cells.cell_nn_i[:,0]])/(cells.nn_len)

        gP_x = gradP_react*cells.cell_nn_tx
        gP_y = gradP_react*cells.cell_nn_ty

        # average the components of the reaction force field at cell centres and get boundary values:
        gPx_cell = np.dot(cells.M_sum_mems,gP_x)/cells.num_mems
        gPy_cell = np.dot(cells.M_sum_mems,gP_y)/cells.num_mems

        self.u_cells_x = u_gj_xo - gPx_cell
        self.u_cells_y = u_gj_yo - gPy_cell

        # enforce the boundary conditions:
        self.u_cells_x[cells.bflags_cells] = 0
        self.u_cells_y[cells.bflags_cells] = 0

    def getFlow(self, cells, p):

        """
        Calculate the electroosmotic fluid flow in the cell and extracellular
        networks using Hagen–Poiseuille "pipe flow" equation.

        """

        # First do extracellular space electroosmotic flow--------------------------------------------------------------

        if p.sim_ECM is True:

            # force of gravity:
            if p.closed_bound is True:

                btag = 'closed'

            else:

                btag = 'open'

            # estimate the inverse viscosity for extracellular flow based on the diffusion constant weighting
            # for the world:
            alpha = ((p.cell_space**2)/(p.mu_water*32))*self.D_env_weight

            if p.deform_electro is True:

                # map charge density to rectangular grid volume:
                Qenv = self.rho_env*cells.ave2ecmV

                Qenv = gaussian_filter(Qenv.reshape(cells.X.shape), 1)

                Fe_x = Qenv*self.E_env_x
                Fe_y = Qenv*self.E_env_y


            else:

                Fe_x = np.zeros(cells.X.shape)
                Fe_y = np.zeros(cells.Y.shape)

            # sum the forces:
            Fx = Fe_x
            Fy = Fe_y

            # calculated the base fluid flow using the Hagen-Poiseuille equation:
            ux_ecm_o = gaussian_filter(Fx*alpha,1)
            uy_ecm_o = gaussian_filter(Fy*alpha,1)

            # calculate the divergence of the flow field as the sum of the two spatial derivatives:
            div_uo = fd.divergence(ux_ecm_o, uy_ecm_o,cells.delta, cells.delta)

            # calculate the alpha-scaled internal pressure from the divergence of the force:
            P = np.dot(cells.lapENV_P_inv, div_uo.ravel())
            P = P.reshape(cells.X.shape)

            # enforce zero normal gradient boundary conditions on P:
            P[:, 0] = P[:, 1]
            P[:, -1] = P[:, -2]
            P[0, :] = P[1, :]
            P[-1, :] = P[-2, :]

            # Take the grid gradient of the scaled internal pressure:
            gPx, gPy = fd.gradient(P, cells.delta)

            # subtract the pressure term from the solution to yield a divergence-free flow field
            u_env_x = ux_ecm_o.reshape(cells.X.shape) - gPx
            u_env_y = uy_ecm_o.reshape(cells.X.shape) - gPy

            # velocities at grid-cell centres:
            self.u_env_x = u_env_x[:]
            self.u_env_y = u_env_y[:]

            # boundary conditions reinforced:
            self.u_env_x[:, 0] = 0
            # right
            self.u_env_x[:, -1] = 0
            # top
            self.u_env_x[-1, :] = 0
            # bottom
            self.u_env_x[0, :] = 0

            # left
            self.u_env_y[:, 0] = 0
            # right
            self.u_env_y[:, -1] = 0
            # top
            self.u_env_y[-1, :] = 0
            # bottom
            self.u_env_y[0, :] = 0

        #-------Next do flow through gap junction connected cells-------------------------------------------------------

        # calculate the inverse viscosity for the cell collection, which is scaled by gj state:
        alpha_gj = (1/(32*p.mu_water))*((self.gjopen*5e-10)**2)
        # alpha_gj = (1/(128*p.mu_water))*((self.gjopen*10e-9)**2)*(1/(cells.mem_sa*p.gj_surface))

        if p.deform_electro is True:

            Fe_cell_x = self.F_gj_x
            Fe_cell_y = self.F_gj_y

        else:
            Fe_cell_x = np.zeros(len(cells.mem_i))
            Fe_cell_y = np.zeros(len(cells.mem_i))

        if p.deform_osmo is True:

            F_osmo_x = self.F_hydro_x_gj
            F_osmo_y = self.F_hydro_y_gj

        else:

            F_osmo_x = np.zeros(len(cells.mem_i))
            F_osmo_y = np.zeros(len(cells.mem_i))

        # net force is the sum of electrostatic and hydrostatic pressure induced body forces:
        F_net_x = Fe_cell_x + F_osmo_x
        F_net_y = Fe_cell_y + F_osmo_y

        # Calculate flow under body forces:
        u_gj_xo = F_net_x*alpha_gj
        u_gj_yo = F_net_y*alpha_gj

        # calculate divergence as the sum of this vector x each surface area, divided by cell volume:
        u_gj = np.sqrt(u_gj_xo**2 + u_gj_yo**2)

        # calculate divergence as the sum of this vector x each surface area, divided by cell volume:
        div_u = (np.dot(cells.M_sum_mems, u_gj*cells.mem_sa*p.gj_surface*self.gjopen)/cells.cell_vol)

        # calculate the reaction pressure required to counter-balance the flow field:
        P_react = np.dot(cells.lapGJ_P_inv, div_u)

        # calculate its gradient:
        gradP_react = (P_react[cells.cell_nn_i[:, 1]] - P_react[cells.cell_nn_i[:, 0]]) / (cells.nn_len)

        gP_x = gradP_react * cells.cell_nn_tx
        gP_y = gradP_react * cells.cell_nn_ty

        u_gj_x = u_gj_xo - gP_x
        u_gj_y = u_gj_yo - gP_y

        # average the components at cell centres:
        self.u_cells_x = np.dot(cells.M_sum_mems, u_gj_x)/cells.num_mems
        self.u_cells_y = np.dot(cells.M_sum_mems, u_gj_y)/cells.num_mems

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

        # components of fluid flow velocity at the membrane:
        if p.fluid_flow is True and p.sim_ECM is True:
            ux_mem = self.u_env_x.ravel()[cells.map_mem2ecm]
            uy_mem = self.u_env_y.ravel()[cells.map_mem2ecm]

        else:
            ux_mem = 0
            uy_mem = 0

        # get the gradient of rho concentration around each membrane:

        rho_pump = np.dot(cells.M_sum_mems, self.rho_pump)/cells.num_mems
        rho_channel = np.dot(cells.M_sum_mems, self.rho_channel)/cells.num_mems

        grad_c = (rho_pump[cells.cell_nn_i[:,0]] - rho_pump[cells.cell_nn_i[:,1]])/cells.nn_len
        grad_c_ch = (rho_channel[cells.cell_nn_i[:,0]] - rho_channel[cells.cell_nn_i[:,1]])/cells.nn_len

        # get the gradient components:
        gcx = grad_c*cells.cell_nn_tx
        gcy = grad_c*cells.cell_nn_ty

        gcx_ch = grad_c_ch*cells.cell_nn_tx
        gcy_ch = grad_c_ch*cells.cell_nn_ty

        # total average electric field at each membrane
        if p.sim_ECM is True:

            Ex = self.E_env_x.ravel()[cells.map_mem2ecm]
            Ey = self.E_env_y.ravel()[cells.map_mem2ecm]

            # Ex = self.E_env_x.ravel()[cells.map_mem2ecm]
            # Ey = self.E_env_y.ravel()[cells.map_mem2ecm]

        else:
            Ex = self.E_gj_x
            Ey = self.E_gj_y

        # calculate the total Nernst-Planck flux at each membrane for rho_pump factor:

        fx_pump, fy_pump = stb.nernst_planck_flux(self.rho_pump, gcx, gcy, Ex, Ey,ux_mem,uy_mem,p.D_membrane,p.z_pump,
            self.T,p)

        # component of total flux in direction of membrane
        # ftot = fx*cells.mem_vects_flat[:,4] + fy*cells.mem_vects_flat[:,5]

        # map the flux to cell centres:
        fx_pump_o = np.dot(cells.M_sum_mems,fx_pump)/cells.num_mems
        fy_pump_o = np.dot(cells.M_sum_mems,fy_pump)/cells.num_mems

        # divergence of the total flux:
        gfx_o = (fx_pump_o[cells.cell_nn_i[:,1]] - fx_pump_o[cells.cell_nn_i[:,0]])/cells.nn_len

        fxx = gfx_o*cells.cell_nn_tx

        gfy_o = (fy_pump_o[cells.cell_nn_i[:,1]] - fy_pump_o[cells.cell_nn_i[:,0]])/cells.nn_len
        fyy = gfy_o*cells.cell_nn_ty

        divF_pump = fxx + fyy

        self.rho_pump = self.rho_pump + divF_pump*p.dt

        #------------------------------------------------
        # calculate the total Nernst-Planck flux at each membrane for rho_channel factor:

        fx_chan, fy_chan = stb.nernst_planck_flux(self.rho_channel, gcx_ch, gcy_ch, Ex, Ey,ux_mem,uy_mem,p.D_membrane,
            p.z_channel,self.T,p)

        # map the flux to cell centres:
        fx_chan_o = np.dot(cells.M_sum_mems,fx_chan)/cells.num_mems
        fy_chan_o = np.dot(cells.M_sum_mems,fy_chan)/cells.num_mems

        # divergence of the total flux:
        gfx_o = (fx_chan_o[cells.cell_nn_i[:,1]] - fx_chan_o[cells.cell_nn_i[:,0]])/cells.nn_len

        fxx = gfx_o*cells.cell_nn_tx

        gfy_o = (fy_chan_o[cells.cell_nn_i[:,1]] - fy_chan_o[cells.cell_nn_i[:,0]])/cells.nn_len
        fyy = gfy_o*cells.cell_nn_ty

        divF_chan = fxx + fyy

        self.rho_channel = self.rho_channel + divF_chan*p.dt

        if p.sim_ECM is False:

            # average to the cell centre:
            self.rho_pump_o = np.dot(cells.M_sum_mems, self.rho_pump)/cells.num_mems
            self.rho_channel_o = np.dot(cells.M_sum_mems, self.rho_channel)/cells.num_mems

        #------------------------------------------------
        # make sure nothing is non-zero:
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
            logs.warning('Oops! Morphogen gated ion channel target not found!')
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

    # FIXME consider moving to a pressure module

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
            self.osmo_P_env = fd.integrator(self.osmo_P_env)
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

        # the inflow of mass adds to base pressure in cells
        # (this format is used to avoid conflict with pressure channels):
        if p.run_sim is True:
            self.P_base = self.P_base + P_react[:]

        else:
            self.P_cells = self.P_cells + P_react[:]

        #------------------------------------------------------------------------------------------------------------
        # actual volume change depends on the mechanical properties (young's modulus) of the
        # tissue:

        self.delta_vol = (1/p.youngMod)*div_u_osmo

        # update concentrations and volume in the cell:
        vo = cells.cell_vol[:]

        v1 = (1 + self.delta_vol)*vo

        self.cc_cells =self.cc_cells*(vo/v1)

        # reassign cell volume:
        cells.cell_vol = v1[:]

        if p.sim_ECM is True:
            vo_ecm = cells.ecm_vol[cells.map_cell2ecm]
            v1_ecm = (1 - self.delta_vol)*vo_ecm

            for i, cc_array in enumerate(self.cc_env):
                self.cc_env[i][cells.map_cell2ecm] = cc_array[cells.map_cell2ecm]*(vo_ecm/v1_ecm)

            cells.ecm_vol[cells.map_cell2ecm] = v1_ecm

        if p.voltage_dye is True:

            self.cDye_cell= self.cDye_cell*(vo/v1)

        if p.scheduled_options['IP3'] != 0 or p.Ca_dyn is True:

            self.cIP3 = self.cIP3*(vo/v1)

    def getHydroF(self,cells,p):
        #----Calculate body forces due to hydrostatic pressure gradients---------------------------------------------

        # determine body force due to hydrostatic pressure gradient between cells:

        gPcells = -(self.P_cells[cells.cell_nn_i[:,1]] - self.P_cells[cells.cell_nn_i[:,0]])/cells.nn_len

        self.F_hydro_x_gj = gPcells*cells.cell_nn_tx
        self.F_hydro_y_gj = gPcells*cells.cell_nn_ty

        # calculate a shear electrostatic body force at the cell centre:
        self.F_hydro_x = np.dot(cells.M_sum_mems, self.F_hydro_x_gj)/cells.num_mems
        self.F_hydro_y = np.dot(cells.M_sum_mems, self.F_hydro_y_gj)/cells.num_mems

        self.F_hydro = np.sqrt(self.F_hydro_x ** 2 + self.F_hydro_y ** 2)

    def electro_P(self,cells,p):
        """
        Calculates electrostatic pressure in collection of cells and
        an electrostatic body force.

        """

        # map charge in cell to the membrane, averaging between two cells:
        Q_mem = (self.rho_cells[cells.cell_nn_i[:,1]] + self.rho_cells[cells.cell_nn_i[:,0]])/2

        # calculate force at each membrane:
        self.F_gj_x = Q_mem*self.Egj*cells.cell_nn_tx
        self.F_gj_y = Q_mem*self.Egj*cells.cell_nn_ty

        # calculate a shear electrostatic body force at the cell centre:
        self.F_electro_x = np.dot(cells.M_sum_mems, self.F_gj_x)/cells.num_mems
        self.F_electro_y = np.dot(cells.M_sum_mems, self.F_gj_y)/cells.num_mems

        self.F_electro = np.sqrt(self.F_electro_x**2 + self.F_electro_y**2)

        P_x = (self.F_electro_x*cells.cell_vol)/cells.cell_sa
        P_y = (self.F_electro_y*cells.cell_vol)/cells.cell_sa

        self.P_electro = np.sqrt(P_x**2 + P_y**2)



    # FIXME consider moving to a deformation module

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

            # enforce boundary conditions on u:
            if p.fixed_cluster_bound is True:
                u_x_o[cells.bflags_cells] = 0
                u_y_o[cells.bflags_cells] = 0
                u_x_o[cells.nn_bound] = 0
                u_y_o[cells.nn_bound] = 0

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

        P_react = np.dot(cells.lapGJ_P_inv,2*div_u)  # trial and error suggest this needs to be 2x the divergence?

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


        #-------------------------------------------------------------------------------------------------

        self.dx_time.append(self.d_cells_x[:]) # append the solution to the time-save vector
        self.dy_time.append(self.d_cells_y[:])


        # Initial value solution--------------------------------------------------------------------------------
        if t == 0.0:

            wave_speed = np.sqrt(p.lame_mu/1000)
            wave_speed = np.float(wave_speed)
            wave_speed = np.round(wave_speed,2)

            logs.log_info(
                'Your wave speed is approximately: ' +
                 str(wave_speed) + ' m/s '
            )

            logs.log_info('Try a world size of at least: ' + str(round((5 / 3) * (wave_speed / 500) * 1e6))
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

        P_react = np.dot(cells.lapGJ_P_inv,2*div_u)

        # self.P_cells = (p.lame_mu/k_const)*P_react[:]

        # calculate its gradient:
        gradP_react = (P_react[cells.cell_nn_i[:,1]] - P_react[cells.cell_nn_i[:,0]])/(cells.nn_len)

        gP_x = gradP_react*cells.cell_nn_tx
        gP_y = gradP_react*cells.cell_nn_ty

        # average the components of the reaction force field at cell centres and get boundary values:
        gPx_cell = np.dot(cells.M_sum_mems,gP_x)/cells.num_mems
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
        stb.check_v(self.d_cells_x)

    def get_mass_flux(self,cells,p):
        """
        Sum up individual trans-membrane and
        trans-gap junction ion fluxes to obtain the
        net flow of mass into a cell.

        Assumes that each ion travels with a hydration
        shell of 6 water molecules.

        """

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

    def implement_deform_timestep(self,cells,t,p):
        # Map individual cell deformations to their membranes. In this case,
        # this is better than interpolation.
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

        # Repackage the structure to include individual cell data.
        cells.ecm_verts = [] # null the original ecm verts data structure...

        # Convert region to a numpy array so it can be sorted.
        for i in range(0, len(cells.cell_to_mems)):
            ecm_nest = ecm_new_flat[cells.cell_to_mems[i]]
            ecm_nest = np.asarray(ecm_nest)
            cells.ecm_verts.append(ecm_nest)

        # Voila! Deformed ecm_verts!
        cells.ecm_verts = np.asarray(cells.ecm_verts)
        cells.deformWorld(p)


    # ..................{ PLOTTERS                           }..................

    def _plot_loop(self, cells: 'Cells', p: 'Parameters') -> (np.ndarray, set):
        '''
        Display and/or save an animation during solving if requested _and_
        calculate data common to solving both with and without extracellular
        spaces.

        Returns
        --------
        tt : np.ndarray
            Time-steps vector appropriate for the current run.
        tsamples : set
            Time-steps vector resampled to save data at substantially fewer
            times. The length of this vector governs the number of frames in
            plotted animations, for example.
        '''

        # Human-readable state of extracellular spaces simulation.
        if p.sim_ECM is False:
            ecm_state_label = ''
        else:
            ecm_state_label = '(with extracellular spaces) '

        # Human-readable type and maximum number of steps of the current run.
        if p.run_sim is False:
            figure_type_label = 'Initializing'
            loop_type_label = 'initialization'
            loop_time_step_max = p.init_tsteps
        else:
            figure_type_label = 'Simulating'
            loop_type_label = 'simulation'
            loop_time_step_max = p.sim_tsteps

        # Maximum number of seconds simulated by the current run.
        loop_seconds_max = loop_time_step_max * p.dt

        # Time-steps vector appropriate for the current run.
        tt = np.linspace(0, loop_seconds_max, loop_time_step_max)

        #FIXME: Refactor into a for loop calling the range() builtin. Sunsets!
        # Resample this vector to save data at substantially fewer times.
        tsamples = set()
        i = 0
        while i < len(tt) - p.t_resample:
            i += p.t_resample
            tsamples.add(tt[i])

        # Log this run. # FIXME can we please say "with x total time steps (y sampled)" for brevity? XOXOXXX!
        logs.log_info(
            'Your {} {}is running from 0 to {:.1f} seconds of in-world time '
            'with {} unsampled and {} sampled time steps.'.format(
            loop_type_label,
            ecm_state_label,
            loop_seconds_max,
            len(tt),
            len(tsamples),
        ))

        # If displaying and/or saving an animation during solving, do so.
        if p.plot_while_solving is True:
            self._anim_cells_while_solving = AnimCellsWhileSolving(
                sim=self, cells=cells, p=p,
                type='Vmem',
                figure_title='Vmem while {}'.format(
                    figure_type_label),
                colorbar_title='Voltage [mV]',
                is_color_autoscaled=p.autoscale_Vmem,
                color_min=p.Vmem_min_clr,
                color_max=p.Vmem_max_clr,
            )

        return tt, tsamples

    def _replot_loop(self, p: 'Parameters') -> None:
        '''
        Update the currently displayed and/or saved animation during solving
        with the results of the most recently solved time step, if requested.
        '''

        # Update this animation for the "last frame" if desired, corresponding
        # to the results of the most recently solved time step.
        if p.plot_while_solving is True:
            self._anim_cells_while_solving.plot_frame(frame_number=-1)

    def _deplot_loop(self) -> None:
        '''
        Explicitly close the previously displayed and/or saved animation if
        any _or_ noop otherwise.

        To conserve memory, this method nullifies and hence garbage
        collects both this animation and this animation's associate figure.
        '''

        if self._anim_cells_while_solving is not None:
            self._anim_cells_while_solving.close()
            self._anim_cells_while_solving = None



#-----------------------------------------------------------------------------------------------------------------------
# WASTELANDS
#-----------------------------------------------------------------------------------------------------------------------
#FIXME: I don't quite grok our usage of "sim.run_sim". This undocumented
#attribute appears to be internally set by the Simulator.run_phase_sans_ecm()
#method. That makes sense; however, what's the parallel "p.run_sim" attribute
#for, then?  Interestingly, the "SimRunner" class sets "p.run_sim" as follows:
#
#* To "False" if an initialization is being performed.
#* To "True" if a simulation is being performed.
#
#This doesn't seem quite ideal, however. Ideally, there would exist one and only
#one attribute whose value is an instance of a multi-state "PhaseEnum" class
#rather than two binary boolean attributes. Possible enum values might include:
#
#* "PhaseEnum.seed" when seeding a new cluster.
#* "PhaseEnum.init" when initializing a seeded cluster.
#* "PhaseEnum.sim" when simulation an initialized cluster.
#
#This attribute would probably exist in the "Simulator" class -- say, as
#"sim.phase". In light of that, consider the following refactoring:
#
#* Define a new "PhaseEnum" class in the "sim" module with the above attributes.
#* Define a new "Simulator.phase" attribute initialized to None.
#* Replace all existing uses of the "p.run_sim" and "sim.run_sim" booleans with
#  "sim.phase" instead. Note that only the:
#  * "SimRunner" class sets "p.run_sim".
#  * "Simulator" class sets "sim.run_sim".
#
#Note also the "plot_type" parameter passed to the pipeline.plot_all() function
#*AND* seemingly duplicate "p.plot_type" attribute, which should probably
#receive similar treatment. Wonder temptress at the speed of light and the
#sound of love!


