#!/usr/bin/env python3
# Copyright 2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.


# FIXME test sped up pumps
# FIXME allow time dependent voltage to be added to global boundaries
# FIXME allow user to specify open or closed concentration boundaries
# FIXME create a proper plot while solving -- create an alternative viz module
# FIXME add this grid computation method as feature to params and config and sim runner...
# FIXME add in currents and electric fields
# FIXME complete Ca2+, dye, IP3, etc
# FIXME add in electroosmosis to ecm and gj
# FIXME add in magnetic vector potential and magnetic field calculation -- does it affect E?
# FIXME check cutting holes and other dynamics

import numpy as np
from scipy import interpolate as interp
from scipy import ndimage
import os, os.path
from random import shuffle
from betse.science import filehandling as fh
from betse.science import visualize as viz
from betse.science import toolbox as tb
from betse.science.dynamics import Dynamics
from betse.science.finitediff import FiniteDiffSolver, laplacian, gradient
import matplotlib.pyplot as plt
from betse.exceptions import BetseExceptionSimulation
from betse.util.io import loggers
import time

class SimGrid(object):
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
        self.D_free = []                 # a list of single-valued free diffusion constants for each ion
        self.D_env = []                 # an array of diffusion constants for each ion defined on env grid
        self.movingIons = []            # moving ions indices
        self.ionlabel = {}              # dictionary to hold ion label names
        self.c_env_bound = []           # moving ion concentration at global boundary

        self.T = p.T                # set the base temperature for the simulation

        self.flx_gj_i = np.zeros(len(cells.gj_i))  # flux matrix across gj for individual ions
        self.fluxes_gj = []
        self.I_gj =np.zeros(len(cells.gj_i))     # total current in the gj network

        # Membrane current data structure initialization
        self.flx_mem_i = np.zeros(len(cells.cell_i))
        self.fluxes_mem = []
        self.I_mem =np.zeros(len(cells.cell_i))     # total current across membranes

        self.flx_env_i = np.zeros(len(cells.xypts))
        self.fluxes_env_x = []
        self.fluxes_env_y = []
        self.I_env =np.zeros(len(cells.xypts))     # total current in environment

        ion_names = list(p.ions_dict.keys())

        #-------------------------------------------------------------------------------------------------------------

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

                    setattr(self, str_Dm, np.zeros(len(cells.cell_i)))
                    vars(self)[str_Dm][:] = p.mem_perms[name]

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

            DmH = np.zeros(len(cells.cell_i))
            DmH[:] = p.Dm_H

            self.zH = np.zeros(len(cells.cell_i))
            self.zH[:] = p.z_H

            self.zH2 = np.zeros(len(cells.xypts))
            self.zH2[:] = p.z_H

            # environmental diffusion for each ion
            DenvH = np.zeros(len(cells.xypts))
            DenvH[:] = p.free_diff['H']

            self.cc_cells.append(self.cH_cells)
            self.cc_env.append(self.cH_env)
            self.zs.append(p.z_H)
            self.z_array.append(self.zH)
            self.z_array_env.append(self.zH2)
            self.Dm_cells.append(DmH)
            self.D_env.append(DenvH)
            self.D_free.append(p.Do_H)

            self.fluxes_gj.append(self.flx_gj_i)
            self.fluxes_mem.append(self.flx_mem_i)
            self.fluxes_env_x.append(self.flx_env_i)
            self.fluxes_env_y.append(self.flx_env_i)

        #-------------------------------------------------------------------------------------------------------
        # set up the ecm diffusion constants properly and filter it to avoid artifacts:

        for i, dmat in enumerate(self.D_env):

            self.D_env[i][cells.map_cell2ecm] = dmat[cells.map_cell2ecm]*p.D_ecm_mult

            D2d = np.zeros(cells.X.shape)

            D2d[cells.map_ij2k[:,0],cells.map_ij2k[:,1]] = self.D_env[i]

            D2d = ndimage.filters.gaussian_filter(D2d, 1.0)

            self.D_env[i] = D2d.ravel()

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
        # scheduled dynamics might varry the values
        self.bound_V = {}
        self.bound_V['T'] = 0
        self.bound_V['B'] = 0
        self.bound_V['L'] = 0
        self.bound_V['R'] = 0

    def tissueInit(self,cells,p):

        self.dyna = Dynamics(self,cells,p)   # create the tissue dynamics object
        self.dyna.tissueProfiles(self,cells,p)  # initialize all tissue profiles

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

        self.Dm_vg = np.copy(Dm_cellsA)
        self.Dm_vg[:] = 0

        self.Dm_cag = np.copy(Dm_cellsA)
        self.Dm_cag[:] = 0

        self.Dm_er_base = np.copy(Dm_cellsER)

        self.Dm_er_CICR = np.copy(Dm_cellsER)
        self.Dm_er_CICR[:] = 0

        self.dcc_ER = []

        if p.global_options['gj_block'] != 0:

            self.gj_block = np.ones(len(cells.gj_i))   # initialize the gap junction blocking vector to ones

        else:

            self.gj_block = 1

        if p.scheduled_options['IP3'] != 0 or p.Ca_dyn == True:

            self.cIP3 = np.zeros(len(cells.cell_i))  # initialize a vector to hold IP3 concentrations
            self.cIP3[:] = p.cIP3_to                 # set the initial concentration of IP3 from params file

            self.cIP3_flux_gj = np.zeros(len(cells.gj_i))
            self.cIP3_flux_mem = np.zeros(len(cells.cell_i))

            self.cIP3_flux_gj_time = []
            self.cIP3_flux_mem_time = []

            self.cIP3_env = np.zeros(len(cells.xypts))     # initialize IP3 concentration of the environment
            self.cIP3_env[:] = p.cIP3_to_env

        if p.voltage_dye == True:

            self.cDye_cell = np.zeros(len(cells.cell_i))   # initialize voltage sensitive dye array for cell and env't
            self.cDye_cell[:] = p.cDye_to_cell

            self.Dye_flux_gj = np.zeros(len(cells.gj_i))
            self.Dye_flux_mem = np.zeros(len(cells.cell_i))

            self.Dye_env = np.zeros(len(cells.xypts))     # initialize Dye concentration in the environment
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
        self.cc_env_time = [] # data array holding extracellular concentrations at time points

        self.vm_time = []  # data array holding trans-membrane voltage at time points
        self.vcell_time = []
        self.venv_time = []

        self.dvm_time = []  # data array holding derivative of voltage at time points
        self.time = []     # time values of the simulation

        self.gjopen_time = []   # stores the fractional gap junction open state at each time
        self.fgj_time = []      # stores the gj fluxes for each ion at each time
        self.Igj_time = []      # total current for each gj at each time

        self.I_gj_time = []    #  gap junction current data storage
        self.I_env_time = []   # initialize environmental matrix data storage
        self.I_mem_time = []   # initialize trans-membrane current matrix data storage

        self.cc_er_time = []   # retains er concentrations as a function of time
        self.cIP3_time = []    # retains cellular ip3 concentrations as a function of time

        self.efield_x_time = []   # matrices storing smooth electric field in gj connected cells
        self.efield_y_time = []

        self.efield_env_x_time = []   # matrices storing smooth electric field in ecm
        self.efield_env_y_time = []

        self.I_gjmem_Matrix_x = []
        self.I_gjmem_Matrix_y = []
        self.I_ecm_Matrix_x = []
        self.I_ecm_Matrix_y = []

        if p.voltage_dye == True:

            self.Dye_flux_env_time = []
            self.Dye_flux_gj_time = []
            self.Dye_flux_mem_time = []
            self.cDye_time = []

        # gap junction specific arrays:
        self.id_gj = np.ones(len(cells.gj_i))  # identity array for gap junction indices...
        self.gjopen = np.ones(len(cells.gj_i))   # holds gap junction open fraction for each gj
        self.gjl = np.zeros(len(cells.gj_i))    # gj length for each gj
        self.gjl[:] = cells.gj_len
        self.gjsa = np.zeros(len(cells.gj_i))        # gj x-sec surface area for each gj
        self.gjsa[:] = p.gjsa

        # get the net, unbalanced charge and corresponding voltage in each cell to initialize values of voltages:
        self.update_V(cells,p)

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

            loggers.log_info('Your simulation (quasi-continuous) is running from '+ str(0) + ' to '+ str(round(p.sim_tsteps*p.dt,3))
                         + ' s of in-world time.')

        else:
            loggers.log_info('Your initialization (quasi-continuous) is running from '+ str(0) + ' to '+ str(round(p.init_tsteps*p.dt,3))
                         + ' s of in-world time.')


        # if p.plot_while_solving == True:
        #
        #     self.checkPlot = viz.PlotWhileSolving(cells,self,p,clrAutoscale = p.autoscale_Vmem, clrMin = p.Vmem_min_clr,
        #         clrMax = p.Vmem_max_clr)

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
            fNa_NaK, fK_NaK = pumpNaKATP(self.cc_cells[self.iNa],self.cc_env[self.iNa][cells.map_cell2ecm],
                self.cc_cells[self.iK],self.cc_env[self.iK][cells.map_cell2ecm],cells.cell_sa,self.vm,
                self.T,p,self.NaKATP_block)

            self.fluxes_mem[self.iNa] = fNa_NaK
            self.fluxes_mem[self.iK] = fK_NaK

            # update the concentrations
            self.update_C(self.iNa,fNa_NaK,cells,p)
            self.update_C(self.iK,fK_NaK,cells,p)

            # recalculate the net, unbalanced charge and voltage in each cell:
            self.update_V(cells,p)

            if p.ions_dict['Ca'] == 1:

                _,_, f_CaATP =\
                    pumpCaATP(self.cc_cells[self.iCa],self.cc_env[self.iCa][cells.map_cell2ecm],
                        cells.cell_sa,self.vm,self.T,p)

                # update calcium concentrations in cell and ecm:
                self.update_C(self.iCa,f_CaATP,cells,p)

                self.fluxes_mem[self.iCa] = f_CaATP

                # recalculate the net, unbalanced charge and voltage in each cell:
                self.update_V(cells,p)

                if p.Ca_dyn ==1:

                    self.cc_er[0],self.cc_cells[self.iCa], _ =\
                        pumpCaER(self.cc_er[0],self.cc_cells[self.iCa],cells.cell_sa,p.ER_vol*cells.cell_vol,
                            cells.cell_vol,self.v_er,self.T,p)

                    # recalculate the net, unbalanced charge and voltage in each cell:
                    self.update_V(cells,p)

                    q_er = get_charge_density(self.cc_er,self.z_array_er,p)
                    # self.v_er = get_volt(q_er,p.ER_sa*cells.cell_sa,p)
                    self.v_er = 0  # FIXME this is hacked up

            if p.ions_dict['H'] == 1:

                self.Hplus_electrofuse_ecm(cells,p,t)

                if p.HKATPase_dyn == 1:

                    self.Hplus_HKATP_ecm(cells,p,t)

                if p.VATPase_dyn == 1:

                    self.Hplus_VATP_ecm(cells,p,t)

            # electro-diffuse all ions (except for proteins, which don't move) across the cell membrane:
            shuffle(cells.gj_i)
            shuffle(self.movingIons)

            for i in self.movingIons:

                f_ED = electroflux(self.cc_env[i][cells.map_cell2ecm],self.cc_cells[i],
                         self.Dm_cells[i], p.tm, cells.cell_sa, self.zs[i], self.vm, self.T, p)
                #
                #
                self.fluxes_mem[i] = self.fluxes_mem[i] + f_ED
                #
                # # update ion concentrations in cell and ecm:
                self.update_C(i,f_ED,cells,p)

                # transport in the extracellular environment:
                self.update_ENV(cells,p,i)

                # update flux between cells due to gap junctions
                self.update_GJ(cells,p,i)

            # # recalculate the net, unbalanced charge and voltage in each cell:
            self.update_V(cells,p)

            if p.scheduled_options['IP3'] != 0 or p.Ca_dyn == True:

                self.update_IP3(cells,p)

            if p.Ca_dyn == 1 and p.ions_dict['Ca'] == 1:

               self.update_ER(cells,p)

            # if p.voltage_dye=1 electrodiffuse voltage sensitive dye between cell and environment
            if p.voltage_dye ==1:

                self.update_Dye(cells,p)

            if p.dynamic_noise == 1 and p.ions_dict['P']==1:
                # add a random walk on protein concentration to generate dynamic noise:
                self.protein_noise_factor = p.dynamic_noise_level*(np.random.random(len(cells.cell_i)) - 0.5)
                self.cc_cells[self.iP] = self.cc_cells[self.iP]*(1+ self.protein_noise_factor)

                # recalculate the net, unbalanced charge and voltage in each cell:
                self.update_V(cells,p)

            check_v(self.vm)


            if t in tsamples:
                # #
                # self.get_current(cells,p)   # get the current in the gj network connection of cells

                # add the new concentration and voltage data to the time-storage matrices:

                # efieldx = self.E_CELL_x[:]
                # self.efield_x_time.append(efieldx)
                # efieldx = None
                #
                # efieldy = self.E_CELL_y[:]
                # self.efield_y_time.append(efieldy)
                # efieldy = None
                #
                # efield_ecm_x = self.E_ECM_x[:]
                # self.efield_env_x_time.append(efield_ecm_x)
                # efield_ecm_x = None
                #
                # efield_ecm_y = self.E_ECM_y[:]
                # self.efield_env_y_time.append(efield_ecm_y)
                # efield_ecm_y = None

                concs = self.cc_cells[:]
                self.cc_time.append(concs)
                concs = None

                ecmsc = self.cc_env[:]
                self.cc_env_time.append(ecmsc)
                ecmsc = None

                vmm = self.vm[:]
                self.vm_time.append(vmm)
                vmm = None

                dvmm = self.dvm[:]
                self.dvm_time.append(dvmm)
                dvmm = None

                ggjopen = self.gjopen[:]
                self.gjopen_time.append(ggjopen)
                ggjopen = None

                vvcell = self.v_cell[:]
                self.vcell_time.append(vvcell)
                vvcell = None

                vvecm = self.v_env[:]
                self.venv_time.append(vvecm)
                vvecm = None

                # Igj = self.I_gj[:]
                # self.I_gj_time.append(Igj)
                # Igj = None
                #
                # Imem = self.I_mem[:]
                # self.I_mem_time.append(Imem)
                # Imem = None

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

                # if p.plot_while_solving == True:
                #     self.checkPlot.updatePlot(self,p)

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

        # self.checkPlot = None

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

        plt.close()
        loggers.log_info('Simulation completed successfully.')

    def update_V(self,cells,p):

        """
        Calculates net charge density and voltages in the cell and local extracellular spaces
        and calculates a trans-membrane voltage self.vm.

        """
        self.rho_cells = get_charge_density(self.cc_cells, self.z_array, p)
        self.rho_env = get_charge_density(self.cc_env, self.z_array_env, p)
        self.v_cell = get_Vcell(self.rho_cells,cells,p)
        self.v_env = get_Venv(self,cells,p)

        self.vm = self.v_cell - self.v_env[cells.map_cell2ecm]  # calculate v_mem

    def update_C(self,ion_i,flux,cells,p):

        self.cc_cells[ion_i] = self.cc_cells[ion_i] + flux*(cells.cell_sa/cells.cell_vol)*p.dt

        self.cc_env[ion_i][cells.map_cell2ecm] = self.cc_env[ion_i][cells.map_cell2ecm] - \
                                                 flux*(cells.cell_sa/cells.ecm_vol)*p.dt

    def update_GJ(self,cells,p,i):
         # calculate voltage difference between cells:
        vmA,vmB = self.v_cell[cells.gap_jun_i][:,0], self.v_cell[cells.gap_jun_i][:,1]
        self.vgj = vmB - vmA

        if p.v_sensitive_gj == True:

            # determine the open state of gap junctions:
            self.gjopen = self.gj_block*((1.0 - tb.step(abs(self.vgj),p.gj_vthresh,p.gj_vgrad) + 0.1))

        # determine flux through gap junctions for this ion:
        fgj = electroflux(self.cc_cells[i][cells.gap_jun_i][:,0],self.cc_cells[i][cells.gap_jun_i][:,1],
            self.D_free[i],cells.gj_len,self.gjopen*p.gjsa,self.zs[i],self.vgj,self.T,p)

        # update cell concentration due to gap junction flux:
        deltac = fgj*p.dt*((self.gjopen*p.gjsa)/cells.cell_vol)
        self.cc_cells[i] = self.cc_cells[i] + np.dot(deltac, cells.gjMatrix)

        # recalculate the net, unbalanced charge and voltage in each cell:
        self.update_V(cells,p)

        self.fluxes_gj[i] = fgj  # store gap junction flux for this ion

    def update_ENV(self,cells,p,i):

        # enforce boundary conditions on environmental voltage: FIXME allow this to be variable -- closed bound, open, applied V
        # self.v_env[cells.bL_k] = 0
        # self.v_env[cells.bR_k] = 0
        # self.v_env[cells.bTop_k] = 0
        # self.v_env[cells.bBot_k] = 0

        # open boundary
        # self.cc_env[i][cells.bL_k] = self.c_env_bound[i]
        # self.cc_env[i][cells.bR_k] = self.c_env_bound[i]
        # self.cc_env[i][cells.bTop_k] = self.c_env_bound[i]
        # self.cc_env[i][cells.bBot_k] = self.c_env_bound[i]

        # make v_env and cc_env into 2d matrices
        cenv = self.cc_env[i][:]
        denv = self.D_env[i][:]

        self.v_env = self.v_env.reshape(cells.X.shape)
        # self.v_env = ndimage.filters.gaussian_filter(self.v_env, 1.0)
        # self.rho_env = self.rho_env.reshape(cells.X.shape)

        cenv = cenv.reshape(cells.X.shape)
        # cenv = ndimage.filters.gaussian_filter(cenv, 1.0)

        denv = denv.reshape(cells.X.shape)

        # calculate gradients in the environment
        grad_V_env_x, grad_V_env_y = gradient(self.v_env, cells.delta)
        grad_cc_env_x, grad_cc_env_y = gradient(cenv, cells.delta)
        grad_Dx, grad_Dy = gradient(denv, cells.delta)

        # calculate laplacians for data in the environment
        dd_cc_env = laplacian(cenv, cells.delta)
        dd_V_env = laplacian(self.v_env,cells.delta)

        f_env_x, f_env_y = nernst_planck_flux(cenv,grad_cc_env_x,grad_cc_env_y,
            grad_V_env_x, grad_V_env_y, denv,self.zs[i],self.T,p)

        delta_c = nernst_planck(cenv,grad_cc_env_x,grad_cc_env_y,dd_cc_env,grad_V_env_x,grad_V_env_y,
            dd_V_env,grad_Dx,grad_Dy,denv,self.zs[i],self.T,p)

        cenv = cenv + delta_c*p.dt

        # Neumann boundary condition (flux at boundary)
        # zero flux boundaries for concentration
        cenv[:,-1] = cenv[:,-2]
        cenv[:,0] = cenv[:,1]
        cenv[0,:] = cenv[1,:]
        cenv[-1,:] = cenv[-2,:]

        # reshape the matrices into vectors
        # self.rho_env = self.rho_env.ravel()
        self.v_env = self.v_env.ravel()
        self.cc_env[i] = cenv.ravel()

         # recalculate the net, unbalanced charge and voltage in each ecm space:
        # self.rho_env = get_charge_density(self.cc_env,self.z_array_env,p)
        # self.v_env = get_Venv(self,cells,p)

        self.fluxes_env_x[i] = f_env_x.ravel()  # store ecm junction flux for this ion
        self.fluxes_env_y[i] = f_env_y.ravel()  # store ecm junction flux for this ion

    def update_ER(self,cells,p):

         # electrodiffusion of ions between cell and endoplasmic reticulum
        self.cc_cells[self.iCa],self.cc_er[0],_ = \
        electroflux(self.cc_cells[self.iCa],self.cc_er[0],self.Dm_er[0],self.tm,p.ER_sa*cells.cell_sa,
            cells.cell_vol,p.ER_vol*cells.cell_vol,self.z_er[0],self.v_er,self.T,p)

        # Electrodiffusion of charge compensation anion
        self.cc_cells[self.iM],self.cc_er[1],_ = \
        electroflux(self.cc_cells[self.iM],self.cc_er[1],self.Dm_er[1],self.tm,p.ER_sa*cells.cell_sa,
            cells.cell_vol,p.ER_vol*cells.cell_vol,self.z_er[1],self.v_er,self.T,p)

        # recalculate the net, unbalanced charge and voltage in each cell:
        self.update_V(cells,p)

        q_er = get_charge_density(self.cc_er,self.z_array_er,p)
        self.v_er = get_volt(q_er,p.ER_sa*cells.cell_sa,p)

    def update_DYE(self,cells,p):


        flux_dye = electroflux(self.cDye_ecm[cells.mem_to_ecm],self.cDye_cell[cells.mem_to_cells],
                            p.Dm_Dye,self.tm[cells.mem_to_cells],cells.mem_sa,
                            p.z_Dye,self.vm,self.T,p)

         # update the dye concentrations in the cell and ecm due to ED fluxes at membrane
        self.cDye_cell = self.cDye_cell + \
                            np.dot((flux_dye/cells.cell_vol[cells.mem_to_cells])*p.dt,cells.cell_UpdateMatrix)

        self.cDye_ecm = self.cDye_ecm - \
                            np.dot((flux_dye/cells.ecm_vol[cells.mem_to_ecm])*p.dt,cells.ecm_UpdateMatrix)

        fDye_gj = electroflux(self.cDye_cell[cells.gap_jun_i][:,0],self.cDye_cell[cells.gap_jun_i][:,1],
            p.Do_Dye,self.gjl,self.gjopen*self.gjsa,p.z_Dye,self.vgj,self.T,p)

        # update cell voltage-sensitive dye concentration due to gap junction flux:
        self.cDye_cell = (self.cDye_cell*cells.cell_vol + np.dot((fDye_gj*p.dt), cells.gjMatrix))/cells.cell_vol

        # electrodiffuse dye through ecm <---> ecm junctions
        flux_ecm_dye = electroflux(self.cDye_ecm[cells.ecm_nn_i[:,0]],self.cDye_ecm[cells.ecm_nn_i[:,1]],
                p.Do_Dye,cells.len_ecm_junc,self.ec2ec_sa,p.z_Dye,self.v_ec2ec,self.T,p)
                    #
        self.cDye_ecm = (self.cDye_ecm*cells.ecm_vol + np.dot(flux_ecm_dye*p.dt,cells.ecmMatrix))/cells.ecm_vol

    def update_IP3(self,cells,p):
         # determine flux through gap junctions for IP3:
        # _,_,fIP3 = electrofuse(self.cIP3[cells.gap_jun_i][:,0],self.cIP3[cells.gap_jun_i][:,1],
        #     self.id_gj*p.Do_IP3,self.gjl,self.gjopen*self.gjsa,cells.cell_vol[cells.gap_jun_i][:,0],
        #     cells.cell_vol[cells.gap_jun_i][:,1],p.z_IP3,self.vgj,self.T,p)
        fIP3 = electroflux(self.cIP3[cells.gap_jun_i][:,0],self.cIP3[cells.gap_jun_i][:,1], p.Do_IP3,self.gjl,
            self.gjopen*self.gjsa,p.z_IP3,self.vgj,self.T,p)

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

    def get_current(self,cells,p):

        self.I_gj = np.zeros(len(cells.gj_i))
        self.I_mem = np.zeros(len(cells.mem_i))

        # calculate current across gap junctions:
        for flux_array, zi in zip(self.fluxes_gj,self.zs):

            # I_i = (flux_array*zi*p.F)/(self.gjopen*self.gjsa)
            I_i = flux_array*zi*p.F

            self.I_gj = self.I_gj + I_i

        # calculate current across cell membranes:
        for flux_array, zi in zip(self.fluxes_mem,self.zs):

            # I_i = (flux_array*zi*p.F)/(self.gjopen*self.gjsa)
            I_i = flux_array*zi*p.F

            self.I_mem = self.I_mem + I_i

        # calculate and store currents components:
        Jmag = np.hstack((self.I_gj,self.I_mem))

        jx = Jmag*cells.nx_Igj
        jy = Jmag*cells.ny_Igj

        self.X_Igj,self.Y_Igj,J_x,J_y = tb.grid_vector_data(cells.xpts_Igj,cells.ypts_Igj,jx,jy,cells,p)

        J_x = np.nan_to_num(J_x)
        J_y = np.nan_to_num(J_y)

        self.I_gjmem_Matrix_x.append(J_x)
        self.I_gjmem_Matrix_y.append(J_y)

        if p.sim_ECM == True:

            self.I_ecm = np.zeros(len(cells.ecm_nn_i))
            self.I_env = np.zeros(len(cells.env_i))

            for flux_array, zi in zip(self.fluxes_ecm,self.zs):

                # I_i = (flux_array*zi*p.F)/(self.gjopen*self.gjsa)
                I_i = flux_array*zi*p.F

                self.I_ecm = self.I_ecm + I_i

            for flux_array, zi in zip(self.fluxes_env,self.zs):

                I_i = flux_array*zi*p.F

                self.I_env = self.I_env + I_i

            # calculate interpolated verts and midpoint data for currents:

            Jmag_ecm = self.I_ecm

            jx_ecm = Jmag_ecm*cells.nx_Iecm
            jy_ecm = Jmag_ecm*cells.ny_Iecm

            # data on environmental currents
            Jmag_env = self.I_env

            # environmental <---> boundary ecm current components:
            jx_env = cells.nx_Ienv*Jmag_env
            jy_env = cells.ny_Ienv*Jmag_env

            # update ecm currents vector components to include environmental current:
            jx_ecm[cells.bflags_ecm] = jx_ecm[cells.bflags_ecm] + jx_env
            jy_ecm[cells.bflags_ecm] = jy_ecm[cells.bflags_ecm] + jy_env

            self.X_Iecm,self.Y_Iecm,J_x,J_y = tb.grid_vector_data(cells.xpts_Iecm,cells.ypts_Iecm,jx_ecm,jy_ecm,cells,p)

            J_x = np.nan_to_num(J_x)
            J_y = np.nan_to_num(J_y)

            self.I_ecm_Matrix_x.append(J_x)
            self.I_ecm_Matrix_y.append(J_y)

    def Hplus_electrofuse_ecm(self,cells,p,t):

        # electrofuse the H+ ion between the cytoplasm and the ecms
        _,_,f_H1 = \
            electroflux(self.cc_ecm[self.iH][cells.mem_to_ecm],self.cc_cells[self.iH][cells.mem_to_cells],
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

def get_Vcell(rho_cell,cells,p):

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

    v_cell = rho_cell*(p.rc**2)/(3*p.eo*80.0)

    return v_cell

def get_Venv(self,cells,p):   # FIXME may need to make this a Poisson solver!

        """
        Calculates the voltage in each extracellular space from Poisson equation charge density

        Parameters
        ---------------
        rho_ecm:        an array listing charge density in each ecm space [C/m3]
        p:              an instance of the Parameters object

        Returns
        ---------
        v_ecm           the voltage in each ecm space   [V]

        """

        # rec = p.rc + p.cell_space
        # v_ecm = (rho_ecm*(rec**3 - p.rc**3))/(3*p.eo*80.0*rec)

        # self.v_env = (self.rho_env*(p.cell_space**2))/(8*p.eo*80.0)

        # Poisson solver----------------------------------------------------------------

        rho_env = self.rho_env[cells.core_k][:]

        # modify the RHS of the equation to incorporate Dirichlet boundary conditions:

        rho_env[cells.bBot_kp] = rho_env[cells.bBot_kp] - (self.bound_V['B']/cells.delta**2)
        rho_env[cells.bTop_kp] = rho_env[cells.bTop_kp] - (self.bound_V['T']/cells.delta**2)
        rho_env[cells.bL_kp] = rho_env[cells.bL_kp] - (self.bound_V['L']/cells.delta**2)
        rho_env[cells.bR_kp] = rho_env[cells.bR_kp] - (self.bound_V['R']/cells.delta**2)

        # Solve Poisson's electrostatic equation:
        V = np.dot(cells.Ainv,-rho_env)

        # Re-pack the solution to the original vector:
        self.v_env[cells.core_k] = V

        # set boundary values on the original vector
        self.v_env[cells.bTop_k] = self.bound_V['T']
        self.v_env[cells.bBot_k] = self.bound_V['B']
        self.v_env[cells.bL_k] = self.bound_V['L']
        self.v_env[cells.bR_k] = self.bound_V['R']

        return self.v_env

def electroflux(cA,cB,Dc,d,sa,zc,vBA,T,p):
    """
    Electro-diffusion between two connected volumes. Note for cell work, 'b' is 'inside', 'a' is outside, with
    a positive flux moving from a to b. The voltage is defined as
    Vb - Va (Vba), which is equivalent to Vmem.

    This function defaults to regular diffusion if Vba == 0.0

    This function takes numpy matrix values as input. All inputs must be matrices of
    the same shape.

    Parameters
    ----------
    cA          concentration in region A [mol/m3] (out)
    cB          concentration in region B [mol/m3] (in)
    Dc          Diffusion constant of c  [m2/s]
    d           Distance between region A and region B [m]
    sa          Surface area separating region A and B [m2]
    zc          valence of ionic species c
    vBA         voltage difference between region B (in) and A (out) = Vmem
    p           an instance of the Parameters class


    Returns
    --------
    flux        Chemical flux magnitude between region A and B [mol/s]

    """

    grad_c = (cB - cA)/d

    c_mid = (cA + cB)/2

    grad_V = vBA/d

    flux = -Dc*grad_c - ((Dc*p.q*zc)/(p.kb*T))*c_mid*grad_V
    # flux = sa*flux

    return flux

def nernst_planck_flux(c, gcx, gcy, gvx, gvy,D,z,T,p):

    alpha = (D*z*p.q)/(p.kb*T)

    fx =  -D*gcx - alpha*gvx*c

    fy =  -D*gcy - alpha*gvy*c

    return fx, fy

def nernst_planck(c,gcx,gcy,ddc,gvx,gvy,ddv,gdx,gdy,D,z,T,p):

    alpha = (D*z*p.q)/(p.kb*T)
    beta = (c*z*p.q)/(p.kb*T)

    delta_c = (gdx*gcx + gdy*gcy) + D*ddc + alpha*(gcx*gvx + gcy*gvy) + beta*(gdx*gvx + gdy*gvy) + c*alpha*ddv

    return delta_c

def pumpNaKATP(cNai,cNao,cKi,cKo,sa,Vm,T,p,block):

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

def pumpCaATP(cCai,cCao,sa,Vm,T,p):

    """
    Parameters
    ----------
    cCai            Concentration of Ca2+ inside the cell
    cCao            Concentration of Ca2+ outside the cell
    sa              Surface area over which flux transfer occurs
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

        alpha = sa*p.alpha_Ca*tb.step(delG,p.halfmax_Ca,p.slope_Ca)

        f_Ca  = -alpha*(cCai)      #flux as [mol/s], scaled to concentration in cell

    elif p.backward_pumps == True:

        alpha = sa*signG*p.alpha_Ca*tb.step(delG,p.halfmax_Ca,p.slope_Ca)

        truth_forwards = signG == 1
        truth_backwards = signG == -1

        inds_forwards = (truth_forwards).nonzero()  # indices of forward-running cells
        inds_backwards = (truth_backwards).nonzero() # indices of backward-running cells

        f_Ca = np.zeros(len(cCai))

        f_Ca[inds_forwards]  = -alpha*(cCai)      #flux as [mol/s], scaled to concentration in cell
        f_Ca[inds_backwards]  = -alpha*(cCao)      #flux as [mol/s], scaled to concentration out of cell

    return f_Ca

def pumpCaER(cCai,cCao,sa,Vm,T,p):


    alpha = sa*p.alpha_CaER

    f_Ca  = alpha*(cCao)*(1.0 - cCai)      #flux as [mol/s]


    return f_Ca

def pumpHKATP(cHi,cHo,cKi,cKo,sa,Vm,T,p,block):

    """
    Parameters
    ----------
    cHi            Concentration of H+ inside the cell
    cHo            Concentration of H+ outside the cell
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

    delG_H = p.R*T*np.log(cHo/cHi) - p.F*Vm
    delG_K = p.R*T*np.log(cKi/cKo) + p.F*Vm

    delG_HKATP = deltaGATP - (delG_H + delG_K)
    delG_pump = (delG_HKATP/1000)
    delG = np.absolute(delG_pump)
    signG = np.sign(delG)

    if p.backward_pumps == False:

        alpha = sa*block*p.alpha_HK*tb.step(delG,p.halfmax_HK,p.slope_HK)
        f_H  = -alpha*cHi*cKo      #flux as [mol/s], scaled by concentrations in and out

    elif p.backward_pumps == True:

        alpha = sa*signG*block*p.alpha_HK*tb.step(delG,p.halfmax_HK,p.slope_HK)

        truth_forwards = signG == 1
        truth_backwards = signG == -1

        inds_forwards = (truth_forwards).nonzero()  # indices of forward-running cells
        inds_backwards = (truth_backwards).nonzero() # indices of backward-running cells

        f_H = np.zeros(len(cHi))

        f_H[inds_forwards]  = -alpha*cHi*cKo      #flux as [mol/s], scaled by concentrations in and out
        f_H[inds_backwards]  = -alpha*cHo*cKi

    f_K = -f_H          # flux as [mol/s]


    return f_H, f_K

def pumpVATP(cHi,cHo,sa,Vm,T,p):

    deltaGATP = 20*p.R*T

    delG_H = p.R*T*np.log(cHo/cHi) - p.F*Vm  # free energy to move H+ out of cell

    delG_VATP = deltaGATP - delG_H   # free energy available to move H+ out of cell
    delG_pump = (delG_VATP/1000)
    delG = np.absolute(delG_pump)
    signG = np.sign(delG)

    if p.backward_pumps == False:

        alpha = sa*p.alpha_V*tb.step(delG,p.halfmax_V,p.slope_V)
        f_H  = -alpha*cHi      #flux as [mol/s], scaled by concentrations in and out

    elif p.backward_pumps == True:

        alpha = sa*signG*p.alpha_V*tb.step(delG,p.halfmax_V,p.slope_V)

        truth_forwards = signG == 1
        truth_backwards = signG == -1

        inds_forwards = (truth_forwards).nonzero()  # indices of forward-running cells
        inds_backwards = (truth_backwards).nonzero() # indices of backward-running cells

        f_H = np.zeros(len(cHi))

        f_H[inds_forwards]  = -alpha*cHi      #flux as [mol/s], scaled by concentrations in and out
        f_H[inds_backwards]  = -alpha*cHo


    return f_H

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

