#!/usr/bin/env python3
# Copyright 2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.


# FIXME implement stability safety threshhold parameter checks and loss-of-stability detection + error message
# FIXME would be nice to have a time estimate for the simulation
# FIXME Initialization function of Dmems, gjs, ect for simulation
# FIXME combine effects of scheduling and dynamic gating into one function so effects can overlap.
# FIXME Calcium dynamics
# FIXME ECM diffusion and discrete membrane domains?
 # FIXME would be nice to track ATP use

import numpy as np
import os, os.path
import copy
from random import shuffle
from betse.science import filehandling as fh


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
        betse_cache_dir = os.path.expanduser(p.cache_path)
        os.makedirs(betse_cache_dir, exist_ok=True)

        # Define data paths for saving an initialization and simulation run:
        self.savedInit = os.path.join(betse_cache_dir, 'saved_init.btse')
        self.savedSim = os.path.join(betse_cache_dir, 'saved_sim.btse')

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
        self.zs = []   # ion valence state initialized
        self.Dm_cells = []              # membrane diffusion constants initialized
        self.D_free = []                 # a list of single-valued free diffusion constants for each ion
        self.movingIons = []            # moving ions indices
        self.ionlabel = {}              # dictionary to hold ion label names

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
            self.ionlabel[self.iCl] = 'chlorine'

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

        if p.ions_dict['H'] == 1:

            i =i+1

            self.iH = i
            self.movingIons.append(self.iH)
            self.ionlabel[self.iH] = 'protons'

            cH_cells = np.zeros(len(cells.cell_i))
            cH_cells[:]=p.cH_cell

            cH_env = np.zeros(len(cells.cell_i))
            cH_env[:]=p.cH_env

            DmH = np.zeros(len(cells.cell_i))
            DmH[:] = p.Dm_H

            self.cc_cells.append(cH_cells)
            self.cc_env.append(cH_env)
            self.zs.append(p.z_H)
            self.Dm_cells.append(DmH)
            self.D_free.append(p.Do_H)

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
            self.Dm_cells.append(DmP)        # FIXME is it a problem that this is here?
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

        self.vm_to = np.zeros(len(cells.cell_i))

        print(self.ionlabel)

    def tissueInit(self,cells,p):


        # Initialize an array structure that will hold user-scheduled changes to membrane permeabilities:
        Dm_cellsA = np.asarray(self.Dm_cells)

        self.Dm_base = np.copy(Dm_cellsA) # make a copy that will serve as the unaffected values base

        self.Dm_scheduled = np.copy(Dm_cellsA)
        self.Dm_scheduled[:] = 0

        # Initialize an array structure that will hold dynamic voltage-gated channel changes to mem permeability:
        self.Dm_vg = np.copy(Dm_cellsA)
        self.Dm_vg[:] = 0

        if p.vg_options['Na_vg'] != 0:

            # Initialization of logic values for voltage gated sodium channel
            self.maxDmNa = p.vg_options['Na_vg'][0]
            self.v_on_Na = p.vg_options['Na_vg'][1]
            self.v_off_Na = p.vg_options['Na_vg'][2]
            self.v_reactivate_Na = p.vg_options['Na_vg'][3]
            self.t_alive_Na = p.vg_options['Na_vg'][4]

            # Initialize matrices defining states of vgNa channels for each cell:
            self.active_Na = np.zeros(len(cells.cell_i))
            self.crossed_inactivate_Na = np.zeros(len(cells.cell_i))
            self.vgNa_state = np.zeros(len(cells.cell_i))

            self.vgNa_OFFtime = np.zeros(len(cells.cell_i)) # sim time at which vgK starts to close

        if p.vg_options['K_vg']  !=0:

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
            self.scheduled_target_inds = np.zeros(len(cells.cell_i))

        elif p.scheduled_targets == 'all':
            self.scheduled_target_inds = np.ones(len(cells.cell_i))

        elif p.scheduled_targets =='random1':
            shuffle(cells.cell_i)
            trgt2 = cells.cell_i[0]
            self.scheduled_target_inds = np.zeros(len(cells.cell_i))
            self.scheduled_target_inds[trgt2] = 1

        elif p.scheduled_targets == 'random50':
            shuffle(cells.cell_i)
            self.scheduled_target_inds = np.random.random(len(cells.cell_i))
            self.scheduled_target_inds = np.rint(self.scheduled_target_inds)

        elif isinstance(p.scheduled_targets, list):
            self.scheduled_target_inds = np.zeros(len(cells.cell_i))
            self.scheduled_target_inds[p.scheduled_targets] = 1

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
        self.cc_time = []  # data array holding the concentrations at time points
        self.vm_time = []  # data array holding voltage at time points
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
        print('Your sim initialization is running for', round((p.init_tsteps*p.dt)/60,2),'minutes of in-world time.')

        for t in tt:   # run through the loop

            # get the net, unbalanced charge in each cell:
            q_cells = get_charge(self.cc_cells,self.zs,cells.cell_vol,p)

            # calculate the voltage in the cell (which is also Vmem as environment is zero):
            self.vm = get_volt(q_cells,cells.cell_sa,p)

            # run the Na-K-ATPase pump:
            self.cc_cells[self.iNa],self.cc_env[self.iNa],self.cc_cells[self.iK],self.cc_env[self.iK], fNa_NaK, fK_NaK =\
                pumpNaKATP(self.cc_cells[self.iNa],self.cc_env[self.iNa],self.cc_cells[self.iK],self.cc_env[self.iK],
                    cells.cell_vol,self.envV,self.vm,p)

            if p.ions_dict['Ca'] == 1:

                self.cc_cells[self.iCa],self.cc_env[self.iCa], _ =\
                    pumpCaATP(self.cc_cells[self.iCa],self.cc_env[self.iCa],cells.cell_vol,self.envV,self.vm,p)

            # electro-diffuse all ions (except for proteins, which don't move!) across the cell membrane:
            shuffle(self.movingIons)  # shuffle the ion indices so it's not the same order every time step

            for i in self.movingIons:

                # recalculate the net, unbalanced charge and voltage in each cell:
                q_cells = get_charge(self.cc_cells,self.zs,cells.cell_vol,p)
                self.vm = get_volt(q_cells,cells.cell_sa,p)

                self.cc_env[i],self.cc_cells[i],fNa = \
                    electrofuse(self.cc_env[i],self.cc_cells[i],self.Dm_cells[i],self.tm,cells.cell_sa,
                        self.envV,cells.cell_vol,self.zs[i],self.vm,p)

            if t in tsamples:
                # add the new concentration and voltage data to the time-storage matrices:
                concs = copy.deepcopy(self.cc_cells)
                self.cc_time.append(concs)
                vmm = copy.deepcopy(self.vm)
                self.vm_time.append(vmm)
                self.time.append(t)

        celf = copy.deepcopy(self)

        datadump = [celf,cells,p]
        fh.saveSim(self.savedInit,datadump)
        message_1 = 'Initialization run saved to' + ' ' + p.cache_path
        print(message_1)

        self.vm_to = copy.deepcopy(self.vm)

        for i in range(0,len(self.ionlabel)):
            endconc = np.round(np.mean(self.cc_time[-1][i]),2)
            label = self.ionlabel[i]
            concmess = 'Final cytoplasmic concentration of'+ ' '+ label + ': '
            print(concmess,endconc,' mmol/L')


        final_vmean = 1000*np.round(np.mean(self.vm_time[-1]),2)
        vmess = 'Final cell Vmem of ' + ': '
        print(vmess,final_vmean, ' mV')

    def runSim(self,cells,p,save=None):
        """
        Drives the actual time-loop iterations for the simulation.
        """

        self.tissueInit(cells,p)   # Initialize all structures used for gap junctions, ion channels, and other dynamics

        # Reinitialize all time-data structures
        self.cc_time = []  # data array holding the concentrations at time points
        self.envcc_time = [] # data array holding environmental concentrations at time points
        self.vm_time = []  # data array holding voltage at time points
        self.dvm_time = []  # data array holding derivative of voltage at time points
        self.time = []     # time values of the simulation

        self.fgj_time = []      # stores the gj fluxes for each ion at each time
        self.Igj_time = []      # current for each gj at each time
        self.vgj_time = []

        self.active_Na_time = []   # stores the activation state of Na and K voltage gated channels at time
        self.active_K_time = []

        # gap junction specific arrays:
        self.id_gj = np.ones(len(cells.gj_i))
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
        print('Your simulation is running from',0,'to',round(p.sim_tsteps*p.dt,3),'seconds of in-world time.')


        for t in tt:   # run through the loop

            # get the net, unbalanced charge in each cell:
            q_cells = get_charge(self.cc_cells,self.zs,cells.cell_vol,p)

            # calculate the voltage in the cell (which is also Vmem as environment is zero):
            self.vm = get_volt(q_cells,cells.cell_sa,p)

            # run the Na-K-ATPase pump:
            self.cc_cells[self.iNa],self.cc_env[self.iNa],self.cc_cells[self.iK],self.cc_env[self.iK], fNa_NaK, fK_NaK =\
                pumpNaKATP(self.cc_cells[self.iNa],self.cc_env[self.iNa],self.cc_cells[self.iK],self.cc_env[self.iK],
                    cells.cell_vol,self.envV,self.vm,p)

            if p.ions_dict['Ca'] == 1:

                self.cc_cells[self.iCa],self.cc_env[self.iCa], _ =\
                    pumpCaATP(self.cc_cells[self.iCa],self.cc_env[self.iCa],cells.cell_vol,self.envV,self.vm,p)

            # electro-diffuse all ions (except for proteins, which don't move!) across the cell membrane:

            shuffle(cells.gj_i)
            shuffle(self.movingIons)

            for i in self.movingIons:

                # recalculate the net, unbalanced charge and voltage in each cell:
                q_cells = get_charge(self.cc_cells,self.zs,cells.cell_vol,p)
                self.vm = get_volt(q_cells,cells.cell_sa,p)

                self.dvm = (self.vm - self.vm_to)/p.dt    # calculate the change in the voltage derivative
                self.vm_to = copy.deepcopy(self.vm)       # reassign the history-saving vm

                # calculate the values of ion channel multipliers:

                self.allDynamics(t,p)  # user-scheduled (forced) interventions

                # electrodiffusion of ion between cell and extracellular matrix
                self.cc_env[i],self.cc_cells[i],_ = \
                    electrofuse(self.cc_env[i],self.cc_cells[i],self.Dm_cells[i],self.tm,cells.cell_sa,
                        self.envV,cells.cell_vol,self.zs[i],self.vm,p)

                # recalculate the net, unbalanced charge and voltage in each cell:
                q_cells = get_charge(self.cc_cells,self.zs,cells.cell_vol,p)
                self.vm = get_volt(q_cells,cells.cell_sa,p)

                # calculate volatge difference between cells:
                vmA,vmB = self.vm[cells.gap_jun_i][:,0], self.vm[cells.gap_jun_i][:,1]
                vgj = vmB - vmA

                # determine the open state of gap junctions:
                self.gjopen = (1.0 - step(abs(vgj),p.gj_vthresh,p.gj_vgrad)) +0.1

                # determine flux through gap junctions for this ion:
                _,_,fgj = electrofuse(self.cc_cells[i][cells.gap_jun_i][:,0],self.cc_cells[i][cells.gap_jun_i][:,1],
                    self.id_gj*self.D_free[i],self.gjl,self.gjopen*self.gjsa,cells.cell_vol[cells.gap_jun_i][:,0],
                    cells.cell_vol[cells.gap_jun_i][:,1],self.zs[i],vgj,p)

                # update cell concentration due to gap junction flux:
                #mole_delta = (fgj*p.dt)

                self.cc_cells[i] = (self.cc_cells[i]*cells.cell_vol + np.dot((fgj*p.dt), cells.gjMatrix))/cells.cell_vol

                self.fluxes_gj[i] = fgj  # store gap junction flux for this ion

                # self.cc_cells[i] = check_c(self.cc_cells[i])

            if t in tsamples:
                # add the new concentration and voltage data to the time-storage matrices:
                concs = copy.deepcopy(self.cc_cells)
                envsc = copy.deepcopy(self.cc_env)
                flxs = copy.deepcopy(self.fluxes_gj)
                vmm = copy.deepcopy(self.vm)
                vgjj = copy.deepcopy(vgj)
                dvmm = copy.deepcopy(self.dvm)
                aNa = copy.deepcopy(self.active_Na)
                aK = copy.deepcopy(self.active_K)
                self.cc_time.append(concs)
                self.envcc_time.append(envsc)
                self.vm_time.append(vmm)
                self.time.append(t)
                self.fgj_time.append(flxs)
                self.gjopen_time.append(self.gjopen)
                self.vgj_time.append(vgjj)
                self.dvm_time.append(dvmm)
                self.active_Na_time.append(aNa)
                self.active_K_time.append(aK)

        # End off by calculating the current through the gap junction network:
        self.Igj_time = []
        for tflux in self.fgj_time:
            igj=0
            for zi, flx in zip(self.zs,tflux):
                igj = igj+ zi*flx

            igj = -p.F*igj
            self.Igj_time.append(igj)

        if save==True:
            celf = copy.deepcopy(self)
            datadump = [celf,cells,p]
            fh.saveSim(self.savedSim,datadump)
            message_2 = 'Simulation run saved to' + ' ' + p.cache_path
            print(message_2)

        print('Simulation completed successfully.')

    def allDynamics(self,t,p):

        target_length = len(self.scheduled_target_inds)

        if p.scheduled_options['Na_mem'] != 0:

            if p.ions_dict['Na'] == 0 or target_length == 0:
                pass

            else:

                t_on = p.scheduled_options['Na_mem'][0]
                t_off = p.scheduled_options['Na_mem'][1]
                t_change = p.scheduled_options['Na_mem'][2]
                mem_mult_Na = p.scheduled_options['Na_mem'][3]

                effector_Na = pulse(t,t_on,t_off,t_change)

                self.Dm_scheduled[self.iNa][self.scheduled_target_inds] = mem_mult_Na*effector_Na*p.Dm_Na

        if p.scheduled_options['K_mem'] != 0:

            if p.ions_dict['K'] == 0 or target_length == 0:
                pass

            else:
                t_on = p.scheduled_options['K_mem'][0]
                t_off = p.scheduled_options['K_mem'][1]
                t_change = p.scheduled_options['K_mem'][2]
                mem_mult_K = p.scheduled_options['K_mem'][3]

                effector_K = pulse(t,t_on,t_off,t_change)

                self.Dm_scheduled[self.iK][self.scheduled_target_inds] = mem_mult_K*effector_K*p.Dm_K

        if p.scheduled_options['Cl_mem'] != 0:

            if p.ions_dict['Cl'] == 0 or target_length == 0:
                pass

            else:
                t_on = p.scheduled_options['Cl_mem'][0]
                t_off = p.scheduled_options['Cl_mem'][1]
                t_change = p.scheduled_options['Cl_mem'][2]
                mem_mult_Cl = p.scheduled_options['Cl_mem'][3]

                effector_Cl = pulse(t,t_on,t_off,t_change)

                self.Dm_scheduled[self.iCl][self.scheduled_target_inds] = mem_mult_Cl*effector_Cl*p.Dm_Cl

        if p.scheduled_options['Ca_mem'] != 0:

            if p.ions_dict['Ca'] == 0 or target_length == 0:
                pass

            else:

                t_on = p.scheduled_options['Ca_mem'][0]
                t_off = p.scheduled_options['Ca_mem'][1]
                t_change = p.scheduled_options['Ca_mem'][2]
                mem_mult_Ca = p.scheduled_options['Ca_mem'][3]

                effector_Ca = pulse(t,t_on,t_off,t_change)

                self.Dm_scheduled[self.iCa][self.scheduled_target_inds] = mem_mult_Ca*effector_Ca*p.Dm_Ca

        if p.scheduled_options['H_mem'] != 0:

            if p.ions_dict['H'] == 0 or target_length == 0:
                pass

            else:

                t_on = p.scheduled_options['H_mem'][0]
                t_off = p.scheduled_options['H_mem'][1]
                t_change = p.scheduled_options['H_mem'][2]
                mem_mult_H = p.scheduled_options['H_mem'][3]

                effector_H = pulse(t,t_on,t_off,t_change)

                self.Dm_scheduled[self.iH][self.scheduled_target_inds] = mem_mult_H*effector_H*p.Dm_H

        if p.scheduled_options['K_env'] != 0:

            t_on = p.scheduled_options['K_env'][0]
            t_off = p.scheduled_options['K_env'][1]
            t_change = p.scheduled_options['K_env'][2]
            mem_mult_Kenv = p.scheduled_options['K_env'][3]

            effector_Kenv = pulse(t,t_on,t_off,t_change)

            K_env_o = np.ones(len(self.cc_env[self.iK]))
            K_env_o[:] = p.cK_env

            self.cc_env[self.iK] = mem_mult_Kenv*effector_Kenv*K_env_o

        # Voltage gated channel effects

        dvsign = np.sign(self.dvm)

        if p.vg_options['Na_vg'] != 0:

            if p.ions_dict['Na'] == 0 or target_length == 0:
                pass

            else:

                # Logic phase 1: find out which cells have activated their vgNa channels
                truth_vmGTvon_Na = self.vm > self.v_on_Na  # returns bools of vm that are bigger than threshhold
                truth_depol_Na = dvsign==1  # returns bools of vm that are bigger than threshhold
                truth_crossed_inactivate_Na = self.crossed_inactivate_Na == 0  # return bools of vm that can activate
                truth_vgNa_Off = self.vgNa_state == 0 # hasn't been turned on yet

                # find the cell indicies that correspond to all statements of logic phase 1:
                inds_activate_Na = (truth_vmGTvon_Na*truth_depol_Na*truth_crossed_inactivate_Na*truth_vgNa_Off*
                                    self.target_cells).nonzero()

                self.vgNa_state[inds_activate_Na] = 1 # set the crossed_activate term to 1
                self.vgNa_OFFtime[inds_activate_Na] = t + self.t_alive_Na # set the timers for the total active state

                # Logic phase 2: find out which cells have closed their gates due to crossing inactivating voltage or timeout:
                truth_vgNa_ON = self.vgNa_state == 1  # channel must be on already
                truth_vmGTvoff_Na = self.vm > self.v_off_Na  # bools of cells that have vm greater than shut-off volts

                inds_shut_Na = (truth_vgNa_ON*truth_vmGTvoff_Na*self.target_cells).nonzero()

                self.vgNa_state[inds_shut_Na] = 0    # close the vg sodium channels
                self.crossed_inactivate_Na[inds_shut_Na] = 1   # switch these so cells do not re-activate
                self.vgNa_OFFtime[inds_shut_Na] = 0            # reset any timers to zero


                # Logic phase 3: find out if cells have passed below threshhold to become capable of activating:
                truth_vmLTvreact_Na = self.vm < self.v_reactivate_Na # voltage is lower than the deactivate voltage
                truth_vgNa_inactive = self.crossed_inactivate_Na == 1  # the cell had previously been inactivated
                inds_reactivate_Na = (truth_vmLTvreact_Na*truth_vgNa_inactive*self.target_cells).nonzero()

                self.crossed_inactivate_Na[inds_reactivate_Na] = 0  # turn the inhibition to activation off
                self.vgNa_state[inds_reactivate_Na] = 0   # shut the Na channel off if it's on
                self.vgNa_OFFtime[inds_reactivate_Na] = 0            # reset the timers to zero

                # Logic phase 4: find out if cells have timed out:
                truth_vgNa_timeout = self.vgNa_OFFtime < t   # find cells that have timed out their vgNa open state
                inds_timeout_Na = (truth_vgNa_timeout*truth_vgNa_ON*self.target_cells).nonzero()
                self.vgNa_state[inds_timeout_Na] = 0             # set the state to closed
                self.vgNa_OFFtime[inds_timeout_Na] = 0            # reset the timers to zero
                self.crossed_inactivate_Na[inds_timeout_Na] = 1    # inactivate the channel

                # Set activity of Na channel:
                inds_open_Na = (self.vgNa_state == 1).nonzero()
                self.active_Na[inds_open_Na] = 1

                inds_closed_Na =(self.vgNa_state == 0).nonzero()
                self.active_Na[inds_closed_Na] = 0

                self.Dm_vg[self.iNa] = self.maxDmNa*self.active_Na

        if p.vg_options['K_vg'] != 0:

            if p.ions_dict['K'] == 0 or target_length == 0:
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




        # if p.vg_options['K_vg'] !=0:
        #
        #     if p.ions_dict['K'] == 0 or target_length == 0:
        #         pass
        #
        #     else:
        #
        #          # Logic phase 1: find out which cells have activated their vgNa channels
        #         truth_vmGTvon_K = self.vm > self.v_on_K  # returns bools of vm that are bigger than threshhold
        #         truth_depol_K = dvsign==1  # returns bools of vm that are bigger than threshhold
        #         # truth_crossed_inactivate_K = self.crossed_inactivate_K == 0  # return bools of vm that can activate
        #
        #          # set the cell indicies that correspond to all statements of logic phase 1:
        #         #inds_activate_K = (truth_vmGTvon_K*truth_depol_K*truth_crossed_inactivate_K*self.target_cells).nonzero()
        #         inds_activate_K = (truth_vmGTvon_K*truth_depol_K*self.target_cells).nonzero()
        #         self.crossed_activate_K[inds_activate_K] = 1 # set the crossed_activate term to 1
        #
        #          # Logic phase 2: find out which cells have closed their gates due to crossing shut-off voltage:
        #         truth_vmLTvoff_K = self.vm < self.v_off_K  # bools of cells that have vm greater than shut-off volts
        #         inds_shut_K = (truth_vmLTvoff_K*self.target_cells).nonzero()
        #         self.crossed_activate_K[inds_shut_K] = 0    # close the vg sodium channels
        #         # self.crossed_inactivate_K[inds_shut_K] = 1   # switch these so cells do not re-activate
        #
        #         # # Logic phase 3: find out which cells can re-activate
        #         # truth_vmGTvreact_K = self.vm > self.v_reactivate_K
        #         # inds_reactivate_K = (truth_vmGTvreact_K*self.target_cells).nonzero()
        #         # self.crossed_inactivate_K[inds_reactivate_K] = 0
        #
        #         # Set activity of K channel:
        #
        #         inds_open_K = (self.crossed_activate_K == 1).nonzero()
        #         self.active_K[inds_open_K] = 1
        #
        #         inds_closed_K =(self.crossed_activate_K == 0).nonzero()
        #         self.active_K[inds_closed_K] = 0
        #
        #         self.Dm_vg[self.iK] = self.maxDmK*self.active_K


        # finally, add together all effects to make change on the cell membrane permeabilities:
        self.Dm_cells = self.Dm_scheduled + self.Dm_vg + self.Dm_base

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

    assert vola>0
    assert volb>0
    assert d >0

    #volab = (vola + volb)/2
    #qualityfactor = (Dc/d)*(sa/volab)*p.dt   # quality factor should be <1.0 for stable simulations

    flux = -sa*Dc*(cB - cA)/d

    if p.method == 0:

        dmol = sa*p.dt*Dc*(cB - cA)/d

        cA2 = cA + dmol/vola
        cB2 = cB - dmol/volb

        cA2 = check_c(cA2)
        cB2 = check_c(cB2)

    elif p.method == 1:

        k1 = sa*Dc*(cB - cA)/d

        k2 = sa*Dc*(cB - (cA + (1/2)*k1*p.dt))/d

        k3 = sa*Dc*(cB - (cA + (1/2)*k2*p.dt))/d

        k4 = sa*Dc*(cB - (cA + k3*p.dt))/d

        dmol = (p.dt/6)*(k1 + 2*k2 + 2*k3 + k4)

        cA2 = cA + dmol/vola
        cB2 = cB - dmol/volb

        cA2 = check_c(cA2)
        cB2 = check_c(cB2)

    return cA2, cB2, flux

def electrofuse(cA,cB,Dc,d,sa,vola,volb,zc,Vba,p):
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

    alpha = (zc*Vba*p.F)/(p.R*p.T)

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

def pumpNaKATP(cNai,cNao,cKi,cKo,voli,volo,Vm,p):

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

    delG_Na = p.R*p.T*np.log(cNao/cNai) - p.F*Vm
    delG_K = p.R*p.T*np.log(cKi/cKo) + p.F*Vm
    delG_NaKATP = p.deltaGATP - (3*delG_Na + 2*delG_K)
    delG = (delG_NaKATP/1000)

    alpha = p.alpha_NaK*step(delG,p.halfmax_NaK,p.slope_NaK)

    f_Na  = -alpha*cNai*cKo      #flux as [mol/s]
    f_K = -(2/3)*f_Na          # flux as [mol/s]

    if p.method == 0:

        dmol = -alpha*cNai*cKo*p.dt

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

def pumpCaATP(cCai,cCao,voli,volo,Vm,p):

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

    delG_Ca = p.R*p.T*np.log(cCao/cCai) - 2*p.F*Vm
    delG_CaATP = p.deltaGATP - (delG_Ca)
    delG = (delG_CaATP/1000)

    alpha = p.alpha_Ca*step(delG,p.halfmax_Ca,p.slope_Ca)

    f_Ca  = -alpha*cCai      #flux as [mol/s]

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

    q = 0

    for conc,z in zip(concentrations,zs):
        q = q+ conc*z

    netcharge = p.F*q*vol

    return netcharge

def check_c(cA):
    """
    Does a quick check on two values (concentrations)
    and sets one to zero if it is below zero.

    """
    # if isinstance(cA,np.float64):  # if we just have a singular value
    #     if cA < 0.0:
    #         cA2 = 0.0
    #
    # elif isinstance(cA,np.ndarray): # if we have matrix data
    isubzeros = (cA<0).nonzero()
    if isubzeros:  # if there's anything in the isubzeros matrix...
        cA[isubzeros] = 0.0

    return cA

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
    assert x.all() > 0

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



