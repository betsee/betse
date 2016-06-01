#!/usr/bin/env python3
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

# FIXME include other channels in morphogen (dye) dynamics

from random import shuffle

import numpy as np
from scipy import interpolate as interp
from scipy import spatial as sps

from betse.science import toolbox as tb
from betse.science.event import modulators as mod
from betse.science.tissue.channels import vg_na as vgna
from betse.science.tissue.channels_o import vgPotassium, cagPotassium, vgCalcium, vgPotassium_init, vgCalcium_init
from betse.util.io.log import logs
from betse.util.type import types


class TissueHandler(object):
    '''
    A high-level handler for user-specified tissue-centric functionality,
    including both tissue profiles _and_ scheduled interventions.

    This handler governs all:

    * Tissue profiles and objects required by these profiles, including all
      geometry-specifying bitmaps.
    * Scheduled interventions, even those _not_ pertaining to tissue profiles
      (e.g., global scheduled interventions).

    Attributes (General)
    ----------------------------
    '''

    def __init__(
        self, sim: 'Simulator', cells: 'Cells', p: 'Parameters') -> None:
        #FIXME: Reduce to the following single line:
        #    self.data_length = len(cells.mem_i if p.sim_ECM else cells.cell_i)
        #Actually, as this duplicates logic below, we should probably just
        #move this logic into a new private method _init_data_length().
        #Actually, the ideal solution is to shift this into the
        #tissueProfiles() method if feasible; doing so would permit us to
        #simply remove this logic both here and below.
        if p.sim_ECM is True:
            self.data_length = len(cells.mem_i)
        else:
            self.data_length = len(cells.mem_i)

    def runAllInit(
        self, sim: 'Simulator', cells: 'Cells', p: 'Parameters') -> None:
        '''
        Initialize all tissue manipulations specified by the passed
        user-specified parameters with the passed tissue simulation and
        cellular world.
        '''

        self._init_events_global(  sim, cells, p)
        self._init_events_tissue(  sim, cells, p)
        self._init_channels_tissue(sim, cells, p)

    def runAllDynamics(self, sim, cells, p, t):
        '''
        Apply all tissue manipulations specified by the passed user-specified
        parameters to the passed tissue simulation and cellular world for the
        passed time step.
        '''

        self._sim_events_global(  sim, cells, p, t)
        self._sim_events_tissue(  sim, cells, p, t)
        self._sim_channels_tissue(sim, cells, p, t)
        self.makeAllChanges(sim)

    def _init_events_global(self,sim,cells,p):
        '''
        Initialize all **global scheduled interventions** (i.e., events globally
        applicable to all cells) specified by the passed user-specified
        parameters with the passed tissue simulation and cellular world.
        '''

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

        if p.global_options['Morph_env'] != 0:

            self.t_on_MorphEnv = p.global_options['Morph_env'][0]
            self.t_off_MorphEnv = p.global_options['Morph_env'][1]
            self.t_change_MorphEnv = p.global_options['Morph_env'][2]
            self.conc_MorphEnv = p.global_options['Morph_env'][3]

        if p.global_options['T_change'] != 0:

            self.tonT = p.global_options['T_change'][0]
            self.toffT = p.global_options['T_change'][1]
            self.trampT = p.global_options['T_change'][2]
            self.multT = p.global_options['T_change'][3]

        if p.global_options['gj_block'] != 0:

            self.tonGJ = p.global_options['gj_block'][0]
            self.toffGJ = p.global_options['gj_block'][1]
            self.trampGJ = p.global_options['gj_block'][2]

            numo = p.global_options['gj_block'][3]

            numo = int(numo)

            if numo > 100:
                numo = 100
            elif numo < 1:
                numo = 1

            data_length = len(cells.mem_i)
            data_fraction = int((numo/100)*data_length)

            mem_i_copy = np.copy(cells.mem_i[:])
            shuffle(mem_i_copy)
            self.targets_gj_block = [cells.mem_i[x] for x in range(0,data_fraction)]


        if p.global_options['NaKATP_block'] != 0:
            self.tonNK = p.global_options['NaKATP_block'][0]
            self.toffNK = p.global_options['NaKATP_block'][1]
            self.trampNK = p.global_options['NaKATP_block'][2]

        if p.global_options['HKATP_block'] != 0:
            self.tonHK = p.global_options['HKATP_block'][0]
            self.toffHK = p.global_options['HKATP_block'][1]
            self.trampHK = p.global_options['HKATP_block'][2]

        if p.global_options['VATP_block'] != 0:
            self.tonV = p.global_options['VATP_block'][0]
            self.toffV = p.global_options['VATP_block'][1]
            self.trampV = p.global_options['VATP_block'][2]

    def _init_events_tissue(self, sim, cells, p):
        '''
        Initialize all **targeted scheduled interventions** (i.e., events only
        applicable to specific tissue profiles) specified by the passed
        user-specified parameters with the passed tissue simulation and cellular
        world.
        '''

        if p.scheduled_options['Na_mem'] != 0:
            self.t_on_Namem = p.scheduled_options['Na_mem'][0]
            self.t_off_Namem = p.scheduled_options['Na_mem'][1]
            self.t_change_Namem = p.scheduled_options['Na_mem'][2]
            self.mem_mult_Namem = p.scheduled_options['Na_mem'][3]
            self.apply_Namem = p.scheduled_options['Na_mem'][4]
            self.function_Namem = p.scheduled_options['Na_mem'][5]

            self.targets_Namem = []
            for profile in self.apply_Namem:
                targets = self.tissue_target_inds[profile]
                self.targets_Namem.append(targets)

            self.targets_Namem = [
                item for sublist in self.targets_Namem for item in sublist]

            self.scalar_Namem = 1
            self.dyna_Namem = lambda t: 1

            # call a special toolbox function to change membrane permeability: spatial grads
            # 'gradient_x', 'gradient_y', 'gradient_r'

            if self.function_Namem != 'None':
                self.scalar_Namem, self.dyna_Namem = getattr(mod, self.function_Namem)(self.targets_Namem,cells,p)

        #----------------------------------------------

        if p.scheduled_options['K_mem'] != 0:
            self.t_on_Kmem = p.scheduled_options['K_mem'][0]
            self.t_off_Kmem = p.scheduled_options['K_mem'][1]
            self.t_change_Kmem = p.scheduled_options['K_mem'][2]
            self.mem_mult_Kmem = p.scheduled_options['K_mem'][3]
            self.apply_Kmem = p.scheduled_options['K_mem'][4]
            self.function_Kmem = p.scheduled_options['K_mem'][5]

            self.targets_Kmem = []
            for profile in self.apply_Kmem:
                targets = self.tissue_target_inds[profile]
                self.targets_Kmem.append(targets)

            self.targets_Kmem = [
                item for sublist in self.targets_Kmem for item in sublist]

            self.scalar_Kmem = 1
            self.dyna_Kmem = lambda t: 1

            if self.function_Kmem != 'None':
                # call a special toolbox function to change membrane permeability: spatial grads
                # 'gradient_x', 'gradient_y', 'gradient_r'

                self.scalar_Kmem, self.dyna_Kmem = getattr(mod,self.function_Kmem)(self.targets_Kmem,cells,p)

        #----------------------------------------------

        if p.scheduled_options['Cl_mem'] != 0:
            self.t_on_Clmem = p.scheduled_options['Cl_mem'][0]
            self.t_off_Clmem = p.scheduled_options['Cl_mem'][1]
            self.t_change_Clmem = p.scheduled_options['Cl_mem'][2]
            self.mem_mult_Clmem = p.scheduled_options['Cl_mem'][3]
            self.apply_Clmem = p.scheduled_options['Cl_mem'][4]
            self.function_Clmem = p.scheduled_options['Cl_mem'][5]

            self.targets_Clmem = []
            for profile in self.apply_Clmem:
                targets = self.tissue_target_inds[profile]
                self.targets_Clmem.append(targets)

            self.targets_Clmem = [
                item for sublist in self.targets_Clmem for item in sublist]

            self.scalar_Clmem = 1
            self.dyna_Clmem = lambda t: 1

            if self.function_Clmem != 'None':

                # call a special toolbox function to change membrane permeability: spatial grads
                # 'gradient_x', 'gradient_y', 'gradient_r'
                self.scalar_Clmem, self.dyna_Clmem = getattr(mod,self.function_Clmem)(self.targets_Clmem,cells,p)

        #----------------------------------------------

        if p.scheduled_options['Ca_mem'] != 0:
            self.t_on_Camem = p.scheduled_options['Ca_mem'][0]
            self.t_off_Camem = p.scheduled_options['Ca_mem'][1]
            self.t_change_Camem = p.scheduled_options['Ca_mem'][2]
            self.mem_mult_Camem = p.scheduled_options['Ca_mem'][3]
            self.apply_Camem = p.scheduled_options['Ca_mem'][4]
            self.function_Camem = p.scheduled_options['Ca_mem'][5]

            self.targets_Camem = []
            for profile in self.apply_Camem:
                targets = self.tissue_target_inds[profile]
                self.targets_Camem.append(targets)

            self.targets_Camem = [
                item for sublist in self.targets_Camem for item in sublist]

            self.scalar_Camem = 1
            self.dyna_Camem = lambda t: 1

            if self.function_Camem != 'None':
                # call a special toolbox function to change membrane permeability: spatial grads
                # 'gradient_x', 'gradient_y', 'gradient_r'

                self.scalar_Camem, self.dyna_Camem = getattr(mod, self.function_Camem)(self.targets_Camem,cells,p)

        #----------------------------------------------

        if p.scheduled_options['IP3'] != 0:
            self.t_onIP3 = p.scheduled_options['IP3'][0]
            self.t_offIP3 = p.scheduled_options['IP3'][1]
            self.t_changeIP3 = p.scheduled_options['IP3'][2]
            self.rate_IP3 = p.scheduled_options['IP3'][3]
            self.apply_IP3 = p.scheduled_options['IP3'][4]
            self.function_IP3 = p.scheduled_options['IP3'][5]

            self.targets_IP3 = []
            for profile in self.apply_IP3:
                targets = self.cell_target_inds[profile]
                self.targets_IP3.append(targets)

            self.targets_IP3 = [
                item for sublist in self.targets_IP3 for item in sublist]

            self.scalar_IP3 = 1
            self.dyna_IP3 = lambda t: 1

            # call a special toolbox function to change membrane permeability: spatial grads
            # 'gradient_x', 'gradient_y', 'gradient_r'
            if self.function_IP3 != 'None':
                self.scalar_IP3, self.dyna_IP3 = getattr(mod, self.function_IP3)(self.targets_IP3,cells, p)

        #----------------------------------------------

        if p.scheduled_options['pressure'] != 0:

            self.t_onP = p.scheduled_options['pressure'][0]
            self.t_offP = p.scheduled_options['pressure'][1]
            self.t_changeP = p.scheduled_options['pressure'][2]
            self.rate_P = p.scheduled_options['pressure'][3]
            self.apply_P = p.scheduled_options['pressure'][4]
            self.function_P = p.scheduled_options['pressure'][5]

            self.targets_P = []
            for profile in self.apply_P:
                targets = self.cell_target_inds[profile]
                self.targets_P.append(targets)

            self.targets_P = [
                item for sublist in self.targets_P for item in sublist]

            self.scalar_P = 1
            self.dyna_P = lambda t: 1

            # call a special toolbox function to change membrane permeability: spatial grads
            # 'gradient_x', 'gradient_y', 'gradient_r'
            if self.function_P != 'None':
                self.scalar_P, self.dyna_P = getattr(mod, self.function_P)(self.targets_P,cells, p)

        #--------------------------------------------------------

        if p.scheduled_options['ecmJ'] != 0 and p.sim_ECM is True:
            self.t_on_ecmJ  = p.scheduled_options['ecmJ'][0]
            self.t_off_ecmJ = p.scheduled_options['ecmJ'][1]
            self.t_change_ecmJ = p.scheduled_options['ecmJ'][2]
            self.apply_ecmJ = p.scheduled_options['ecmJ'][3]

            self.targets_ecmJ = []
            for profile in self.apply_ecmJ:
                targets = self.env_target_inds[profile]
                self.targets_ecmJ.append(targets)

            self.targets_ecmJ = [
                item for sublist in self.targets_ecmJ for item in sublist]

    def _init_channels_tissue(self, sim, cells, p):
        '''
        Initialize all **targeted ion channels** (i.e., ion channels only
        applicable to specific tissue profiles) specified by the passed
        user-specified parameters with the passed tissue simulation and
        cellular world.
        '''

        if p.vgNa_bool:

            # Initialization of logic values for voltage gated sodium channel
            self.maxDmNa = p.vgNa_max
            # self.v_activate_Na = p.vg_options['Na_vg'][1]
            # self.v_inactivate_Na = p.vg_options['Na_vg'][2]
            # self.v_deactivate_Na = p.vg_options['Na_vg'][3]
            # self.t_alive_Na = p.vg_options['Na_vg'][4]
            # self.t_dead_Na = p.vg_options['Na_vg'][5]
            self.apply_vgNa = p.vgNa_apply

            # Initialize matrices defining states of vgNa channels for each cell membrane:
            # self.inactivated_Na = np.zeros(self.data_length)
            # self.vgNa_state = np.zeros(self.data_length)
            #
            # self.vgNa_aliveTimer = np.zeros(self.data_length) # sim time at which vgNa starts to close if activated
            # self.vgNa_deadTimer = np.zeros(self.data_length) # sim time at which vgNa reactivates after inactivation

            self.targets_vgNa = []
            for profile in self.apply_vgNa:
                targets = self.tissue_target_inds[profile]
                self.targets_vgNa.append(targets)

            self.targets_vgNa = [item for sublist in self.targets_vgNa for item in sublist]
            self.targets_vgNa = np.asarray(self.targets_vgNa)

            # self.target_mask_vgNa = np.zeros(self.data_length)
            # self.target_mask_vgNa[self.targets_vgNa] = 1

            # self.m_Na = np.zeros(len(self.targets_vgNa))
            # self.h_Na = np.zeros(len(self.targets_vgNa))

            # create the desired voltage gated sodium channel instance:

            Na_class_ = getattr(vgna,p.vgNa_type,'vgNa_Default')
            self.vgNa_object = Na_class_()

            if p.run_sim is True:
                # initialize the voltage-gated sodium object
                self.vgNa_object.init(self, sim, p)
                # vgSodium_init(self,sim,p)


            # self.m_Na = self.mInf_o
            # self.h_Na = self.hInf_o


        if p.vg_options['K_vg'] !=0:
            # Initialization of logic values forr voltage gated potassium channel
            self.maxDmK = p.vg_options['K_vg'][0]
            self.v_on_K = p.vg_options['K_vg'][1]
            self.v_off_K = p.vg_options['K_vg'][2]
            self.t_alive_K = p.vg_options['K_vg'][3]
            self.apply_vgK = p.vg_options['K_vg'][4]

            # Initialize matrices defining states of vgK channels for each cell:
            self.active_K = np.zeros(self.data_length)
            self.crossed_activate_K = np.zeros(self.data_length)
            self.crossed_inactivate_K = np.zeros(self.data_length)

            # Initialize other matrices for vgK timing logic: NEW!
            self.vgK_state = np.zeros(self.data_length)   # state can be 0 = off, 1 = open
            self.vgK_OFFtime = np.zeros(self.data_length) # sim time at which vgK starts to close

            self.targets_vgK = []
            for profile in self.apply_vgK:
                targets = self.tissue_target_inds[profile]
                self.targets_vgK.append(targets)

            self.targets_vgK = [item for sublist in self.targets_vgK for item in sublist]

            self.targets_vgK = np.asarray(self.targets_vgK)

            self.target_mask_vgK = np.zeros(self.data_length)
            self.target_mask_vgK[self.targets_vgK] = 1

            # self.m_K = np.zeros(len(self.targets_vgK))
            # self.h_K = np.zeros(len(self.targets_vgK))
            # self.n_K = np.zeros(len(self.targets_vgK))

            if p.run_sim is True:

                vgPotassium_init(self,sim,p)


        if p.vg_options['Ca_vg'] !=0:
            # Initialization of logic values for voltage gated calcium channel
            self.maxDmCa = p.vg_options['Ca_vg'][0]
            self.v_on_Ca = p.vg_options['Ca_vg'][1]
            self.v_off_Ca = p.vg_options['Ca_vg'][2]
            self.ca_upper_ca = p.vg_options['Ca_vg'][3]
            self.ca_lower_ca = p.vg_options['Ca_vg'][4]
            self.apply_vgCa = p.vg_options['Ca_vg'][5]

            # Initialize matrices defining states of vgK channels for each cell membrane:
            self.active_Ca = np.zeros(self.data_length)

            self.vgCa_state = np.zeros(self.data_length)   # state can be 0 = off, 1 = open

            self.targets_vgCa = []
            for profile in self.apply_vgCa:
                targets = self.tissue_target_inds[profile]
                self.targets_vgCa.append(targets)

            self.targets_vgCa = [item for sublist in self.targets_vgCa for item in sublist]

            self.targets_vgCa = np.asarray(self.targets_vgCa)

            self.target_mask_vgCa = np.zeros(self.data_length)
            self.target_mask_vgCa[self.targets_vgCa] = 1

            # self.m_Ca = np.zeros(self.data_length)
            # self.h_Ca = np.zeros(self.data_length)
            #
            # self.m_Ca[:] = 1.0000/(1+ np.exp(-((sim.vm*1e3 - 10) + 30.000)/6))
            # self.h_Ca[:] = 1.0000/(1+ np.exp(((sim.vm*1e3 - 10) + 80.000)/6.4))

            if p.run_sim is True:

                vgCalcium_init(self,sim,p)

        if p.vg_options['K_cag'] != 0:
            self.maxDmKcag = p.vg_options['K_cag'][0]
            self.Kcag_halfmax = p.vg_options['K_cag'][1]
            self.Kcag_n = p.vg_options['K_cag'][2]
            self.apply_cagK = p.vg_options['K_cag'][3]

            # Initialize matrices defining states of cag K channels for each cell membrane:
            self.active_cagK = np.zeros(self.data_length)

            self.targets_cagK = []
            for profile in self.apply_cagK:
                targets = self.tissue_target_inds[profile]
                self.targets_cagK.append(targets)

            self.targets_cagK = [item for sublist in self.targets_cagK for item in sublist]

            self.targets_cagK = np.asarray(self.targets_cagK)

            self.target_mask_cagK = np.zeros(self.data_length)
            self.target_mask_cagK[self.targets_cagK] = 1

        # stretch activated ion channels:
        if p.vg_options['Na_stretch'] != 0:
            self.maxDmNaStretch = p.vg_options['Na_stretch'][0]
            self.NaStretch_halfmax = p.vg_options['Na_stretch'][1]
            self.NaStretch_n = p.vg_options['Na_stretch'][2]
            self.apply_NaStretch = p.vg_options['Na_stretch'][3]

             # Initialize matrices defining states of cag K channels for each cell membrane:
            self.active_NaStretch = np.zeros(self.data_length)

            self.targets_NaStretch = []
            for profile in self.apply_NaStretch:
                targets = self.tissue_target_inds[profile]
                self.targets_NaStretch.append(targets)

            self.targets_NaStretch = [item for sublist in self.targets_NaStretch for item in sublist]

            self.targets_NaStretch = np.asarray(self.targets_NaStretch)

            self.target_mask_NaStretch = np.zeros(self.data_length)
            self.target_mask_NaStretch[self.targets_NaStretch] = 1

        # calcium dynamics
        if p.Ca_dyn_options['CICR'] != 0:
            self.stateER = np.zeros(len(cells.cell_i))   # state of ER membrane Ca permeability

            self.maxDmCaER = p.Ca_dyn_options['CICR'][0][0]
            self.topCa = p.Ca_dyn_options['CICR'][0][1]
            self.bottomCa =  p.Ca_dyn_options['CICR'][0][2]

            if len(p.Ca_dyn_options['CICR'][1])!=0:
                self.midCaR = p.Ca_dyn_options['CICR'][1][0]
                self.widthCaR = p.Ca_dyn_options['CICR'][1][1]

            if len(p.Ca_dyn_options['CICR'][2])!=0:
                self.KhmIP3 = p.Ca_dyn_options['CICR'][2][0]
                self.n_IP3 = p.Ca_dyn_options['CICR'][2][1]

            self.apply_Ca = p.Ca_dyn_options['CICR'][3]

            #FIXME: The following five lines appear to be reducible to just:
            #
            # self.targets_Ca = list(itertools.chain.from_iterable(
            #     self.cell_target_inds[profile]
            #     for profile in self.apply_Ca))
            #
            #If that's a bit too much Python magic, you can also just replace
            #append() by extend() and remove the final list comprehension:
            #
            # self.targets_Ca = []
            # for profile in self.apply_Ca:
            #     self.targets_Ca.extend(self.cell_target_inds[profile])
            #
            #This probably applies everywhere above as well. Shimmy up, Jimbo!
            self.targets_Ca = []
            for profile in self.apply_Ca:
                targets = self.cell_target_inds[profile]
                self.targets_Ca.append(targets)
            self.targets_Ca = [item for sublist in self.targets_Ca for item in sublist]

            self.targets_Ca = np.asarray(self.targets_Ca)

            self.target_mask_Ca = np.zeros(len(cells.cell_i))
            self.target_mask_Ca[self.targets_Ca] = 1

    def _sim_events_global(self, sim, cells, p, t):
        '''
        Apply all **global scheduled interventions** (i.e., events globally
        applicable to all cells) specified by the passed user-specified
        parameters to the passed tissue simulation and cellular world for the
        passed time step.
        '''

        if p.global_options['K_env'] != 0:
            effector_Kenv = tb.pulse(t,self.t_on_Kenv,self.t_off_Kenv,self.t_change_Kenv)

            if p.sim_ECM is False:
                sim.cc_env[sim.iK][:] = self.mem_mult_Kenv*effector_Kenv*p.cK_env + p.cK_env

            elif p.sim_ECM is True: # simulate addition of potassium salt to remain charge neutral
                sim.c_env_bound[sim.iK] = self.mem_mult_Kenv*effector_Kenv*p.env_concs['K'] + p.env_concs['K']
                sim.c_env_bound[sim.iM] = self.mem_mult_Kenv*effector_Kenv*p.env_concs['K'] + p.env_concs['M']

        if p.global_options['Cl_env'] != 0 and p.ions_dict['Cl'] == 1:
            effector_Clenv = tb.pulse(t,self.t_on_Clenv,self.t_off_Clenv,self.t_change_Clenv)

            if p.sim_ECM is False:
                sim.cc_env[sim.iCl][:] = self.mem_mult_Clenv*effector_Clenv*p.cCl_env + p.cCl_env

            elif p.sim_ECM is True:  # simulate addition of sodium chloride to remain charge neutral
                sim.c_env_bound[sim.iCl] = self.mem_mult_Clenv*effector_Clenv*p.env_concs['Cl'] + p.env_concs['Cl']
                sim.c_env_bound[sim.iNa] = self.mem_mult_Clenv*effector_Clenv*p.env_concs['Cl'] + p.env_concs['Na']

        if p.global_options['Na_env'] != 0:
            effector_Naenv = tb.pulse(t,self.t_on_Naenv,self.t_off_Naenv,self.t_change_Naenv)

            if p.sim_ECM is False:
                sim.cc_env[sim.iNa][:] = self.mem_mult_Naenv*effector_Naenv*p.cNa_env + p.cNa_env

            elif p.sim_ECM is True: # simulate addition of sodium salt to remain charge neutral
                sim.c_env_bound[sim.iNa] = self.mem_mult_Naenv*effector_Naenv*p.env_concs['Na'] + p.env_concs['Na']
                sim.c_env_bound[sim.iM] = self.mem_mult_Naenv*effector_Naenv*p.env_concs['Na'] + p.env_concs['M']

        if p.global_options['Morph_env'] != 0 and p.voltage_dye is True:
            effector_MorphEnv = tb.pulse(t,self.t_on_MorphEnv,self.t_off_MorphEnv,self.t_change_MorphEnv)

            if p.sim_ECM is False:
                sim.cDye_env[:] = self.conc_MorphEnv*effector_MorphEnv + sim.cDye_env*(1-effector_MorphEnv)

            elif p.sim_ECM is True: # simulate addition of counter salt to maintain charge neutrality:
                sim.c_dye_bound = self.conc_MorphEnv*effector_MorphEnv + p.cDye_to*(1-effector_MorphEnv)

        if p.global_options['T_change'] != 0:
            sim.T = self.multT*tb.pulse(t,self.tonT,self.toffT,self.trampT)*p.T + p.T

        if p.global_options['gj_block'] != 0:
            sim.gj_block[self.targets_gj_block] = (1.0 - tb.pulse(t,self.tonGJ,self.toffGJ,self.trampGJ))
            # sim.gj_block[:] = (1.0 - tb.pulse(t,self.tonGJ,self.toffGJ,self.trampGJ))

        if p.global_options['NaKATP_block'] != 0:
            sim.NaKATP_block = (1.0 - tb.pulse(t,self.tonNK,self.toffNK,self.trampNK))

        if p.global_options['HKATP_block'] != 0:
            sim.HKATP_block = (1.0 - tb.pulse(t,self.tonHK,self.toffHK,self.trampHK))

        if p.global_options['VATP_block'] != 0:
            sim.VATP_block = (1.0 - tb.pulse(t,self.tonV,self.toffV,self.trampV))

    def _sim_events_tissue(self, sim, cells, p, t):
        '''
        Apply all **targeted scheduled interventions** (i.e., events only
        applicable to specific tissue profiles) specified by the passed
        user-specified parameters to the passed tissue simulation and cellular
        world for the passed time step.
        '''

        if p.scheduled_options['Na_mem'] != 0:
            effector_Na = self.scalar_Namem*self.dyna_Namem(t)*\
                          tb.pulse(t,self.t_on_Namem,self.t_off_Namem,self.t_change_Namem)



            sim.Dm_scheduled[sim.iNa][self.targets_Namem] = self.mem_mult_Namem*effector_Na*p.Dm_Na



        if p.scheduled_options['K_mem'] != 0:
            effector_K = self.scalar_Kmem*self.dyna_Kmem(t)*\
                         tb.pulse(t,self.t_on_Kmem,self.t_off_Kmem,self.t_change_Kmem)

            sim.Dm_scheduled[sim.iK][self.targets_Kmem] = self.mem_mult_Kmem*effector_K*p.Dm_K

        if p.scheduled_options['Cl_mem'] != 0 and p.ions_dict['Cl'] != 0:
            effector_Cl = self.scalar_Clmem*self.dyna_Clmem(t)*\
                          tb.pulse(t,self.t_on_Clmem,self.t_off_Clmem,self.t_change_Clmem)

            sim.Dm_scheduled[sim.iCl][self.targets_Clmem] = self.mem_mult_Clmem*effector_Cl*p.Dm_Cl

        if p.scheduled_options['Ca_mem'] != 0 and p.ions_dict['Ca'] != 0:
            effector_Ca = self.scalar_Camem*self.dyna_Camem(t)*\
                          tb.pulse(t,self.t_on_Camem,self.t_off_Camem,self.t_change_Camem)

            sim.Dm_scheduled[sim.iCa][self.targets_Camem] = self.mem_mult_Camem*effector_Ca*p.Dm_Ca

        if p.scheduled_options['IP3'] != 0:
            sim.cIP3[self.targets_IP3] = sim.cIP3[self.targets_IP3] + \
                                         self.scalar_IP3*self.dyna_IP3(t)*self.rate_IP3*tb.pulse(t,self.t_onIP3,
                                         self.t_offIP3,self.t_changeIP3)

        if p.scheduled_options['pressure'] != 0:

            sim.P_mod[self.targets_P] =  self.scalar_P*self.dyna_P(t)*self.rate_P*tb.pulse(t,self.t_onP,
                                         self.t_offP,self.t_changeP)

        if p.sim_ECM is True:

            if p.scheduled_options['ecmJ'] != 0:
                for i, dmat in enumerate(sim.D_env):
                    effector_ecmJ = tb.pulse(t,self.t_on_ecmJ,self.t_off_ecmJ,self.t_change_ecmJ)
                    sim.D_env[i][self.targets_ecmJ] = sim.D_env_base[i][self.targets_ecmJ]*(1 - effector_ecmJ) \
                                                      + effector_ecmJ*sim.D_free[i]

                    sim.D_env_weight = sim.D_env_weight.ravel()
                    sim.D_env_weight_base = sim.D_env_weight_base.ravel()

                    sim.D_env_weight[self.targets_ecmJ] = \
                        sim.D_env_weight_base[self.targets_ecmJ]*(1-effector_ecmJ) + \
                                                          effector_ecmJ

                    sim.D_env_weight = sim.D_env_weight.reshape(cells.X.shape)
                    sim.D_env_weight_base = sim.D_env_weight_base.reshape(cells.X.shape)

                for i, dmat in enumerate(sim.D_env):
                    if p.env_type is True:
                        sim.D_env_u[i] = interp.griddata((cells.xypts[:,0],cells.xypts[:,1]),dmat.ravel(),
                            (cells.grid_obj.u_X,cells.grid_obj.u_Y),method='nearest',fill_value = sim.D_free[i])
                        sim.D_env_v[i] = interp.griddata((cells.xypts[:,0],cells.xypts[:,1]),dmat.ravel(),
                            (cells.grid_obj.v_X,cells.grid_obj.v_Y),method='nearest',fill_value=sim.D_free[i])

                    else:
                        sim.D_env_u[i] = interp.griddata((cells.xypts[:,0],cells.xypts[:,1]),dmat.ravel(),
                            (cells.grid_obj.u_X,cells.grid_obj.u_Y),method='nearest',fill_value = 0)

                        sim.D_env_v[i] = interp.griddata((cells.xypts[:,0],cells.xypts[:,1]),dmat.ravel(),
                            (cells.grid_obj.v_X,cells.grid_obj.v_Y),method='nearest',fill_value = 0)

                sim.D_env_weight_u = sim.D_env_u[sim.iP]/sim.D_env_u[sim.iP].max()
                sim.D_env_weight_v = sim.D_env_v[sim.iP]/sim.D_env_v[sim.iP].max()

                if p.closed_bound is True:  # set full no slip boundary condition at exterior bounds
                    sim.D_env_weight_u[:,0] = 0
                    sim.D_env_weight_u[:,-1] = 0
                    sim.D_env_weight_u[0,:] = 0
                    sim.D_env_weight_u[-1,:] = 0

                    sim.D_env_weight_v[:,0] = 0
                    sim.D_env_weight_v[:,-1] = 0
                    sim.D_env_weight_v[0,:] = 0
                    sim.D_env_weight_v[-1,:] = 0

        soc = p.scheduled_options['cuts']
        if soc is not None and not soc._is_fired and t >= soc.time:
            for cut_profile_name in soc.profile_names:
                logs.log_info(
                    'Cutting cell cluster via cut profile "%s"...',
                    cut_profile_name)
                removeCells(p.profiles[cut_profile_name].picker, sim, cells, p)

            logs.log_info("Cutting event successful! Resuming simulation...")

            #FIXME: Duplicate logic. See above. The snow bear dances at noon.

            # redo main data length variable for this dynamics module with updated world:
            if p.sim_ECM is True:
                self.data_length = len(cells.mem_i)
            else:
                self.data_length = len(cells.mem_i)

            self.tissueProfiles(sim, cells, p)
            cells.redo_gj(self, p, savecells=False)
            self.runAllInit(sim, cells, p)

            #FIXME: Excise. And the chill vermouth of winter did exhale!
            # Clear and recreate the currently displayed and/or saved
            # animation. Cutting requires sufficiently "heavy" modifications to
            # plot data that starting over from scratch is the safest and
            # simplest approach.
            # sim._dereplot_loop(p)

            # Avoid repeating this cutting event at subsequent time steps.
            soc._is_fired = True

        # If the voltage event is enabled, adjust the voltage accordingly.
        if p.scheduled_options['extV'] is not None:
           p.scheduled_options['extV'].fire(sim, t)

    def _sim_channels_tissue(self, sim, cells, p, t):
        '''
        Handle all **targeted ion channels** (i.e., ion channels only applicable
        to specific tissue profiles) specified by the passed user-specified
        parameters to the passed tissue simulation and cellular world for the
        passed time step.
        '''

        self.dvsign = np.sign(sim.dvm)

        if p.vgNa_bool:

            # update the voltage-gated sodium object
            self.vgNa_object.run(self, sim, p)

            # vgSodium(self,sim,cells,p)

        if p.vg_options['K_vg'] !=0:

            vgPotassium(self,sim,cells,p)

        if p.vg_options['Ca_vg'] !=0 and p.ions_dict['Ca'] != 0:

            vgCalcium(self,sim,cells,p)

        if p.vg_options['K_cag'] != 0 and p.ions_dict['Ca'] != 0:

            cagPotassium(self,sim,cells,p)

        if p.vg_options['Na_stretch'] != 0 and p.deformation is True:

            self.stretchChannel(sim,cells,p,t)

        if p.Ca_dyn_options['CICR'] != 0 and p.ions_dict['Ca'] != 0:

            self.calciumDynamics(sim,cells,p)

        # update membrane permeability if dye targets an ion channel:
        if p.voltage_dye is True and sim.dye_target is not None:

            if p.Dye_acts_extracell is False:

                sim.Dm_mod_dye = p.Dye_peak_channel*tb.hill(sim.cDye_cell,p.Dye_Hill_K,p.Dye_Hill_exp)

                if p.sim_ECM is True:
                    sim.Dm_morpho[sim.dye_target] = sim.Dm_mod_dye[cells.mem_to_cells]

                else:
                    sim.Dm_morpho[sim.dye_target] = sim.Dm_mod_dye

            elif p.Dye_acts_extracell is True and p.sim_ECM is True:

                sim.Dm_mod_dye = p.Dye_peak_channel*tb.hill(sim.cDye_env,p.Dye_Hill_K,p.Dye_Hill_exp)

                sim.Dm_morpho[sim.dye_target] = sim.Dm_mod_dye[cells.map_mem2ecm]

    def stretchChannel(self,sim,cells,p,t):

        dd = np.sqrt(sim.d_cells_x**2 + sim.d_cells_y**2)

        eta = (dd/cells.R)

        self.active_NaStretch[self.targets_NaStretch] = tb.hill(eta[self.targets_NaStretch],
                self.NaStretch_halfmax,self.NaStretch_n)

        sim.Dm_stretch[sim.iNa] = self.maxDmNaStretch*self.active_NaStretch

    def calciumDynamics(self,sim,cells,p):

        if p.Ca_dyn_options['CICR'] != 0:

            if len(p.Ca_dyn_options['CICR'][1])==0:
                term_Ca_reg = 1.0

            else:
                term_Ca_reg = (np.exp(-((sim.cc_mems[sim.iCa]-self.midCaR)**2)/((2*self.widthCaR)**2)))

            if len(p.Ca_dyn_options['CICR'][2]) == 0:
                term_IP3_reg = 1.0

            else:
                term_IP3_reg = tb.hill(sim.cIP3,self.KhmIP3,self.n_IP3)

            sim.Dm_er_CICR[0] = self.maxDmCaER*term_IP3_reg*term_Ca_reg
            sim.Dm_er = sim.Dm_er_CICR + sim.Dm_er_base

    def tissueProfiles(self, sim, cells, p):
        '''
        Create cell-specific (and if simulating extracellular spaces, membrane-
        specific as well) index sets for all user-defined tissue profiles.
        '''

        profile_names = list(p.profiles.keys())
        self.tissue_target_inds = {}
        self.cell_target_inds = {}
        self.env_target_inds = {}
        self.tissue_profile_names = []

        # Go through again and do traditional tissue profiles:
        for profile_name in profile_names:
            profile = p.profiles[profile_name]

            #FIXME: A pretty horrible hack. It will go away. Praise the flower!
            profile_type = (
                profile['type'] if isinstance(profile, dict) else 'cut')

            # If this is a tissue profile...
            if profile_type == 'tissue':
                profile_picker = profile['picker']
                dmem_list = profile['diffusion constants']

                self.tissue_profile_names.append(profile_name)

                #FIXME: Somewhat ambiguous attribute names. Ideally:
                #
                #* "tissue_target_inds" should probably be renamed to something
                #  like "cell_inds_sans_ecm".
                #* "cell_target_inds" should probably be renamed to something
                #  like "cell_inds_with_ecm".
                #
                #Burning oil of hubris by the dawn's yearning pangs of joy!

                self.tissue_target_inds[profile_name] = (
                    profile_picker.get_cell_indices(cells, p, ignoreECM=False))
                self.cell_target_inds[profile_name] = (
                    profile_picker.get_cell_indices(cells, p, ignoreECM=True))

                if len(self.cell_target_inds[profile_name]):
                    # Get ECM targets.
                    if p.sim_ECM is True:
                        ecm_targs_cell = list(
                            cells.map_cell2ecm[self.cell_target_inds[profile_name]])
                        ecm_targs_mem = list(
                            cells.map_mem2ecm[self.tissue_target_inds[profile_name]])

                        #FIXME: The following six lines are reducible to merely:
                        #    self.env_target_inds[profile_name] = ecm_targs_cell + ecm_targs_mem
                        ecm_targs = []
                        for v in ecm_targs_cell:
                            ecm_targs.append(v)
                        for v in ecm_targs_mem:
                            ecm_targs.append(v)
                        self.env_target_inds[profile_name] = ecm_targs

                    # Set the values of Dmems and ECM diffusion based on the
                    # identified target indices.
                    if p.ions_dict['Na'] == 1:
                        dNa = dmem_list['Dm_Na']
                        sim.Dm_cells[sim.iNa][self.tissue_target_inds[profile_name]] = dNa

                    if p.ions_dict['K'] == 1:
                        dK = dmem_list['Dm_K']
                        sim.Dm_cells[sim.iK][self.tissue_target_inds[profile_name]] = dK

                    if p.ions_dict['Cl'] == 1:
                        dCl = dmem_list['Dm_Cl']
                        sim.Dm_cells[sim.iCl][self.tissue_target_inds[profile_name]] = dCl

                    if p.ions_dict['Ca'] == 1:
                        dCa = dmem_list['Dm_Ca']
                        sim.Dm_cells[sim.iCa][self.tissue_target_inds[profile_name]] = dCa

                    if p.ions_dict['H'] == 1:
                        dH = dmem_list['Dm_H']
                        sim.Dm_cells[sim.iH][self.tissue_target_inds[profile_name]] = dH

                    if p.ions_dict['M'] == 1:
                        dM = dmem_list['Dm_M']
                        sim.Dm_cells[sim.iM][self.tissue_target_inds[profile_name]] = dM

                    if p.ions_dict['P'] == 1:
                        dP = dmem_list['Dm_P']
                        sim.Dm_cells[sim.iP][self.tissue_target_inds[profile_name]] = dP

            # Else if this is a cut profile, do nothing. It feeeels good.
            elif profile_type == 'cut':
                pass

            # Else this is a bad profile.
            else:
                TypeError('Profile type {} unrecognized.'.format(profile_type))

    def makeAllChanges(self, sim):
        '''
        Add together all effects to finalize changes to cell membrane
        permeabilities.
        '''

        sim.Dm_cells = \
            sim.Dm_scheduled + \
            sim.Dm_vg + \
            sim.Dm_cag + \
            sim.Dm_morpho + \
            sim.Dm_stretch +\
            sim.Dm_base

        sim.P_cells = sim.P_mod + sim.P_base


#FIXME: Replace "p.scheduled_options['cuts']" everywhere below by "self".
def removeCells(
    tissue_picker,
    sim: 'Simulation', cells: 'Cells', p: 'Parameters') -> None:
    '''
    Permanently remove all cells selected by the passed tissue picker.

    Parameters
    ---------------------------------
    tissue_picker : TissuePicker
        Object matching all cells to be removed.
    sim : Simulator
        Instance of the `Simulator` class.
    cells : Cells
        Instance of the `Cells` class.
    p : Parameters
        Instance of the `Parameters` class.
    '''
    assert types.is_simulator(sim), types.assert_not_simulator(sim)
    assert types.is_cells(cells),   types.assert_not_cells(cells)
    assert types.is_parameters(p),  types.assert_not_parameters(p)

    # Redo environmental diffusion matrices by setting the environmental spaces
    # around cut world to the free value (True) or not (False)?
    open_TJ = True

    # Subtract this bitmap's clipping mask from the global cluster mask.
    bitmap_mask = tissue_picker.get_bitmapper(cells).clipping_matrix
    cells.cluster_mask = cells.cluster_mask - bitmap_mask

    # Indices of all cells to be removed, ignoring extracellular spaces.
    target_inds_cell = tissue_picker.get_cell_indices(cells, p, ignoreECM=True)

    # get the corresponding flags to membrane entities
    target_inds_mem = cells.cell_to_mems[target_inds_cell]
    target_inds_mem,_,_ = tb.flatten(target_inds_mem)
    # target_inds_gj,_,_ = tb.flatten(cells.cell_to_nn_full[target_inds_cell])

    if p.sim_ECM is True:
        # get environmental targets around each removed cell:
        ecm_targs_cell = list(cells.map_cell2ecm[target_inds_cell])
        ecm_targs_mem = list(cells.map_mem2ecm[target_inds_mem])

        #FIXME: The following five lines are reducible to this single line:
        #    ecm_targs = ecm_targs_cell + ecm_targs_mem
        #Multifoliate rose gallavanting obscenely in a midnight haze!
        ecm_targs = []
        for v in ecm_targs_cell:
            ecm_targs.append(v)
        for v in ecm_targs_mem:
            ecm_targs.append(v)

        # redo environmental diffusion matrices by
        # setting the environmental spaces around cut world to the free value -- if desired!:
        if open_TJ is True:
            # save the x,y coordinates of the original boundary cell and membrane points:
            old_bflag_cellxy = cells.cell_centres[cells.bflags_cells]
            # old_bflag_memxy = cells.mem_mids_flat[cells.bflags_mems]

    # set up the situation to make world joined to cut world have more permeable membranes:
    # hurt_cells = np.zeros(len(cells.cell_i))

    # # If we're creating a dangling gap junction situation...
    # if p.scheduled_options['cuts'].is_hurt_cells_leaky:
    #     target_inds_gj_unique = np.unique(target_inds_gj)
    #
    #     for i, inds in enumerate(cells.cell_to_nn_full): # for all the nn inds to a cell...
    #         inds_array = np.asarray(inds)
    #         inds_in_target = np.intersect1d(inds_array,target_inds_gj_unique)
    #
    #         if len(inds_in_target):
    #             hurt_cells[i] = 1  # flag the cell as a "hurt" cell
    #
    #     hurt_inds = (hurt_cells == 1).nonzero()

        #FIXME: The following two conditional branches are suspiciously
        #similar. They appear to reduce to just:
        #
        # cut_indices = hurt_inds if not p.sim_ECM else (
        #     tb.flatten(cells.cell_to_mems[hurt_inds])[0])
        # for i, dmat_a in enumerate(sim.Dm_cells):
        #     sim.Dm_cells[i][cut_indices] = (
        #         p.scheduled_options['cuts'].hurt_cell_leakage * sim.D_free[i])
        #
        #Loamy-eyed surf and sandy-haired sun!
        #
        # if p.sim_ECM is True:
        #     mem_flags,_,_ = tb.flatten(cells.cell_to_mems[hurt_inds])  # get the flags to the memrbanes
        #
        #     for i,dmat_a in enumerate(sim.Dm_cells):
        #         sim.Dm_cells[i][mem_flags] = (
        #             p.scheduled_options['cuts'].hurt_cell_leakage * sim.D_free[i])
        #
        # else:
        #     for i, dmat_a in enumerate(sim.Dm_cells):
        #         sim.Dm_cells[i][hurt_inds] = (
        #             p.scheduled_options['cuts'].hurt_cell_leakage * sim.D_free[i])
        #
        # # copy the Dm to the base:
        # sim.Dm_base = np.copy(sim.Dm_cells)

    # Names of all attributes in the current "Simulation" object.
    sim_names = list(sim.__dict__.keys())

    #FIXME: Redeclare as a "set" object, rename to "specials_names", and remove
    #the duplicate declaration of that variable below. Jumpin' jallopies!

    # Names of all such attributes to be repaired due to this cutting event.
    specials_list = [
        'cc_cells',
        'cc_mems',
        'cc_env',
        'cc_er',
        'z_array',
        'z_array_er',
        'fluxes_gj_x',
        'fluxes_gj_y',
        'fluxes_mem',
        'Dm_base',
        'Dm_cells',
        'Dm_scheduled',
        'Dm_vg',
        'Dm_cag',
        'Dm_morpho',
        'Dm_er_base',
        'Dm_er_CICR',
        'Dm_stretch',
        'D_gj',
    ]

    if p.sim_ECM is True:
        specials_list.remove('cc_env')
        extra = ['z_array_cells']
        for ent in extra:
            specials_list.append(ent)

    special_names = set(specials_list)

    #FIXME: This is awesome-sauce in a BETSE jar, but I don't quite grok what's
    #going on. It'd be swell if this could recieve some documentation
    #massaging. The tidal waves of time recede, inch by minute inch!
    for name in sim_names:
        if name in special_names: # if this is a nested data structure...
            super_data = getattr(sim, name)
            super_data2 = []

            for i, data in enumerate(super_data):
                if isinstance(data, np.ndarray):
                    if len(data) == len(cells.cell_i):
                        data2 = np.delete(data,target_inds_cell)

                    elif len(data) == len(cells.mem_i):
                        data2 = np.delete(data,target_inds_mem)

                    elif len(data) == len(cells.nn_i):
                        data2 = np.delete(data,target_inds_gj)

                elif isinstance(data, list):
                    data2 = []
                    if len(data) == len(cells.cell_i):
                        for index in sorted(target_inds_cell, reverse=True):
                            del data[index]
                        data2.append(data[index])

                    elif len(data) == len(cells.mem_i):
                        for index in sorted(target_inds_mem, reverse=True):
                            del data[index]
                        data2.append(data[index])

                    elif len(data) == len(cells.nn_i):
                        for index in sorted(target_inds_gj, reverse=True):
                            del data[index]
                        data2.append(data[index])

                super_data2.append(data2)

            if type(super_data) == np.ndarray:
                super_data2 = np.asarray(super_data2)

            setattr(sim, name, super_data2)

        #FIXME: What does this "else" branch signify? Logarithmic love indulge!
        else:
            data = getattr(sim, name)

            if isinstance(data, np.ndarray):
                if len(data) == len(cells.cell_i):
                    data2 = np.delete(data, target_inds_cell)
                    setattr(sim, name, data2)

                elif len(data) == len(cells.mem_i):
                    data2 = np.delete(data, target_inds_mem)
                    setattr(sim, name, data2)

            elif isinstance(data, list):
                data2 = []

                if len(data) == len(cells.cell_i):
                    for index in sorted(target_inds_cell, reverse=True):
                        del data[index]
                    data2.append(data[index])
                    setattr(sim, name, data2)

                elif len(data) == len(cells.mem_i):
                    for index in sorted(target_inds_mem, reverse=True):
                        del data[index]
                    data2.append(data[index])
                    setattr(sim, name, data2)


#-------------------------------Fix-up cell world ----------------------------------------------------------------------
    new_cell_centres = []
    new_ecm_verts = []
    removal_flags = np.zeros(len(cells.cell_i))
    removal_flags[target_inds_cell] = 1

    for i,flag in enumerate(removal_flags):
        if flag == 0:
            new_cell_centres.append(cells.cell_centres[i])
            new_ecm_verts.append(cells.ecm_verts[i])

    cells.cell_centres = np.asarray(new_cell_centres)
    cells.ecm_verts = np.asarray(new_ecm_verts)

    # recalculate ecm_verts_unique:
    ecm_verts_flat,_,_ = tb.flatten(cells.ecm_verts)
    ecm_verts_set = set()

    for vert in ecm_verts_flat:
        ptx = vert[0]
        pty = vert[1]
        ecm_verts_set.add((ptx,pty))

    cells.ecm_verts_unique = [list(verts) for verts in list(ecm_verts_set)]
    cells.ecm_verts_unique = np.asarray(cells.ecm_verts_unique)  # convert to numpy array

    #-----------------------------------------------------------------
    logs.log_info('Recalculating cluster variables for new configuration...')

    cells.cellVerts(p)   # create individual cell polygon vertices and other essential data structures
    cells.cellMatrices(p)  # creates a variety of matrices used in routine cells calculations
    cells.mem_processing(p)  # calculates membrane nearest neighbours, ecm interaction, boundary tags, etc
    cells.near_neigh(p)  # Calculate the nn array for each cell
    cells.voronoiGrid(p)
    cells.calc_gj_vects(p)
    cells.environment(p)  # define features of the ecm grid
    cells.make_maskM(p)
    cells.grid_len = len(cells.xypts)

    if p.sim_ECM is True:

        logs.log_info('Recalculating Capacitance Matrix for new configuration...')

        cells.maxwellCapMatrix(p)  # create Maxwell Capacitance Matrix solver for voltages

        if open_TJ is True:
            # if desire for cut away space to lack tight junctions, remove new bflags from set:
            searchTree = sps.KDTree(cells.cell_centres)
            original_pt_inds = list(searchTree.query(old_bflag_cellxy))[1]
            cells.bflags_cells = original_pt_inds[:]

        sim.initDenv(cells,p)

    if p.fluid_flow is True or p.deformation is True or p.calc_J is True:
        # make a laplacian and solver for discrete transfers on closed, irregular cell network:
        logs.log_info('Re-creating cell network Poisson solver...')
        cells.graphLaplacian(p)
        logs.log_info('Completed major world-building computations.')

    if p.sim_eosmosis is True:

        cells.eosmo_tools(p)
