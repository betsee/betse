#!/usr/bin/env python3
# Copyright 2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

# FIXME cut out all the grey areas -- don't need 'em anymore


import numpy as np
from scipy import spatial as sps
from scipy import interpolate as interp
from random import shuffle
from betse.science import toolbox as tb
from betse.exceptions import BetseExceptionSimulation
from betse.science.bitmapper import Bitmapper
from betse.util.io import loggers

class Dynamics(object):

    def __init__(self, sim, cells, p):

        if p.sim_ECM == True:
            self.data_length = len(cells.mem_i)

        elif p.sim_ECM == False:
            self.data_length = len(cells.cell_i)

    def runAllInit(self,sim,cells,p):
        self.globalInit(sim,cells,p)
        self.scheduledInit(sim,cells,p)
        self.dynamicInit(sim,cells,p)

    def runAllDynamics(self,sim,cells,p,t):
        self.globalDyn(sim,cells,p,t)
        self.scheduledDyn(sim,cells,p,t)
        self.dynamicDyn(sim,cells,p,t)
        self.makeAllChanges(sim)

    def globalInit(self,sim,cells,p):

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

            numo = p.global_options['gj_block'][3]

            numo = int(numo)

            if numo > 100:
                numo = 100
            elif numo < 1:
                numo = 1

            data_length = len(cells.nn_index)
            data_fraction = int((numo/100)*data_length)

            shuffle(cells.nn_index)

            self.targets_gj_block = [cells.nn_index[x] for x in range(0,data_fraction)]

        if p.global_options['ecm_change'] != 0:

            self.t_on_ecm = p.global_options['ecm_change'][0]
            self.t_off_ecm = p.global_options['ecm_change'][1]
            self.t_change_ecm = p.global_options['ecm_change'][2]
            self.mult_ecm = p.global_options['ecm_change'][3]

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

    def scheduledInit(self,sim,cells,p):

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

            self.targets_Namem = [item for sublist in self.targets_Namem for item in sublist]

            self.scalar_Namem = 1

            if self.function_Namem != 'None':

                self.scalar_Namem = getattr(tb,self.function_Namem)(cells,sim,p)

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

            self.targets_Kmem = [item for sublist in self.targets_Kmem for item in sublist]

            self.scalar_Kmem = 1

            if self.function_Kmem != 'None':

                self.scalar_Kmem = getattr(tb,self.function_Kmem)(cells,sim,p)

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

            self.targets_Clmem = [item for sublist in self.targets_Clmem for item in sublist]

            self.scalar_Clmem = 1

            if self.function_Clmem != 'None':

                self.scalar_Clmem = getattr(tb,self.function_Clmem)(cells,sim,p)

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

            self.targets_Camem = [item for sublist in self.targets_Camem for item in sublist]

            self.scalar_Camem = 1

            if self.function_Camem != 'None':

                self.scalar_Camem = getattr(tb,self.function_Camem)(cells,sim,p)

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

            self.targets_IP3 = [item for sublist in self.targets_IP3 for item in sublist]

            self.scalar_IP3 = 1

            if self.function_IP3 != 'None':

                self.scalar_IP3 = getattr(tb,self.function_IP3)(cells,sim,p)

        if p.scheduled_options['extV'] != 0 and p.sim_ECM == True:

            self.t_on_extV = p.scheduled_options['extV'][0]
            self.t_off_extV = p.scheduled_options['extV'][1]
            self.t_change_extV = p.scheduled_options['extV'][2]
            self.peak_val_extV = p.scheduled_options['extV'][3]
            self.apply_extV = p.scheduled_options['extV'][4]

            name_positive = self.apply_extV[0]
            name_negative = self.apply_extV[1]

            self.targets_extV_positive = self.env_target_label[name_positive]
            self.targets_extV_negative = self.env_target_label[name_negative]

        if p.scheduled_options['cuts'] != 0 and cells.do_once_cuts == True:

            self.t_cuts = p.scheduled_options['cuts'][0]
            self.apply_cuts = p.scheduled_options['cuts'][1]
            self.cut_hole = p.scheduled_options['cuts'][2]
            self.targets_cuts = []
            targets = self.cuts_target_inds[self.apply_cuts]
            self.targets_cuts.append(targets)

            self.targets_cuts = [item for sublist in self.targets_cuts for item in sublist]

    def dynamicInit(self,sim,cells,p):

        if p.vg_options['Na_vg'] != 0:

            # Initialization of logic values for voltage gated sodium channel
            self.maxDmNa = p.vg_options['Na_vg'][0]
            self.v_activate_Na = p.vg_options['Na_vg'][1]
            self.v_inactivate_Na = p.vg_options['Na_vg'][2]
            self.v_deactivate_Na = p.vg_options['Na_vg'][3]
            self.t_alive_Na = p.vg_options['Na_vg'][4]
            self.t_dead_Na = p.vg_options['Na_vg'][5]
            self.apply_vgNa = p.vg_options['Na_vg'][6]

            # Initialize matrices defining states of vgNa channels for each cell membrane:
            self.inactivated_Na = np.zeros(self.data_length)
            self.vgNa_state = np.zeros(self.data_length)

            self.vgNa_aliveTimer = np.zeros(self.data_length) # sim time at which vgNa starts to close if activated
            self.vgNa_deadTimer = np.zeros(self.data_length) # sim time at which vgNa reactivates after inactivation

            self.targets_vgNa = []
            for profile in self.apply_vgNa:
                targets = self.tissue_target_inds[profile]
                self.targets_vgNa.append(targets)

            self.targets_vgNa = [item for sublist in self.targets_vgNa for item in sublist]
            self.targets_vgNa = np.asarray(self.targets_vgNa)

            self.target_mask_vgNa = np.zeros(self.data_length)
            self.target_mask_vgNa[self.targets_vgNa] = 1

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

            self.targets_Ca = []
            for profile in self.apply_Ca:
                targets = self.tissue_target_inds[profile]
                self.targets_Ca.append(targets)

            self.targets_Ca = [item for sublist in self.targets_Ca for item in sublist]

            self.targets_Ca = np.asarray(self.targets_Ca)

            self.target_mask_Ca = np.zeros(self.data_length)
            self.target_mask_Ca[self.targets_Ca] = 1

    def globalDyn(self,sim,cells,p,t):

        if p.global_options['K_env'] != 0:

            effector_Kenv = tb.pulse(t,self.t_on_Kenv,self.t_off_Kenv,self.t_change_Kenv)

            if p.sim_ECM == False:

                sim.cc_env[sim.iK][:] = self.mem_mult_Kenv*effector_Kenv*p.cK_env + p.cK_env

            elif p.sim_ECM == True: # simulate addition of potassium salt to remain charge neutral

                sim.c_env_bound[sim.iK] = self.mem_mult_Kenv*effector_Kenv*p.env_concs['K'] + p.env_concs['K']
                sim.c_env_bound[sim.iM] = self.mem_mult_Kenv*effector_Kenv*p.env_concs['K'] + p.env_concs['M']

        if p.global_options['Cl_env'] != 0 and p.ions_dict['Cl'] == 1:

            effector_Clenv = tb.pulse(t,self.t_on_Clenv,self.t_off_Clenv,self.t_change_Clenv)

            if p.sim_ECM == False:

                sim.cc_env[sim.iCl][:] = self.mem_mult_Clenv*effector_Clenv*p.cCl_env + p.cCl_env

            elif p.sim_ECM == True:  # simulate addition of sodium chloride to remain charge neutral

                sim.c_env_bound[sim.iCl] = self.mem_mult_Clenv*effector_Clenv*p.env_concs['Cl'] + p.env_concs['Cl']
                sim.c_env_bound[sim.iNa] = self.mem_mult_Clenv*effector_Clenv*p.env_concs['Cl'] + p.env_concs['Na']


        if p.global_options['Na_env'] != 0:

            effector_Naenv = tb.pulse(t,self.t_on_Naenv,self.t_off_Naenv,self.t_change_Naenv)

            if p.sim_ECM == False:

                sim.cc_env[sim.iNa][:] = self.mem_mult_Naenv*effector_Naenv*p.cNa_env + p.cNa_env

            elif p.sim_ECM == True: # simulate addition of sodium salt to remain charge neutral

                sim.c_env_bound[sim.iNa] = self.mem_mult_Naenv*effector_Naenv*p.env_concs['Na'] + p.env_concs['Na']
                sim.c_env_bound[sim.iM] = self.mem_mult_Naenv*effector_Naenv*p.env_concs['Na'] + p.env_concs['M']

        if p.global_options['T_change'] != 0:

            sim.T = self.multT*tb.pulse(t,self.tonT,self.toffT,self.trampT)*p.T + p.T

        if p.global_options['gj_block'] != 0:

            sim.gj_block[self.targets_gj_block] = (1.0 - tb.pulse(t,self.tonGJ,self.toffGJ,self.trampGJ))

        if p.global_options['ecm_change'] != 0:

            effector_ecm = tb.pulse(t,self.t_on_ecm,self.t_off_ecm,self.t_change_ecm)

            for i, subD in enumerate(sim.D_env_base):

                sim.D_env[i][cells.map_cell2ecm] = self.mult_ecm*effector_ecm*subD[cells.map_cell2ecm] + \
                                                   sim.D_env_base[i][cells.map_cell2ecm]

                sim.D_env[i][cells.map_mem2ecm] = self.mult_ecm*effector_ecm*subD[cells.map_mem2ecm] + \
                                                  sim.D_env_base[i][cells.map_mem2ecm]

                # ensure cells on the cluster boundary aren't involved in the change:
                sim.D_env[i][cells.ecm_bound_k] = sim.D_env_base[i][cells.ecm_bound_k]


        if p.global_options['NaKATP_block'] != 0:

            sim.NaKATP_block = (1.0 - tb.pulse(t,self.tonNK,self.toffNK,self.trampNK))

        if p.global_options['HKATP_block'] != 0:

            sim.HKATP_block = (1.0 - tb.pulse(t,self.tonHK,self.toffHK,self.trampHK))

        if p.global_options['VATP_block'] != 0:

            sim.VATP_block = (1.0 - tb.pulse(t,self.tonV,self.toffV,self.trampV))

    def scheduledDyn(self,sim,cells,p,t):

        if p.scheduled_options['Na_mem'] != 0:

            effector_Na = self.scalar_Namem*tb.pulse(t,self.t_on_Namem,self.t_off_Namem,self.t_change_Namem)

            sim.Dm_scheduled[sim.iNa][self.targets_Namem] = self.mem_mult_Namem*effector_Na*p.Dm_Na

        if p.scheduled_options['K_mem'] != 0:

            effector_K = self.scalar_Kmem*tb.pulse(t,self.t_on_Kmem,self.t_off_Kmem,self.t_change_Kmem)

            sim.Dm_scheduled[sim.iK][self.targets_Kmem] = self.mem_mult_Kmem*effector_K*p.Dm_K

        if p.scheduled_options['Cl_mem'] != 0 and p.ions_dict['Cl'] != 0:

            effector_Cl = self.scalar_Clmem*tb.pulse(t,self.t_on_Clmem,self.t_off_Clmem,self.t_change_Clmem)

            sim.Dm_scheduled[sim.iCl][self.targets_Clmem] = self.mem_mult_Clmem*effector_Cl*p.Dm_Cl

        if p.scheduled_options['Ca_mem'] != 0 and p.ions_dict['Ca'] != 0:

            effector_Ca = self.scalar_Camem*tb.pulse(t,self.t_on_Camem,self.t_off_Camem,self.t_change_Camem)

            sim.Dm_scheduled[sim.iCa][self.targets_Camem] = self.mem_mult_Camem*effector_Ca*p.Dm_Ca

        if p.scheduled_options['IP3'] != 0:

            sim.cIP3[self.targets_IP3] = sim.cIP3[self.targets_IP3] + self.scalar_IP3*self.rate_IP3*tb.pulse(t,self.t_onIP3,
                self.t_offIP3,self.t_changeIP3)

        if p.scheduled_options['cuts'] != 0 and cells.do_once_cuts == True and t>self.t_cuts:

            target_method = p.tissue_profiles[self.apply_cuts]['target method']

            removeCells(self.apply_cuts,target_method,sim,cells,p,cavity_volume=self.cut_hole,simMod = True)

            # redo main data length variable for this dynamics module with updated world:
            if p.sim_ECM == True:

                self.data_length = len(cells.mem_i)

            elif p.sim_ECM == False:

                self.data_length = len(cells.cell_i)

            self.tissueProfiles(sim,cells,p)

            cells.redo_gj(self,p,savecells =False)

            self.runAllInit(sim,cells,p)

            if p.plot_while_solving == True:

                sim.checkPlot.resetData(cells,sim,p)

            cells.do_once_cuts = False  # set the cells' do_once field to prevent attempted repeats

        if p.scheduled_options['extV'] != 0 and p.sim_ECM == True:

            effector_extV = tb.pulse(t,self.t_on_extV,self.t_off_extV,self.t_change_extV)

            sim.bound_V[self.targets_extV_positive] = self.peak_val_extV*effector_extV
            sim.bound_V[self.targets_extV_negative] = -self.peak_val_extV*effector_extV

    def dynamicDyn(self,sim,cells,p,t):

        self.dvsign = np.sign(sim.dvm)

        if p.vg_options['Na_vg'] != 0:

            self.vgSodium(sim,cells,p,t)

        if p.vg_options['K_vg'] !=0:

            self.vgPotassium(sim,cells,p,t)

        if p.vg_options['Ca_vg'] !=0 and p.ions_dict['Ca'] != 0:

            self.vgCalcium(sim,cells,p,t)

        if p.vg_options['K_cag'] != 0 and p.ions_dict['Ca'] != 0:

            self.cagPotassium(sim,cells,p,t)

        if p.Ca_dyn_options['CICR'] != 0 and p.ions_dict['Ca'] != 0:

            self.calciumDynamics(sim,cells,p)

    def vgSodium(self,sim,cells,p,t):

        # Logic phase 1: find out which cells have activated their vgNa channels
        truth_vmGTvon_Na = sim.vm > self.v_activate_Na  # returns bools of vm that are bigger than threshhold
        #truth_depol_Na = dvsign==1  # returns bools of vm that are bigger than threshhold
        truth_not_inactivated_Na = self.inactivated_Na == 0  # return bools of vm that can activate
        truth_vgNa_Off = self.vgNa_state == 0 # hasn't been turned on yet

        # find the cell indicies that correspond to all statements of logic phase 1:
        inds_activate_Na = (truth_vmGTvon_Na*truth_not_inactivated_Na*truth_vgNa_Off*self.target_mask_vgNa).nonzero()

        self.vgNa_state[inds_activate_Na] = 1 # open the channel
        self.vgNa_aliveTimer[inds_activate_Na] = t + self.t_alive_Na # set the timers for the total active state
        self.vgNa_deadTimer[inds_activate_Na] = 0  # reset any timers for an inactivated state to zero

        # Logic phase 2: find out which cells have closed their gates due to crossing inactivating voltage:
        truth_vgNa_On = self.vgNa_state == 1  # channel must be on already
        truth_vmGTvoff_Na = sim.vm > self.v_inactivate_Na  # bools of cells that have vm greater than shut-off volts

        inds_inactivate_Na = (truth_vgNa_On*truth_vmGTvoff_Na*self.target_mask_vgNa).nonzero()

        self.vgNa_state[inds_inactivate_Na] = 0    # close the vg sodium channels
        self.inactivated_Na[inds_inactivate_Na] = 1   # switch these so cells do not re-activate
        self.vgNa_aliveTimer[inds_inactivate_Na] = 0            # reset any alive timers to zero
        self.vgNa_deadTimer[inds_inactivate_Na] = t + self.t_dead_Na # set the timer of the inactivated state

         # Logic phase 3: find out if cell activation state has timed out, also rendering inactivated state:

        truth_vgNa_act_timeout = self.vgNa_aliveTimer < t   # find cells that have timed out their vgNa open state
        truth_vgNa_On = self.vgNa_state == 1 # ensure the vgNa is indeed open
        inds_timeout_Na_act = (truth_vgNa_act_timeout*truth_vgNa_On*self.target_mask_vgNa).nonzero()

        self.vgNa_state[inds_timeout_Na_act] = 0             # set the state to closed
        self.vgNa_aliveTimer[inds_timeout_Na_act] = 0            # reset the timers to zero
        self.inactivated_Na[inds_timeout_Na_act] = 1    # inactivate the channel so it can't reactivate
        self.vgNa_deadTimer[inds_timeout_Na_act] = t + self.t_dead_Na # set the timer of the inactivated state

        # Logic phase 4: find out if inactivation timers have timed out:
        truth_vgNa_inact_timeout = self.vgNa_deadTimer <t  # find cells that have timed out their vgNa inact state
        truth_vgNa_Off = self.vgNa_state == 0 # check to make sure these channels are indeed closed
        inds_timeout_Na_inact = (truth_vgNa_inact_timeout*truth_vgNa_Off*self.target_mask_vgNa).nonzero()

        self.vgNa_deadTimer[inds_timeout_Na_inact] = 0    # reset the inactivation timer
        self.inactivated_Na[inds_timeout_Na_inact] = 0    # remove inhibition to activation

        # Logic phase 5: find out if cells have passed below threshhold to become deactivated:
        truth_vmLTvreact_Na = sim.vm < self.v_deactivate_Na # voltage is lower than the deactivate voltage

        inds_deactivate_Na = (truth_vmLTvreact_Na*self.target_mask_vgNa).nonzero()

        self.inactivated_Na[inds_deactivate_Na] = 0  # turn any inhibition to activation off
        self.vgNa_state[inds_deactivate_Na] = 0   # shut the Na channel off if it's on
        self.vgNa_aliveTimer[inds_deactivate_Na] = 0       # reset any alive-timers to zero
        self.vgNa_deadTimer[inds_deactivate_Na] = 0   # reset any dead-timers to zero

        # Define ultimate activity of the vgNa channel:
        sim.Dm_vg[sim.iNa] = self.maxDmNa*self.vgNa_state

    def vgPotassium(self,sim,cells,p,t):
         # detecting channels to turn on:

        truth_vmGTvon_K = sim.vm > self.v_on_K  # bools for cells with vm greater than the on threshold for vgK
        truth_depol_K = self.dvsign == 1  # bools matrix for cells that are depolarizing
        truth_vgK_OFF = self.vgK_state == 0   # bools matrix for cells that are in the off state

        # cells at these indices will become activated in this time step:
        inds_activate_K = (truth_vmGTvon_K*truth_depol_K*truth_vgK_OFF*self.target_mask_vgK).nonzero()
        self.vgK_state[inds_activate_K] = 1  # set the state of these channels to "open"
        self.vgK_OFFtime[inds_activate_K] = self.t_alive_K + t  # set the time at which these channels will close

        #  detecting channels to turn off:
        truth_vgK_ON = self.vgK_state == 1  # detect cells that are in their on state
        truth_vgK_timeout = self.vgK_OFFtime < t     # detect the cells that have expired off timers
        inds_deactivate_K = (truth_vgK_ON*truth_vgK_timeout*self.target_mask_vgK).nonzero()
        self.vgK_state[inds_deactivate_K] = 0 # turn off the channels to closed
        self.vgK_OFFtime[inds_deactivate_K] = 0

        inds_open_K = (self.vgK_state == 1).nonzero()
        self.active_K[inds_open_K] = 1

        inds_closed_K =(self.vgK_state == 0).nonzero()
        self.active_K[inds_closed_K] = 0

        sim.Dm_vg[sim.iK] = self.maxDmK*self.active_K

    def vgCalcium(self,sim,cells,p,t):
         # detect condition to turn vg_Ca channel on:
        truth_vmGTvon_Ca = sim.vm > self.v_on_Ca  # bools for cells with vm greater than the on threshold for vgK
        truth_caLTcaOff = sim.cc_cells[sim.iCa] < self.ca_lower_ca # check that cellular calcium is below inactivating Ca
        truth_depol_Ca = self.dvsign == 1  # bools matrix for cells that are depolarizing
        truth_vgCa_OFF = self.vgCa_state == 0   # bools matrix for cells that are in the off state

        # cells at these indices will become activated in this time step:
        inds_activate_Ca = (truth_vmGTvon_Ca*truth_depol_Ca*truth_caLTcaOff*truth_vgCa_OFF*self.target_mask_vgCa).nonzero()
        self.vgCa_state[inds_activate_Ca] = 1  # set the state of these channels to "open"

        # detect condition to turn off vg_Ca channel:
        truth_caGTcaOff = sim.cc_cells[sim.iCa] > self.ca_upper_ca   # check that calcium exceeds maximum
        truth_vgCa_ON = self.vgCa_state == 1 # check that the channel is on
        inds_inactivate_Ca = (truth_caGTcaOff*truth_vgCa_ON*self.target_mask_vgCa).nonzero()
        self.vgCa_state[inds_inactivate_Ca] = 0

        # additional condition to turn off vg_Ca via depolarizing voltage:
        truth_vmGTvcaOff = sim.vm > self.v_off_Ca
        inds_inactivate_Ca_2 = (truth_vmGTvcaOff*self.target_mask_vgCa*truth_vgCa_ON).nonzero()
        self.vgCa_state[inds_inactivate_Ca_2] = 0

        inds_open_Ca = (self.vgCa_state == 1).nonzero()
        self.active_Ca[inds_open_Ca] = 1

        inds_closed_Ca =(self.vgCa_state == 0).nonzero()
        self.active_Ca[inds_closed_Ca] = 0

        sim.Dm_vg[sim.iCa] = self.maxDmCa*self.active_Ca

    def cagPotassium(self,sim,cells,p,t):

        # inds_cagK_targets = (self.targets_cagK).nonzero()
        #
        # self.active_cagK[inds_cagK_targets] = tb.hill(sim.cc_cells[sim.iCa][inds_cagK_targets],
        #     self.Kcag_halfmax,self.Kcag_n)

        self.active_cagK[self.targets_cagK] = tb.hill(sim.cc_cells[sim.iCa][self.targets_cagK],
            self.Kcag_halfmax,self.Kcag_n)

        sim.Dm_cag[sim.iK] = self.maxDmKcag*self.active_cagK

    def calciumDynamics(self,sim,cells,p):

        if p.Ca_dyn_options['CICR'] != 0:

            dcc_CaER_sign = np.sign(sim.dcc_ER[0])

            if len(p.Ca_dyn_options['CICR'][1])==0:
                term_Ca_reg = 1.0

            else:
                term_Ca_reg = (np.exp(-((sim.cc_cells[sim.iCa]-self.midCaR)**2)/((2*self.widthCaR)**2)))

            if len(p.Ca_dyn_options['CICR'][2]) == 0:
                term_IP3_reg = 1.0

            else:
                term_IP3_reg = tb.hill(sim.cIP3,self.KhmIP3,self.n_IP3)

            if p.FMmod == 1:
                span = self.topCa - self.bottomCa
                FMmod = p.ip3FM*span
                topCa = self.topCa - FMmod*term_IP3_reg
            else:
                topCa = self.topCa

            truth_overHighCa = sim.cc_er[0] >=  topCa
            truth_increasingCa = dcc_CaER_sign == 1
            truth_alreadyClosed = self.stateER == 0.0
            inds_open_ER = (truth_overHighCa*truth_increasingCa*truth_alreadyClosed).nonzero()

            truth_underBottomCa = sim.cc_er[0]< self.bottomCa
            truth_decreasingCa = dcc_CaER_sign == -1
            truth_alreadyOpen = self.stateER == 1.0
            inds_close_ER = (truth_underBottomCa*truth_alreadyOpen).nonzero()

            self.stateER[inds_open_ER] = 1.0
            self.stateER[inds_close_ER] = 0.0

            sim.Dm_er_CICR[0] = self.maxDmCaER*self.stateER*term_IP3_reg*term_Ca_reg

            sim.Dm_er = sim.Dm_er_CICR + sim.Dm_er_base

    def tissueProfiles(self,sim,cells,p):

        """
        Reads in parameters data to build cell and membrane specific (if p.sim_ECM == True) index
        sets for each user-defined tissue profile.

        """
        profile_names = list(p.tissue_profiles.keys())
        self.tissue_target_inds = {}
        self.cell_target_inds = {}
        self.cuts_target_inds = {}
        self.env_target_inds = {}
        self.tissue_profile_names = []


        # Do a first search to find requests for any cavities as they modify the world structure:

        if cells.do_once_cavity == True:  # if we haven't done this already...
            for name in profile_names:

                data_stream = p.tissue_profiles[name]
                target_method = data_stream['target method']
                dmem_list = data_stream['diffusion constants']
                ecm_val = data_stream['ecm multiplier']
                designation = data_stream['designation']

                if designation == 'cavity':

                    removeCells(name,target_method,sim,cells,p,cavity_volume=True)
                    cells.do_once_cavity = False  # set the cells' do_once field to prevent attempted repeats

        # Go through again and do traditional tissue profiles:
        for name in profile_names:

            data_stream = p.tissue_profiles[name]
            target_method = data_stream['target method']
            dmem_list = data_stream['diffusion constants']
            ecm_val = data_stream['ecm multiplier']
            designation = data_stream['designation']

            if designation == 'None':

                self.tissue_profile_names.append(name)

                self.tissue_target_inds[name] = getCellTargets(name,target_method, cells, p)
                self.cell_target_inds[name] = getCellTargets(name,target_method, cells, p, ignoreECM=True)

                if p.sim_ECM == True:
                    #get ecm targets

                    ecm_targs_cell = list(cells.map_cell2ecm[self.cell_target_inds[name]])
                    ecm_targs_mem = list(cells.map_mem2ecm[self.tissue_target_inds[name]])

                    ecm_targs = []

                    for v in ecm_targs_cell:
                        ecm_targs.append(v)

                    for v in ecm_targs_mem:
                        ecm_targs.append(v)

                    self.env_target_inds[name] = ecm_targs

                # set the values of Dmems and ecm diffusion based on the identified target indices
                if p.ions_dict['Na'] == 1:
                    dNa = dmem_list['Dm_Na']
                    sim.Dm_cells[sim.iNa][self.tissue_target_inds[name]] = dNa

                    if p.sim_ECM == True:
                        sim.D_env[sim.iNa][self.env_target_inds[name]] = p.Do_Na*ecm_val

                if p.ions_dict['K'] == 1:
                    dK = dmem_list['Dm_K']
                    sim.Dm_cells[sim.iK][self.tissue_target_inds[name]] = dK

                    if p.sim_ECM == True:
                        sim.D_env[sim.iK][self.env_target_inds[name]] = p.Do_K*ecm_val

                if p.ions_dict['Cl'] == 1:
                    dCl = dmem_list['Dm_Cl']
                    sim.Dm_cells[sim.iCl][self.tissue_target_inds[name]] = dCl

                    if p.sim_ECM == True:
                        sim.D_env[sim.iCl][self.env_target_inds[name]] = p.Do_Cl*ecm_val

                if p.ions_dict['Ca'] == 1:
                    dCa = dmem_list['Dm_Ca']
                    sim.Dm_cells[sim.iCa][self.tissue_target_inds[name]] = dCa

                    if p.sim_ECM == True:

                        sim.D_env[sim.iCa][self.env_target_inds[name]] = p.Do_Ca*ecm_val

                if p.ions_dict['H'] == 1:
                    dH = dmem_list['Dm_H']
                    sim.Dm_cells[sim.iH][self.tissue_target_inds[name]] = dH

                    if p.sim_ECM == True:

                        sim.D_env[sim.iH][self.env_target_inds[name]] = p.Do_H*ecm_val

                if p.ions_dict['M'] == 1:
                    dM = dmem_list['Dm_M']
                    sim.Dm_cells[sim.iM][self.tissue_target_inds[name]] = dM

                    if p.sim_ECM == True:

                        sim.D_env[sim.iM][self.env_target_inds[name]] = p.Do_M*ecm_val
                        # sim.D_env[sim.iM][ecm_targs_mem] = p.Do_M*ecm_val
                        # sim.D_env[sim.iM][ecm_targs_cell] = p.Do_M*ecm_val

                if p.ions_dict['P'] == 1:
                    dP = dmem_list['Dm_P']
                    sim.Dm_cells[sim.iP][self.tissue_target_inds[name]] = dP

                    if p.sim_ECM == True:

                        sim.D_env[sim.iP][self.env_target_inds[name]] = p.Do_P*ecm_val
                        # sim.D_env[sim.iP][ecm_targs_mem] = p.Do_P*ecm_val
                        # sim.D_env[sim.iP][ecm_targs_cell] = p.Do_P*ecm_val

            elif designation == 'cuts':
                # if the user wants to use this as a region to be cut, define cuts target inds:
                self.cuts_target_inds[name] = getCellTargets(name,target_method, cells, p,ignoreECM=True)

    def makeAllChanges(self,sim):
        # Add together all effects to make change on the cell membrane permeabilities:
        sim.Dm_cells = sim.Dm_scheduled + sim.Dm_vg + sim.Dm_cag + sim.Dm_base

def getCellTargets(profile_key,targets_description,cells,p,ignoreECM = False):

    """
    Using an input description flag, which is a string in the format of
    'random40', a list of integers corresponding to cell indices,
    [4,5,7], this returns the cell or membrane indices to used to define
    tissue profiles.The string format targets the specified random fraction of total
    indices, for instance 'random20' would randomly select 20% of the cell population.

    Parameters
    ---------------------------------
    targets_description                  a string in the format 'random50', 'all', or a list of indices to cell_i
    cells                                an instance of the world module object
    p                                    an instance of the parameters module object
    ignoreECM                            a flag telling the function to ignore p.sim_ECM

    Returns
    ---------------------------------
    target_inds                          a list of integers corresponding to targeted cell or membrane indices

    """

    if isinstance(targets_description,str):

        meter = len(targets_description)

        if meter >= 6:

            chaff = targets_description[0:6]
            numo = targets_description[6:len(targets_description)]

            if chaff == 'bitmap':

                if p.use_bitmaps == True:

                    bitmask = Bitmapper(p,profile_key,cells.xmin,cells.xmax,cells.ymin,cells.ymax)
                    bitmask.clipPoints(cells.cell_centres[:,0],cells.cell_centres[:,1])
                    target_inds = bitmask.good_inds   # get the cell_i indicies falling within the bitmap mask

                    if p.sim_ECM == True and ignoreECM == False:

                        target_inds = cells.cell_to_mems[target_inds]
                        target_inds,_,_ = tb.flatten(target_inds)

                else:

                    target_inds = []

            elif chaff == 'bounda':

                if p.sim_ECM == False or ignoreECM == True:
                    target_inds = cells.bflags_cells

                elif p.sim_ECM == True and ignoreECM == False:
                    target_inds = cells.cell_to_mems[cells.bflags_cells]
                    target_inds,_,_ = tb.flatten(target_inds)


            elif chaff == 'random':

                numo = int(numo)

                if numo > 100:
                    numo = 100
                elif numo < 1:
                    numo = 1

                data_length = len(cells.cell_i)
                data_fraction = int((numo/100)*data_length)

                shuffle(cells.cell_i)

                target_inds_cell = [cells.cell_i[x] for x in range(0,data_fraction)]

                if p.sim_ECM == False or ignoreECM == True:
                    target_inds = target_inds_cell

                elif p.sim_ECM == True and ignoreECM == False:
                    target_inds = cells.cell_to_mems[target_inds_cell]
                    target_inds,_,_ = tb.flatten(target_inds)

        elif targets_description == 'all':

            if p.sim_ECM == False or ignoreECM == True:
                target_inds = cells.cell_i

            elif p.sim_ECM == True and ignoreECM == False:
                target_inds = cells.cell_to_mems[cells.cell_i]
                target_inds,_,_ = tb.flatten(target_inds)

        else:
            raise BetseExceptionSimulation("Error in specifying cell targets for tissue profile."
            "String must be in an appropriate format: 'random10', 'all', 'boundary', 'bitmap'")


    elif isinstance(targets_description, list):

        target_inds_cell = targets_description

        if p.sim_ECM == False or ignoreECM == True:
            target_inds = target_inds_cell

        elif p.sim_ECM == True and ignoreECM == False:
            target_inds = cells.cell_to_mems[target_inds_cell]
            target_inds,_,_ = tb.flatten(target_inds)

    return target_inds

def getEcmTargets(profile_key,targets_description,cells,p,boundaryOnly = True):

    """
    Using an input description flag, which is a string in the format of
    'random40', or a list of integers corresponding to ecm indices,
    [4,5,7], this returns the ecm indices to used to define
    tissue profiles.The string format targets the specified random fraction of total
    indices, for instance 'random20' would randomly select 20% of the ecm spaces.

    Parameters
    ---------------------------------
    targets_description                  a list [8,9,10] of indices to bflags_ecm or ecm_i or the string 'all'
    cells                                an instance of the world module object
    p                                    an instance of the parameters module object
    boundaryOnly                         a flag telling the function we're only interested in bflags_ecm

    Returns
    ---------------------------------
    target_inds                          a list of integers corresponding to targeted ecm indices

    """

    target_inds = []

    if isinstance(targets_description,str):

        inds_env = np.asarray(len(cells.xypts))

        if targets_description == 'bitmap':

            if p.use_bitmaps == True:
                bitmask = Bitmapper(p,profile_key,cells.xmin,cells.xmax,cells.ymin,cells.ymax)
                bitmask.clipPoints(cells.xypts[:,0],cells.xypts[:,1])
                target_inds = bitmask.good_inds   # get the cell_i indicies falling within the bitmap mask

    #     elif targets_description == 'all':
    #
    #         target_inds = inds_env
    #
    # if isinstance(targets_description, list):
    #
    #     inds_cell = cells.map_cell2ecm[targets_description]
    #     inds_mem = cells.map_mem2ecm[targets_description]
    #
    #     target_inds = []
    #
    #     target_inds.append(inds_cell)
    #     target_inds.append(inds_mem)

    return target_inds

def removeCells(profile_name,targets_description,sim,cells,p,cavity_volume = False, simMod = False):

    loggers.log_info('Cutting hole in cell cluster! Removing cells...')

    if isinstance(targets_description,str):

        if p.use_bitmaps == True and targets_description == 'bitmap':

            bitmask = Bitmapper(p,profile_name,cells.xmin,cells.xmax,cells.ymin,cells.ymax)
            bitmask.clipPoints(cells.cell_centres[:,0],cells.cell_centres[:,1])
            target_inds_cell = bitmask.good_inds   # get the cell_i indices falling within the bitmap mask

            cells.cluster_mask = cells.cluster_mask - bitmask.clippingMatrix  # update the cluster mask by subtracting deleted region

    if isinstance(targets_description,list):   # otherwise, if it's a list, then take the targets literally...

        target_inds_cell = targets_description

    # get the corresponding flags to membrane entities
    target_inds_mem = cells.cell_to_mems[target_inds_cell]
    target_inds_mem,_,_ = tb.flatten(target_inds_mem)
    target_inds_gj,_,_ = tb.flatten(cells.cell_to_nn[target_inds_cell])


    if simMod == True:

        sim_names = list(sim.__dict__.keys())
        specials_list = ['cc_cells','cc_env','z_array','z_array_er','Dm_cells','fluxes_gj_x','fluxes_gj_y',
            'fluxes_mem','Dm_base','Dm_scheduled','Dm_vg','Dm_cag','Dm_er_base','Dm_er_CICR','D_gj','cc_er']

        if p.sim_ECM == True:
            specials_list.remove('cc_env')
            extra = ['z_array_cells']
            for ent in extra:
                specials_list.append(ent)


        special_names = set(specials_list)

        for name in sim_names:

            if name in special_names: # if this is a nested data structure...

                super_data = getattr(sim,name)

                super_data2 = []

                for i, data in enumerate(super_data):

                    if isinstance(data,np.ndarray):
                        if len(data) == len(cells.cell_i):
                            data2 = np.delete(data,target_inds_cell)

                        elif len(data) == len(cells.mem_i):
                            data2 = np.delete(data,target_inds_mem)

                        elif len(data) == len(cells.nn_index):
                            data2 = np.delete(data,target_inds_gj)


                    if isinstance(data,list):
                        data2 = []
                        if len(data) == len(cells.cell_i):
                            for index in sorted(target_inds_cell, reverse=True):
                                del data[index]
                            data2.append(data[index])

                        elif len(data) == len(cells.mem_i):
                            for index in sorted(target_inds_mem, reverse=True):
                                del data[index]
                            data2.append(data[index])

                        elif len(data) == len(cells.nn_index):
                            for index in sorted(target_inds_gj, reverse=True):
                                del data[index]
                            data2.append(data[index])


                    super_data2.append(data2)

                if type(super_data) == np.ndarray:
                    super_data2 = np.asarray(super_data2)

                setattr(sim,name,super_data2)

            else:

                data = getattr(sim,name)

                if isinstance(data,np.ndarray):
                    if len(data) == len(cells.cell_i):
                        data2 = np.delete(data,target_inds_cell)
                        setattr(sim,name,data2)

                    elif len(data) == len(cells.mem_i):
                        data2 = np.delete(data,target_inds_mem)
                        setattr(sim,name,data2)

                    elif len(data) == len(cells.nn_index):
                        data2 = np.delete(data,target_inds_gj)
                        setattr(sim,name,data2)

                if isinstance(data,list):
                    data2 = []

                    if len(data) == len(cells.cell_i):
                        for index in sorted(target_inds_cell, reverse=True):
                            del data[index]
                        data2.append(data[index])
                        setattr(sim,name,data2)

                    elif len(data) == len(cells.mem_i):
                        for index in sorted(target_inds_mem, reverse=True):
                            del data[index]
                        data2.append(data[index])
                        setattr(sim,name,data2)

                    elif len(data) == len(cells.nn_index):
                        for index in sorted(target_inds_gj, reverse=True):
                            del data[index]
                        data2.append(data[index])
                        setattr(sim,name,data2)


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

    loggers.log_info('Recalculating cluster variables for new configuration...')

    if p.sim_ECM == True:

        cells.cellVerts(p)   # create individual cell polygon vertices
        cells.bflags_mems,_ = cells.boundTag(cells.mem_mids_flat,p,alpha=0.8)  # flag membranes on the cluster bound
        cells.near_neigh(p)    # Calculate the nn array for each cell
        cells.cleanUp(p)       # Free up memory...
        # cells.makeECM(p)       # create the ecm grid
        cells.environment(p)   # define features of the ecm grid
        cells.grid_len =len(cells.xypts)

        # make a laplacian and solver for discrete transfers on closed, irregular cell network:
        loggers.log_info('Creating cell network Poisson solver...')
        cells.graphLaplacian()
        loggers.log_info('Completed major world-building computations.')

    else:

        cells.cellVerts(p)   # create individual cell polygon vertices and membrane specific data structures
        cells.bflags_cells,_ = cells.boundTag(cells.cell_centres,p,alpha=0.8)
        cells.near_neigh(p)    # Calculate the nn array for each cell
        cells.cleanUp(p)      # Free up memory...














