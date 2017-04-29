#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

# FIXME include other channels in morphogen (dye) dynamics

# ....................{ IMPORTS                            }....................
import copy
import numpy as np
from betse.science.channels import vg_ca as vgca
from betse.science.channels import vg_funny as vgfun
from betse.science.channels import vg_k as vgk
from betse.science.channels import vg_kir as vgkir
from betse.science.channels import vg_na as vgna
from betse.science.channels import vg_nap as vgnap
from betse.science.channels import wound_channel as w
from betse.science.math import modulate as mod
from betse.science.math import toolbox as tb
from betse.science.tissue.channels_o import cagPotassium
from betse.util.io.log import logs
from betse.util.type import types
from betse.util.type.types import type_check
from random import shuffle
from scipy import interpolate as interp
from scipy import spatial as sps

# ....................{ CLASSES                            }....................
class TissueHandler(object):
    '''
    High-level handler for user-specified tissue-centric functionality,
    including both tissue profiles *and* scheduled interventions.

    This handler governs all:

    * Tissue profiles and objects required by these profiles, including all
      geometry-specifying bitmaps.
    * Scheduled interventions, even those *not* pertaining to tissue profiles
      (e.g., global scheduled interventions).

    Attributes (General)
    ----------------------------
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(
        self,
        sim:   'betse.science.sim.Simulator',
        cells: 'betse.science.cells.Cells',
        p:     'betse.science.parameters.Parameters',
    ) -> None:

        if p.sim_ECM:
            self.data_length = len(cells.mem_i)
        else:
            self.data_length = len(cells.mem_i)

        self.wound_channel_used = False

    # ..................{ RUNNERS ~ init                     }..................
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

    # ..................{ RUNNERS ~ apply                    }..................
    def runAllDynamics(self, sim, cells, p, t: float):
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
                self.scalar_Namem, self.dyna_Namem = getattr(
                    mod, self.function_Namem)(self.targets_Namem,cells,p)

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

                self.scalar_Kmem, self.dyna_Kmem = getattr(
                    mod, self.function_Kmem)(self.targets_Kmem, cells, p)

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
                self.scalar_P, self.dyna_P = getattr(mod, self.function_P)(
                    self.targets_P,cells, p)

        #--------------------------------------------------------

        if p.scheduled_options['ecmJ'] != 0 and p.sim_ECM is True:
            self.t_on_ecmJ  = p.scheduled_options['ecmJ'][0]
            self.t_off_ecmJ = p.scheduled_options['ecmJ'][1]
            self.t_change_ecmJ = p.scheduled_options['ecmJ'][2]
            self.apply_ecmJ = p.scheduled_options['ecmJ'][3]
            self.mult_ecmJ = p.scheduled_options['ecmJ'][4]

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

        # voltage gated sodium channel --------------------------------------------------------------------------------
        if p.vgNa_bool:

            # Initialization of logic values for voltage gated sodium channel
            self.maxDmNa = p.vgNa_max

            self.apply_vgNa = p.vgNa_apply

            self.targets_vgNa = []

            for profile in self.apply_vgNa:
                targets = self.tissue_target_inds[profile]
                self.targets_vgNa.append(targets)

            self.targets_vgNa = [item for sublist in self.targets_vgNa for item in sublist]
            self.targets_vgNa = np.asarray(self.targets_vgNa)

            # create the desired voltage gated sodium channel instance:
            for ind, type_string in enumerate(p.vgNa_type):

                Na_class = getattr(vgna,type_string,'Nav1p2')

                object_name = 'vgNa_object' + str(ind)

                setattr(self,object_name,Na_class())

                if p.run_sim is True:
                    # initialize the voltage-gated sodium object
                    init_funk = getattr(self,object_name)
                    init_funk.init(self, sim, cells, p)

        # persistent voltage gated sodium channel:---------------------------------------------------------------------
        if p.vgNaP_bool:

            # Initialization of logic values for persistent voltage gated sodium channel
            self.maxDmNaP = p.vgNaP_max

            self.apply_vgNaP = p.vgNaP_apply

            self.targets_vgNaP = []

            for profile in self.apply_vgNaP:
                targets = self.tissue_target_inds[profile]
                self.targets_vgNaP.append(targets)

            self.targets_vgNaP = [item for sublist in self.targets_vgNaP for item in sublist]
            self.targets_vgNaP = np.asarray(self.targets_vgNaP)

            # create the desired voltage gated sodium channel instances:

            for ind, type_string in enumerate(p.vgNaP_type):

                NaP_class = getattr(vgnap,type_string,'Nav1p6')

                object_name = 'vgNaP_object' + str(ind)

                setattr(self,object_name,NaP_class())

                if p.run_sim is True:
                    # initialize the voltage-gated sodium object
                    init_funk = getattr(self,object_name)
                    init_funk.init(self, sim, cells, p)


            # NaP_class_ = getattr(vgnap,p.vgNaP_type,'Nav1p6')
            # self.vgNaP_object = NaP_class_()
            #
            # if p.run_sim is True:
            #     # initialize the voltage-gated sodium object
            #     self.vgNaP_object.init(self, sim, cells, p)

        # voltage gated potassium channel ----------------------------------------------------------------------------
        if p.vgK_bool:

            # Initialization of logic values forr voltage gated potassium channel
            self.maxDmK = p.vgK_max

            self.apply_vgK = p.vgK_apply

            self.targets_vgK = []

            for profile in self.apply_vgK:
                targets = self.tissue_target_inds[profile]
                self.targets_vgK.append(targets)

            self.targets_vgK = [item for sublist in self.targets_vgK for item in sublist]
            self.targets_vgK = np.asarray(self.targets_vgK)

            # create the desired voltage gated potassium channel instances:
            for ind, type_string in enumerate(p.vgK_type):

                K_class = getattr(vgk,type_string,'Kv1p2')

                object_name = 'vgK_object' + str(ind)

                setattr(self,object_name,K_class())

                if p.run_sim is True:
                    # initialize the voltage-gated sodium object
                    init_funk = getattr(self,object_name)
                    init_funk.init(self, sim, cells, p)


            # K_class_ = getattr(vgk,p.vgK_type,'Kv1p2')
            # self.vgK_object = K_class_()
            #
            # if p.run_sim is True:
            #     # initialize the voltage-gated sodium object
            #     self.vgK_object.init(self, sim, cells, p)

        # Inward rectifying voltage gated potassium channel -----------------------------------------------------------
        if p.vgKir_bool:

            # Initialization of logic values for voltage gated potassium channel
            self.maxDmKir = p.vgKir_max

            self.apply_vgKir = p.vgKir_apply

            self.targets_vgKir = []

            for profile in self.apply_vgKir:
                targets = self.tissue_target_inds[profile]
                self.targets_vgKir.append(targets)

            self.targets_vgKir = [item for sublist in self.targets_vgKir for item in sublist]
            self.targets_vgKir = np.asarray(self.targets_vgKir)

            # create the desired voltage gated potassium channel instances:
            for ind, type_string in enumerate(p.vgKir_type):

                Kir_class = getattr(vgkir,type_string,'Kir2p1')

                object_name = 'vgKir_object' + str(ind)

                setattr(self,object_name,Kir_class())

                if p.run_sim is True:
                    # initialize the voltage-gated potassium object
                    init_funk = getattr(self,object_name)
                    init_funk.init(self, sim, cells, p)

            # Kir_class_ = getattr(vgkir,p.vgKir_type,'Kir2p1')
            # self.vgKir_object = Kir_class_()
            #
            # if p.run_sim is True:
            #     # initialize the voltage-gated potassium object
            #     self.vgKir_object.init(self, sim, cells, p)

        #-------Funny Current------------------------------------------------------------------------------------------
        if p.vgFun_bool:

            # Initialization of logic values for voltage gated potassium channel
            self.maxDmFun = p.vgFun_max

            self.apply_vgFun = p.vgFun_apply

            self.targets_vgFun = []

            for profile in self.apply_vgFun:
                targets = self.tissue_target_inds[profile]
                self.targets_vgFun.append(targets)

            self.targets_vgFun = [item for sublist in self.targets_vgFun for item in sublist]
            self.targets_vgFun = np.asarray(self.targets_vgFun)

            # create the desired funny current channel instances:
            for ind, type_string in enumerate(p.vgFun_type):

                Fun_class = getattr(vgfun,type_string,'HCN4')

                object_name = 'vgFun_object' + str(ind)

                setattr(self,object_name,Fun_class())

                if p.run_sim is True:
                    # initialize the voltage-gated sodium object
                    init_funk = getattr(self,object_name)
                    init_funk.init(self, sim, cells, p)


            # Fun_class_ = getattr(vgfun,p.vgFun_type,'HCN4')
            # self.vgFun_object = Fun_class_()
            #
            # if p.run_sim is True:
            #     # initialize the voltage-gated sodium/potassium object
            #     self.vgFun_object.init(self, sim, cells, p)

        #-------Calcium Channels--------------------------------------------------------------------------------------

        if p.vgCa_bool and p.ions_dict['Ca'] != 0:

            # Initialization of logic values for voltage gated sodium channel
            self.maxDmCa = p.vgCa_max

            self.apply_vgCa = p.vgCa_apply

            self.targets_vgCa = []

            for profile in self.apply_vgCa:
                targets = self.tissue_target_inds[profile]
                self.targets_vgCa.append(targets)

            self.targets_vgCa = [item for sublist in self.targets_vgCa for item in sublist]
            self.targets_vgCa = np.asarray(self.targets_vgCa)

            # create the desired voltage gated calcium channel instances:
            for ind, type_string in enumerate(p.vgCa_type):

                Ca_class = getattr(vgca,type_string,'Ca_L')

                object_name = 'vgCa_object' + str(ind)

                setattr(self,object_name,Ca_class())

                if p.run_sim is True:
                    # initialize the voltage-gated sodium object
                    init_funk = getattr(self,object_name)
                    init_funk.init(self, sim, cells, p)

            # Ca_class_ = getattr(vgca, p.vgCa_type, 'Ca_L')
            # self.vgCa_object = Ca_class_()
            #
            # if p.run_sim is True:
            #     # initialize the voltage-gated sodium object
            #     self.vgCa_object.init(self, sim, cells, p)

        #--------------------------------------------------------------------------------------------------------------
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

    def _sim_events_global(self, sim, cells, p, t) -> None:
        '''
        Apply all **global scheduled interventions** (i.e., events globally
        applicable to all cells) specified by the passed user-specified
        parameters to the passed tissue simulation and cellular world for the
        passed time step.
        '''

        if p.global_options['K_env'] != 0:
            effector_Kenv = tb.pulse(
                t, self.t_on_Kenv, self.t_off_Kenv, self.t_change_Kenv)

            if p.sim_ECM: # simulate addition of potassium salt to remain charge neutral
                sim.c_env_bound[sim.iK] = (
                    self.mem_mult_Kenv*effector_Kenv*p.env_concs['K'] +
                    p.env_concs['K']
                )
                sim.c_env_bound[sim.iM] = (
                    self.mem_mult_Kenv*effector_Kenv*p.env_concs['K'] +
                    p.env_concs['M']
                )
            else:
                sim.cc_env[sim.iK][:] = (
                    self.mem_mult_Kenv*effector_Kenv*p.cK_env + p.cK_env)

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

        if p.global_options['T_change'] != 0:
            sim.T = self.multT*tb.pulse(t,self.tonT,self.toffT,self.trampT)*p.T + p.T

        if p.global_options['gj_block'] != 0:
            sim.gj_block[self.targets_gj_block] = (1.0 - tb.pulse(t,self.tonGJ,self.toffGJ,self.trampGJ))

        if p.global_options['NaKATP_block'] != 0:
            sim.NaKATP_block = (1.0 - tb.pulse(t,self.tonNK,self.toffNK,self.trampNK))

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

        if p.scheduled_options['pressure'] != 0:
            # logs.log_debug('Applying pressure event...')
            sim.P_mod[self.targets_P] = (
                self.scalar_P*self.dyna_P(t)*self.rate_P*tb.pulse(
                    t,self.t_onP, self.t_offP, self.t_changeP))

        if p.sim_ECM is True:

            if p.scheduled_options['ecmJ'] != 0:
                for i, dmat in enumerate(sim.D_env):

                    effector_ecmJ = self.mult_ecmJ*tb.pulse(
                        t,self.t_on_ecmJ,self.t_off_ecmJ,self.t_change_ecmJ)

                    sim.D_env[i][self.targets_ecmJ] = (
                        sim.D_env_base[i][self.targets_ecmJ]*(
                            1 - effector_ecmJ) + effector_ecmJ*sim.D_free[i])

                    sim.D_env_weight = sim.D_env_weight.ravel()
                    sim.D_env_weight_base = sim.D_env_weight_base.ravel()

                    sim.D_env_weight[self.targets_ecmJ] = \
                        sim.D_env_weight_base[self.targets_ecmJ]*(1-effector_ecmJ) + \
                                                          effector_ecmJ

                    sim.D_env_weight = sim.D_env_weight.reshape(cells.X.shape)
                    sim.D_env_weight_base = sim.D_env_weight_base.reshape(cells.X.shape)

                # for i, dmat in enumerate(sim.D_env):
                #     if p.env_type is True:
                #         sim.D_env_u[i] = interp.griddata((cells.xypts[:,0],cells.xypts[:,1]),dmat.ravel(),
                #             (cells.grid_obj.u_X,cells.grid_obj.u_Y),method='nearest',fill_value = sim.D_free[i])
                #         sim.D_env_v[i] = interp.griddata((cells.xypts[:,0],cells.xypts[:,1]),dmat.ravel(),
                #             (cells.grid_obj.v_X,cells.grid_obj.v_Y),method='nearest',fill_value=sim.D_free[i])
                #
                #     else:
                #         sim.D_env_u[i] = interp.griddata((cells.xypts[:,0],cells.xypts[:,1]),dmat.ravel(),
                #             (cells.grid_obj.u_X,cells.grid_obj.u_Y),method='nearest',fill_value = 0)
                #
                #         sim.D_env_v[i] = interp.griddata((cells.xypts[:,0],cells.xypts[:,1]),dmat.ravel(),
                #             (cells.grid_obj.v_X,cells.grid_obj.v_Y),method='nearest',fill_value = 0)
                #
                # sim.D_env_weight_u = sim.D_env_u[sim.iP]/sim.D_env_u[sim.iP].max()
                # sim.D_env_weight_v = sim.D_env_v[sim.iP]/sim.D_env_v[sim.iP].max()

                # if p.closed_bound is True:  # set full no slip boundary condition at exterior bounds
                #     sim.D_env_weight_u[:,0] = 0
                #     sim.D_env_weight_u[:,-1] = 0
                #     sim.D_env_weight_u[0,:] = 0
                #     sim.D_env_weight_u[-1,:] = 0
                #
                #     sim.D_env_weight_v[:,0] = 0
                #     sim.D_env_weight_v[:,-1] = 0
                #     sim.D_env_weight_v[0,:] = 0
                #     sim.D_env_weight_v[-1,:] = 0

        soc = p.scheduled_options['cuts']
        if soc is not None and not soc._is_fired and t >= soc.time:
            for cut_profile_name in soc.profile_names:
                logs.log_info(
                    'Cutting cell cluster via cut profile "%s"...',
                    cut_profile_name)
                self.removeCells(p.profiles[cut_profile_name].picker, sim, cells, p)

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

        # self.dvsign = np.sign(sim.dvm)

        if p.vgNa_bool:

            for ind, type_string in enumerate(p.vgNa_type):

                object_name = 'vgNa_object' + str(ind)

                init_funk = getattr(self, object_name)
                init_funk.run(self, sim, cells, p)


            # update the voltage-gated sodium object
            # self.vgNa_object.run(self, sim, cells, p)

        if p.vgNaP_bool:
            # update the voltage-gated sodium objects
            for ind, type_string in enumerate(p.vgNaP_type):

                object_name = 'vgNaP_object' + str(ind)

                init_funk = getattr(self, object_name)
                init_funk.run(self, sim, cells, p)


            # self.vgNaP_object.run(self, sim, cells, p)

        if p.vgK_bool:
            # update the voltage-gated potassium object
            for ind, type_string in enumerate(p.vgK_type):

                object_name = 'vgK_object' + str(ind)

                init_funk = getattr(self, object_name)
                init_funk.run(self, sim, cells, p)


            # self.vgK_object.run(self, sim, cells, p)

        if p.vgKir_bool:
            # update the voltage-gated potassium object
            for ind, type_string in enumerate(p.vgKir_type):

                object_name = 'vgKir_object' + str(ind)

                init_funk = getattr(self, object_name)
                init_funk.run(self, sim, cells, p)


            # self.vgKir_object.run(self, sim, cells, p)

        if p.vgFun_bool:
            # update the voltage-gated funny current object
            for ind, type_string in enumerate(p.vgFun_type):

                object_name = 'vgFun_object' + str(ind)

                init_funk = getattr(self, object_name)
                init_funk.run(self, sim, cells, p)


            # self.vgFun_object.run(self, sim, cells, p)

        if p.vgCa_bool and p.ions_dict['Ca'] != 0:
            # update the voltage-gated calcium object

            for ind, type_string in enumerate(p.vgCa_type):

                object_name = 'vgCa_object' + str(ind)

                init_funk = getattr(self, object_name)
                init_funk.run(self, sim, cells, p)
            # self.vgCa_object.run(self, sim, cells, p)

        if p.vg_options['K_cag'] != 0 and p.ions_dict['Ca'] != 0:

            cagPotassium(self,sim,cells,p)

        if p.vg_options['Na_stretch'] != 0:

            self.stretchChannel(sim,cells,p,t)

        if self.wound_channel_used is True:

            self.wound_channel.run(self, sim, cells, p)

    def stretchChannel(self,sim,cells,p,t):

        dd = np.sqrt(sim.d_cells_x**2 + sim.d_cells_y**2)

        # calculate strain from displacement
        eta = (dd/cells.R)

        # # create a smooth bivariate spline to interpolate deformation data from cells:
        # cellinterp_x = SmoothBivariateSpline(cells.cell_centres[:, 0], cells.cell_centres[:, 1], sim.d_cells_x, kx=1,
        #                                      ky=1)
        # cellinterp_y = SmoothBivariateSpline(cells.cell_centres[:, 0], cells.cell_centres[:, 1], sim.d_cells_y, kx=1,
        #                                      ky=1)
        #
        # # calculate deformations wrt the ecm using the smooth bivariate spline:
        # dmem_x = cellinterp_x.ev(cells.mem_mids_flat[:, 0], cells.mem_mids_flat[:, 1])
        # dmem_y = cellinterp_y.ev(cells.mem_mids_flat[:, 0], cells.mem_mids_flat[:, 1])
        #
        # # obtain normal component to membrane
        # dd = dmem_x*cells.mem_vects_flat[:,2] + dmem_y*cells.mem_vects_flat[:,3]
        #
        # # strain is the divergence of the displacement:
        # eta = np.dot(cells.M_sum_mems, dd*cells.mem_sa)/cells.cell_vol

        # self.active_NaStretch[self.targets_NaStretch] = tb.hill(sim.P_cells[cells.mem_to_cells][self.targets_NaStretch],
        #         self.NaStretch_halfmax,self.NaStretch_n)

        self.active_NaStretch[self.targets_NaStretch] = tb.hill(eta[cells.mem_to_cells][self.targets_NaStretch],
                self.NaStretch_halfmax,self.NaStretch_n)

        sim.Dm_stretch[sim.iNa] = self.maxDmNaStretch*self.active_NaStretch
        sim.Dm_stretch[sim.iK] = self.maxDmNaStretch*self.active_NaStretch

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

            profile_type = (
                profile['type'] if isinstance(profile, dict) else 'cut')

            # If this is a tissue profile...
            if profile_type == 'tissue':
                profile_picker = profile['picker']
                dmem_list = profile['diffusion constants']

                self.tissue_profile_names.append(profile_name)

                self.tissue_target_inds[profile_name] = (
                    profile_picker.get_cell_indices(cells, p, ignoreECM=False))
                self.cell_target_inds[profile_name] = (
                    profile_picker.get_cell_indices(cells, p, ignoreECM=True))

                if len(self.cell_target_inds[profile_name]):
                    # Get ECM targets.
                    if p.sim_ECM is True:
                        ecm_targs_mem = list(
                            cells.map_mem2ecm[self.tissue_target_inds[profile_name]])

                        self.env_target_inds[profile_name] = ecm_targs_mem

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

    def makeAllChanges(self, sim) -> None:
        '''
        Add together all effects to finalize changes to cell membrane
        permeabilities for the current time step.
        '''

        sim.Dm_cells = (
            sim.Dm_scheduled +
            sim.Dm_cag +
            sim.Dm_morpho +
            sim.Dm_stretch +
            sim.Dm_custom +
            sim.Dm_base
        )

        sim.P_cells = sim.P_mod + sim.P_base

    #FIXME: Replace "p.scheduled_options['cuts']" everywhere below by "self".
    def removeCells(
        self, tissue_picker, sim, cells, p, hole_tag = False) -> None:
        '''
        Permanently remove all cells selected by the passed tissue picker.

        Parameters
        ---------------------------------
        tissue_picker : TissuePickerABC
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
        # open_TJ = True

        # Subtract this bitmap's clipping mask from the global cluster mask.
        bitmap_mask = tissue_picker.get_bitmapper(cells).clipping_matrix
        cells.cluster_mask = cells.cluster_mask - bitmap_mask

        # FIXME, if deformation is too much, the following line will crash as the "target_inds_cell" is null. Rayse
        # an uman weedable eggseption.

        # Indices of all cells to be removed, ignoring extracellular spaces.
        target_inds_cell = tissue_picker.get_cell_indices(cells, p, ignoreECM=True)

        # get the corresponding flags to membrane entities
        target_inds_mem = cells.cell_to_mems[target_inds_cell]
        target_inds_mem,_,_ = tb.flatten(target_inds_mem)
        target_inds_gj,_,_ = tb.flatten(cells.cell_to_nn_full[target_inds_cell])

        if p.sim_ECM is True:
            # get environmental targets around each removed cell:
            ecm_targs_cell = list(cells.map_cell2ecm[target_inds_cell])
            ecm_targs_mem = list(cells.map_mem2ecm[target_inds_mem])

            ecm_targs = []
            for v in ecm_targs_cell:
                ecm_targs.append(v)
            for v in ecm_targs_mem:
                ecm_targs.append(v)

            # if running metabolism, make a situation where the cells burst out their ATP:
            if sim.metabo is not None:
                # get concentration of ATP in cells to be removed:
                cell_ATP = sim.metabo.core.cell_concs['ATP'][target_inds_cell]*cells.cell_vol[target_inds_cell]
                # move this entire concentration to the extracellular spaces (assumed upon cell bursting)
                sim.metabo.core.env_concs['ATP'][ecm_targs_cell] = cell_ATP/(p.cell_height*cells.delta**2)

            elif sim.molecules is not None and 'ATP' in sim.molecules.core.molecules:

                # get concentration of ATP in cells to be removed:
                cell_ATP = sim.molecules.core.cell_concs['ATP'][target_inds_cell] * cells.cell_vol[target_inds_cell]
                # move this entire concentration to the extracellular spaces (assumed upon cell bursting)
                sim.molecules.core.env_concs['ATP'][ecm_targs_cell] = cell_ATP / (p.cell_height * cells.delta ** 2)

            elif sim.grn is not None and 'ATP' in sim.grn.core.molecules:

                # get concentration of ATP in cells to be removed:
                cell_ATP = sim.grn.core.cell_concs['ATP'][target_inds_cell] * cells.cell_vol[target_inds_cell]
                # move this entire concentration to the extracellular spaces (assumed upon cell bursting)
                sim.molecules.core.env_concs['ATP'][ecm_targs_cell] = cell_ATP / (p.cell_height * cells.delta ** 2)

            # redo environmental diffusion matrices by
            # setting the environmental spaces around cut world to the free value -- if desired!:
            # if open_TJ is True:
        # save the x,y coordinates of the original boundary cell and membrane points:
        old_bflag_cellxy = np.copy(cells.cell_centres[cells.bflags_cells])
        old_bflag_memxy = np.copy(cells.mem_mids_flat[cells.bflags_mems])

        # set up the situation to make world joined to cut world have more permeable membranes:
        hurt_cells = np.zeros(len(cells.cell_i))

        # If we're creating a wound-induced channel ------------------------------------------------
        target_inds_gj_unique = np.unique(target_inds_gj)

        for i, inds in enumerate(cells.cell_to_nn_full): # for all the nn inds to a cell...
            inds_array = np.asarray(inds)
            inds_in_target = np.intersect1d(inds_array,target_inds_gj_unique)

            if len(inds_in_target):
                hurt_cells[i] = 1  # flag the cell as a "hurt" cell

        hurt_inds = (hurt_cells == 1).nonzero()
        sim.hurt_mask = np.zeros(sim.cdl)
        sim.hurt_mask[hurt_inds] = 1.0

        #----------------------------------------------------------------------------------

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
            'fluxes_gj',
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
            'Dm_custom',
            'D_gj',
            'fluxes_intra',
            'cc_at_mem',
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

                        # elif len(data) == len(cells.nn_i):
                        #     data2 = np.delete(data,target_inds_gj)

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

                        # elif len(data) == len(cells.nn_i):
                        #     for index in sorted(target_inds_gj, reverse=True):
                        #         del data[index]
                        #     data2.append(data[index])

                    super_data2.append(data2)

                if type(super_data) == np.ndarray:
                    super_data2 = np.asarray(super_data2)

                setattr(sim, name, super_data2)


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


        if p.Ca_dyn is True and sim.endo_retic is not None:

            sim.endo_retic.remove_ers(sim, target_inds_cell)


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
        cells.intra_updater(p)  # creates matrix used for finite volume integration on cell patch
        cells.cell_vols(p)  # calculate the volume of cell and its internal regions
        cells.mem_processing(p)  # calculates membrane nearest neighbours, ecm interaction, boundary tags, etc
        cells.near_neigh(p)  # Calculate the nn array for each cell
        cells.voronoiGrid(p)
        cells.calc_gj_vects(p)
        cells.environment(p)  # define features of the ecm grid
        cells.make_maskM(p)
        cells.grid_len = len(cells.xypts)

        # update sim data lengths with new world values:
        sim.mdl = len(cells.mem_mids_flat)  # mems-data-length
        sim.cdl = len(cells.cell_centres)  # cells-data-length

        if p.sim_ECM is True:  # set environnment data length
            sim.edl = len(cells.xypts)

        else:
            sim.edl = len(cells.mem_mids_flat)

        logs.log_info('Re-creating cell network Poisson solver...')
        cells.graphLaplacian(p)

        # if microtubules
        if sim.mtubes is not None:
            sim.mtubes.remove_mtubes(target_inds_mem, target_inds_cell, cells, sim, p)

        # if running voltage gated gap junctions, reinnitialize them:
        if p.v_sensitive_gj and sim.gj_funk is not None:
            sim.gj_funk.init(sim, cells, p)

        # delete data from molecules objects:
        if p.molecules_enabled and sim.molecules is not None:

            sim.molecules.core.mod_after_cut_event(target_inds_cell, target_inds_mem, sim, cells, p)

        if p.metabolism_enabled and sim.metabo is not None:

            sim.metabo.core.mod_after_cut_event(target_inds_cell, target_inds_mem, sim, cells, p, met_tag = True)

        if p.grn_enabled and sim.grn is not None:

            sim.grn.core.mod_after_cut_event(target_inds_cell, target_inds_mem, sim, cells, p)

        if p.sim_ECM is True:

            # if desire for cut away space to lack tight junctions, remove new bflags from set:
            new_bcells = np.copy(cells.bflags_cells)
            new_bmems = np.copy(cells.bflags_mems)

            searchTree = sps.KDTree(cells.cell_centres)
            original_pt_inds = list(searchTree.query(old_bflag_cellxy))[1]
            cells.bflags_cells = original_pt_inds

            searchTree_m = sps.KDTree(cells.mem_mids_flat)
            original_pt_inds_m = list(searchTree_m.query(old_bflag_memxy))[1]
            cells.bflags_mems = original_pt_inds_m

            # calculate indices to tag TJ at boundary in terms of original boundary flags
            neigh_to_bcells, _, _ = tb.flatten(cells.cell_nn[cells.bflags_cells])
            all_bound_mem_inds_o = cells.cell_to_mems[cells.bflags_cells]
            interior_bound_mem_inds_o = cells.cell_to_mems[neigh_to_bcells]
            interior_bound_mem_inds_o, _, _ = tb.flatten(interior_bound_mem_inds_o)
            all_bound_mem_inds_o, _, _ = tb.flatten(all_bound_mem_inds_o)

            cells.all_bound_mem_inds = cells.map_mem2ecm[all_bound_mem_inds_o]
            cells.interior_bound_mem_inds = cells.map_mem2ecm[interior_bound_mem_inds_o]
            cells.inds_outmem = cells.map_mem2ecm[cells.bflags_mems]
            cells.ecm_inds_bound_cell = cells.map_cell2ecm[cells.bflags_cells]


            # if hole_tag is False: # if we're not defining a hole at the beginning, reassign to new bflags

            sim.initDenv(cells,p)

            # re-assign the boundary flags to the new configuration:
            cells.bflags_cells = new_bcells
            cells.bflags_mems = new_bmems

            sim.conc_J_x = np.zeros(len(cells.xypts))
            sim.conc_J_y = np.zeros(len(cells.xypts))

        if p.fluid_flow is True or p.deformation is True:
            # make a laplacian and solver for discrete transfers on closed, irregular cell network:

            cells.deform_tools(p)

            if p.deformation is True:  # if user desires deformation:

                # create a copy of cells world, to apply deformations to for visualization purposes only:
                sim.cellso = copy.deepcopy(cells)

                if p.td_deform is True:
                    # make a laplacian and solver for discrete transfers on closed, irregular cell network
                    logs.log_info('Creating cell network Poisson solver...')
                    cells.graphLaplacian(p)

        if p.sim_eosmosis is True:

            if sim.move_pumps_channels is not None:

                sim.move_pumps_channels.remove_data(target_inds_cell)

        # WOUND CHANNEL FINALIZATION-----------------------------------------

        if p.use_wound_channel is True:

            # recalculate targets for wound channel:
            match_inds = (sim.hurt_mask == 1.0).nonzero()

            mem_match_inds = tb.flatten(cells.cell_to_mems[match_inds[0]])[0]

            self.targets_vgWound = mem_match_inds

            self.maxDmWound = p.wound_Dmax   # FIXME add to params and config!

            self.wound_channel_used = True
            self.wound_channel = w.TRP()
            self.wound_channel.init(self, sim, cells, p)

        # need to also re-do tissue profiles and GJ
        self.tissueProfiles(sim, cells, p)
        cells.redo_gj(self, p)  # redo gap junctions to isolate different tissue types
