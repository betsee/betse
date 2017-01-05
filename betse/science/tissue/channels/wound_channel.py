#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Wound-induced transient channel response classes.
'''

# .................... IMPORTS                            ....................
from abc import ABCMeta, abstractmethod

import numpy as np
from betse.science.tissue.channels.channels_abc import ChannelsABC
from betse.util.io.log import logs
from betse.science import toolbox as tb
from betse.science import sim_toolbox as stb
from betse.exceptions import BetseSimConfigException
# from betse.science.chemistry.molecule import get_influencers


# .................... BASE                               ....................
class WoundABC(ChannelsABC, metaclass=ABCMeta):
    '''
    Abstract base class of all wound-induced (TRP based) channel classes.

    Attributes
    ----------
    _mInf :
        Equillibrium value function for m-gates of channel
    _mTau :
        Time-constant function for m-gates of channel
    _hInf :
        Equillibrium value function for h-gates of channel
    _hTau :
        Time-constant function for h-gates of channel

    '''

    def init(self, dyna, sim, cells, p):
        '''
        Initialize wound-induced (TRP based) channel at the point of wounding.

        '''

        self.modulator = 1.0

        V = sim.vm*1000

        self._init_state(V=V, dyna=dyna, sim=sim, p=p)

        self.W_factor = 2.0   # initial concentration of the 'wound factor' (which is basically pressure)
        self.W_decay = p.wound_close_factor


    def run(self, dyna, sim, cells, p):
        '''
        Simulate wound-induced (TRP based) channel activity after wounding.

        '''

        V = sim.vm*1000

        self._calculate_state(V, dyna, sim, p)

        self.W_factor = self.W_factor - self.W_factor*self.W_decay*p.dt

        self._implement_state(V, dyna, sim, cells, p)

    def _implement_state(self, V, dyna, sim, cells, p):
        # calculate m and h channel states using RK4:
        dmWound = tb.RK4(lambda m: (self._mInf - m) / self._mTau)
        dhWound = tb.RK4(lambda h: (self._hInf - h) / self._hTau)

        dyna.m_Wound = dmWound(dyna.m_Wound, p.dt * 1e3) + dyna.m_Wound
        dyna.h_Wound = dhWound(dyna.h_Wound, p.dt * 1e3) + dyna.h_Wound

        # calculate the open-probability of the channel:
        P = (dyna.m_Wound ** self._mpower) * (dyna.h_Wound ** self._hpower)

        # get modulation coefficients by any activating/inhibiting substances:
        # # FIXME won't this crash if sim.molecules is None???
        # activator_alpha, inhibitor_alpha = get_influencers(sim, sim.molecules, p.wound_channel_activators_list,
        #                                                    p.wound_channel_activators_Km, p.wound_channel_activators_n,
        #                                                    p.wound_channel_inhibitors_list,
        #                                                    p.wound_channel_inhibitors_Km, p.wound_channel_inhibitors_n,
        #                                                    reaction_zone='mems')


        # make use of activators and inhibitors to modulate open probability:
        # P = P*activator_alpha*inhibitor_alpha

        # calculate the change of charge described for this channel, as a trans-membrane flux (+ into cell):

        # if type(P) == float:
        #     delta_Q = - (dyna.maxDmWound * P * (V - self.vrev)) * (self.W_factor / (1 + self.W_factor))
        #
        # else:
        #
        #     delta_Q = - (dyna.maxDmWound*P[dyna.targets_vgWound]*(V - self.vrev))*(self.W_factor/(1 + self.W_factor))


        # obtain concentration of ion inside and out of the cell, as well as its charge z:
        c_mem_Na = sim.cc_cells[sim.iNa][cells.mem_to_cells]
        c_mem_K = sim.cc_cells[sim.iK][cells.mem_to_cells]

        if p.sim_ECM is True:
            c_env_Na = sim.cc_env[sim.iNa][cells.map_mem2ecm]
            c_env_K = sim.cc_env[sim.iK][cells.map_mem2ecm]

        else:
            c_env_Na = sim.cc_env[sim.iNa]
            c_env_K = sim.cc_env[sim.iK]

        IdM = np.ones(sim.mdl)

        z_Na = sim.zs[sim.iNa] * IdM
        z_K = sim.zs[sim.iK] * IdM

        # membrane diffusion constant of the channel:
        Dchan = dyna.maxDmWound*P*1.0e-9*(self.W_factor/(1 + self.W_factor))

        # calculate specific ion flux contribution for this channel:
        delta_Q_Na = stb.electroflux(c_env_Na, c_mem_Na, Dchan, p.tm * IdM, z_Na, sim.vm, sim.T, p, rho=sim.rho_channel)
        delta_Q_K = stb.electroflux(c_env_K, c_mem_K, Dchan, p.tm * IdM, z_K, sim.vm, sim.T, p, rho=sim.rho_channel)

        self.clip_flux(delta_Q_Na, threshold=p.flux_threshold)
        self.clip_flux(delta_Q_K, threshold=p.flux_threshold)

        self.update_charge(sim.iNa, delta_Q_Na, dyna.targets_vgWound, sim, cells, p)
        self.update_charge(sim.iK, delta_Q_K, dyna.targets_vgWound, sim, cells, p)

        # if p.ions_dict['Ca'] == 1.0:
        #     self.update_charge(sim.iCa, 0.1*delta_Q, dyna.targets_vgWound, sim, cells, p)


    @abstractmethod
    def _init_state(self, V, dyna, sim, p):
        '''
        Initializes values of the m and h gates of the channel.
        '''
        pass


    @abstractmethod
    def _calculate_state(self, V, dyna, sim, p):
        '''
        Calculates time-dependent values of the m and h gates of the channel.
        '''
        pass

# ....................{ SUBCLASS                           }....................
class TRP(WoundABC):
    '''
    TRP model default.

    '''

    def _init_state(self, V, dyna, sim, p):
        """

        Run initialization calculation for m and h gates of the channel at starting Vmem value.

        """

        logs.log_info('You are using the wound-induced channel: TRP')

        self.v_corr = 0

        # initialize values of the m and h gates of the channel on m_inf and h_inf:
        dyna.m_Wound = np.ones(sim.mdl)
        dyna.h_Wound = np.ones(sim.mdl)

        # define the power of m and h gates used in the final channel state equation:
        self._mpower = 0
        self._hpower = 0


    def _calculate_state(self, V, dyna, sim, p):
        """

        Update the state of m and h gates of the channel given their present value and present
        simulation Vmem.

        """

        self.vrev = 0  # reversal voltage used in model [mV]

        self._mInf = 1.0
        self._mTau = 1.0
        self._hInf = 1.0
        self._hTau = 1.0

def get_influencers(sim, sim_metabo, a_list, Km_a_list, n_a_list, i_list, Km_i_list,
                            n_i_list, reaction_zone='cell'):


    """
    Get coefficients representing the net effect of all activators and inhibitors on a particular reaction.

    Parameters
    ------------
    sim                 Instance of BETSE simulator
    sim_metabo:         Instance of MasterOfMetabolism
    a_list:             activator names list
    Km_a_list:          activator half-max constants
    n_a_list:           activator Hill exponents
    i_list:             inhibitor names list
    Km_i_list:          inhibitor half-max constants
    n_i_list:           inhibitor Hill exponents
    reaction_zone:      Reaction occurring in 'cell' or 'mitochondria'

    Returns
    ------------
    activator_alpha         Coefficient of net effect of activators
    inhibitor_alpha         Coefficient of net effect of inhibitors
    """

    if reaction_zone == 'cell':
        type_self = 'c_cells'
        type_sim = 'cc_cells'

    elif reaction_zone == 'mems':
        type_self = 'c_mems'
        type_sim = 'cc_mems'

    elif reaction_zone == 'mitochondria':

        type_self = 'c_mit'
        type_sim = 'cc_mit'

    elif reaction_zone == 'env':

        type_self = 'c_env'
        type_sim = 'cc_env'

    # initialize a blank list
    activator_terms = []

    if a_list is not None and a_list != 'None' and len(a_list) > 0:  # if user specified activators for growth/decay

        # get reaction zone for data type:

        # get the activator concentration for the substance, and
        # create a term based on Hill form:
        for i, activator_name in enumerate(a_list):

            label = 'i' + activator_name
            ion_check = getattr(sim, label, None)

            if ion_check is None:

                try:
                    obj_activator = getattr(sim_metabo, activator_name)
                    c_act = getattr(obj_activator, type_self)

                except KeyError:

                    raise BetseSimConfigException('Name of reaction activator is not a defined chemical, '
                                                   'or is not an ion currently included in the ion profile '
                                                   'being used.'
                                                   'Please check biomolecule definitions and ion profile'
                                                   'settings of your config(s) and try again.')

            else:
                # define the reactant as the ion concentration from the cell concentrations object in sim:
                sim_conco = getattr(sim, type_sim)
                c_act = sim_conco[ion_check]

            Km_act = Km_a_list[i]
            n_act = n_a_list[i]

            cs = (c_act / Km_act) ** n_act

            act_term = cs / (1 + cs)

            activator_terms.append(act_term)

        activator_terms = np.asarray(activator_terms)

        # calculate the net effect of all activator terms:
        activator_alpha = np.prod(activator_terms, axis=0)


    else:

        activator_alpha = 1

    # initialize a blank list
    inhibitor_terms = []

    if i_list is not None and i_list != 'None' and len(i_list) > 0:  # if user specified inhibitors for growth/decay

        # get the inhibitor concentration for the substance, and
        # create a term based on Hill form:
        for j, inhibitor_name in enumerate(i_list):

            label = 'i' + inhibitor_name
            ion_check = getattr(sim, label, None)

            if ion_check is None:

                try:
                    obj_inhibitor = getattr(sim_metabo, inhibitor_name)
                    c_inh = getattr(obj_inhibitor, type_self)

                except KeyError:

                    raise BetseSimConfigException('Name of substance is not a defined chemical, '
                                                   'or is not an ion currently included in the ion profile '
                                                   'being used.'
                                                   'Please check biomolecule definitions and ion profile'
                                                   'settings of your config(s) and try again.')

            else:
                # define the reactant as the ion concentration from the cell concentrations object in sim:
                sim_conco = getattr(sim, type_sim)
                c_inh = sim_conco[ion_check]

            # print(i_list, Km_i_list, n_i_list)

            Km_inh = Km_i_list[j]
            n_inh = n_i_list[j]

            cs = (c_inh / Km_inh) ** n_inh

            inh_term = 1 / (1 + cs)

            inhibitor_terms.append(inh_term)

        inhibitor_terms = np.asarray(inhibitor_terms)

        # calculate the net effect of all activator terms:
        inhibitor_alpha = np.prod(inhibitor_terms, axis=0)

    else:
        inhibitor_alpha = 1

    return activator_alpha, inhibitor_alpha