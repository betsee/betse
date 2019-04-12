#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Endoplasmic reticulum (ER) functionality.
'''

# ....................{ IMPORTS                           }....................
import copy
import numpy as np
import matplotlib.pyplot as plt
from betse.science import sim_toolbox as stb
from betse.util.path import dirs, pathnames

# ....................{ CLASSES                           }....................
class EndoRetic(object):
    '''
    Endoplasmic reticulum (ER) modeller.

    This object models all physical phenomena pertaining to the ER, including
    ER-specific pumps, channels, and calcium dynamics (e.g., calcium-induced
    calcium release controlled by inositol-triphosphate). This modeller also
    contains the facilities to initialize, define the core computations for a
    simulation loop, remove the ER during a cutting event, and save and report
    on relevant data.

    Attributes
    ----------
    Ca_er_time : list
        Two-dimensional list of all cellular endoplasmic reticulum calcium ion
        concentrations for all time steps, whose:

        #. First dimension indexes each sampled time step.
        #. Second dimension indexes each cell such that each item is the
           concentration of calcium ions in that cell's endoplasmic reticulum.
    ver_time : list
        Two-dimensional list of all cellular endoplasmic reticulum Vmems
        (transmembrane voltages) for all time steps, whose:

        #. First dimension indexes each sampled time step.
        #. Second dimension indexes each cell such that each item is the
           Vmem (transmembrane voltage) across that cell's endoplasmic
           reticulum.
    '''

    # ..................{ INITIALIZERS                      }..................
    def __init__(self, sim, cells, p):

        # init basic fields
        self.er_vol = 0.1*cells.cell_vol     # er volume
        self.er_sa = 1.0*cells.cell_sa      # er surface areas
        self.Ver = np.zeros(sim.cdl)   # initial transmembrane voltage for er
        self.Q = np.zeros(sim.cdl)     # total charge in ER
        self.cm_er = p.cm    # ER membrane capacitance

        sim.cc_er = copy.deepcopy(sim.cc_cells)    # ion concentrations
        sim.cc_er[sim.iCa][:] = 0.2                # initial concentration in the ER
        self.Dm_er = copy.deepcopy(sim.cc_cells)    # membrane permeability

        # Membrane permeability altered so all are minimal.
        for arr in self.Dm_er:
            arr[:] = 1.0e-18

        self.Dm_er_base = copy.deepcopy(self.Dm_er)  # copies of Dm for ion channel dynamics
        self.Dm_channels = copy.deepcopy(self.Dm_er)

        self.zer = copy.deepcopy(sim.cc_cells)

        for i, arr in enumerate(self.zer):
            arr[:] = sim.zs[i]

    # ..................{ GETTERS                           }..................
    def get_v(self, sim, p):

        self.Q = np.sum(self.zer*p.F*sim.cc_er, axis = 0)
        self.Ver = (1/self.cm_er)*self.Q*(self.er_vol/self.er_sa)

    # ..................{ CACHERS                           }..................
    def clear_cache(self):

        self.ver_time = []
        self.Ca_er_time = []


    def write_cache(self, sim):

        self.ver_time.append(1*self.Ver)
        self.Ca_er_time.append(1*sim.cc_er[sim.iCa][:])

    # ..................{ UPDATERS                          }..................
    def update(self, sim, cells, p):

        # if p.run_sim:

        self.channels(sim, cells, p)

        # run SERCA pump:
        f_CaATP = stb.pumpCaER(sim.cc_er[sim.iCa], sim.cc_cells[sim.iCa], self.Ver, sim.T, p)

        # update with flux
        sim.cc_cells[sim.iCa] = sim.cc_cells[sim.iCa] - f_CaATP * (self.er_sa / cells.cell_vol) * p.dt
        sim.cc_er[sim.iCa] = sim.cc_er[sim.iCa] + f_CaATP * (self.er_sa / self.er_vol) * p.dt

        for i in sim.movingIons:

            IdCM = np.ones(sim.cdl)

            f_ED = stb.electroflux(sim.cc_cells[i], sim.cc_er[i], self.Dm_er[i], p.tm*IdCM, sim.zs[i]*IdCM,
                self.Ver, sim.T, p, rho=1)

            # update with flux
            sim.cc_cells[i] = sim.cc_cells[i] - f_ED*(self.er_sa/cells.cell_vol)*p.dt
            sim.cc_er[i] = sim.cc_er[i] + f_ED*(self.er_sa/self.er_vol)*p.dt

        self.get_v(sim, p)

        # print(1e3*self.Ver.mean(), sim.cc_er[sim.iCa].mean(), sim.cc_cells[sim.iCa].mean())


    def channels(self, sim, cells, p):

        # Dm_mod_mol = self.gating_max_val * tb.hill(sim.cc_cells[sim.iCa], self.gating_Hill_K, self.gating_Hill_n)
        cCa_act = (sim.cc_cells[sim.iCa]/p.act_Km_Ca)**p.act_n_Ca
        cCa_inh = (sim.cc_cells[sim.iCa] /p.inh_Km_Ca)**p.inh_n_Ca

        Dm_mod_mol = (cCa_act/(1+cCa_act))*(1/(1+cCa_inh))


        if sim.molecules is not None or sim.metabo is not None or sim.grn is not None:

            if sim.molecules is not None and 'IP3' in sim.molecules.molecules:
                cIP3_act = (sim.molecules.cell_concs['IP3']/p.act_Km_IP3)**p.act_n_IP3

            elif sim.metabo is not None and 'IP3' in sim.metabo.core.cell_concs:
                cIP3_act = (sim.metabo.core.cell_concs['IP3'] / p.act_Km_IP3) ** p.act_n_IP3

            elif sim.grn is not None and 'IP3' in sim.grn.core.cell_concs:
                cIP3_act = (sim.grn.core.cell_concs['IP3'] / p.act_Km_IP3) ** p.act_n_IP3

            else:
                cIP3_act = np.zeros(sim.cdl)

            Dm_mod_mol = (cCa_act/(1 + cCa_act))*(1/(1 + cCa_inh))*(cIP3_act/(1+cIP3_act))

        self.Dm_channels[sim.iCa] = p.max_er*Dm_mod_mol

        self.Dm_er = self.Dm_er_base + self.Dm_channels

        # print(sim.cc_er[sim.iCa][p.visual.single_cell_index], sim.cc_cells[sim.iCa][p.visual.single_cell_index],
        # sim.molecules.IP3.c_cells[p.visual.single_cell_index])


    def remove_ers(self, sim, target_inds_cell):

        # remove cells from the mit voltage list:
        ver2 = np.delete(self.Ver, target_inds_cell)
        # reassign the new data vector to the object:
        self.Ver = ver2

        erv2 = np.delete(self.er_vol, target_inds_cell)
        self.er_vol = erv2

        erca2 = np.delete(self.er_sa, target_inds_cell)
        self.er_sa = erca2

        Q2 = np.delete(self.Q, target_inds_cell)
        self.Q = Q2

        cm2 = np.delete(self.cm_er, target_inds_cell)
        self.cm_er = cm2

        # cc_er2 = []
        #
        # for i, arr in enumerate(sim.cc_er):
        #
        #     # remove cells from the mit ion array in sim:
        #     arr2 = np.delete(arr, target_inds_cell)
        #     cc_er2.append(arr2)
        #
        # sim.cc_er = np.asarray(cc_er2)

        der1 = []
        der2 = []
        der3 = []
        zer = []

        for i, arr in enumerate(self.Dm_er):

            arr2 = np.delete(arr, target_inds_cell)
            der1.append(arr2)

        self.Dm_er = np.asarray(der1)

        for i, arr in enumerate(self.Dm_er_base):

            arr2 = np.delete(arr, target_inds_cell)
            der2.append(arr2)

        self.Dm_er_base = np.asarray(der2)

        for i, arr in enumerate(self.Dm_channels):

            arr2 = np.delete(arr, target_inds_cell)
            der3.append(arr2)

        self.Dm_channels = np.asarray(der3)

        for i, arr in enumerate(self.zer):

            arr2 = np.delete(arr, target_inds_cell)
            zer.append(arr2)

        self.zer = np.asarray(zer)
