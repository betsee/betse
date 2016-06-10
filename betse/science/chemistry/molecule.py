#!/usr/bin/env python3
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

"""

Create and electrodiffuse a general molecule in the BETSE ecosystem, including functionality to
pump the molecule, use it as a gating ligand, produce and consume it, and use it an enzyme to
facilitate another reaction. The molecule is assumed to be at low concentrations and to not
have a significant effect on system voltages or currents.

"""

import copy
import os
import os.path
import time
import numpy as np
from betse.science import filehandling as fh
from betse.science import toolbox as tb
from betse.science import sim_toolbox as stb

class MasterOfMolecules(object):

    def __init__(self, sim, config_settings, p):

        """
         Initializes the MasterOfMolecules object.

        sim:                An instance of simulator
        config_settings     List of dictionaries storing key settings (p.molecules_config)
        p                   An instance of params
        """

        self.init(sim, config_settings, p)

    def init(self, sim, config_settings, p):
        """
            Initialize all of the main variables for all molecules included in the simulation.
            config settings is typically = p.molecules_config

        """
        # Initialize a list that will keep track of molecule names in the simulation
        self.molecule_names = []

        for mol_dic in config_settings:
            # get each user-defined name-filed in the dictionary:
            name = str(mol_dic['name'])

            # add the name to the name catalogue:
            self.molecule_names.append(name)

            # add a field to the MasterOfMolecules corresponding to Molecule object with that name
            setattr(self, name, Molecule())

            # now set the attributes of that Molecule object with the cornucopia of variables:

            Dm = mol_dic['Dm']  # membrane diffusion coefficient [m2/s]
            Do = mol_dic['Do']  # free diffusion constant in extra and intracellular spaces [m2/s]
            z = mol_dic['z']  # charge (oxidation state)
            c_env = mol_dic['env conc']  # initial concentration in the environment [mmol/L]
            c_cell = mol_dic['cell conc']  # initial concentration in the cytoplasm [mmol/L]
            r_production = mol_dic['production rate']  # rate at which molecule is produced in the cytoplasm
            r_decay = mol_dic['decay rate']  # rate at which molecule is consumed in the cytoplasm

            # get MasterOfMolecules.name
            obj = getattr(self, name)

            # assign all properties
            obj.Dm = Dm
            obj.Do = Do
            obj.z = z
            obj.c_envo = c_env
            obj.c_cello = c_cell
            obj.r_production = r_production
            obj.r_decay = r_decay

            # assign ion channel gating properties
            icg = mol_dic['ion channel gating']

            gating_ion_o = icg['ion channel target']  # get a target ion label to gate membrane to (or 'None')

            if gating_ion_o is not None:
                obj.use_gating_ligand = True
                self.gating_ion = sim.get_ion(gating_ion_o)

            else:
                obj.use_gating_ligand = False
                self.gating_ion = []

            obj.gating_Hill_K = icg['target Hill coefficient']
            obj.gating_Hill_n = icg['target Hill exponent']
            obj.gating_max_val = icg['peak channel opening']
            obj.gating_extracell = icg['acts extracellularly']

            # assign active pumping properties
            ap = mol_dic['active pumping']
            obj.use_pumping = ap['turn on']
            obj.pump_to_cell = ap['pump to cell']
            obj.pump_max_val = ap['maximum rate']
            obj.pump_Km = ap['pump Km']

            # assign boundary change event properties
            cab = mol_dic['change at bounds']
            obj.change_at_bounds = cab['event happens']
            obj.change_bounds_start = cab['change start']
            obj.change_bounds_end = cab['change finish']
            obj.change_bounds_rate = cab['change rate']
            obj.change_bounds_target = cab['concentration']

            # assign plotting properties
            pd = mol_dic['plotting']
            obj.make_plots = pd['plot 2D']
            obj.make_ani = pd['animate']
            obj.plot_autoscale = pd['autoscale colorbar']
            obj.plot_max = pd['max val']
            obj.plot_min = pd['min val']

            # create data structures to use with sim --------
            # initialize concentrations in cells:
            obj.c_cells = np.ones(sim.cdl)*obj.c_cello
            obj.c_mems = np.ones(sim.mdl) * obj.c_cello

            # initialize concentration in the environment:
            if p.sim_ECM is False:
                obj.c_env = np.ones(sim.mdl) * obj.c_envo
            else:
                obj.c_env = np.ones(sim.edl) * obj.c_envo

            # initialize concentration at the boundary
            obj.c_bound = obj.c_envo

    def run_loop(self, sim, cells, p):

        # get the name of the specific substance:
        for name in self.molecule_names:

            obj = getattr(self,name)

            obj.pump(sim, cells, p)
            obj.transport(sim, cells, p)
            obj.updateIntra(sim, cells, p)

    def write_data(self, sim, p):

        pass

    def report(self, sim, p):
        pass

    def export_data(self, sim, p):
        pass

    def plot(self, sim, p):

        pass

    def anim(self, sim, p):

        pass


class Molecule(object):


    def transport(self, sim, cells, p):
        """
        Transports the molecule across the membrane,
        through gap junctions, and if p.sim_ECM is true,
        through extracellular spaces and the environment.

        Parameters
        -----------
        cells               Instance of Cells
        p                   Instance of Parameters


        """

        self.c_mems, self.c_env, _, _, _, _ = stb.molecule_mover(sim, self.c_mems,
                                                                self.c_env, cells, p, z=self.z, Dm = self.Dm,
                                                                Do = self.Do, c_bound = self.c_bound,
                                                                ignoreECM = True)


    def updateC(self, flux, sim, cells, p):
        """

        General updater for a flux defined on membranes and updating concentrations in
        cells and environment.

        """
        self.c_mems, self.c_env = stb.update_Co(sim, self.c_mems, self.c_cells, flux, cells, p, ignoreECM=False)


    def updateIntra(self, sim, cells, p):

        self.c_mems, self.c_cells, _ = stb.update_intra(sim, cells, self.c_mems, self.c_cells, self.Do, self.z, p)


    def pump(self, sim, cells, p):

        """
        Defines a generic active transport pump that can be used to move
        a general molecule (such as serotonin or glutamate)
        into or out of the cell by active transport.

        Works on the basic premise of enzymatic pumps defined elsewhere:

        pump_out is True:

        cX_cell + cATP  -------> cX_env + cADP + cPi

        pump_out is False:

        cX_env + cATP  <------- cX_cell + cADP + cPi

        """

        if self.use_pumping:

            self.c_mems, self.c_env, _ = stb.molecule_pump(sim, self.c_mems, self.c_env,
                                                                 cells, p, Df=self.Do, z=self.z,
                                                                 pump_into_cell=self.pump_to_cell,
                                                                 alpha_max=self.pump_max_val, Km_X=self.pump_Km,
                                                                 Km_ATP=1.0)


    def gating(self, sim, cells, p):

        # update membrane permeability if dye targets an ion channel:
        if self.use_gating_ligand:

            if self.gating_extracell is False:

                sim.Dm_mod_dye = sim.rho_channel*self.gating_max_val*tb.hill(self.c_mems,
                                                                        self.gating_Hill_K,self.gating_Hill_n)

                sim.Dm_morpho[self.gating_ion] = sim.Dm_mod_dye


            elif self.gating_extracell is True and p.sim_ECM is True:

                sim.Dm_mod_dye = self.gating_max_val*tb.hill(self.c_env,self.gating_Hill_K,self.gating_Hill_n)

                sim.Dm_morpho[self.gating_ion] = sim.rho_channel*sim.Dm_mod_dye[cells.map_mem2ecm]

    def reaction(self, sim, cells, p):

        pass


    def remove_cells(self, dyna, sim, cells, p):

        pass

    def update_boundary(self, t, p):
        """
        Run a dynamic event in the sim, which alters concentration at the global boundary.

        t:          simulation world time
        p:          parameters instance

        """

        if self.change_at_bounds:

            effector_MorphEnv = tb.pulse(t,self.change_bounds_start,self.change_bounds_end,self.change_bounds_rate)

            if p.sim_ECM is False:
                self.c_env[:] = self.conc_MorphEnv*effector_MorphEnv + self.c_envo*(1-effector_MorphEnv)

            elif p.sim_ECM is True: # simulate addition of counter salt to maintain charge neutrality:
                self.c_bound = self.conc_MorphEnv*effector_MorphEnv + self.c_envo*(1-effector_MorphEnv)

    def plot_cells(self, t, p):

        pass

    def plot_env(self, t, p):

        pass

    def anim(self, t, p):

        pass
