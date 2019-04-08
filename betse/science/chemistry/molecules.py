#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Controls a gene regulatory network.

Creates and electrodiffuses a suite of customizable general gene products in
the BETSE ecosystem, where the gene products are assumed to activate and/or
inhibit the expression of other genes (and therefore the production of other
gene products) in the gene regulatory network (GRN).
'''

#FIXME: Unify the large amount of code shared in common between this and the
#"gene" submodule. See the latter for details.

# ....................{ IMPORTS                           }....................
import numpy as np
from betse.science import filehandling as fh
from betse.science.chemistry.networks import MasterOfNetworks
from betse.science.chemistry.netplot import set_net_opts
from betse.science.phase.phasecls import SimPhase
from betse.util.io.log import logs
from betse.util.path import pathnames
from betse.util.type.types import type_check

# ....................{ CLASSES                           }....................
class MasterOfMolecules(object):

    # ..................{ INITIALIZERS                      }..................
    def __init__(self, p):

        #FIXME: Extract the "GeneralNetwork.betse" basename into a new
        #configuration option, perhaps in the "init file saving" section.
        #FIXME: Replace "GeneralNetwork.betse" with "GeneralNetwork.betse.gz" to
        #compress this pickled file.

        # Define data paths for saving an initialization and simulation run:
        self.savedMoM = pathnames.join(p.init_pickle_dirname, 'GeneralNetwork.betse')


    #FIXME: Most of this method appears to have been copy-and-pasted
    #from the read_mol_config() method below. That's... not the best. Let's
    #extract the code shared in common between these two methods into a new
    #_init_mols() method internally called by these two methods.
    @type_check
    def reinitialize(self, phase: SimPhase) -> None:

        # Localize high-level phase objects for convenience.
        cells = phase.cells
        p     = phase.p
        sim   = phase.sim

        # read the config file into a dictionary:
        config_dic = p.network_config

        # Time dilation:
        self.core.time_dila = float(config_dic.get('time dilation factor', 1.0))

        # obtain specific sub-dictionaries from the config file:
        substances_config = config_dic['biomolecules']
        reactions_config = config_dic.get('reactions', None)
        transporters_config = config_dic.get('transporters', None)
        channels_config = config_dic.get('channels', None)
        modulators_config = config_dic.get('modulators', None)

        self.core.tissue_init(phase, substances_config)

        if reactions_config is not None:
            # initialize the reactions of metabolism:
            self.core.read_reactions(reactions_config, sim, cells, p)
            self.core.write_reactions()
            self.core.create_reaction_matrix()
            self.core.write_reactions_env()
            self.core.create_reaction_matrix_env()

            self.reactions = True

        else:
            self.core.create_reaction_matrix()
            self.core.create_reaction_matrix_env()
            self.reactions = False

        # initialize transporters, if defined:
        if transporters_config is not None:
            self.core.read_transporters(transporters_config, phase)
            self.core.write_transporters(sim, cells, p)

            self.transporters = True

        else:
            self.transporters = False

        # initialize channels, if desired:
        if channels_config is not None:
            self.core.read_channels(channels_config, phase)
            self.channels = True

        else:
            self.channels = False

        # initialize modulators, if desired:
        if modulators_config is not None:
            self.core.read_modulators(modulators_config, sim, cells, p)
            self.modulators = True

        else:
            self.modulators = False

    # ..................{ READERS                           }..................
    @type_check
    def read_mol_config(self, phase: SimPhase) -> None:

        # Localize high-level phase objects for convenience.
        cells = phase.cells
        p     = phase.p
        sim   = phase.sim

        # read the config file into a dictionary:
        config_dic = p.network_config

        # obtain specific sub-dictionaries from the config file:
        substances_config = config_dic['biomolecules']
        reactions_config = config_dic.get('reactions', None)
        transporters_config = config_dic.get('transporters', None)
        channels_config = config_dic.get('channels', None)
        modulators_config = config_dic.get('modulators', None)

        # initialize the substances of metabolism in a core field encapsulating
        # Master of Molecules:
        self.core = MasterOfNetworks(sim, cells, substances_config, p)

        # Time dilation:
        self.core.time_dila = float(config_dic.get('time dilation factor', 1.0))

        # read in substance properties from the config file, and initialize basic properties:
        self.core.read_substances(sim, cells, substances_config, p)
        self.core.tissue_init(phase, substances_config)

        if reactions_config is not None:
            # initialize the reactions of metabolism:
            self.core.read_reactions(reactions_config, sim, cells, p)
            self.core.write_reactions()
            self.core.create_reaction_matrix()
            self.core.write_reactions_env()
            self.core.create_reaction_matrix_env()

            self.reactions = True
        else:
            self.core.create_reaction_matrix()
            self.core.create_reaction_matrix_env()
            self.reactions = False

        # initialize transporters, if defined:
        if transporters_config is not None:
            self.core.read_transporters(transporters_config, phase)
            self.core.write_transporters(sim, cells, p)

            self.transporters = True
        else:
            self.transporters = False

        # initialize channels, if desired:
        if channels_config is not None:
            self.core.read_channels(channels_config, phase)
            self.channels = True
        else:
            self.channels = False

        # initialize modulators, if desired:
        if modulators_config is not None:
            self.core.read_modulators(modulators_config, sim, cells, p)
            self.modulators = True
        else:
            self.modulators = False

        # read in network plotting options:
        self.core.net_plot_opts = config_dic.get('network plotting', None)

        # set plotting options for the network:
        set_net_opts(self.core, self.core.net_plot_opts, p)

        # after primary initialization, check and see if optimization required:
        optim_exists = config_dic.get('optimization', None)
        if optim_exists is not None:
            opti = config_dic['optimization']['optimize network']
            self.core.opti_N = config_dic['optimization']['optimization steps']
            self.core.opti_method = config_dic['optimization']['optimization method']
            self.core.target_vmem = float(config_dic['optimization']['target Vmem'])
            self.core.opti_T = float(config_dic['optimization']['optimization T'])
            self.core.opti_step = float(config_dic['optimization']['optimization step'])

            if opti:
                logs.log_info(
                    'Analyzing the general network for optimal rates...')
                self.core.optimizer(sim, cells, p)
                self.reinitialize(phase)
