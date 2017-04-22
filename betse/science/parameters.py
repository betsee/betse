#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

# ....................{ IMPORTS                            }....................
import numpy as np
from betse.exceptions import BetseSimConfigException, BetseSimPhaseException
from betse.lib.matplotlib import mplutil
from betse.science.config import confio #, confalias
from betse.science.config.confabc import conf_alias
from betse.science.config.event import eventcut
from betse.science.config.event import eventvoltage
from betse.science.config.export.confanim import SimConfAnimAll
from betse.science.config.export.confplot import SimConfPlotAll
from betse.science.simulate.simphase import SimPhaseKind
from betse.science.tissue import tissuecls
from betse.science.tissue.tissuepick import TissuePickerBitmap
from betse.util.io.log import logs
from betse.util.path import paths
from betse.util.type.types import type_check, NumericTypes, SequenceTypes
from collections import OrderedDict

# ....................{ CLASSES                            }....................
#FIXME: Rename the "I_overlay" attribute to "is_plot_current_overlay".
class Parameters(object):
    '''
    Storage for all user-defined parameters used in world-building,
    simulation, and plotting.

    These parameters are *deserialized* (i.e., read, loaded, and converted)
    from the user-defined YAML configuration file passed to this object on
    initialization.

    Attributes (Private)
    ----------
    _config : MappingType
        Low-level dictionary containing this entire simulation configuration,
        deserialized from this configuration's YAML-formatted file and safely
        deserializable back to the same file. This private dictionary should
        *never* be accessed directly. Note that the settings encapsulated by
        this dictionary are safely retrievable and modifiable by callers *only*
        via the public :func:`conf_alias` data descriptors in this class.

    Attributes (General: Path)
    ----------
    config_dirname : str
        Absolute path of the directory containing the source YAML configuration
        file from which this object was first deserialized. This directory
        typically also contains example resources for use by BETSE's default
        configuration file (e.g., geometry-defining bitmaps).
    config_filename : str
        Absolute path of the source YAML configuration file from which this
        object was first deserialized.

    Attributes (General: Boolean)
    ----------
    sim_ECM : bool
        ``True`` only if this simulation enables the extracellular matrix (ECM),
        commonly referred to as extracellular spaces.

    Attributes (General: Scalars)
    ----------
    cell_polarizability : NumericTypes
        Constant defining the rate of cell polarizability change in electric
        fields, typically ranging ``[0.0, 1.0e-3]``.

    Attributes (Phase: Time)
    ----------
    dt : float
        Duration in seconds of each time step for the current simulation phase.
    total_time : float
        Duration in seconds of the current simulation phase, *not* accelerated
        by the current gap junction acceleration factor.
    total_time_accelerated : float
        Duration in seconds of the current simulation phase, accelerated by the
        current gap junction acceleration factor.

    Attributes (Ions)
    ----------
    ions_dict : dict
        Dictionary mapping:
        * From each of the following names of a supported core ion:
          * ``Na``, the sodium ion Na+.
          * ``K``, the potassium ion K+.
          * ``Ca``, the calcium ion Ca2+.
          * ``Cl``, the chloride ion Cl-.
          * ``H``, the hydrogen ion H+.
          * ``M``, the M anion bicarbonate HCO3-.
          * ``P``, the anionic protein P-.
        * To either:
          * 0, if that ion is disabled by the ion profile in the current
            simulation configuration.
          * 1, if that ion is enabled by that ion profile.

    Attributes (Exports)
    ----------
    anim : SimConfAnimAll
        Subconfiguration encapsulating exported simulation animations.
    plot : SimConfPlotAll
        Subconfiguration encapsulating exported simulation plots.
    plot_cell : int
        0-based index of the cell to isolate all single-cell time plots to.
        Defaults to 0, the index assigned to the first cell guaranteed to exist.
        Note that cell indices are seed-specific and may be visualized via the
        :attr:`enumerate_cells` boolean.
    I_overlay : bool
        ``True`` only if overlaying either electric current or concentration
        flux streamlines on appropriate plots and animations.

    Attributes (Tissue)
    ----------
    clipping_bitmap_matcher : TissuePickerBitmap
        Object encapsulating the bitmap whose colored pixel area specifies the
        global geometry mask to which all tissue profile bitmaps will be
        clipped.
    closed_bound : bool
        `True` if environmental boundaries are closed (i.e., _not_ open).
    tissue_profiles : list
        List of ordered dictionaries, each describing a **tissue profile**
        (i.e., instance of the :class:`ProfileABC` class identifying cells
        to be associated with particular simulation constants and parameters).
    '''

    # ..................{ ALIASES ~ path                     }..................
    world_filename = conf_alias("['init file saving']['worldfile']", str)
    init_filename  = conf_alias("['init file saving']['file']", str)
    sim_filename   = conf_alias("['sim file saving']['file']", str)

    # ..................{ ALIASES ~ bool                     }..................
    sim_ECM = conf_alias(
        "['general options']['simulate extracellular spaces']", bool)

    # ..................{ ALIASES ~ scalar                   }..................
    cell_polarizability = conf_alias(
        "['internal parameters']['cell polarizability']", NumericTypes)

    # ..................{ INITIALIZERS                       }..................
    #FIXME: Convert all or most of the variables parsed in the __init__() method
    #below should be converted into aliases of the above form. Brainy rainbows!

    @type_check
    def __init__(self, config_filename: str) -> None:
        '''
        Parse all settings from the passed YAML-formatted simulation
        configuration file.

        Parameters
        ----------------------------
        config_filename : str
            Absolute or relative path of this file.
        '''

        # Unique absolute path of the passed file and directory containing this
        # file. Since the latter uses the former, this dirname is guaranteed to
        # be non-empty and hence *NOT* raise an exception.
        self.config_filename = paths.canonicalize(config_filename)
        self.config_dirname = paths.get_dirname(self.config_filename)

        # Dictionary loaded from this YAML file.
        self._conf = confio.read(self.config_filename)

        # Preserve backward compatibility with prior configuration formats.
        self._init_backward_compatibility()

        #---------------------------------------------------------------------------------------------------------------
        # FILE HANDLING
        #---------------------------------------------------------------------------------------------------------------

        #FIXME: Convert the following variables into typical @property-style
        #properties, as they require dynamic logic and hence *CANNOT* be
        #encapsulated via conf_alias() above. When doing so, note that it would
        #be quite nice if:
        #
        #* We cached rather than recomputed each such path using
        #  @property_cached. To do so, we'll need to redefine that decorator to
        #  return a new descriptor rather than merely a @property, which will
        #  then permit us to redefine the setter() method of that descriptor to
        #  invalidate the currently cached value on being set. Pretty sweet
        #  functionality, which we *REALLY* want elsewhere as well.
        #* We set each such path using:
        #
        #     # ...this:
        #     os.path.expanduser(paths.join(self.config_dirname, ...))
        #
        #     # ...rather than merely this:
        #     paths.join(self.config_dirname, ...)
        #
        #  Since all logic that accesses these variables in the codebase
        #  *ALWAYS* wraps these variables with os.path.expanduser(), let's just
        #  do so once for each variable.

         # Define paths for saving initialization runs, simulation runs, and results:
        self.init_path = paths.join(
            self.config_dirname, self._conf['init file saving']['directory'])  # world, inits, and sims are saved and read to/from this directory.
        self.sim_path = paths.join(
            self.config_dirname, self._conf['sim file saving']['directory']) # folder to save unique simulation and data linked to init
        self.sim_results = paths.join(
            self.config_dirname, self._conf['results file saving']['sim directory']) # folder to auto-save results (graphs, images, animations)
        self.init_results = paths.join(
            self.config_dirname, self._conf['results file saving']['init directory']) # folder to auto-save results (graphs, images, ani

        #---------------------------------------------------------------------------------------------------------------
        # INIT & SIM SETTINGS
        #---------------------------------------------------------------------------------------------------------------

        #FIXME: Redundant. The set_time_profile() method already sets the more
        #appropriately named "self.init_end", "self.sim_end", and
        #"self.total_time" attributes to these values. Remove these, please.
        self.time4init = self._conf['init time settings']['total time']      # set the time for the initialization sim [s]
        self.time4sim = self._conf['sim time settings']['total time']        # set total time for simulation [s]

        #---------------------------------------------------------------------------------------------------------------
        # GENERAL OPTIONS
        #---------------------------------------------------------------------------------------------------------------

        self.autoInit = self._conf['automatically run initialization']

        self.grid_size = int(self._conf['general options']['comp grid size'])
        self.plot_grid_size = int(self._conf['general options']['plot grid size'])

       # set ion profile to be used: 'basic', 'basic_Ca', 'animal', 'xenopus', 'scratch'
        self.ion_profile = self._conf['general options']['ion profile']

        #---------------------------------------------------------------------------------------------------------------
        # WORLD OPTIONS
        #---------------------------------------------------------------------------------------------------------------


        # Geometric constants and factors
        self.wsx = float(self._conf['world options']['world size'])  # the x-dimension of the world space
        self.wsy = self.wsx  # the y-dimension of the world space [m]
        self.rc = float(self._conf['world options']['cell radius'])  # radius of single cell
        self.cell_height = float(self._conf['world options']['cell height'])  # the height of a cell in the z-direction
        self.lattice_type = self._conf['world options']['lattice type']  # hex or rect lattice base
        self.cell_space = float(self._conf['world options']['cell spacing'])  # the true cell-cell spacing
        self.nl = float(self._conf['world options']['lattice disorder'])  # noise level for the lattice

        volmult = float(self._conf['internal parameters']['environment volume multiplier'])

        self.vol_env = volmult*self.wsx*self.wsy*self.cell_height  # environmental volume for "no ECM" simulation

        #---------------------------------------------------------------------------------------------------------------
        # TISSUE PROFILES
        #---------------------------------------------------------------------------------------------------------------

        self._init_tissue_and_cut_profiles()

        #---------------------------------------------------------------------------------------------------------------
        # TARGETED INTERVENTIONS
        #---------------------------------------------------------------------------------------------------------------

        # initialize dictionary keeping track of targeted scheduled options for the sim:
        self.scheduled_options = {}

        bool_Namem = bool(self._conf['change Na mem']['event happens'])
        bool_Kmem = bool(self._conf['change K mem']['event happens'])
        bool_Clmem = bool(self._conf['change Cl mem']['event happens'])
        bool_Camem = bool(self._conf['change Ca mem']['event happens'])
        bool_ecmj = bool(self._conf['break ecm junctions']['event happens'])
        bool_pressure = bool(self._conf['apply pressure']['event happens'])

        if bool_Namem is False:
            self.scheduled_options['Na_mem'] = 0
        elif bool_Namem is True:
            on_Namem = float(self._conf['change Na mem']['change start'])
            off_Namem = float(self._conf['change Na mem']['change finish'])
            rate_Namem = float(self._conf['change Na mem']['change rate'])
            multi_Namem = float(self._conf['change Na mem']['multiplier'])
            apply_Namem = self._conf['change Na mem']['apply to']
            function = self._conf['change Na mem']['modulator function']
            Namem = [on_Namem, off_Namem, rate_Namem, multi_Namem, apply_Namem,function]
            self.scheduled_options['Na_mem'] = Namem

        if bool_Kmem is False:
            self.scheduled_options['K_mem'] = 0
        elif bool_Kmem is True:
            on_Kmem = float(self._conf['change K mem']['change start'])
            off_Kmem = float(self._conf['change K mem']['change finish'])
            rate_Kmem = float(self._conf['change K mem']['change rate'])
            multi_Kmem = float(self._conf['change K mem']['multiplier'])
            apply_Kmem = self._conf['change K mem']['apply to']
            function = self._conf['change K mem']['modulator function']
            Kmem = [on_Kmem, off_Kmem, rate_Kmem, multi_Kmem, apply_Kmem,function]
            self.scheduled_options['K_mem'] = Kmem

        if bool_Clmem is False:
            self.scheduled_options['Cl_mem'] = 0
        elif bool_Clmem is True:
            on_Clmem = float(self._conf['change Cl mem']['change start'])
            off_Clmem = float(self._conf['change Cl mem']['change finish'])
            rate_Clmem = float(self._conf['change Cl mem']['change rate'])
            multi_Clmem = float(self._conf['change Cl mem']['multiplier'])
            apply_Clmem = self._conf['change Cl mem']['apply to']
            function = self._conf['change Cl mem']['modulator function']
            Clmem = [on_Clmem, off_Clmem, rate_Clmem, multi_Clmem, apply_Clmem, function]
            self.scheduled_options['Cl_mem'] = Clmem

        if bool_Camem is False:
            self.scheduled_options['Ca_mem'] = 0
        elif bool_Camem is True:
            on_Camem = float(self._conf['change Ca mem']['change start'])
            off_Camem = float(self._conf['change Ca mem']['change finish'])
            rate_Camem = float(self._conf['change Ca mem']['change rate'])
            multi_Camem = float(self._conf['change Ca mem']['multiplier'])
            apply_Camem = self._conf['change Ca mem']['apply to']
            function = self._conf['change Ca mem']['modulator function']
            Camem = [on_Camem, off_Camem, rate_Camem, multi_Camem, apply_Camem,function]
            self.scheduled_options['Ca_mem'] = Camem

        if bool_pressure:
            on_p = float(self._conf['apply pressure']['change start'])
            off_p = float(self._conf['apply pressure']['change finish'])
            rate_p = float(self._conf['apply pressure']['change rate'])
            multi_p = float(self._conf['apply pressure']['multiplier'])
            apply_p = self._conf['apply pressure']['apply to']
            function = self._conf['apply pressure']['modulator function']
            pressure_ops = [on_p, off_p, rate_p, multi_p, apply_p,function]

            self.scheduled_options['pressure'] = pressure_ops
        else:
            self.scheduled_options['pressure'] = 0

        #FIXME: Rename this dictionary key from "extV" to "external voltage".
        #Thus spake Sessums!

        # Parameterize the voltage event if enabled.
        self.scheduled_options['extV'] = eventvoltage.make(p=self)

        if bool_ecmj is False:
            self.scheduled_options['ecmJ'] = 0
        elif bool_ecmj is True:
            on_ecmj = float(self._conf['break ecm junctions']['change start'])
            off_ecmj = float(self._conf['break ecm junctions']['change finish'])
            rate_ecmj = float(self._conf['break ecm junctions']['change rate'])
            apply_ecmj = self._conf['break ecm junctions']['apply to']
            mult_ecmj = 1.0 - float(self._conf['break ecm junctions'].get('multiplier', 0.0))
            ecmj = [on_ecmj, off_ecmj, rate_ecmj, apply_ecmj, mult_ecmj]

            self.scheduled_options['ecmJ'] = ecmj

        # Parameterize the cutting event if enabled.
        wc = self._conf['cutting event']['wound channel']
        self.use_wound_channel = wc['use channel']
        self.wound_Dmax = float(wc['max conductivity'])
        self.wound_close_factor = float(wc['closure delay'])
        self.wound_channel_activators_list = wc.get('activators', None)
        self.wound_channel_activators_Km = wc.get('Km activators', None)
        self.wound_channel_activators_n = wc.get('n activators', None)
        self.wound_channel_inhibitors_list = wc.get('inhibitors', None)
        self.wound_channel_inhibitors_Km = wc.get('Km inhibitors', None)
        self.wound_channel_inhibitors_n = wc.get('n inhibitors', None)

        self.scheduled_options['cuts'] = eventcut.make(p=self)

        #---------------------------------------------------------------------------------------------------------------
        # GLOBAL INTERVENTIONS
        #---------------------------------------------------------------------------------------------------------------

        # initialize dictionary keeping track of global scheduled options for the sim:
        self.global_options = {}

        bool_Naenv = bool(self._conf['change Na env']['event happens'])
        bool_Kenv = bool(self._conf['change K env']['event happens'])
        bool_Clenv = bool(self._conf['change Cl env']['event happens'])
        bool_gjblock = bool(self._conf['block gap junctions']['event happens'])
        bool_temp =  bool(self._conf['change temperature']['event happens'])
        bool_NaKblock = bool(self._conf['block NaKATP pump']['event happens'])

        if bool_Kenv is False:
            self.global_options['K_env'] = 0
        elif bool_Kenv is True:
            on_Kenv = float(self._conf['change K env']['change start'])
            off_Kenv = float(self._conf['change K env']['change finish'])
            rate_Kenv = float(self._conf['change K env']['change rate'])
            multi_Kenv = float(self._conf['change K env']['multiplier'])
            kenv = [on_Kenv, off_Kenv, rate_Kenv, multi_Kenv]
            self.global_options['K_env'] = kenv

        if bool_Clenv is False:
            self.global_options['Cl_env'] = 0
        elif bool_Clenv is True:
            on_Clenv = float(self._conf['change Cl env']['change start'])
            off_Clenv = float(self._conf['change Cl env']['change finish'])
            rate_Clenv = float(self._conf['change Cl env']['change rate'])
            multi_Clenv = float(self._conf['change Cl env']['multiplier'])
            Clenv = [on_Clenv, off_Clenv, rate_Clenv, multi_Clenv]
            self.global_options['Cl_env'] = Clenv

        if bool_Naenv is False:
            self.global_options['Na_env'] = 0
        elif bool_Naenv is True:
            on_Naenv = float(self._conf['change Na env']['change start'])
            off_Naenv = float(self._conf['change Na env']['change finish'])
            rate_Naenv = float(self._conf['change Na env']['change rate'])
            multi_Naenv = float(self._conf['change Na env']['multiplier'])
            Naenv = [on_Naenv, off_Naenv, rate_Naenv, multi_Naenv]
            self.global_options['Na_env'] = Naenv

        if bool_gjblock is False:
            self.global_options['gj_block'] = 0
        elif bool_gjblock is True:
            on_gj = float(self._conf['block gap junctions']['change start'])
            off_gj = float(self._conf['block gap junctions']['change finish'])
            rate_gj = float(self._conf['block gap junctions']['change rate'])
            fraction_gj = float(self._conf['block gap junctions']['random fraction'])
            gjb = [on_gj,off_gj,rate_gj,fraction_gj]
            self.global_options['gj_block'] = gjb

        if bool_temp is False:
            self.global_options['T_change'] = 0
        elif bool_temp is True:
            on_T = float(self._conf['change temperature']['change start'])
            off_T = float(self._conf['change temperature']['change finish'])
            rate_T = float(self._conf['change temperature']['change rate'])
            multi_T = float(self._conf['change temperature']['multiplier'])
            temper = [on_T, off_T, rate_T, multi_T]
            self.global_options['T_change'] = temper

        if bool_NaKblock is False:
            self.global_options['NaKATP_block'] = 0
        elif bool_NaKblock is True:
            on_nak = float(self._conf['block NaKATP pump']['change start'])
            off_nak = float(self._conf['block NaKATP pump']['change finish'])
            rate_nak = float(self._conf['block NaKATP pump']['change rate'])
            nak = [on_nak,off_nak,rate_nak]
            self.global_options['NaKATP_block'] = nak


        #--------------------------------------------------------------------------------------------------------------
        # DYNAMIC CHANNELS
        #--------------------------------------------------------------------------------------------------------------

        # cells to effect with voltage gated channels: (choices = 'none','all','random1','random50', [1,2,3])
        # self.gated_targets = self._config['ion channel target cells']

        bool_cagK = bool(self._conf['calcium gated K+']['turn on'])
        bool_stretch = bool(self._conf['stretch gated Na+']['turn on'])

        # set specific character of gated ion channel dynamics:
        opNa = self._conf['voltage gated Na+']
        opNaP = self._conf['persistent vg Na+']
        opK = self._conf['voltage gated K+']
        opKir = self._conf['additional voltage gated K+']
        opFun = self._conf['funny current']
        opCa = self._conf['voltage gated Ca2+']
        # opcK = self._config['calcium gated K+']

        opStretch = self._conf['gated ion channel options']['stretch gated Na']

        self.flux_threshold = self._conf['gated ion channel options']['flux threshold']

        # voltage gated sodium:
        self.vgNa_bool = opNa['turn on']
        self.vgNa_type = opNa['channel type']
        self.vgNa_max = opNa['max value']
        self.vgNa_apply = opNa['apply to']

        # persistent voltage gated sodium:
        self.vgNaP_bool = opNaP['turn on']
        self.vgNaP_type = opNaP['channel type']
        self.vgNaP_max = opNaP['max value']
        self.vgNaP_apply = opNaP['apply to']

        # voltage gated potassium:
        self.vgK_bool = opK['turn on']
        self.vgK_type = opK['channel type']
        self.vgK_max = opK['max value']
        self.vgK_apply = opK['apply to']

        # inward rectifying voltage gated potassium:
        self.vgKir_bool = opKir['turn on']
        self.vgKir_type = opKir['channel type']
        self.vgKir_max = opKir['max value']
        self.vgKir_apply = opKir['apply to']

        # funny current voltage gated potassium:
        self.vgFun_bool = opFun['turn on']
        self.vgFun_type = opFun['channel type']
        self.vgFun_max = opFun['max value']
        self.vgFun_apply = opFun['apply to']

        # voltage gated Calcium:
        self.vgCa_bool = opCa['turn on']
        self.vgCa_type = opCa['channel type']
        self.vgCa_max = opCa['max value']
        self.vgCa_apply = opCa['apply to']


        # calcium gated K

        max_NaSt = float(self._conf['stretch gated Na+']['max value']) * 1.0e-9
        max_cagK = float(self._conf['calcium gated K+']['max value']) * 1.0e-9

        cagK = [max_cagK, 1.0e-3, 3.0]

        apply_cagK = self._conf['calcium gated K+']['apply to']

        cagK.append(apply_cagK)

        # stretch gated Na

        stNa = [max_NaSt,float(opStretch['hill K_half']),float(opStretch['hill n'])]

        apply_stNa = self._conf['stretch gated Na+']['apply to']

        stNa.append(apply_stNa)

        # initialize dictionary holding options for dynamic channels:
        self.vg_options = {}

        if bool_cagK is False:
            self.vg_options['K_cag'] = 0
        elif bool_cagK is True:
            self.vg_options['K_cag'] = cagK

        if bool_stretch is False:
            self.vg_options['Na_stretch'] = 0
        elif bool_stretch is True:
            self.vg_options['Na_stretch'] = stNa

        # Calcium TissueHandler: Calcium Induced Calcium Release (CICR).....................................................

        # include full calcium dynamics in the situation (i.e. endoplasmic reticulum, etc)?
        self.Ca_dyn = self._conf['Ca dynamics']['turn on']

        if self.Ca_dyn:
            self.serca_max = self._conf['Ca dynamics']['serca pump max'] # maximum rate of the ER membrane Ca++ ATPase (serca)
            self.max_er = self._conf['Ca dynamics']['max er']    # maximum conductivity of ER membrane Ca++ channel
            self.act_Km_Ca = self._conf['Ca dynamics']['act Km Ca'] # concentration at which calcium activates ER opening
            self.act_n_Ca = self._conf['Ca dynamics']['act n Ca']    # exponent for Ca activation
            self.inh_Km_Ca = self._conf['Ca dynamics']['inh Km Ca']  # concentration at which calcium inhibits ER opening
            self.inh_n_Ca = self._conf['Ca dynamics']['inh n Ca']     # exponent for Ca inhibition
            self.act_Km_IP3 = self._conf['Ca dynamics']['act Km IP3']   # concentration at which IP3 (if included in biomolecules) activates opening
            self.act_n_IP3 = self._conf['Ca dynamics']['act n IP3']     # exponent for IP3 (if included in biomolecules) activates opening

        #--------------------------------------------------------------------------------------------------------------
        #  CUSTOM GENERAL NETWORK (defined in main config file)
        #--------------------------------------------------------------------------------------------------------------

        self.network_config = self._conf.get('general network', None)

        if self.network_config is not None:
            self.molecules_enabled = self.network_config['implement network']
            self.mol_mit_enabled = self.network_config['enable mitochondria']

        else:
            self.mol_mit_enabled = False

        #---------------------------------------------------------------------------------------------------------------
        #  METABOLISM
        #---------------------------------------------------------------------------------------------------------------

        self.metabolism_enabled = self._conf['metabolism settings']['metabolism simulated']

        self.metabo_config_filename = self._conf['metabolism settings']['metabolism config']

        #---------------------------------------------------------------------------------------------------------------
        #  GENE REGULATORY NETWORKS
        #---------------------------------------------------------------------------------------------------------------

        self.grn_enabled = self._conf['gene regulatory network settings']['gene regulatory network simulated']

        self.grn_config_filename = self._conf['gene regulatory network settings']['gene regulatory network config']

        #--------------------------------------------------------------------------------------------------------------
        # VARIABLE SETTINGS
        #--------------------------------------------------------------------------------------------------------------

        self.gravity = False

        self.T = float(self._conf['variable settings']['temperature'])  # system temperature

        # electroosmotic fluid flow-----------------------------------------------------
        self.fluid_flow = self._conf['variable settings']['fluid flow']['include fluid flow']
        self.mu_water = float(self._conf['variable settings']['fluid flow']['water viscocity']) # visc water [Pa.s]

        # electrodiffusive movement pumps and channels -----------------------------------
        self.sim_eosmosis = self._conf['variable settings']['channel electroosmosis']['turn on']
        self.D_membrane = float(self._conf['variable settings']['channel electroosmosis']['membrane mobility'])
        self.z_channel = float(self._conf['variable settings']['channel electroosmosis']['channel charge'])
        self.z_pump = float(self._conf['variable settings']['channel electroosmosis']['pump charge'])

        # mechanical deformation ----------------------------------------------------------
        self.deformation = self._conf['variable settings']['deformation']['turn on']

        self.galvanotropism = float(self._conf['variable settings']['deformation']['galvanotropism'])
        self.td_deform = False # this has been disabled due to ongoing technical difficulties
        self.fixed_cluster_bound = self._conf['variable settings']['deformation']['fixed cluster boundary']
        self.youngMod = float(self._conf['variable settings']['deformation']['young modulus'])
        self.mu_tissue = float(self._conf['variable settings']['deformation']['viscous damping'])

        # osmotic and electrostatic pressures --------------------------------
        self.deform_osmo = self._conf['variable settings']['pressures']['include osmotic pressure']
        self.aquaporins = float(self._conf['variable settings']['pressures']['membrane water conductivity'])

        # calculate lame's parameters from young mod and the poisson ratio:
        self.poi = 0.49 # Poisson's ratio for the biological medium
        self.lame_mu = self.youngMod/(2*(1+self.poi))
        self.lame_lamb = (self.youngMod*self.poi)/((1+self.poi)*(1-2*self.poi))
        self.mu_membrane = 1.0 # membrane viscocity
        self.zeta = -70e-3  # zeta potential of cell membrane [V]

        # Gap junction parameters ------------------
        self.gj_surface = float(self._conf['variable settings']['gap junctions']['gap junction surface area'])
        self.gj_flux_sensitive = False
        self.gj_vthresh = float(self._conf['variable settings']['gap junctions']['gj voltage threshold'])
        self.gj_vgrad  = float(self._conf['variable settings']['gap junctions']['gj voltage window'])
        self.gj_min = float(self._conf['variable settings']['gap junctions']['gj minimum'])
        self.gj_respond_flow = False # (feature currently unsupported)
        self.v_sensitive_gj = self._conf['variable settings']['gap junctions']['voltage sensitive gj']

        # Environmental features and tight junctions ---------------------------------------------------
        self.env_type = True # for now, can't handle air boundaries
        self.cluster_open = True
        self.closed_bound = False
        self.D_tj = float(self._conf['variable settings']['tight junction scaling'])
        self.D_adh = float(self._conf['variable settings']['adherens junction scaling'])
        # tight junction relative ion movement properties:
        self.Dtj_rel = {}  # use a dictionary to hold the tj values:

        self.Dtj_rel['Na']=float(self._conf['variable settings']['tight junction relative diffusion']['Na'])
        self.Dtj_rel['K']=float(self._conf['variable settings']['tight junction relative diffusion']['K'])
        self.Dtj_rel['Cl']=float(self._conf['variable settings']['tight junction relative diffusion']['Cl'])
        self.Dtj_rel['Ca']=float(self._conf['variable settings']['tight junction relative diffusion']['Ca'])
        self.Dtj_rel['M']=float(self._conf['variable settings']['tight junction relative diffusion']['M'])
        self.Dtj_rel['P']=float(self._conf['variable settings']['tight junction relative diffusion']['P'])
        self.Dtj_rel['H']=float(self._conf['variable settings']['tight junction relative diffusion']['H'])

        # default membrane diffusion constants: easy control of cell's base resting potential
        self.Dm_Na = float(self._conf['variable settings']['default tissue properties']['Dm_Na'])     # sodium [m2/s]
        self.Dm_K = float(self._conf['variable settings']['default tissue properties']['Dm_K'])     #  potassium [m2/s]
        self.Dm_Cl = float(self._conf['variable settings']['default tissue properties']['Dm_Cl'])    # chloride [m2/s]
        self.Dm_Ca = float(self._conf['variable settings']['default tissue properties']['Dm_Ca'])   #  calcium [m2/s]
        self.Dm_H = float(self._conf['variable settings']['default tissue properties']['Dm_H'])    #  hydrogen [m2/s]
        self.Dm_M = float(self._conf['variable settings']['default tissue properties']['Dm_M'])    #  anchor ion [m2/s]
        self.Dm_P = float(self._conf['variable settings']['default tissue properties']['Dm_P'])     #  proteins [m2/s]

        # environmental (global) boundary concentrations:
        self.cbnd = self._conf['variable settings']['env boundary concentrations']

        # include noise in the simulation?
        self.channel_noise_level = float(self._conf['variable settings']['noise']['static noise level'])

        self.dynamic_noise = self._conf['variable settings']['noise']['dynamic noise']
        self.dynamic_noise_level = float(self._conf['variable settings']['noise']['dynamic noise level'])

        # Modulator functions ------------------------------------------------------------------------------------------
        self.gradient_x_properties = {}
        self.gradient_y_properties = {}
        self.gradient_r_properties = {}

        self.periodic_properties = {}
        self.f_scan_properties = {}

        self.gradient_x_properties['slope'] =float(self._conf['modulator function properties']['gradient_x']['slope'])
        self.gradient_x_properties['offset'] =float(self._conf['modulator function properties']['gradient_x']['offset'])
        self.gradient_x_properties['exponent'] = float(
                                        self._conf['modulator function properties']['gradient_x'].get('exponent', 1))

        self.gradient_y_properties['slope'] =float(self._conf['modulator function properties']['gradient_y']['slope'])
        self.gradient_y_properties['offset'] = float(self._conf['modulator function properties']['gradient_y']['offset'])
        self.gradient_y_properties['exponent'] = float(
                                    self._conf['modulator function properties']['gradient_y'].get('exponent', 1))

        self.gradient_r_properties['slope'] = float(self._conf['modulator function properties']['gradient_r']['slope'])
        self.gradient_r_properties['offset'] = float(self._conf['modulator function properties']['gradient_r']['offset'])
        self.gradient_r_properties['exponent'] = float(
            self._conf['modulator function properties']['gradient_r'].get('exponent', 1))

        self.periodic_properties['frequency'] = float(self._conf['modulator function properties']['periodic']['frequency'])
        self.periodic_properties['phase'] = float(self._conf['modulator function properties']['periodic']['phase'])

        self.f_scan_properties['f start'] = \
                                float(self._conf['modulator function properties']['f_sweep']['start frequency'])

        self.f_scan_properties['f stop'] = \
                                float(self._conf['modulator function properties']['f_sweep']['end frequency'])

        #initialize the f vect field to None as it's set depending on the sim timestep:

        self.f_scan_properties['f slope'] = None

        # ................{ EXPORTS                            }................
        ro = self._conf['results options']

        # Animation subconfiguration.
        self.anim = SimConfAnimAll(conf=self._conf)

        # CSV subconfiguration.
        self.exportData = ro['save']['data']['all']['enabled']     # export all stored data for the plot_cell to a csv text file
        self.exportData2D = ro['save']['data']['vmem']['enabled']

        # use the GHK equation to calculate alt Vmem from params?
        self.GHK_calc = self._conf['variable settings']['use Goldman calculator']

        # ................{ PLOTS                              }................
        # Plot subconfiguration.
        self.plot = SimConfPlotAll(conf=self._conf)

        #FIXME: Replace all instances of "p.turn_all_plots_off" in the codebase
        #by "not p.plot.is_after_sim_show" and remove this attribute entirely.
        self.turn_all_plots_off = not self.plot.is_after_sim_show

        #FIXME: Replace all instances of "p.autosave" in the codebase
        #by "p.plot.is_after_sim_save" and remove this attribute entirely.
        self.autosave = self.plot.is_after_sim_save  # autosave all still images to a results directory

        self.plot_cutlines = ro['plot cutlines']

        # Colormaps.
        self.default_cm = mplutil.get_colormap(ro['default colormap'])
        self.background_cm = mplutil.get_colormap(ro['background colormap'])

        # Colormap for plotting gj currents on top of default colormap.
        self.gj_cm = mplutil.get_colormap(ro['gj colormap'])

        # new options for plotting reaction network graphs:
        self.plot_network = ro.get('plot networks', False)
        self.network_cm = mplutil.get_colormap(ro.get('network colormap', 'coolwarm'))

        # Colors.
        self.vcolor = ro['vector and stream color']  # color of vector and streamlines

        # True if numbering cells in plots and animations.
        self.enumerate_cells = ro['enumerate cells']
        self.plot_cell = ro['plot cell index']             # State the cell index to use for single-cell time plots
        self.plot_networks_single_cell = ro['plot networks single cell']
        self.showCells = ro['show cells']     # True = polygon patch plots, False = trimesh
        self.I_overlay = ro['overlay currents']
        self.stream_density = ro['streamline density']
        self.IecmPlot = ro['plot total current']    # True = plot extracellular currents, false plot gj
        self.plotMask = ro['plot masked geometry']

        # Plot seed options.
        self.plot_cell_cluster = ro.get('plot cell cluster', True)
        self.plot_cell_connectivity = ro.get('plot cell connectivity diagram', True)
        self.plot_cluster_mask = ro.get('plot cluster mask', True)

        #--------------------------------------------------------------------------------------------------------------
        # INTERNAL USE ONLY
        #--------------------------------------------------------------------------------------------------------------

        iu = self._conf['internal parameters']

        self.interp_type = 'nearest'

        self.smooth_level = float(iu['gaussian smoothing'])

        self.smooth_concs = iu['smooth concentrations']

        # self.media_rho = float(iu['media resistivity'])
        # self.tissue_rho = float(iu['tissue resistivity'])

        self.substances_affect_charge = iu['substances affect Vmem']

         # default free diffusion constants (cytoplasmic)
        self.Do_Na = float(iu['Do_Na'])      # free diffusion constant sodium [m2/s]
        self.Do_K = float(iu['Do_K'])      # free diffusion constant potassium [m2/s]
        self.Do_Cl = float(iu['Do_Cl'])     # free diffusion constant chloride [m2/s]
        self.Do_Ca = float(iu['Do_Ca'])     # free diffusion constant calcium [m2/s]
        self.Do_H = float(iu['Do_H'])      # free diffusion constant hydrogen [m2/s]
        self.Do_M = float(iu['Do_M'])     # free diffusion constant mystery anchor ion [m2/s]
        self.Do_P = float(iu['Do_P'])      # free diffusion constant protein [m2/s]

        # gap junction acceleration for molecular substances:+
        # self.gj_acceleration = float(iu['time acceleration'])
        self.gj_acceleration = 1.0

        # ATP charge in the cell (for metabolism mode off)
        # FIXME add these as options to the config
        self.cATP = 1.5
        self.cADP = 0.1
        self.cPi = 0.1

        # pump parameters
        self.alpha_NaK = float(iu['alpha_NaK']) # maximum rate constant sodium-potassium ATPase per unit surface area

        # FIXME add these as options to the config
        self.KmNK_Na = 12.0   # NaKATPase enzyme ext Na half-max sat value (alpha1 = 12, alpha2 = 20, alpha3 = 60)
        self.KmNK_K = 0.2     # NaKATPase enzyme ext K half-max sat value (alpha1 = 0.2, alpha2 = 0.20, alpha3 = 0.09)
        self.KmNK_ATP = 0.5   # NaKATPase enzyme ATP half-max sat value

        self.alpha_Ca = float(iu['alpha_Ca']) # pump rate for calcium ATPase in membrane [1/mol*s] 2.0e-15

        # FIXME add these as options to the config:
        self.KmCa_Ca = 3.0e-3   # CaATPase enzyme Ca half-max sat value (1.7 - 2.8 for vascular, 0.25 for platlets)
        self.KmCa_ATP = 0.5    # CaATPase enzyme ATP half-max sat value

        # partial pressure dissolved CO2
        self.CO2 = 40.0   # [mmHg]
        self.cCO2 = 1.2

        # charge states of ions
        self.z_Na = 1
        self.z_K = 1
        self.z_Cl = -1
        self.z_Ca = 2
        self.z_H = 1
        self.z_P = -1
        self.z_M = -1

        # molar mass of ions (kg/mol):
        self.M_Na = 23e-3
        self.M_K = 39e-3
        self.M_Cl = 35e-3
        self.M_Ca = 40e-3
        self.M_H = 1e-3
        self.M_P = 500e-3
        self.M_M = 60e-3

        # fundamental constants
        self.F = 96485 # Faraday constant [J/V*mol]
        self.R = 8.314  # Gas constant [J/K*mol]
        self.eo = 8.854e-12 # permittivity of free space [F/m]
        self.kb = 1.3806e-23  # Boltzmann constant [m2 kg/ s2 K1]
        self.q = 1.602e-19    # electron charge [C]
        self.mu = 1.275e-6   # magnetic permeability [H/m or N/A2]
        self.NAv = 6.022e23     # Avagadro's Number
        self.er = 80.0          # relative dielectric constant of electrolyte
        self.eedl = 1.0         # relative dielectric constant of screening layer

        self.deltaGATP = -37000    # free energy released in ATP hydrolysis under standard phys conditions [J/mol]

        self.ac = 1.0e-6  # cell-cell separation for drawing
        self.scale_cell = 0.90          # the amount to scale cell membranes in from ecm edges (only affects drawing)
        self.cm = float(iu['membrane capacitance'])           # patch capacitance of cell membrane 0.022 [F/m2]

        self.tm = 7.5e-9           # thickness of cell membrane [m]
        self.cell_sides = 4      # minimum number of membrane domains per cell (must be >2)
        self.scale_alpha = 1.4   # the amount to scale (1/d_cell) when calculating the concave hull (boundary search)
        self.merge_cut_off = (1/50)  # the fraction of nominal cell perimeter at which nearby ecm points are merged

        self.d_cell = self.rc * 2  # diameter of single cell
        self.nx = int(self.wsx / self.d_cell)  # number of lattice sites in world x index
        self.ny = int(self.wsy / self.d_cell)  # number of lattice sites in world y index
        self.wsx = self.wsx + 5 * self.nl * self.d_cell  # readjust the world size for noise
        self.wsy = self.wsy + 5 * self.nl * self.d_cell

        self.gjl = 2*self.tm + self.cell_space     # gap junction length
        self.gj_radius = 1.0e-9              # effective radius of gap junctions connecting cells [m] (range 0 to 5.0 e-9 m)

        self.um = 1e6    # multiplication factor to convert m to um

        # self.self_cap_cell = (8 + 4.1*((self.cell_height/self.rc)**0.76))*self.eo*80*self.rc

        self.isamples = 40.0  # sampling of vector data for currents

        self.mem_const = (80.0*self.eo/self.mu_membrane)*self.zeta    # electroosmosis constant

        self.water_const = (80.0*self.eo/self.mu_water)*self.zeta    # electroosmosis constant

        self.rho = 1050 # mass density of system [kg/m3]

        self.Keqm_ph = 7.94e-4          # equilibrium constant for bicarbonate buffer

        self.vm_ph = 0.1             # rate constant for bicarbonate buffer [mol/s] 5.0e-5 originally

        self.cytoplasm_viscocity = float(iu.get('cytoplasm viscocity', 1.0))

        self.u_mtube = float(iu.get('microtubule eosmo', 0.0))

        self.D_mtube = float(iu.get('microtubule diffusion', 1.0e-12))

        # simplest ion ion_profile giving realistic results with minimal ions (Na+ & K+ focus):
        if self.ion_profile == 'basic':

            self.cNa_env = 145.0
            self.cK_env = 5.0
            self.cP_env = 10.0

            zs = [self.z_Na, self.z_K, self.z_P]

            conc_env = [self.cNa_env,self.cK_env, self.cP_env]
            self.cM_env, self.z_M_env = bal_charge(conc_env,zs)

            assert self.z_M_env == -1

            self.cNa_cell = 8.0
            self.cK_cell = 125.0
            self.cP_cell = 100.0

            conc_cell = [self.cNa_cell,self.cK_cell, self.cP_cell]

            self.cM_cell, self.z_M_cell = bal_charge(conc_cell,zs)

            assert self.z_M_cell == -1

            self.ions_dict = {'Na':1,'K':1,'Cl':0,'Ca':0,'H':0,'P':1,'M':1}

            self.cell_concs ={'Na':self.cNa_cell,'K':self.cK_cell,'P':self.cP_cell,'M':self.cM_cell}
            self.env_concs ={'Na':self.cNa_env,'K':self.cK_env,'P':self.cP_env,'M':self.cM_env}
            self.mem_perms = {'Na':self.Dm_Na,'K':self.Dm_K,'P':self.Dm_P,'M':self.Dm_M}
            self.ion_charge = {'Na':self.z_Na,'K':self.z_K,'P':self.z_P,'M':self.z_M}
            self.free_diff = {'Na':self.Do_Na,'K':self.Do_K,'P':self.Do_P,'M':self.Do_M}
            self.molar_mass = {'Na':self.M_Na,'K':self.M_K,'P':self.M_P,'M':self.M_M}
            self.ion_long_name = {'Na':'sodium','K':'potassium','P':'proteins','M':'anion'}


        elif self.ion_profile == 'basic_Ca':

            self.cNa_env = 145.0
            self.cK_env = 5.0
            self.cCa_env = 1.0
            self.cP_env = 10.0

            zs = [self.z_Na, self.z_K, self.z_Ca, self.z_P]

            conc_env = [self.cNa_env,self.cK_env, self.cCa_env, self.cP_env]
            self.cM_env, self.z_M_env = bal_charge(conc_env,zs)

            assert self.z_M_env == -1

            self.cNa_cell = 8.0
            self.cK_cell = 125.0
            self.cCa_cell = 1.0e-4
            self.cP_cell = 100.0

            conc_cell = [self.cNa_cell,self.cK_cell, self.cCa_cell, self.cP_cell]

            self.cM_cell, self.z_M_cell = bal_charge(conc_cell,zs)

            assert self.z_M_cell == -1

            self.cCa_er = 0.5
            self.cM_er = self.cCa_er

            self.ions_dict = {'Na':1,'K':1,'Cl':0,'Ca':1,'H':0,'P':1,'M':1}

            self.cell_concs ={'Na':self.cNa_cell,'K':self.cK_cell,'Ca':self.cCa_cell,'P':self.cP_cell,'M':self.cM_cell}
            self.env_concs ={'Na':self.cNa_env,'K':self.cK_env,'Ca':self.cCa_env,'P':self.cP_env,'M':self.cM_env}
            self.mem_perms = {'Na':self.Dm_Na,'K':self.Dm_K,'Ca':self.Dm_Ca, 'P':self.Dm_P,'M':self.Dm_M}
            self.ion_charge = {'Na':self.z_Na,'K':self.z_K,'Ca':self.z_Ca, 'P':self.z_P,'M':self.z_M}
            self.free_diff = {'Na':self.Do_Na,'K':self.Do_K,'Ca':self.Do_Ca, 'P':self.Do_P,'M':self.Do_M}
            self.molar_mass = {'Na':self.M_Na,'K':self.M_K,'Ca':self.M_Ca, 'P':self.M_P,'M':self.M_M}
            self.ion_long_name = {'Na':'sodium','K':'potassium','Ca':'calcium','P':'proteins','M':'anion'}

        # default environmental and cytoplasmic initial values mammalian cells
        elif self.ion_profile == 'animal':

            # initialize proton concentrations to "None" placeholders
            self.cH_cell = None
            self.cH_env = None

            self.cNa_env = 145.0
            self.cK_env = 5.0
            self.cCl_env = 105.0
            self.cCa_env = 1.0
            # self.cH_env = 3.98e-5
            self.cP_env = 10.0

            zs = [self.z_Na, self.z_K, self.z_Cl, self.z_Ca, self.z_P]

            conc_env = [self.cNa_env,self.cK_env, self.cCl_env, self.cCa_env, self.cP_env]
            self.cM_env, self.z_M_env = bal_charge(conc_env,zs)

            assert self.z_M_env == -1

            self.cNa_cell = 8.0
            self.cK_cell = 125.0
            self.cCl_cell = 20.0
            self.cCa_cell = 1.0e-4
            self.cP_cell = 100.0

            conc_cell = [self.cNa_cell,self.cK_cell, self.cCl_cell, self.cCa_cell, self.cP_cell]
            self.cM_cell, self.z_M_cell = bal_charge(conc_cell,zs)

            assert self.z_M_cell == -1

            self.cCa_er = 0.5
            self.cM_er = self.cCa_er

            self.ions_dict = {'Na':1,'K':1,'Cl':1,'Ca':1,'H':1,'P':1,'M':1}

            self.cell_concs ={'Na':self.cNa_cell,'K':self.cK_cell,'Ca':self.cCa_cell,'Cl':self.cCl_cell,
                              'H':self.cH_cell,'P':self.cP_cell,'M':self.cM_cell}

            self.env_concs ={'Na':self.cNa_env,'K':self.cK_env,'Ca':self.cCa_env,'Cl':self.cCl_env,
                             'H':self.cH_env,'P':self.cP_env,'M':self.cM_env}

            self.mem_perms = {'Na':self.Dm_Na,'K':self.Dm_K,'Ca':self.Dm_Ca,'Cl':self.Dm_Cl,
                              'H':self.Dm_H,'P':self.Dm_P,'M':self.Dm_M}

            self.ion_charge = {'Na':self.z_Na,'K':self.z_K,'Ca':self.z_Ca,'Cl':self.z_Cl,
                               'H':self.z_H,'P':self.z_P,'M':self.z_M}

            self.free_diff = {'Na':self.Do_Na,'K':self.Do_K,'Ca':self.Do_Ca,'Cl':self.Do_Cl,
                              'H':self.Do_H,'P':self.Do_P,'M':self.Do_M}

            self.molar_mass = {'Na':self.M_Na,'K':self.M_K,'Ca':self.M_Ca,'Cl':self.M_Cl,
                               'H':self.M_H,'P':self.M_P,'M':self.M_M}

            self.ion_long_name = {'Na':'sodium','K':'potassium','Ca':'calcium','Cl':'chloride',
                                  'H':'protons','P':'proteins','M':'anion'}

         # default environmental and cytoplasm values invertebrate cells
        elif self.ion_profile == 'xenopus':

            # initialize proton concentrations to "None" placeholders
            self.cH_cell = None
            self.cH_env = None

            self.cNa_env = 14.50
            self.cK_env = 0.5
            self.cCl_env = 10.50
            self.cCa_env = 0.2
            self.cP_env = 0.0

            zs = [self.z_Na, self.z_K, self.z_Cl, self.z_Ca, self.z_P]

            conc_env = [self.cNa_env,self.cK_env, self.cCl_env, self.cCa_env, self.cP_env]
            self.cM_env, self.z_M_env = bal_charge(conc_env,zs)

            assert self.z_M_env == -1

            self.cNa_cell = 8.0
            self.cK_cell = 125.0
            self.cCl_cell = 20.0
            self.cCa_cell = 1.0e-3
            self.cP_cell = 100.0

            conc_cell = [self.cNa_cell,self.cK_cell, self.cCl_cell, self.cCa_cell, self.cP_cell]
            self.cM_cell, self.z_M_cell = bal_charge(conc_cell,zs)

            assert self.z_M_cell == -1

            self.cCa_er = 0.5
            self.cM_er = self.cCa_er

            self.ions_dict = {'Na':1,'K':1,'Cl':1,'Ca':1,'H':1,'P':1,'M':1}

            self.cell_concs ={'Na':self.cNa_cell,'K':self.cK_cell,'Ca':self.cCa_cell,'Cl':self.cCl_cell,
                              'H':self.cH_cell,'P':self.cP_cell,'M':self.cM_cell}
            self.env_concs ={'Na':self.cNa_env,'K':self.cK_env,'Ca':self.cCa_env,'Cl':self.cCl_env,
                             'H':self.cH_env,'P':self.cP_env,'M':self.cM_env}
            self.mem_perms = {'Na':self.Dm_Na,'K':self.Dm_K,'Ca':self.Dm_Ca,'Cl':self.Dm_Cl,
                              'H':self.Dm_H,'P':self.Dm_P,'M':self.Dm_M}
            self.ion_charge = {'Na':self.z_Na,'K':self.z_K,'Ca':self.z_Ca,'Cl':self.z_Cl,
                               'H':self.z_H,'P':self.z_P,'M':self.z_M}
            self.free_diff = {'Na':self.Do_Na,'K':self.Do_K,'Ca':self.Do_Ca,'Cl':self.Do_Cl,
                              'H':self.Do_H,'P':self.Do_P,'M':self.Do_M}

            self.molar_mass = {'Na':self.M_Na,'K':self.M_K,'Ca':self.M_Ca,'Cl':self.M_Cl,
                               'H':self.M_H,'P':self.M_P,'M':self.M_M}

            self.ion_long_name = {'Na':'sodium','K':'potassium','Ca':'calcium','Cl':'chloride',
                                  'H':'protons','P':'proteins','M':'anion'}

        elif self.ion_profile == 'scratch':
            self.cNa_env = 145.0
            self.cK_env = 5.0
            self.cCl_env = 105.0
            self.cP_env = 10.0

            zs = [self.z_Na, self.z_K, self.z_Cl, self.z_P]

            conc_env = [self.cNa_env,self.cK_env, self.cCl_env, self.cP_env]
            self.cM_env, self.z_M_env = bal_charge(conc_env,zs)

            assert self.z_M_env == -1

            self.cNa_cell = 120.0
            self.cCl_cell = 60.0
            self.cK_cell = 30.0
            self.cP_cell = 80.0

            conc_cell = [self.cNa_cell,self.cK_cell, self.cCl_cell,self.cP_cell]

            self.cM_cell, self.z_M_cell = bal_charge(conc_cell,zs)

            assert self.z_M_cell == -1

            self.ions_dict = {'Na':1,'K':1,'Cl':1,'Ca':0,'H':0,'P':1,'M':1}

            self.cell_concs ={'Na':self.cNa_cell,'K':self.cK_cell,'Cl':self.cCl_cell,'P':self.cP_cell,'M':self.cM_cell}
            self.env_concs ={'Na':self.cNa_env,'K':self.cK_env,'Cl':self.cCl_env,'P':self.cP_env,'M':self.cM_env}
            self.mem_perms = {'Na':self.Dm_Na,'K':self.Dm_K,'Cl':self.Dm_Cl ,'P':self.Dm_P,'M':self.Dm_M}
            self.ion_charge = {'Na':self.z_Na,'K':self.z_K,'Cl':self.z_Cl,'P':self.z_P,'M':self.z_M}
            self.free_diff = {'Na':self.Do_Na,'K':self.Do_K,'Cl':self.Do_Cl,'P':self.Do_P,'M':self.Do_M}
            self.molar_mass = {'Na':self.M_Na,'K':self.M_K,'Cl':self.M_Cl,'P':self.M_P,'M':self.M_M}
            self.ion_long_name = {'Na':'sodium','K':'potassium','Cl':'chloride','P':'proteins','M':'anion'}

        # user-specified environmental and cytoplasm values (customized)
        elif self.ion_profile == 'customized':

            # initialize proton concentrations to "None" placeholders
            self.cH_cell = None
            self.cH_env = None

            cip = self._conf['general options']['customized ion profile']

            self.cNa_env = float(cip['extracellular Na+ concentration'])
            self.cK_env = float(cip['extracellular K+ concentration'])
            self.cCl_env = float(cip['extracellular Cl- concentration'])
            self.cCa_env = float(cip['extracellular Ca2+ concentration'])
            self.cM_env = float(cip['extracellular HCO3- concentration'])   # FIXME change in config to be 'anion'
            self.cP_env = float(cip['extracellular protein- concentration'])

            zs = [self.z_Na, self.z_K, self.z_Cl, self.z_Ca, self.z_P]

            conc_env = [self.cNa_env, self.cK_env, self.cCl_env, self.cCa_env, self.cP_env]
            self.cM_env, self.z_M_env = bal_charge(conc_env, zs)

            self.cNa_cell = float(cip['cytosolic Na+ concentration'])
            self.cK_cell = float(cip['cytosolic K+ concentration'])
            self.cCl_cell = float(cip['cytosolic Cl- concentration'])
            self.cCa_cell = float(cip['cytosolic Ca2+ concentration'])
            self.cM_cell = float(cip['cytosolic HCO3- concentration'])
            self.cP_cell = float(cip['cytosolic protein- concentration'])

            conc_cell = [self.cNa_cell, self.cK_cell, self.cCl_cell, self.cCa_cell, self.cP_cell]
            self.cM_cell, self.z_M_cell = bal_charge(conc_cell, zs)

            assert self.z_M_cell == -1

            self.cCa_er = float(cip['endoplasmic reticulum Ca2+'])
            self.cM_er = self.cCa_er

            self.ions_dict = {'Na':1,'K':1,'Cl':1,'Ca':1,'H':0,'P':1,'M':1}

            self.cell_concs ={'Na':self.cNa_cell,'K':self.cK_cell,'Ca':self.cCa_cell,'Cl':self.cCl_cell, 'P':self.cP_cell,'M':self.cM_cell}

            self.env_concs ={'Na':self.cNa_env,'K':self.cK_env,'Ca':self.cCa_env,'Cl':self.cCl_env, 'P':self.cP_env,'M':self.cM_env}

            self.mem_perms = {'Na':self.Dm_Na,'K':self.Dm_K,'Ca':self.Dm_Ca,'Cl':self.Dm_Cl, 'P':self.Dm_P,'M':self.Dm_M}

            self.ion_charge = {'Na':self.z_Na,'K':self.z_K,'Ca':self.z_Ca,'Cl':self.z_Cl, 'P':self.z_P,'M':self.z_M}

            self.free_diff = {'Na':self.Do_Na,'K':self.Do_K,'Ca':self.Do_Ca,'Cl':self.Do_Cl, 'P':self.Do_P,'M':self.Do_M}

            self.molar_mass = {'Na':self.M_Na,'K':self.M_K,'Ca':self.M_Ca,'Cl':self.M_Cl, 'P':self.M_P,'M':self.M_M}

            self.ion_long_name = {'Na':'sodium','K':'potassium','Ca':'calcium','Cl':'chloride',
                                  'P':'proteins','M':'anion'}

        # Else, this ion profile name is unrecognized. Raise an exception.
        else:
            raise BetseSimConfigException(
                'Ion profile name "{}" unrecognized.'.format(self.ion_profile))

    # ..................{ PRIVATE ~ initializers             }..................
    def _init_backward_compatibility(self) -> None:
        '''
        Attempt to preserve backward compatibility with prior configuration
        file formats, converting all obsolete key-value pairs of this
        configuration into their modern equivalents.
        '''

        #FIXME: Excise this hack *AFTER* refactoring the codebase to use the
        #preferable "self._config" attribute defined above. Since the
        #conf_alias() data descriptor expects that private attribute rather
        #than this public attribute, this public attribute is unhelpful now.
        self.config = self._conf

        # For convenience, localize configuration subdictionaries.
        results = self._conf['results options']

        # For backward compatibility, convert the prior into the current
        # configuration format.
        if not (
            'while solving' in results and
            'after solving' in results and
            'save' in results
        ):
            # Log a non-fatal warning.
            logs.log_warning(
                'Config file results options '
                '"while solving", "after solving", and/or "save" not found. '
                'Repairing to preserve backward compatibility. '
                'Consider upgrading to the newest config file format!',
            )

            # For convenience, localize configuration subdictionaries.
            anim_save = results['save animations']
            anim_save_frames = anim_save['frames']

            # Convert the prior into the current configuration format.
            results['while solving'] = {
                'animations': {
                    'enabled': (
                               results['plot while solving'] or
                               results['save solving plot']
                    ),
                    'show':    results['plot while solving'],
                    'save':    results['save solving plot'],
                },
            }
            results['after solving'] = {
                'plots': {
                    'enabled': (
                               results['display plots'] or
                               results['automatically save plots']
                    ),
                    'show':    results['display plots'],
                    'save':    results['automatically save plots'],
                },
                'animations': {
                    'enabled': results['create all animations'],
                    'show':    results['display plots'],
                    'save':    anim_save_frames['enabled'],
                },
            }
            results['save'] = {
                'plots': {
                    'filetype': anim_save_frames['filetype'],
                    'dpi':      anim_save_frames['dpi'],
                },
                'animations': {
                    'images': {
                        'enabled':  anim_save_frames['enabled'],
                        'filetype': anim_save_frames['filetype'],
                        'dpi':      anim_save_frames['dpi'],
                    },
                    'video': {
                        'enabled':  False,
                        'filetype': 'mkv',
                        'dpi': 300,
                        'bitrate': 1500,
                        'framerate': 5,
                        'metadata': {
                            'artist':  'BETSE',
                            'genre':   'Bioinformatics',
                            'subject': 'Bioinformatics',
                            'comment': 'Produced by BETSE.',
                        },
                        'writers': [
                            'ffmpeg', 'avconv', 'mencoder', 'imagemagick'],
                        'codecs': ['auto'],
                    },
                },
                'data': {
                    'all': {
                        'enabled': results['export data to file'],
                        'filetype': 'csv',
                    },
                    'vmem': {
                        'enabled': results['export 2D data to file'],
                        'filetype': 'csv',
                    },
                }
            }

        after_solving_anims = results['after solving']['animations']
        after_solving_plots = results['after solving']['plots']

        if 'pipeline' not in after_solving_anims:
            # Log a non-fatal warning.
            logs.log_warning(
                'Config file setting "results options" -> "after solving" -> '
                '"animations" -> "pipeline" not found. '
                'Repairing to preserve backward compatibility. '
                'Consider upgrading to the newest config file format!',
            )

            # Default the value for this dictionary key to the empty list.
            after_solving_anims['pipeline'] = []

        while_solving_anims = results['while solving']['animations']

        if 'colorbar' not in while_solving_anims:
            # Log a non-fatal warning.
            logs.log_warning(
                'Config file setting "results options" -> "while solving" -> '
                '"animations" -> "colorbar" not found. '
                'Repairing to preserve backward compatibility. '
                'Consider upgrading to the newest config file format!',
            )

            # Default the value for this dictionary key to the typical settings.
            while_solving_anims['colorbar'] = {
                'autoscale': True,
                'minimum': -70.0,
                'maximum':  10.0,
            }

        if 'single cell pipeline' not in after_solving_plots:
            # Log a non-fatal warning.
            if 'single cell' not in after_solving_plots:
                logs.log_warning(
                    'Config file setting "results options" -> "after solving" -> '
                    '"plots" -> "single cell pipeline" not found. '
                    'Repairing to preserve backward compatibility. '
                    'Consider upgrading to the newest config file format!',
                )

            # Default the value for this dictionary key to the empty list.
            after_solving_plots['single cell pipeline'] = (
                after_solving_plots['single cell']['pipeline']
                if 'single cell' in after_solving_plots else [])

        if 'cell cluster pipeline' not in after_solving_plots:
            # Log a non-fatal warning.
            if 'cell cluster' not in after_solving_plots:
                logs.log_warning(
                    'Config file setting "results options" -> "after solving" -> '
                    '"plots" -> "cell cluster pipeline" not found. '
                    'Repairing to preserve backward compatibility. '
                    'Consider upgrading to the newest config file format!',
                )

            after_solving_plots['cell cluster pipeline'] = (
                after_solving_plots['cell cluster']['pipeline']
                if 'cell cluster' in after_solving_plots else [])

        if 'plot networks single cell' not in results:
            # Log a non-fatal warning.
            logs.log_warning(
                'Config file setting "results options" -> '
                '"plot networks single cell" not found. '
                'Repairing to preserve backward compatibility. '
                'Consider upgrading to the newest config file format!',
            )

            # Default the value for this dictionary key to the empty list.
            results['plot networks single cell'] = results[
                'plot single cell graphs']


    def _init_tissue_and_cut_profiles(self) -> None:
        '''
        Parse tissue and cut profile-specific parameters from the current YAML
        configuration file.
        '''

        tpd = self._conf['tissue profile definition']

        self.default_tissue_name = (
            self._conf['variable settings']['default tissue name'])
        self.clipping_bitmap_matcher = TissuePickerBitmap(
            tpd['clipping']['bitmap']['file'], self.config_dirname)

        model_hole = tpd.get('internal pore file', None)

        if model_hole is not None:
            self.clipping_bitmap_hole = TissuePickerBitmap(
                tpd['internal pore file'], self.config_dirname)
        else:
            self.clipping_bitmap_hole = None

        # If tissue profiles are currently enabled, parse all profiles.
        self.profiles = OrderedDict()
        if tpd['profiles enabled']:
            for i, profile_config in enumerate(tpd['profiles']):
                self.profiles[profile_config['name']] = tissuecls.make(
                    p=self,
                    conf=profile_config,
                    # Convert from 0-based list indices to 1-based z order.
                    z_order=i + 1,
                )

    # ..................{ EXCEPTIONS                         }..................
    def die_unless_ecm(self) -> None:
        '''
        Raise an exception unless this configuration has enabled simulation of
        the extracellular matrix and hence environmental grid spaces.
        '''

        if not self.sim_ECM:
            raise BetseSimConfigException(
                'Extracellular spaces disabled by '
                'this simulation configuration.')

    # ..................{ SETTERS                            }..................
    #FIXME: Refactor to accept an enumeration value rather than raw string --
    #or, better yet, to accept no parameter and leverage an enumeration value
    #identifying the current phase already set as an attribute of this object.
    #Ideally, use the new "betse.science.sim.SimPhaseKind" enumeration value.
    #FIXME: Actually, this method makes little sense here; clearly, this method
    #sets phase-specific attributes. All phase-specific attributes should reside
    #directly in the "SimPhase" API rather than in this class, which is
    #supposed to be phase-agnostic and hence reusable between phases.

    @type_check
    def set_time_profile(self, phase_kind: SimPhaseKind) -> None:
        '''
        Set temporal attributes specific to the passed simulation phase type.

        These attributes include:

        * ``dt``, the duration in seconds of each time step for this phase.
        * ``total_time``, the duration in seconds of this phase.

        Parameters
        ----------
        phase_kind : SimPhaseKind
            Current simulation phase type.
        '''

        # If this simulation phase is an initialization...
        if phase_kind is SimPhaseKind.INIT:
            self.dt = float(self._conf['init time settings']['time step'])
            self.init_end = float(self._conf['init time settings']['total time'])
            self.init_tsteps = self.init_end/self.dt
            self.resample = float(self._conf['init time settings']['sampling rate'])
            self.t_resample = self.resample/self.dt

            # Duration in seconds of the current simulation phase.
            self.total_time = self.init_end
        # Else if this simulation phase is a simulation.
        elif phase_kind is SimPhaseKind.SIM:
            self.dt = float(self._conf['sim time settings']['time step'])
            self.sim_end = float(self._conf['sim time settings']['total time'])
            self.sim_tsteps = self.sim_end/self.dt
            self.resample = float(self._conf['sim time settings']['sampling rate'])
            self.t_resample = self.resample/self.dt

            # Duration in seconds of the current simulation phase.
            self.total_time = self.sim_end
        else:
            raise BetseSimPhaseException(
                'Simulation phase "{}" unsupported.'.format(phase_kind.name))

        # Duration in seconds of the current simulation phase accelerated by
        # the current gap junction acceleration factor.
        self.total_time_accelerated = self.total_time * self.gj_acceleration

# ....................{ HELPERS                            }....................
#FIXME: Shift this into a more appropriate math-oriented module. Funny sunning!
def bal_charge(concentrations: SequenceTypes, zs: SequenceTypes) -> tuple:
    '''
    Sum the concentrations of profile ions with their charge state to
    determine how much net positive charge exists, returning the concentration
    of the charge compensation anion M- needed to have zero net charge.

    Parameters
    -------------
    concentrations : SequenceType
        Array defining the concentrations of all ions in a space.
    zs : SequenceType
        Array (in complementary order to `concentrations`) of ion valence
        state.

    Returns
    ---------
    (float, float)
        2-tuple `(bal_conc, valence)`, where:
        * `bal_conc` is the concentration of anion M- to create zero net
          charge.
        * `valence` is the charge of the `bal_conc`. Ideally, this should
          _always_ be -1.
    '''

    q = 0

    for conc,z in zip(concentrations, zs):
        q = q+ conc*z

        to_zero = -q
        bal_conc = abs(to_zero)
        valance = np.sign(to_zero)

        assert bal_conc >= 0

    return bal_conc, valance
