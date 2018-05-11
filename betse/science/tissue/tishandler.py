#!/usr/bin/env python3
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

# FIXME include other channels in morphogen (dye) dynamics

# ....................{ IMPORTS                            }....................
import copy
import numpy as np
from betse.exceptions import BetseSimTissueException
from betse.science.config.confenum import CellsPickerType, SolverType
from betse.science.math import modulate as mod
from betse.science.math import toolbox as tb
from betse.science.tissue.tisprofile import CutProfile, TissueProfile
from betse.science.tissue.event.tisevecut import SimEventCut
from betse.science.tissue.picker.tispickcls import (
    TissuePickerAll, TissuePickerIndices, TissuePickerPercent)
from betse.science.tissue.picker.tispickimage import TissuePickerImage
from betse.util.io.log import logs
# from betse.util.type import types
from betse.util.type.types import type_check
from collections import OrderedDict
from random import shuffle

# ....................{ CLASSES                            }....................
#FIXME: Since this class handles both tissue profiles *AND* scheduled
#interventions, consider renaming this class to a less ambiguous name -- say, to
#"TissueEventHandler".
#FIXME: Document all instance variables defined within this class.
class TissueHandler(object):
    '''
    High-level tissue handler, managing all tissue-centric functionality
    including both tissue and cut profiles *and* scheduled interventions.

    This handler governs all:

    * Tissue profiles and objects required by these profiles, including all
      geometry-specifying bitmaps.
    * Scheduled interventions, even those *not* pertaining to tissue profiles
      (e.g., global scheduled interventions).

    Attributes (Profile: Tissue)
    ----------
    cell_target_inds : dict
        Dictionary mapping from the name of each tissue profile enabled by the
        current simulation configuration to a one-dimensional NumPy array of the
        indices of all cells in the cluster belonging to this tissue.
    tissue_default : TissueProfile
        **Default tissue profile** (i.e., profile applicable to all cells *not*
        already assigned to a non-default tissue profile). For convenience, this
        profile is accessible via both this instance variable *and* as a
        standard tissue profile (i.e., as a value of a key-value pair of the
        ordered :attr:`tissue_name_to_profile` dictionary whose key is the name
        of this profile).
    tissue_name_to_profile : OrderedDict
        Ordered dictionary mapping from the name of each tissue profile enabled
        by the current simulation configuration (in YAML list order) to the
        :class:`TissueProfile` instance encapsulating this profile.
    tissue_target_inds : dict
        Dictionary mapping from the name of each tissue profile enabled by the
        current simulation configuration to a one-dimensional Numpy array of the
        indices of all cell membranes in the cluster belonging to this tissue.

    Attributes (Profile: Cut)
    ----------
    cut_name_to_profile : OrderedDict
        Ordered dictionary mapping from the name of each cut profile enabled by
        the current simulation configuration (in YAML list order) to the
        :class:`CutProfile` instance encapsulating this profile.
    event_cut : SimEventCut
        **Cutting event** (i.e., event removing a region of the current cluster
        at some time step during the simulation phase) if enabled by the
        current simulation configuration *or* ``None`` otherwise.
    _wound_channel : {TRP, NoneType}
        Wound-induced TRP-based channel if the following conditions are all
        satisfied *or* ``None`` otherwise:
        * The current simulation configuration has enabled simulation of the
          mechanosensitive Na channel associated with wounding.
        * The cutting event has already been performed, thus generating a wound
          and activating this channel.
    '''

    # ..................{ INITIALIZERS                       }..................
    # Avoid circular import dependencies.
    @type_check
    def __init__(self, p: 'betse.science.parameters.Parameters') -> None:
        '''
        Initialize this tissue handler.

        Parameters
        ----------
        p : betse.science.parameters.Parameters
            Current simulation configuration.
        '''

        # Initialize tissue and cut profiles.
        self._init_tissues(p)
        self._init_cuts(p)

        # Initialize scheduled interventions *AFTER* tissue and cut profiles, as
        # the former requires the latter.
        self._init_events(p)


    @type_check
    def _init_tissues(self, p: 'betse.science.parameters.Parameters') -> None:
        '''
        Internally initialize all tissue profiles defined by the passed
        simulation configuration.

        Specifically, this method encapsulates each low-level YAML-formatted
        list item defining a tissue profile in this configuration with a
        high-level instance of the :class:`CellsProfileABC` class.
        '''

        # Ordered dictionaries mapping from the name of each each tissue profile
        # to the corresponding profile.
        self.tissue_name_to_profile = OrderedDict()

        # Nullify all remaining tissue-centric attributes for safety.
        self.cell_target_inds = None
        self.env_target_inds = None
        self.tissue_target_inds = None

        # If tissue profiles are disabled, silently noop.
        if not p.is_tissue_profiles:
            return

        # 1-based index of the default tissue profile.
        tissue_z_order = 1

        # Object assigning a cell cluster region to the default tissue profile.
        tissue_picker = TissuePickerAll()

        #FIXME: Set this picker for the "tissue_default" object, while still
        #setting the above tissue picker for the
        #"tissue_name_to_profile["p.tissue_default.name]" object -- implying
        #that these two objects should be different objects. Viola!
        # tissue_picker = TissuePickerImage(
        #     filename=p.tissue_default.picker_image_filename,
        #     dirname=p.conf_dirname)

        # For convenience, encapsulate the default tissue profile as a
        # high-level object accessible both as a separate instance variable
        # *AND* as a standard tissue profile (i.e., as a value of a key-value
        # pair whose key is the name of the default tissue profile).
        self.tissue_default = self.tissue_name_to_profile[
            p.tissue_default.name] = TissueProfile(
                name=p.tissue_default.name,
                z_order=tissue_z_order,
                picker=tissue_picker,

                # By definition, the cell cluster shares no gap junctions with cells
                # in the environment -- since, of course, there *ARE* no cells
                # in the environment.
                is_gj_insular=True,

                Dm_Na=p.tissue_default.Dm_Na,
                Dm_K=p.tissue_default.Dm_K,
                Dm_Cl=p.tissue_default.Dm_Cl,
                Dm_Ca=p.tissue_default.Dm_Ca,
                Dm_M=p.tissue_default.Dm_M,
                Dm_P=p.tissue_default.Dm_P,
            )

        # For each non-default tissue profile...
        for tissue_profile in p.tissue_profiles:
            # If a prior profile collides with this profile's name, this profile
            # is non-unique. In this case, raise an exception.
            if tissue_profile.name in self.tissue_name_to_profile:
                raise BetseSimTissueException(
                    'Tissue profile "{0}" non-unique '
                    '(i.e., two or more tissue profiles named "{0}").'.format(
                        tissue_profile.name))

            # 1-based index of this tissue profile.
            tissue_z_order += 1

            # Object assigning a cell cluster region to this tissue profile.
            tissue_picker = None

            # Conditionally define this object.
            if tissue_profile.picker_type is CellsPickerType.ALL:
                tissue_picker = TissuePickerAll()
            elif tissue_profile.picker_type is CellsPickerType.IMAGE:
                tissue_picker = TissuePickerImage(
                    filename=tissue_profile.picker_image_filename,
                    dirname=p.conf_dirname)
            elif tissue_profile.picker_type is CellsPickerType.INDICES:
                tissue_picker = TissuePickerIndices(
                    cells_index=tissue_profile.picker_cells_index)
            elif tissue_profile.picker_type is CellsPickerType.PERCENT:
                tissue_picker = TissuePickerPercent(
                    cells_percent=tissue_profile.picker_cells_percent)
            else:
                raise BetseSimTissueException(
                    'Tissue profile picker type "{}" unrecognized.'.format(
                        tissue_profile.picker_type))

            # Map this profile's name to a high-level tissue profile object.
            self.tissue_name_to_profile[tissue_profile.name] = TissueProfile(
                name=tissue_profile.name,
                z_order=tissue_z_order,
                picker=tissue_picker,
                is_gj_insular=tissue_profile.is_gj_insular,
                Dm_Na=tissue_profile.Dm_Na,
                Dm_K=tissue_profile.Dm_K,
                Dm_Cl=tissue_profile.Dm_Cl,
                Dm_Ca=tissue_profile.Dm_Ca,
                Dm_M=tissue_profile.Dm_M,
                Dm_P=tissue_profile.Dm_P,
            )


    @type_check
    def _init_cuts(self, p: 'betse.science.parameters.Parameters') -> None:
        '''
        Internally initialize all cut profiles defined by the passed simulation
        configuration.

        Specifically, this method encapsulates each low-level YAML-formatted
        list item defining a cut profile in this configuration with a high-level
        instance of the :class:`CellsProfileABC` class.
        '''

        # Ordered dictionaries mapping from the name of each each cut profile
        # to the corresponding profile.
        self.cut_name_to_profile = OrderedDict()

        # If cut profiles are disabled, silently noop.
        if not p.is_tissue_profiles:
            return

        # For each low-level cut profile...
        for cut_index, cut_profile in enumerate(p.cut_profiles):
            # If this profile is non-unique, raise an exception.
            if cut_profile.name in self.cut_name_to_profile:
                raise BetseSimTissueException(
                    'Cut profile "{0}" non-unique '
                    '(i.e., two or more cut profiles named "{0}").'.format(
                        cut_profile.name))

            # Map this profile's name to a high-level profile object.
            self.cut_name_to_profile[cut_profile.name] = CutProfile(
                name=cut_profile.name,
                z_order=cut_index + 1,
                picker=TissuePickerImage(
                    filename=cut_profile.picker_image_filename,
                    dirname=p.conf_dirname)
            )


    @type_check
    def _init_events(self, p: 'betse.science.parameters.Parameters') -> None:
        '''
        Initialize all scheduled interventions defined by the passed simulation
        configuration.
        '''

        # Cutting event if enabled by this configuration or "None" otherwise.
        self.event_cut = None

        # Wound-induced TRP-based channel if enabled by this configuration and a
        # cutting event generating such a wound has already been performed or
        # "None" otherwise.
        self._wound_channel = None

        #FIXME: Access only public "Parameters" attributes.
        # If this event is enabled...
        ce = p._conf['cutting event']
        if bool(ce['event happens']):
            # If cut profiles are enabled...
            if self.cut_name_to_profile:
                # List of the names of all cut profiles performed by this event.
                cut_profile_names = ce['apply to']

                # For each such name...
                for cut_profile_name in cut_profile_names:
                    # If this cut profile does *NOT* exist, raise an exception.
                    if cut_profile_name not in self.cut_name_to_profile:
                        raise BetseSimTissueException(
                            'Cut profile "{}" referenced by '
                            'cutting event not found.'.format(
                                cut_profile_name))

                # Define this event.
                self.event_cut = SimEventCut(
                    # Time step at which to cut. For simplicity, this is coerced
                    # to be the start of the simulation.
                    time_step=0.0,
                    profile_names=cut_profile_names,
                )
            # Else, log a non-fatal warning.
            else:
                logs.log_warning(
                    'Ignoring cutting event, as cut profiles are disabled.')


    #FIXME: Rename to _map_tissue_profiles_to_cells(). See the
    #SimPhase.__init__() method for further discussion.
    @type_check
    def tissueProfiles(
        self,
        sim:   'betse.science.sim.Simulator',
        cells: 'betse.science.cells.Cells',
        p:     'betse.science.parameters.Parameters',
    ) -> None:
        '''
        Define all dictionaries mapping from tissue profile names to
        one-dimensional NumPy arrays of the indices of various cellular objects
        (e.g., cells, cell membranes) in the cluster belonging to these tissues.

        Specifically, this method defines the :attr:`cell_target_inds`,
        :attr:`env_target_inds`, and :attr:`tissue_target_inds` dictionaries.
        '''

        # Log this initialization.
        logs.log_debug('Mapping tissue, cut, and boundary profiles...')

        #FIXME: These three dictionaries are all indexed by tissue profile names
        #and hence should simply be folded into the "TissueProfile" class as
        #instance variables whose values are the corresponding one-dimensional
        #Numpy lists for that specific tissue. This eliminates awkward lookups.
        self.cell_target_inds = {}
        self.env_target_inds = {}

        #FIXME: Let's rename this to something less ambiguous. Since it indexes
        #cell membranes, "mem_target_inds" might be a reasonable name? Lo-ho-ho!
        self.tissue_target_inds = {}

        # For each tissue profiles defined by this simulation configuration...
        for tissue_name, tissue_profile in self.tissue_name_to_profile.items():
            # One-dimensional Numpy arrays of the indices of all cells and cell
            # membranes comprising this tissue.
            tissue_cells_index, tissue_mems_index = (
                tissue_profile.picker.pick_cells_and_mems(cells=cells, p=p))

            # Persist these arrays for this tissue.
            self.cell_target_inds[tissue_name] = tissue_cells_index
            self.tissue_target_inds[tissue_name] = tissue_mems_index

            # If this tissue is non-empty (i.e., contains at least one cell)...
            if len(self.cell_target_inds[tissue_name]):
                # Get ECM targets.
                if p.is_ecm:
                    ecm_targs_mem = list(cells.map_mem2ecm[tissue_mems_index])
                    self.env_target_inds[tissue_name] = ecm_targs_mem

                # Set the values of Dmems and ECM diffusion based on the
                # identified target indices.
                if p.ions_dict['Na'] == 1:
                    sim.Dm_cells[sim.iNa][tissue_mems_index] = (
                        tissue_profile.Dm_Na)

                if p.ions_dict['K'] == 1:
                    sim.Dm_cells[sim.iK][tissue_mems_index] = (
                        tissue_profile.Dm_K)

                if p.ions_dict['Cl'] == 1:
                    sim.Dm_cells[sim.iCl][tissue_mems_index] = (
                        tissue_profile.Dm_Cl)

                if p.ions_dict['Ca'] == 1:
                    sim.Dm_cells[sim.iCa][tissue_mems_index] = (
                        tissue_profile.Dm_Ca)

                if p.ions_dict['M'] == 1:
                    sim.Dm_cells[sim.iM][tissue_mems_index] = (
                        tissue_profile.Dm_M)

                if p.ions_dict['P'] == 1:
                    sim.Dm_cells[sim.iP][tissue_mems_index] = (
                        tissue_profile.Dm_P)

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

    # ..................{ RUNNERS ~ apply                    }..................
    def runAllDynamics(self, sim, cells, p, t: float):
        '''
        Apply all tissue manipulations specified by the passed user-specified
        parameters to the passed tissue simulation and cellular world for the
        passed time step.
        '''

        self._sim_events_global(  sim, cells, p, t)
        self._sim_events_tissue(  sim, cells, p, t)
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

            #FIXME: Isn't the exact same "data_length" already available as the
            #"self.data_length" variable initialized by the __init__() method?
            #Leaping unicorns deny the inevitability of love!
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

        if p.scheduled_options['ecmJ'] != 0 and p.is_ecm:
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

            if p.is_ecm: # simulate addition of potassium salt to remain charge neutral
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
                    self.mem_mult_Kenv*effector_Kenv*p.conc_env_k + p.conc_env_k)

        if p.global_options['Cl_env'] != 0 and p.ions_dict['Cl'] == 1:
            effector_Clenv = tb.pulse(t,self.t_on_Clenv,self.t_off_Clenv,self.t_change_Clenv)

            if not p.is_ecm:
                sim.cc_env[sim.iCl][:] = self.mem_mult_Clenv*effector_Clenv*p.conc_env_cl + p.conc_env_cl

            else:  # simulate addition of sodium chloride to remain charge neutral
                sim.c_env_bound[sim.iCl] = self.mem_mult_Clenv*effector_Clenv*p.env_concs['Cl'] + p.env_concs['Cl']
                sim.c_env_bound[sim.iNa] = self.mem_mult_Clenv*effector_Clenv*p.env_concs['Cl'] + p.env_concs['Na']

        if p.global_options['Na_env'] != 0:
            effector_Naenv = tb.pulse(t,self.t_on_Naenv,self.t_off_Naenv,self.t_change_Naenv)

            if not p.is_ecm:
                sim.cc_env[sim.iNa][:] = self.mem_mult_Naenv*effector_Naenv*p.conc_env_na + p.conc_env_na

            else: # simulate addition of sodium salt to remain charge neutral
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
            effector_Na = (
                self.scalar_Namem*self.dyna_Namem(t)*
                tb.pulse(t,self.t_on_Namem,self.t_off_Namem,self.t_change_Namem))
            sim.Dm_scheduled[sim.iNa][self.targets_Namem] = (
                self.mem_mult_Namem*effector_Na*self.tissue_default.Dm_Na)

        if p.scheduled_options['K_mem'] != 0:
            effector_K = (
                self.scalar_Kmem*self.dyna_Kmem(t)*
                tb.pulse(t,self.t_on_Kmem,self.t_off_Kmem,self.t_change_Kmem))
            sim.Dm_scheduled[sim.iK][self.targets_Kmem] = (
                self.mem_mult_Kmem*effector_K*self.tissue_default.Dm_K)

        if p.scheduled_options['Cl_mem'] != 0 and p.ions_dict['Cl'] != 0:
            effector_Cl = (
                self.scalar_Clmem*self.dyna_Clmem(t)*
                tb.pulse(t,self.t_on_Clmem,self.t_off_Clmem,self.t_change_Clmem))
            sim.Dm_scheduled[sim.iCl][self.targets_Clmem] = (
                self.mem_mult_Clmem*effector_Cl*self.tissue_default.Dm_Cl)

        if p.scheduled_options['Ca_mem'] != 0 and p.ions_dict['Ca'] != 0:
            effector_Ca = (
                self.scalar_Camem*self.dyna_Camem(t)*
                tb.pulse(t,self.t_on_Camem,self.t_off_Camem,self.t_change_Camem))
            sim.Dm_scheduled[sim.iCa][self.targets_Camem] = (
                self.mem_mult_Camem*effector_Ca*self.tissue_default.Dm_Ca)

        if p.scheduled_options['pressure'] != 0:
            # logs.log_debug('Applying pressure event...')
            sim.P_mod[self.targets_P] = (
                self.scalar_P*self.dyna_P(t)*self.rate_P*tb.pulse(
                    t,self.t_onP, self.t_offP, self.t_changeP))

        if p.scheduled_options['ecmJ'] != 0:
            if p.is_ecm:
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

            else:

                effector_ecmJ = self.mult_ecmJ * tb.pulse(
                    t, self.t_on_ecmJ, self.t_off_ecmJ, self.t_change_ecmJ)

                sim.D_env_weight = sim.D_env_weight.ravel()
                sim.D_env_weight_base = sim.D_env_weight_base.ravel()

                sim.D_env_weight[self.targets_ecmJ] = \
                    sim.D_env_weight_base[self.targets_ecmJ] * (1 - effector_ecmJ) + \
                    effector_ecmJ

                sim.D_env_weight = sim.D_env_weight.reshape(cells.X.shape)
                sim.D_env_weight_base = sim.D_env_weight_base.reshape(cells.X.shape)

        # If the cutting event is enabled but has yet to be performed, do so.
        if self.event_cut is not None and not self.event_cut.is_fired and t > p.cut_time:
            for cut_profile_name in self.event_cut.profile_names:
                logs.log_info(
                    'Cutting cell cluster via cut profile "%s"...',
                    cut_profile_name)

                cut_profile_picker = (
                    self.cut_name_to_profile[cut_profile_name].picker)
                self.removeCells(cut_profile_picker, sim, cells, p)

            logs.log_info("Cutting event successful! Resuming simulation...")

            # Redo main data length variable for this dynamics module with
            # updated world.
            self.data_length = len(cells.mem_i)

            self.tissueProfiles(sim, cells, p)
            cells.redo_gj(self, p)
            self.runAllInit(sim, cells, p)

            # Avoid repeating this cutting event at subsequent time steps.
            self.event_cut.fire()

        # If the voltage event is enabled, adjust the voltage accordingly.
        if p.scheduled_options['extV'] is not None:
            p.scheduled_options['extV'].fire(sim, t)


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

    def makeAllChanges(self, sim) -> None:
        '''
        Add together all effects to finalize changes to cell membrane
        permeabilities for the current time step.
        '''

        sim.Dm_cells = (
            sim.Dm_scheduled +
            sim.Dm_base
        )

        sim.P_cells = sim.P_mod + sim.P_base


    @type_check
    def removeCells(
        self,
        tissue_picker,
        sim:   'betse.science.sim.Simulator',
        cells: 'betse.science.cells.Cells',
        p:     'betse.science.parameters.Parameters',
    ) -> None:
        '''
        Permanently remove all cells selected by the passed tissue picker.

        Parameters
        ---------------------------------
        tissue_picker : TissuePickerABC
            Object matching all cells to be removed.
        sim : betse.science.sim.Simulation
            Current simulation.
        cells : betse.science.cells.Cells
            Current cell cluster.
        p : betse.science.parameters.Parameters
            Current simulation configuration.
        '''

        # Redo environmental diffusion matrices by setting the environmental spaces
        # around cut world to the free value (True) or not (False)?
        # open_TJ = True

        # Subtract this bitmap's clipping mask from the global cluster mask.
        # bitmap_mask = tissue_picker.get_image_mask(cells).clipping_matrix
        # cells.cluster_mask = cells.cluster_mask - bitmap_mask

        # FIXME, if deformation is too much, the following line will crash as
        # the "target_inds_cell" is null. Rayse an uman weedable eggseption.

        # Indices of all cells and cell membranes to be removed.
        target_inds_cell, target_inds_mem = tissue_picker.pick_cells_and_mems(
            cells=cells, p=p)
            # dmem_list = tissue_profile['diffusion constants']

        # get the corresponding flags to membrane entities
        target_inds_gj,_,_ = tb.flatten(cells.cell_to_nn_full[target_inds_cell])

        if p.is_ecm:
            # get environmental targets around each removed cell:
            ecm_targs_cell = list(cells.map_cell2ecm[target_inds_cell])
            ecm_targs_mem = list(cells.map_mem2ecm[target_inds_mem])

            ecm_targs = []
            for v in ecm_targs_cell:
                ecm_targs.append(v)
            for v in ecm_targs_mem:
                ecm_targs.append(v)

            # if sim.molecules is not None and 'ATP' in sim.molecules.core.molecules:
            #
            #     # get concentration of ATP in cells to be removed:
            #     cell_ATP = sim.molecules.core.cell_concs['ATP'][target_inds_cell] * cells.cell_vol[target_inds_cell]
            #     # move this entire concentration to the extracellular spaces (assumed upon cell bursting)
            #     sim.molecules.core.env_concs['ATP'][ecm_targs_cell] = cell_ATP / (p.cell_height * cells.delta ** 2)
            #
            # elif sim.grn is not None and 'ATP' in sim.grn.core.molecules:
            #
            #     # get concentration of ATP in cells to be removed:
            #     cell_ATP = sim.grn.core.cell_concs['ATP'][target_inds_cell] * cells.cell_vol[target_inds_cell]
            #     # move this entire concentration to the extracellular spaces (assumed upon cell bursting)
            #     sim.grn.core.env_concs['ATP'][ecm_targs_cell] = cell_ATP / (p.cell_height * cells.delta ** 2)

            # redo environmental diffusion matrices by
            # setting the environmental spaces around cut world to the free value -- if desired!:
            # if open_TJ is True:
        # save the x,y coordinates of the original boundary cell and membrane points:
        # old_bflag_cellxy = np.copy(cells.cell_centres[cells.bflags_cells])
        # old_bflag_memxy = np.copy(cells.mem_mids_flat[cells.bflags_mems])

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

        if p.is_ecm:
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
                    index = None
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

        if p.is_ecm is True:  # set environnment data length
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

        if p.grn_enabled and sim.grn is not None:

            sim.grn.core.mod_after_cut_event(target_inds_cell, target_inds_mem, sim, cells, p)

        # Save target inds so they can be used outside of the core simulator (i.e. dynamically in sim-grn)
        sim.target_inds_cell_o = target_inds_cell
        sim.target_inds_mem_o = target_inds_mem

        # if hole_tag is False: # if we're not defining a hole at the beginning, reassign to new bflags
        sim.initDenv(cells,p)

        sim.conc_J_x = np.zeros(len(cells.xypts))
        sim.conc_J_y = np.zeros(len(cells.xypts))

        if p.fluid_flow or p.deformation:
            # make a laplacian and solver for discrete transfers on closed, irregular cell network:

            cells.deform_tools(p)

            if p.deformation:  # if user desires deformation:

                # create a copy of cells world, to apply deformations to for visualization purposes only:
                sim.cellso = copy.deepcopy(cells)

                if p.td_deform:
                    # make a laplacian and solver for discrete transfers on closed, irregular cell network
                    logs.log_info('Creating cell network Poisson solver...')
                    cells.graphLaplacian(p)

        # if p.sim_eosmosis is True:
        #
        #     if sim.move_pumps_channels is not None:
        #
        #         sim.move_pumps_channels.remove_data(target_inds_cell)

        # calculate targets for wound channel:
        match_inds = (sim.hurt_mask == 1.0).nonzero()

        mem_match_inds = tb.flatten(cells.cell_to_mems[match_inds[0]])[0]

        self.targets_vgWound = mem_match_inds

        if p.break_TJ: # we need to redo the TJ barrier at the wound site:

            # first identify indices of exterior and interior extracellular grid sites corresponding to the wound TJ:
            cell_wound_inds = (sim.hurt_mask == 1.0).nonzero()
            neigh_to_wound_cells, _, _ = tb.flatten(cells.cell_nn[cell_wound_inds])

            interior_wound_inds = cells.cell_to_mems[neigh_to_wound_cells]
            interior_wound_inds, _, _ = tb.flatten(interior_wound_inds)

            ecmtarg_a = cells.map_mem2ecm[self.targets_vgWound]
            ecmtarg_b = cells.map_mem2ecm[interior_wound_inds]

            ecmtarg = np.hstack((ecmtarg_a, ecmtarg_b))
            # ecmtarg = ecmtarg_a


            Dw = np.copy(sim.D_env_weight.ravel())

            # assign the wound TJ to the user's requested value:
            Dw[ecmtarg] = 1.0*p.wound_TJ

            # re-assign the environmental weight matrix:
            sim.D_env_weight = Dw.reshape(cells.X.shape)*1

            # furthermore, recalculate the individual ion diffusion maps:
            if p.is_ecm:

                for i, dmat in enumerate(sim.D_env):
                    Denv_o = np.copy(dmat)

                    Denv_o[ecmtarg] = sim.D_free[i]*p.wound_TJ*sim.Dtj_rel[i]

                    # re-assign the ecm diffusion grid filled with the environmental values at wound site:
                    sim.D_env[i] = Denv_o*1.0

            # # re-calculate specific maps of conductivity in cells and environment
            # sim.sigma_cell = np.asarray([((z ** 2) * p.q * p.F * cc * D) / (p.kb * p.T) for (z, cc, D) in
            #                               zip(sim.zs, sim.cc_cells, sim.D_free)]).mean(axis=0)
            #
            # sim.sigma_env = np.asarray(
            #     [((z ** 2) * p.q * p.F * cc * D) / (p.kb * p.T) for (z, cc, D) in
            #      zip(sim.zs, sim.cc_env, sim.D_env)]).mean(
            #     axis=0)
            #
            # if p.is_ecm:
            #
            #     # environmental conductivity matrix needs to be smoothed to assured simulation stability:
            #     sim.sigma_env = gaussian_filter(sim.sigma_env.reshape(cells.X.shape), 1).ravel()

        # need to also re-do tissue profiles and GJ
        self.tissueProfiles(sim, cells, p)
        cells.redo_gj(self, p)  # redo gap junctions to isolate different tissue types

        # If this is the fast BETSE solver, re-initialize this solver.
        if p.solver_type is SolverType.FAST:
            sim.fast_sim_init(cells, p)
