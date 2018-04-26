#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Facilities guaranteeing backward compatibility with prior file formats for
simulation configurations.
'''

# ....................{ IMPORTS                            }....................
from betse.science.parameters import Parameters
from betse.util.io.log import logs
from betse.util.type import iterables
from betse.util.type.types import type_check, MappingType

# ....................{ UPGRADERS                          }....................
#FIXME: If the current third-party YAML dependency is "ruamel.yaml" rather than
#PyYAML, improve this function to preserve all in-memory changes back to disk --
#albeit, for safety, presumably in a separate file preserving the existing
#filetype and basename prefix of the original file (e.g., writing a new
#"sim_conf-auto_upgraded.yaml" file given an input "sim_conf.yaml" file). If
#this file already exists, silently overwrite it; there is absolutely no chance
#and hence concern of a user intentionally embedding the substring
#"-auto_upgraded" in a simulation configuration filename.
#
#Actually, as a safety check, ensure that the current input filename does *NOT*
#already contain this substring. Trivially accomplished, thankfully.

@type_check
def upgrade_sim_conf(p: Parameters) -> None:
    '''
    Upgrade the in-memory contents of the passed simulation configuration to
    reflect the newest structure of these contents expected by the current
    version of this application.

    This function preserves backward compatibility with all prior supported
    configuration file formats, converting all obsolete key-value pairs of this
    configuration into their modern equivalents. Specifically:

    * Any configuration file produced by any version of this application no
      older than (i.e., at least as new as) the version specified by the
      :attr:`betse.metadata.GIT_TAG_OLDEST_BACKWARD_COMPATIBILITY` string global
      is explicitly supported by this function and hence guaranteed to be safely
      loadable with the current version of this application.
    * No configuration files produced by any older version of this application
      is explicitly supported by this function. Indeed, these files are unlikely
      to be safely loadable with the current version of this application. These
      files *must* be manually upgraded by end users to conform with the current
      configuration format.
    '''

    #FIXME: Excise this hack *AFTER* refactoring the codebase to use the
    #preferable "p._config" attribute defined above. Since the
    #yaml_alias() data descriptor expects that private attribute rather
    #than this public attribute, this public attribute is unhelpful now.
    p.config = p._conf

    # Upgrade this configuration to each successive format. For safety, each
    # upgrade is performed in strict chronological order.
    _upgrade_sim_conf_to_0_5_0(p)
    _upgrade_sim_conf_to_0_5_2(p)
    _upgrade_sim_conf_to_0_6_0(p)
    _upgrade_sim_conf_to_0_7_1(p)

# ....................{ UPGRADERS ~ 0.5.0                  }....................
@type_check
def _upgrade_sim_conf_to_0_5_0(p: Parameters) -> None:
    '''
    Upgrade the in-memory contents of the passed simulation configuration to
    reflect the newest structure of these contents expected by version 0.5.0
    (i.e., "Happy Hodgkin") of this application.
    '''

    # Log this upgrade attempt.
    logs.log_debug('Upgrading simulation configuration to 0.5.0 format...')

    # Localize configuration subdictionaries for convenience.
    results_dict = p._conf['results options']

    # For backward compatibility, convert the prior into the current
    # configuration format.
    if not (
        'while solving' in results_dict and
        'after solving' in results_dict and
        'save' in results_dict
    ):
        # Log a non-fatal warning.
        logs.log_warning(
            'Config file results options '
            '"while solving", "after solving", and/or "save" not found. '
            'Repairing to preserve backward compatibility. '
            'Consider upgrading to the newest config file format!',
        )

        # For convenience, localize configuration subdictionaries.
        anim_save = results_dict['save animations']
        anim_save_frames = anim_save['frames']

        # Convert the prior into the current configuration format.
        results_dict['while solving'] = {
            'animations': {
                'enabled': (
                            results_dict['plot while solving'] or
                            results_dict['save solving plot']
                ),
                'show':    results_dict['plot while solving'],
                'save':    results_dict['save solving plot'],
            },
        }
        results_dict['after solving'] = {
            'plots': {
                'enabled': (
                            results_dict['display plots'] or
                            results_dict['automatically save plots']
                ),
                'show':    results_dict['display plots'],
                'save':    results_dict['automatically save plots'],
            },
            'animations': {
                'enabled': results_dict['create all animations'],
                'show':    results_dict['display plots'],
                'save':    anim_save_frames['enabled'],
            },
        }
        results_dict['save'] = {
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
                    'enabled': results_dict['export data to file'],
                    'filetype': 'csv',
                },
                'vmem': {
                    'enabled': results_dict['export 2D data to file'],
                    'filetype': 'csv',
                },
            }
        }

    after_solving_anims = results_dict['after solving']['animations']
    after_solving_plots = results_dict['after solving']['plots']

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

    while_solving_anims = results_dict['while solving']['animations']

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

    if 'plot networks single cell' not in results_dict:
        # Log a non-fatal warning.
        logs.log_warning(
            'Config file setting "results options" -> '
            '"plot networks single cell" not found. '
            'Repairing to preserve backward compatibility. '
            'Consider upgrading to the newest config file format!',
        )

        # Default the value for this dictionary key to the empty list.
        results_dict['plot networks single cell'] = results_dict[
            'plot single cell graphs']

    # For each pipelined animation and cell cluster plot...
    for anim_conf in iterables.iter_items(
        results_dict['after solving']['animations']['pipeline'],
        results_dict['after solving']['plots']['cell cluster pipeline'],
    ):
        # Add the "enabled" boolean.
        if 'enabled' not in anim_conf:
            anim_conf['enabled'] = True

        # For disambiguity, rename:
        #
        # * "polarization" to "voltage_polarity".
        # * "junction_state" to "gj_permeability".
        if anim_conf['type'] == 'polarization':
            anim_conf['type'] = 'voltage_polarity'
        elif anim_conf['type'] == 'junction_state':
            anim_conf['type'] = 'gj_permeability'

# ....................{ UPGRADERS ~ 0.5.2                  }....................
@type_check
def _upgrade_sim_conf_to_0_5_2(p: Parameters) -> None:
    '''
    Upgrade the in-memory contents of the passed simulation configuration to
    reflect the newest structure of these contents expected by version 0.5.2
    (i.e., "Happiest Hodgkin") of this application.
    '''

    # Log this upgrade attempt.
    logs.log_debug('Upgrading simulation configuration to 0.5.2 format...')

    # Localize configuration subdictionaries for convenience.
    general_dict = p._conf['general options']
    results_dict = p._conf['results options']
    world_dict = p._conf['world options']

    # Support newly introduced networks plotting.
    if 'plot networks' not in results_dict:
        results_dict['plot networks'] = False
    if 'network colormap' not in results_dict:
        results_dict['network colormap'] = 'coolwarm'

    # Patch old- to new-style cell lattice types.
    if world_dict['lattice type'] == 'hexagonal':
        world_dict['lattice type'] = 'hex'
    elif world_dict['lattice type'] == 'rect':
        world_dict['lattice type'] = 'square'

    # Patch old- to new-style ion profile types.
    if general_dict['ion profile'] == 'basic_Ca':
        general_dict['ion profile'] = 'basic_ca'
    elif general_dict['ion profile'] == 'animal':
        general_dict['ion profile'] = 'mammal'
    elif general_dict['ion profile'] == 'xenopus':
        general_dict['ion profile'] = 'amphibian'
    elif general_dict['ion profile'] == 'customized':
        general_dict['ion profile'] = 'custom'

# ....................{ UPGRADERS ~ 0.6.0                  }....................
@type_check
def _upgrade_sim_conf_to_0_6_0(p: Parameters) -> None:
    '''
    Upgrade the in-memory contents of the passed simulation configuration to
    reflect the newest structure of these contents expected by version 0.6.0
    of this application.
    '''

    # Log this upgrade attempt.
    logs.log_debug('Upgrading simulation configuration to 0.6.0 format...')

    # Localize configuration subdictionaries for convenience.
    tissue_dict = p._conf['tissue profile definition']
    vars_dict = p._conf['variable settings']

    # Split ambiguously unified tissue and cut profiles into unambiguous lists.
    if 'profiles' in tissue_dict:
        tissue_dict['cut profiles'] = []
        tissue_dict['tissue profiles'] = []

        for profile in tissue_dict['profiles']:
            if profile['type'] == 'cut':
                tissue_dict['cut profiles'].append(profile)
            else:
                tissue_dict['tissue profiles'].append(profile)

    # Shift the tissue profiles list into a nested dictionary key.
    if 'tissue profiles' in tissue_dict:
        tissue_dict['tissue'] = {
            'profiles': tissue_dict['tissue profiles']
        }

    # Shift the default tissue profile into the tissue profile subsection.
    if 'default' not in tissue_dict['tissue']:
        tissue_dict['tissue']['default'] = {
            'name': vars_dict['default tissue name'],
            'insular': False,
            'diffusion constants': vars_dict['default tissue properties'],
        }

    # If double layer permittivity is undefined, define a sensible default.
    if 'double layer permittivity' not in p._conf['internal parameters']:
        p._conf['internal parameters']['double layer permittivity'] = 80.0

    # Shift the clipping bitmap mask into the default tissue profile.
    if 'clipping' in tissue_dict:
        tissue_dict['tissue']['default']['cell targets'] = tissue_dict[
            'clipping']

    # For disambiguity, upgrade tissue and cut profiles as follows:
    #
    # * Rename "bitmap" to "image" and reduce "image" from a dictionary to
    #   string by promoting its sole key-value pair to itself.
    # * Rename "random" to "percent".
    if 'image' not in tissue_dict['tissue']['default']:
        tissue_dict['tissue']['default']['image'] = tissue_dict[
            'tissue']['default']['cell targets']['bitmap']
    if isinstance(tissue_dict['tissue']['default']['image'], MappingType):
        tissue_dict['tissue']['default']['image'] = tissue_dict[
            'tissue']['default']['image']['file']
    # logs.log_debug('Default tissue: %r', tissue_dict['tissue']['default'])

    for profile in tissue_dict['tissue']['profiles']:
        if profile['cell targets']['type'] == 'bitmap':
            profile['cell targets']['type'] = 'image'
        elif profile['cell targets']['type'] == 'random':
            profile['cell targets']['type'] = 'percent'

        if 'image' not in profile['cell targets']:
            profile['cell targets']['image'] = profile['cell targets']['bitmap']
        if 'percent' not in profile['cell targets']:
            profile['cell targets']['percent'] = profile[
                'cell targets']['random']
        if isinstance(profile['cell targets']['image'], MappingType):
            profile['cell targets']['image'] = profile[
                'cell targets']['image']['file']

    for profile in tissue_dict['cut profiles']:
        if 'image' not in profile:
            profile['image'] = profile['bitmap']
        if isinstance(profile['image'], MappingType):
            profile['image'] = profile['image']['file']


@type_check
def _upgrade_sim_conf_to_0_7_1(p: Parameters) -> None:
    '''
    Upgrade the in-memory contents of the passed simulation configuration to
    reflect the newest structure of these contents expected by version 0.7.1
    of this application.
    '''

    # Log this upgrade attempt.
    logs.log_debug('Upgrading simulation configuration to 0.7.1 format...')

    # Localize configuration subdictionaries for convenience.
    grn_dict     = p._conf['gene regulatory network settings']
    results_dict = p._conf['results options']

    # If the solver type is undefined, default to the complete BETSE solver.
    p._conf.setdefault('fast solver', False)

    # If solver settings are undefined, synthesize from existing defaults.
    p._conf.setdefault('solver options', {
        'type': 'fast' if p._conf['fast solver'] else 'full'
    })

    # If the CSV file pipeline is undefined, synthesize from existing defaults.
    if 'csvs' not in results_dict['after solving']:
        results_dict['after solving']['csvs'] = {
            'save': True,
            'pipeline': [],
        }

        # If the single-cell time series CSV file is enabled, preserve.
        if results_dict['save']['data']['all']['enabled']:
            results_dict['after solving']['csvs']['pipeline'].append({
                'type': 'cell_series', 'enabled': True,})

        # If the all-cells Vmem time series CSV files are enabled, preserve.
        if results_dict['save']['data']['vmem']['enabled']:
            results_dict['after solving']['csvs']['pipeline'].append({
                'type': 'cells_vmem', 'enabled': True,})

    # If CSV save settings are undefined, synthesize from existing defaults.
    results_dict['save'].setdefault('csvs', {'filetype': 'csv',})

    # If "sim-grn" settings are undefined, synthesize from existing defaults.
    grn_dict.setdefault('sim-grn settings', {
        'run network on': 'seed',
        'save to directory': 'RESULTS/GRN',
        'save to file': 'GRN_1.betse.gz',
        'load from': None,
    })

    # If newer "sim-grn" settings are undefined, synthesize from existing
    # defaults.
    grn_dict['sim-grn settings'].setdefault('time step', 0.1)
    grn_dict['sim-grn settings'].setdefault('total time', 1.0e2)
    grn_dict['sim-grn settings'].setdefault('sampling rate', 1.0e1)
    grn_dict['sim-grn settings'].setdefault('run as sim', False)

    # If a "sim-grn" setting previously erroneously defaulting to the string
    # "None" rather than the singleton "None" is set, coerce this setting from
    # the former to the latter.
    if grn_dict['sim-grn settings']['load from'] == 'None':
        grn_dict['sim-grn settings']['load from'] = None
