#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Simulation configuration backward compatibility facilities.
'''

# ....................{ IMPORTS                            }....................
from betse.util.io.log import logs
from betse.util.type.types import type_check

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
def upgrade_sim_conf(p: 'betse.science.parameters.Parameters') -> None:
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

# ....................{ UPGRADERS ~ 0.5.0                  }....................
@type_check
def _upgrade_sim_conf_to_0_5_0(
    p: 'betse.science.parameters.Parameters') -> None:
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

# ....................{ UPGRADERS ~ 0.5.2                  }....................
@type_check
def _upgrade_sim_conf_to_0_5_2(
    p: 'betse.science.parameters.Parameters') -> None:
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
    if world_dict['lattice type'] == 'hex':
        world_dict['lattice type'] = 'hexagonal'
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
def _upgrade_sim_conf_to_0_6_0(
    p: 'betse.science.parameters.Parameters') -> None:
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

            del profile['type']

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

    # For disambiguity, rename "bitmap" to "image" in tissue and cut profiles.
    if 'image' not in tissue_dict['tissue']['default']:
        tissue_dict['tissue']['default']['image'] = tissue_dict[
            'tissue']['default']['cell targets']['bitmap']
    for profile in tissue_dict['tissue']['profiles']:
        if profile['cell targets']['type'] == 'bitmap':
            profile['cell targets']['type'] = 'image'
        if 'image' not in profile['cell targets']:
            profile['cell targets']['image'] = profile['cell targets']['bitmap']
    for profile in tissue_dict['cut profiles']:
        if 'image' not in profile:
            profile['image'] = profile['bitmap']
