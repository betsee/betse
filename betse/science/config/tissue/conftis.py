#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
YAML-backed simulation subconfiguration classes for tissue and cut profiles.
'''

# ....................{ IMPORTS                            }....................
from abc import ABCMeta  #, abstractproperty
from betse.lib.yaml.yamlalias import (
    yaml_alias,
    yaml_alias_float_nonnegative,
    yaml_alias_float_percent,
    yaml_enum_alias,
)
from betse.lib.yaml.abc import yamllistabc
from betse.lib.yaml.abc.yamlabc import YamlABC
from betse.lib.yaml.abc.yamllistabc import YamlList, YamlListItemABC
# from betse.util.io.log import logs
from betse.util.type.enums import make_enum
from betse.util.type.mapping import maputil
from betse.util.type.types import type_check, MappingType, SequenceTypes

# ....................{ ENUMS                              }....................
TissueProfilePickerType = make_enum(
    class_name='TissueProfilePickerType',
    member_names=('ALL', 'IMAGE', 'INDICES', 'PERCENT',))
'''
Enumeration of all supported types of **tissue profile pickers** (i.e., objects
assigning a subset of all cells matching some criteria to the corresponding
tissue profile).

Attributes
----------
ALL : enum
    All-inclusive tissue picker, unconditionally matching *all* cells.
IMAGE : enum
    Image-based tissue picker, matching all cells residing inside the colored
    pixel area defined by an associated on-disk image mask file.
INDICES : enum
    Cell indices-based tissue picker, matching all cells whose indices are
    defined by a given sequence.
PERCENT : enum
    Randomized cell picker, randomly matching a given percentage of all cells.
'''

# ....................{ DEFAULTS                           }....................
def _get_tissue_conf_default() -> MappingType:
    '''
    Low-level dictionary suitable for use as the default value for the ``conf``
    parameter of YAML-centric methods (e.g., :meth:`YamlListItemABC.__init__`)
    generally pertaining to tissue profiles.

    For orthogonality, the contents of this dictionary exactly duplicate the
    default tissue profile listed first in our default YAML-formatted
    simulation configuration file.
    '''

    return {
        'name': 'base',
        'insular': True,
        'diffusion constants': {
            'Dm_Na': 1.0e-18,
            'Dm_K': 15.0e-18,
            'Dm_Cl': 2.0e-18,
            'Dm_Ca': 1.0e-18,
            'Dm_H': 1.0e-18,
            'Dm_M': 1.0e-18,
            'Dm_P': 0.0,
        },
        'cell targets': {
            'type': 'all',
            'bitmap': {'file': 'geo/circle/circle_base.png'},
            'indices': [3, 14, 15, 9, 265],
            'random': 50,
        },
    }

# ....................{ SUPERCLASSES                       }....................
class SimConfTissueABC(object, metaclass=ABCMeta):
    '''
    Abstract mixin generalizing implementation common to all YAML-backed tissue
    profile subconfiguration subclasses.

    This mixin encapsulates configuration of a single tissue profile parsed from
    the current YAML-formatted simulation configuration file. For generality,
    this mixin provides no support for the YAML ``name`` key or the suite of
    YAML keys pertaining to tissue pickers.

    Attributes
    ----------
    is_gj_insular : bool
        ``True`` only if gap junctions originating at cells in this tissue are
        **insular** (i.e., prevented from connecting to cells in other tissues),
        implying these gap junctions to be strictly intra-tissue.

    Attributes (Constant)
    ----------
    Dm_Na : float
        Sodium (Na+) membrane diffusion constant in m2/s.
    Dm_K : float
        Potassium (K+) membrane diffusion constant in m2/s.
    Dm_Cl : float
        Chloride (Cl-) membrane diffusion constant in m2/s.
    Dm_Ca : float
        Calcium (Ca2+) membrane diffusion constant in m2/s.
    Dm_H : float
        Hydrogen (H+) membrane diffusion constant in m2/s.
    Dm_M : float
        Charge balance anion (M-) membrane diffusion constant in m2/s.
    Dm_P : float
        Protein (P-) membrane diffusion constant in m2/s.
    '''

    # ..................{ ALIASES                            }..................
    is_gj_insular = yaml_alias("['insular']", bool)

    # ..................{ ALIASES ~ diffusion                }..................
    Dm_Na = yaml_alias_float_nonnegative("['diffusion constants']['Dm_Na']")
    Dm_K  = yaml_alias_float_nonnegative("['diffusion constants']['Dm_K']")
    Dm_Cl = yaml_alias_float_nonnegative("['diffusion constants']['Dm_Cl']")
    Dm_Ca = yaml_alias_float_nonnegative("['diffusion constants']['Dm_Ca']")
    Dm_H  = yaml_alias_float_nonnegative("['diffusion constants']['Dm_H']")
    Dm_M  = yaml_alias_float_nonnegative("['diffusion constants']['Dm_M']")
    Dm_P  = yaml_alias_float_nonnegative("['diffusion constants']['Dm_P']")

# ....................{ SUBCLASSES                         }....................
class SimConfTissueDefault(SimConfTissueABC, YamlABC):
    '''
    YAML-backed default tissue profile subconfiguration, encapsulating the
    configuration of a single tissue profile applicable to all cells parsed from
    a dictionary configuring at least this profile in the current YAML-formatted
    simulation configuration file.
    '''

    # ..................{ CLASS                              }..................
    @classmethod
    @type_check
    def make_default(cls, yaml_list: YamlList) -> YamlListItemABC:

        # Create and return the default YAML-backed tissue profile.
        return SimConfTissueListItem(conf=_get_tissue_conf_default())

# ....................{ SUBCLASSES                         }....................
#FIXME: Define a similar "SimConfCutListItem" class as well.
#FIXME: Actually leverage this in "Parameters".
class SimConfTissueListItem(SimConfTissueABC, YamlListItemABC):
    '''
    YAML-backed tissue profile list item subconfiguration, encapsulating the
    configuration of a single tissue profile parsed from a list of these
    profiles in the current YAML-formatted simulation configuration file.

    Attributes
    ----------
    name : str
        Arbitrary string uniquely identifying this tissue profile in this list.

    Attributes (Picker)
    ----------
    picker_type : TissueProfilePickerType
        Type of **tissue profile picker** (i.e., object assigning a subset of
        all cells matching some criteria to this tissue profile).
    picker_cells_index : SequenceTypes
        Ignored unless :attr:`picker_type` is
        :attr:`TissueProfilePickerType.INDICES`.
    picker_cells_percent : float
        **Percentage** (i.e., floating point number in the range ``[0.0,
        100.0]``) of the total cell population to randomly match. Ignored unless
        :attr:`picker_type` is :attr:`TissueProfilePickerType.PERCENT`.
    picker_image_filename : str
        Ignored unless :attr:`picker_type` is
        :attr:`TissueProfilePickerType.IMAGE`.
    '''

    # ..................{ ALIASES                            }..................
    name = yaml_alias("['name']", str)

    # ..................{ ALIASES ~ picker                   }..................
    picker_type = yaml_enum_alias(
        "['cell targets']['type']", TissueProfilePickerType)
    picker_cells_index = yaml_alias(
        "['cell targets']['indices']", SequenceTypes)
    picker_cells_percent = yaml_alias_float_percent(
        "['cell targets']['random']")
    picker_image_filename = yaml_alias(
        "['cell targets']['bitmap']['file']", str)

    # ..................{ CLASS                              }..................
    @classmethod
    @type_check
    def make_default(cls, yaml_list: YamlList) -> YamlListItemABC:

        # Name of this tissue profile unique to this list.
        tissue_name = yamllistabc.get_list_item_name_unique(
            yaml_list=yaml_list, name_format='tissue ({{}})')

        # Low-level YAML-centric dictionary configuring this tissue profile,
        # produced by merging (in order):
        tissue_conf = maputil.merge(
            # The default dictionary for tissue profiles.
            _get_tissue_conf_default(),

            # The default dictionary for tissue profile list items, duplicating
            # the first such item in our default YAML file. The list-specific
            # "name" key of this dictionary takes precedence over the generic
            # "name" key of the prior dictionary. This dictionary *MUST* be
            # listed last.
            {
                'name': tissue_name,
                'cell targets': {
                    'type': 'all',
                    'bitmap': {'file': 'geo/circle/circle_base.png'},
                    'indices': [3, 14, 15, 9, 265],
                    'random': 50,
                },
            },
        )

        # Create and return the equivalent YAML-backed tissue profile list item.
        return SimConfTissueListItem(conf=tissue_conf)
