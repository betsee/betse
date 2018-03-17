#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
YAML-backed simulation subconfiguration classes for tissue and cut profiles.
'''

# ....................{ IMPORTS                            }....................
from abc import ABCMeta
from betse.lib.yaml.yamlalias import (
    yaml_alias,
    yaml_alias_float_nonnegative,
    yaml_alias_float_percent,
    yaml_enum_alias,
)
from betse.lib.yaml.abc.yamlabc import YamlABC
from betse.lib.yaml.abc.yamllistabc import YamlList, YamlListItemABC
from betse.science.config.confenum import CellsPickerType
# from betse.util.io.log import logs
from betse.util.type.types import type_check, SequenceTypes

# ....................{ SUPERCLASSES ~ tissue              }....................
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
    name : str
        Arbitrary string uniquely identifying this tissue profile in the list of
        all tissue profiles for this simulation.

    Attributes (Membrane Diffusion)
    ----------
    Dm_Na : float
        Sodium (Na+) membrane diffusion constant in m2/s.
    Dm_K : float
        Potassium (K+) membrane diffusion constant in m2/s.
    Dm_Cl : float
        Chloride (Cl-) membrane diffusion constant in m2/s.
    Dm_Ca : float
        Calcium (Ca2+) membrane diffusion constant in m2/s.
    Dm_M : float
        Charge balance anion (M-) membrane diffusion constant in m2/s.
    Dm_P : float
        Protein (P-) membrane diffusion constant in m2/s.
    '''

    # ..................{ ALIASES                            }..................
    name = yaml_alias("['name']", str)

    # ..................{ ALIASES ~ diffusion                }..................
    Dm_Na = yaml_alias_float_nonnegative("['diffusion constants']['Dm_Na']")
    Dm_K  = yaml_alias_float_nonnegative("['diffusion constants']['Dm_K']")
    Dm_Cl = yaml_alias_float_nonnegative("['diffusion constants']['Dm_Cl']")
    Dm_Ca = yaml_alias_float_nonnegative("['diffusion constants']['Dm_Ca']")
    # Dm_H  = yaml_alias_float_nonnegative("['diffusion constants']['Dm_H']")
    Dm_M  = yaml_alias_float_nonnegative("['diffusion constants']['Dm_M']")
    Dm_P  = yaml_alias_float_nonnegative("['diffusion constants']['Dm_P']")

# ....................{ SUBCLASSES ~ tissue                }....................
class SimConfTissueDefault(SimConfTissueABC, YamlABC):
    '''
    YAML-backed default tissue profile subconfiguration, encapsulating the
    configuration of a single tissue profile unconditionally applicable to all
    cells parsed from a dictionary configuring at least this profile in the
    current YAML-formatted simulation configuration file.

    Attributes (Cell Picker)
    ----------
    picker_image_filename : str
        Absolute or relative filename of the image mask whose pure-black pixels
        define the shape of the cell cluster to be populated with cells. See
        the :attr:`SimConfTissueListItem.picker_image_filename` variable for
        details.
    '''

    # ..................{ ALIASES ~ picker                   }..................
    picker_image_filename = yaml_alias("['image']", str)


class SimConfTissueListItem(SimConfTissueABC, YamlListItemABC):
    '''
    YAML-backed tissue profile list item subconfiguration, encapsulating the
    configuration of a single tissue profile parsed from a list of these
    profiles in the current YAML-formatted simulation configuration file.

    Attributes
    ----------
    is_gj_insular : bool
        ``True`` only if gap junctions originating at cells in this tissue are
        **insular** (i.e., prevented from connecting to cells in other tissues),
        implying these gap junctions to be strictly intra-tissue.

    Attributes (Cell Picker)
    ----------
    picker_type : CellsPickerType
        Type of **tissue profile picker** (i.e., object assigning a subset of
        all cells matching some criteria to this tissue profile).
    picker_cells_index : SequenceTypes
        Sequence of the indices of all cells to be assigned to this tissue.
        Ignored unless :attr:`picker_type` is :attr:`CellsPickerType.INDICES`.
    picker_cells_percent : float
        **Percentage** (i.e., floating point number in the range ``[0.0,
        100.0]``) of the total cell population to be randomly assigned to this
        tissue. Ignored unless :attr:`picker_type` is
        :attr:`CellsPickerType.PERCENT`.
    picker_image_filename : str
        Absolute or relative filename of the image mask whose pure-black pixels
        (i.e., pixels whose red, green, and blue color components are all 0)
        define the region of the cell cluster whose cells are all to be assigned
        to this tissue. This image *must*:
        * Be square (i.e., have equal width and height).
        * Contain no alpha transparency layer.
        Ignored unless :attr:`picker_type` is :attr:`CellsPickerType.IMAGE`.
    '''

    # ..................{ ALIASES                            }..................
    is_gj_insular = yaml_alias("['insular']", bool)

    # ..................{ ALIASES ~ picker                   }..................
    picker_type = yaml_enum_alias("['cell targets']['type']", CellsPickerType)
    picker_cells_index = yaml_alias(
        "['cell targets']['indices']", SequenceTypes)
    picker_cells_percent = yaml_alias_float_percent(
        "['cell targets']['percent']")
    picker_image_filename = yaml_alias(
        "['cell targets']['image']", str)

    # ..................{ CLASS                              }..................
    @classmethod
    @type_check
    def make_default(cls, yaml_list: YamlList) -> YamlListItemABC:

        # Name of this tissue profile unique to this list.
        tissue_name = yaml_list.get_item_name_unique(
            name_format='tissue ({{}})')

        # Create and return the equivalent YAML-backed tissue profile list item,
        # duplicating the first such item in our default YAML file.
        return SimConfTissueListItem().load(conf={
            'name': tissue_name,
            'insular': True,
            'diffusion constants': {
                'Dm_Na': 1.0e-18,
                'Dm_K': 15.0e-18,
                'Dm_Cl': 2.0e-18,
                'Dm_Ca': 1.0e-18,
                # 'Dm_H': 1.0e-18,
                'Dm_M': 1.0e-18,
                'Dm_P': 0.0,
            },
            'cell targets': {
                'type': 'all',
                'image': 'geo/circle/circle_base.png',
                'indices': [3, 14, 15, 9, 265],
                'percent': 50,
            },
        })

# ....................{ SUBCLASSES ~ cut                   }....................
class SimConfCutListItem(YamlListItemABC):
    '''
    YAML-backed cut profile list item subconfiguration, encapsulating the
    configuration of a single cut profile parsed from a list of these profiles
    in the current YAML-formatted simulation configuration file.

    Attributes
    ----------
    name : str
        Arbitrary string uniquely identifying this cut profile in this list.
    picker_image_filename : str
        Absolute or relative filename of the image mask whose pure-black pixels
        define the region of the cell cluster whose cells are all to be removed
        by this cut profile. See the
        :attr:`SimConfTissueListItem.picker_image_filename` variable for
        details.
    '''

    # ..................{ ALIASES                            }..................
    name = yaml_alias("['name']", str)
    picker_image_filename = yaml_alias("['image']", str)

    # ..................{ CLASS                              }..................
    @classmethod
    @type_check
    def make_default(cls, yaml_list: YamlList) -> YamlListItemABC:

        # Name of this cut profile unique to this list.
        cut_name = yaml_list.get_item_name_unique(name_format='cut ({{}})')

        # Create and return the equivalent YAML-backed cut profile list item,
        # duplicating the first such item in our default YAML file.
        return SimConfTissueListItem().load(conf={
            'name': cut_name,
            'image': 'geo/circle/wedge.png',
        })
