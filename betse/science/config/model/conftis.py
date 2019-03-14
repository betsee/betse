#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
YAML-backed simulation subconfiguration classes for tissue and cut profiles.
'''

# ....................{ IMPORTS                           }....................
from abc import ABCMeta
from betse.lib.yaml.yamlalias import (
    yaml_alias,
    yaml_alias_float_nonnegative,
    yaml_alias_float_percent,
    yaml_enum_alias,
)
from betse.lib.yaml.abc.yamlabc import YamlABC
from betse.lib.yaml.abc.yamllistabc import YamlList, YamlListItemABC
from betse.lib.yaml.abc.yamlmixin import YamlNamedMixin
from betse.science.enum.enumconf import CellsPickerType
# from betse.util.io.log import logs
from betse.util.type.types import type_check, SequenceTypes

# ....................{ SUPERCLASSES ~ tissue             }....................
class SimConfTissueABC(YamlNamedMixin, metaclass=ABCMeta):
    '''
    Mixin of all **YAML-backed tissue profile subconfiguration** (i.e.,
    configuration of a single tissue profile parsed from the current
    YAML-formatted simulation configuration file) subclasses.

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

    # ..................{ ALIASES ~ diffusion               }..................
    Dm_Na = yaml_alias_float_nonnegative("['diffusion constants']['Dm_Na']")
    Dm_K  = yaml_alias_float_nonnegative("['diffusion constants']['Dm_K']")
    Dm_Cl = yaml_alias_float_nonnegative("['diffusion constants']['Dm_Cl']")
    Dm_Ca = yaml_alias_float_nonnegative("['diffusion constants']['Dm_Ca']")
    Dm_M  = yaml_alias_float_nonnegative("['diffusion constants']['Dm_M']")
    Dm_P  = yaml_alias_float_nonnegative("['diffusion constants']['Dm_P']")

# ....................{ SUBCLASSES ~ tissue               }....................
class SimConfTissueDefault(SimConfTissueABC, YamlABC):
    '''
    YAML-backed **default tissue profile** (i.e., profile applied to all cells
    *not* already targeted by another tissue profile) subconfiguration,
    encapsulating the configuration of a single tissue profile unconditionally
    applicable to all cells parsed from a dictionary configuring at least this
    profile in the current YAML-formatted simulation configuration file.

    Attributes (Cell Picker)
    ----------
    picker_image_filename : str
        Absolute or relative filename of the image mask whose pure-black pixels
        define the shape of the cell cluster to be populated with cells. See
        the :attr:`SimConfTissueListItem.picker_image_filename` variable for
        details.
    '''

    # ..................{ ALIASES ~ picker                  }..................
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
        **insular** (i.e., prevented from connecting to cells in other
        tissues), implying these gap junctions to be strictly intra-tissue.

    Attributes (Cell Picker)
    ----------
    picker_type : CellsPickerType
        Type of **tissue profile picker** (i.e., object assigning a subset of
        all cells matching some criteria to this tissue profile).
    picker_cells_color : str
        **Hexadecimal-formatted color** (i.e., string of six hexadecimal digits
        specifying this color's red, green, and blue components) of all circles
        within the vector image (defined by the ``cells from svg`` setting in
        the current simulation configuration) to be assigned to this tissue.
        Ignored unless :attr:`picker_type` is :attr:`CellsPickerType.COLOR`.
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
        define the region of the cell cluster whose cells are all to be
        assigned to this tissue. This image *must*:

        * Be square (i.e., have equal width and height).
        * Contain no alpha transparency layer.

        Ignored unless :attr:`picker_type` is :attr:`CellsPickerType.IMAGE`.
    '''

    # ..................{ ALIASES                           }..................
    is_gj_insular = yaml_alias("['insular']", bool)

    # ..................{ ALIASES ~ picker                  }..................
    picker_type = yaml_enum_alias("['cell targets']['type']", CellsPickerType)

    #FIXME: Create a new yaml_alias_color() data descriptor validating this
    #string to be a valid hexadecimal-formatted color.
    #FIXME: In the YAML, if color specifications don't have quotations around
    #them (e.g. '808080') the yaml sometimes thinks all numerical hex codes are
    #strings (e.g. 008080), and sometimes integers (e.g. 800080), in which case
    #the whole thing crashed with a type error! I'm not quite sure how to force
    #it to read in a string... For now, I've put quotes around the hex code
    #numeral...
    picker_cells_color = yaml_alias("['cell targets']['color']", str)

    picker_cells_index = yaml_alias(
        "['cell targets']['indices']", SequenceTypes)
    picker_cells_percent = yaml_alias_float_percent(
        "['cell targets']['percent']")
    picker_image_filename = yaml_alias(
        "['cell targets']['image']", str)

    # ..................{ SUPERCLASS                        }..................
    @classmethod
    @type_check
    def make_default(cls, yaml_list: YamlList) -> YamlListItemABC:

        # Duplicate the first tissue profile in our default configuration file.
        return cls._make_loaded(conf={
            'name': yaml_list.get_item_name_uniquified('Tissue ({})'),
            'insular': True,
            'diffusion constants': {
                'Dm_Na': 1.0e-18,
                'Dm_K':  15.0e-18,
                'Dm_Cl': 2.0e-18,
                'Dm_Ca': 1.0e-18,
                'Dm_M':  1.0e-18,
                'Dm_P':  0.0,
            },
            'cell targets': {
                'type': 'all',
                'color': 'ff0000',
                'image': 'geo/circle/circle_base.png',
                'indices': [3, 14, 15, 9, 265],
                'percent': 50,
            },
        })

# ....................{ SUBCLASSES ~ cut                  }....................
class SimConfCutListItem(YamlNamedMixin, YamlListItemABC):
    '''
    YAML-backed cut profile list item subconfiguration, encapsulating the
    configuration of a single cut profile parsed from a list of these profiles
    in the current YAML-formatted simulation configuration file.

    Attributes
    ----------
    picker_image_filename : str
        Absolute or relative filename of the image mask whose pure-black pixels
        define the region of the cell cluster whose cells are all to be removed
        by this cut profile. See the
        :attr:`SimConfTissueListItem.picker_image_filename` variable for
        details.
    '''

    # ..................{ ALIASES                           }..................
    picker_image_filename = yaml_alias("['image']", str)

    # ..................{ SUPERCLASS                        }..................
    @classmethod
    @type_check
    def make_default(cls, yaml_list: YamlList) -> YamlListItemABC:

        # Duplicate the first tissue profile in our default configuration file.
        return cls._make_loaded(conf={
            'name': yaml_list.get_item_name_uniquified('Cut ({})'),
            'image': 'geo/circle/wedge.png',
        })
