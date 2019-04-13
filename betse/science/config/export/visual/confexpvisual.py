#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
**Visual export subconfiguration** (i.e., YAML-backed wrapper exposing settings
globally applicable to all visual exports) functionality.
'''

#FIXME: Refactor all remaining visual export-centric settings currently defined
#as "Parameters" instance variables into data descriptors defined below.

# ....................{ IMPORTS                           }....................
from betse.lib.yaml.yamlalias import yaml_alias
from betse.lib.yaml.abc.yamlabc import YamlABC
# from betse.util.type.types import type_check

# ....................{ SUBCLASSES                        }....................
class SimConfExportVisual(YamlABC):
    '''
    **Visual export subconfiguration** (i.e., YAML-backed wrapper exposing
    settings globally applicable to all visual exports).

    Attributes (Cell: Indices)
    ----------
    is_show_cell_indices : bool
        ``True`` only if the 0-based indices of all cells are to be displayed
        over all cell cluster visuals (e.g., as integers situated at cell
        centres).
    single_cell_index : int
        0-based index of the cell to be visualized for all single cell visuals.
        Defaults to 0, the index assigned to the first cell guaranteed to
        exist. Note that cell indices are seed-specific and may be visualized
        by enabling the :attr:`is_show_cell_indices` boolean.
    '''

    # ..................{ ALIASES ~ cell : index            }..................
    is_show_cell_indices = yaml_alias(
        "['results options']['visuals']['cell indices']['show']", bool)

    #FIXME, Consider generalizing this into a list.
    single_cell_index = yaml_alias(
        "['results options']['visuals']['cell indices']['single cell']", int)
