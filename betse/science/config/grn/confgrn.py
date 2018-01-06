#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
YAML-backed gene regulatory network (GRN) subconfigurations.
'''

# ....................{ IMPORTS                            }....................
from betse.lib.yaml.yamlalias import yaml_alias
from betse.lib.yaml.abc.yamlabc import YamlFileABC
# from betse.science.config.grn.confgrnyadda import SimConfGrnYadda
from betse.util.type.types import type_check  #, MappingType, SequenceTypes

# ....................{ SUBCLASSES                         }....................
class SimConfGrnFile(YamlFileABC):
    '''
    YAML-backed in-memory and on-disk gene regulatory network (GRN)
    subconfiguration, encapsulating a low-level container of *all* GRN-related
    configuration settings both loaded from and saved back to a YAML-formatted
    configuration file.

    Attributes (General)
    ----------
    '''

    # ..................{ ALIASES                            }..................
    # is_overlay_current = yaml_alias(
    #     "['results options']['overlay currents']", bool)

    # ..................{ INITIALIZERS                       }..................
    def __init__(self, *args, **kwargs) -> None:

        # Initialize our superclass with all passed parameters.
        super().__init__(*args, **kwargs)

        # Encapsulate low-level dictionaries with high-level wrappers.
        # self.anim_while_sim = SimConfVisualCellsEmbedded()

        # Encapsulate low-level lists of dictionaries with high-level wrappers.
        # self.anims_after_sim = SimConfVisualCellsListItem.make_list()

    # ..................{ LOADERS                            }..................
    def load(self, *args, **kwargs) -> None:

        # Load our superclass with all passed arguments.
        super().load(*args, **kwargs)

        # Load all subconfigurations of this configuration.
        # self.anim_while_sim.load(conf=self._conf[
        #     'results options']['while solving']['animations'])


    def unload(self) -> None:

        # Unload our superclass.
        super().unload()

        # Unload all subconfigurations of this configuration.
        # self.anim_while_sim.unload()
