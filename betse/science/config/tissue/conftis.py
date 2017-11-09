#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
YAML-backed simulation subconfiguration classes for tissue and cut profiles.
'''

# ....................{ IMPORTS                            }....................
from betse.lib.yaml.yamlalias import yaml_alias
from betse.lib.yaml.abc.yamllistabc import YamlListItemABC
# from betse.util.io.log import logs
# from betse.util.type.types import type_check

# ....................{ SUBCLASSES                         }....................
#FIXME: Define a "SimConfTissueProfileListItem" class as well.
#FIXME: Actually leverage this in "Parameters".
class SimConfTissueProfileListItem(YamlListItemABC):
    '''
    YAML-backed tissue profile subconfiguration, encapsulating the configuration
    of a single tissue profile parsed from a list of these profiles in the
    current YAML-formatted simulation configuration file.

    Attributes
    ----------
    is_gj_insular : bool
        ``True`` only if gap junctions originating at cells in this tissue are
        **insular** (i.e., prevented from connecting to cells in other tissues),
        implying these gap junctions to be strictly intra-tissue.
    name : str
        Arbitrary string uniquely identifying this tissue profile.

    Attributes (Membrane Diffusion Constants)
    ----------

    Attributes (Cell Selection)
    ----------
    '''

    # ..................{ ALIASES                            }..................
    #FIXME: Alias all remaining tissue profile settings.
    is_gj_insular = yaml_alias("['insular']", bool)
    name = yaml_alias("['name']", str)

    # ..................{ CLASS                              }..................
    @classmethod
    def make_default(cls) -> YamlListItemABC:

        # Duplicate the default tissue listed first in our default YAML file.
        return SimConfTissueProfileListItem(conf={
            #FIXME: We'll need to leverage a closure to dynamically fabricate a
            #new tissue profile *GUARANTEED* to be unique. Actually, wait...
            #Since this method does *NOT* accept the current parent "YamlList"
            #object, it has no means of actually guaranteeing uniqueness. Ergo,
            #we'll need to refactor the superclass method signature as follows:
            #
            #     @classmethod
            #     def make_default(cls, parent_list: 'YamlList') -> YamlListItemABC:
            'name': 'tissue (1)',
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
        })

    # ..................{ INITIALIZERS                       }..................
    def __init__(self, *args, **kwargs) -> None:

        # Initialize our superclass with all passed parameters.
        super().__init__(*args, **kwargs)
