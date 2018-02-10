#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
YAML-backed simulation subconfigurations for exporting comma-separated value
(CSV) files.
'''

# ....................{ IMPORTS                            }....................
from betse.lib.yaml.yamlalias import yaml_alias
from betse.lib.yaml.abc.yamlabc import YamlABC
from betse.lib.yaml.abc.yamllistabc import (
    YamlList, YamlListItemABC, YamlListItemTypedABC)
from betse.util.type.types import type_check
# from betse.util.type.types import type_check, MappingType, SequenceTypes

# ....................{ SUBCLASSES                         }....................
class SimConfExportCSVs(YamlABC):
    '''
    YAML-backed subconfiguration for exporting *all* comma-separated value (CSV)
    files enabled by the current YAML-formatted simulation configuration file.

    Attributes (After Solving)
    ----------
    is_after_sim_save : bool
        ``True`` only if this configuration saves post-simulation CSV files.
    anims_after_sim : YamlList
        YAML-backed list of all post-simulation CSV files to be animated.
        Ignored if :attr:`is_after_sim_save` is ``False``.

    Attributes (Save)
    ----------
    csv_filetype : str
        Filetype of all CSV files saved by this configuration. Ignored if
        :attr:`is_after_sim_save` is ``False``.
    '''

    # ..................{ ALIASES ~ after                    }..................
    is_after_sim_save = yaml_alias(
        "['results options']['after solving']['csvs']['save']", bool)

    # ..................{ ALIASES ~ save                     }..................
    filetype = yaml_alias(
        "['results options']['save']['csvs']['filetype']", str)

    # ..................{ INITIALIZERS                       }..................
    def __init__(self, *args, **kwargs) -> None:

        # Initialize our superclass with all passed parameters.
        super().__init__(*args, **kwargs)

        # Encapsulate low-level lists of dictionaries with high-level wrappers.
        self.csvs_after_sim = SimConfExportCSV.make_list()

    # ..................{ LOADERS                            }..................
    def load(self, *args, **kwargs) -> None:

        # Load our superclass with all passed arguments.
        super().load(*args, **kwargs)

        # Load all subconfigurations of this configuration.
        self.csvs_after_sim.load(conf=self._conf[
            'results options']['after solving']['csvs']['pipeline'])


    def unload(self) -> None:

        # Unload our superclass.
        super().unload()

        # Unload all subconfigurations of this configuration.
        self.csvs_after_sim.unload()

# ....................{ SUBCLASSES : item                  }....................
class SimConfExportCSV(YamlListItemTypedABC):
    '''
    YAML-backed subconfiguration list item for exporting a comma-separated value
    (CSV) file enabled by a list of these exports in the current YAML-formatted
    simulation configuration file.
    '''

    # ..................{ MAKERS                             }..................
    @classmethod
    @type_check
    def make_default(cls, yaml_list: YamlList) -> YamlListItemABC:

        # Duplicate the default CSV file listed first in our default YAML file.
        return SimConfExportCSV().load(conf={
            'type': 'cell_Series',
            'enabled': True,
        })
