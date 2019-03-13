#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
**Simulation export subconfiguration** (i.e., YAML-backed list item configuring
the exportation of one or more external files from the ``plot`` subcommand when
passed the simulation configuration file containing this item) superclasses.
'''

# ....................{ IMPORTS                           }....................
from betse.lib.yaml.abc.yamllistabc import YamlListItemABC
from betse.lib.yaml.abc.yamlmixin import (
    YamlBooledMixin, YamlNamedMixin, YamlTypedMixin)
# from betse.util.type.types import type_check

# ....................{ SUPERCLASSES                      }....................
class SimConfExportABC(
    YamlBooledMixin,
    YamlNamedMixin,
    YamlTypedMixin,
    YamlListItemABC,
):
    '''
    Abstract base class of all **simulation export subconfiguration** (i.e.,
    YAML-backed list item configuring the exportation of one or more external
    files from the ``plot`` subcommand when passed the simulation configuration
    file containing this item) subclasses.
    '''

    # That's all she wrote.
    pass
