#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
**Simulation export pipeline** (i.e., container of related, albeit isolated,
simulation export actions to be run iteratively) hierarchy.
'''

# ....................{ IMPORTS                           }....................
from betse.science.pipe.pipeabc import SimPipeABC
from betse.util.type.descriptor.descs import classproperty_readonly

# ....................{ SUPERCLASSES                      }....................
class SimPipeExportABC(SimPipeABC):
    '''
    Abstract base class of all **simulation export pipeline** (i.e., object
    iteratively exporting *all* variations on a single type of simulation
    export) subclasses.
    '''

    # ..................{ SUPERCLASS                        }..................
    @classproperty_readonly
    def _RUNNER_METHOD_NAME_PREFIX(self) -> str:
        return  'export_'


    @classproperty_readonly
    def _VERB_CONTINUOUS(cls) -> str:
        return 'Exporting'
