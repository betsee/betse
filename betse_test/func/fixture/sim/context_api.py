#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Simulation-specific fixture context classes.

Simulation-specific fixtures typically return instances of these classes to
caller fixtures and tests, formalizing communication between functional tests.
'''

# ....................{ IMPORTS                            }....................
from betse.science.config import sim_config

# ....................{ CLASSES                            }....................
class SimTestContext(object):
    '''
    Simulation-specific test context encapsulating simulation configuration,
    state, and metadata.

    Attributes
    ----------
    config : dict
        Dictionary of all configuration data deserialized from the YAML-
        formatted file with path `config_filename`. Note that the contents of
        this in-memory dictionary differ from that of this on-disk file. For
        efficiency, callers are expected to additionally modify this dictionary
        to suite test requirements before finally reserializing this dictionary
        to this file.
    config_filename : py.path.local
        Absolute path of a temporary simulation configuration file specific to
        the parent fixture as a `py.path.local` instance, defining an
        object-oriented superset of the non-object-oriented `os.path` module.

    See Also
    ----------
    https://py.readthedocs.org/en/latest/path.html
        Official `py.path` class documentation.
    '''

    def __init__(self, config_filename : 'py.path.local') -> None:
        '''
        Initialize this test configuration.

        Parameters
        ----------
        config_filename : py.path.local
            Absolute path to which this method copies BETSE's default simulation
            configuration file. If this file already exists, an exception is
            raised.
        '''

        # Configuration filename.
        self.config_filename = config_filename

        # Configuration deserialized from this file, reducing this filename from
        # a high-level "py.path.local" instance to a low-level string.
        self.config = sim_config.read(str(self.config_filename))


    #FIXME: Implement me as a convenience for fixtures!
    def write(self) -> None:
        '''
        Write the current configuration to the current configuration file.

        If this file already exists, this file will be silently overwritten.
        '''

        pass
