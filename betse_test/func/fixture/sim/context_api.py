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
from betse.science.config.wrapper import SimConfigWrapper

# ....................{ CLASSES                            }....................
class SimTestContext(object):
    '''
    Simulation-specific test context encapsulating simulation configuration,
    state, and metadata.

    Simulation-specific fixtures typically return instances of this class as a
    means of communicating this context to other fixtures and tests.

    Attributes
    ----------
    config : SimConfigWrapper
        Simulation configuration wrapper wrapping the low-level dictionary
        deserialized from the YAML-formatted simulation configuration file with
        path `config_filename`. Note that the contents of this in-memory
        dictionary be desynchronized from those of this file. For efficiency,
        callers may modify this dictionary to suite test requirements _before_
        reserializing this dictionary back to this file.
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
        Initialize this test context.

        Parameters
        ----------
        config_filename : py.path.local
            Absolute path to which this method copies BETSE's default simulation
            configuration file. If this file already exists, an exception is
            raised.
        '''

        # Configuration filename. While the "self.config" object classified
        # below also provides this filename as a low-level string, classify this
        # higher-level "py.path.local" instance for use in fixtures and tests.
        self._config_filename = config_filename

        # Configuration deserialized from this file, reducing this filename from
        # a high-level "py.path.local" instance to a low-level string.
        self.config = SimConfigWrapper(filename=str(self.config_filename))

        # Disable configuration options either requiring interactive input *OR*
        # displaying interactive output.
        self.config.disable_interaction()
