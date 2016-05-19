#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Medium-level simulation configuration classes wrapping low-level YAML
dictionaries deserialized from disk.
'''

#FIXME: Refactor all functions defined by "betse.science.config.sim_config"
#into methods of the "SimConfigWrapper" class defined below.

# ....................{ IMPORTS                            }....................
# import yaml
# from betse import pathtree
from betse.science.config import sim_config
# from betse.util.type import types
# from betse.util.io.log import logs
# from betse.util.path import dirs, files, paths

# ....................{ CLASSES                            }....................
class SimConfigWrapper(object):
    '''
    Medium-level simulation configuration wrapper wrapping a low-level
    dictionary deserialized from a YAML-formatted simulation configuration file.

    This wrapper superficially wraps this dictionary with convenience methods
    safely modifying this dictionary. This wrapper is lower level than the
    high-level `Parameters` simulation configuration, which transforms this
    dictionary into numerous high-level objects rather than merely wrapping this
    dictionary. While the `Parameters` configuration is principally used by
    backend simulation modelling in the `betse.science` package, this wrapper is
    principally used by frontend logic modifying simulation configurations on
    behalf of either interactive users (e.g., BETSE's GUI) or automated tests.

    Attributes
    ----------
    _config : dict
        Low-level dictionary deserialized from `_filename`.
    _filename : str
        Absolute or relative path of the YAML-formatted simulation configuration
        file deserialized into `_config`.
    '''

    # ..................{ INITIALIZE                         }..................
    def __init__(self, filename: str) -> None:
        '''
        Wrap the low-level dictionary deserialized from the passed
        YAML-formatted simulation configuration file.

        Parameters
        ----------
        filename : str
            Absolute or relative path of this file.
        '''

        # Classify the passed parameters.
        self._filename = filename

        # Deserialize this YAML file into a dictionary.
        self._config = sim_config.read(filename)


    # ..................{ GETTERS                            }..................

    # ..................{ SETTERS                            }..................

    # ..................{ WRITERS                            }..................
    def overwrite(self) -> None:
        '''
        Reserialize the current low-level configuration dictionary to the
        current configuration file.

        This method silently overwrites the contents of this file with the
        (possibly modified) contents of this dictionary.
        '''

        #FIXME: We'll probably need to explicitly delete this file if it
        #currently exists first (in a safe manner hopefully avoiding race
        #conditions).
        self.write(self._filename)


    #FIXME: Implement us up.
    def write(self, filename: str) -> None:
        '''
        Serialize the current low-level configuration dictionary to the passed
        simulation configuration file in YAML format.

        If this file already exists, an exception is raised.
        '''

        pass

    # ..................{ DISABLERS                          }..................
    def disable_interaction(self) -> None:
        '''
        Disable all simulation configuration options either requiring
        interactive user input _or_ displaying graphical output intended for
        interactive user consumption (e.g., plots, animations).

        This method is intended to be called by non-interactive automation
        (e.g., tests, scripts) expecting simulations to behave silently.
        '''

        # Disable display of all mid- and post-simulation plots and animations.
        self._config['results options']['plot while solving'] = False
        self._config['results options']['plot after solving'] = False
