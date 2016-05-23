#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level classes returned by simulation configuration fixtures to other
fixtures and tests requiring a simulation configuration.
'''

# ....................{ IMPORTS                            }....................
from betse.science.config.wrapper import SimConfigWrapper
from betse.util.path import dirs
from betse.util.type import types

# ....................{ CLASSES                            }....................
class SimTestConfig(object):
    '''
    Simulation configuration context encapsulating simulation configuration,
    state, and metadata.

    Simulation configuration fixtures typically return instances of this class
    as a means of communicating this context to other fixtures and tests.

    Attributes
    ----------
    config : SimConfigWrapper
        Simulation configuration wrapper wrapping the low-level dictionary
        deserialized from the YAML-formatted simulation configuration file with
        path `config_filepath`. Note the contents of this in-memory dictionary
        may be desynchronized from those of this file. For efficiency, callers
        may modify this dictionary to suite test requirements _before_
        reserializing this dictionary back to this file.
    config_filepath : py.path.local
        Absolute path of a temporary simulation configuration file specific to
        the parent fixture as a `py.path.local` instance, defining an
        object-oriented superset of the non-object-oriented `os.path` module.

    See Also
    ----------
    https://py.readthedocs.org/en/latest/path.html
        Official `py.path` class documentation.
    '''


    def __init__(self, config_filepath: 'py.path.local') -> None:
        '''
        Initialize this simulation configuration context.

        This method copies BETSE's default simulation configuration file,
        complete with all external assets (e.g., geometry masks) referenced and
        required by this file, to the passed path.

        Parameters
        ----------
        config_filepath : py.path.local
            Absolute path to which this method copies BETSE's default simulation
            configuration file. If this file already exists, an exception is
            raised.
        '''
        assert types.is_py_path_local(config_filepath), (
            types.assert_not_py_path_local(config_filepath))

        # Classify the passed parameters. While the "self.config" object
        # classified below provides this filename as a low-level string, this
        # high-level "py.path.local" instance is useful in fixtures and tests.
        self.config_filepath = config_filepath

        # Configuration deserialized from this file, reducing this filename from
        # a high-level "py.path.local" instance to a low-level string.
        self.config = SimConfigWrapper.wrap_new_default(
            filename=str(config_filepath))

        # Unconditionally disable configuration options either requiring
        # interactive input *OR* displaying interactive output for all fixtures
        # and tests.
        self.config.disable_interaction()


    def get_command_context(self) -> 'contextlib.contextmanager':
        '''
        Context manager changing the current working directory (CWD) of the
        current test to the directory containing this configuration file for the
        duration of this context.
        '''

        return dirs.current(self.config.dirname)
