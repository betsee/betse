#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level support facilities for Yet Another Markup Language (YAML), the
file format encapsulating most input and output data for this application.
'''

#FIXME: Add support for "ruamel.yaml" and "ruamel_yaml". First, note that the
#currently selected YAML implementation is given by
#"betse.metadeps.RUNTIME_MANDATORY_YAML_PROJECT_NAME". Next, note that
#supporting the new object-oriented "ruamel.yaml" API will require considerably
#different logic from that of PyYAML. To support both, consider:
#
#* Defining the following new functions:
#  * _load_pyyaml(), implementing PyYAML-specific load functionality.
#  * _load_ruamel(), implementing both "ruamel.yaml" and "ruamel_yaml"-specific
#    load functionality. For space efficiency, don't both preserving the object
#    required by this API; just discard it when the function returns.
#* Ditto for save functionality.
#* The top-level load() and save() functions should then internally call the
#  appropriate private implementations manually switched on the value of the
#  "betse.metadeps.RUNTIME_MANDATORY_YAML_PROJECT_NAME" global. *DO NOT BOTHER
#  WITH A DICTIONARY MAPPING.* We are serious. There are only two available
#  choices. Do *NOT* go overkill here, please.
#FIXME: Actually, note that there is currently no point in supporting
#"ruamel_yaml" -- which is painfully obsolete and unlikely to ever catch up.

#FIXME: Consider contributing various portions of this submodule back to
#"ruamel.yaml" -- particularly the Numpy-type-to-YAML-native-type conversions.

# ....................{ IMPORTS                            }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: To permit YAML implementations to be conditionally imported at
# application startup, no implementation (e.g., the top-level "yaml" package
# corresponding to PyYAML) are importable here.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

import numpy
from betse.util.io import iofiles
from betse.util.io.log import logs
from betse.util.path import pathnames
from betse.util.type.types import type_check, MappingOrSequenceTypes

# ....................{ GLOBALS                            }....................
FILETYPES = {'yaml', 'yml',}
'''
Set of all YAML-compliant filetypes.
'''

# ....................{ GLOBALS ~ private                  }....................
_TAG_NUMPY_NDARRAY = '!ndarray'
'''
YAML-formatted tag identifying each YAML sequence converted from a Numpy array.
'''

# ....................{ WARNERS                            }....................
@type_check
def warn_unless_filetype_yaml(filename: str) -> None:
    '''
    Log a non-fatal warning unless the passed filename has a YAML-compliant
    filetype (i.e., is suffixed by either ``.yaml`` or ``.yml``).

    Parameters
    ----------
    filename : str
        Absolute or relative path of this file.
    '''

    # Filetype of this file.
    filetype = pathnames.get_filetype_undotted_or_none(filename)

    # If this filetype is *NOT* YAML-compliant...
    if filetype not in FILETYPES:
        # If this file has *NO* filetype, log an appropriate warning.
        if filetype is None:
            logs.log_warning('YAML file "%s" has no filetype.', filename)
        # Else, this file has a filetype. Log an appropriate warning.
        else:
            logs.log_warning(
                'YAML file "%s" filetype "%s" neither "yaml" nor "yml".',
                filename, filetype)

# ....................{ LOADERS                            }....................
@type_check
def load(filename: str) -> MappingOrSequenceTypes:
    '''
    Load (i.e., read, deserialize) and return the contents of the YAML-formatted
    file with the passed path as either a dictionary or list.

    Parameters
    ----------
    filename : str
        Absolute or relative path of this file.

    Returns
    ----------
    MappingOrSequenceTypes
        Dictionary or list corresponding to the contents of this file.
    '''

    # Delay importation of the desired YAML implementation.
    import yaml

    # If this filename has no YAML-compliant filetype, log a warning.
    warn_unless_filetype_yaml(filename)

    # Load and return the contents of this YAML file.
    with iofiles.reading_chars(filename) as yaml_file:
        return yaml.load(stream=yaml_file)

# ....................{ SAVERS                             }....................
@type_check
def save(container: MappingOrSequenceTypes, filename: str) -> None:
    '''
    Save (i.e., write, serialize) the passed dictionary or list to the YAML-
    formatted file with the passed path.

    Parameters
    ----------
    container: MappingOrSequenceTypes
        Dictionary or list to be written as the contents of this file.
    filename : str
        Absolute or relative path of this file.
    '''

    # Delay importation of the desired YAML implementation.
    import yaml

    # If this filename has no YAML-compliant filetype, log a warning.
    warn_unless_filetype_yaml(filename)

    # Save this container to this YAML file.
    with iofiles.writing_chars(filename) as yaml_file:
        yaml.dump(
            data=container,
            stream=yaml_file,
            allow_unicode=True,
            default_flow_style=False,
            encoding=None,
        )

# ....................{ INITIALIZERS                       }....................
def init() -> None:
    '''
    Select and initialize the best available third-party YAML implementation
    (e.g., :mod:`ruamel.yaml`, PyYAML).

    This function selects the first such implementation importable under the
    active Python interpreter from the following list (in descending order of
    preference): :mod:`ruamel.yaml`, PyYAML.
    '''

    # Delay importation of the desired YAML implementation.
    import yaml
    from yaml.representer import SafeRepresenter

    # When saving arbitrary objects to YAML, stringify all:
    #
    # * Numpy arrays as YAML sequences.
    # * Numpy booleans as YAML booleans.
    # * Numpy complex numbers as YAML complex numbers.
    # * Numpy floats as YAML floats.
    # * Numpy integers as YAML integers.
    #
    # Since most YAML implementations lack explicit support for Numpy types,
    # attempting to saving Numpy objects to YAML results in non-human-readable
    # output defeating the whole purpose of YAML. This output resembles:
    #
    #     >>> print(yaml.dump({'a':numpy.array([3.1415,])}))
    #     /usr/lib64/python3.4/site-packages/yaml/representer.py:135:
    #     FutureWarning: comparison to `None` will result in an elementwise
    #     object comparison in the future.
    #       if data in [None, ()]:
    #     a: !!python/object/apply:numpy.core.multiarray._reconstruct
    #       args:
    #       - !!python/name:numpy.ndarray ''
    #       - !!python/tuple [0]
    #       - !!binary |
    #         Yg==
    #       state: !!python/tuple
    #       - 1
    #       - !!python/tuple [1]
    #       - !!python/object/apply:numpy.dtype
    #         args: [f8, 0, 1]
    #         state: !!python/tuple [3, <, null, null, null, -1, -1, 0]
    #       - false
    #       - !!binary |
    #         bxKDwMohCUA=
    yaml.add_representer(numpy.ndarray,  _represent_numpy_ndarray)
    yaml.add_representer(numpy.bool_,    _represent_numpy_bool)
    yaml.add_representer(numpy.complex_, _represent_numpy_complex)
    yaml.add_representer(numpy.float_,   _represent_numpy_float)
    yaml.add_representer(numpy.int_,     _represent_numpy_int)

    #FIXME: Perform this *ONLY* if the selected YAML implementation is PyYAML.
    #ruamel.yaml already fixes all Numpy-specific warnings.

    # Monkeypatch this PyYAML method to eliminate Numpy-specific warnings.
    SafeRepresenter.ignore_aliases = _ignore_aliases_monkeypatch

# ....................{ REPRESENTERS                       }....................
def _represent_numpy_ndarray(dumper, ndarray: numpy.ndarray) -> str:
    '''
    Convert the passed Numpy array into a YAML-formatted string.

    Specifically, this function returns a stringified YAML mapping whose:

    * Keys are the names of parameters accepted by the `numpy.array()` function.
    * Values are the values of these parameters required to losslessly restore
      the contents of this array on YAML deserialization of this string.

    Parameters
    ----------
    dumper: yaml.Dumper
        Object converting arbitrary Python objects to YAML-formatted strings.
    ndarray: numpy.ndarray
        Numpy array to be converted into a YAML-formatted string.

    Returns
    ----------
    str
        YAML-formatted string representing this Numpy array.
    '''

    # Dictionary mapping from numpy.array() parameter names to values.
    array_params = {
        'object': ndarray.tolist(),
        'dtype': ndarray.dtype.name,
    }

    # Convert this dictionary into a stringified YAML mapping.
    return dumper.represent_mapping(_TAG_NUMPY_NDARRAY, array_params)


def _represent_numpy_bool(dumper, npbool: numpy.bool_) -> str:
    '''
    Convert the passed Numpy boolean into a YAML-formatted string.

    Specifically, this function returns a stringified YAML boolean equivalent
    to the passed boolean.

    Parameters
    ----------
    dumper: yaml.Dumper
        Object converting arbitrary Python objects to YAML-formatted strings.
    npbool: numpy.ndarray
        Numpy boolean to be converted into a YAML-formatted string.

    Returns
    ----------
    str
        YAML-formatted string representing this Numpy bool.
    '''

    return dumper.represent_bool(bool(npbool))


def _represent_numpy_complex(dumper, npcomplex: numpy.complex_) -> str:
    '''
    Convert the passed Numpy complex number into a YAML-formatted string.

    Specifically, this function returns a stringified YAML complex number
    equivalent to the passed complex number.

    Parameters
    ----------
    dumper: yaml.Dumper
        Object converting arbitrary Python objects to YAML-formatted strings.
    npcomplex: numpy.ndarray
        Numpy complex number to be converted into a YAML-formatted string.

    Returns
    ----------
    str
        YAML-formatted string representing this Numpy complex number.
    '''

    return dumper.represent_complex(complex(npcomplex))


def _represent_numpy_float(dumper, npfloat: numpy.float_) -> str:
    '''
    Convert the passed Numpy float into a YAML-formatted string.

    Specifically, this function returns a stringified YAML float equivalent to
    the passed float.

    Parameters
    ----------
    dumper: yaml.Dumper
        Object converting arbitrary Python objects to YAML-formatted strings.
    npfloat: numpy.ndarray
        Numpy float to be converted into a YAML-formatted string.

    Returns
    ----------
    str
        YAML-formatted string representing this Numpy float.
    '''

    return dumper.represent_float(float(npfloat))


def _represent_numpy_int(dumper, npint: numpy.int_) -> str:
    '''
    Convert the passed Numpy integer into a YAML-formatted string.

    Specifically, this function returns a stringified YAML integer equivalent to
    the passed integer.

    Parameters
    ----------
    dumper: yaml.Dumper
        Object converting arbitrary Python objects to YAML-formatted strings.
    npint: numpy.ndarray
        Numpy integer to be converted into a YAML-formatted string.

    Returns
    ----------
    str
        YAML-formatted string representing this Numpy integer.
    '''

    return dumper.represent_int(int(npint))

# ....................{ MONKEYPATCHES                      }....................
def _ignore_aliases_monkeypatch(self, data):
    '''
    :meth:`yaml.representer.SafeRepresenter.ignore_aliases` method monkeypatched
    to eliminate Numpy-specific warnings.

    This method has been refactored to eliminate future warnings resembling:

        /usr/lib64/python3.4/site-packages/yaml/representer.py:135:
        FutureWarning: comparison to `None` will result in an elementwise object
        comparison in the future.
          if data in [None, ()]:
    '''

    # Refactor the default test of "if data in [None, ()]:" in the stock
    # ignore_aliases() method to avoid the future warning detailed above. To do
    # so safely, note that testing the passed "data" parameter against:
    #
    # * The singleton "None" object with the "is" operator is required.
    # * The non-singleton () object with the "==" rather than "is" operator is
    #   required. Bizarrely, despite the empty tuple being frozen, Python
    #   documentation notes that "...two occurrences of the empty tuple may or
    #   may not yield the same object". Hence, the comparison "data is ()" could
    #   unsafely return "False" even when "data" is the empty tuple!
    if data is None or data == ():
        return True
    if isinstance(data, (str, bytes, bool, int, float)):
        return True
