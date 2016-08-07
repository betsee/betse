#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level support facilities for Yet Another Markup Language (YAML), the
file format encapsulating most input and output data for this application.
'''

# ....................{ IMPORTS                            }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: To permit YAML implementations to be conditionally imported at
# application startup, no implementation (e.g., the top-level "yaml" package,
# corresponding to PyYAML) may be imported here.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

from betse.util.type.types import type_check, MappingType, SequenceTypes

# ....................{ CONSTANTS                          }....................
#FIXME: Initialize and leverage this below.
_YAML_PACKAGE_NAME = None
'''
'''

# ....................{ SAVERS                             }....................
@type_check
def save(container: (MappingType,) + SequenceTypes, filename: str) -> None:
    '''
    Save (i.e., write, serialize) the passed dictionary or list to the YAML-
    formatted file with the passed path.

    Parameters
    ----------
    container: (MappingType,) + SequenceTypes
        Dictionary or list to be written as the contents of this file.
    filename : str
        Absolute or relative path of this file.
    '''

    # Delay importation of the desired YAML implementation.
    import yaml

    # Save this container to this YAML file.
    with open(filename, 'w') as yaml_file:
        yaml.dump(container, yaml_file)

# ....................{ LOADERS                            }....................
@type_check
def load(filename: str) -> (MappingType,) + SequenceTypes:
    '''
    Load (i.e., read, deserialize) and return the contents of the YAML-formatted
    file with the passed path as either a dictionary or list.

    Parameters
    ----------
    filename : str
        Absolute or relative path of this file.

    Returns
    ----------
    (MappingType,) + SequenceTypes
        Dictionary or list corresponding to the contents of this file.
    '''

    # Delay importation of the desired YAML implementation.
    import yaml

    # Load and return the contents of this YAML file.
    with open(filename, 'r') as yaml_file:
        return yaml.load(yaml_file)