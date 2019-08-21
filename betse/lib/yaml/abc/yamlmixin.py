#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Abstract mixins of YAML-backed configuration subclasses, standardizing property
nomenclature for common YAML design patterns.
'''

# ....................{ IMPORTS                           }....................
from betse.lib.yaml.yamlalias import yaml_alias
# from betse.util.type.types import type_check, ClassType, SequenceTypes

# ....................{ MIXINS ~ alias                    }....................
# Mixins standardizing various general-purpose YAML aliases.

class YamlBooledMixin(object):
    '''
    Mixin of all **YAML-backed booled configuration** (i.e., backed by a YAML
    dictionary with top-level key ``enabled`` whose value is a boolean
    specifying whether this configuration is enabled or disabled) subclasses.

    This class is suitable for use as a multiple-inheritance mixin. To preserve
    the expected method resolution order (MRO) semantics, this class should
    typically be inherited *first* rather than *last* in subclasses.

    Attributes
    ----------
    is_enabled : bool
        ``True`` only if this list item is enabled.
    '''

    # ..................{ ALIASES                           }..................
    is_enabled = yaml_alias("['enabled']", bool)


class YamlNamedMixin(object):
    '''
    Mixin of all **YAML-backed named configuration** (i.e., backed by a YAML
    dictionary with top-level key ``name`` whose value is a human- and/or
    machine-readable string identifying this configuration's presumably unique
    name) subclasses.

    This class is suitable for use as a multiple-inheritance mixin. To preserve
    the expected method resolution order (MRO) semantics, this class should
    typically be inherited *first* rather than *last* in subclasses.

    Attributes
    ----------
    name : str
        Arbitrary string typically uniquely identifying this configuration
        (e.g., ``spot``, naming a spatially isolated tissue profile).
    '''

    # ..................{ ALIASES                           }..................
    name = yaml_alias("['name']", str)


class YamlTypedMixin(object):
    '''
    Mixin of all **YAML-backed typed configuration** (i.e., backed by a YAML
    dictionary with top-level key ``type`` whose value is a machine-readable
    string identifying this configuration's type) subclasses.

    This class is suitable for use as a multiple-inheritance mixin. To preserve
    the expected method resolution order (MRO) semantics, this class should
    typically be inherited *first* rather than *last* in subclasses.

    Attributes
    ----------
    kind : str
        Lowercase alphanumeric string uniquely identifying the type of this
        configuration (e.g., ``voltage_membrane``, identifying a transmembrane
        voltage). See each ``type`` key of the corresponding list in the
        default simulation configuration file for real-world examples.
    '''

    # ..................{ ALIASES                           }..................
    kind = yaml_alias("['type']", str)
