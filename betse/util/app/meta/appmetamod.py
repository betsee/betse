#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level **application dependency metadata** (i.e., lists of version-pinned
dependencies synopsizing application requirements) functionality.
'''

# ....................{ IMPORTS                           }....................
# from betse.util.io.log import logs
from betse.util.type.types import (
    type_check, MappingType, ModuleType, IterableTypes)

# ....................{ GLOBALS                           }....................
MERGE_MODULE_METADEPS_DICTS_NAME = (
    'RUNTIME_MANDATORY',
    'RUNTIME_OPTIONAL',
    'TESTING_MANDATORY',
    'REQUIREMENT_NAME_TO_COMMANDS',
)
'''
Tuple of the names of all global dictionaries required to be defined by each
application dependency metadata module passed as input to the
:func:`merge_module_metadeps` function via the ``modules_metadeps`` parameter.
'''

# ....................{ MAKERS                            }....................
#FIXME: Improve documentation with respect to key collisions. Specifically,
#provide concrete examples without code (which would probably be a bit too
#much, frankly) of how this function resolves the following similar edge cases:
#
#* Two or more dictionaries to be merged contain key-value pairs containing the
#  same key but differing values.
#* Two or more dictionaries to be merged contain the same key-value pairs.

@type_check
def merge_module_metadeps(
    module_name: str, modules_metadeps: IterableTypes) -> ModuleType:
    '''
    Dynamically create and return a new **application dependency metadata
    module** (i.e., module satisfying the informal protocol defined by the
    :mod:`betse.metadeps` module) with the passed name, iteratively merging the
    contents of all passed application dependency metadata modules.

    Specifically, this function (in order):

    #. Dynamically creates a new module with the passed name.
    #. For each of the global attributes assumed by the application dependency
       metadata protocol (i.e., dictionaries named ``RUNTIME_MANDATORY``,
       ``RUNTIME_OPTIONAL``, ``TESTING_MANDATORY``, and
       ``REQUIREMENT_NAME_TO_COMMANDS``), creates a new global attribute of the
       same name in this new module whose value is a dictionary merging the
       contents of all dictionaries of the same name defined by all passed
       modules. For example, if passed the :mod:`betse.metadeps` and
       :mod:`betsee.guimetadeps` modules, this function creates a new module
       defining a global dictionary named ``RUNTIME_MANDATORY`` merging the
       contents of the :attr:`betse.metadeps.RUNTIME_MANDATORY` and
       :mod:`betsee.guimetadeps.RUNTIME_MANDATORY` dictionaries.
    #. Returns this new module.

    Motivation
    ----------
    This utility function is *only* intended to be called by subclass
    implementations of the abstract private
    :meth:`betse.util.app.meta.appmetaabc.AppMetaABC._module_metadeps` property
    in downstream consumers. Notably, BETSEE calls this function from its
    concrete implementation of this property to merge its own dependency
    requirements with those of BETSE; doing so guarantees that calls to
    dependency-based functions within the BETSE codebase (e.g., of the
    :func:`betse.lib.libs.import_runtime_optional` function dynamically
    importing one or more optional dependencies) behave as expected.

    Caveats
    ----------
    **Order is insignificant.** If any of the requisite global dictionaries
    defined by any of the passed modules contain one or more key-value pairs
    contained in any other dictionary of the same name defined by any other
    passed module, this function raises an exception. Since this prevents any
    collisions between key-value pairs, *no* implicit precedence exists between
    these modules.

    Parameters
    ----------
    module_name : str
        Fully-qualified name of the module to be created.
    modules_metadeps : IterableTypes
        Iterable of all application dependency metadata modules to be merged.
        Order is insignificant. See above.

    Raises
    ----------
    BetseAttrException
        If any of these modules fail to define a requisite attribute.
    BetseMappingException
        If any two global dictionaries of the same name defined by any two of
        these modules **item-collide** (i.e., if any key-value pair in any such
        dictionary is also a key-value pair in any other such dictionary).
    BetseModuleException
        If a module with this target module name already exists.
    BetseTypeException
        If any of the requisite attributes defined by any of these modules are
        *not* dictionaries.

    Returns
    ----------
    ModuleType
        Output application dependency metadata module with this name, merging
        the contents of these input application dependency metadata modules.
    '''

    # Avoid circular import dependencies.
    from betse.util.py.module import pymodname, pymodule
    from betse.util.type.iterable import itertest, sequences
    from betse.util.type.iterable.mapping import mapmerge
    from betse.util.type.obj import objects
    from betse.util.type.text.string import strjoin

    # Sequence of modules converted from this iterable of modules.
    modules_metadeps = sequences.to_sequence(modules_metadeps)

    # If less than two modules were passed, raise an exception.
    sequences.die_if_length_less_than(sequence=modules_metadeps, length=2)
    # Else, at least two modules were passed.

    # If any of the passed modules is *NOT* a module, raise an exception.
    itertest.die_unless_items_instance_of(
        iterable=modules_metadeps, cls=ModuleType)

    # Dictionary mapping from the name to value of each module-scoped
    # attribute to be declared in the module to be created and returned,
    # defaulting to the "RequirementCommand" class globally defined by the
    # first such module.
    trg_module_attr_name_to_value = {
        'RequirementCommand': objects.get_attr(
            obj=modules_metadeps[0], attr_name='RequirementCommand'),
    }

    # For the name of each such global dictionary...
    for module_dict_name in MERGE_MODULE_METADEPS_DICTS_NAME:
        # Generator comprehension aggregating all of the global dictionaries
        # defined by all of these input modules, raising exceptions if any such
        # module fails to define such a dictionary.
        src_modules_dict = (
            objects.get_attr(
                obj=modules_metadep,
                attr_name=module_dict_name,
                attr_type=MappingType,
            )
            for modules_metadep in modules_metadeps
        )

        # Merge these dictionaries into the dictionary to be returned, raising
        # exceptions if any requirement defined by any such dictionary
        # item-collides (i.e., if any key-value pair in any such dictionary is
        # also a key-value pair in any other such dictionary).
        trg_module_attr_name_to_value[module_dict_name] = mapmerge.merge_maps(
            src_modules_dict)

    # Double-quoted conjunction of the fully-qualified names of all passed
    # input application dependency metadata modules.
    src_modules_metadeps_name = (
        strjoin.join_iterable_as_conjunction_double_quoted(
            pymodule.get_name_qualified(module_metadep)
            for module_metadep in modules_metadeps
        ))

    # Output application dependency metadata module to be returned.
    trg_module_metadeps = pymodname.make_module(
        module_name=module_name,
        module_attr_name_to_value=trg_module_attr_name_to_value,
        module_doc='''
**Application dependency metadata** (i.e., sequences of version-pinned
dependencies synopsizing application requirements), dynamically merged from the
{} modules.'''.format(src_modules_metadeps_name),
    )

    # Return this module.
    return trg_module_metadeps
