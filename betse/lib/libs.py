#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level **application dependency** (i.e., both mandatory and optional
third-party Python packages required by this application) facilities.

This low-level submodule defines functions intended to be called by high-level
submodules (e.g., :mod:`betse.util.cli.cliabc`) *before* attempting to import
any such dependencies.
'''

#FIXME: Refactor *ALL* functions defined by this submodule into methods of a
#new "LibResolverABC" superclass -- except the following, which are called by
#setuptools at installation-time and hence *MUST* reside outside this API for
#copy-and-paste-style reuse by downstream consumers:
#
#* get_runtime_mandatory_tuple().
#* get_runtime_optional_tuple().
#* get_testing_mandatory_tuple().
#FIXME: Ideally eliminate the boilerplate repeated across the equivalent of
#this submodule in BETSE, BETSEE, and so on by generalizing this functionality.
#To do so, we might consider:
#
#* Define a new "betse.util.py.abc.pylibabc" submodule (...or some such).
#* Define a new "LibResolverABC" abstract base class in this submodule.
#* Define a new "_metadeps" abstract property of this class, returning an
#  actual imported module of type "ModuleType".
#* For each top-level public function defined by this "nimlib" submodule,
#  define a corresponding *CONCRETE* public method of the same name in this
#  class. The implementation of each such method should defer to the
#  aforementioned "_metadeps" abstract property of this class.
#* Define a new "NimmeLibResolver(LibResolverABC)" subclass in this submodule.
#  Ideally, the only attribute defined by this subclass should be the
#  "_metadeps" attribute.
#* Instantiate a public singleton instance of this class in this submodule.
#* Rewrite all external import statements referencing this submodule to
#  instead reference that instance.
#* Remove all top-level public functions defined from this "nimlib" submodule.
#
#As a succinct example:
#
#    # In "betse.util.py.abc.pylibabc".
#    class LibResolverABC(metaclass=ABCMeta):
#        ...
#
#        @abstractproperty
#        def _metadeps(self) -> ModuleType:
#            pass
#
#
#        def die_unless_runtime_mandatory_all(self) -> None:
#            '''
#            Raise an exception unless all mandatory runtime dependencies of this
#            application are **satisfiable** (i.e., both importable and of a
#            satisfactory version) *and* all external commands required by these
#            dependencies (e.g., GraphViz's ``dot`` command) reside in the current
#            ``${PATH}``.
#
#            Raises
#            ----------
#            BetseLibException
#                If at least one mandatory runtime dependency is unsatisfiable.
#
#            See Also
#            ----------
#            :func:`betse_libs.die_unless_runtime_mandatory_all`
#                Further details.
#            '''
#
#            betse_libs.die_unless_requirements_dict(
#                self._metadeps.RUNTIME_MANDATORY)
#
#        ...
#
#
#    # In this submodule.
#    from nimme import nimmetadeps
#
#    class NimmeLibResolver(LibResolverABC):
#        @property
#        def _metadeps(self) -> ModuleType:
#            return nimmetadeps
#
#    lib_resolver = NimmeLibResolver()
#
#Ergo, given the above structure, we'd replace all import statements of the
#form "from nimme.lib import nimlib" to
#"from nimme.lib.nimlib import lib_resolver".

# ....................{ IMPORTS                           }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: To raise human-readable exceptions on missing mandatory
# dependencies, the top-level of this module may import *ONLY* from packages
# guaranteed to exist at initial runtime (i.e., standard Python and application
# packages).
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

from betse import metadeps
from betse.exceptions import BetseLibException
# from betse.util.io.log import logs
from betse.util.py.module import pymodname
from betse.util.type.types import (
    type_check, MappingType, ModuleOrSequenceTypes)

# ....................{ EXCEPTIONS                        }....................
def die_unless_runtime_mandatory_all() -> None:
    '''
    Raise an exception unless all mandatory runtime dependencies of this
    application are **satisfiable** (i.e., both importable and of a
    satisfactory version) *and* all external commands required by these
    dependencies (e.g., GraphViz's ``dot`` command) reside in the current
    ``${PATH}``.

    Equivalently, this function raises an exception if at least one such
    dependency is unsatisfied. For importable unsatisfied dependencies with
    :mod:`setuptools`-specific metadata (e.g., ``.egg-info/``-suffixed
    subdirectories of the ``site-packages/`` directory for the active Python 3
    interpreter, typically created by :mod:`setuptools` at install time), this
    function additionally validates the versions of these dependencies to
    satisfy all application requirements.

    Raises
    ----------
    BetseLibException
        If at least one mandatory runtime dependency is unsatisfiable.
    '''

    die_unless_requirements_dict(metadeps.RUNTIME_MANDATORY)


@type_check
def die_unless_runtime_optional(*requirement_names: str) -> None:
    '''
    Raise an exception unless all optional runtime dependencies of this
    application with the passed :mod:`setuptools`-specific project names are
    **satisfiable** (i.e., both importable and of a satisfactory version)
    *and* all external commands required by these dependencies (e.g.,
    GraphViz's ``dot`` command) reside in the current ``${PATH}``.

    Parameters
    ----------
    requirement_names : Tuple[str]
        Tuple of the names of all :mod:`setuptools`-specific projects
        implementing these dependencies (e.g., ``NetworkX``). If any such name
        is unrecognized (i.e., is *not* a key of the
        :data:`metadeps.RUNTIME_OPTIONAL` dictionary), an exception is raised.

    Raises
    ----------
    BetseLibException
        If at least one such dependency is unsatisfiable.

    See Also
    ----------
    :func:`die_unless_runtime_mandatory_all`
        Further details.
    '''

    die_unless_requirements_dict_keys(
        metadeps.RUNTIME_OPTIONAL, *requirement_names)

# ....................{ EXCEPTIONS ~ dict                 }....................
@type_check
def die_unless_requirements_dict(requirements_dict: MappingType) -> None:
    '''
    Raise an exception unless all dependencies described by the passed
    dictionary are **satisfiable** (i.e., both importable and of a satisfactory
    version) *and* all external commands required by these dependencies (e.g.,
    GraphViz's ``dot`` command) reside in the current ``${PATH}``.

    Raises
    ----------
    BetseLibException
        If at least passed dependency is unsatisfiable.
    '''

    # If the "pkg_resources" setuptools dependency is missing, raise an
    # exception *BEFORE* importing this dependency below.
    pymodname.die_unless_module(
        module_name='pkg_resources',
        exception_message='Mandatory dependency "pkg_resources" not found.',
    )

    # Import the following submodule, which globally imports "pkg_resources",
    # *AFTER* validating "pkg_resources" to be importable.
    from betse.lib.setuptools import setuptool

    # Validate these dependencies.
    setuptool.die_unless_requirements_dict(requirements_dict)

    # Validate all external commands required by these dependencies.
    die_unless_command(*requirements_dict.keys())


@type_check
def die_unless_requirements_dict_keys(
    requirements_dict: MappingType, *requirement_names: str) -> None:
    '''
    Raise an exception unless all dependencies with the passed
    :mod:`setuptools`-specific project names described by the passed dictionary
    are **satisfiable** (i.e., both importable and of a satisfactory version)
    *and* all external commands required by these dependencies (e.g.,
    GraphViz's ``dot`` command) reside in the current ``${PATH}``.

    Parameters
    ----------
    requirements_dict : MappingType
        Dictionary mapping from the names of all :mod:`setuptools`-specific
        projects implementing these dependencies to the requirements strings
        constraining these dependencies.
    requirement_names : tuple[str]
        Tuple of the names of all :mod:`setuptools`-specific projects
        implementing these dependencies (e.g., ``NetworkX``). If any such
        name is *not* a key of this dictionary, an exception is raised.

    Raises
    ----------
    BetseLibException
        If at least passed dependency is unsatisfiable.
    '''

    # Avoid circular import dependencies.
    from betse.lib.setuptools import setuptool

    # Validate these dependencies.
    setuptool.die_unless_requirements_dict_keys(
        requirements_dict, *requirement_names)

    # Validate all external commands required by these dependencies.
    die_unless_command(*requirement_names)

# ....................{ EXCEPTIONS ~ commands             }....................
@type_check
def die_unless_command(*requirement_names: str) -> None:
    '''
    Raise an exception unless all external commands required by all application
    dependencies (of any type, including optional, mandatory, runtime, testing,
    or otherwise) with the passed :mod:`setuptools`-specific project names are
    **installed** (i.e., are executable files in the current ``${PATH}``).

    Parameters
    ----------
    requirement_names : Tuple[str]
        Tuple of the names of all :mod:`setuptools`-specific projects
        corresponding to these dependencies (e.g., ``NetworkX``).

    Raises
    ----------
    BetseLibException
        If any external command required by any such dependency is *not* found.
    '''

    # Avoid circular import dependencies.
    from betse.util.path.command import cmdpath

    # If any external command required by any such dependency is *NOT* found,
    # iteratively search for the first such missing command and raise a
    # human-readable exception synopsizing this command. For efficiency, this
    # inefficient iteration is performed *ONLY* as required.
    if not is_command(*requirement_names):
        # For the name of each such dependency...
        for requirement_name in requirement_names:
            # For each "betse.metadata.DependencyCommand" instance
            # describing each external command required by this dependency if
            # any *OR* the empty tuple otherwise...
            for dependency_command in metadeps.EXTERNAL_COMMANDS.get(
                requirement_name, ()):
                # If this command is *NOT* in the ${PATH}, raise an exception.
                if not cmdpath.is_pathable(dependency_command.basename):
                    raise BetseLibException(
                        'Dependency "{}" unsatisfied, as '
                        '{} not installed '
                        '(i.e., command "{}" not found in ${{PATH}}).'.format(
                            requirement_name,
                            dependency_command.name,
                            dependency_command.basename,
                        ))

# ....................{ TESTERS                           }....................
@type_check
def is_command(*requirement_names: str) -> bool:
    '''
    ``True`` only if all external commands required by all application
    dependencies (of any type, including optional, mandatory, runtime, testing,
    or otherwise) with the passed :mod:`setuptools`-specific project names are
    **installed** (i.e., are executable files in the current ``${PATH}``).

    Parameters
    ----------
    requirement_names : Tuple[str]
        Tuple of the names of all :mod:`setuptools`-specific projects
        corresponding to these dependencies (e.g., ``NetworkX``).
    '''

    # Avoid circular import dependencies.
    from betse.util.path.command import cmdpath

    # Return True only if...
    return all(
        # Each external command required by each dependency is in the ${PATH}.
        cmdpath.is_pathable(dependency_command.basename)
        # For the name of each passed dependency...
        for requirement_name in requirement_names
        # For the tuple of all "betse.metadata.DependencyCommand" instances
        # describing all external commands required by this dependency if any
        # *OR* the empty tuple otherwise...
        for dependency_command in metadeps.EXTERNAL_COMMANDS.get(
            requirement_name, ())
    )


@type_check
def is_runtime_optional(*requirement_names: str) -> bool:
    '''
    ``True`` only if all optional runtime dependencies of this application with
    the passed :mod:`setuptools`-specific project names are **satisfiable**
    (i.e., both importable and of a satisfactory version) *and* all external
    commands required by these dependencies (e.g., GraphViz's ``dot`` command)
    reside in the current ``${PATH}``.

    Parameters
    ----------
    requirement_names : tuple[str]
        Tuple of the names of all :mod:`setuptools`-specific projects
        implementing these dependencies (e.g., ``NetworkX``). If any such
        name is *not* a key of the :data:`metadeps.RUNTIME_OPTIONAL` dictionary
        and is thus unrecognized, an exception is raised.

    Returns
    ----------
    bool
        ``True`` only if these optional runtime dependencies are satisfiable.

    See Also
    ----------
    :func:`die_unless_runtime_mandatory_all`
        Further details.
    '''

    return is_requirements_dict_keys(
        metadeps.RUNTIME_OPTIONAL, *requirement_names)


@type_check
def is_requirements_dict_keys(
    requirements_dict: MappingType, *requirement_names: str) -> bool:
    '''
    ``True`` only if all dependencies with the passed
    :mod:`setuptools`-specific project names are **satisfiable** (i.e., both
    importable and of a satisfactory version) *and* all external commands
    required by these dependencies (e.g., GraphViz's ``dot`` command) reside in
    the current ``${PATH}``.

    Parameters
    ----------
    requirements_dict : MappingType
        Dictionary mapping from the names of all :mod:`setuptools`-specific
        projects implementing these dependencies to the requirements strings
        constraining these dependencies.
    requirement_names : tuple[str]
        Tuple of the names of all :mod:`setuptools`-specific projects
        implementing these dependencies (e.g., ``NetworkX``). If any such
        name is *not* a key of this dictionary, an exception is raised.

    Returns
    ----------
    bool
        ``True`` only if these dependencies are all satisfiable.

    See Also
    ----------
    :func:`die_unless_runtime_mandatory_all`
        Further details.
    '''

    # Avoid circular import dependencies.
    from betse.lib.setuptools import setuptool

    # Tuple of requirements strings converted from this subset of key-value
    # pairs of these requirements.
    requirements_tuple = setuptool.get_requirements_str_from_dict_keys(
        requirements_dict, *requirement_names)

    # Return True only if...
    return (
        # These dependencies are all satisfied.
        setuptool.is_requirement_str(*requirements_tuple) and
        # All external commands required by these dependencies are installed.
        is_command(*requirement_names)
    )

# ....................{ GETTERS ~ runtime                 }....................
def get_runtime_mandatory_tuple() -> tuple:
    '''
    Tuple listing the :mod:`setuptools`-specific requirement string containing
    the mandatory name and optional version and extras constraints of each
    mandatory runtime dependency for this application, dynamically converted
    from the :data:`metadata.RUNTIME_MANDATORY` dictionary.
    '''

    # Avoid circular import dependencies.
    from betse.lib.setuptools import setuptool

    # Convert this dictionary into a tuple.
    return setuptool.get_requirements_str_from_dict(
        metadeps.RUNTIME_MANDATORY)


def get_runtime_optional_tuple() -> tuple:
    '''
    Tuple listing the :mod:`setuptools`-specific requirement string containing
    the mandatory name and optional version and extras constraints of each
    optional runtime dependency for this application, dynamically converted
    from the :data:`metadata.RUNTIME_OPTIONAL` dictionary.
    '''

    # Avoid circular import dependencies.
    from betse.lib.setuptools import setuptool

    # Convert this dictionary into a tuple.
    return setuptool.get_requirements_str_from_dict(
        metadeps.RUNTIME_OPTIONAL)

# ....................{ GETTERS ~ testing                 }....................
def get_testing_mandatory_tuple() -> tuple:
    '''
    Tuple listing the :mod:`setuptools`-specific requirement string containing
    the mandatory name and optional version and extras constraints of each
    mandatory testing dependency for this application, dynamically converted
    from the :data:`metadata.RUNTIME_OPTIONAL` dictionary.
    '''

    # Avoid circular import dependencies.
    from betse.lib.setuptools import setuptool

    # Convert this dictionary into a tuple.
    return setuptool.get_requirements_str_from_dict(
        metadeps.TESTING_MANDATORY)

# ....................{ GETTERS ~ metadata                }....................
def get_metadatas() -> tuple:
    '''
    Tuple of 2-tuples `(metedata_name, metadata_value`), describing all
    currently installed optional and mandatory third-party dependencies
    required at both runtime and testing time.
    '''

    # Defer heavyweight imports.
    from betse.lib.matplotlib.matplotlibs import mpl_config
    from betse.lib.numpy import numpys
    from betse.lib.setuptools import setuptool

    # Tuple of all dependency versions.
    LIB_VERSION_METADATA = (
        # Dependencies metadata.
        ('runtime dependencies (mandatory)',
         setuptool.get_requirements_dict_synopsis(metadeps.RUNTIME_MANDATORY)),
        ('runtime dependencies (optional)',
         setuptool.get_requirements_dict_synopsis(metadeps.RUNTIME_OPTIONAL)),
        ('testing dependencies (mandatory)',
         setuptool.get_requirements_dict_synopsis(metadeps.TESTING_MANDATORY)),
    )

    # Tuple of all dependency-specific metadata.
    LIB_METADATA = (
        # matplotlib metadata.
        ('matplotlib', mpl_config.get_metadata()),
    )

    # Return a tuple aggregating these tuples.
    return (
        LIB_VERSION_METADATA +
        LIB_METADATA +
        numpys.get_metadatas()
    )

# ....................{ IMPORTERS                         }....................
@type_check
def import_runtime_optional(*requirements_name: str) -> ModuleOrSequenceTypes:
    '''
    Import and return either a single top-level module object if passed one
    :mod:`setuptools`-specific requirement name *or* sequence of top-level
    module objects otherwise (i.e., if passed multiple such names).

    Each returned module object is guaranteed to satisfy the optional runtime
    dependency of this application with the corresponding requirement name.

    Parameters
    ----------
    requirements_name : tuple[str]
        Tuple of the names of all :mod:`setuptools`-specific projects
        implementing these dependencies (e.g., ``NetworkX``). If any such name
        is unrecognized (i.e., is *not* a key of the
        :data:`metadeps.RUNTIME_OPTIONAL` dictionary), an exception is raised.
    '''

    return import_requirements_dict_keys(
        metadeps.RUNTIME_OPTIONAL, *requirements_name)


@type_check
def import_requirements_dict_keys(
    requirements_dict: MappingType, *requirements_name: str) -> (
    ModuleOrSequenceTypes):
    '''
    Import and return either a single top-level module object if passed one
    :mod:`setuptools`-specific requirement name *or* sequence of top-level
    module objects otherwise (i.e., if passed multiple such names).

    Each returned module object is guaranteed to satisfy the dependency of this
    application with the corresponding requirement name as a key of the passed
    dictionary.
    '''

    # Avoid circular import dependencies.
    from betse.lib.setuptools import setuptool

    # Validate these dependencies.
    die_unless_requirements_dict_keys(requirements_dict, *requirements_name)

    # Validate all external commands required by these dependencies.
    return setuptool.import_requirements_dict_keys(
        requirements_dict, *requirements_name)
