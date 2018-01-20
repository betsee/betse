#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level **application dependency** (i.e., both mandatory and optional
third-party Python packages required by this application) facilities.

This low-level submodule defines functions intended to be called by high-level
submodules (e.g., :mod:`betse.util.cli.cliabc`) *before* attempting to import any
such dependencies.
'''

# ....................{ IMPORTS                            }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: To raise human-readable exceptions on missing mandatory dependencies,
# the top-level of this module may import *ONLY* from packages guaranteed to
# exist at initial runtime (i.e., standard Python and application packages).
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

from betse import metadata, metadeps
from betse.exceptions import BetseLibException
from betse.util.io.log import logs
from betse.util.type.types import type_check, MappingType, StrOrNoneTypes

# ....................{ GLOBALS                            }....................
_IS_INITTED = False
'''
``True`` only if the :func:`init` function has already been called.

That function uses this private boolean to guard against repeated invocations of
the :func:`init` function from multiple modules in the same Python process
(e.g., :mod:`betse.science.__init__`, :mod:`betse.util.cli.cliabc`). While that
function does technically support repeated calls, each additional call after the
first inefficiently performs no meaningful work and is thus safely ignorable.
'''

# ....................{ EXCEPTIONS                         }....................
def die_unless_runtime_mandatory_all() -> None:
    '''
    Raise an exception unless all mandatory runtime dependencies of this
    application are **satisfiable** (i.e., both importable and of a satisfactory
    version) *and* all external commands required by these dependencies (e.g.,
    GraphViz's ``dot`` command) reside in the current ``${PATH}``.

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

# ....................{ EXCEPTIONS ~ dict                  }....................
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

    # Avoid circular import dependencies.
    from betse.util.type import modules

    # If the "pkg_resources" setuptools dependency is missing, raise an
    # exception *BEFORE* importing this dependency below.
    modules.die_unless_module(
        module_name='pkg_resources',
        exception_message='Mandatory dependency "pkg_resources" not found.',
    )

    # Import this submodule *AFTER* validating "pkg_resources" to exist above.
    from betse.lib.setuptools import setuptool

    # Validate these dependencies.
    setuptool.die_unless_requirements_dict(requirements_dict)

    # Validate all external commands required by these dependencies.
    die_unless_commands(*requirements_dict.keys())


@type_check
def die_unless_requirements_dict_keys(
    requirements_dict: MappingType, *requirement_names: str) -> None:
    '''
    Raise an exception unless all dependencies with the passed
    :mod:`setuptools`-specific project names described by the passed dictionary
    are **satisfiable** (i.e., both importable and of a satisfactory version)
    *and* all external commands required by these dependencies (e.g., GraphViz's
    ``dot`` command) reside in the current ``${PATH}``.

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
    die_unless_commands(*requirement_names)

# ....................{ EXCEPTIONS ~ commands              }....................
@type_check
def die_unless_commands(*requirement_names: str) -> None:
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
    if not is_commands(*requirement_names):
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

# ....................{ TESTERS                            }....................
@type_check
def is_commands(*requirement_names: str) -> bool:
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
    ``True`` only if all dependencies with the passed :mod:`setuptools`-specific
    project names are **satisfiable** (i.e., both importable and of a
    satisfactory version) *and* all external commands required by these
    dependencies (e.g., GraphViz's ``dot`` command) reside in the current
    ``${PATH}``.

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
    requirements_tuple = setuptool.convert_requirements_dict_keys_to_tuple(
        requirements_dict, *requirement_names)

    # Return True only if...
    return (
        # These dependencies are all satisfied.
        setuptool.is_requirement_str(*requirements_tuple) and
        # All external commands required by these dependencies are installed.
        is_commands(*requirement_names)
    )

# ....................{ GETTERS ~ runtime                  }....................
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
    return setuptool.convert_requirements_dict_to_tuple(
        metadeps.RUNTIME_MANDATORY)


def get_runtime_optional_tuple() -> tuple:
    '''
    Tuple listing the :mod:`setuptools`-specific requirement string containing
    the mandatory name and optional version and extras constraints of each
    optional runtime dependency for this application, dynamically converted from
    the :data:`metadata.RUNTIME_OPTIONAL` dictionary.
    '''

    # Avoid circular import dependencies.
    from betse.lib.setuptools import setuptool

    # Convert this dictionary into a tuple.
    return setuptool.convert_requirements_dict_to_tuple(
        metadeps.RUNTIME_OPTIONAL)

# ....................{ GETTERS ~ testing                  }....................
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
    return setuptool.convert_requirements_dict_to_tuple(
        metadeps.TESTING_MANDATORY)

# ....................{ GETTERS ~ metadata                 }....................
def get_metadatas() -> tuple:
    '''
    Tuple of 2-tuples `(metedata_name, metadata_value`), describing all
    currently installed optional and mandatory third-party dependencies required
    at both runtime and testing time.
    '''

    # Defer heavyweight imports.
    from betse.lib.matplotlib.matplotlibs import mpl_config
    from betse.lib.numpy import numpys
    from betse.lib.setuptools import setuptool

    # Tuple of all dependency versions.
    LIB_VERSION_METADATA = (
        # Dependencies metadata.
        ('runtime dependencies (mandatory)',
         setuptool.get_requirements_dict_metadata(metadeps.RUNTIME_MANDATORY)),
        ('runtime dependencies (optional)',
         setuptool.get_requirements_dict_metadata(metadeps.RUNTIME_OPTIONAL)),
        ('testing dependencies (mandatory)',
         setuptool.get_requirements_dict_metadata(metadeps.TESTING_MANDATORY)),
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

# ....................{ IMPORTERS                          }....................
@type_check
def import_runtime_optional(*requirement_names: str) -> object:
    '''
    Import and return the top-level module object satisfying each optional
    runtime dependency of this application with the passed name.

    Parameters
    ----------
    requirement_names : tuple[str]
        Tuple of the names of all :mod:`setuptools`-specific projects
        implementing these dependencies (e.g., ``NetworkX``). If any such name
        is unrecognized (i.e., is *not* a key of the
        :data:`metadeps.RUNTIME_OPTIONAL` dictionary), an exception is raised.

    See Also
    ----------
    :func:`import_requirements_dict_keys`
        Further details.
    '''

    return import_requirements_dict_keys(
        metadeps.RUNTIME_OPTIONAL, *requirement_names)


@type_check
def import_requirements_dict_keys(
    requirements_dict: MappingType, *requirement_names: str) -> object:
    '''
    Import and return the top-level module object satisfying each dependency
    with the passed name described by the passed dictionary.

    See Also
    ----------
    :func:`setuptool.import_requirements_dict_keys`
        Further details.
    '''

    # Avoid circular import dependencies.
    from betse.lib.setuptools import setuptool

    # Validate these dependencies.
    die_unless_requirements_dict_keys(requirements_dict, *requirement_names)

    # Validate all external commands required by these dependencies.
    return setuptool.import_requirements_dict_keys(
        requirements_dict, *requirement_names)

# ....................{ INITIALIZERS                       }....................
def reinit(*args, **kwargs) -> None:
    '''
    (Re-)initialize all mandatory runtime dependencies of this application with
    the passed parameters.

    Parameters
    ----------
    All passed parameters are passed to the :func:`init` function as is.
    '''

    # Force the init() function to reinitialize this application.
    global _IS_INITTED
    _IS_INITTED = False

    # Reinitialize these dependencies with these parameters.
    init(*args, **kwargs)


@type_check
def init(matplotlib_backend_name: StrOrNoneTypes = None) -> None:
    '''
    Initialize all mandatory runtime dependencies of this application.

    Specifically, this function:

    * Reconfigures matplotlib with sane defaults specific to the current
      platform.

    Parameters
    ----------
    matplotlib_backend_name : StrOrNoneTypes
        Name of the matplotlib backend to explicitly enable. Defaults to
        ``None``, in which case this method implicitly enables the first
        importable backend known to be both usable and supported by this
        application (in descending order of preference).
    '''

    # If this function has already been called, noop.
    global _IS_INITTED
    if     _IS_INITTED:
        return

    # Defer heavyweight imports.
    from betse.lib.matplotlib.matplotlibs import mpl_config
    from betse.lib.numpy import numpys
    from betse.lib.pil import pils
    from betse.lib.yaml import yamls

    # Log this initialization. Since initializing heavyweight third-party
    # dependencies (especially matplotlib) consumes non-trivial time, this
    # message is intentionally exposed to all users by default.
    logs.log_info('Loading third-party %s dependencies...', metadata.NAME)

    # Initialize these dependencies.
    mpl_config.init(backend_name=matplotlib_backend_name)
    numpys.init()
    pils.init()
    yamls.init()

    # Record this function as having been called *AFTER* successfully doing so.
    _IS_INITTED = True
