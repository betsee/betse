#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level **resource** (i.e., arbitrary file or directory contained within an
arbitrary package, regardless of whether that file or directory physically
exists *or* is compressed within an EGG-like archive file and therefore only
exists in the abstract) facilities.

See Also
----------
https://setuptools.readthedocs.io/en/latest/pkg_resources.html#resourcemanager-api
    Official :mod:setuptools documentation for the "ResourceManager API."
'''

# ....................{ IMPORTS                            }....................
import pkg_resources
from betse.util.type.types import type_check
from pkg_resources import Requirement

# ....................{ TESTERS                            }....................
@type_check
def is_dir(module_name: (str, Requirement), dirname: str) -> bool:
    '''
    ``True`` only if the resource whose pathname is the concatenation of the
    following strings is an existing directory (*in order*):

    #. The absolute pathname of the top-level directory containing the module or
       requirement with the passed name (e.g.,
       ``/home/hiroprotagonist/betse/betse/``).
    #. The directory separator for the current platform (e.g., ``/`` on
       POSIX-compatible platforms).
    #. The passed relative pathname (e.g., ``data/yaml/sim_config.yaml``).

    Parameters
    ----------
    module_name : str, Requirement
        Either the fully-qualified name of the module *or* the
        :mod:`setuptools`- specific requirement this resource is relative to.
    dirname : str
        ``/``-separated pathname relative to this module or requirement. This
        pathname *must* be POSIX-compliant and hence separated with the POSIX-
        specific ``/`` directory separator rather than the Windows-specific ``\``
        directory separator (e.g., by calling ``'\\'.join()`` rather than either
        :func:`betse.util.path.paths.join` or :func:`os.path.join`).

    Returns
    ----------
    bool
        ``True`` only if this resource is a directory.
    '''

    return pkg_resources.resource_isdir(module_name, dirname)

# ....................{ GETTERS                            }....................
@type_check
def get_pathname(module_name: (str, Requirement), pathname: str) -> str:
    '''
    Absolute path of the resource whose path the pathname of the module or
    requirement with passed name joined with the passed relative pathname is an
    existing directory.

    Caveats
    ----------
    If this resource resides in an EGG-like archive file, this function silently
    extracts this resource and all resources transitively required by this
    resource to a temporary :mod:`setuptools`-specific directory *before*
    returning the absolute path of this resource within this directory.

    To quote the `ResourceManager API`_ documentation:

         If the resource is in an archive distribution (such as a zipped egg),
         it will be extracted to a cache directory, and the filename within the
         cache will be returned. If the named resource is a directory, then all
         resources within that directory (including subdirectories) are also
         extracted. If the named resource is a C extension or “eager resource”
         (see the setuptools documentation for details), then all C extensions
         and eager resources are extracted at the same time.

    :: _`ResourceManager API`:
       https://setuptools.readthedocs.io/en/latest/pkg_resources.html#resourcemanager-api

    Parameters
    ----------
    module_name : str, Requirement
        Either the fully-qualified name of the module *or* the
        :mod:`setuptools`-specific requirement this resource is relative to.
    pathname : str
        ``/``-separated pathname relative to this module or requirement. This
        pathname *must* be POSIX-compliant and hence separated with the POSIX-
        specific ``/`` directory separator rather than the Windows-specific ``\``
        directory separator (e.g., by calling ``'\\'.join()`` rather than either
        :func:`betse.util.path.paths.join` or :func:`os.path.join`).

    Returns
    ----------
    str
        Absolute path of this resource.
    '''

    return pkg_resources.resource_filename(module_name, pathname)
