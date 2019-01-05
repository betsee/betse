#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level package-specific facilities.

All functions defined by this submodule accept at least a previously imported
package object; most also accept the fully-qualified name of a package.

See Also
----------
:mod:`betse.util.py.module.pyname`
    Related submodule whose functions accept only fully-qualified names.
'''

# ....................{ IMPORTS                           }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: To raise human-readable exceptions on missing mandatory dependencies,
# the top-level of this module may import *ONLY* from packages guaranteed to
# exist at installation time -- which typically means *ONLY* BETSE packages and
# stock Python packages.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# from betse.util.io.log import logs
from betse.util.type.types import type_check, ModuleType  #, ModuleOrStrTypes

# ....................{ GETTERS                           }....................
@type_check
def get_object_type_package_root(obj: object) -> ModuleType:
    '''
    **Root package** (i.e., topmost package whose package name contains no
    ``.`` delimiters) transitively defining the class of the passed object.

    Parameters
    ----------
    obj : object
        Object to retrieve the root package name of.

    Returns
    ----------
    ModuleType
        Root package transitively defining this object's class.
    '''

    # Avoid circular import dependencies.
    from betse.util.py.module import pymodname

    # Fully-qualified name of the root package of this object's class.
    class_package_root_name = get_object_type_package_root_name(obj)

    # Return the previously imported instance of this root package.
    return pymodname.import_module(module_name=class_package_root_name)


@type_check
def get_object_type_package_root_name(obj: object) -> str:
    '''
    Name of the **root package** (i.e., topmost package whose package name
    contains no ``.`` delimiters) transitively defining the class of the passed
    object if any *or* raise an exception otherwise (i.e., if that class is
    defined by no module, as is the case for standard C-based classes).

    Design
    ----------
    The name of this function is intentionally suffixed by neither
    ``_qualified`` nor ``_unqualified``. Since root package names contain no
    ``.`` delimiters (e.g., :mod:`betse`), the fully-qualified and unqualified
    names for any root package are necessarily identical.

    Parameters
    ----------
    obj : object
        Object to retrieve the root package name of.

    Returns
    ----------
    str
        Name of the root package transitively defining this object's class.

    Raises
    ----------
    BetseTypeException
        If this object's class is *not* defined by a module (e.g., :mod:`str`).
        See the :func:`betse.util.type.cls.classes.get_module_name_qualified`
        function for further details.
    '''

    # Avoid circular import dependencies.
    from betse.util.type.cls import classes
    from betse.util.type.obj import objects
    from betse.util.type.text import strs

    # Class of this object.
    cls = objects.get_class(obj)

    # Fully-qualified name of the module defining this class if any *OR* raise
    # the "BetseTypeException" otherwise.
    module_name = classes.get_module_name_qualified(cls)

    # Name of the root package defining this module. Note that:
    #
    # * If this module name is itself this package name, this module name
    #   contains no "." delimiter, in which case this call correctly returns
    #   this module name as is.
    # * Else, this module name contains at least one "." delimiter, in which
    #   case this call correctly returns the prefix of this module name
    #   preceding the first such ".".
    package_root_name = strs.get_prefix_preceding_char_or_text(
        text=module_name, char='.')

    # Return this name.
    return package_root_name
