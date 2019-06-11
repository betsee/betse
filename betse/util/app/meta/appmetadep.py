#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level **application dependency metadata** (i.e., lists of version-pinned
dependencies synopsizing application requirements) functionality.
'''

# ....................{ IMPORTS                           }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: To avoid race conditions during setuptools-based installation, this
# module may import *ONLY* from modules guaranteed to exist at the start of
# installation. This includes all standard Python and application modules but
# *NOT* third-party dependencies, which if currently uninstalled will only be
# installed at some later time in the installation. Likewise, to avoid circular
# import dependencies, the top-level of this module should avoid importing
# application modules where feasible.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# from betse.util.io.log import logs
from betse.util.type.types import (
    type_check, MappingOrNoneTypes, ModuleType, StrOrNoneTypes)

# ....................{ MAKERS                            }....................
# @type_check
# def make_module(
#     module_name: str,
#     module_doc: StrOrNoneTypes = None,
#     module_attr_name_to_value: MappingOrNoneTypes = None,
#     is_importable: bool = False,
# ) -> ModuleType:
#     '''
#     Dynamically create and return a new module with the passed name and
#     optional passed docstring and attributes, optionally added to the standard
#     :attr:`sys.modules` dictionary for external importation elsewhere.
#
#     Parameters
#     ----------
#     module_name : str
#         Fully-qualified name of the module to be created.
#     module_doc: StrOrNoneTypes
#         **Docstring** (i.e., human-readable documentation in reStructuredText
#         (reST) format) to be associated with this module. Defaults to ``None``,
#         in which case this module is undocumented by default.
#     module_attr_name_to_value : MappingOrNoneTypes
#         Dictionary mapping from the name to value of each module-scoped
#         attribute to be declared in this module. Defaults to ``None``, in which
#         case this module is empty (i.e., contains *no* attributes) by default.
#     is_importable : bool
#         ``True`` only if external callers may trivially import this module via
#         this module name elsewhere in this codebase. Defaults to ``False`` for
#         safety, in which case the returned object is the *only* initial
#         reference to this module.
#
#     Raises
#     ----------
#     ????
#         If the ``is_importable`` parameter is ``True`` *and* one or more parent
#         packages of this module are unimportable.
#
#     See Also
#     ----------
#     https://stackoverflow.com/questions/2931950/dynamic-module-creation
#         StackOverflow *question* strongly inspiring this implementation. Note
#         that, against all expectations, this question is substantially more
#         informative than its answers.
#     '''
#
#     pass
