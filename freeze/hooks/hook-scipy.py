#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Include all `scipy`-specific shared libraries.
'''

# ....................{ IMPORTS                            }....................
from PyInstaller.utils.hooks import hookutils
from PyInstaller.utils import misc
import os

import PyInstaller.log as logging
logger = logging.getLogger(__name__)

#FIXME: The following should probably be pushed to the PyInstaller codebase
#itself as official "scipy" hooks. For a roughly analogous hook, see
#"PyInstaller/hooks/hook-PyQt4.QtGui.py". Perhaps add such hooks by adding a
#similar hook() method to a new "PyInstaller/hooks/hook-scipy.py" script.
# hook_scipy(a)

def hook(mod):
    '''
    Include all `scipy`-specific shared libraries.
    '''
    mod.add_binary(collect_package_binaries('scipy'))
    return mod

#FIXME: This is sufficiently useful that it should probably be generalized
#into a new "hookutils" utility function.

def collect_package_binaries(package_name):
    '''
    Get a `mod.binaries`-formatted list of 3-tuples installing all dynamic
    libraries recursively found in the source directory containing the package
    with the passed name to corresponding subdirectories of a target directory
    with the same name.
    '''
    # Absolute path of such package's "__init__.py" script.
    package_init_file = hookutils.exec_statement(
        'import {package_name}; print({package_name}.__file__)'.format(
            package_name = package_name))

    return collect_binaries(
        os.path.dirname(package_init_file), package_name)

#FIXME: This is sufficiently useful that it should probably be generalized
#into a new "hookutils" utility function.
#FIXME: qt5_qml_plugins_binaries() should be refactored to call this instead.

def collect_binaries(src_root_dir, trg_root_dir):
    '''
    Get a `mod.binaries`-formatted list of 3-tuples installing all dynamic
    libraries recursively found in the passed source directory to corresponding
    subdirectories of the passed target directory.

    Such source directory should be absolute; such target directory should be
    relative. For example, to list all shared libraries for `scipy`:

        >>> from PyInstaller.utils.hooks.hookutils import collect_binaries
        >>> collect_binaries(
        ...     '/usr/lib64/python3.3/site-packages/scipy', 'scipy')
    '''
    binaries = []
    src_dlls = misc.dlls_in_subdirs(src_root_dir)

    for src_dll in src_dlls:
        src_dll_relative = os.path.relpath(src_dll, src_root_dir)
        src_dll_relative_dir, src_dll_base = os.path.split(src_dll_relative)
        trg_dll = os.path.join(trg_root_dir, src_dll_relative_dir, src_dll_base)

        binaries.append((trg_dll, src_dll, 'BINARY'))
        logger.debug(
            'Dynamic library "%s" installing to "%s".', src_dll, trg_dll)

    return binaries

# --------------------( WASTELANDS                         )--------------------
    # print('FILE: '+package_init_file)
#FUXME: The following code is unused. PyInstaller appears to be unable to
#tolerate functions in spec files. (We're not kidding.) That said, the following
