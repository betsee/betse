# -*- mode: python -*-
# ====================[ betse.core.spec                    ]====================
#
# --------------------( SYNOPSIS                           )--------------------
# Low-level dependency of all high-level PyInstaller spec files freezing the CLI
# script "betse" into an executable file or directory.
#
# --------------------( CAVEATS                            )--------------------
# This file expects "betse" to have been editably installed via the custom
# setuptools command "symlink" rather than non-editably installed via the core
# setuptools command "install". For details on how to adjust this file to
# support the latter, see:
#     https://github.com/pyinstaller/pyinstaller/wiki/Recipe-Setuptools-Entry-Point

#FIXME: The produced executables are rather large, at the moment. The reason why
#appears to be indiscriminate inclusion of *ALL* possible matplotlib backends.
#We should probably correct this explicitly excluding the largest such backends
#-- particularly "PyQt4" and subsidiaries. For details on how to exclude modules
#and libraries, see this great stackoverflow answer:
#    https://stackoverflow.com/questions/4890159/python-excluding-modules-pyinstaller/17595149#17595149

#FIXME: When we begin work on freezing "betse-qt", note that the correct
#solution is almost certainly to leverage a so-called "multipackage bundle."
#This prevents us from duplicating embedded libraries in "betse-qt" that are
#already bundled with "betse", dramatically reducing executable filesize. See:
#    http://pythonhosted.org/PyInstaller/#id67

# ....................{ IMPORTS                            }....................
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: PyInstaller spec files are *NOT* valid Python code. Specifically,
# attempting to import modules in such files does nothing at beast and silently
# corrupts the current global namespace at worst, preventing further processing.
# The following modules are "imported" only to avoid lint errors during editing.
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

from os import path
import os, platform

# ....................{ CONSTANTS ~ os                     }....................
# True if the current operating system is Linux.
IS_OS_LINUX = platform.system() == 'Linux'

# True if the current operating system is OS X.
IS_OS_OS_X = platform.system() == 'Darwin'

# True if the current operating system is Windows.
IS_OS_WINDOWS = platform.system() == 'Windows'

# ....................{ CONSTANTS ~ path                   }....................
# Absolute path of the current directory.
CURRENT_DIRNAME = os.getcwd()

# Absolute path of the root Python module running such script. For portability,
# such path references a module in the Python package tree for this application
# rather than a script wrapper previously installed by setuptools (e.g.,
# "/usr/bin/betse").
MODULE_ROOT_FILENAME = path.join(
    CURRENT_DIRNAME, 'betse', 'cli', '__main__.py')

# Absolute path of the top-level directory containing all non-Python data files
# to be bundled with such executable.
DATA_ROOT_DIRNAME = path.join(CURRENT_DIRNAME, 'betse', 'data')

# ....................{ CONSTANTS ~ exe                    }....................
# Basename of the executable to be output.
EXE_basename = 'betse'
if IS_OS_WINDOWS:
    EXE_basename += '.exe'

# Keyword arguments to be subsequently passed to EXE() by parent files.
EXE_options = {
    'name': EXE_basename,
    'debug': False,
    'strip': None,
    'upx': True,
    'console': True,

    #FIXME: Uncomment *AFTER* adding such file. Note that this is probably
    #Linux- and Windows-specific, however; OS X appears to require a completely
    #different icon file format. (Naturally. Why not?)
    # 'icon': path.join(APPPATH, "..", "res", "Icon.ico"),
}

# ....................{ MAIN                               }....................
# Analyze such top-level Python module.
a = Analysis(
    [MODULE_ROOT_FILENAME],
    pathex = [CURRENT_DIRNAME],

    # List of the names of all imported and hence required modules *NOT*
    # automatically detectable by PyInstaller.
    hiddenimports = [],

    # List of the relative paths of all directories containing non-runtime
    # PyInstaller hooks.
    hookspath = ['freeze/hooks'],

    # List of the relative paths of all directories containing runtime
    # PyInstaller hooks.
    runtime_hooks = None,
)

# ....................{ MAIN ~ datas                       }....................
# Record all non-Python data files to be bundled with such executable. For
# further details, see http://pythonhosted.org/PyInstaller/#id38.
a.datas += Tree(
    DATA_ROOT_DIRNAME,

    # Relative path of the output top-level directory containing such files.
    prefix = 'data',
)

# ....................{ MAIN ~ pyz                         }....................
# Create an intermediate archive containing only Python modules and scripts.
pyz = PYZ(a.pure)

# --------------------( WASTELANDS                         )--------------------
        # 'scipy.sparse.sparsetools._csr',
# ....................{ MAIN ~ binaries : scipy            }....................
# Include all `scipy`-specific shared libraries.

# import PyInstaller.log as logging
# from PyInstaller.utils.hooks import hookutils
# from PyInstaller.utils import misc
#
# #FUXME: The following code is unused. PyInstaller appears to be unable to
# #tolerate functions in spec files. (We're not kidding.) That said, the following
# #
#
# #FUXME: The following should probably be pushed to the PyInstaller codebase
# #itself as official "scipy" hooks. For a roughly analogous hook, see
# #"PyInstaller/hooks/hook-PyQt4.QtGui.py". Perhaps add such hooks by adding a
# #similar hook() method to a new "PyInstaller/hooks/hook-scipy.py" script.
# # hook_scipy(a)
#
# def hook_scipy(a):
#     '''
#     Include all `scipy`-specific shared libraries.
#     '''
#     a.binaries.extend(collect_package_binaries('scipy'))
#
# #FUXME: This is sufficiently useful that it should probably be generalized
# #into a new "hookutils" utility function.
# def collect_package_binaries(package_name):
#     '''
#     Get a `mod.binaries`-formatted list of 3-tuples installing all dynamic
#     libraries recursively found in the source directory containing the package
#     with the passed name to corresponding subdirectories of a target directory
#     with the same name.
#     '''
#     # Absolute path of such package's "__init__.py" script.
#     package_init_file = hookutils.eval_statement('''
#         import {package_name}
#         print({package_name}.__file__)
#     '''.format(package_name = package_name))
#
#     return collect_binaries(
#         os.path.dirname(package_init_file), package_name)
#
# #FUXME: This is sufficiently useful that it should probably be generalized
# #into a new "hookutils" utility function.
# #FUXME: qt5_qml_plugins_binaries() should be refactored to call this instead.
#
# def collect_binaries(src_root_dir, trg_root_dir):
#     '''
#     Get a `mod.binaries`-formatted list of 3-tuples installing all dynamic
#     libraries recursively found in the passed source directory to corresponding
#     subdirectories of the passed target directory.
#
#     Such source directory should be absolute; such target directory should be
#     relative. For example, to list all shared libraries for `scipy`:
#
#         >>> from PyInstaller.utils.hooks.hookutils import collect_binaries
#         >>> collect_binaries(
#         ...     '/usr/lib64/python3.3/site-packages/scipy', 'scipy')
#     '''
#     logger = logging.getLogger(__name__)
#     binaries = []
#     src_dlls = misc.dlls_in_subdirs(src_root_dir)
#
#     for src_dll in src_dlls:
#         src_dll_relative = os.path.relpath(src_dll, src_root_dir)
#         src_dll_relative_dir, src_dll_base = os.path.split(src_dll_relative)
#         trg_dll = os.path.join(trg_root_dir, src_dll_relative_dir, src_dll_base)
#
#         binaries.append((trg_dll, src_dll, 'BINARY'))
#         logger.info(
#             'Dynamic library "%s" installing to "%s".', src_dll, trg_dll)
#
#     return binaries

    # import scipy
    # a.binaries.extend(
    #     collect_binaries(os.path.dirname(scipy.__file__), 'scipy'))
    # Absolute paths of shared libraries in such directory.
        # Relative path of such library relative to such top-level directory.
            #FUXME: Just use "src_dll_relative", no?
    # # Relative path of the top-level "site-packages/scipy" directory.
    # trg_root_dir = 'scipy'
    #
    # # Absolute path of the top-level "site-packages/scipy" directory.
    # src_root_dir = os.path.dirname(scipy.__file__)

    # assert isinstance(module_object, module), \
    #     '"{}" not a module.'.format(module_object)
    # List of TOC-style tuples describing such libraries.
    # scipy_binaries = []

    # excludes = ['PyQt4', 'PyQt4.QtCore', 'PyQt4.QtGui'],
# PyInstaller-specific collection of such files, which may be passed as is to
# the EXE(), COLLECT(), and BUNDLE() functions. For further details, see:
    # Absolute path of the input top-level directory containing such files.
    # root = path.join(CURRENT_DIRNAME, 'betse', 'data'),
# a = Analysis(['/usr/bin/betse'],
# Absolute path of the root Python package for this application.
# PACKAGE_ROOT_DIRNAME = path.join(os.getcwd(), ')
